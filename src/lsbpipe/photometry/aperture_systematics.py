"""Stage 04.3: read-only conditional systematics and blank-sky diagnostics.

All results describe the CURRENT native coadd, draft masks and fixed aperture
hypotheses. They are not adopted sky corrections, missing-light reconstructions,
or calibrated uncertainties of a total galaxy magnitude.
"""
from __future__ import annotations

import json
import math

import numpy as np

from .aperture_flux import (
    ApertureGeometryICRS, ApertureProvenance, ApertureSkyScenario,
    measure_native_aperture, project_icrs_ellipse_to_native,
)
from .contracts import NativeBandFrame


_SCHEMA = "desidera.photometry.aperture-systematics/1"
_BLANK_SCHEMA = "desidera.photometry.blank-annulus-covariance/1"


def _scene_masks(frame, source_mask, manual_exclude):
    """Snapshot user-controlled masks and retain immutable hard invalidity."""
    arrays = []
    for name, raw in (("source", source_mask), ("manual", manual_exclude)):
        array = np.asarray(raw)
        if array.dtype.kind != "b" or array.shape != frame.shape:
            raise ValueError(name + " mask must be native-grid Boolean")
        arrays.append(np.array(array, dtype=bool, copy=True))
    return arrays[0], arrays[1]


def _ring_aperture(x, y, projected, radii, index):
    """Exact annulus (including 0..r0), translated without resampling."""
    from photutils.aperture import EllipticalAnnulus, EllipticalAperture

    scale = radii[index] / radii[-1]
    if index == 0:
        return EllipticalAperture(
            (x, y), projected.a_pix * scale,
            projected.b_pix * scale,
            theta=projected.theta_rad_pixel_ccw)
    return EllipticalAnnulus(
        (x, y),
        a_in=projected.a_pix * radii[index - 1] / radii[-1],
        a_out=projected.a_pix * scale,
        b_in=projected.b_pix * radii[index - 1] / radii[-1],
        b_out=projected.b_pix * scale,
        theta=projected.theta_rad_pixel_ccw)


def estimate_blank_annulus_covariance(
        frame: NativeBandFrame, outer: ApertureGeometryICRS,
        radii_arcsec, *, effective_source_mask, manual_exclude,
        min_placements=12, max_placements=48):
    """Empirical covariance of signed fluxes in non-overlapping blank rings.

    Candidates have the same native-pixel shape as the galaxy rings; they
    cannot overlap each other or the target aperture. A placement with any
    masked/invalid geometric area is rejected. No assumption of pixel-white
    noise, no sky smoothing, and no inverse covariance are introduced.

    The target exclusion extends only to its specified outer ellipse;
    undetected galaxy outskirts, background gradients and residual
    contaminants remain possible biases. Hence even successful sampling
    remains DIAGNOSTIC and cannot certify a science covariance.
    """
    if not isinstance(frame, NativeBandFrame) or not isinstance(
            outer, ApertureGeometryICRS):
        raise TypeError("Validated native frame and ICRS geometry required")
    radii = tuple(float(v) for v in radii_arcsec)
    if (not 2 <= len(radii) <= 40 or
            any(not math.isfinite(v) or v <= 0 for v in radii) or
            any(b <= a for a, b in zip(radii[:-1], radii[1:])) or
            not math.isclose(radii[-1], outer.semimajor_arcsec,
                             rel_tol=1e-9, abs_tol=1e-10)):
        raise ValueError("Increasing radii must end at the outer ellipse")
    if (isinstance(min_placements, bool) or
            not isinstance(min_placements, int) or
            min_placements < 3 or min_placements > 256 or
            isinstance(max_placements, bool) or
            not isinstance(max_placements, int) or
            max_placements < min_placements or max_placements > 256):
        raise ValueError("Invalid bounded blank-aperture sampling limits")

    source, manual = _scene_masks(
        frame, effective_source_mask, manual_exclude)
    invalid = ~frame.valid | source | manual
    projected = project_icrs_ellipse_to_native(outer, frame.wcs)
    margin = int(math.ceil(projected.a_pix + 3))
    step = max(2 * margin + 2, 8)
    ny, nx = frame.shape
    centres, samples = [], []
    examined, rejected = 0, 0
    required = max(min_placements, 2 * len(radii) + 2)
    for row, y in enumerate(range(margin, ny - margin, step)):
        start = margin + (step // 2 if row % 2 else 0)
        for x in range(start, nx - margin, step):
            if len(centres) >= max_placements:
                break
            examined += 1
            # At least the full outer aperture separates each sky placement
            # from the putative target. Extended LSB light is NOT excluded.
            if math.hypot(
                    x - projected.x_center,
                    y - projected.y_center) <= (
                    2. * projected.a_pix + 3.):
                rejected += 1
                continue
            vector = []
            for i in range(len(radii)):
                mask = _ring_aperture(x, y, projected, radii, i).to_mask(
                    method="exact")
                large, small = mask.get_overlap_slices(frame.shape)
                if large is None:
                    vector = []
                    break
                weights = np.asarray(mask.data[small], dtype=float)
                if (not np.isfinite(weights).all() or
                        np.any((weights < -1e-10) | (weights > 1 + 1e-10))):
                    raise ValueError("Invalid blank-aperture exact weights")
                weights = np.clip(weights, 0., 1.)
                # This checks all partially covered boundary pixels too.
                bad_area = float(np.sum(weights * invalid[large]))
                area = float(np.sum(weights))
                if area <= 0. or bad_area > 1e-6 * area:
                    vector = []
                    break
                vector.append(float(np.sum(
                    weights * frame.image_nMgy[large],
                    dtype=np.float64)))
            if len(vector) != len(radii):
                rejected += 1
                continue
            centres.append([int(x), int(y)])
            samples.append(vector)
        if len(centres) >= max_placements:
            break

    report = {
        "schema": _BLANK_SCHEMA,
        "dataset_id": frame.dataset_id,
        "radii_arcsec": list(radii),
        "placement_policy": (
            "DETERMINISTIC_NONOVERLAPPING_NATIVE_ELLIPTICAL_RINGS_V1"),
        "n_examined": int(examined),
        "n_rejected": int(rejected),
        "n_accepted": len(centres),
        "n_required": int(required),
        "centres_native_xy": centres,
        "annulus_mean_blank_flux_nMgy": None,
        "annulus_covariance_nMgy2": None,
        "annulus_sigma_empirical_nMgy": None,
        "cumulative_covariance_nMgy2": None,
        "cumulative_sigma_empirical_nMgy": None,
        "status": "INSUFFICIENT_INDEPENDENT_BLANK_POSITIONS",
        "covariance_calibrated_for_science": False,
        "warnings": [
            "TARGET_LIGHT_OUTSIDE_TESTED_ELLIPSE_NOT_EXCLUDED",
            "UNDETECTED_SOURCES_OR_BACKGROUND_GRADIENTS_MAY_REMAIN",
            "BLANK_POSITIONS_FROM_ONE_COADD_NOT_AN_INDEPENDENT_EXPOSURE_TEST",
        ],
    }
    if len(centres) >= required:
        observations = np.asarray(samples, dtype=float)
        covariance = np.atleast_2d(np.cov(observations, rowvar=False, ddof=1))
        diagonal = np.diag(covariance)
        if (np.isfinite(covariance).all() and
                np.all(diagonal >= 0) and np.any(diagonal > 0)):
            report["status"] = "EMPIRICAL_BLANK_COVARIANCE_DIAGNOSTIC_ONLY"
            report["annulus_mean_blank_flux_nMgy"] = (
                np.mean(observations, axis=0).tolist())
            report["annulus_covariance_nMgy2"] = covariance.tolist()
            report["annulus_sigma_empirical_nMgy"] = np.sqrt(
                diagonal).tolist()
            # Nested apertures share their inner rings: C_CoG = T C_ring T^T.
            # No independent-error quadrature and no covariance inversion.
            transfer = np.tril(np.ones((len(radii), len(radii))))
            cumulative = transfer @ covariance @ transfer.T
            report["cumulative_covariance_nMgy2"] = cumulative.tolist()
            report["cumulative_sigma_empirical_nMgy"] = np.sqrt(
                np.maximum(0., np.diag(cumulative))).tolist()
        else:
            report["status"] = "DEGENERATE_OR_NONFINITE_BLANK_VARIANCE"
    json.dumps(report, allow_nan=False)
    return report


def diagnose_aperture_systematics(
        frame: NativeBandFrame, *, anchor_icrs, provenance: ApertureProvenance,
        effective_source_mask, manual_exclude, curves):
    """Conditional sensitivity at the largest common, observed aperture.

    - sky response = exact derivative of observed valid-area flux to a
      residual constant or tilted native-pixel sky at the target centre;
    - source-mask counterfactual = recompute with ONLY provisional catalogue
      exclusions disabled; never restores hard-invalid or manual-excluded data;
    - geometry span = difference of distinct observed regions, NOT 1-sigma;
    - blank covariance = non-overlapping empty-field diagnostic, if feasible.
    """
    if not isinstance(frame, NativeBandFrame) or not isinstance(
            provenance, ApertureProvenance):
        raise TypeError("Validated frame and project provenance required")
    source, manual = _scene_masks(
        frame, effective_source_mask, manual_exclude)
    if not curves:
        raise ValueError("At least one frozen curve required")
    first = curves[0]
    last = first["cumulative"][-1]
    radius = float(first["radii_arcsec"][-1])
    hypothesis = first["geometry_hypothesis"]
    outer = ApertureGeometryICRS(
        ra_deg=float(anchor_icrs[0]), dec_deg=float(anchor_icrs[1]),
        semimajor_arcsec=radius, axis_ratio=float(hypothesis["q"]),
        pa_deg_east_of_north=float(hypothesis["pa_deg"]))
    projected = project_icrs_ellipse_to_native(outer, frame.wcs)

    from photutils.aperture import EllipticalAperture
    aperture = EllipticalAperture(
        (projected.x_center, projected.y_center),
        projected.a_pix, projected.b_pix,
        theta=projected.theta_rad_pixel_ccw)
    mask = aperture.to_mask(method="exact")
    large, small = mask.get_overlap_slices(frame.shape)
    area = dx_sum = dy_sum = 0.
    if large is not None:
        weights = np.clip(np.asarray(mask.data[small], dtype=float), 0., 1.)
        usable = frame.valid[large] & ~source[large] & ~manual[large]
        yy, xx = np.indices(weights.shape, dtype=float)
        xx += large[1].start
        yy += large[0].start
        effective = weights * usable
        area = float(np.sum(effective))
        dx_sum = float(np.sum(effective * (xx - projected.x_center)))
        dy_sum = float(np.sum(effective * (yy - projected.y_center)))
    if not math.isclose(
            area, last["area_native_pixel2"]["valid_native_pixel2"],
            rel_tol=1e-10, abs_tol=1e-6):
        raise ValueError("Stage 04.3 masks disagree with frozen Stage 04.2")

    relaxed = measure_native_aperture(
        frame, outer, provenance,
        effective_source_mask=np.zeros(frame.shape, dtype=bool),
        manual_exclude=manual, sky=ApertureSkyScenario())
    original_flux = last["signed_valid_flux_nMgy"]
    relaxed_flux = relaxed.signed_valid_flux_nMgy
    candidate_fluxes = [
        float(curve["cumulative"][-1]["signed_valid_flux_nMgy"])
        for curve in curves
        if (math.isclose(float(curve["radii_arcsec"][-1]), radius,
                         rel_tol=1e-10, abs_tol=1e-10) and
            curve["cumulative"][-1]["signed_valid_flux_nMgy"] is not None)
    ]
    blanks = estimate_blank_annulus_covariance(
        frame, outer, first["radii_arcsec"],
        effective_source_mask=source, manual_exclude=manual)
    report = {
        "schema": _SCHEMA,
        "project_id": provenance.project_id,
        "dataset_id": frame.dataset_id,
        "draft_sha256": provenance.scene_draft_sha256,
        "catalogue_sha256": provenance.catalogue_sha256,
        "source_policy_id": provenance.source_policy_id,
        "baseline_geometry": hypothesis,
        "reference_semimajor_arcsec": radius,
        "reference_valid_flux_nMgy": original_flux,
        "sky_response": {
            "reference_native_xy": [
                projected.x_center, projected.y_center],
            "dflux_db0_native_pixel2": -area,
            "dflux_dbx_native_pixel3": -dx_sum,
            "dflux_dby_native_pixel3": -dy_sum,
            "sky_level_uncertainty_assumed": None,
            "sky_plane_uncertainty_assumed": None,
            "sky_science_adopted": False,
        },
        "catalogue_mask_counterfactual": {
            "flux_with_provisional_source_mask_nMgy": original_flux,
            "flux_without_provisional_source_mask_nMgy": relaxed_flux,
            "conditional_delta_nMgy": (
                float(relaxed_flux - original_flux)
                if relaxed_flux is not None and original_flux is not None
                else None),
            "manual_and_hard_exclusions_retained": True,
            "delta_is_not_missing_galaxy_light": True,
        },
        "geometry_hypothesis_span": {
            "n_valid_hypotheses": len(candidate_fluxes),
            "min_flux_nMgy": min(candidate_fluxes) if candidate_fluxes else None,
            "max_flux_nMgy": max(candidate_fluxes) if candidate_fluxes else None,
            "span_nMgy": (
                max(candidate_fluxes) - min(candidate_fluxes)
                if candidate_fluxes else None),
            "not_a_statistical_uncertainty": True,
        },
        "blank_apertures": blanks,
        "adopted_for_science": False,
        "covariance_calibrated_for_science": False,
        "status": "CONDITIONAL_SYSTEMATICS_DIAGNOSTIC_ONLY",
        "warnings": [
            "SKY_RESPONSE_IS_A_DERIVATIVE_NOT_A_SKY_ERROR_BAR",
            "RELAXED_CATALOGUE_MASK_INCLUDES_CONTAMINANTS_AND_TARGET_MEMBERS",
            "GEOMETRY_SPAN_COMPARES_DIFFERENT_SKY_FOOTPRINTS",
            "NO_TOTAL_FLUX_OR_PSF_CORRECTION_INFERRED",
        ],
    }
    json.dumps(report, allow_nan=False)
    return report
