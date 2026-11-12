"""Stage 04.3: read-only conditional systematics and blank-sky diagnostics.

All results describe the CURRENT native coadd, draft masks and fixed aperture
hypotheses. They are not adopted sky corrections, missing-light reconstructions,
or calibrated uncertainties of a total galaxy magnitude.
"""
from __future__ import annotations

from dataclasses import asdict
from hashlib import sha256
import json
import math
from bisect import bisect_right

import numpy as np

from .aperture_flux import (
    ApertureGeometryICRS, ApertureProvenance, ApertureSkyScenario,
    measure_native_aperture, project_icrs_ellipse_to_native,
)
from .contracts import NativeBandFrame
from lsbpipe.scene import grid_fingerprint


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


def _bounded_quadrant_blank_candidates(shape, margin, step, *, seed=20261112):
    """Predeclare a staggered four-quadrant lattice with O(rows) memory.

    A fully enumerated `(x, y)` list can be millions of entries on a
    large native survey image, defeating the accepted/evaluation budgets.
    Instead, each image row contributes at most TWO contiguous
    arithmetic-progression blocks. A seeded coprime stride visits each
    block's flattened quadrant index exactly once without constructing
    any list of pixel centres. No fluxes or masks enter candidate order.

    Returns (candidate_counts_by_quadrant, candidate_iterator).
    """
    ny, nx = shape
    if not (isinstance(ny, (int, np.integer)) and
            isinstance(nx, (int, np.integer)) and ny > 0 and nx > 0
            and isinstance(margin, int) and margin >= 0
            and isinstance(step, int) and step > 0):
        raise ValueError("Invalid bounded blank candidate geometry")
    blocks = [[], [], [], []]
    totals = [0, 0, 0, 0]
    for row, y in enumerate(range(margin, ny - margin, step)):
        x_first = margin + (step // 2 if row % 2 else 0)
        n_cols = len(range(x_first, nx - margin, step))
        if n_cols == 0:
            continue
        # The two halves are defined solely in original native pixels.
        n_left = max(0, min(
            n_cols, math.ceil((nx / 2. - x_first) / step)))
        choices = (
            (0, x_first, n_left),
            (1, x_first + n_left * step, n_cols - n_left),
        )
        for side, first, count in choices:
            if count == 0:
                continue
            quadrant = side + 2 * int(y >= ny / 2.)
            totals[quadrant] += count
            # prefix end, original y row, first x column
            blocks[quadrant].append(
                (totals[quadrant], int(y), int(first)))
    # Seeded permutation of a quadrant's flattenable integer indices,
    # no permutation array and no sampling with replacement.
    rng = np.random.default_rng(seed)
    offsets, strides, ends = [], [], []
    for quadrant, count in enumerate(totals):
        if count == 0:
            offsets.append(0)
            strides.append(1)
        else:
            offsets.append(int(rng.integers(0, count)))
            stride = int(rng.integers(1, count + 1))
            while math.gcd(stride, count) != 1:
                stride = 1 if stride == count else stride + 1
            strides.append(stride)
        ends.append([record[0] for record in blocks[quadrant]])

    def iterator():
        for i in range(max(totals, default=0)):
            for quadrant, count in enumerate(totals):
                if i >= count:
                    continue
                flat_index = (offsets[quadrant] +
                              i * strides[quadrant]) % count
                block_index = bisect_right(ends[quadrant], flat_index)
                preceding = (ends[quadrant][block_index - 1]
                             if block_index else 0)
                _, y, x_first = blocks[quadrant][block_index]
                yield (quadrant,
                       x_first + step * (flat_index - preceding),
                       y)
    return totals, iterator()


def estimate_blank_annulus_covariance(
        frame: NativeBandFrame, outer: ApertureGeometryICRS,
        radii_arcsec, *, effective_source_mask, manual_exclude,
        min_placements=12, max_placements=None,
        full_target_exclusion_outer: ApertureGeometryICRS | None = None):
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
            (max_placements is not None and (
                isinstance(max_placements, bool) or
                not isinstance(max_placements, int) or
                max_placements < min_placements or max_placements > 256))):
        raise ValueError("Invalid bounded blank-aperture sampling limits")

    source, manual = _scene_masks(
        frame, effective_source_mask, manual_exclude)
    invalid = ~frame.valid | source | manual
    projected = project_icrs_ellipse_to_native(outer, frame.wcs)
    # A smaller-scale blank probe must STILL avoid the full galaxy aperture.
    # The centre-distance test is deliberately conservative for ellipses.
    if full_target_exclusion_outer is None:
        full_radius_pix = projected.a_pix
    else:
        if not isinstance(full_target_exclusion_outer, ApertureGeometryICRS):
            raise TypeError("Full target exclusion must be ICRS ellipse")
        full_projection = project_icrs_ellipse_to_native(
            full_target_exclusion_outer, frame.wcs)
        if (math.hypot(
                full_projection.x_center - projected.x_center,
                full_projection.y_center - projected.y_center) > 1e-3):
            raise ValueError("Full target exclusion must share project anchor")
        full_radius_pix = full_projection.a_pix
        if full_radius_pix < projected.a_pix:
            raise ValueError("Full target exclusion cannot be smaller than blank rings")
    exclusion_distance_pix = full_radius_pix + projected.a_pix + 3.
    margin = int(math.ceil(projected.a_pix + 3))
    step = max(2 * margin + 2, 8)
    ny, nx = frame.shape
    centres, samples = [], []
    examined, rejected = 0, 0
    required = max(min_placements, 2 * len(radii) + 2)
    # A heavily masked large CCD must not trigger an unbounded search.
    # For long radial vectors, a hardcoded 48-position default can be
    # LESS than the required 2*n_rings+2. Auto-scale the default
    # sampling budget, but never silently override an explicit user cap.
    requested_max_placements = max_placements
    if max_placements is None:
        max_placements = min(256, max(48, required))
    candidate_budget = max(256, 32 * max_placements)
    # Enumerate at most two integer-progression blocks per native row.
    # Never allocate all possible (x, y) positions: an image can have
    # millions of lattice nodes even when at most 256 blanks are used.
    available_by_quadrant, candidate_positions = (
        _bounded_quadrant_blank_candidates(frame.shape, margin, step))
    accepted_by_quadrant = [0, 0, 0, 0]
    for quadrant, x, y in candidate_positions:
        if len(centres) >= max_placements or examined >= candidate_budget:
            break
        examined += 1
        # At least the full outer aperture separates each sky placement
        # from the putative target. Extended LSB light is NOT excluded.
        if math.hypot(
                x - projected.x_center,
                y - projected.y_center) <= exclusion_distance_pix:
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
            # Reject every positive-area contact with a bad pixel. A
            # tiny tolerated masked fraction is not an unmasked blank.
            # Crucially, even a zero fractional weight times NaN is NaN
            # in NumPy: never sum the entire rectangular mask bounding
            # box, which may contain invalid pixels OUTSIDE the ring.
            support = weights > 0.
            if (not np.any(support) or
                    np.any(invalid[large] & support)):
                vector = []
                break
            vector.append(float(np.sum(
                weights[support] * frame.image_nMgy[large][support],
                dtype=np.float64)))
        if len(vector) != len(radii):
            rejected += 1
            continue
        centres.append([int(x), int(y)])
        accepted_by_quadrant[quadrant] += 1
        samples.append(vector)

    report = {
        "schema": _BLANK_SCHEMA,
        "dataset_id": frame.dataset_id,
        "radii_arcsec": list(radii),
        "placement_policy": (
            "STRATIFIED_QUADRANT_BOUNDED_MEMORY_NATIVE_RINGS_V3"),
        "candidate_order_independent_of_image_brightness": True,
        "candidate_order_random_seed": 20261112,
        "candidate_storage_policy": "ROW_BLOCKS_AND_COPRIME_INDEX_STRIDE_NO_CENTER_ARRAY",
        "candidate_grid_nodes_by_quadrant": available_by_quadrant,
        "accepted_centres_by_quadrant": accepted_by_quadrant,
        "n_quadrants_with_accepted_blanks": sum(
            n > 0 for n in accepted_by_quadrant),
        "n_examined": int(examined),
        "n_rejected": int(rejected),
        "n_accepted": len(centres),
        "n_required": int(required),
        "n_max_placements": int(max_placements),
        "placement_cap_explicit": requested_max_placements is not None,
        "n_candidate_budget": int(candidate_budget),
        "target_exclusion_radius_native_pix": float(full_radius_pix),
        "blank_outer_radius_native_pix": float(projected.a_pix),
        "target_plus_blank_exclusion_distance_native_pix": float(
            exclusion_distance_pix),
        "candidate_search_budget_exhausted": examined >= candidate_budget,
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
            "NONOVERLAPPING_PLACEMENTS_MAY_STILL_SHARE_LARGE_SCALE_NOISE",
        ],
    }
    # Quasi-independent positions across a single coadd are not
    # automatically representative. Capture that failure explicitly
    # even if the nominal sample covariance matrix is finite.
    if len(centres) > 0:
        if report["n_quadrants_with_accepted_blanks"] < 3:
            report["warnings"].append(
                "BLANK_APERTURE_SPATIAL_QUADRANT_COVERAGE_LIMITED")
        if max(accepted_by_quadrant) / len(centres) > .5:
            report["warnings"].append(
                "BLANK_APERTURE_SAMPLES_DOMINATED_BY_ONE_QUADRANT")
    if required > max_placements:
        report["status"] = "EXPLICIT_PLACEMENT_CAP_BELOW_REQUIRED_COVARIANCE_SAMPLE"
        report["warnings"].append(
            "REQUESTED_CAP_CANNOT_MEET_MINIMUM_SAMPLES_FOR_RADIAL_DIMENSION")
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


def _smaller_blank_scales(frame, full_outer, full_radii, source, manual, *,
                          full_scale_available):
    """Fixed prefix tests, never represent them as full-size covariance.

    A 105-arcsec blank ellipse cannot be tiled in a 4-arcmin cutout. Smaller
    apertures can still constrain native-coadd fluctuations at shorter scales
    IF they remain outside the *full* target footprint and are clean.
    """
    result = {
        "status": "NOT_REQUIRED_FULL_SCALE_DIAGNOSTIC_AVAILABLE"
        if full_scale_available else "SMALLER_SCALES_NOT_YET_ESTABLISHED",
        "full_target_semimajor_arcsec": full_outer.semimajor_arcsec,
        "full_scale_covariance_available": bool(full_scale_available),
        "probes": [],
        "first_available_probe_outer_arcsec": None,
        "must_not_replace_full_scale_covariance": True,
        "no_target_mask_or_fits_modified": True,
        "covariance_calibrated_for_science": False,
        "adopted_for_science": False,
    }
    if full_scale_available or len(full_radii) < 3:
        return result
    # Predeclared, decreasing radii: no data-dependent sampling selection.
    indices = []
    for value in (
            len(full_radii) - 4, len(full_radii) // 2 - 1, 2, 1):
        if 1 <= value < len(full_radii) - 1 and value not in indices:
            indices.append(value)
    indices.sort(reverse=True)
    for index in indices:
        radius = float(full_radii[index])
        small_outer = ApertureGeometryICRS(
            ra_deg=full_outer.ra_deg, dec_deg=full_outer.dec_deg,
            semimajor_arcsec=radius, axis_ratio=full_outer.axis_ratio,
            pa_deg_east_of_north=full_outer.pa_deg_east_of_north)
        probe = estimate_blank_annulus_covariance(
            frame, small_outer, full_radii[:index+1],
            effective_source_mask=source, manual_exclude=manual,
            full_target_exclusion_outer=full_outer)
        result["probes"].append(probe)
        if probe["annulus_covariance_nMgy2"] is not None:
            result["status"] = "SMALL_SCALE_COVARIANCE_DIAGNOSTIC_ONLY"
            result["first_available_probe_outer_arcsec"] = radius
            break
    if result["first_available_probe_outer_arcsec"] is None:
        result["status"] = "NO_CLEAN_SMALL_SCALE_COVARIANCE_AVAILABLE"
    return result



def _verify_frozen_growth_curves(frame, anchor_icrs, provenance,
                                 source, manual, curves, sky):
    """Fail closed on stale Stage 04.2 inputs, including same-area mask edits.

    Comparing only valid pixel area cannot identify changed source pixels,
    image values, inverse variances or a different conditional sky model.
    Verify each frozen receipt and recompute its *outermost* direct aperture
    against the original native frame and the exact current mask/sky recipe.
    This is bounded by the number of declared geometry hypotheses; no
    photometry, FITS or user draft is written.
    """
    if not isinstance(curves, (tuple, list)) or not curves:
        raise ValueError("At least one frozen Stage 04.2 curve required")
    expected_grid = grid_fingerprint(frame.wcs, frame.shape)
    reference_radii = tuple(curves[0]["radii_arcsec"])
    anchor = tuple(map(float, anchor_icrs))
    for curve in curves:
        if not isinstance(curve, dict):
            raise ValueError("Stage 04.2 curve must be a frozen QA mapping")
        receipt = sha256(json.dumps(
            {key: value for key, value in curve.items()
             if key != "receipt_sha256"},
            sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode("utf-8")).hexdigest()
        if curve.get("receipt_sha256") != receipt:
            raise ValueError("Stage 04.2 growth receipt integrity mismatch")
        expected = {
            "project_id": provenance.project_id,
            "dataset_id": frame.dataset_id,
            "grid_fingerprint": expected_grid,
            "draft_sha256": provenance.scene_draft_sha256,
            "catalogue_sha256": provenance.catalogue_sha256,
            "source_policy_id": provenance.source_policy_id,
            "sky_method": sky.method,
            "sky_evidence_id": sky.evidence_id,
        }
        if any(curve.get(k) != v for k, v in expected.items()):
            raise ValueError("Stage 04.2 lineage or sky scenario mismatch")
        if (tuple(curve.get("fixed_project_anchor_icrs", ())) != anchor or
                tuple(curve.get("radii_arcsec", ())) != reference_radii):
            raise ValueError("Stage 04.2 anchor or radial grid mismatch")
        final = curve["cumulative"][-1]
        hypothesis = curve["geometry_hypothesis"]
        outer = ApertureGeometryICRS(
            ra_deg=anchor[0], dec_deg=anchor[1],
            semimajor_arcsec=float(reference_radii[-1]),
            axis_ratio=float(hypothesis["q"]),
            pa_deg_east_of_north=float(hypothesis["pa_deg"]))
        fresh = measure_native_aperture(
            frame, outer, provenance,
            effective_source_mask=source, manual_exclude=manual, sky=sky)
        if fresh.recipe_sha256 != final.get("recipe_sha256"):
            raise ValueError(
                "Stage 04.2 aperture recipe mismatch: mask pixels, "
                "geometry, calibration lineage or sky changed")
        for name, recorded, current in (
            ("valid area", final["area_native_pixel2"]["valid_native_pixel2"],
             fresh.area_native_pixel2["valid_native_pixel2"]),
            ("signed flux", final["signed_valid_flux_nMgy"],
             fresh.signed_valid_flux_nMgy),
            ("formal variance", final["formal_diagonal_sigma_nMgy"],
             fresh.formal_diagonal_sigma_nMgy),
        ):
            if recorded is None and current is None:
                continue
            if (recorded is None or current is None or
                    not math.isclose(float(recorded), float(current),
                                     rel_tol=1e-11, abs_tol=1e-8)):
                raise ValueError(
                    "Stage 04.2 original native pixels or weights changed: "
                    + name)

def diagnose_aperture_systematics(
        frame: NativeBandFrame, *, anchor_icrs, provenance: ApertureProvenance,
        effective_source_mask, manual_exclude, curves,
        sky: ApertureSkyScenario | None = None):
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
    sky = ApertureSkyScenario() if sky is None else sky
    if not isinstance(sky, ApertureSkyScenario):
        raise TypeError("Versioned frozen sky scenario required")
    _verify_frozen_growth_curves(
        frame, anchor_icrs, provenance, source, manual, curves, sky)
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
        manual_exclude=manual, sky=sky)
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
    small_scale = _smaller_blank_scales(
        frame, outer, first["radii_arcsec"], source, manual,
        full_scale_available=(blanks["annulus_covariance_nMgy2"] is not None))
    report = {
        "schema": _SCHEMA,
        "project_id": provenance.project_id,
        "dataset_id": frame.dataset_id,
        "draft_sha256": provenance.scene_draft_sha256,
        "catalogue_sha256": provenance.catalogue_sha256,
        "source_policy_id": provenance.source_policy_id,
        "stage_04_2_receipts_and_outer_apertures_verified": len(curves),
        "conditional_sky_scenario": asdict(sky),
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
        "smaller_scale_blank_diagnostics": small_scale,
        "adopted_for_science": False,
        "covariance_calibrated_for_science": False,
        "status": "CONDITIONAL_SYSTEMATICS_DIAGNOSTIC_ONLY",
        "warnings": [
            "SKY_RESPONSE_IS_A_DERIVATIVE_NOT_A_SKY_ERROR_BAR",
            "RELAXED_CATALOGUE_MASK_INCLUDES_CONTAMINANTS_AND_TARGET_MEMBERS",
            "GEOMETRY_SPAN_COMPARES_DIFFERENT_SKY_FOOTPRINTS",
            "NO_TOTAL_FLUX_OR_PSF_CORRECTION_INFERRED",
        ] + ([
            "BLANK_COVARIANCE_USES_UNCORRECTED_NATIVE_COADD"
        ] if sky.method != "zero_additional_correction" else []),
    }
    json.dumps(report, allow_nan=False)
    return report
