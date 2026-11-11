"""Read-only native coadd joint scene decomposition around fixed PSF-0 stars.

This pilot *holds the already built empirical PSF kernel fixed* and solves
simultaneously for the central-star amplitude, compatible catalogued point
source amplitudes and a planar sky on a larger region of ORIGINAL coadd
pixels. It does not recover PSF flux beyond the empirical kernel support
or account for extended catalogued objects or correlated coadd covariance.

Neither a better residual nor a nonnegative fitted neighbor demonstrates
scientific PSF validation. No fit output feeds back into training, holdouts,
star rejection, sky subtraction, the target, or the adopted PSF state.
"""
from __future__ import annotations

from hashlib import sha256
import math

import numpy as np
from scipy.ndimage import map_coordinates
from scipy.optimize import nnls
from scipy.spatial import cKDTree


SCHEMA = "desidera.psf0-fixed-psf-joint-native-scene/3"
MAX_NUISANCE_SOURCES = 18  # engineering compute budget, NOT science calibration


def _finite_float(value):
    try:
        x = float(value)
    except (ValueError, TypeError, OverflowError):
        return None
    return x if math.isfinite(x) else None


def _shifted_native_kernel(kernel, xx, yy, x, y):
    """Sample frozen native PSF at exact subpixel Tractor/native coordinates."""
    kh, kw = kernel.shape
    coords = np.array([yy-y+(kh-1)/2., xx-x+(kw-1)/2.])
    # Cubic interpolation can introduce small negative numerical lobes.
    # The delivered PSF0 kernel is nonnegative; this diagnostic explicitly
    # uses a nonnegative basis without modifying the original kernel.
    return np.maximum(map_coordinates(
        kernel, coords, order=3, mode="constant", cval=0.,
        prefilter=True), 0.)


def _projected_nnls(y, B, P, w):
    """Solve an unpenalized free sky and nonnegative point-source fluxes.

    Project away the *unconstrained* three sky coefficients by QR, then
    solve NNLS for ALL point sources. This avoids arbitrary Tikhonov
    regularization or fitting sky to a pre-subtracted stellar image.
    """
    y = np.asarray(y, dtype=float)
    w = np.asarray(w, dtype=float)
    Bw = B*w[:, None]
    Pw = P*w[:, None]
    yw = y*w
    if len(y) <= P.shape[1]+3:
        raise ValueError("INSUFFICIENT_SCENE_PIXELS_FOR_JOINT_FIT")
    q, r = np.linalg.qr(Bw, mode="reduced")
    if np.linalg.matrix_rank(r) < 3:
        raise ValueError("DEGENERATE_SCENE_SKY_GEOMETRY")
    projected = Pw-q@(q.T@Pw)
    data_projected = yw-q@(q.T@yw)
    singular = np.linalg.svd(projected, compute_uv=False)
    if (not np.isfinite(singular).all()
            or not len(singular) or singular[-1] <=
            np.finfo(float).eps*max(projected.shape)*singular[0]):
        raise ValueError("NON_IDENTIFIABLE_SCENE_SOURCE_TEMPLATES")
    amps, _ = nnls(projected, data_projected)
    sky_coeff = np.linalg.lstsq(Bw, yw-Pw@amps, rcond=None)[0]
    if not np.isfinite(amps).all() or not np.isfinite(sky_coeff).all():
        raise ValueError("NONFINITE_JOINT_SCENE_SOLUTION")
    resid = yw-Bw@sky_coeff-Pw@amps
    return amps, sky_coeff, float(resid@resid), float(singular[0]/singular[-1])


def fit_fixed_psf_scene(
    original_image, original_invvar, original_valid, *,
    kernel, target_xy, neighbors, guessed_fwhm_pix, original_sky_nMgy_pix,
    scene_size_pix=None,
):
    """Fit star+plane+catalog PSF neighbors on actual, unmodified pixels.

    The patch width defaults to exactly twice the empirical PSF's linear
    pixel support minus one (81 -> 161), a declared *engineering probe*
    geometry, not a measured sky-free boundary. It is bounded by the
    input native patch; callers must not secretly pad or inpaint edges.
    """
    data = np.asarray(original_image, dtype=float)
    ivar = np.asarray(original_invvar, dtype=float)
    valid = np.asarray(original_valid, dtype=bool)
    k = np.asarray(kernel, dtype=float)
    if (data.ndim != 2 or data.shape != ivar.shape or data.shape != valid.shape
            or k.ndim != 2 or k.shape[0] != k.shape[1]
            or k.shape[0] % 2 != 1 or k.shape[0] < 15
            or not np.isfinite(k).all() or np.any(k < 0)
            or not math.isclose(float(np.sum(k)), 1., rel_tol=1e-8, abs_tol=1e-8)):
        raise ValueError("Invalid original native arrays or frozen PSF")
    width = 2*k.shape[0]-1 if scene_size_pix is None else scene_size_pix
    if (type(width) is not int or width < k.shape[0]
            or width > 2*k.shape[0]-1 or width % 2 == 0):
        raise ValueError("Invalid bounded native scene footprint")
    if (len(target_xy) != 2 or
            any(_finite_float(v) is None for v in target_xy) or
            _finite_float(guessed_fwhm_pix) is None or guessed_fwhm_pix <= 0 or
            _finite_float(original_sky_nMgy_pix) is None):
        raise ValueError("Invalid scene coordinates or native sky reference")
    if len(neighbors) > MAX_NUISANCE_SOURCES:
        return {
            "status": "TOO_MANY_PSF_NEIGHBORS_FOR_BOUNDED_FIT",
            "n_catalog_point_neighbors": len(neighbors),
            "science_adopted": False,
        }
    cx, cy = map(float, target_xy)
    half = width//2
    x0, y0 = int(np.rint(cx))-half, int(np.rint(cy))-half
    if x0 < 0 or y0 < 0 or x0+width > data.shape[1] or y0+width > data.shape[0]:
        return {"status": "FULL_EXPANDED_NATIVE_STAMP_OUTSIDE_IMAGE",
                "scene_width_pix": width, "science_adopted": False}
    rows = slice(y0, y0+width)
    cols = slice(x0, x0+width)
    z = data[rows, cols]
    inv = ivar[rows, cols]
    m = valid[rows, cols] & np.isfinite(z) & np.isfinite(inv) & (inv > 0)
    if not m[half, half]:
        return {"status": "INVALID_CENTRAL_NATIVE_PIXEL",
                "science_adopted": False}
    n_valid = int(np.count_nonzero(m))
    # Only impose the mathematical number-of-unknowns requirement;
    # no uncalibrated 25%-valid-stamp acceptance threshold is invented.
    # The QR/SVD tests below check the actual numerical rank.
    if n_valid <= len(neighbors) + 3:
        return {"status": "INSUFFICIENT_NATIVE_SCENE_PIXEL_DOF",
                "n_valid": n_valid, "n_total": width*width,
                "n_point_source_coefficients": len(neighbors)+1,
                "n_planar_sky_coefficients": 3,
                "science_adopted": False}
    gy, gx = np.indices(z.shape, dtype=float)
    gx += x0
    gy += y0
    extent = max(1., float(half))
    B = np.column_stack((
        np.ones(n_valid),
        ((gx-cx)/extent)[m],
        ((gy-cy)/extent)[m],
    ))
    center_template = _shifted_native_kernel(k, gx, gy, cx, cy)
    templates = [center_template]
    active_neighbors = []
    zero_support_neighbors = []
    # A catalog position can be inside the *approximate* patch+kernel
    # bounding rectangle but its interpolated finite-support template can
    # have EXACTLY zero usable native pixels (especially near corners and
    # masks). This is an identically zero design-matrix column, NOT a
    # mathematical degeneracy of the central star. Exclude it explicitly
    # with evidence; never silently discard a source with any model support.
    for row in neighbors:
        x, y = float(row["x_pix"]), float(row["y_pix"])
        if not math.isfinite(x) or not math.isfinite(y):
            raise ValueError("Nonfinite neighbor catalogue coordinate")
        template = _shifted_native_kernel(k, gx, gy, x, y)
        if not np.any(template[m] > 0.):
            zero_support_neighbors.append({
                "tractor_objid": int(row["tractor_objid"]),
                "separation_from_star_pix": float(math.hypot(x-cx, y-cy)),
                "status": "EXACT_ZERO_FROZEN_KERNEL_ON_ALL_VALID_SCENE_PIXELS",
            })
            continue
        active_neighbors.append(row)
        templates.append(template)
    columns = []
    template_peak_scales = []
    for temp in templates:
        t = temp[m]
        maximum = float(np.max(t)) if len(t) else 0.
        if maximum <= 0 or not math.isfinite(maximum):
            return {"status": "CENTRAL_OR_ACTIVE_TEMPLATE_HAS_NO_VALID_SUPPORT",
                    "science_adopted": False}
        columns.append(t/maximum)
        template_peak_scales.append(maximum)
    P = np.column_stack(columns)
    weight = np.sqrt(inv[m])
    try:
        amplitudes, sky, full_rss, condition = _projected_nnls(
            z[m], B, P, weight)
        null_amplitudes, null_sky, null_rss, null_condition = _projected_nnls(
            z[m], B, P[:, :1], weight)
    except (ValueError, RuntimeError) as exc:
        return {
            "status": "JOINT_SCENE_NUMERICALLY_UNIDENTIFIABLE",
            "reason": str(exc), "n_valid": n_valid,
            "n_catalog_point_neighbors": len(neighbors),
            "science_adopted": False,
        }
    if amplitudes[0] <= 0 or not math.isfinite(float(amplitudes[0])):
        return {
            "status": "NONPOSITIVE_FITTED_CENTRAL_STAR_AMPLITUDE",
            "n_valid": n_valid, "n_catalog_point_neighbors": len(neighbors),
            "science_adopted": False,
        }
    # Full and null fits use the *identical original native pixels*.
    # A nested fit is guaranteed to be at least as good. A reduced RSS
    # is NOT a significance, p-value, evidence ratio or true-PSF metric.
    r2 = np.hypot(gx-cx, gy-cy)
    within_core = m & (r2 <= guessed_fwhm_pix)
    within_outer = m & (r2 <= 4*guessed_fwhm_pix)
    # Global residuals can be dominated by other bright stars near the
    # edge of this 161px scene. Separately compare the IDENTICAL central
    # apertures, still using amplitudes fitted globally in both models.
    # Local residual improvements may be NEGATIVE: never clamp them.
    full_residual = (z[m] - B@sky - P@amplitudes) * weight
    null_residual = (z[m] - B@null_sky - P[:, :1]@null_amplitudes)*weight
    native_radius = r2[m]
    # The full-minus-null image model is an exact sum of THREE changes:
    # change in the central star amplitude, change in the free sky plane,
    # and added catalog point-neighbor templates. This is a descriptive
    # decomposition of fitted model changes, NOT astrophysical causation.
    model_changes = {
        "central_star_amplitude_reallocation": (
            weight * P[:, 0] * (amplitudes[0]-null_amplitudes[0])),
        "planar_sky_reallocation": (
            weight * (B @ (sky-null_sky))),
        "direct_catalog_point_neighbors": (
            weight * (P[:, 1:] @ amplitudes[1:])),
    }
    total_model_change = sum(model_changes.values())
    def radial_residual_comparison(radius_pix):
        aperture = (np.ones(len(native_radius), dtype=bool)
                    if radius_pix is None else native_radius <= radius_pix)
        if not np.any(aperture):
            return {
                "status": "NO_VALID_PIXELS_IN_APERTURE",
                "n_native_pixels": 0,
                "radius_pix": float(radius_pix),
                "weighted_squared_residual_full": None,
                "weighted_squared_residual_star_plane_only": None,
                "signed_improvement_fraction": None,
            }
        sum_full = float(np.sum(full_residual[aperture]**2))
        sum_null = float(np.sum(null_residual[aperture]**2))
        # Since RSS(null)-RSS(full) = 2<e_null,d> - <d,d>,
        # allocate each cross-term symmetrically to its two components:
        # G_i = 2<e_null,d_i> - <d_i,d>. Their sum reproduces EXACTLY
        # the signed RSS gain, including a local worsening. This is
        # an algebraic partition, not a causal source-identification
        # or a significance test. The fitted sky/stellar amplitudes are
        # mutually coupled and depend on the entire expanded field.
        change = total_model_change[aperture]
        residual = null_residual[aperture]
        gains = {
            name: float(2*np.dot(residual, values[aperture])
                        -np.dot(values[aperture], change))
            for name, values in model_changes.items()
        }
        neighbors = model_changes["direct_catalog_point_neighbors"][aperture]
        return {
            "status": "SAME_NATIVE_PIXELS_WITH_GLOBAL_PARAMETERS",
            "n_native_pixels": int(np.count_nonzero(aperture)),
            "radius_pix": (float(radius_pix)
                           if radius_pix is not None else None),
            "weighted_squared_residual_full": sum_full,
            "weighted_squared_residual_star_plane_only": sum_null,
            "signed_improvement_fraction": (
                float((sum_null-sum_full)/sum_null) if sum_null > 0 else None),
            "not_a_covariance_calibrated_or_independent_validation": True,
            "retains_negative_improvement_when_central_fit_worsens": True,
            "weighted_rss_gain_partition_by_fitted_model_change": gains,
            "weighted_rss_gain_partition_closure_error": float(
                sum(gains.values())-(sum_null-sum_full)),
            "direct_catalog_point_neighbors_model_power_weighted":
                float(neighbors@neighbors),
            "gain_partition_is_algebraic_noncausal_with_symmetric_cross_terms": True,
            "background_and_central_amplitude_changes_are_globally_coupled": True,
            "no_new_local_refit_or_science_adoption": True,
        }
    def neighbor_fraction(aperture):
        if not np.any(aperture):
            return None
        central_flux = float(
            (amplitudes[0]/template_peak_scales[0])*
            np.sum(templates[0][aperture]))
        other_flux = float(sum(
            (amplitudes[i]/template_peak_scales[i])*
            np.sum(templates[i][aperture])
            for i in range(1, len(templates))))
        return (other_flux/central_flux if central_flux > 0 else None)

    result = {
        "schema": SCHEMA,
        "status": "FIXED_PSF_JOINT_SCENE_DIAGNOSTIC",
        "scene_width_pix": width,
        "scene_geometry": "2*EMPIRICAL_KERNEL_WIDTH_MINUS_ONE",
        "scene_size_is_not_a_measured_wing_boundary": True,
        "n_original_native_pixels": width*width,
        "n_valid_native_pixels": n_valid,
        "valid_native_pixel_fraction": float(n_valid/(width*width)),
        "no_arbitrary_valid_pixel_fraction_cut_imposed": True,
        "n_catalog_point_neighbors": len(active_neighbors),
        "n_catalog_point_neighbors_in_expanded_search": len(neighbors),
        "n_catalog_neighbors_with_exact_zero_finite_psf_support":
            len(zero_support_neighbors),
        "catalog_neighbors_with_exact_zero_finite_psf_support":
            zero_support_neighbors,
        "zero_support_exclusion_is_only_relative_to_frozen_kernel": True,
        "neighbors": [{
            "tractor_objid": int(row["tractor_objid"]),
            "tractor_type": row["tractor_type"],
            "separation_from_star_pix": float(
                math.hypot(row["x_pix"]-cx, row["y_pix"]-cy)),
            "catalog_band_flux_nMgy": float(row["flux_nMgy"]),
            "fitted_peak_nMgy_per_pixel": float(amplitudes[i+1]),
            "fitted_flux_within_model_support_nMgy": float(
                amplitudes[i+1]/template_peak_scales[i+1] *
                np.sum(templates[i+1])),
        } for i, row in enumerate(active_neighbors)],
        "target_fitted_peak_nMgy_per_pixel": float(amplitudes[0]),
        "star_only_null_fitted_peak_nMgy_per_pixel":
            float(null_amplitudes[0]),
        "full_minus_null_fitted_peak_nMgy_per_pixel":
            float(amplitudes[0]-null_amplitudes[0]),
        "star_only_null_fitted_sky_at_star_nMgy_per_native_pixel":
            float(null_sky[0]),
        "full_minus_null_fitted_sky_at_star_nMgy_per_native_pixel":
            float(sky[0]-null_sky[0]),
        "target_fitted_model_support_flux_nMgy": float(
            amplitudes[0]/template_peak_scales[0] * np.sum(templates[0])),
        "fitted_sky_at_star_nMgy_per_native_pixel": float(sky[0]),
        "fitted_sky_gradient_x_nMgy_per_native_pixel_per_pix":
            float(sky[1]/extent),
        "fitted_sky_gradient_y_nMgy_per_native_pixel_per_pix":
            float(sky[2]/extent),
        "original_annulus_sky_at_star_nMgy_per_native_pixel":
            float(original_sky_nMgy_pix),
        "joint_minus_annulus_sky_nMgy_per_native_pixel":
            float(sky[0]-original_sky_nMgy_pix),
        "neighbor_to_central_core_model_flux_ratio":
            neighbor_fraction(within_core),
        "neighbor_to_central_4fwhm_model_flux_ratio":
            neighbor_fraction(within_outer),
        "weighted_squared_residual_star_plane_only": float(null_rss),
        "weighted_squared_residual_star_neighbors_plane": float(full_rss),
        "central_aperture_residual_comparison": {
            "entire_native_scene": radial_residual_comparison(None),
            "one_guessed_fwhm": radial_residual_comparison(
                float(guessed_fwhm_pix)),
            "four_guessed_fwhm": radial_residual_comparison(
                float(4*guessed_fwhm_pix)),
            "same_global_fit_coefficients_no_local_refit": True,
            "not_an_independent_psf_accuracy_test": True,
        },
        "weighted_residual_improvement_fraction": (
            float(max(0., null_rss-full_rss)/null_rss)
            if null_rss > 0 else None),
        "joint_source_matrix_condition_number": condition,
        "null_source_matrix_condition_number": null_condition,
        "joint_fit_uses_original_unsmoothed_unsubtracted_pixels": True,
        "star_and_neighbor_templates_fixed_to_existing_empirical_kernel": True,
        "point_source_fits_nonnegative": True,
        "sky_fit_unconstrained_planar": True,
        "kernel_support_beyond_observed_ePSF_unmodeled": True,
        "catalog_extended_neighbors_not_forward_modeled": True,
        "correlated_noise_and_fit_parameter_covariance_unmodeled": True,
        "comparison_fit_has_more_parameters_and_cannot_be_science_gate": True,
        "science_adopted": False,
        "changes_psf_or_source_selection_or_original_pixels": False,
    }
    return result


def audit_fixed_psf_joint_scenes(
    frame, kernel, stars, catalog_matcher, *,
    protected_mask, guessed_fwhm_pix,
):
    """Evaluate the SAME preselected stars. Never resplit or retune sample."""
    if catalog_matcher is None:
        return {
            "schema": SCHEMA, "status": "NO_NATIVE_TRACTOR_CATALOG",
            "n_star_rows": len(stars), "science_adopted": False,
        }
    if catalog_matcher.dataset_id != frame.dataset_id:
        raise ValueError("Joint scene catalog dataset does not match native FITS")
    kernel_arr = np.asarray(kernel, dtype=float)
    half = kernel_arr.shape[0]-1  # half-width of (2*N - 1) patch
    kh = kernel_arr.shape[0]//2
    max_extent = half+kh+1
    index = cKDTree(catalog_matcher.pixel_positions)
    if protected_mask.shape != frame.shape:
        raise ValueError("Protected target mask differs from native coadd")
    # Compute the original native acceptance mask ONCE per analysis band;
    # avoid copying a full 3600x3600 frame separately for every star.
    valid = np.asarray(frame.valid, dtype=bool)
    for key in (f"SATUR_{frame.band.upper()}",
                f"ALLMASK_{frame.band.upper()}"):
        if key in frame.mask_reasons:
            valid = valid & ~frame.mask_reasons[key]
    valid = valid & ~protected_mask
    outputs = []
    for star in stars:
        x, y = float(star.x_pix), float(star.y_pix)
        item = {"source_id": str(star.star_id), "role": str(star.role),
                "gaia_source_id": star.gaia_source_id}
        matches = index.query_ball_point((x, y),
                                         r=float(catalog_matcher.match_radius_pix))
        if len(matches) != 1:
            item["status"] = "CENTRAL_TRACTOR_MATCH_NOT_UNIQUE"
            item["n_central_matches"] = len(matches)
            outputs.append(item)
            continue
        own = int(matches[0])
        own_row = catalog_matcher.rows[own]
        if (star.gaia_source_id is not None and
                int(own_row.get("ref_id", -1)) != star.gaia_source_id):
            item["status"] = "CENTRAL_GAIA_ID_MISMATCH"
            outputs.append(item)
            continue
        if own_row.get("type") != "PSF":
            item["status"] = "CENTRAL_TRACTOR_SOURCE_NOT_PSF"
            outputs.append(item)
            continue
        neighboring = []
        extended = []
        # Query source positions over the patch+full ePSF support, NOT
        # only within the existing 50px significant-neighbor preflight.
        candidates = index.query_ball_point((x,y), r=math.sqrt(2)*max_extent)
        for j in candidates:
            if j == own:
                continue
            px, py = map(float, catalog_matcher.pixel_positions[j])
            if abs(px-x) > max_extent or abs(py-y) > max_extent:
                continue
            row = catalog_matcher.rows[j]
            flux = _finite_float(row.get("flux"))
            record = {
                "x_pix":px, "y_pix":py,
                "tractor_objid": int(row["objid"]),
                "tractor_type": row["type"],
                "flux_nMgy": flux,
                # These MUST originate from the original DR11 Tractor
                # table; missing morphology remains None, never guessed.
                "shape_r_arcsec": row.get("shape_r_arcsec"),
                "shape_e1": row.get("shape_e1"),
                "shape_e2": row.get("shape_e2"),
                "sersic_n": row.get("sersic_n"),
                "shape_r_ivar_arcsec_m2":
                    row.get("shape_r_ivar_arcsec_m2"),
                "shape_e1_ivar": row.get("shape_e1_ivar"),
                "shape_e2_ivar": row.get("shape_e2_ivar"),
            }
            if row["type"] == "PSF" and flux is not None and flux > 0:
                neighboring.append(record)
            else:
                extended.append(record)
        neighboring.sort(key=lambda r: (
            math.hypot(r["x_pix"]-x,r["y_pix"]-y),r["tractor_objid"]))
        extended.sort(key=lambda r: (
            math.hypot(r["x_pix"]-x,r["y_pix"]-y),r["tractor_objid"]))
        item["n_unmodeled_extended_or_nonpositive_catalog_sources"] = len(extended)
        item["nearest_unmodeled_source_distance_pix"] = (
            float(math.hypot(extended[0]["x_pix"]-x,
                             extended[0]["y_pix"]-y))
            if extended else None)
        item["nearest_unmodeled_catalog_source"] = ({
            "tractor_objid": extended[0]["tractor_objid"],
            "tractor_type": extended[0]["tractor_type"],
            "catalog_band_flux_nMgy": extended[0]["flux_nMgy"],
            "separation_from_star_pix":
                item["nearest_unmodeled_source_distance_pix"],
            "reason": (
                "TRACTOR_NON_PSF_MODEL" if
                extended[0]["tractor_type"] != "PSF"
                else "PSF_NONPOSITIVE_OR_INVALID_CATALOG_FLUX"),
            "source_not_forward_modeled": True,
        } if extended else None)
        item["catalog_point_neighbor_count_full_support"] = len(neighboring)
        item["catalog_scene_fully_modelled"] = len(extended) == 0
        item["known_unmodeled_catalog_sources_prevent_science_certification"] = (
            len(extended) > 0)
        item["nearest_catalog_point_neighbor_distance_pix"] = (
            float(math.hypot(neighboring[0]["x_pix"]-x,
                             neighboring[0]["y_pix"]-y))
            if neighboring else None)
        fitted = fit_fixed_psf_scene(
            frame.image_nMgy, frame.invvar_nMgy2, valid,
            kernel=kernel_arr, target_xy=(x,y),
            neighbors=neighboring,
            guessed_fwhm_pix=guessed_fwhm_pix,
            original_sky_nMgy_pix=star.cutout_sky_nMgy)
        item.update(fitted)
        # P15: one separate, non-adopting observational comparison of
        # actual DR11 REX/EXP/DEV morphology, never fed back into PSF.
        # A local import avoids a cycle: the extended renderer reuses
        # the already tested projected-NNLS primitives in this module.
        if fitted["status"] != "FIXED_PSF_JOINT_SCENE_DIAGNOSTIC":
            item["p15_extended_tractor_scene"] = {
                "status":"NOT_RUN_NO_IDENTIFIABLE_P14_POINT_SCENE",
                "science_adopted":False,
            }
        elif getattr(frame,"wcs",None) is None:
            item["p15_extended_tractor_scene"] = {
                "status":"NOT_RUN_NO_VERIFIED_ORIGINAL_WCS",
                "science_adopted":False,
            }
        else:
            from lsbpipe.photometry.psf_extended_scene import (
                pilot_extended_scene,
            )
            item["p15_extended_tractor_scene"] = pilot_extended_scene(
                frame.image_nMgy, frame.invvar_nMgy2, valid,
                wcs=frame.wcs, kernel=kernel_arr,
                target_xy=(x,y), point_sources=neighboring,
                extended_sources=extended,
                guessed_fwhm_pix=guessed_fwhm_pix)
        outputs.append(item)
    ok = [row for row in outputs
          if row["status"] == "FIXED_PSF_JOINT_SCENE_DIAGNOSTIC"]
    return {
        "schema": SCHEMA,
        "status": "FIXED_PSF_ORIGINAL_NATIVE_SCENE_AUDIT_NO_SCIENCE_ADOPTION",
        "dataset_id": frame.dataset_id,
        "native_tractor_catalog_sha256": catalog_matcher.source_sha256,
        "frozen_input_kernel_sha256":
            sha256(kernel_arr.tobytes(order="C")).hexdigest(),
        "n_preselected_stars": len(stars),
        "n_completed_joint_scene_fits": len(ok),
        "n_failed_joint_scene_fits": len(stars)-len(ok),
        "n_completed_scenes_with_unmodeled_catalog_sources": sum(
            row.get("known_unmodeled_catalog_sources_prevent_science_certification",
                    False)
            for row in ok),
        "n_completed_scenes_with_zero_support_catalog_psf_neighbors": sum(
            row.get("n_catalog_neighbors_with_exact_zero_finite_psf_support",
                    0) > 0 for row in ok),
        "n_total_zero_support_catalog_psf_neighbors_in_completed_fits": sum(
            row.get("n_catalog_neighbors_with_exact_zero_finite_psf_support",
                    0) for row in ok),
        "n_completed_fixed_holdouts": sum(
            row["role"] == "holdout" for row in ok),
        "n_holdouts_expected": sum(star.role == "holdout" for star in stars),
        "per_star": outputs,
        "p15_extended_scene_pilot_summary": {
            "schema": "desidera.psf0-native-tractor-extended-scene-summary/1",
            "n_original_stars":len(stars),
            "n_attempted_after_p14":sum(
                row.get("p15_extended_tractor_scene", {}).get("status")
                not in ("NOT_RUN_NO_IDENTIFIABLE_P14_POINT_SCENE",
                        "NOT_RUN_NO_VERIFIED_ORIGINAL_WCS")
                for row in outputs),
            "n_numerical_partial_extended_fits":sum(
                row.get("p15_extended_tractor_scene", {}).get("status")
                == "PARTIAL_APPROXIMATE_EXTENDED_SCENE_DIAGNOSTIC"
                for row in outputs),
            "n_fits_with_unmodelled_extended_sources":sum(
                row.get("p15_extended_tractor_scene", {}).get(
                    "n_extended_catalog_unmodeled",0)>0 for row in outputs),
            "no_psf_adoption_or_star_split_changes":True,
            "not_exact_tractor_profile_or_noise_covariance":True,
            "science_adopted":False,
        },
        "all_stars_audited_without_resplitting": True,
        "catalog_only_nuisance_sources_are_not_new_training_stars": True,
        "frozen_psf_still_has_no_wide_angle_wings": True,
        "unmodeled_extended_sources_not_certified_clean": True,
        "coadd_covariance_not_calibrated": True,
        "nuisance_parameter_fit_to_holdouts_not_blind_validation": True,
        "not_used_to_select_psf_recipe_or_science_thresholds": True,
        "science_adopted": False,
    }



def audit_finite_psf0_foreground(
        image, invvar, native_hard, *, kernel, source_xy, target_xy,
        measured_fwhm_pix, source_id, catalog_flux_nMgy=None,
        other_catalog_centres=(), component_indices=None):
    """Read-only local *finite* PSF-0 star probe, never wing subtraction.

    A source-anchored fit uses two predeclared sky hypotheses. Training
    consists of 0.8..2.5 FWHM stellar pixels and 4..4.5 FWHM sky-control
    pixels; 2.5..4 FWHM is reserved as a held-out check. The locally
    normalized coadd ePSF is immutable. Neither hypothesis estimates
    flux outside that kernel, true stellar wings or target missing light.
    """
    schema = "desidera.scene.psf0-foreground-finite-probe/1"
    base = {
        "schema": schema, "source_id": str(source_id),
        "status": "UNASSESSED",
        "science_adopted": False, "psf_wings_validated": False,
        "background_validated": False, "correlated_noise_validated": False,
        "source_light_subtracted": False, "original_pixels_modified": False,
        "physical_target_flux_loss_measured": False,
        "model_support_beyond_stored_psf": "UNKNOWN_NOT_EXTRAPOLATED",
        "training_fwhm_intervals": [[0.8, 2.5], [4., 4.5]],
        "holdout_fwhm_interval": [2.5, 4.],
        "sky_hypotheses": ["constant", "planar"],
        "warning": ("An experimental empirical coadd core is not a calibrated "
                    "bright-star wide-angle PSF. No science-mask, sky or "
                    "photometry update follows from this diagnostic."),
    }
    a, iv, hard = (np.asarray(image, dtype=float),
                   np.asarray(invvar, dtype=float),
                   np.asarray(native_hard, dtype=bool))
    k = np.asarray(kernel, dtype=float)
    if (a.ndim != 2 or iv.shape != a.shape or hard.shape != a.shape or
            k.ndim != 2 or k.shape[0] != k.shape[1] or
            k.shape[0] < 15 or k.shape[0] % 2 != 1 or
            not np.isfinite(k).all() or (k < 0).any() or
            not np.isclose(k.sum(), 1., rtol=1e-8, atol=1e-8)):
        raise ValueError("Native arrays or restored normalized PSF0 invalid")
    try:
        sx, sy = map(float, source_xy)
        tx, ty = map(float, target_xy)
        fwhm = float(measured_fwhm_pix)
    except (ValueError, TypeError) as exc:
        raise ValueError("Invalid source, target or PSF0 geometry") from exc
    if not (all(math.isfinite(z) for z in (sx, sy, tx, ty, fwhm))
            and 2. <= fwhm <= 25. and bool(source_id)):
        raise ValueError("Invalid measured PSF0 FWHM or coordinates")
    half_kernel = (k.shape[0]-1)//2
    offset = float(math.hypot(sx-tx, sy-ty))
    base.update(
        measured_fwhm_native_pix=fwhm,
        finite_kernel_width_pix=int(k.shape[0]),
        finite_kernel_half_width_native_pix=int(half_kernel),
        star_to_target_centre_native_pix=offset,
        source_xy_native_pix=[sx, sy],
        target_xy_native_pix=[tx, ty],
        catalog_flux_nMgy=(float(catalog_flux_nMgy)
                           if _finite_float(catalog_flux_nMgy) is not None
                           else None),
        galaxy_centre_outside_finite_kernel_square=bool(
            abs(sx-tx) > half_kernel or abs(sy-ty) > half_kernel),
        galaxy_centre_stellar_wing_flux_nMgy_pixel=None,
    )
    if half_kernel < math.ceil(4.5*fwhm)+2:
        return dict(base, status="FINITE_KERNEL_TOO_SMALL_FOR_HELDOUT_TEST")
    half = int(math.ceil(4.5*fwhm)+2)
    x0, y0 = int(np.rint(sx))-half, int(np.rint(sy))-half
    n = 2*half+1
    if x0 < 0 or y0 < 0 or x0+n > a.shape[1] or y0+n > a.shape[0]:
        return dict(base, status="FULL_NATIVE_STAR_PATCH_OUTSIDE_VIEW")
    base["probe_width_native_pix"] = n
    nearby = []
    for entry in other_catalog_centres:
        objid, ox, oy = entry
        if all(math.isfinite(float(z)) for z in (ox,oy)):
            separation = math.hypot(sx-float(ox), sy-float(oy))
            if separation < 6*fwhm:
                nearby.append({
                    "objid":int(objid),
                    "distance_native_pix":float(separation),
                })
    nearby.sort(key=lambda e:e["distance_native_pix"])
    base["other_tractor_centres_within_six_fwhm"] = nearby[:20]
    if nearby and nearby[0]["distance_native_pix"] < 2.8*fwhm:
        return dict(base, status="CLOSE_CATALOG_NEIGHBOR_PREVENTS_ISOLATED_FIT")
    yy, xx = np.indices((n,n), dtype=float)
    xx += x0
    yy += y0
    radial = np.hypot(xx-sx, yy-sy)/fwhm
    stamp = a[y0:y0+n,x0:x0+n]
    ivs = iv[y0:y0+n,x0:x0+n]
    usable = (~hard[y0:y0+n,x0:x0+n] &
              np.isfinite(stamp) & np.isfinite(ivs) & (ivs > 0))
    train = usable & (((radial >= .8) & (radial < 2.5)) |
                      ((radial >= 4.) & (radial <= 4.5)))
    held = usable & (radial >= 2.5) & (radial < 4.)
    base["native_training_pixels"] = int(train.sum())
    base["native_heldout_pixels"] = int(held.sum())
    base["native_training_valid_fraction"] = float(
        train.sum()/max(1, np.count_nonzero(
            ((radial >= .8) & (radial < 2.5)) |
            ((radial >= 4.) & (radial <= 4.5)))))
    base["native_heldout_valid_fraction"] = float(
        held.sum()/max(1,np.count_nonzero(
            (radial >= 2.5) & (radial < 4.))))
    if train.sum() < 80 or held.sum() < 80:
        return dict(base, status="INSUFFICIENT_NATIVE_TRAIN_OR_HELDOUT_PIXELS")
    template = _shifted_native_kernel(k, xx, yy, sx, sy)
    if not np.any(template[train] > 0):
        return dict(base, status="NO_EMPIRICAL_PSF_SUPPORT_IN_TRAINING_RING")
    # Two fits use IDENTICAL original pixels and a strictly held-out ring.
    results = {}
    for sky_name in ("constant", "planar"):
        if sky_name == "constant":
            sky_matrix = np.ones((n,n,1),dtype=float)
        else:
            sky_matrix = np.stack((
                np.ones((n,n),dtype=float),
                (xx-sx)/half, (yy-sy)/half),axis=-1)
        B = sky_matrix[train]
        w = np.sqrt(ivs[train])
        w /= np.median(w)
        Bw = B*w[:,None]
        yw = stamp[train]*w
        pw = template[train]*w
        q, rank = np.linalg.qr(Bw,mode="reduced")
        if np.linalg.matrix_rank(rank) != B.shape[1]:
            return dict(base,status="DEGENERATE_FOREGROUND_SKY_GEOMETRY")
        projected = pw-q@(q.T@pw)
        response = yw-q@(q.T@yw)
        denom = float(projected@projected)
        if not denom > np.finfo(float).eps*max(1.,float(pw@pw)):
            return dict(base,status="NON_IDENTIFIABLE_FOREGROUND_PSF_AMPLITUDE")
        amplitude = max(0.,float((projected@response)/denom))
        sky = np.linalg.lstsq(Bw,yw-amplitude*pw,rcond=None)[0]
        sky_null = np.linalg.lstsq(Bw,yw,rcond=None)[0]
        Bh = sky_matrix[held]
        observed = stamp[held]
        model = Bh@sky + amplitude*template[held]
        model_null = Bh@sky_null
        held_w = np.sqrt(ivs[held])
        residual = observed-model
        residual_null = observed-model_null
        angle = np.mod(np.arctan2((yy-sy)[held],
                                  (xx-sx)[held]),2*np.pi)
        sectors = np.minimum((angle/(np.pi/4)).astype(int),7)
        octants = [
            (float(np.median(residual[sectors == octant]))
             if np.any(sectors == octant) else None)
            for octant in range(8)]
        results[sky_name] = {
            "frozen_finite_psf_amplitude_nMgy":float(amplitude),
            "sky_at_source_nMgy_pixel":float(sky[0]),
            "sky_gradient_x_nMgy_pixel_per_pix":(
                float(sky[1]/half) if len(sky)>1 else None),
            "sky_gradient_y_nMgy_pixel_per_pix":(
                float(sky[2]/half) if len(sky)>1 else None),
            "heldout_weighted_rms_empirical_not_significance":float(
                np.sqrt(np.mean((residual*held_w)**2))),
            "heldout_sky_only_weighted_rms_empirical":float(
                np.sqrt(np.mean((residual_null*held_w)**2))),
            "heldout_median_residual_nMgy_pixel_by_octant":octants,
            "same_heldout_pixels_as_other_sky_hypothesis":True,
        }
    if not all(math.isfinite(r["frozen_finite_psf_amplitude_nMgy"])
               for r in results.values()):
        return dict(base,status="NONFINITE_FOREGROUND_PSF_PROBE")
    amplitude = results["planar"]["frozen_finite_psf_amplitude_nMgy"]
    if component_indices is not None:
        inds = np.asarray(component_indices)
        if (inds.ndim != 1 or not np.issubdtype(inds.dtype,np.integer)
                or np.any(inds < 0) or np.any(inds >= a.size)):
            raise ValueError("Candidate component indices not on native image")
        cx = inds%a.shape[1]
        cy = inds//a.shape[1]
        covered = ((np.abs(cx-sx) <= half_kernel) &
                   (np.abs(cy-sy) <= half_kernel))
        if np.any(covered):
            finite_template = _shifted_native_kernel(
                k,cx[covered].astype(float),cy[covered].astype(float),sx,sy)
            in_support = int(np.count_nonzero(finite_template > 0))
            finite_flux = float(amplitude*np.sum(finite_template))
        else:
            in_support = 0
            finite_flux = 0.
        base["candidate_component_finite_psf0_probe"] = {
            "component_pixels":int(inds.size),
            "pixels_with_positive_finite_kernel":in_support,
            "model_flux_nMgy_within_positive_finite_kernel":finite_flux,
            "missing_flux_and_true_physical_star_wings_unconstrained":True,
            "not_an_estimate_of_galaxy_flux_or_flux_loss":True,
        }
    base.update(
        status="FINITE_PSF0_LOCAL_PROBE_DIAGNOSTIC_NOT_ADOPTED",
        sky_fit_comparison=results,
        nearest_catalog_neighbour_distance_native_pix=(
            nearby[0]["distance_native_pix"] if nearby else None),
        near_neighbour_wing_risk_unquantified=True,
        target_light_in_star_fit_region_not_independently_excluded=True,
        model_subpixel_interpolation_not_exactly_flux_conserving=True,
    )
    return base



# P21: a SOURCE-BOUND DIRECTIONAL sensitivity test on original pixels only.
# It is neither a new PSF, a fitted halo, nor a science-ready sky estimate.
P21_SECTOR_RINGS_FWHM = ((6., 8.), (8., 12.), (12., 16.),
                         (16., 20.), (20., 24.))
P21_OUTER_REFERENCE_RING_FWHM = (26., 30.)
# Conservatively label annuli with sparse or strongly unequal mask support.
# These are QA warnings, never significance/wing certification thresholds.
P21_MIN_CONE_GEOMETRIC_COVERAGE = 0.70
P21_MAX_CONE_COVERAGE_DIFFERENCE = 0.20
P21_MIN_ANTIPODAL_PAIRS = 30


def audit_same_project_foreground_sector_profile(
        image, invvar, native_hard, *, source_xy, target_xy,
        measured_fwhm_pix, other_catalog_centres=()):
    """Read-only directional star/sky comparison strictly on the loaded image.

    Use fixed 90-degree cones facing / opposing the target, five radial
    shells and two *conditional* outer-sky rings. Compare the angular
    difference with 4x4 fixed-tile catalogue-only controls from the SAME
    native image. A control is not asserted to be galaxy-free. No star or
    sky flux, external brick, additional source, PSF wing or mask is fitted.

    Control locations are specified before reading any control pixels.
    Failed or catalogue-crowded tiles are not replaced using their values.
    """
    base = {
        "schema": "desidera.scene.project-only-foreground-sector/2",
        "status": "NOT_MEASURED",
        "science_adopted": False,
        "only_current_project_loaded_native_image": True,
        "additional_bricks_or_catalogues_used": False,
        "original_pixels_modified": False,
        "absolute_stellar_wings_identified": False,
        "galaxy_light_loss_measured": False,
        "background_galaxy_free_certified": False,
        "noise_covariance_calibrated": False,
        "fixed_radial_rings_measured_fwhm": [
            list(v) for v in P21_SECTOR_RINGS_FWHM],
        "conditional_outer_sky_rings_measured_fwhm":
            [[24., 26.], list(P21_OUTER_REFERENCE_RING_FWHM)],
        "control_selection": "FIXED_4x4_NATIVE_CROP_TILES_NO_BACKFILL",
        "direction": "TOWARD_VS_AWAY_FROM_PROJECT_TARGET_90_DEG_CONES",
        "matched_opposite_pixel_method":
            "NEAREST_NATIVE_PIXEL_REFLECTED_ABOUT_SUBPIXEL_CENTER",
        "min_cone_geometric_coverage_for_descriptive_comparison":
            P21_MIN_CONE_GEOMETRIC_COVERAGE,
        "max_cone_coverage_fraction_difference":
            P21_MAX_CONE_COVERAGE_DIFFERENCE,
        "min_antipodal_matched_pairs": P21_MIN_ANTIPODAL_PAIRS,
        "warning": ("Galaxy surface brightness and unmodelled stellar wings "
                    "can both contribute in the target-facing cone. "
                    "Sky-ring and blank controls are conditional, not "
                    "independent calibration or a halo upper bound."),
    }
    a = np.asarray(image, dtype=float)
    iv = np.asarray(invvar, dtype=float)
    hard = np.asarray(native_hard, dtype=bool)
    if a.ndim != 2 or iv.shape != a.shape or hard.shape != a.shape:
        raise ValueError("Project-only probe expects matching native image arrays")
    try:
        sx, sy = map(float, source_xy)
        tx, ty = map(float, target_xy)
        fwhm = float(measured_fwhm_pix)
    except (TypeError, ValueError) as exc:
        raise ValueError("Invalid native source/target geometry") from exc
    if (not all(math.isfinite(z) for z in (sx, sy, tx, ty, fwhm))
            or not 2. <= fwhm <= 25.):
        raise ValueError("Nonfinite source geometry or unmeasured FWHM")
    dist = float(math.hypot(tx-sx, ty-sy))
    if dist <= fwhm:
        return dict(base, status="SOURCE_TARGET_DIRECTION_NOT_IDENTIFIABLE")
    ux, uy = (tx-sx)/dist, (ty-sy)/dist
    neighbours = []
    for objid, x, y in other_catalog_centres:
        x, y = float(x), float(y)
        if math.isfinite(x) and math.isfinite(y):
            neighbours.append((int(objid), x, y))
    half = int(math.ceil(P21_OUTER_REFERENCE_RING_FWHM[1]*fwhm))+2
    rings = P21_SECTOR_RINGS_FWHM
    cos_limit = math.cos(math.pi/4.)
    rows, columns = a.shape

    def measure(cx, cy, catalog_centres):
        x0, y0 = int(round(cx))-half, int(round(cy))-half
        x1, y1 = x0+2*half+1, y0+2*half+1
        if x0 < 0 or y0 < 0 or x1 > columns or y1 > rows:
            return {"status": "OUTER_30_FWHM_PATCH_CROSSES_NATIVE_IMAGE_EDGE"}
        patch = a[y0:y1, x0:x1]
        inverse = iv[y0:y1, x0:x1]
        usable = (~hard[y0:y1, x0:x1] & np.isfinite(patch)
                  & np.isfinite(inverse) & (inverse > 0))
        yy, xx = np.indices(patch.shape, dtype=float)
        xx += x0
        yy += y0
        dx, dy = xx-cx, yy-cy
        radius = np.hypot(dx, dy)/fwhm
        # The 3-FWHM cores only remove identified catalogue neighbours.
        # Their unknown outer wings are NOT declared removed.
        catalogue_masked = 0
        for _, nx, ny in catalog_centres:
            if abs(nx-cx) > half+3*fwhm or abs(ny-cy) > half+3*fwhm:
                continue
            hit = (xx-nx)**2+(yy-ny)**2 < (3.*fwhm)**2
            catalogue_masked += int((hit & usable).sum())
            usable &= ~hit
        projection = (dx*ux+dy*uy)/np.maximum(np.hypot(dx,dy), 1.e-12)
        toward = projection >= cos_limit
        away = projection <= -cos_limit
        # Pair a surviving toward-side pixel with the nearest original
        # pixel at its antipodal location, using exactly the same native
        # neighbour-core and hard-pixel masks. Reflection of the integer
        # pixel lattice is one-to-one, including at subpixel star centres.
        reflected_x = np.rint(2.*cx-xx).astype(np.int64)-x0
        reflected_y = np.rint(2.*cy-yy).astype(np.int64)-y0
        on_patch = ((reflected_x >= 0) & (reflected_y >= 0)
                    & (reflected_x < patch.shape[1])
                    & (reflected_y < patch.shape[0]))
        mirrored_x = np.clip(reflected_x, 0, patch.shape[1]-1)
        mirrored_y = np.clip(reflected_y, 0, patch.shape[0]-1)
        mirrored_usable = usable[mirrored_y, mirrored_x] & on_patch
        mirrored_radius = radius[mirrored_y, mirrored_x]
        mirrored_away = away[mirrored_y, mirrored_x]
        reflected_pixels = patch[mirrored_y, mirrored_x]

        def median_in(low, high, direction):
            selected = (usable & direction & (radius >= low)
                        & (radius < high))
            n = int(selected.sum())
            return (float(np.median(patch[selected])) if n >= 30 else None,
                    n)

        sky_24, n_24 = median_in(24., 26., away)
        sky_26, n_26 = median_in(26., 30., away)
        report = []
        for low, high in rings:
            facing, n_facing = median_in(low, high, toward)
            opposite, n_opposite = median_in(low, high, away)
            difference = (float(facing-opposite)
                          if facing is not None and opposite is not None
                          else None)
            in_ring = (radius >= low) & (radius < high)
            total_facing = int((in_ring & toward).sum())
            total_away = int((in_ring & away).sum())
            coverage_facing = (float(n_facing/total_facing)
                               if total_facing else None)
            coverage_away = (float(n_opposite/total_away)
                             if total_away else None)
            mask_asymmetric = (
                coverage_facing is None or coverage_away is None
                or min(coverage_facing, coverage_away) <
                P21_MIN_CONE_GEOMETRIC_COVERAGE
                or abs(coverage_facing-coverage_away) >
                P21_MAX_CONE_COVERAGE_DIFFERENCE)
            matched = (in_ring & toward & usable & mirrored_usable
                       & mirrored_away
                       & (mirrored_radius >= low)
                       & (mirrored_radius < high))
            n_matched = int(matched.sum())
            paired_difference = (
                float(np.median((patch-reflected_pixels)[matched]))
                if n_matched >= P21_MIN_ANTIPODAL_PAIRS else None)
            report.append({
                "radius_in_measured_fwhm": [low,high],
                "n_geometrically_available_facing_pixels": total_facing,
                "n_geometrically_available_opposite_pixels": total_away,
                "n_original_valid_facing_pixels": n_facing,
                "n_original_valid_opposite_pixels": n_opposite,
                "facing_cone_valid_geometric_fraction": coverage_facing,
                "opposite_cone_valid_geometric_fraction": coverage_away,
                "coverage_difference_facing_minus_opposite": (
                    float(coverage_facing-coverage_away)
                    if coverage_facing is not None and coverage_away is not None
                    else None),
                "mask_coverage_interpretation": (
                    "SPARSE_OR_ASYMMETRIC_RAW_SECTOR_MEDIANS_NOT_COMPARABLE"
                    if mask_asymmetric else
                    "MASK_GEOMETRY_COMPARABLE_NOT_SCIENCE_VALIDATION"),
                "n_valid_antipodal_facing_to_opposite_pixel_pairs": n_matched,
                "antipodal_paired_difference_median_nMgy_per_pixel":
                    paired_difference,
                "antipodal_paired_native_geometry_not_calibrated": True,
                "toward_median_nMgy_per_pixel": facing,
                "away_median_nMgy_per_pixel": opposite,
                "toward_minus_away_nMgy_per_pixel": difference,
                "toward_relative_to_outer_sky_24_26_nMgy_per_pixel": (
                    float(facing-sky_24)
                    if facing is not None and sky_24 is not None else None),
                "toward_relative_to_outer_sky_26_30_nMgy_per_pixel": (
                    float(facing-sky_26)
                    if facing is not None and sky_26 is not None else None),
                "sky_reference_difference_not_statistical_error": (
                    float(sky_24-sky_26)
                    if sky_24 is not None and sky_26 is not None else None),
            })
        return {
            "status": "FIXED_DIRECTIONAL_ANNULI_DESCRIPTIVE_ONLY",
            "n_catalog_neighbour_core_pixels_excluded": catalogue_masked,
            "outer_reference_away_24_26_nMgy_per_pixel": sky_24,
            "outer_reference_away_26_30_nMgy_per_pixel": sky_26,
            "outer_reference_away_pixels": [n_24,n_26],
            "n_rings_measured_with_both_cones": sum(
                row["toward_minus_away_nMgy_per_pixel"] is not None
                for row in report),
            "per_ring": report,
        }

    star = measure(sx, sy, neighbours)
    if star["status"] != "FIXED_DIRECTIONAL_ANNULI_DESCRIPTIVE_ONLY":
        return dict(base,status=star["status"],source=star)
    controls = []
    for tile_y in range(4):
        for tile_x in range(4):
            cx = int(round((tile_x+.5)*columns/4.))
            cy = int(round((tile_y+.5)*rows/4.))
            entry = {"tile": [tile_x,tile_y],
                     "center_native_xy_pix": [cx,cy]}
            if (math.hypot(cx-tx,cy-ty) < 60.*fwhm or
                    math.hypot(cx-sx,cy-sy) < 60.*fwhm):
                entry["status"] = "PROJECT_TARGET_OR_STAR_GEOMETRY_EXCLUDED"
            elif any(math.hypot(cx-nx,cy-ny) < 3.*fwhm
                     for _,nx,ny in neighbours):
                entry["status"] = "CATALOG_CENTRE_AT_CONTROL_ORIGIN"
            else:
                # Treat the selected star as an additional contaminant
                # at control locations; never mask the source itself
                # when measuring its own radial profile above.
                entry.update(measure(
                    cx,cy,neighbours+[(-1,sx,sy)]))
            controls.append(entry)
    comparison = []
    for ring_index, (low,high) in enumerate(rings):
        star_row = star["per_ring"][ring_index]
        controls_here = [
            row["per_ring"][ring_index]["toward_minus_away_nMgy_per_pixel"]
            for row in controls
            if row["status"] == "FIXED_DIRECTIONAL_ANNULI_DESCRIPTIVE_ONLY"
            and row["per_ring"][ring_index][
                "toward_minus_away_nMgy_per_pixel"] is not None
        ]
        paired_controls_here = [
            row["per_ring"][ring_index][
                "antipodal_paired_difference_median_nMgy_per_pixel"]
            for row in controls
            if row["status"] == "FIXED_DIRECTIONAL_ANNULI_DESCRIPTIVE_ONLY"
            and row["per_ring"][ring_index][
                "antipodal_paired_difference_median_nMgy_per_pixel"] is not None
        ]
        comparison.append({
            "radius_in_measured_fwhm": [low,high],
            "source_toward_minus_away_nMgy_per_pixel":
                star_row["toward_minus_away_nMgy_per_pixel"],
            "source_mask_coverage_interpretation":
                star_row["mask_coverage_interpretation"],
            "source_antipodal_paired_difference_median_nMgy_per_pixel":
                star_row["antipodal_paired_difference_median_nMgy_per_pixel"],
            "source_n_antipodal_pixel_pairs":
                star_row["n_valid_antipodal_facing_to_opposite_pixel_pairs"],
            "n_catalog_geometry_fixed_controls_with_valid_cones":
                len(controls_here),
            "median_abs_control_asymmetry_nMgy_per_pixel": (
                float(np.median(np.abs(controls_here)))
                if controls_here else None),
            "max_abs_control_asymmetry_nMgy_per_pixel": (
                float(np.max(np.abs(controls_here)))
                if controls_here else None),
            "n_fixed_controls_with_valid_antipodal_pairs":
                len(paired_controls_here),
            "median_abs_paired_control_difference_nMgy_per_pixel": (
                float(np.median(np.abs(paired_controls_here)))
                if paired_controls_here else None),
            "max_abs_paired_control_difference_nMgy_per_pixel": (
                float(np.max(np.abs(paired_controls_here)))
                if paired_controls_here else None),
            "not_a_wing_detection_or_significance": True,
        })
    return dict(
        base,
        status="SAME_PROJECT_DIRECTIONAL_SKY_SENSITIVITY_ONLY",
        source_to_target_centre_measured_fwhm=float(dist/fwhm),
        source=star,
        n_source_rings_with_mask_geometry_comparable=sum(
            row["mask_coverage_interpretation"] ==
            "MASK_GEOMETRY_COMPARABLE_NOT_SCIENCE_VALIDATION"
            for row in star["per_ring"]),
        n_source_rings_with_antipodal_matched_pixels=sum(
            row["antipodal_paired_difference_median_nMgy_per_pixel"]
            is not None for row in star["per_ring"]),
        controls=controls,
        n_fixed_control_tiles=16,
        n_controls_with_at_least_one_measured_ring=sum(
            r.get("n_rings_measured_with_both_cones",0)>0
            for r in controls),
        comparison_by_ring=comparison,
        fixed_control_geometry_not_selected_on_residuals=True,
        all_control_centres_from_current_image_only=True,
        coadd_controls_share_sky_and_pixel_correlations=True,
        no_empirical_psf_wing_extrapolation=True,
    )
