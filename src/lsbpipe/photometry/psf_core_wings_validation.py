"""Finite PSF-0 core and wing-readiness checks; never adopt a science PSF.

PSF half-maximum crossings use native-grid axial sections through the
original peak, and encircled fractions refer ONLY to the truncated kernel.
These are diagnostics, not calibrated uncertainty estimates or evidence of
true wide-angle wings. No observation pixels or PSF kernels are modified.
"""
from __future__ import annotations

import math

import numpy as np

SCHEMA = "desidera.psf0-core-wings-validation/1"
RADII_IN_MEASURED_FWHM = (1., 2., 3., 4.)


def _half_max_crossing(profile, peak_index):
    """First outward crossing, linearly interpolated in native pixels."""
    center = float(profile[peak_index])
    half = 0.5 * center
    distances = []
    for direction in (-1, +1):
        crossing = None
        for step in range(1, len(profile)):
            j = peak_index + direction * step
            if j < 0 or j >= len(profile):
                break
            below = float(profile[j])
            if below <= half:
                inside = float(profile[j-direction])
                if not (inside > half and inside > below):
                    break
                crossing = float(step - 1 + (inside-half)/(inside-below))
                break
        distances.append(crossing)
    return distances


def measure_finite_native_psf_core(kernel, *, guessed_fwhm_pix):
    """Measure core-axis half maxima and finite-kernel EEF without fit priors."""
    result = {
        "schema": "desidera.psf0-native-core-geometry/1",
        "status": "UNMEASURED",
        "method": "NATIVE_PEAK_AXIS_HALF_MAX_LINEAR_INTERPOLATION",
        "measured_from_final_nonnegative_normalized_kernel": True,
        "fwhm_is_kernel_axis_intercept_not_a_2d_model_fit": True,
        "subpixel_psf_variation_not_measured": True,
        "uncertainty_not_calibrated": True,
        "finite_kernel_encircled_fraction_not_true_infinite_total": True,
        "no_psf_model_change_or_science_adoption": True,
        "guessed_detection_fwhm_pix_not_a_measurement": float(
            guessed_fwhm_pix),
    }
    image = np.asarray(kernel, dtype=float)
    if (image.ndim != 2 or any(v < 15 or v % 2 == 0 for v in image.shape)
            or not np.all(np.isfinite(image)) or np.any(image < 0)
            or not np.isclose(float(image.sum()), 1., rtol=1e-8,
                              atol=1e-9)):
        return dict(result, status="INVALID_FINITE_NORMALIZED_KERNEL")
    ny, nx = image.shape
    iy, ix = np.unravel_index(np.argmax(image), image.shape)
    center = ((ny-1)//2, (nx-1)//2)
    result.update(
        kernel_shape_native_pix=[int(ny), int(nx)],
        peak_native_pixel_yx=[int(iy), int(ix)],
        grid_center_native_pixel_yx=[int(center[0]), int(center[1])],
        peak_offset_from_native_grid_center_pix=float(math.hypot(
            ix-center[1], iy-center[0])),
        peak_value_normalized_native_pixel=float(image[iy, ix]),
    )
    if (iy == 0 or ix == 0 or iy == ny-1 or ix == nx-1
            or image[iy, ix] <= 0):
        return dict(result, status="PEAK_AT_EDGE_OR_NONPOSITIVE")
    x_cross = _half_max_crossing(image[iy, :], ix)
    y_cross = _half_max_crossing(image[:, ix], iy)
    result["axial_half_peak_crossing_distances_pix"] = {
        "x_negative": x_cross[0], "x_positive": x_cross[1],
        "y_negative": y_cross[0], "y_positive": y_cross[1],
    }
    if any(v is None for v in x_cross+y_cross):
        return dict(result, status="UNBRACKETED_CORE_HALF_MAXIMUM")
    width_x = float(sum(x_cross))
    width_y = float(sum(y_cross))
    if not (width_x > 0 and width_y > 0
            and math.isfinite(width_x) and math.isfinite(width_y)):
        return dict(result, status="INVALID_CORE_HALF_MAXIMUM_WIDTH")
    measured = float(math.sqrt(width_x*width_y))
    # Center the radial apertures at the actual kernel peak. Because
    # the native PSF grid is finite, even all-included-patch fractions
    # cannot establish a true infinite-radius encircled energy curve.
    yy, xx = np.indices(image.shape)
    radius = np.hypot(xx-ix, yy-iy)
    enclosed = []
    full_radius = float(min(ix, nx-1-ix, iy, ny-1-iy))
    for multiple in RADII_IN_MEASURED_FWHM:
        r = float(multiple*measured)
        fully_contained = bool(r <= full_radius)
        # Sum in float64, divide by the *actual* finite-kernel
        # normalization and remove only machine-roundoff overshoots.
        # The [0,1] clamp is a numerical probability bound, NOT a
        # correction for unknown astrophysical power in outer wings.
        included_fraction = None
        if fully_contained:
            raw = float(image[radius <= r].sum()/image.sum())
            included_fraction = float(min(1., max(0., raw)))
        enclosed.append({
            "radius_in_measured_axial_geometric_fwhm":multiple,
            "radius_native_pix":r,
            "enclosed_fraction_of_finite_kernel":included_fraction,
            "full_circular_aperture_inside_kernel":fully_contained,
            "only_float_roundoff_clamped_to_unit_interval":True,
            "not_total_physical_psf_flux":True,
        })
    return dict(
        result,
        status="FINITE_NATIVE_CORE_GEOMETRY_MEASURED",
        fwhm_x_native_pix=width_x,
        fwhm_y_native_pix=width_y,
        measured_core_fwhm_axis_geometric_mean_native_pix=measured,
        measured_over_guessed_fwhm_ratio=(
            measured / float(guessed_fwhm_pix)
            if math.isfinite(guessed_fwhm_pix)
            and guessed_fwhm_pix>0 else None),
        fwhm_axis_ratio_max_over_min=float(
            max(width_x,width_y)/min(width_x,width_y)),
        fully_contained_radial_radius_limit_native_pix=full_radius,
        finite_kernel_encircled_flux=enclosed,
    )


def evaluate_core_wings_readiness(
    *, core_geometry, n_train, n_holdout, coverage,
    holdout_growth_summary, preclip_negative_mass_fraction,
    verified_native_source_hashes,
):
    """PASS means evidence *exists*, never that science accuracy passed.

    Calibration limits are deliberately NOT chosen by looking at the
    current Leo P holdouts. Science readiness requires independently
    justified tolerance for source bias, noise, coadd covariance, wings,
    spatial representativity and prospective source injections.
    """
    coverage = coverage if isinstance(coverage, dict) else {}
    growth = (
        holdout_growth_summary if isinstance(holdout_growth_summary,dict)
        else {})
    geom = core_geometry if isinstance(core_geometry,dict) else {}
    complete = (
        coverage.get("all_fixed_holdouts_and_annuli_scored") is True
        and coverage.get("n_expected_holdout_stars")==n_holdout
        and coverage.get("n_successfully_scored_stars")==n_holdout
        and n_holdout>=1
    )
    def row(status, evidence, requirement):
        return {
            "status": status,
            "current_evidence": evidence,
            "requirement_before_science_pass": requirement,
        }
    stages = {
        "original_native_input_provenance": row(
            "PASS" if verified_native_source_hashes else "HOLD",
            "Original local FITS hash quartet verified"
            if verified_native_source_hashes else "No verified original FITS quartet",
            "Local original survey image, inverse-variance, NEXP and maskbits with SHA256 lineage"),
        "finite_kernel_numerical_geometry": row(
            "PASS" if geom.get("status")=="FINITE_NATIVE_CORE_GEOMETRY_MEASURED"
            else "HOLD",
            geom.get("status","NO_CORE_GEOMETRY"),
            "Finite normalized kernel with measurable, bracketed native-axis core FWHM"),
        "predeclared_holdout_coverage": row(
            "PASS" if complete else "HOLD",
            f"{coverage.get('n_successfully_scored_stars','?')}/{n_holdout} unchanged heldout stars",
            "All originally reserved holdout stars scored on original masks and predeclared annuli"),
        "core_accuracy": row(
            "HOLD",
            "Core geometry/holdout residuals are descriptive; no prospective bias tolerance",
            "Externally justified core flux/shape error budget, stable across brightness, position and independent stars"),
        "wide_angle_wing_measurement": row(
            "HOLD",
            "Only a finite 81px-class coadd kernel; measured EEF denominator is truncated",
            "Independent bright-star large-radius measurements and calibrated outer-wing model with sky separation"),
        "sky_versus_psf_wing_identifiability": row(
            "HOLD",
            "Stellar annulus and coadd sky can absorb outer-wing flux",
            "Independent sky/PSF-wing degeneracy stress and bias bounds"),
        "spatial_representativity": row(
            "HOLD",
            "Single native brick global PSF not checked at independent positions",
            "Independent brightness- and position-matched stars plus survey PSF-map comparison"),
        "resampled_coadd_noise_covariance": row(
            "HOLD",
            "Diagonal invvar alone cannot validate correlated noise",
            "Measured/propagated covariance or empirically validated correlated-noise uncertainty model"),
        "galaxy_injection_recovery": row(
            "HOLD",
            "No independent known-input LSB-galaxy bias budget",
            "Predeclared surface-brightness, radius, ellipticity and total-flux bias tolerances tested on injections"),
    }
    return {
        "schema": SCHEMA,
        "status": "SCIENCE_GATE_HOLD",
        "science_adopted": False,
        "n_train_stars": int(n_train),
        "n_fixed_holdout_stars": int(n_holdout),
        "core_geometry": geom,
        "preclipping_negative_absolute_mass_fraction":(
            float(preclip_negative_mass_fraction)
            if preclip_negative_mass_fraction is not None else None),
        "preclipping_negative_mass_not_a_physical_psf_probability":True,
        "holdout_growth_summary_existing_descriptive_only": {
            "n_fixed_holdouts": growth.get("n_fixed_holdouts"),
            "radii": [
                {
                    "radius_in_guessed_fwhm": item.get(
                        "radius_in_guessed_fwhm"),
                    "median_absolute_growth_excess_over_core":
                        item.get("median_absolute_growth_excess_over_core"),
                    "n_full_aperture_holdouts":
                        item.get("n_full_aperture_holdouts"),
                }
                for item in growth.get("radii", [])
            ],
        },
        "stages": stages,
        "passes_are_evidence_completeness_not_scientific_accuracy":True,
        "no_threshold_optimized_on_original_holdouts":True,
        "no_automatic_recipe_selection":True,
        "no_physical_wing_flux_inference_from_finite_kernel":True,
        "no_psf_kernel_change":True,
        "no_science_adoption_without_photometric_bias_calibration":True,
    }
