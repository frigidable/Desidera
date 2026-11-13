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
    # The FWHM finder prior must be finite/positive and not a Boolean.
    # No invalid prior is allowed in a completed PSF QA receipt.
    guess = None
    if not isinstance(guessed_fwhm_pix, (bool, np.bool_, str)):
        try:
            candidate = float(guessed_fwhm_pix)
            if math.isfinite(candidate) and candidate > 0:
                guess = candidate
        except (TypeError, ValueError, OverflowError):
            pass
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
        "guessed_detection_fwhm_pix_not_a_measurement": guess,
    }
    if guess is None:
        return dict(result, status="INVALID_DETECTION_FWHM_PRIOR")
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
        measured_over_guessed_fwhm_ratio=float(measured / guess),
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
            "Predeclared error budget tested with fixed brightness/position-stratified withheld stars and injections from this brick"),
        "wide_angle_wing_measurement": row(
            "HOLD",
            "Only a finite 81px-class coadd kernel; measured EEF denominator is truncated",
            "Disjoint bright stars within the SAME native brick, wide-radius coadd profiles with conditional sky separation and explicit wing-flux uncertainty or unresolved bound"),
        "sky_versus_psf_wing_identifiability": row(
            "HOLD",
            "Stellar annulus and coadd sky can absorb outer-wing flux",
            "Within-brick blank-sky and neighbor-control tests, using original pixels; quantify sky/wing degeneracy and photometric bias without changing brick"),
        "spatial_representativity": row(
            "HOLD",
            "The original training/holdout split is too small for brightness-matched within-brick spatial representativity",
            "Position-separated and brightness-matched withheld stars within the SAME native brick; map spatial residual variation without requiring other bricks"),
        "resampled_coadd_noise_covariance": row(
            "HOLD",
            "Diagonal invvar alone cannot validate correlated noise",
            "Native-brick blank-sky empirical covariance or a demonstrably calibrated correlated-noise uncertainty model"),
        "galaxy_injection_recovery": row(
            "HOLD",
            "No independent known-input LSB-galaxy bias budget",
            "Predeclared surface-brightness, radius, ellipticity and total-flux bias tolerances tested on injections into this native brick"),
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



# The engineering stage can be finalized without promoting an unmeasured
# outer halo to a physical point-spread kernel or demanding data from
# another brick. This is a read-only terminal characterization of the
# *same* brick, not a new wing-fit experiment or a science gate override.
SINGLE_BRICK_CLOSURE_SCHEMA = "desidera.psf0-single-brick-empirical-closure/1"
SINGLE_BRICK_RINGS = (
    (4., 6.), (6., 8.), (8., 12.), (12., 16.),
    (16., 20.), (20., 22.),
)
SINGLE_BRICK_COHORT_ROLES = (
    "P16C_LEGACY_EXPLORATORY",
    "NEW_CATALOG_EXPLORATORY",
    "NEW_CATALOG_HOLDOUT",
)
SINGLE_BRICK_EXPECTED_COHORT = (16, 24, 8)


def _brick_finite_median(values, *, absolute=False):
    """Return a finite descriptive median; missing values never become zero."""
    finite = [
        abs(float(v)) if absolute else float(v)
        for v in values if isinstance(v, (int, float, np.number))
        and not isinstance(v, (bool, np.bool_))
        and math.isfinite(float(v))
    ]
    return float(np.median(finite)) if finite else None


def finalize_single_brick_empirical_psf(
    *, dataset_id, image_shape, core_validation,
    p17_qa, p18_qa, verified_native_sources,
):
    """Close PSF-0 engineering for a single native brick without adoption.

    Combines already measured P17 and P18 *matched identities* only.
    The returned radial medians are contaminated sky-conditional
    diagnostics, NOT a normalized physical wing PSF or flux upper limits.
    Completeness of this inventory differs from photometric science
    accuracy; independent bricks, new datasets and selected ring
    thresholds are never requested or invoked.
    """
    result = {
        "schema": SINGLE_BRICK_CLOSURE_SCHEMA,
        "status": "INCOMPLETE",
        "dataset_id": dataset_id,
        "fixed_native_brick_only": True,
        "engineering_stage_finished": False,
        "science_adopted": False,
        "wide_wing_psf_kernel_constructed": False,
        "wide_wing_detection_significance_claimed": False,
        "wide_wing_photometric_correction_permitted": False,
        "no_other_brick_required_or_used": True,
        "outer_profiles_sky_conditional_not_physical_ee_fraction": True,
        "observed_coadd_psf_not_raw_optical_instrument_psf": True,
        "finite_epsf_normalization_not_true_total_stellar_flux": True,
        "stellar_identities_not_reranked_or_replaced": True,
        "statistics_correlated_within_the_same_coadd": True,
        "coverage_thresholds_are_descriptive_not_science_gates": True,
    }
    if not isinstance(dataset_id, str) or dataset_id.count("/") != 3:
        return dict(result, status="INVALID_SINGLE_BRICK_DATASET_ID")
    if (not isinstance(image_shape, (tuple, list))
            or len(image_shape) != 2
            or any(type(x) is not int or x <= 0 for x in image_shape)):
        return dict(result, status="INVALID_ORIGINAL_BRICK_SHAPE")
    if not verified_native_sources:
        return dict(result, status="ORIGINAL_NATIVE_SOURCE_HASHES_UNVERIFIED")
    if (not isinstance(core_validation, dict)
            or core_validation.get("status") != "SCIENCE_GATE_HOLD"):
        return dict(result, status="MISSING_FROZEN_CORE_SCIENCE_GATE")
    geom = core_validation.get("core_geometry", {})
    if geom.get("status") != "FINITE_NATIVE_CORE_GEOMETRY_MEASURED":
        return dict(result, status="NO_MEASURED_FINITE_CORE")
    fwhm = geom.get("measured_core_fwhm_axis_geometric_mean_native_pix")
    full_radius = geom.get("fully_contained_radial_radius_limit_native_pix")
    if (not isinstance(fwhm, (int, float))
            or not math.isfinite(fwhm) or fwhm <= 0
            or not isinstance(full_radius, (int, float))
            or not math.isfinite(full_radius) or full_radius <= 0):
        return dict(result, status="INVALID_FROZEN_CORE_GEOMETRY")
    if (not isinstance(p17_qa, dict)
            or p17_qa.get("status") !=
            "FIXED_EXPANDED_COHORT_MASK_FEASIBILITY_ONLY"
            or not isinstance(p18_qa, dict)
            or p18_qa.get("status") !=
            "FIXED_P17_COHORT_MORPHOLOGY_SKY_DIAGNOSTIC_ONLY"):
        return dict(result, status="FIXED_P17_P18_EVIDENCE_MISSING")
    if (not isinstance(p17_qa.get("native_catalog_sha256"), str)
            or p17_qa["native_catalog_sha256"] !=
            p18_qa.get("native_catalog_sha256")):
        return dict(result, status="NATIVE_TRACTOR_CATALOG_HASH_MISMATCH")
    p17 = p17_qa.get("per_star", [])
    p18 = p18_qa.get("per_star", [])
    if not isinstance(p17, list) or not isinstance(p18, list):
        return dict(result, status="INVALID_FROZEN_STAR_LISTS")
    if len(p17) != 48 or len(p18) != 48:
        return dict(result, status="INCOMPLETE_48_STAR_COHORT")
    by_rank = {}
    for old in p17:
        rank = old.get("catalog_rank")
        if (type(rank) is not int or rank < 1 or rank > 48
                or rank in by_rank):
            return dict(result, status="DUPLICATE_OR_INVALID_P17_RANK")
        by_rank[rank] = old
    if set(by_rank) != set(range(1, 49)):
        return dict(result, status="INCOMPLETE_OR_NONCONTIGUOUS_RANKS")
    matched = []
    encountered = set()
    for new in p18:
        rank = new.get("catalog_rank")
        if type(rank) is not int or rank not in by_rank or rank in encountered:
            return dict(result, status="DUPLICATE_OR_UNKNOWN_P18_RANK")
        encountered.add(rank)
        old = by_rank[rank]
        if any(new.get(k) != old.get(k) for k in (
                "gaia_source_id", "source_id", "cohort_role")):
            return dict(result, status="P18_P17_STAR_IDENTITY_MISMATCH")
        if len(old.get("mask_variants", [])) != 4:
            return dict(result, status="P17_NOT_ALL_FOUR_MASKS")
        if (len(new.get("radial_profile", [])) != len(SINGLE_BRICK_RINGS)
                or len(old["mask_variants"][0].get("radial_profile", []))
                != len(SINGLE_BRICK_RINGS)):
            return dict(result, status="MISSING_FIXED_RADIAL_BINS")
        for index, ring in enumerate(SINGLE_BRICK_RINGS):
            for entry in (new["radial_profile"][index],
                          old["mask_variants"][0]["radial_profile"][index]):
                if (entry.get("radius_inner_in_measured_fwhm"),
                        entry.get("radius_outer_in_measured_fwhm")) != ring:
                    return dict(result, status="RADIAL_RING_BOUNDARIES_CHANGED")
        xy = old.get("native_xy_pix")
        if (not isinstance(xy, (tuple, list)) or len(xy) != 2
                or not all(isinstance(x, (int, float)) and
                           math.isfinite(x) for x in xy)
                or not (0 <= xy[0] < image_shape[1]
                        and 0 <= xy[1] < image_shape[0])):
            return dict(result, status="INVALID_NATIVE_STAR_LOCATION")
        matched.append((old,new))
    if len({old["gaia_source_id"] for old, _ in matched}) != 48:
        return dict(result, status="DUPLICATE_GAIA_IDENTITIES")
    if tuple(sum(old["cohort_role"] == role for old,_ in matched)
             for role in SINGLE_BRICK_COHORT_ROLES) != (
             SINGLE_BRICK_EXPECTED_COHORT):
        return dict(result, status="FROZEN_COHORT_ROLES_CHANGED")
    if (p17_qa.get("n_legacy_baseline_reproduced") != 16
            or p17_qa.get("n_multimask_stars_measured") != 48):
        return dict(result, status="P16C_BASELINE_OR_P17_RUN_INCOMPLETE")
    if any(old.get("status") != "MULTIMASK_DESCRIPTIVE_ONLY"
           for old, _ in matched):
        return dict(result, status="P17_FIXED_STAR_NOT_MEASURED")

    subsets = {"all_48": matched}
    for label in SINGLE_BRICK_COHORT_ROLES:
        subsets[label] = [
            (old, new) for old,new in matched
            if old["cohort_role"] == label
        ]
    population = {}
    radial_nulls = 0
    radial_no_pixels = 0
    for label, rows in subsets.items():
        radial = []
        for index,(lo,hi) in enumerate(SINGLE_BRICK_RINGS):
            old_r = [old["mask_variants"][0]["radial_profile"][index]
                     for old,new in rows]
            new_r = [new["radial_profile"][index] for old,new in rows]
            old_signal = [r.get("sky_b_subtracted_median_over_peak")
                          for r in old_r]
            new_signal = [r.get("sky_b_subtracted_median_over_peak")
                          for r in new_r]
            paired = [
                abs(float(a)-float(b))
                for a,b in zip(old_signal, new_signal)
                if a is not None and b is not None
                and math.isfinite(float(a)) and math.isfinite(float(b))
            ]
            n_available = sum(
                isinstance(x, (int, float)) and math.isfinite(x)
                for x in new_signal
            )
            if label == "all_48":
                radial_nulls += len(rows)-n_available
                radial_no_pixels += sum(
                    r.get("n_unmasked_native_pixels") == 0 for r in new_r)
            qualified = sum(
                r.get("unmasked_fraction") is not None and
                r["unmasked_fraction"] >= .5 and
                r.get("n_octants_with_at_least_25_pixels",0) >= 6
                for r in new_r
            )
            spread = []
            for r in new_r:
                octants = [
                    float(v) for v in
                    r.get("sky_b_subtracted_octant_medians_over_peak", [])
                    if isinstance(v, (int,float)) and math.isfinite(v)
                ]
                spread.append(
                    max(octants)-min(octants) if len(octants) >= 4 else None
                )
            radial.append({
                "radius_range_in_measured_fwhm": [lo,hi],
                "n_stars": len(rows),
                "n_nonnull_p18_radial_residuals": n_available,
                "n_geometric_coverage_illustration_only": qualified,
                "p18_median_unmasked_annular_fraction": _brick_finite_median(
                    [r.get("unmasked_fraction") for r in new_r]),
                "p18_median_signed_residual_over_native_peak":
                    _brick_finite_median(new_signal),
                "p18_median_abs_residual_over_native_peak":
                    _brick_finite_median(new_signal, absolute=True),
                "p17_3fwhm_median_signed_residual_over_native_peak":
                    _brick_finite_median(old_signal),
                "median_abs_p18_vs_p17_3fwhm_shift_over_native_peak":
                    _brick_finite_median(paired),
                "median_p18_max_minus_min_octant_over_native_peak":
                    _brick_finite_median(spread),
                "no_significance_or_true_psf_flux_inferred": True,
            })
        sky = {
            "median_abs_sky_a_minus_sky_b_over_native_peak":
                _brick_finite_median(
                    [new.get("sky_a_minus_sky_b_over_peak")
                     for _,new in rows], absolute=True),
            "median_max_abs_sky_b_octant_deletion_over_native_peak":
                _brick_finite_median(
                    [new.get(
                        "sky_b_max_abs_leave_one_octant_change_over_peak")
                     for _,new in rows]),
            "eight_deletions_not_independent_error_estimates": True,
        }
        population[label] = {
            "n_stars": len(rows),
            "radial_profiles": radial,
            "sky_sensitivity": sky,
        }
    # Fixed x/y half-brick geometry is catalog-only. Sparse quadrants are
    # reported rather than pooled into fictitious location-independent QA.
    quadrants = []
    for east in (False, True):
        for north in (False, True):
            rows = [(old,new) for old,new in matched
                    if (old["native_xy_pix"][0] >= image_shape[1]/2.) == east
                    and (old["native_xy_pix"][1] >= image_shape[0]/2.)
                    == north]
            bright = [(old,new) for old,new in rows
                      if old["cohort_role"] ==
                      "P16C_LEGACY_EXPLORATORY"]
            quadrants.append({
                "brick_pixel_x_upper_half": east,
                "brick_pixel_y_upper_half": north,
                "n_fixed_stars": len(rows),
                "n_legacy_bright_stars": len(bright),
                "median_bright_p18_4to6_fwhm_residual_over_native_peak":
                    _brick_finite_median([
                        new["radial_profile"][0].get(
                            "sky_b_subtracted_median_over_peak")
                        for _,new in bright]),
                "brightness_matched_spatial_qa_validated": False,
            })

    # Explicitly bound what the result does and does not mean. A PSF
    # "engineering stage done" is not a science-adopted physical halo.
    return dict(
        result,
        status="BRICK_ONLY_EMPIRICAL_CORE_AND_WING_CONSTRAINTS_COMPLETE",
        engineering_stage_finished=True,
        native_catalog_sha256=p17_qa["native_catalog_sha256"],
        n_frozen_stars=len(matched),
        n_original_kernel_training_stars=core_validation.get(
            "n_train_stars"),
        n_original_kernel_holdout_stars=core_validation.get(
            "n_fixed_holdout_stars"),
        finite_core_measured_fwhm_native_pix=float(fwhm),
        finite_core_kernel_full_radial_support_native_pix=float(full_radius),
        physical_core_flux_normalization_known=False,
        outer_diagnostic_radii_fwhm=[4.,22.],
        outer_conditioning_sky_radii_fwhm=[22.,30.],
        n_p18_radial_annuli_missing_residual=radial_nulls,
        n_p18_radial_annuli_with_no_pixels=radial_no_pixels,
        cohorts=population,
        same_brick_catalog_pixel_quadrants=quadrants,
        scientifically_validated_extended_wings=False,
        final_outcome=(
            "Finite native 2D empirical core + sky-conditional stellar "
            "profiles and masking sensitivity. No wide-angle physical "
            "PSF inferred and no science photometric correction."),
    )
