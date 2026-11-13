"""PSF Core & Wings Validation: numeric evidence PASS != science adoption."""
import numpy as np
import pytest

from lsbpipe.photometry.psf_core_wings_validation import (
    measure_finite_native_psf_core, evaluate_core_wings_readiness,
)


def gaussian_kernel(size=81, sigma_x=3., sigma_y=3.):
    yy, xx = np.indices((size, size), dtype=float)
    r = (size-1)/2
    image = np.exp(-.5*(((xx-r)/sigma_x)**2 +
                        ((yy-r)/sigma_y)**2))
    return image/image.sum()


def complete_coverage():
    return {
        "all_fixed_holdouts_and_annuli_scored": True,
        "n_expected_holdout_stars": 3,
        "n_successfully_scored_stars": 3,
    }


def growth():
    return {
        "n_fixed_holdouts":3,
        "radii":[{
            "radius_in_guessed_fwhm":4.,
            "median_absolute_growth_excess_over_core":.095869648,
            "n_full_aperture_holdouts":3,
        }]
    }


def test_symmetric_native_gaussian_fwhm_is_measured_not_guessed():
    kernel=gaussian_kernel()
    original=kernel.copy()
    result=measure_finite_native_psf_core(
        kernel,guessed_fwhm_pix=6.57)
    exact=2*np.sqrt(2*np.log(2))*3.
    assert result["status"]=="FINITE_NATIVE_CORE_GEOMETRY_MEASURED"
    assert result["measured_core_fwhm_axis_geometric_mean_native_pix"]==(
        pytest.approx(exact,rel=.015))
    assert result["fwhm_x_native_pix"]==pytest.approx(
        result["fwhm_y_native_pix"],rel=1e-12)
    assert result["measured_over_guessed_fwhm_ratio"]==pytest.approx(
        result["measured_core_fwhm_axis_geometric_mean_native_pix"]/6.57)
    assert result["peak_offset_from_native_grid_center_pix"]==0
    assert result["fwhm_axis_ratio_max_over_min"]==pytest.approx(1.)
    assert result["uncertainty_not_calibrated"]
    assert result["fwhm_is_kernel_axis_intercept_not_a_2d_model_fit"]
    assert result["no_psf_model_change_or_science_adoption"]
    rows=result["finite_kernel_encircled_flux"]
    assert len(rows)==4
    fractions=[row["enclosed_fraction_of_finite_kernel"] for row in rows]
    assert all(row["full_circular_aperture_inside_kernel"] for row in rows)
    assert all(0<=x<=1 for x in fractions)
    assert fractions==sorted(fractions)
    assert fractions[-1]<=1.  # finite-precision Gaussian tail can underflow
    np.testing.assert_array_equal(kernel,original)


def test_elliptical_core_width_and_censored_outer_growth():
    image=gaussian_kernel(size=81,sigma_x=8.,sigma_y=4.)
    result=measure_finite_native_psf_core(image,guessed_fwhm_pix=6.57)
    assert result["status"]=="FINITE_NATIVE_CORE_GEOMETRY_MEASURED"
    assert result["fwhm_x_native_pix"]>result["fwhm_y_native_pix"]
    assert result["fwhm_axis_ratio_max_over_min"]==pytest.approx(
        2.,rel=.03)
    # 4 x geometric FWHM is beyond the smallest edge distance.
    assert result["finite_kernel_encircled_flux"][-1][
        "full_circular_aperture_inside_kernel"] is False
    assert result["finite_kernel_encircled_flux"][-1][
        "enclosed_fraction_of_finite_kernel"] is None
    assert result["finite_kernel_encircled_fraction_not_true_infinite_total"]


@pytest.mark.parametrize("bad", [
    np.zeros((20,20)),
    np.ones((81,81))/80.**2,
    -gaussian_kernel(),
    np.full((81,81),np.nan),
])
def test_malformed_or_unnormalized_psf_fails_closed(bad):
    result=measure_finite_native_psf_core(bad,guessed_fwhm_pix=6.57)
    assert result["status"]=="INVALID_FINITE_NORMALIZED_KERNEL"
    assert result["no_psf_model_change_or_science_adoption"]


def test_nonbracketed_half_maximum_never_fakes_measured_fwhm():
    flat=np.ones((81,81))
    flat[40,40]=1.1  # center is the max but edges stay above half peak
    flat/=flat.sum()
    result=measure_finite_native_psf_core(flat,guessed_fwhm_pix=6.57)
    assert result["status"]=="UNBRACKETED_CORE_HALF_MAXIMUM"
    assert "measured_core_fwhm_axis_geometric_mean_native_pix" not in result


def test_offcenter_normalized_psf_keeps_measured_offset():
    source=gaussian_kernel(size=81)
    moved=np.roll(source,3,axis=1)
    result=measure_finite_native_psf_core(
        moved,guessed_fwhm_pix=6.57)
    assert result["status"]=="FINITE_NATIVE_CORE_GEOMETRY_MEASURED"
    assert result["peak_offset_from_native_grid_center_pix"]==3.
    assert result["grid_center_native_pixel_yx"]==[40,40]
    assert result["peak_native_pixel_yx"]==[40,43]


def test_evidence_passes_do_not_promote_to_science():
    geometry=measure_finite_native_psf_core(
        gaussian_kernel(),guessed_fwhm_pix=6.57)
    out=evaluate_core_wings_readiness(
        core_geometry=geometry,n_train=10,n_holdout=3,
        coverage=complete_coverage(),holdout_growth_summary=growth(),
        preclip_negative_mass_fraction=.016537987,
        verified_native_source_hashes=True)
    assert out["status"]=="SCIENCE_GATE_HOLD"
    assert out["science_adopted"] is False
    for name in (
        "original_native_input_provenance",
        "finite_kernel_numerical_geometry",
        "predeclared_holdout_coverage"):
        assert out["stages"][name]["status"]=="PASS"
    for name in (
        "core_accuracy",
        "wide_angle_wing_measurement",
        "sky_versus_psf_wing_identifiability",
        "spatial_representativity",
        "resampled_coadd_noise_covariance",
        "galaxy_injection_recovery"):
        assert out["stages"][name]["status"]=="HOLD"
    assert out["preclipping_negative_absolute_mass_fraction"]==(
        pytest.approx(.016537987))
    assert out["holdout_growth_summary_existing_descriptive_only"][
        "radii"][0]["median_absolute_growth_excess_over_core"]==(
            pytest.approx(.095869648))
    assert out["no_threshold_optimized_on_original_holdouts"]
    assert out["no_physical_wing_flux_inference_from_finite_kernel"]
    assert out["no_science_adoption_without_photometric_bias_calibration"]


def test_failed_input_or_incomplete_holdout_cannot_pass_gate():
    geometry=measure_finite_native_psf_core(
        gaussian_kernel(),guessed_fwhm_pix=6.57)
    out=evaluate_core_wings_readiness(
        core_geometry=geometry,n_train=10,n_holdout=3,
        coverage={
            "all_fixed_holdouts_and_annuli_scored":False,
            "n_expected_holdout_stars":3,
            "n_successfully_scored_stars":2,
        },holdout_growth_summary={},
        preclip_negative_mass_fraction=None,
        verified_native_source_hashes=False)
    assert out["stages"]["original_native_input_provenance"]["status"]=="HOLD"
    assert out["stages"]["predeclared_holdout_coverage"]["status"]=="HOLD"
    assert out["preclipping_negative_absolute_mass_fraction"] is None
    assert out["status"]=="SCIENCE_GATE_HOLD"
    assert out["science_adopted"] is False



def test_psf_core_wings_gui_renders_measured_fwhm_and_closed_science_gates():
    pytest.importorskip("PySide6")
    from lsbpipe.gui.psf0_panel import PSF0Panel
    core = measure_finite_native_psf_core(
        gaussian_kernel(), guessed_fwhm_pix=6.57)
    readiness = evaluate_core_wings_readiness(
        core_geometry=core, n_train=10, n_holdout=3,
        coverage=complete_coverage(), holdout_growth_summary=growth(),
        preclip_negative_mass_fraction=.016537987,
        verified_native_source_hashes=True)
    summary = PSF0Panel._calibration_diagnostic_summary({
        "core_wings_validation_qa":readiness,
    })
    assert "measured native-kernel" in summary
    assert "SCIENCE_GATE_HOLD" in summary
    assert "finite_kernel_numerical_geometry: PASS" in summary
    assert "predeclared_holdout_coverage: PASS" in summary
    assert "wide_angle_wing_measurement: HOLD" in summary
    assert "galaxy_injection_recovery: HOLD" in summary
    assert "PASS is diagnostic completeness only" in summary
    assert "LSB injection-recovery are unvalidated" in summary


@pytest.mark.parametrize("guess", (0., -1., np.nan, np.inf, True, "4.8", None))
def test_invalid_finder_fwhm_prior_cannot_produce_measured_psf_core(guess):
    result = measure_finite_native_psf_core(
        gaussian_kernel(), guessed_fwhm_pix=guess)
    assert result["status"] == "INVALID_DETECTION_FWHM_PRIOR"
    assert result["guessed_detection_fwhm_pix_not_a_measurement"] is None
    assert result["no_psf_model_change_or_science_adoption"] is True
