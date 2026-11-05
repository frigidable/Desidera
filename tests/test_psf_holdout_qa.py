"""Independent holdout signed-PSF residual diagnostics are not PSF adoption."""
import numpy as np
import pytest

from lsbpipe.photometry.psf_holdout_qa import (
    heldout_signed_residual, heldout_final_kernel_residual,
    summarize_signed_holdouts, audit_signed_holdout_coverage,
)


def _source():
    yy, xx = np.indices((31,31))
    kernel = np.exp(-.5*((xx-15)**2+(yy-15)**2)/2.**2)
    kernel[15,26] = -.005
    return kernel


def test_signed_model_holdout_uses_core_amplitude_and_does_not_clip_negative():
    kernel = _source()
    original = kernel.tobytes()
    qa = heldout_signed_residual(
        source_id="holdout-1", image=85.*kernel,
        sqrt_ivar=np.ones(kernel.shape), signed_kernel=kernel,
        local_center_xy=(15.,15.), guessed_fwhm_pix=4.8,
        measured_peak=85.)
    assert qa["source_id"] == "holdout-1"
    assert qa["fitted_amplitude_model_units"] == pytest.approx(85.,rel=1e-7)
    assert qa["signed_kernel_unchanged"] is True
    assert qa["correlated_coadd_noise_not_whitened"] is True
    for region in qa["radial_residuals"].values():
        assert region["median_absolute_residual_over_peak"] < 1e-7
    assert kernel.tobytes() == original


def test_outer_residual_is_reported_by_independent_star_sample():
    kernel = _source()
    yy, xx = np.indices(kernel.shape)
    r = np.hypot(xx-15, yy-15)
    measured = 90.*kernel + np.where(r>9.6, .05, 0.)
    qa = heldout_signed_residual(
        source_id="holdout-2", image=measured,
        sqrt_ivar=np.ones(kernel.shape), signed_kernel=kernel,
        local_center_xy=(15.,15.), guessed_fwhm_pix=4.8,
        measured_peak=90.)
    outer = qa["radial_residuals"]["outer_r_gt_2x_guessed_fwhm"]
    inner = qa["radial_residuals"]["core_r_le_1x_guessed_fwhm"]
    assert outer["median_signed_residual_over_peak"] > 0.
    assert outer["median_absolute_residual_over_peak"] > inner[
        "median_absolute_residual_over_peak"]
    another = heldout_signed_residual(
        source_id="holdout-3", image=95.*kernel + np.where(r>9.6, .07, 0.),
        sqrt_ivar=np.ones(kernel.shape), signed_kernel=kernel,
        local_center_xy=(15., 15.), guessed_fwhm_pix=4.8,
        measured_peak=95.)
    summary=summarize_signed_holdouts([qa, another])
    assert summary["n_usable_holdout_stars"] == 2
    assert summary["regions"]["outer_r_gt_2x_guessed_fwhm"]["n_stars"] == 2
    assert summary["regions"]["outer_r_gt_2x_guessed_fwhm"][
        "median_star_abs_residual_over_peak"] > 0.


def test_empty_or_bad_holdout_does_not_fabricate_qa():
    empty = summarize_signed_holdouts([])
    assert empty["n_usable_holdout_stars"] == 0
    assert empty["regions"]["outer_r_gt_2x_guessed_fwhm"][
        "median_star_abs_residual_over_peak"] is None
    kernel = _source()
    with pytest.raises(ValueError, match="Too few"):
        heldout_signed_residual(
            source_id="bad", image=kernel,
            sqrt_ivar=np.zeros(kernel.shape), signed_kernel=kernel,
            local_center_xy=(15.,15.), guessed_fwhm_pix=4.8,
            measured_peak=1.)


def test_repeated_source_id_cannot_create_fake_independent_holdout_evidence():
    kernel = _source()
    qa = heldout_signed_residual(
        source_id="gaia-physical-star-42", image=90.*kernel,
        sqrt_ivar=np.ones(kernel.shape), signed_kernel=kernel,
        local_center_xy=(15., 15.), guessed_fwhm_pix=4.8,
        measured_peak=90.)
    with pytest.raises(ValueError, match="Duplicate holdout source_id"):
        summarize_signed_holdouts([qa, qa])
    with pytest.raises(ValueError, match="Duplicate holdout source_id"):
        summarize_signed_holdouts([qa, {**qa, "source_id": " gaia-physical-star-42 "}])
    with pytest.raises(ValueError, match="nonempty source_id"):
        summarize_signed_holdouts([{**qa, "source_id": "  "}])
    with pytest.raises(TypeError, match="mapping"):
        summarize_signed_holdouts([None])

def test_final_kernel_holdout_uses_identical_radii_after_clipping():
    signed = _source()
    yy, xx = np.indices(signed.shape)
    # A scientifically unphysical signed wing influences the preclip
    # residual. Its delivered model must be checked separately.
    image = 85.*signed
    final = np.maximum(signed, 0.)
    final /= final.sum()
    orig = final.tobytes()
    kwargs = dict(source_id="holdout-fixed", image=image,
                  sqrt_ivar=np.ones(signed.shape),
                  local_center_xy=(15., 15.), guessed_fwhm_pix=4.8,
                  measured_peak=85.)
    signed_qa = heldout_signed_residual(signed_kernel=signed, **kwargs)
    final_qa = heldout_final_kernel_residual(final_kernel=final, **kwargs)
    assert final_qa["kernel_representation"] == "FINAL_NONNEGATIVE_NORMALIZED"
    assert final_qa["final_kernel_nonnegative_and_unit_sum"] is True
    assert "signed_kernel_unchanged" not in final_qa
    assert final.tobytes() == orig
    assert set(signed_qa["radial_residuals"]) == set(final_qa["radial_residuals"])
    for zone in final_qa["radial_residuals"]:
        assert (final_qa["radial_residuals"][zone]["n_valid_pixels"]
                == signed_qa["radial_residuals"][zone]["n_valid_pixels"])
    # Pre-clipping signed residual is zero by construction; the negative
    # wing is absent from the delivered positive PSF.
    assert signed_qa["radial_residuals"]["outer_r_gt_2x_guessed_fwhm"][
        "median_absolute_residual_over_peak"] == pytest.approx(0., abs=1e-12)
    assert final_qa["radial_residuals"]["outer_r_gt_2x_guessed_fwhm"][
        "median_absolute_residual_over_peak"] >= 0.
    coverage = audit_signed_holdout_coverage(["holdout-fixed"],
                                              [final_qa], [])
    assert coverage["status"] == "COMPLETE"


@pytest.mark.parametrize("bad", ["signed", "unnormalized", "nan"])
def test_final_holdout_rejects_nonphysical_kernel(bad):
    kernel = _source()
    if bad == "unnormalized":
        kernel = np.maximum(kernel, 0.)
    elif bad == "nan":
        kernel[0,0] = np.nan
    with pytest.raises(ValueError, match="nonnegative and normalized"):
        heldout_final_kernel_residual(
            source_id="holdout-1", image=np.ones((31,31)),
            sqrt_ivar=np.ones((31,31)), final_kernel=kernel,
            local_center_xy=(15.,15.), guessed_fwhm_pix=4.8,
            measured_peak=1.)
