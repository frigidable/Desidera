"""Independent holdout signed-PSF residual diagnostics are not PSF adoption."""
import numpy as np
import pytest

from lsbpipe.photometry.psf_holdout_qa import (
    heldout_signed_residual, summarize_signed_holdouts,
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
    summary=summarize_signed_holdouts([qa,qa])
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
