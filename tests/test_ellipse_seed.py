"""Preliminary threshold-isophote geometry, not adopted photometry."""
import numpy as np
import pytest

from lsbpipe.photometry.ellipse_seed import (
    SEED_STATUS, propose_isophotal_seed, seed_as_sky_geometry)


def _image(shape=(141, 141), center=(70, 70), a=22., q=.45, theta=32.):
    yy, xx = np.indices(shape)
    th = np.deg2rad(theta)
    dx, dy = xx-center[0], yy-center[1]
    u = np.cos(th)*dx + np.sin(th)*dy
    v = -np.sin(th)*dx + np.cos(th)*dy
    return 20. * np.exp(-.5*((u/a)**2+(v/(a*q))**2))


def test_seed_recovers_ellipse_and_preserves_target_center():
    image = _image()
    good = np.ones_like(image, bool)
    output = propose_isophotal_seed(
        image, .5, good, (70., 70.),
        threshold_rms=3., smooth_sigma_pix=1.5, min_pixels=30)
    assert output.state == SEED_STATUS
    assert output.x_pix == 70. and output.y_pix == 70.
    assert 0.35 < output.q < .60
    assert min(abs(output.theta_native_deg-32),
               abs(output.theta_native_deg-212)) < 10
    assert output.sma_pix > 10.
    assert "SMOOTHED_CONTRAST_NOT_CALIBRATED_SIGNIFICANCE" in output.flags


def test_nearby_bright_masked_contaminant_does_not_move_center():
    image = _image()
    image[20:30, 20:30] = 100000.
    valid = np.ones_like(image, bool)
    valid[20:30, 20:30] = False
    output = propose_isophotal_seed(image, 1., valid, (70., 70.))
    assert 60 < output.sma_pix + 70 < 130
    assert output.x_pix == 70.


def test_fails_gracefully_on_non_detected_target():
    im = np.zeros((51, 51), float)
    with pytest.raises(ValueError, match="No isophotal component"):
        propose_isophotal_seed(im, 1., np.ones_like(im, bool), (25., 25.))
    with pytest.raises(ValueError, match="outside"):
        propose_isophotal_seed(im, 1., np.ones_like(im, bool), (200., 25.))
    with pytest.raises(ValueError, match="positive RMS"):
        propose_isophotal_seed(im, 0., np.ones_like(im, bool), (25., 25.))


def test_local_celestial_pa_conversion_on_inverted_ra_axis():
    from astropy.wcs import WCS
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [71., 71.]
    wcs.wcs.crval = [155., 18.]
    wcs.wcs.cd = np.array([[-1, 0], [0, 1]])/3600.
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    seed = propose_isophotal_seed(
        _image(theta=90., q=.5), 1.,
        np.ones((141, 141), bool), (70., 70.))
    geometry = seed_as_sky_geometry(seed, wcs)
    assert geometry["source"] == "LOCAL_WCS_ISOPHOTE_MOMENT_SEED"
    assert geometry["q_sky"] == pytest.approx(seed.q, rel=.005)
    assert geometry["pa_north_through_east_deg"] == pytest.approx(
        0., abs=10.) or geometry["pa_north_through_east_deg"] == pytest.approx(180., abs=10.)
    assert geometry["sma_arcsec"] > 10.
