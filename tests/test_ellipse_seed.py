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
    with pytest.raises(ValueError, match="finite and positive"):
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


def test_local_isophote_returns_connected_contour_and_equivalent_ellipse():
    from lsbpipe.photometry.ellipse_seed import propose_local_isophote
    rng = np.random.default_rng(65)
    ny = nx = 241
    yy, xx = np.indices((ny, nx))
    ang = np.deg2rad(34.)
    dx, dy = xx-120., yy-120.
    u = np.cos(ang)*dx + np.sin(ang)*dy
    v = -np.sin(ang)*dx + np.cos(ang)*dy
    image = (3.7 + rng.normal(0.,1.,(ny,nx))
             + 2.2*np.exp(-.5*((u/28.)**2+(v/14.)**2)))
    original = image.copy()
    # A much brighter unrelated object away from ICRS target must not
    # attract the connected region's centre.
    image[19:27, 210:218] += 1000.
    invalid = np.zeros_like(image, bool)
    invalid[125,125] = True
    output = propose_local_isophote(
        image, ~invalid, (120.,120.),
        threshold_contrast=2.2, smooth_sigma_pix=4.)
    assert output.state == "LOCAL_ISOPHOTE_PROPOSAL_NOT_ADOPTED"
    assert output.seed.x_pix == 120. and output.seed.y_pix == 120.
    assert output.seed.component_pixels > 100
    assert .30 < output.seed.q < .8
    assert 0. < output.raw_sky_rms < 2.
    assert output.smoothed_sky_rms < output.raw_sky_rms
    assert len(output.contour_xy) <= 192
    assert len(output.contour_xy) >= 24
    x, y = np.asarray(output.contour_xy).T
    assert np.isfinite(x).all() and np.isfinite(y).all()
    assert max(abs(x-120.)) < 110. and max(abs(y-120.)) < 110.
    assert np.array_equal(original[:18,:18], image[:18,:18])
    assert "BACKGROUND_CONSTANT_NOT_VALIDATED" in output.flags


def test_local_isophote_fails_closed_on_blank_or_unanchored_light():
    from lsbpipe.photometry.ellipse_seed import propose_local_isophote
    blank = np.zeros((241,241))
    with pytest.raises(ValueError, match="background dispersion"):
        propose_local_isophote(
            blank, np.ones_like(blank, bool), (120.,120.))
    rng = np.random.default_rng(19)
    image = rng.normal(0.,1.,(241,241))
    yy, xx = np.indices(image.shape)
    image += 50*np.exp(-.5*(((xx-208.)/6.)**2+((yy-40.)/6.)**2))
    with pytest.raises(ValueError, match="No connected"):
        propose_local_isophote(
            image, np.ones_like(image,bool), (120.,120.),
            threshold_contrast=7., smooth_sigma_pix=3.)
    with pytest.raises(ValueError, match="outside"):
        propose_local_isophote(
            image, np.ones_like(image,bool), (500.,120.))
