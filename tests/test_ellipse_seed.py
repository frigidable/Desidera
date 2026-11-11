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


def test_multiscale_sensitivity_reports_nearby_thresholds_and_low_confidence():
    from dataclasses import replace
    from lsbpipe.photometry.ellipse_seed import (
        propose_local_isophote, assess_local_isophote_sensitivity,
    )
    rng=np.random.default_rng(106)
    yy,xx=np.indices((241,241))
    image=(2.+rng.normal(0.,.4,(241,241))
           + 4.*np.exp(-.5*(((xx-120.)/24.)**2+
                            ((yy-120.)/13.)**2)))
    good=np.ones_like(image,bool)
    model=propose_local_isophote(
        image,good,(120.,120.),threshold_contrast=2.2,
        smooth_sigma_pix=3.,bridge_guard_pix=2)
    assert "THIN_BRIDGE_OPENING_2PX" in model.flags
    qa=assess_local_isophote_sensitivity(
        image,good,(120.,120.),model,
        threshold_contrast=2.2,smooth_sigma_pix=3.,bridge_guard_pix=2)
    assert qa["trial_count"]==4
    assert qa["sensitivity_status"]=="MULTISCALE_DIAGNOSTIC_NOT_SCIENCE_VALIDATED"
    assert all("state" in trial for trial in qa["trials"])
    assert qa["minimum_spatial_iou"] is not None
    assert 0. <= qa["minimum_spatial_iou"] <= 1.
    # A reference mask in the wrong spatial location is not allowed to
    # pass merely because its area and axis ratio look plausible.
    offset_pixels = np.array(model.component_indices + image.size,copy=True)
    offset_pixels.setflags(write=False)
    shifted=replace(model,component_indices=offset_pixels)
    spatial=assess_local_isophote_sensitivity(
        image,good,(120.,120.),shifted,
        threshold_contrast=2.2,smooth_sigma_pix=3.,bridge_guard_pix=2)
    assert "LOW_SPATIAL_OVERLAP_BETWEEN_ISOPHOTES" in spatial["flags"]
    assert spatial["confidence"]=="LOW_CONFIDENCE_REQUIRES_EXPLICIT_OVERRIDE"
    # Force a known displaced reference to confirm the safety gate does
    # not silently accept large offsets even when other trials succeed.
    displaced=replace(model,seed=replace(
        model.seed,centroid_offset_pix=65.))
    warning=assess_local_isophote_sensitivity(
        image,good,(120.,120.),displaced,
        threshold_contrast=2.2,smooth_sigma_pix=3.,bridge_guard_pix=2)
    assert "UNANCHORED_CENTROID_DRIFT" in warning["flags"]
    assert warning["confidence"]=="LOW_CONFIDENCE_REQUIRES_EXPLICIT_OVERRIDE"


def test_bridge_guard_is_bounded_and_does_not_mutate_native_pixels():
    from lsbpipe.photometry.ellipse_seed import propose_local_isophote
    rng=np.random.default_rng(114)
    yy,xx=np.indices((241,241))
    image=(3.+rng.normal(0.,.3,(241,241))
           +2.*np.exp(-.5*(((xx-120.)/20.)**2+
                            ((yy-120.)/13.)**2)))
    raw=image.copy()
    good=np.ones_like(image,bool)
    model=propose_local_isophote(
        image,good,(120.,120.),bridge_guard_pix=2)
    assert model.seed.component_pixels>100
    assert np.array_equal(raw,image)
    with pytest.raises(ValueError,match="bridge settings"):
        propose_local_isophote(
            image,good,(120.,120.),bridge_guard_pix=9)


def test_sga_guided_radial_extent_follows_native_light_not_parent_diameter():
    from astropy.wcs import WCS
    from lsbpipe.photometry.sky_masks import SkyMaskShape
    from lsbpipe.photometry.ellipse_seed import sga_informed_radial_seed
    w=WCS(naxis=2)
    w.wcs.crpix=[201.,201.]
    w.wcs.crval=[155.,18.]
    w.wcs.cd=np.diag([-1.,1.])/3600.
    w.wcs.ctype=["RA---TAN","DEC--TAN"]
    y,x=np.indices((401,401))
    east=-(x-200.)*np.cos(np.deg2rad(18.))
    north=y-200.
    angle=np.deg2rad(150.)
    major=east*np.sin(angle)+north*np.cos(angle)
    minor=east*np.cos(angle)-north*np.sin(angle)
    a=35.
    q=.55
    gaussian=3.*np.exp(-.5*((major/a)**2+(minor/(a*q))**2))
    noise=np.random.default_rng(115).normal(0.,.08,(401,401))
    image=.3+noise+gaussian
    original=image.copy()
    valid=np.ones_like(image,bool)
    valid[180:206,244:264]=False
    prior=SkyMaskShape("ellipse",center=(155.,18.),
                       a_arcsec=14.,q=q,pa_deg=150.)
    result=sga_informed_radial_seed(
        image,valid,w,prior,
        local_sky=.3,smoothed_rms=.05,smooth_sigma_pix=3.)
    assert result.status=="SGA_INFORMED_RADIAL_SEED_NOT_ADOPTED"
    assert result.reference_a_arcsec==pytest.approx(14.)
    assert result.detected_a_arcsec is not None
    assert 26. < result.detected_a_arcsec < 135.
    assert result.radial_step_arcsec is not None
    assert result.radial_step_arcsec > 0.
    # The reported radial endpoint is the outer edge of an integer
    # annulus, not a continuously resolved isophotal semi-major axis.
    ring_number=result.detected_a_arcsec/result.radial_step_arcsec
    assert ring_number==pytest.approx(round(ring_number),abs=1e-9)
    assert result.q_sky==pytest.approx(q)
    assert result.pa_deg==pytest.approx(150.)
    assert 0 < result.masked_fraction < .30
    assert len(result.bins)>10
    # New diagnostic octants are position-preserving: missing sectors
    # remain null rather than shifting the identity of later sectors.
    for ring in result.bins:
        octants=ring["sector_contrast_by_octant"]
        assert len(octants)==8
        assert sum(v is not None for v in octants)==ring["sectors_with_data"]
        assert all(np.isfinite(v) for v in octants if v is not None)
        assert sum(v >= .4*.65 for v in octants if v is not None)==(
            ring["sectors_positive"])
    assert np.array_equal(image,original)


def test_sga_guided_radial_seed_does_not_invent_extent_on_sky_noise():
    from astropy.wcs import WCS
    from lsbpipe.photometry.sky_masks import SkyMaskShape
    from lsbpipe.photometry.ellipse_seed import sga_informed_radial_seed
    w=WCS(naxis=2)
    w.wcs.crpix=[201.,201.]
    w.wcs.crval=[155.,18.]
    w.wcs.cd=np.diag([-1.,1.])/3600.
    w.wcs.ctype=["RA---TAN","DEC--TAN"]
    image=np.random.default_rng(214).normal(0.,.1,(401,401))
    prior=SkyMaskShape("ellipse",center=(155.,18.),
                       a_arcsec=14.,q=.53,pa_deg=151.)
    result=sga_informed_radial_seed(
        image,np.ones_like(image,bool),w,prior,
        local_sky=0.,smoothed_rms=.1)
    assert result.detected_a_arcsec is None
    assert "NO_DETECTED_OUTER_RADIUS" in result.flags
    assert all(not b["supported"] for b in result.bins[:5])
    assert result.status=="SGA_INFORMED_RADIAL_SEED_NOT_ADOPTED"
    with pytest.raises(ValueError,match="valid selected SGA ellipse"):
        sga_informed_radial_seed(image,np.ones_like(image,bool),w,None,
                                 local_sky=0.,smoothed_rms=.1)


def test_sga_radial_matches_target_protection_radius_flagged_not_confirmed():
    """A coincident 60-arcsec screening circle cannot validate Leo P extent."""
    from astropy.wcs import WCS
    from lsbpipe.photometry.sky_masks import SkyMaskShape
    from lsbpipe.photometry.ellipse_seed import sga_informed_radial_seed
    w = WCS(naxis=2)
    w.wcs.crpix = [201., 201.]
    w.wcs.crval = [155., 18.]
    w.wcs.cd = np.diag([-1., 1.])/3600.
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    y, x = np.indices((401,401))
    q = .55
    theta = np.deg2rad(150.)
    east, north = 200.-x, y-200.
    major = east*np.sin(theta)+north*np.cos(theta)
    minor = east*np.cos(theta)-north*np.sin(theta)
    image = (.3 + 3.*np.exp(-.5*((major/35.)**2+
                                  (minor/(35.*q))**2))
             + np.random.default_rng(115).normal(0.,.08,(401,401)))
    shape = SkyMaskShape("ellipse",center=(155.,18.),
                         a_arcsec=14.,q=q,pa_deg=150.)
    good = np.ones_like(image,bool)
    without = sga_informed_radial_seed(
        image, good, w, shape,local_sky=.3,smoothed_rms=.05,
        smooth_sigma_pix=3.)
    assert without.detected_a_arcsec is not None
    # The coincident radius is an INPUT SCREENING CHOICE, not a
    # photometric measurement, so coincidence is only diagnostic QA.
    near = sga_informed_radial_seed(
        image,good,w,shape,local_sky=.3,smoothed_rms=.05,
        smooth_sigma_pix=3.,
        target_protection_radius_arcsec=without.detected_a_arcsec)
    assert near.detected_a_arcsec == without.detected_a_arcsec
    assert near.target_protection_radius_arcsec == pytest.approx(
        without.detected_a_arcsec)
    assert "RADIAL_EXTENT_NEAR_PROVISIONAL_TARGET_PROTECTION_RADIUS" in near.flags
    assert "MASK_BOUNDARY_SENSITIVITY_NOT_TESTED" in near.flags
    far = sga_informed_radial_seed(
        image,good,w,shape,local_sky=.3,smoothed_rms=.05,
        smooth_sigma_pix=3.,target_protection_radius_arcsec=12.)
    assert "RADIAL_EXTENT_NEAR_PROVISIONAL_TARGET_PROTECTION_RADIUS" not in far.flags
    with pytest.raises(ValueError,match="protection radius invalid"):
        sga_informed_radial_seed(
            image,good,w,shape,local_sky=.3,smoothed_rms=.05,
            target_protection_radius_arcsec=0.)


def test_outer_boundary_keeps_edges_at_diagonal_void_contacts():
    """A valid 4-connected galaxy region may have diagonally touching holes.

    Previously an edge dictionary silently overwrote one outgoing edge
    at the shared pixel-grid vertex and raised "Disconnected isophotal
    pixel boundary", although the outer contour was fully connected.
    """
    from lsbpipe.photometry.ellipse_seed import _outer_boundary

    region = np.zeros((42, 42), dtype=bool)
    region[4:38, 4:38] = True
    region[14, 14] = False
    region[15, 15] = False
    # The source region is still 4-connected; only the two holes touch
    # diagonally. The exterior must enclose the *whole* source, not a hole.
    from scipy.ndimage import label
    components, count = label(region)
    assert count == 1

    outline = np.asarray(_outer_boundary(region, 100, 200))
    assert 24 <= len(outline) <= 192
    assert np.isfinite(outline).all()
    assert np.ptp(outline[:, 0]) > 30.
    assert np.ptp(outline[:, 1]) > 30.
    assert 101. <= outline[:, 0].min() < 108.
    assert 201. <= outline[:, 1].min() < 208.


def test_outer_boundary_diagonal_hole_chain_is_deterministic():
    """Three corner-touching holes do not change the chosen outer ring."""
    from lsbpipe.photometry.ellipse_seed import _outer_boundary

    region = np.zeros((42, 42), dtype=bool)
    region[4:38, 4:38] = True
    for row in (14, 15, 16):
        region[row, row] = False
    a = _outer_boundary(region, 0, 0)
    b = _outer_boundary(region.copy(), 0, 0)
    assert a == b
    assert len(a) >= 24
