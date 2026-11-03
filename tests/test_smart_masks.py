"""Synthetic, no-network acceptance for researcher-guided star mask proposals."""
import numpy as np
import pytest

pytest.importorskip("astropy")
from astropy.wcs import WCS

from lsbpipe.photometry.contracts import NativeBandFrame
from lsbpipe.photometry.smart_masks import (
    SMART_MASK_STATUS, StarMaskSeed, propose_seeded_star_mask)
from lsbpipe.scene import LayerKind, SceneState


def _wcs(*, mirrored=False):
    w = WCS(naxis=2)
    w.wcs.crval = [359.999, 18.]
    w.wcs.crpix = [61., 61.]
    w.wcs.cd = np.array([[1 if mirrored else -1, 0], [0, 1]]) / 3600.
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    return w


def _star(shape, center, sigma=2.5, flux=80.):
    y, x = np.indices(shape)
    return flux * np.exp(-.5 * ((x-center[0])**2+(y-center[1])**2) / sigma**2)


def _frame(wcs, *, band="r", center_xy=(60., 60.), sigma=2.5, ivar=None,
           flux=80.):
    shape = (121, 121)
    original = 4. + _star(shape, center_xy, sigma, flux)
    return NativeBandFrame(
        image_nMgy=original,
        invvar_nMgy2=np.ones(shape) if ivar is None else ivar,
        nexp=np.ones(shape), maskbits=np.zeros(shape, np.int16),
        native_wcs=wcs, release="dr9", region="south",
        brick="0401p012", band=band,
        source_evidence={"synthetic": True})


def _seed(wcs, at_xy=(60., 60.), radius=12.):
    ra, dec = wcs.all_pix2world(at_xy[0], at_xy[1], 0)
    return StarMaskSeed("foreground-star-1", float(ra) % 360,
                        float(dec), radius)


def test_single_click_proposes_connected_star_contour_but_never_adopts():
    w = _wcs()
    frame = _frame(w)
    originals = frame.image_nMgy.tobytes()
    result = propose_seeded_star_mask(
        frame, _seed(w), residual_sky_nMgy_per_pixel=4.,
        local_rms_nMgy_per_pixel=1., smoothing_sigma_pix=.7)
    assert result.status == SMART_MASK_STATUS
    assert result.dataset_id == "dr9/south/0401p012/r"
    assert result.mask[60, 60] and not result.mask[0, 0]
    assert 20 < result.n_pixels < 300
    assert "EXTENDED_PSF_WINGS_NOT_VALIDATED" in result.flags
    assert not result.summary()["mask_is_validated_science"]
    assert not result.mask.flags.writeable
    with pytest.raises(ValueError):
        result.mask.setflags(write=True)
    assert frame.image_nMgy.tobytes() == originals
    with pytest.raises(ValueError, match="explicitly confirm"):
        result.to_reviewed_scene_layer(name="masked_star", confirmed=False)
    layer = result.to_reviewed_scene_layer(name="masked_star", confirmed=True)
    assert layer.kind == LayerKind.MANUAL_EXCLUDE
    assert layer.band == "r"
    state = SceneState(target_id="dwarf", grid_id=result.grid_id,
                       band="r", shape=frame.shape)
    updated = state.with_layer(layer)
    assert updated.layers[0].name == "masked_star"
    assert state.layers == ()


def test_one_icrs_click_reprojects_and_refits_separately_across_bands():
    w_r, w_g = _wcs(), _wcs(mirrored=True)
    # World-space click 10 original pixels east of frame center.
    anchor = _seed(w_r, at_xy=(50., 60.), radius=15.)
    gx, gy = w_g.all_world2pix(anchor.ra_deg, anchor.dec_deg, 0)
    r = _frame(w_r, center_xy=(50., 60.), sigma=2., flux=80.)
    g = _frame(w_g, band="g", center_xy=(gx, gy), sigma=4., flux=80.)
    r_proposal = propose_seeded_star_mask(
        r, anchor, residual_sky_nMgy_per_pixel=4.,
        local_rms_nMgy_per_pixel=1.)
    g_proposal = propose_seeded_star_mask(
        g, anchor, residual_sky_nMgy_per_pixel=4.,
        local_rms_nMgy_per_pixel=1.)
    assert r_proposal.band == "r" and g_proposal.band == "g"
    assert abs(r_proposal.peak_x_pix - 50) <= 1
    assert abs(g_proposal.peak_x_pix - gx) <= 1
    assert abs(g_proposal.peak_y_pix - gy) <= 1
    assert r_proposal.peak_x_pix != g_proposal.peak_x_pix
    assert g_proposal.n_pixels > r_proposal.n_pixels
    assert r_proposal.grid_id != g_proposal.grid_id


def test_bad_native_pixel_cannot_be_invented_by_candidate():
    w = _wcs()
    ivar = np.ones((121, 121))
    ivar[60, 60] = 0.
    f = _frame(w, ivar=ivar)
    proposal = propose_seeded_star_mask(
        f, _seed(w), residual_sky_nMgy_per_pixel=4.,
        local_rms_nMgy_per_pixel=1., dilation_pixels=2)
    assert not f.valid[60, 60]
    assert not proposal.mask[60, 60]
    assert "INVALID_NATIVE_PIXELS_NEAR_SOURCE" in proposal.flags


def test_target_overlap_must_be_reviewed_and_never_automatically_blocked():
    w = _wcs()
    target = np.zeros((121, 121), dtype=bool)
    target[55:65, 55:65] = True
    candidate = propose_seeded_star_mask(
        _frame(w), _seed(w), residual_sky_nMgy_per_pixel=4.,
        local_rms_nMgy_per_pixel=1., target_protection=target)
    assert candidate.mask[60, 60]
    assert "TARGET_PROTECTION_OVERLAP_REVIEW_REQUIRED" in candidate.flags


def test_galaxy_filling_search_region_is_refused():
    w = _wcs()
    f = _frame(w, sigma=30., flux=50.)
    with pytest.raises(ValueError, match="fills most"):
        propose_seeded_star_mask(
            f, _seed(w), residual_sky_nMgy_per_pixel=4.,
            local_rms_nMgy_per_pixel=1.)


def test_low_contrast_false_seed_does_not_produce_mask():
    w = _wcs()
    f = _frame(w, flux=0.)
    with pytest.raises(ValueError, match="No connected source"):
        propose_seeded_star_mask(
            f, _seed(w), residual_sky_nMgy_per_pixel=4.,
            local_rms_nMgy_per_pixel=1.)
    with pytest.raises(ValueError, match="RMS"):
        propose_seeded_star_mask(
            f, _seed(w), residual_sky_nMgy_per_pixel=4.,
            local_rms_nMgy_per_pixel=0.)


def test_seed_and_parameters_fail_closed():
    w = _wcs()
    f = _frame(w)
    with pytest.raises(ValueError):
        StarMaskSeed("star", 12., 91., 10.)
    with pytest.raises(ValueError):
        StarMaskSeed("star", 12., 18., -5.)
    with pytest.raises(ValueError, match="outside"):
        propose_seeded_star_mask(
            f, StarMaskSeed("elsewhere", 122., -20., 10.),
            residual_sky_nMgy_per_pixel=4.,
            local_rms_nMgy_per_pixel=1.)
    with pytest.raises(ValueError, match="Target protection"):
        propose_seeded_star_mask(
            f, _seed(w), residual_sky_nMgy_per_pixel=4.,
            local_rms_nMgy_per_pixel=1.,
            target_protection=np.ones((20, 20), dtype=bool))
