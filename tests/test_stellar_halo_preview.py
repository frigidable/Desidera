import numpy as np
import pytest
from lsbpipe.photometry.stellar_halo_preview import conditional_halo_images


def run(a,mask=None):
    zero=np.zeros(a.shape,dtype=bool)
    c=dict(training_r_range_arcsec=[12.,80.],radial_knots_arcsec=[12.,30.,70.],
        halo_at_knots_nMgy_pix=[.02,.005,0.],sky_constant_nMgy_pix=-.001,
        sky_x_nMgy_pix_per_arcsec=0.,sky_y_nMgy_pix_per_arcsec=0.,
        side_median_abs_residual_nMgy_pix=.001,side_radial_residuals=[
        dict(sector=2,r_inner_arcsec=20.,n=4,observed_minus_model_median_nMgy_pix=.002),
        dict(sector=6,r_inner_arcsec=20.,n=4,observed_minus_model_median_nMgy_pix=-.002)])
    d=dict(c,halo_at_knots_nMgy_pix=[.025,.006,0.],side_median_abs_residual_nMgy_pix=.002)
    return conditional_halo_images(a,np.ones_like(a),zero,zero if mask is None else mask,
        source_xy=(100.,50.),target_xy=(100.,100.),pixel_scale_arcsec=.5,
        comparison={"candidates":[c,d]},radial_edges=[12.,20.,40.,80.])


def test_signed_apertures_support_and_immutable_input():
    a=np.full((201,201),-.001); before=a.copy()
    mask=np.zeros(a.shape,dtype=bool); mask[100,101]=True
    p=run(a,mask)
    assert np.array_equal(a,before)
    assert not p["science_adopted"]
    assert np.isnan(p["layers"]["residual"][50,100])
    assert np.isnan(p["layers"]["residual"][100,101])
    assert not p["layers"]["residual"].flags.writeable
    for row in p["apertures"]:
        assert row["raw_matched_flux_nMgy"]<0
        for m in row["model_fluxes"]:
            assert m["halo_corrected_flux_nMgy"]==pytest.approx(row["raw_matched_flux_nMgy"]-m["halo_flux_nMgy"])
            assert m["halo_sky_residual_flux_nMgy"]==pytest.approx(m["halo_corrected_flux_nMgy"]-m["sky_flux_nMgy"])
            assert m["stress_supported_pixels"]<=row["common_support_pixels"]
    assert p["apertures"][-1]["unsupported_valid_pixels"]>0
    assert not p["apertures"][-1]["full_valid_pixel_model_support"]


def test_fixed_model_injection_algebra_only():
    a=np.zeros((201,201)); b=a.copy(); b[95:106,95:106]=.003
    p=run(a); q=run(b)
    good=np.isfinite(p["layers"]["residual"])
    np.testing.assert_allclose((q["layers"]["residual"]-p["layers"]["residual"])[good],b[good],atol=1e-15)
    assert p["apertures"][0]["model_fluxes"][0]["stress_supported_pixels"]>0
