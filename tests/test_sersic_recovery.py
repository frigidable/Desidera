import numpy as np
import pytest
pytest.importorskip("astropy")
from astropy.wcs import WCS
from scipy.special import gammainc,gammaincinv
from lsbpipe.photometry.sersic_recovery import (
    total_flux, radius_pixels, convolved_profile, strict_convolve, fit_profile, coverage_sites)
from lsbpipe.photometry.galaxy_geometry_diagnostic import local_grid


def wcs(n):
    w=WCS(naxis=2); w.wcs.crpix=[(n+1)/2]*2
    w.wcs.crval=[96.,78.]; w.wcs.cdelt=[-.5/3600,.5/3600]
    w.wcs.ctype=["RA---TAN","DEC--TAN"]
    return w


def test_sersic_total_and_half_light():
    assert gammainc(2*.8,gammaincinv(1.6,.5))==pytest.approx(.5)
    assert total_flux(.1,15.,1.,.6)==pytest.approx(2*np.pi*.6*.1*(15/1.67834699001666)**2)


def test_strict_cross_convolution_masks_and_flux():
    a=np.ones((61,61)); m=np.ones_like(a,bool); m[30,30]=False
    k=np.array([[0,1,0],[1,4,1],[0,1,0]],float)/8
    out,valid=strict_convolve(a,m,k)
    assert not valid[29:32,29:32].any() and not valid[0].any()
    assert np.allclose(out[valid],1.)
    assert np.isnan(out[~valid]).all()


def test_forward_profile_recovers_signed_observed_ring_fit():
    n=201; w=wcs(n); raw=np.zeros((n,n))
    grid=local_grid(dict(target_xy_native_pix=[100.,100.],layers={"original":raw}),w)
    geom=dict(east_offset_arcsec=0.,north_offset_arcsec=0.,q=.56,pa_deg=44.5)
    k=np.ones((5,5))/25
    rho=radius_pixels(grid["x"].shape,(100.,100.),grid["k"],.56,44.5,
        origin=(grid["x"][0,0],grid["y"][0,0]),pad=(2,2))
    image=.12*grid["area"]*convolved_profile(rho,27.,.8,k)
    valid=np.ones_like(image,bool);valid[:,95:100]=False
    f=fit_profile(image,valid,grid,geom,k)
    assert f["status"]=="SERSIC_CONDITIONAL_NOT_ADOPTED"
    assert f["re_arcsec"]==pytest.approx(27.,rel=.02)
    assert f["n"]==pytest.approx(.8,rel=.03)
    assert f["analytic_total_nMgy"]==pytest.approx(total_flux(.12,27.,.8,.56),rel=.03)
    assert f["apertures"][-1]["observed_signed_flux_nMgy"]==pytest.approx(
        f["apertures"][-1]["model_on_observed_pixels_nMgy"],rel=.02)
    blocked=fit_profile(image,np.zeros_like(valid),grid,geom,k)
    assert blocked["status"]=="SERSIC_INSUFFICIENT_RADIAL_SUPPORT"


def test_plateau_screen_keeps_signed_growth_and_never_adopts():
    from lsbpipe.photometry.sersic_recovery import growth_saturation
    def rows(flux):
        return [dict(sector="ALL",outer_arcsec=a,scenarios={"r_plane":{"cumulative_flux_nMgy":f}})
                for a,f in zip((30.,35.,40.),flux)]
    assert growth_saturation(rows([216.,254.,288.]),"r_plane")["status"]=="COG_NO_EMPIRICAL_PLATEAU_AT_OUTER_RADIUS"
    p=growth_saturation(rows([100.,99.,100.]),"r_plane")
    assert p["status"]=="COG_SMALL_LAST_INCREMENTS_UNCALIBRATED"
    assert p["last_two_signed_increments_over_outer_flux"]==[-.01,.01]
    assert not p["science_adopted"]
