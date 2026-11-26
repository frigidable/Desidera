import numpy as np
import pytest
from astropy.wcs import WCS
from lsbpipe.photometry.galaxy_geometry_diagnostic import (
    local_grid,block_design,fit_geometry,radial_profiles,geometry_diagnostics,paired_geometry)


def fixture():
    w=WCS(naxis=2); w.wcs.crpix=[181.,181.]; w.wcs.crval=[96.,78.]
    w.wcs.cdelt=[-.262/3600,.262/3600]; w.wcs.ctype=["RA---TAN","DEC--TAN"]
    a=np.zeros((361,361))
    p=dict(layers={"original":a,"residual":a},target_xy_native_pix=[180.,180.],
           source_xy_native_pix=[188.,-35.],pixel_scale_arcsec=.262)
    g=local_grid(p,w)
    e=g["e"]-1.5;n=g["n"]+2.
    t=np.deg2rad(35.)
    x=e*np.sin(t)+n*np.cos(t); y=e*np.cos(t)-n*np.sin(t)
    a[:]=.012*np.exp(-np.log(2)*np.hypot(x,y/.62)/8.)-.0008
    return w,p,g,a


def test_robust_geometry_with_mask_and_signed_background():
    w,p,g,a=fixture()
    v=a.copy(); v[:,174:187]=np.nan
    v[170:176,150:156]+=.05
    old=v.copy()
    d=block_design(g,np.isfinite(v),35.,3.)
    f=fit_geometry(d,v,g)
    assert f["status"]=="GEO1_CONDITIONAL_FIT_NOT_ADOPTED"
    assert f["east_offset_arcsec"]==pytest.approx(1.5,abs=.7)
    assert f["north_offset_arcsec"]==pytest.approx(-2.,abs=.7)
    assert f["q"]==pytest.approx(.62,abs=.08)
    assert f["pa_deg"]==pytest.approx(35.,abs=8)
    np.testing.assert_array_equal(v,old)
    rows=radial_profiles(g,{"a":v,"b":2*v},np.isfinite(v),f)
    for row in rows:
        if row["used_pixels"]:
            assert row["scenarios"]["b"]["ring_flux_nMgy"]==pytest.approx(2*row["scenarios"]["a"]["ring_flux_nMgy"])


def test_coverage_gate_and_pair_identity():
    w,p,g,a=fixture()
    d=block_design(g,np.zeros(a.shape,dtype=bool),35.)
    assert fit_geometry(d,a,g)["status"]=="GEO1_INSUFFICIENT_COVERAGE"
    with pytest.raises(ValueError,match="pair mismatch"):
        paired_geometry({}, {},w)


def test_constant_field_has_no_geometry():
    w,p,g,a=fixture()
    d=block_design(g,np.ones(a.shape,dtype=bool),35.)
    f=fit_geometry(d,np.full_like(a,-.001),g)
    assert f["status"]=="GEO1_NO_MEASURABLE_PROFILE_CONTRAST"
