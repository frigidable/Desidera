import numpy as np
import pytest
pytest.importorskip("astropy")
from astropy.wcs import WCS
from lsbpipe.photometry.sky_clean_injection import saturations,compact_peaks,select_sites,galaxy_image,inject


def test_saturation_is_independent_of_flux_ivar():
    points=np.array([[50.,50.],[100.,50.]])
    rows=[dict(objid=9,maskbits=4,flux=500.,flux_ivar=0.,brick_primary=False),
          dict(objid=10,maskbits=8,flux=1000.,flux_ivar=0.)]
    m,s=saturations((151,151),points,rows,"g",.5)
    assert m[50,50] and not m[50,100]
    assert [r["objid"] for r in s]==[9]
    assert s[0]["flux_ivar"]==0.


def test_compact_detection_and_negative_control():
    rng=np.random.default_rng(42)
    a=rng.normal(0,.001,(201,201))
    a[60:65,60:65]+=.2
    a[140:145,140:145]-=.2
    mask=np.ones_like(a,bool)
    saved=a.copy()
    grown,qa=compact_peaks(a,mask,.5)
    assert grown[62,62] and not grown[142,142]
    assert qa["positive_components"]>=1 and qa["negative_components"]>=1
    np.testing.assert_array_equal(a,saved)


def test_sites_and_exponential_truth():
    m=np.ones((601,601),bool)
    sites=select_sites(m,300.,300.,.5)
    assert len(sites)==2
    assert all(np.hypot(s["x_pix"]-300,s["y_pix"]-300)*.5>=100 for s in sites)
    w=WCS(naxis=2); w.wcs.crpix=[301.,301.]; w.wcs.crval=[96.,78.]
    w.wcs.cdelt=[-.5/3600,.5/3600]; w.wcs.ctype=["RA---TAN","DEC--TAN"]
    site=dict(x_pix=300,y_pix=300)
    im,rho=galaxy_image(m.shape,w,site,(15.,24.5,.56,44.5),.5)
    assert im[300,300]==pytest.approx(.25*10**(.4*(22.5-24.5)),rel=1e-3)
    assert np.sum(im[rho<=15])/im.sum()==pytest.approx(.5,abs=.015)


def test_zero_sky_differential_recovery_no_noise_and_refit():
    shape=(701,701); y,x=np.indices(shape)
    w=WCS(naxis=2); w.wcs.crpix=[351.,351.]; w.wcs.crval=[96.,78.]
    w.wcs.cdelt=[-.5/3600,.5/3600]; w.wcs.ctype=["RA---TAN","DEC--TAN"]
    z=np.zeros(shape); m=np.ones(shape,bool)
    site=dict(x_pix=145,y_pix=350)
    f=inject(z,m,w,site,(15.,24.5,.56,44.5),350.,350.,.5,48.)
    a=f["models"]["zero"]
    assert a["geometry"]["status"]=="GEO1_CONDITIONAL_FIT_NOT_ADOPTED"
    assert a["geometry_errors"]["delta_q"]==pytest.approx(0,abs=.04)
    assert a["geometry_errors"]["delta_pa_deg"]==pytest.approx(0,abs=3)
    assert all(abs(v["differential_fractional_error"])<1e-12 for v in a["apertures"])
    assert f["models"]["constant"]["sky_after"] is not None
