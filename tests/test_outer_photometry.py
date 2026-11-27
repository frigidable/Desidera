import numpy as np
import pytest
pytest.importorskip("astropy")
from lsbpipe.photometry.outer_photometry import outer_rows,plateau,hull_support
from lsbpipe.photometry.sersic_recovery import fit_profile, radius_pixels, convolved_profile
from lsbpipe.photometry.galaxy_geometry_diagnostic import local_grid
from astropy.wcs import WCS


def test_outer_flux_signed_and_missing_sector_blocks_plateau():
    y,x=np.mgrid[-80:81,-80:81]
    g=dict(x=x,e=x,n=y,area=1.)
    geom=dict(east_offset_arcsec=0.,north_offset_arcsec=0.,q=1.,pa_deg=0.)
    mask=np.ones_like(x,bool)
    a=np.where(np.hypot(x,y)<40,1.,-.1)
    rows=outer_rows(g,mask,{"g":a},geom,0.,1.)
    assert rows[-1]["scenarios"]["g"]["ring_flux_nMgy"]<0
    assert plateau(rows,"g")["supported_outer_arcsec"]==70
    mask[(x>0)&(y>0)&(np.hypot(x,y)>40)]=False
    rows=outer_rows(g,mask,{"g":a},geom,0.,1.)
    assert plateau(rows,"g")["supported_outer_arcsec"]==40
    assert not rows[-1]["support_screen"]


def test_hull_requires_entire_disk():
    rows=[dict(x_pix=x,y_pix=y) for x,y in [(-80,-80),(-80,80),(80,-80),(80,80)]]
    assert hull_support(rows,[(0,0)],[40,70,90],1.)==[40.,70.]
    assert hull_support(rows,[(60,0)],[40],1.)==[]


def test_profile_extends_observed_aperture_to_outer_radius():
    w=WCS(naxis=2); w.wcs.crpix=[181.,181.]; w.wcs.crval=[96.,78.]
    w.wcs.cdelt=[-.5/3600,.5/3600]; w.wcs.ctype=["RA---TAN","DEC--TAN"]
    a=np.zeros((361,361))
    g=local_grid(dict(target_xy_native_pix=[180.,180.],layers={"original":a}),w,margin_arcsec=80.)
    geom=dict(east_offset_arcsec=0.,north_offset_arcsec=0.,q=.56,pa_deg=44.5)
    k=np.ones((3,3))/9
    r=radius_pixels(g["x"].shape,(180.,180.),g["k"],.56,44.5,origin=(g["x"][0,0],g["y"][0,0]),pad=(1,1))
    im=.12*g["area"]*convolved_profile(r,26.,.5,k)
    fit=fit_profile(im,np.ones_like(im,bool),g,geom,k,outer=60.)
    assert fit["status"]=="SERSIC_CONDITIONAL_NOT_ADOPTED"
    assert fit["apertures"][-1]["a_arcsec"]==60.
    assert fit["re_arcsec"]==pytest.approx(26.,rel=.02)
