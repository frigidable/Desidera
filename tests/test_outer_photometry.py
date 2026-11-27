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
    fit=fit_profile(im,np.ones_like(im,bool),g,geom,k,outer=60.,validation_outer=70.)
    assert fit["status"]=="SERSIC_CONDITIONAL_NOT_ADOPTED"
    assert fit["apertures"][-1]["a_arcsec"]==60.
    assert fit["re_arcsec"]==pytest.approx(26.,rel=.02)

    assert len(fit["heldout_outer_rings"])==2
    assert all(not row["used_in_profile_fit"] for row in fit["heldout_outer_rings"])
    assert max(abs(row["observed_minus_model_nMgy"]) for row in fit["heldout_outer_rings"])<.01


def test_native_overlap_compares_only_shared_pixels():
    from lsbpipe.photometry.outer_photometry import native_overlap
    y,x=np.mgrid[0:21,0:21]
    a=np.ones_like(x,dtype=float)
    old=a.copy(); old[0,0]=np.nan
    mask=np.ones_like(x,bool); mask[1,1]=False
    grid=dict(x=x,y=y,e=x-10,n=y-10)
    geom=dict(east_offset_arcsec=0.,north_offset_arcsec=0.,q=1.,pa_deg=0.)
    p=dict(target_xy_native_pix=[10.,10.],layers={"residual":old})
    z=native_overlap(grid,mask,a,p,geom,10.,10.)
    assert z["identical_pixels"]==439
    assert z["status"]=="OUTER_NATIVE_OVERLAP_AGREES"
    a[2,2]+=1.
    assert native_overlap(grid,mask,a,p,geom,10.,10.)["status"]=="OUTER_NATIVE_OVERLAP_REVIEW"


def test_halo_support_is_per_scenario_not_global():
    from lsbpipe.photometry.outer_photometry import scenario_photo_masks
    y,x=np.mgrid[-80:81,-80:81]
    base=np.ones_like(x,bool)
    base[80,80]=False  # native hard/source exclusion must survive all trials
    ref=np.ones_like(x,float)
    g=ref.copy(); r=ref.copy()
    g[np.hypot(x,y+56)<38]=np.nan
    r[np.hypot(x,y+56)<36]=np.nan
    masks=scenario_photo_masks(base,{"g":{"reference":ref,"protect_90":g},
                                     "r":{"reference":ref,"protect_90":r}})
    assert np.array_equal(masks["reference"],base)
    assert np.array_equal(masks["protect_90"],base&np.isfinite(g)&np.isfinite(r))
    assert np.array_equal(base,masks["reference"])
    grid=dict(x=x,e=x,n=y,area=1.)
    geom=dict(east_offset_arcsec=0.,north_offset_arcsec=0.,q=1.,pa_deg=0.)
    rows=outer_rows(grid,masks["reference"],{"g":ref,"r":ref},geom,0.,1.)
    assert plateau(rows,"g")["supported_outer_arcsec"]==70.
    rows=outer_rows(grid,masks["protect_90"],{"g":g,"r":r},geom,0.,1.)
    assert plateau(rows,"g")["supported_outer_arcsec"]<40.
    assert np.isfinite([row["scenarios"]["g"]["ring_flux_nMgy"] for row in rows]).all()
    assert not masks["reference"][80,80] and not masks["protect_90"][80,80]


def test_native_overlap_reports_support_loss_separately_from_residual_delta():
    from lsbpipe.photometry.outer_photometry import native_overlap
    y,x=np.mgrid[0:21,0:21]
    a=np.ones_like(x,float); mask=np.ones_like(x,bool)
    mask[1,1]=False
    grid=dict(x=x+7,y=y+11,e=x-10,n=y-10)
    geom=dict(east_offset_arcsec=0.,north_offset_arcsec=0.,q=1.,pa_deg=0.)
    p=dict(target_xy_native_pix=[10.,10.],layers={"residual":a})
    z=native_overlap(grid,mask,a,p,geom,17.,21.)
    assert z["max_abs_pixel_delta_nMgy"]==0.
    assert z["native_only_pixels"]==1
    assert not z["support_agrees_in_native_frame"]
    mask[1,1]=True
    assert native_overlap(grid,mask,a,p,geom,17.,21.)["support_agrees_in_native_frame"]
