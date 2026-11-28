import numpy as np
import pytest
from lsbpipe.photometry.stellar_halo import measure_saturated_halo_sectors


def test_saturated_sector_profiles_preserve_signed_flux_and_target_veto():
    n=321
    yy,xx=np.indices((n,n))
    sx,sy=160.,110.
    r=np.hypot(xx-sx,yy-sy)
    image=12/(1+r)**1.5-.1
    hard=(r<10)|(abs(yy-sy)<2)
    image[hard]=1e6
    protected=(xx-160)**2+(yy-220)**2<60**2
    image[protected]+=20
    source=np.zeros_like(hard)
    source[30:50,150:170]=True
    veto=np.zeros_like(hard)
    before=image.copy()
    result=measure_saturated_halo_sectors(image,np.ones_like(image),hard,
        source_xy=(sx,sy),target_xy=(160,220),pixel_scale_arcsec=.262,
        source_mask=source,protected=protected,unresolved_sky_veto=veto,
        neighbor_centres=[(160,40)])
    assert result["status"]=="SH1_SECTORS_MEASURED_NOT_ADOPTED"
    assert not result["sky_subtracted"] and not result["halo_fitted"]
    assert not result["compact_psf_used"] and not result["catalog_flux_used"]
    assert np.array_equal(image,before)
    rows=result["rows"]
    facing=next(x for x in rows if x["sector_index"]==0 and x["r_inner_arcsec"]==25.)
    assert facing["descriptive"]["median_nMgy_pix"]>10
    assert facing["outside_protection_and_veto"]["n"]==0
    away=next(x for x in rows if x["sector_index"]==4 and x["r_inner_arcsec"]==16.)
    assert away["outside_protection_and_veto"]["median_nMgy_pix"]<0
    assert away["adequate_for_descriptive_comparison"]
    edge=next(x for x in rows if x["sector_index"]==4 and x["r_inner_arcsec"]==100.)
    assert edge["in_frame_pixels"]==0 and not edge["adequate_for_descriptive_comparison"]


def test_saturated_sector_rejects_unaligned_masks():
    a=np.zeros((10,10))
    with pytest.raises(ValueError,match="aligned"):
        measure_saturated_halo_sectors(a,a,a,source_xy=(1,1),target_xy=(5,5),
            pixel_scale_arcsec=.262,source_mask=a,protected=a[:2],
            unresolved_sky_veto=a)


def test_sh1b_block_hypotheses_exclude_target_light_and_keep_disjoint_tiles():
    n=801
    yy,xx=np.indices((n,n))
    sx,sy=400.,400.
    r=np.hypot(xx-sx,yy-sy)*.262
    halo=.012*np.clip(1-r/70.,0,1)
    image=halo-.0003+2e-6*(xx-sx)
    hard=r<10
    protected=(xx-400)**2+(yy-610)**2<100**2
    source=np.zeros_like(hard)
    def run(a):
        return measure_saturated_halo_sectors(a,np.ones_like(a),hard,
            source_xy=(sx,sy),target_xy=(400,610),pixel_scale_arcsec=.262,
            source_mask=source,protected=protected,unresolved_sky_veto=source,
            compare_models=True)["sh1b_model_comparison"]
    result=run(image)
    injected=image.copy(); injected[protected]+=10.
    other=run(injected)
    assert result["candidates"]
    assert result["candidates"]==other["candidates"]
    assert not result["science_adopted"] and not result["source_light_subtracted"]
    for group in result["block_sets"]:
        train={tuple(b["tile"]) for b in group["blocks"] if b["role"]=="TRAIN"}
        side={tuple(b["tile"]) for b in group["blocks"] if b["role"]=="SIDE_HOLDOUT"}
        assert train.isdisjoint(side)
        assert all(b["sector"] in (3,4,5) for b in group["blocks"] if b["role"]=="TRAIN")
    for c in result["candidates"]:
        h=c["halo_at_knots_nMgy_pix"]
        assert all(x>=0 for x in h)
        assert all(a>=b for a,b in zip(h,h[1:]))
        assert h[-1]==0
    assert any(c["sky_model"]=="plane" for c in result["candidates"])
    assert result["target_distance_hypothesis_range_nMgy_pix"] is not None
    c=result["candidates"][0]
    key=(c["tile_size_pix"],c["assumed_zero_halo_radius_arcsec"],c["sky_model"])
    fixed=measure_saturated_halo_sectors(image,np.ones_like(image),hard,
        source_xy=(sx,sy),target_xy=(400,610),pixel_scale_arcsec=.262,
        source_mask=source,protected=protected,unresolved_sky_veto=source,
        compare_models=True,fixed_hypothesis=key)["sh1b_model_comparison"]
    assert fixed["candidates"]==[c]
    assert len(fixed["block_sets"])==1


def test_sh1b2_extends_radial_support_and_emits_nonadopted_tile_map():
    n=741
    yy,xx=np.indices((n,n))
    sx,sy=370.,370.
    r=np.hypot(xx-sx,yy-sy)*.5
    image=.01*np.clip(1-r/120.,0,1)-.0002
    zero=np.zeros_like(image,dtype=bool)
    protected=(xx-370)**2+(yy-480)**2<60**2
    result=measure_saturated_halo_sectors(image,np.ones_like(image),r<12.,
        source_xy=(sx,sy),target_xy=(370,480),pixel_scale_arcsec=.5,
        source_mask=zero,protected=protected,unresolved_sky_veto=zero,
        compare_models=True,expanded=True)
    assert result["radial_edges_arcsec"][-1]==180.
    assert len(result["rows"])==120
    comp=result["sh1b_model_comparison"]
    assert any(c["assumed_zero_halo_radius_arcsec"]==150. for c in comp["candidates"])
    assert all(c["training_r_range_arcsec"][1]>160. for c in comp["candidates"])
    m=comp["diagnostic_tile_map"]
    assert m["cells"] and not m["science_adopted"]
    assert any(c["protected_fraction_of_used"]>0 for c in m["cells"])
    assert all(c["halo_median_range_nMgy_pix"][0]>=0 for c in m["cells"])


def test_expanded_halo_crop_clips_without_padding_and_preserves_wcs():
    from astropy.wcs import WCS
    from types import SimpleNamespace
    from lsbpipe.photometry.stellar_halo_scene import expanded_halo_cutout
    from lsbpipe.scene import grid_fingerprint
    wcs=WCS(naxis=2)
    wcs.wcs.crpix=[501.,501.]
    wcs.wcs.cdelt=[-.262/3600,.262/3600]
    wcs.wcs.crval=[96.,78.]
    wcs.wcs.ctype=["RA---TAN","DEC--TAN"]
    a=np.ones((1000,1000))
    hashes={k:"a"*64 for k in ("image","invvar","nexp","maskbits")}
    hdus={k:1 for k in hashes}
    frame=SimpleNamespace(dataset_id="dr11/north/0965p780/r",wcs=wcs,
        shape=a.shape,image_nMgy=a,invvar_nMgy2=a,nexp=a,maskbits=np.zeros(a.shape,dtype=np.uint32),
        source_evidence={"local_hashes":hashes,"hdu_indices":hdus})
    cut=wcs.slice(np.s_[200:800,400:900])
    scene={"source_quartet_key":["dr11","north","0965p780","r"],
           "original_fits_sha256_local":hashes,"hdu_indices":hdus,
           "wcs":cut,"pixel_shape":(600,500),"source":"native",
           "project_id":"p","grid_fingerprint":grid_fingerprint(cut,(600,500))}
    star=tuple(wcs.all_pix2world(850.,500.,0))
    expanded,qa=expanded_halo_cutout(frame,scene,star)
    assert qa["clipped_by_brick"] and qa["actual_bounds_xyxy_exclusive"][2]==1000
    assert qa["scene_origin_full_brick_xy"]==[400,200]
    assert qa["original_scene_grid_fingerprint"]==scene["grid_fingerprint"]
    assert expanded["data"]["image"].shape==tuple(qa["native_shape"])
    assert not expanded["data"]["image"].flags.writeable
    assert not expanded["native_hard"].any()
    assert not np.shares_memory(a,expanded["data"]["image"])
    with pytest.raises(ValueError,match="hashes differ"):
        expanded_halo_cutout(frame,dict(scene,original_fits_sha256_local={}),star)
