import numpy as np
import pytest
from lsbpipe.photometry.residual_sky_diagnostic import tile_rows,empty_samples,candidate,predict


def test_plane_recovered_without_protected_galaxy_and_heldout_apertures():
    shape=(601,601); tx=ty=300.; scale=.5
    y,x=np.indices(shape); dx=(x-tx)*scale; dy=(y-ty)*scale
    plane=.0005+.000003*dx-.000004*dy
    safe=np.hypot(dx,dy)>=60.
    a=plane.copy(); a[~safe]+=10.
    saved=a.copy()
    rows=tile_rows(a,safe,tx,ty,scale,32)
    samples=empty_samples(a,safe,tx,ty,scale)
    f=candidate(rows,"plane",samples,tx,ty,scale)
    assert f["status"]=="SKY1_CONDITIONAL_NOT_ADOPTED"
    assert f["fit"]["coefficients_nMgy_pix"]==pytest.approx([.0005,.0003,-.0004],abs=1e-10)
    assert len(f["heldout_quadrants"])==4
    assert f["cv_equal_quadrant_mae_nMgy_pix"]<1e-10
    assert any(v["count"]>=5 for v in f["empty_apertures"])
    assert max(v["rms_nMgy"] for v in f["empty_apertures"] if v["count"])<1e-7
    np.testing.assert_array_equal(a,saved)


def test_quadratic_and_missing_quadrants():
    y,x=np.indices((601,601)); dx=(x-300.)*.5; dy=(y-300.)*.5
    a=.0003+1e-8*dx**2-2e-8*dx*dy
    safe=np.hypot(dx,dy)>=60.
    rows=tile_rows(a,safe,300.,300.,.5,32)
    f=candidate(rows,"quadratic",[],300.,300.,.5)
    assert f["status"]=="SKY1_CONDITIONAL_NOT_ADOPTED"
    assert predict(f["fit"],np.array([0.,25.]),np.array([0.,10.]))==pytest.approx([.0003,.00030125],abs=1e-10)
    blocked=candidate([r for r in rows if r["quadrant"]!=0],"plane",[],300.,300.,.5)
    assert blocked["status"]=="SKY1_INSUFFICIENT_SPATIAL_SUPPORT"


def test_larger_required_target_disk_blocks_extrapolation():
    y,x=np.indices((401,401)); dx=(x-200.)*.5; dy=(y-200.)*.5
    a=.0003+1e-6*dx
    safe=np.hypot(dx,dy)>=55.
    rows=tile_rows(a,safe,200.,200.,.5,32)
    f=candidate(rows,"plane",[],200.,200.,.5,required_radius_arcsec=110.)
    assert f["status"]=="SKY1_VALIDATION_OR_HULL_BLOCKED"
    assert not f["target_region_inside_training_hull"]


def test_pair_uses_identical_blocks_and_apertures():
    from lsbpipe.photometry.residual_sky_diagnostic import paired_sky_diagnostics
    shape=(401,401); y,x=np.indices(shape); dx=(x-200)*.5; dy=(y-200)*.5
    g=.0003+1e-6*dx; r=.0006+2e-6*dy
    base=dict(project_id="p",draft_sha256="d",catalog_sha256="c",foreground_source_id="star",
        source_triplet=["dr11","north","brick"],original_grid_fingerprint="native",sh1_policy={},
        original_fits_sha256_local={"image":"hash"})
    def band(name,a,mask):
        return dict(base,band=name,expanded_sky=dict(expanded_grid_fingerprint="expanded",
            expanded_shape=list(shape),expanded_target_xy=[200.,200.],expanded_source_xy=[200.,100.],
            pixel_scale_arcsec=.5,expanded_layers={"SH1C":a},expanded_safe_mask=mask,
            required_target_support_radius_arcsec=48.))
    one=band("g",g,x>10); two=band("r",r,y>15)
    result=paired_sky_diagnostics(one,two)
    assert result["status"]=="SKY2_PAIRED_NOT_ADOPTED"
    assert result["common_pixels"]==int(((x>10)&(y>15)&(np.hypot(dx,dy)>=60.)).sum())
    assert result["bands"]["g"]["empty_aperture_centres"]==result["bands"]["r"]["empty_aperture_centres"]
    two["expanded_sky"]["expanded_grid_fingerprint"]="wrong"
    with pytest.raises(ValueError,match="grid mismatch"):
        paired_sky_diagnostics(one,two)
