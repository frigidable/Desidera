import numpy as np
import pytest
from lsbpipe.photometry.diagnostic_cache import cached,read_cache,write_cache
from lsbpipe.photometry.project_diagnostics import stage_keys


def test_cache_roundtrip_and_selective_invalidation(tmp_path):
    p=tmp_path/"cache.zip"; calls=[]
    def calculate():
        calls.append(1)
        return {"a":np.array([1.,np.nan,-2.]),"t":(4,5),"x":float("inf"),"map":{1:"a"}}
    x,s,h=cached(p,"a",calculate)
    assert s=="COMPUTED"
    y,s,h2=cached(p,"a",calculate)
    assert s=="HIT" and h==h2 and len(calls)==1
    np.testing.assert_array_equal(x["a"],y["a"])
    assert not y["a"].flags.writeable and y["t"]==(4,5) and y["map"]=={1:"a"}
    assert read_cache(p,"other") is None
    cached(p,"b",calculate)
    assert len(calls)==2
    p.write_bytes(b"broken zip")
    assert read_cache(p,"b") is None
    assert cached(p,"b",calculate)[1]=="COMPUTED"


def test_key_changes_sky_only_and_mask_changes():
    ds=dict(band="g",original_fits_sha256_local={"image":"x"},hdu_indices={"image":0},grid_fingerprint="grid")
    m=np.zeros((4,4),bool)
    c=dict(catalog_sha256="cat",mask=m.copy(),protected=m,hard=m,unresolved_target_overlap_sky_veto_mask=m)
    e=dict(source_id="s",objid=1)
    a,b=stage_keys({"p":1},ds,c,e,{"sh1":"a","sky":"b"},"artifact")
    aa,bb=stage_keys({"p":1},ds,c,e,{"sh1":"a","sky":"c"},"artifact")
    assert a==aa and b!=bb
    c["mask"][0,0]=True
    assert stage_keys({"p":1},ds,c,e,{"sh1":"a","sky":"b"},"artifact")[0]!=a


def test_failed_compute_does_not_replace_cache(tmp_path):
    p=tmp_path/"cache.zip"
    write_cache(p,"old",{"a":1})
    def fail(): raise ValueError("test")
    with pytest.raises(ValueError): cached(p,"new",fail)
    assert read_cache(p,"old")=={"a":1}


def test_headless_runner_reuses_both_stages_and_exports(monkeypatch,tmp_path):
    pytest.importorskip("PySide6")
    pytest.importorskip("astropy")
    from types import SimpleNamespace
    import json
    import zipfile
    from lsbpipe.photometry import project_diagnostics as runner
    from lsbpipe.surveys import object_projects
    from lsbpipe.photometry import mask_drafts, stellar_halo_sectors, galaxy_geometry_diagnostic, residual_sky_diagnostic
    from lsbpipe.gui import scene_editor, sky_workspace
    project=SimpleNamespace(project_id="p",root=tmp_path,data_root=tmp_path,
        selected_dataset=("dr11","north","brick","g"),target_radius_arcsec=70.,
        field=SimpleNamespace(ra_deg=96.,dec_deg=78.,width_arcmin=2.6))
    draft=SimpleNamespace(digest="d")
    monkeypatch.setattr(object_projects,"read_object_project",lambda root:project)
    monkeypatch.setattr(mask_drafts,"read_active_mask_draft",lambda root:(draft,"d"))
    a=np.arange(64,dtype=float).reshape(8,8)*.001
    mask=np.zeros(a.shape,bool)
    calls={"sh1":0,"sky":0}
    code=dict(sh1="one",sky="two",versions={})
    monkeypatch.setattr(runner,"code_identity",lambda:dict(code))
    monkeypatch.setattr(runner,"verify_sources",lambda *args:None)
    entry=dict(objid=360,source_id="s",tractor_type="PSF",gaia_dr3_source_id="1",
        foreground_evidence_assessment="GAIA_ASTROMETRY_FOREGROUND_CANDIDATE_REVIEW",
        maskbits_names=["SATUR_G","SATUR_R"])
    def load(root,band,width):
        return dict(band=band,wcs=object(),original_fits_sha256_local={"image":band},
                    hdu_indices={"image":0},grid_fingerprint="grid")
    def masks(*args,**kw):
        return dict(catalog_sha256="cat",mask=mask,protected=mask,hard=mask,
            unresolved_target_overlap_sky_veto_mask=mask,source_inventory={"entries":[entry]})
    def probe(*args,**kw):
        calls["sh1"]+=1
        return dict(status="SH1_SECTORS_MEASURED_NOT_ADOPTED",sh1c_preview=dict(
            sky1_inputs={"image":a},layers={"original":a,"residual":a}))
    def sky(*args):
        calls["sky"]+=1
        return dict(status="SKY1_NOT_ADOPTED",expanded_sky=dict(expanded_shape=[8,8],
            expanded_layers={"original":a,"SH1C":a},layers={},trials=[]))
    monkeypatch.setattr(scene_editor,"load_scene_original_cutout",load)
    monkeypatch.setattr(scene_editor,"prepare_scene_neighbour_masks",masks)
    monkeypatch.setattr(scene_editor,"run_scene_psf0_foreground_probe",probe)
    monkeypatch.setattr(sky_workspace,"run_expanded_sky",sky)
    monkeypatch.setattr(stellar_halo_sectors,"paired_sector_diagnostics",lambda *a:{"status":"TEST"})
    monkeypatch.setattr(galaxy_geometry_diagnostic,"paired_geometry",lambda *a:{"status":"TEST"})
    monkeypatch.setattr(residual_sky_diagnostic,"paired_sky_diagnostics",lambda *a:{"status":"TEST"})
    run=lambda:runner.run_project(tmp_path,objid=360,progress=lambda msg:None)
    path=run()
    assert calls=={"sh1":2,"sky":2}
    assert run()==path and calls=={"sh1":2,"sky":2}
    with zipfile.ZipFile(path) as z:
        manifest=json.loads(z.read("manifest.json"))
        assert manifest["stages"]["g"]["sh1"]["cache"]=="HIT"
        assert "g-native-sh1c.png" in z.namelist()
        assert "summary.json" in z.namelist()
        assert z.read("r-native-original.png").startswith(b"\x89PNG")
    code["sky"]="changed"
    run()
    assert calls=={"sh1":2,"sky":4}
    code["sh1"]="changed"
    run()
    assert calls=={"sh1":4,"sky":6}
