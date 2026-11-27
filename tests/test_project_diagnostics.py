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
