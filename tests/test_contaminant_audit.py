import importlib.util
import json
from pathlib import Path
import zipfile

import pytest

spec=importlib.util.spec_from_file_location("contaminant_audit",Path(__file__).resolve().parents[1]/"scripts"/"contaminant_audit.py")
audit=importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def test_elliptical_major_minor():
    geom=dict(centre_icrs_deg=[96.,0.],q=.5,pa_deg=0.)
    assert audit.ellipse_coordinate(96.,0.,geom)==pytest.approx(0.,abs=1e-9)
    assert audit.ellipse_coordinate(96.,10/3600,geom)==pytest.approx(10.,rel=1e-6)
    assert audit.ellipse_coordinate(96.+10/3600,0.,geom)==pytest.approx(20.,rel=1e-6)


def test_union_and_no_adoption(tmp_path):
    path=tmp_path/"report.zip"
    geom=dict(centre_icrs_deg=[96.,0.],q=.5,pa_deg=0.)
    e=dict(source_id="dr11/north/b/OBJID=1",objid=1,ra_deg=96.,dec_deg=0.,
           screening_disposition="PROTECTED_OVERLAP",foreground_evidence_assessment="INSUFFICIENT_EVIDENCE",
           tractor_type="SER",model_fluxes_nMgy={"g":1.,"r":2.})
    with zipfile.ZipFile(path,"w") as z:
        for name,value in {"manifest.json":{"foreground_objid":360},
            "outer-qa.json":{"reference_geometry":geom},
            "g-qa.json":{"scene":{"catalog_sha256":"same","source_inventory":{"truncated":0,"entries":[e]}}},
            "r-qa.json":{"scene":{"catalog_sha256":"same","source_inventory":{"truncated":0,"entries":[]}}}}.items():
            z.writestr(name,json.dumps(value))
    rows=audit.audit(path)
    assert len(rows)==1
    assert rows[0]["inventory_bands"]=="g"
    assert rows[0]["decision"]=="REVIEW_NOT_ADOPTED"
    assert rows[0]["review_class"]=="EXTENDED_TARGET_OR_BACKGROUND_UNRESOLVED"
    assert rows[0]["r_action"] is None
