"""P0.6: metadata-only CCD lineage, release overlap and evidence tests."""
from __future__ import annotations

import json

import pytest

from lsbpipe.surveys.dataset_registry import (
    import_legacy_fits, read_dataset_index, registry_path,
)
from lsbpipe.surveys.legacy import BrickHit, SkyField, product_for
from lsbpipe.surveys.object_projects import create_object_project
from lsbpipe.surveys.observation_lineage import (
    CCDIdentity, UNVERIFIED, compare_coadd_lineage, lineage_path,
    read_observation_lineage, register_ccd_assertions,
)


@pytest.fixture
def field_project(tmp_path, monkeypatch):
    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path / "sessions"))
    return create_object_project(
        "Leo P", tmp_path / "Leo_P", field=SkyField(155.438012, 18.088025, 5.0))


def _brick_fits(root, *, release, region="south", band="r",
                brick="1554p180"):
    hit = BrickHit(release, region, brick, 0, 0, 0, 0, {band: 3})
    paths = {}
    for kind in ("image", "invvar", "nexp", "maskbits"):
        rel = product_for(hit, kind,
                          None if kind == "maskbits" else band).relative_path
        path = root / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(
            f"ORIGINAL_UNMODIFIED:{release}/{brick}/{kind}/{band}".encode())
        paths[kind] = path
    return paths


def _two_images(project):
    root = project.root / "data"
    originals = {
        "dr9": _brick_fits(root, release="dr9"),
        "dr11": _brick_fits(root, release="dr11"),
    }
    import_legacy_fits(project.root, root)
    index = read_dataset_index(project.root)
    releases = {s["source_id"]: s["release"] for s in index["sources"]}
    images = {releases[a["source_id"]]: a["asset_id"]
              for a in index["assets"] if a["product_kind"] == "image"}
    return images, originals


def _row(exposure="0001234", detector="N12", locator="row:19",
         band="r"):
    return {"archive": "NOIRLab", "instrument": "DECam",
            "exposure_id": exposure, "detector_id": detector,
            "band": band, "record_locator": locator}


def _assertion(project, image_id, row=None):
    return register_ccd_assertions(
        project.root, coadd_asset_id=image_id,
        catalogue_reference="Legacy brick-CCDs, row source unverified",
        records=[_row() if row is None else row])


def test_no_sidecar_is_unknown_and_pure_read(field_project):
    path = lineage_path(field_project.root)
    assert not path.exists()
    state = read_observation_lineage(field_project.root)
    assert state["project_id"] == field_project.project_id
    assert not state["observations"] and not state["contributions"]
    assert state["scientific_status"] == "NOT_SCIENTIFICALLY_ADOPTED"
    assert not path.exists()


def test_immutable_ccd_identity_release_independent_and_preserve_leading_zero():
    a = CCDIdentity("NOIRLab", "DECam", "0001234", "n12")
    b = CCDIdentity("noirlab", "decam", "0001234", "N12")
    c = CCDIdentity("noirlab", "decam", "1234", "N12")
    assert a == b and a.ccd_id == b.ccd_id
    assert a.ccd_id != c.ccd_id
    assert a.exposure_id_stable == b.exposure_id_stable
    assert c.exposure_id_stable != a.exposure_id_stable
    with pytest.raises(ValueError, match="preserve leading zeros"):
        CCDIdentity("NOIRLab", "DECam", 1234, "N12")


def test_dr9_dr11_same_ccd_preserves_origins_and_is_idempotent(field_project):
    images, files = _two_images(field_project)
    before = {str(path): path.read_bytes()
              for group in files.values() for path in group.values()}
    project_bytes = (field_project.root / "lsbpipe-object.json").read_bytes()
    original_index = registry_path(field_project.root).read_bytes()
    a = _assertion(field_project, images["dr9"])
    b = _assertion(field_project, images["dr11"])
    assert a["observations_added"] == 1
    assert b["observations_added"] == 0
    assert a["evidence_state"] == UNVERIFIED
    result = compare_coadd_lineage(
        field_project.root, images["dr9"], images["dr11"])
    assert result["status"] == "SHARED_CCD_IDS_ASSERTED"
    assert len(result["shared_ccd_ids"]) == 1
    assert len(result["shared_exposure_ids"]) == 1
    assert result["independence_certified"] is False
    assert result["joint_scientific_combination_allowed"] is False
    assert len(read_observation_lineage(field_project.root)["contributions"]) == 2

    sidecar_before = lineage_path(field_project.root).read_bytes()
    repeated = _assertion(field_project, images["dr9"])
    assert repeated["observations_added"] == 0
    assert repeated["contributions_added"] == 0
    assert lineage_path(field_project.root).read_bytes() == sidecar_before
    assert registry_path(field_project.root).read_bytes() == original_index
    assert (field_project.root / "lsbpipe-object.json").read_bytes() == project_bytes
    assert {str(path): path.read_bytes()
            for group in files.values() for path in group.values()} == before


def test_shared_exposure_but_different_detector_is_not_same_ccd(field_project):
    images, _ = _two_images(field_project)
    _assertion(field_project, images["dr9"], _row(detector="N12"))
    _assertion(field_project, images["dr11"], _row(detector="N13"))
    result = compare_coadd_lineage(
        field_project.root, images["dr9"], images["dr11"])
    assert result["status"] == "SHARED_EXPOSURE_IDS_ASSERTED"
    assert result["shared_ccd_ids"] == []
    assert len(result["shared_exposure_ids"]) == 1
    assert not result["independence_certified"]


def test_missing_and_distinct_lineage_never_certify_independence(field_project):
    images, _ = _two_images(field_project)
    missing = compare_coadd_lineage(
        field_project.root, images["dr9"], images["dr11"])
    assert missing["status"] == "UNKNOWN_LINEAGE_NOT_INDEXED"
    _assertion(field_project, images["dr9"], _row(exposure="001"))
    still_missing = compare_coadd_lineage(
        field_project.root, images["dr9"], images["dr11"])
    assert still_missing["status"] == "UNKNOWN_LINEAGE_NOT_INDEXED"
    _assertion(field_project, images["dr11"], _row(exposure="002"))
    distinct = compare_coadd_lineage(
        field_project.root, images["dr9"], images["dr11"])
    assert distinct["status"] == "NO_SHARED_REGISTERED_IDS_NOT_CERTIFIED"
    for result in (missing, still_missing, distinct):
        assert result["independence_certified"] is False
        assert result["joint_scientific_combination_allowed"] is False


def test_invalid_batch_has_no_partial_sidecar_or_fits_edits(field_project):
    images, files = _two_images(field_project)
    with pytest.raises(ValueError, match="band conflicts"):
        register_ccd_assertions(
            field_project.root, coadd_asset_id=images["dr9"],
            catalogue_reference="test", records=[_row(), _row(band="z")])
    assert not lineage_path(field_project.root).exists()
    with pytest.raises(ValueError, match="non-coadd"):
        index = read_dataset_index(field_project.root)
        invvar = next(a["asset_id"] for a in index["assets"]
                      if a["product_kind"] == "invvar")
        register_ccd_assertions(
            field_project.root, coadd_asset_id=invvar,
            catalogue_reference="test", records=[_row()])
    assert files["dr9"]["image"].is_file()


def test_conflicting_locator_and_tampered_sidecar_fail_closed(field_project):
    images, _ = _two_images(field_project)
    _assertion(field_project, images["dr9"])
    sidecar = lineage_path(field_project.root)
    before = sidecar.read_bytes()
    with pytest.raises(ValueError, match="Conflicting CCD"):
        _assertion(field_project, images["dr9"], _row(exposure="different"))
    assert sidecar.read_bytes() == before
    raw = json.loads(sidecar.read_text(encoding="utf-8"))
    raw["project_id"] = "different-project"
    sidecar.write_text(json.dumps(raw), encoding="utf-8")
    with pytest.raises(ValueError, match="different galaxy"):
        read_observation_lineage(field_project.root)
    raw["project_id"] = field_project.project_id
    raw["observations"][0]["identity"]["exposure_id"] = "forged"
    sidecar.write_text(json.dumps(raw), encoding="utf-8")
    with pytest.raises(ValueError, match="mismatch"):
        read_observation_lineage(field_project.root)


def test_same_asset_is_not_independent_measurement(field_project):
    images, _ = _two_images(field_project)
    result = compare_coadd_lineage(
        field_project.root, images["dr9"], images["dr9"])
    assert result["status"] == "SAME_ASSET"
    assert result["independence_certified"] is False
