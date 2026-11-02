"""Real offline Legacy brick-CCDs FITS ingestion with conservative provenance."""
from __future__ import annotations

import json
import hashlib

import numpy as np
import pytest
from astropy.io import fits

from lsbpipe.surveys.dataset_registry import (
    import_legacy_fits, read_dataset_index, registry_path,
    verify_asset_checksum,
)
from lsbpipe.surveys.legacy import BrickHit, SkyField, product_for
from lsbpipe.surveys.object_projects import create_object_project
from lsbpipe.surveys.observation_lineage import (
    compare_coadd_lineage, lineage_path, read_observation_lineage,
)
from lsbpipe.surveys.legacy_ccd_catalog import (
    inspect_brick_ccd_catalog, import_brick_ccd_catalog,
)


@pytest.fixture
def project(tmp_path, monkeypatch):
    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path / "sessions"))
    return create_object_project("Leo P", tmp_path / "Leo_P",
             field=SkyField(155.438012, 18.088025, 5))


def add_original_coadds(project):
    data = project.root / "data"
    for release in ("dr9", "dr11"):
        hit = BrickHit(release, "south", "1554p180", 0, 0, 0, 0, {"r": 3})
        for kind in ("image", "invvar", "nexp", "maskbits"):
            path = data / product_for(
                hit, kind, None if kind == "maskbits" else "r").relative_path
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(f"TEST_SOURCE:{release}:{kind}".encode())
    import_legacy_fits(project.root, data)
    index = read_dataset_index(project.root)
    sources = {s["source_id"]: s["release"] for s in index["sources"]}
    images = {
        sources[a["source_id"]]: a["asset_id"] for a in index["assets"]
        if a["product_kind"] == "image"}
    for aid in images.values():
        verify_asset_checksum(project.root, aid)
    return images


def make_catalog(tmp_path, *, brick="1554p180", rows=None,
                 missing=(), extra_hdu=False, bad_header=False):
    if rows is None:
        rows = [
            ("decam", 574299, "N10", "r", 0, 20, 4, 25),
            ("decam", 574300, "N12", "g", 0, 20, 4, 25),
            ("decam", 574301, "N13", "z", 0, 20, 4, 25),
        ]
    rows = list(rows)
    names = ("camera", "expnum", "ccdname", "filter", "brick_x0",
             "brick_x1", "brick_y0", "brick_y1")
    types = ("6A", "K", "6A", "1A", "J", "J", "J", "J")
    cols = []
    for pos, (name, kind) in enumerate(zip(names, types)):
        if name not in missing:
            cols.append(fits.Column(
                name=name.upper(), format=kind,
                array=np.array([r[pos] for r in rows])))
    hdus = [fits.PrimaryHDU(), fits.BinTableHDU.from_columns(cols)]
    if extra_hdu:
        hdus.append(fits.BinTableHDU.from_columns(cols))
    if bad_header:
        hdus[0].header["CHECKSUM"] = "wrong-hash"
    path = tmp_path / f"legacysurvey-{brick}-ccds.fits"
    fits.HDUList(hdus).writeto(path, overwrite=True)
    return path


def test_real_fits_import_cross_release_receipts_are_idempotent(project, tmp_path):
    images = add_original_coadds(project)
    file = make_catalog(tmp_path)
    before_catalogue = file.read_bytes()
    before_object = (project.root / "lsbpipe-object.json").read_bytes()
    before_dataset = registry_path(project.root).read_bytes()
    assert inspect_brick_ccd_catalog(
        project.root, coadd_asset_id=images["dr9"],
        catalogue_file=file)["selected_ccd_count"] == 1
    assert not lineage_path(project.root).exists()

    a = import_brick_ccd_catalog(
        project.root, coadd_asset_id=images["dr9"], catalogue_file=file)
    b = import_brick_ccd_catalog(
        project.root, coadd_asset_id=images["dr11"], catalogue_file=file)
    assert a["observations_added"] == 1
    assert b["observations_added"] == 0
    assert a["catalogue_receipt_added"]
    assert b["catalogue_receipt_added"]
    assert a["catalogue_sha256"] == hashlib.sha256(before_catalogue).hexdigest()
    state = read_observation_lineage(project.root)
    assert len(state["catalogue_receipts"]) == 2
    assert all(r["status"] == "CATALOGUE_ASSERTION_NOT_VALIDATED"
               for r in state["catalogue_receipts"])
    res = compare_coadd_lineage(
        project.root, images["dr9"], images["dr11"])
    assert res["status"] == "SHARED_CCD_IDS_ASSERTED"
    assert not res["independence_certified"]
    assert not res["joint_scientific_combination_allowed"]
    snapshot = lineage_path(project.root).read_bytes()
    again = import_brick_ccd_catalog(
        project.root, coadd_asset_id=images["dr9"], catalogue_file=file)
    assert again["observations_added"] == 0
    assert again["contributions_added"] == 0
    assert not again["catalogue_receipt_added"]
    assert lineage_path(project.root).read_bytes() == snapshot
    assert file.read_bytes() == before_catalogue
    assert registry_path(project.root).read_bytes() == before_dataset
    assert (project.root / "lsbpipe-object.json").read_bytes() == before_object


@pytest.mark.parametrize("kwargs,match", [
    ({"brick": "1554p182"}, "filename brick"),
    ({"missing": ("camera",)}, "mandatory columns"),
    ({"extra_hdu": True}, "exactly one"),
    ({"rows": []}, "Empty"),
    ({"rows": [("decam", 1, "N10", "r", 0, 0, 3, 4)]},
     "nonpositive brick overlap"),
    ({"rows": [("decam", 1, "N10", "r", 0, 10, 3, 4),
               ("decam", 1, "N10", "r", 0, 10, 3, 4)]},
     "Duplicate physical CCD"),
    ({"rows": [("decam", 1, "N10", "g", 0, 10, 3, 4)]},
     "No CCD catalogue rows"),
    ({"rows": [("UNKNOWN", 1, "N10", "r", 0, 10, 3, 4)]},
     "Unrecognized CAMERA"),
])
def test_fits_fail_closed_and_leave_no_partial_metadata(
        project, tmp_path, kwargs, match):
    images = add_original_coadds(project)
    file = make_catalog(tmp_path, **kwargs)
    with pytest.raises(ValueError, match=match):
        import_brick_ccd_catalog(
            project.root, coadd_asset_id=images["dr9"],
            catalogue_file=file)
    assert not lineage_path(project.root).exists()


def test_modified_catalogue_does_not_replace_old_lineage(project, tmp_path):
    images = add_original_coadds(project)
    file = make_catalog(tmp_path)
    import_brick_ccd_catalog(
        project.root, coadd_asset_id=images["dr9"], catalogue_file=file)
    old = lineage_path(project.root).read_bytes()
    make_catalog(tmp_path, rows=[("decam", 999, "N11", "r", 0, 10, 0, 10)])
    with pytest.raises(ValueError, match="catalogue changed"):
        import_brick_ccd_catalog(
            project.root, coadd_asset_id=images["dr9"], catalogue_file=file)
    assert lineage_path(project.root).read_bytes() == old


def test_corrupted_catalogue_receipt_is_detected(project, tmp_path):
    images = add_original_coadds(project)
    file = make_catalog(tmp_path)
    import_brick_ccd_catalog(
        project.root, coadd_asset_id=images["dr9"], catalogue_file=file)
    manifest = lineage_path(project.root)
    raw = json.loads(manifest.read_text(encoding="utf-8"))
    raw["catalogue_receipts"][0]["catalogue_sha256"] = "0" * 64
    manifest.write_text(json.dumps(raw), encoding="utf-8")
    with pytest.raises(ValueError, match="reference does not match"):
        read_observation_lineage(project.root)


def test_only_verified_original_image_is_accepted(project, tmp_path):
    images = add_original_coadds(project)
    file = make_catalog(tmp_path)
    index = read_dataset_index(project.root)
    invvar = next(a["asset_id"] for a in index["assets"]
                  if a["product_kind"] == "invvar")
    with pytest.raises(ValueError, match="coadd IMAGE"):
        import_brick_ccd_catalog(
            project.root, coadd_asset_id=invvar, catalogue_file=file)
    raw = json.loads(registry_path(project.root).read_text(encoding="utf-8"))
    a = next(x for x in raw["assets"] if x["asset_id"] == images["dr9"])
    a["checksum_status"] = "NOT_HASHED"
    registry_path(project.root).write_text(json.dumps(raw), encoding="utf-8")
    with pytest.raises(ValueError, match="SHA-256 must be verified"):
        import_brick_ccd_catalog(
            project.root, coadd_asset_id=images["dr9"], catalogue_file=file)
    assert not lineage_path(project.root).exists()
