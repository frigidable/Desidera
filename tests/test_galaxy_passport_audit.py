"""P10h saved dataset QA: specific, historical and never silently adopted."""
import hashlib
import json
import os
from pathlib import Path

import pytest

from lsbpipe.surveys.galaxy_passport_audit import inspect_saved_dataset_audit
from lsbpipe.surveys.dataset_integrity import SCHEMA
from lsbpipe.surveys.legacy import BrickHit, SkyField, product_for
from lsbpipe.surveys.object_projects import create_object_project


def _fixture(tmp_path, *, selected=("dr11", "south", "1554p180", "r")):
    project = create_object_project(
        "Leo P", tmp_path / "LeoP", field=SkyField(155.438012, 18.088025, 4))
    data = project.root / "data"
    data.mkdir()
    project.data_root = data
    project.selected_dataset = selected
    return project


def _audit(project, *, record_status="FITS_WCS_VALIDATED_DIAGNOSTIC"):
    key = project.selected_dataset
    hit = BrickHit(key[0], key[1], key[2], 0, 0, 0, 0, {key[3]: None})
    products = {}
    for kind in ("image", "invvar", "nexp", "maskbits"):
        product = product_for(hit, kind, key[3] if kind != "maskbits" else None)
        file = project.data_root / product.relative_path
        file.parent.mkdir(parents=True, exist_ok=True)
        file.write_bytes((kind + " - test original contents").encode("ascii"))
        st = file.stat()
        products[kind] = {
            "relative_path": product.relative_path.as_posix(),
            "size_bytes": st.st_size, "mtime_ns": st.st_mtime_ns,
            "sha256_local": hashlib.sha256(file.read_bytes()).hexdigest(),
        }
    dataset = {
        "release": key[0], "region": key[1], "brick": key[2],
        "band": key[3], "status": record_status,
        "products": products,
        "fractions": {"instrumentally_usable": 0.9975},
    }
    stamp = "2026-11-03T11:00:00+00:00"
    receipt = {
        "schema": SCHEMA, "project_id": project.project_id,
        "science_status": "DIAGNOSTIC_NOT_ADOPTED",
        "created_utc": stamp, "data_root": str(project.data_root),
        "field": {
            "ra_deg": project.field.ra_deg,
            "dec_deg": project.field.dec_deg,
            "width_arcmin": project.field.width_arcmin,
        },
        "datasets": [dataset],
    }
    folder = project.root / "dataset-audits" / "run123"
    folder.mkdir(parents=True)
    out = folder / "audit.json"
    out.write_text(json.dumps(receipt), encoding="utf-8")
    project.last_dataset_audit = {
        "relative_path": "dataset-audits/run123/audit.json",
        "created_utc": stamp, "data_root": str(project.data_root),
    }
    return receipt, out


def test_no_saved_audit_never_means_fits_validated(tmp_path):
    p = _fixture(tmp_path)
    result = inspect_saved_dataset_audit(p)
    assert result.state == "NO_SAVED_AUDIT"
    assert "No saved FITS/WCS" in result.label


def test_past_validated_quartet_matches_stats_not_rehashed(tmp_path, monkeypatch):
    p = _fixture(tmp_path)
    report, path = _audit(p)
    original = path.read_bytes()
    import astropy.io.fits
    monkeypatch.setattr(astropy.io.fits, "open", lambda *args, **kwargs:
                        pytest.fail("The passport must never decode FITS"))
    result = inspect_saved_dataset_audit(p)
    assert result.state == "AUDIT_RECORDED_STAT_MATCH"
    assert "metadata unchanged" in result.label
    assert "NOT rechecked" in result.explanation
    assert "NOT certified" in result.explanation
    assert result.instrumentally_usable_fraction == 0.9975
    assert result.receipt_relative_path.endswith("audit.json")
    assert path.read_bytes() == original
    assert result.created_utc == report["created_utc"]


def test_stat_equal_does_not_claim_content_hash_was_rechecked(tmp_path):
    p = _fixture(tmp_path)
    report, _ = _audit(p)
    image = p.data_root / report["datasets"][0]["products"]["image"]["relative_path"]
    old = image.stat()
    image.write_bytes(b"Z" * old.st_size)
    os.utime(image, ns=(old.st_atime_ns, old.st_mtime_ns))
    result = inspect_saved_dataset_audit(p)
    assert result.state == "AUDIT_RECORDED_STAT_MATCH"
    assert "hashes were NOT rechecked" in result.explanation
    assert "scientific suitability" in result.explanation


def test_changed_original_is_not_current_validated_qa(tmp_path):
    p = _fixture(tmp_path)
    report, _ = _audit(p)
    image = p.data_root / report["datasets"][0]["products"]["image"]["relative_path"]
    image.write_bytes(b"changed size invalidated old report")
    result = inspect_saved_dataset_audit(p)
    assert result.state == "AUDIT_INPUTS_CHANGED"
    assert "needs rerun" in result.label


@pytest.mark.parametrize("change", [
    "project_id", "schema", "science_status", "timestamp",
    "source", "input_hash", "product_path", "audit_path",
])
def test_corrupt_or_unbound_evidence_fails_closed(tmp_path, change):
    p = _fixture(tmp_path)
    report, path = _audit(p)
    if change == "project_id":
        report["project_id"] = "other-galaxy"
    elif change == "schema":
        report["schema"] = "unsupported"
    elif change == "science_status":
        report["science_status"] = "SCIENCE_CERTIFIED"
    elif change == "timestamp":
        report["created_utc"] = "changed"
    elif change == "source":
        report["data_root"] = str(tmp_path / "other-source")
    elif change == "input_hash":
        report["datasets"][0]["products"]["image"]["sha256_local"] = "fake"
    elif change == "product_path":
        report["datasets"][0]["products"]["image"]["relative_path"] = "../x"
    elif change == "audit_path":
        p.last_dataset_audit["relative_path"] = "dataset-audits/../audit.json"
    if change != "audit_path":
        path.write_text(json.dumps(report), encoding="utf-8")
    result = inspect_saved_dataset_audit(p)
    assert result.state in ("AUDIT_EVIDENCE_UNAVAILABLE", "AUDIT_INPUTS_CHANGED")
    assert result.state != "AUDIT_RECORDED_STAT_MATCH"


def test_field_change_or_dataset_switch_does_not_reuse_qa(tmp_path):
    p = _fixture(tmp_path)
    _audit(p)
    p.field = SkyField(155.438, 18.085, 4)
    assert inspect_saved_dataset_audit(p).state == "AUDIT_SELECTION_STALE"
    p.field = SkyField(155.438012, 18.088025, 4)
    p.selected_dataset = ("dr11", "south", "1554p180", "g")
    assert inspect_saved_dataset_audit(p).state == "AUDIT_SELECTION_STALE"


def test_failed_report_is_historical_not_current_error(tmp_path):
    p = _fixture(tmp_path)
    _audit(p, record_status="INVALID_FITS")
    result = inspect_saved_dataset_audit(p)
    assert result.state == "AUDIT_REPORTED_ISSUE"
    assert "historical" in result.explanation.lower()
    assert "current-data failure" in result.explanation


def test_missing_receipt_and_external_source_are_not_adopted(tmp_path):
    p = _fixture(tmp_path)
    _, path = _audit(p)
    path.unlink()
    assert inspect_saved_dataset_audit(p).state == "AUDIT_EVIDENCE_UNAVAILABLE"


def test_project_data_root_can_point_to_parent_of_data_folder(tmp_path):
    p = _fixture(tmp_path)
    _audit(p)
    p.data_root = p.root
    assert inspect_saved_dataset_audit(p).state == "AUDIT_RECORDED_STAT_MATCH"


def test_receipt_symlink_fails_closed(tmp_path):
    p = _fixture(tmp_path)
    _, path = _audit(p)
    duplicate = path.parent / "elsewhere.json"
    duplicate.write_bytes(path.read_bytes())
    path.unlink()
    try:
        path.symlink_to(duplicate)
    except (OSError, NotImplementedError):
        pytest.skip("Symlink privileges unavailable")
    assert inspect_saved_dataset_audit(p).state == "AUDIT_EVIDENCE_UNAVAILABLE"
