"""P10h: saved FITS/WCS QA appears in Project Hub without science adoption."""
import hashlib
import json

import pytest

pytest.importorskip("pytestqt")
pytest.importorskip("PySide6.QtWebEngineWidgets")

from scripts.offline_gui_fixture import offline_legacy_viewer
from lsbpipe.gui.project_hub import ProjectHub
from lsbpipe.surveys.dataset_integrity import SCHEMA
from lsbpipe.surveys.legacy import BrickHit, SkyField, product_for
from lsbpipe.surveys.object_projects import (
    create_object_project, save_object_project, read_object_project,
)


def _make_project(tmp_path, *, with_receipt):
    p = create_object_project(
        "Leo P", tmp_path / "Leo_P",
        field=SkyField(155.438012, 18.088025, 4.))
    p.selected_dataset = ("dr11", "south", "1554p180", "r")
    p.data_root = p.root / "data"
    p.data_root.mkdir()
    hit = BrickHit("dr11", "south", "1554p180", 0, 0, 0, 0, {"r": None})
    products = {}
    for kind in ("image", "invvar", "nexp", "maskbits"):
        rel = product_for(
            hit, kind, "r" if kind != "maskbits" else None).relative_path
        file = p.data_root / rel
        file.parent.mkdir(parents=True, exist_ok=True)
        file.write_bytes(("source:" + kind).encode())
        st = file.stat()
        products[kind] = {
            "relative_path": rel, "size_bytes": st.st_size,
            "mtime_ns": st.st_mtime_ns,
            "sha256_local": hashlib.sha256(file.read_bytes()).hexdigest()}
    if with_receipt:
        stamp = "2026-11-03T12:00:00+00:00"
        report = {
            "schema": SCHEMA, "project_id": p.project_id,
            "created_utc": stamp, "data_root": str(p.data_root),
            "science_status": "DIAGNOSTIC_NOT_ADOPTED",
            "field": {"ra_deg": p.field.ra_deg,
                      "dec_deg": p.field.dec_deg,
                      "width_arcmin": p.field.width_arcmin},
            "datasets": [{
                "release": "dr11", "region": "south", "brick": "1554p180",
                "band": "r", "status": "FITS_WCS_VALIDATED_DIAGNOSTIC",
                "products": products,
                "fractions": {"instrumentally_usable": 0.9995},
            }],
        }
        folder = p.root / "dataset-audits" / "run123"
        folder.mkdir(parents=True)
        (folder / "audit.json").write_text(json.dumps(report), encoding="utf8")
        p.last_dataset_audit = {
            "relative_path": "dataset-audits/run123/audit.json",
            "created_utc": stamp, "data_root": str(p.data_root),
        }
    save_object_project(p)
    return p, products


def _hub(qtbot, monkeypatch):
    monkeypatch.setenv("DESIDERA_DISABLE_AUTOMATIC_RGB", "1")
    # Deliberately fake FITS, only their stat metadata is used by the passport.
    monkeypatch.setattr(
        "lsbpipe.gui.project_hub.has_local_image_source",
        lambda project: False)
    with offline_legacy_viewer():
        hub = ProjectHub()
    qtbot.addWidget(hub)
    hub.show()
    return hub


def test_no_saved_audit_is_not_qa_approval(qtbot, tmp_path, monkeypatch):
    p, _ = _make_project(tmp_path, with_receipt=False)
    manifest = (p.root / "lsbpipe-object.json").read_bytes()
    hub = _hub(qtbot, monkeypatch)
    card = hub.galaxy_card
    assert "4/4" in card.fits.text()
    assert "No saved FITS/WCS audit" in card.qa_status.text()
    assert "QA pending" in hub.projects.item(0).text()
    assert not card.next_step_title.isVisible()
    assert not card.next_step.isVisible()
    assert card.open_button.text() == "Open project"
    assert "inspect FITS" in card.open_button.toolTip()
    assert (p.root / "lsbpipe-object.json").read_bytes() == manifest


def test_saved_qa_visible_but_not_scientific_certification(
        qtbot, tmp_path, monkeypatch):
    p, products = _make_project(tmp_path, with_receipt=True)
    manifest = (p.root / "lsbpipe-object.json").read_bytes()
    original = {
        kind: (p.data_root / item["relative_path"]).read_bytes()
        for kind, item in products.items()
    }
    hub = _hub(qtbot, monkeypatch)
    card = hub.galaxy_card
    assert card.snapshot.audit.state == "AUDIT_RECORDED_STAT_MATCH"
    assert "metadata unchanged" in card.qa_status.text()
    assert "NOT rechecked" in card.qa_note.text()
    assert "NOT certified" in card.qa_note.text()
    assert "review saved FITS/WCS QA" in card.open_button.toolTip()
    assert card.open_button.text() == "Open project"
    assert not card.qa_note.isVisible()
    assert not card.observation_explanation.isVisible()
    assert "prior WCS QA" in hub.projects.item(0).text()
    assert not card.audit_receipt_note.isVisible()
    card.details_toggle.click()
    assert card.qa_note.isVisible()
    assert "NOT rechecked" in card.qa_note.text()
    assert card.observation_explanation.isVisible()
    assert "dataset-audits/run123/audit.json" in card.audit_receipt_note.text()
    assert (p.root / "lsbpipe-object.json").read_bytes() == manifest
    assert all((p.data_root / products[k]["relative_path"]).read_bytes()
               == original[k] for k in products)
    assert read_object_project(p.root).field == p.field


def test_changed_original_marks_prior_qa_stale_on_refresh(
        qtbot, tmp_path, monkeypatch):
    p, products = _make_project(tmp_path, with_receipt=True)
    hub = _hub(qtbot, monkeypatch)
    image = p.data_root / products["image"]["relative_path"]
    image.write_bytes(b"changed content and size")
    hub.refresh()
    card = hub.galaxy_card
    assert card.snapshot.audit.state == "AUDIT_INPUTS_CHANGED"
    assert "needs rerun" in card.qa_status.text()
    assert "QA needs rerun" in hub.projects.item(0).text()
    assert card.open_button.text() == "Open project"
    assert "rerun dataset QA" in card.open_button.toolTip()
