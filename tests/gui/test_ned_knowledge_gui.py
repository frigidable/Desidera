"""P10b Qt acceptance: manual-only NED lookup and offline project passport."""
import pytest

pytest.importorskip("pytestqt")
pytest.importorskip("PySide6.QtWebEngineWidgets")

from PySide6.QtWidgets import QMessageBox

from scripts.offline_gui_fixture import offline_legacy_viewer
from lsbpipe.gui.project_hub import ProjectHub
from lsbpipe.gui import ned_knowledge_dialog as ned_gui
from lsbpipe.surveys.object_projects import create_object_project
from lsbpipe.surveys.legacy import SkyField


def test_ned_dialog_manual_only_then_confirm_and_offline_reopen(
        qtbot, tmp_path, monkeypatch):
    project = create_object_project(
        "Leo P", tmp_path / "Leo_P",
        field=SkyField(155.438012, 18.088025, 4.))
    manifest = project.root / "lsbpipe-object.json"
    before = manifest.read_bytes()
    calls = []

    def fake_remote(name):
        calls.append(name)
        return {
            "schema": "desidera.knowledge.ned-objectlookup/1",
            "provider": "NED",
            "api": "https://ned.ipac.caltech.edu/srs/ObjectLookup",
            "fetched_at_utc": "2026-11-03T15:00:00+00:00",
            "query_name": name, "api_version": "2.0", "result_code": 3,
            "status": "KNOWN_OBJECT",
            "candidate": {
                "main_id": "Leo P", "ra_deg": 155.4382,
                "dec_deg": 18.088, "object_type": "G",
                "position_refcode": "2021ApJ...000A..01B",
                "redshift": {"value": 0.0009, "uncertainty": 0.0001,
                             "refcode": "2021AJ...000...1A",
                             "quality_flag": "", "frame": "heliocentric",
                             "adopted": False}},
            "aliases_unadopted": ["Leo P", "AGC 208583"],
            "ambiguous_names": [],
            "source_response": {"Version": "2.0", "StatusCode": 100,
                                "ResultCode": 3},
            "science_status": "NOT_SCIENTIFICALLY_ADOPTED",
        }

    monkeypatch.setattr(ned_gui, "query_ned_identity", fake_remote)
    dialog = ned_gui.NedKnowledgeDialog(project.root)
    qtbot.addWidget(dialog)
    dialog.show()
    assert not calls
    assert not dialog.confirm_btn.isEnabled()
    dialog.search_btn.click()
    qtbot.waitUntil(lambda: dialog._future is None, timeout=5000)
    assert calls == ["Leo P"]
    assert "Object type: G" in dialog.response.text()
    assert dialog.confirm_btn.isEnabled()
    monkeypatch.setattr(QMessageBox, "question", lambda *a, **k: QMessageBox.Yes)
    dialog.confirm_btn.click()
    assert "Leo P" in dialog.confirmed.text()
    assert 'href="https://ned.ipac.caltech.edu/byname?objname=Leo+P"' in (
        dialog.confirmed.text())
    assert dialog.confirmed.openExternalLinks()

    reopened = ned_gui.NedKnowledgeDialog(project.root)
    qtbot.addWidget(reopened)
    assert reopened.confirm_btn.isEnabled()
    assert "Leo P" in reopened.confirmed.text()
    assert 'ned.ipac.caltech.edu/byname?objname=Leo+P' in reopened.confirmed.text()
    assert len(calls) == 1

    with offline_legacy_viewer():
        hub = ProjectHub()
    qtbot.addWidget(hub)
    hub.show()
    assert hub.catalogues_menu_btn.isVisible()
    assert hub.ned_action.isEnabled()
    assert "NED: " in hub.galaxy_card.catalogue.text()
    assert "Leo P" in hub.galaxy_card.catalogue.text()
    assert 'href="https://ned.ipac.caltech.edu/byname?objname=Leo+P"' in (
        hub.galaxy_card.catalogue.text())
    assert "unadopted" in hub.galaxy_card.catalogue_provenance.text()
    assert manifest.read_bytes() == before


def test_ambiguous_ned_does_not_enable_confirmation(qtbot, tmp_path, monkeypatch):
    project = create_object_project("Leo", tmp_path / "Leo")
    def ambiguous(name):
        return {
            "schema": "desidera.knowledge.ned-objectlookup/1",
            "provider": "NED",
            "api": "https://ned.ipac.caltech.edu/srs/ObjectLookup",
            "fetched_at_utc": "2026-11-03T15:00:00+00:00",
            "query_name": name, "api_version": "2.0", "result_code": 1,
            "status": "AMBIGUOUS_NAME", "candidate": None,
            "aliases_unadopted": [], "ambiguous_names": ["Leo I", "Leo II"],
            "source_response": {"Version": "2.0", "ResultCode": 1,
                                "StatusCode": 100},
            "science_status": "NOT_SCIENTIFICALLY_ADOPTED",
        }
    monkeypatch.setattr(ned_gui, "query_ned_identity", ambiguous)
    dialog = ned_gui.NedKnowledgeDialog(project.root)
    qtbot.addWidget(dialog)
    dialog.search_btn.click()
    qtbot.waitUntil(lambda: dialog._future is None, timeout=5000)
    assert "ambiguous" in dialog.response.text()
    assert not dialog.confirm_btn.isEnabled()
