"""Offline Qt test: SIMBAD is opt-in and identities require confirmation."""
import pytest

pytest.importorskip("pytestqt")
pytest.importorskip("PySide6.QtWidgets")

from PySide6.QtWidgets import QMessageBox

from lsbpipe.surveys.legacy import SkyField
from lsbpipe.surveys.object_projects import create_object_project
from lsbpipe.surveys.galaxy_knowledge import latest_simbad_confirmation
from lsbpipe.gui import galaxy_knowledge_dialog as gui


def test_knowledge_dialog_offline_then_manual_query_and_confirm(
        qtbot, tmp_path, monkeypatch):
    project = create_object_project(
        "UGC 2162", tmp_path / "UGC_2162",
        field=SkyField(40.096208, 1.22925, 5.))
    manifest = (project.root / "lsbpipe-object.json").read_bytes()
    calls = []

    def fake_query(**kwargs):
        calls.append(kwargs)
        return {
            "schema": "desidera.knowledge.simbad-query/1",
            "source": "SIMBAD", "fetched_at_utc": "2026-11-03T12:00:00+00:00",
            "query": {"name": kwargs["name"], "radius_arcsec": 60},
            "candidates": [
                {"main_id": "UGC 2162", "ra_deg": 40.096208,
                 "dec_deg": 1.22925, "separation_arcsec": 0.,
                 "coordinate_bibcode": "", "matched_id": "UGC 2162",
                 "discovered_via": ["name", "cone"]},
                {"main_id": "Nearby object", "ra_deg": 40.097,
                 "dec_deg": 1.23, "separation_arcsec": 4.,
                 "coordinate_bibcode": "", "matched_id": "",
                 "discovered_via": ["cone"]}],
            "truncated": False, "status": "CANDIDATES_UNVERIFIED",
            "science_status": "NOT_SCIENTIFICALLY_ADOPTED",
        }

    monkeypatch.setattr(gui, "query_simbad_candidates", fake_query)
    dialog = gui.GalaxyKnowledgeDialog(project.root)
    qtbot.addWidget(dialog)
    dialog.show()
    assert not calls
    assert dialog.table.rowCount() == 0
    assert not dialog.confirm_btn.isEnabled()
    dialog.search_btn.click()
    qtbot.waitUntil(lambda: dialog._future is None, timeout=5000)
    assert len(calls) == 1
    assert dialog.table.rowCount() == 2
    assert not dialog.confirm_btn.isEnabled()
    dialog.table.selectRow(1)
    assert dialog.confirm_btn.isEnabled()
    monkeypatch.setattr(QMessageBox, "question", lambda *a, **k: QMessageBox.Yes)
    dialog.confirm_btn.click()
    saved = latest_simbad_confirmation(project.root)
    assert saved["candidate"]["main_id"] == "Nearby object"
    assert "Nearby object" in dialog.confirmed.text()
    assert saved["coordinate_adoption"] is False
    assert (project.root / "lsbpipe-object.json").read_bytes() == manifest

    # New dialogue opens offline from the immutable saved evidence.
    reopened = gui.GalaxyKnowledgeDialog(project.root)
    qtbot.addWidget(reopened)
    assert reopened.table.rowCount() == 2
    assert "Nearby object" in reopened.confirmed.text()
    assert len(calls) == 1
