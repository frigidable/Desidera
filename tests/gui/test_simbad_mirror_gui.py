"""Offline Qt regression for SIMBAD mirror choice and provenance."""
import pytest

pytest.importorskip("pytestqt")
pytest.importorskip("PySide6.QtWidgets")

from lsbpipe.gui import galaxy_knowledge_dialog as gui
from lsbpipe.surveys.object_projects import create_object_project
from lsbpipe.surveys.legacy import SkyField


def test_mirror_policy_is_explicit_and_server_origin_remains_visible(
        qtbot, tmp_path, monkeypatch):
    project = create_object_project(
        "UGC 2162", tmp_path / "UGC_2162",
        field=SkyField(40.096208, 1.22925, 5))
    seen = []

    def fake_lookup(**kwargs):
        seen.append(kwargs)
        return {
            "schema": "desidera.knowledge.simbad-query/1",
            "source": "SIMBAD",
            "fetched_at_utc": "2026-11-03T14:00:00+00:00",
            "query": {"name": kwargs["name"], "radius_arcsec": 60},
            "candidates": [],
            "truncated": False,
            "endpoint_evidence": {
                "provider": "SIMBAD",
                "endpoint_host": "simbad.harvard.edu",
                "endpoint_policy": kwargs["mirror_policy"],
                "attempted_hosts": [
                    "simbad.cds.unistra.fr", "simbad.harvard.edu"],
                "failed_attempts": [{
                    "host": "simbad.cds.unistra.fr",
                    "error_type": "TimeoutError"}],
                "fallback_used": True,
            },
            "status": "CANDIDATES_UNVERIFIED",
            "science_status": "NOT_SCIENTIFICALLY_ADOPTED",
        }

    monkeypatch.setattr(gui, "query_simbad_candidates", fake_lookup)
    dialog = gui.GalaxyKnowledgeDialog(project.root)
    qtbot.addWidget(dialog)
    dialog.show()
    assert not seen  # opening project remains offline
    assert dialog.endpoint_policy.currentData() == "auto"
    dialog.search_btn.click()
    qtbot.waitUntil(lambda: dialog._future is None, timeout=5000)
    assert seen[0]["mirror_policy"] == "auto"
    assert "simbad.harvard.edu" in dialog.endpoint_note.text()
    assert "fallback" in dialog.endpoint_note.text()
    assert "simbad.cds.unistra.fr" in dialog.endpoint_note.text()

    reopened = gui.GalaxyKnowledgeDialog(project.root)
    qtbot.addWidget(reopened)
    assert "simbad.harvard.edu" in reopened.endpoint_note.text()
    assert len(seen) == 1
