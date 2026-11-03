"""Consolidation gate: compact passport retains science and knowledge actions."""
import pytest

pytest.importorskip("pytestqt")
pytest.importorskip("PySide6.QtWebEngineWidgets")

from lsbpipe.gui.project_hub import ProjectHub
from lsbpipe.surveys.object_projects import create_object_project
from lsbpipe.surveys.legacy import SkyField
from scripts.offline_gui_fixture import offline_legacy_viewer


def test_compact_passport_preserves_simbad_and_immutable_setup_controls(
        qtbot, tmp_path, monkeypatch):
    original = create_object_project(
        "UGC 2162", tmp_path / "UGC_2162",
        field=SkyField(40.096208, 1.22925, 5.))
    manifest = original.root / "lsbpipe-object.json"
    unchanged = manifest.read_bytes()
    with offline_legacy_viewer():
        hub = ProjectHub()
    qtbot.addWidget(hub)
    hub.show()
    assert hub.galaxy_card.snapshot.project_id == original.project_id
    assert hub.galaxy_card.open_button.isEnabled()
    assert hub.catalogues_menu_btn.isVisible()
    assert hub.simbad_action.isEnabled()
    assert hub.tools_menu_btn.isVisible()
    assert hub.config_action.isEnabled()
    assert not hub.continue_btn.isVisible()
    assert not hub.remove_btn.isVisible()
    assert hub.remove_action.isEnabled()
    assert manifest.read_bytes() == unchanged
    hub.filter.setText("unknown-no-match")
    assert hub.galaxy_card.snapshot is None
    assert not hub.simbad_action.isEnabled()
    assert not hub.config_action.isEnabled()
    hub.filter.clear()
    assert hub.galaxy_card.snapshot.project_id == original.project_id
