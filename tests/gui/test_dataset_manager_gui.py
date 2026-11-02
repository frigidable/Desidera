"""Qt offline acceptance: project-bound Dataset Manager, no downloads."""
from __future__ import annotations

import pytest

pytest.importorskip("pytestqt")
pytest.importorskip("PySide6.QtWebEngineWidgets")

from lsbpipe.surveys.object_projects import create_object_project
from lsbpipe.surveys.dataset_registry import (
    register_external_fits, read_dataset_index)
from lsbpipe.surveys.legacy import SkyField
from lsbpipe.gui.dataset_manager import DatasetManager
from lsbpipe.gui.survey_explorer import SurveyExplorer
from scripts.offline_gui_fixture import offline_legacy_viewer


@pytest.fixture
def project(tmp_path, monkeypatch):
    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path / "session"))
    return create_object_project(
        "UGC 2162", tmp_path / "UGC_2162",
        field=SkyField(40.096208, 1.22925, 5.))


def test_workspace_navigation_to_dataset_manager_keeps_old_fits_tabs(
        qtbot, project):
    with offline_legacy_viewer():
        workspace = SurveyExplorer(restore_session=False, launch_remote=False)
    qtbot.addWidget(workspace)
    workspace.attach_object_project(project)
    workspace.show()
    assert workspace.dataset_manager_btn.isVisible()
    assert workspace.preview_tabs.count() == 2
    workspace.dataset_manager_btn.click()
    assert workspace.preview_tabs.currentWidget() is workspace.dataset_manager
    assert workspace.preview_tabs.count() == 3
    assert workspace.dataset_manager.tree.topLevelItemCount() == 0
    assert "No indexed sources" in workspace.dataset_manager.source_summary.text()
    assert "not inspected" in workspace.dataset_manager.map_status.text().lower()
    assert not workspace.dataset_manager.map_btn.isEnabled()
    workspace.dataset_manager.back_btn.click()
    assert workspace.preview_tabs.currentWidget() is workspace.sky
    assert workspace.preview_tabs.indexOf(workspace.fits_explorer) == 1
    assert not (project.root / "datasets" / "index.json").exists()
    workspace._returning_to_hub = True
    workspace.close()


def test_dataset_manager_shows_actual_external_source_and_hash_action(
        qtbot, project, tmp_path):
    image = tmp_path / "hsc_external_i.fits"
    image.write_bytes(b"USER ORIGINAL BYTES")
    register_external_fits(project.root, image, provider="HSC", band="i")
    widget = DatasetManager()
    qtbot.addWidget(widget)
    widget.bind_project(project, read_dataset_index(project.root))
    widget.resize(1200, 700)
    widget.show()
    assert widget.tree.topLevelItemCount() == 1
    parent = widget.tree.topLevelItem(0)
    assert "HSC" in parent.text(0)
    assert parent.text(1) == "i"
    assert parent.childCount() == 1
    assert not widget.hash_btn.isEnabled()
    asset = parent.child(0)
    widget.tree.setCurrentItem(asset)
    assert widget.hash_btn.isEnabled()
    requested = []
    widget.checksumRequested.connect(requested.append)
    widget.hash_btn.click()
    assert requested == [asset.data(0, 0x0100)[1]]  # Qt.UserRole
    assert "not yet verified" in widget.details.text()
    assert widget.map_btn.isEnabled()
    assert image.read_bytes() == b"USER ORIGINAL BYTES"


def test_manager_real_wcs_footprint_display_is_read_only(
        qtbot, project):
    widget = DatasetManager()
    qtbot.addWidget(widget)
    widget.bind_project(project, read_dataset_index(project.root))
    original = widget.scene.itemsBoundingRect()
    assert original.width() >= 300  # target search square, arcseconds
    widget.set_coverage({
        "status": "GEOMETRY_ONLY_NOT_ADOPTED",
        "footprints": [{
            "asset_id": "asset-1", "source_id": "source-1",
            "provider": "Synthetic", "release": None, "region": None,
            "brick": None, "band": "r",
            "points": [[-100, -100], [100, -100],
                       [100, 100], [-100, 100]],
        }],
        "rejected": [],
        "omitted_count": 0,
    })
    assert "1 footprints" in widget.map_status.text()
    assert "Geometry only" in widget.map_status.text()
    assert widget.scene.itemsBoundingRect().width() >= original.width()
