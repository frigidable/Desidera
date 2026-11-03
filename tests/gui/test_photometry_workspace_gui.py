"""Offline Qt transition: Projects -> Workspace -> Photometry screen."""
from __future__ import annotations

import pytest

pytest.importorskip("pytestqt")
pytest.importorskip("PySide6.QtWebEngineWidgets")

from lsbpipe.gui.survey_explorer import SurveyExplorer
from lsbpipe.surveys.object_projects import create_object_project
from lsbpipe.surveys.legacy import SkyField
from scripts.offline_gui_fixture import offline_legacy_viewer


@pytest.fixture
def workspace(qtbot, tmp_path, monkeypatch):
    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path / "session"))
    project = create_object_project(
        "Leo P", tmp_path / "Leo_P",
        field=SkyField(155.438012, 18.088025, 5.))
    with offline_legacy_viewer():
        window = SurveyExplorer(restore_session=False, launch_remote=False)
    qtbot.addWidget(window)
    window.attach_object_project(project)
    window.show()
    yield window, project
    window._returning_to_hub = True
    window.close()


def test_screen_is_isolated_until_open_and_returns_to_workspace(workspace):
    window, project = workspace
    assert window.photometry_btn.isVisible()
    assert window.preview_tabs.count() == 2
    assert window.preview_tabs.indexOf(window.photometry_workspace) == -1
    assert not window.photometry_workspace.isVisible()
    window.photometry_btn.click()
    assert window.preview_tabs.currentWidget() is window.photometry_workspace
    assert window.preview_tabs.count() == 3
    assert window.workflow_stage.text() == "03 / 04 · Photometry"
    assert "No dataset selected" in window.photometry_workspace.source_identity.text()
    assert "NO_SELECTED_DATASET" in window.photometry_workspace.source_state.text()
    assert "NOT_MEASURED_NOT_ADOPTED" in window.photometry_workspace.blocker.text()
    assert not window.band_compare_btn.isVisible()
    assert not window.comparison_btn.isVisible()
    assert not window.data_toggle.isEnabled()
    assert not (project.root / "photometry-runs").exists()
    window.photometry_workspace.back_btn.click()
    assert window.preview_tabs.currentWidget() is window.sky
    assert window.data_toggle.isEnabled()


def test_manage_sources_from_photometry_opens_real_dataset_manager(workspace):
    window, project = workspace
    window.photometry_btn.click()
    window.photometry_workspace.sources_btn.click()
    assert window.preview_tabs.currentWidget() is window.dataset_manager
    assert window.dataset_manager.tree.topLevelItemCount() == 0
    assert not (project.root / "datasets" / "index.json").exists()
