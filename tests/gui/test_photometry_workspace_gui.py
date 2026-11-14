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



def test_photometry_initial_status_does_not_repeat_copy_qa_instructions(qtbot):
    from lsbpipe.gui.photometry_workspace import PhotometryWorkspace
    w = PhotometryWorkspace()
    qtbot.addWidget(w)
    w.show()
    assert "no adopted galaxy photometry" in w.subtitle.text()
    assert not w.project_qa_status.isVisible()
    assert not w.batch_status.isVisible() or not w.batch_status.text()
    w.copy_project_qa_btn.click()
    assert w.project_qa_status.isVisible()
    assert "Open a saved galaxy project" in w.project_qa_status.text()



def test_empty_photometry_does_not_show_disabled_multiband_toolbar(qtbot):
    from lsbpipe.gui.photometry_workspace import PhotometryWorkspace
    w = PhotometryWorkspace()
    qtbot.addWidget(w)
    w.resize(1366, 768)
    w.show()
    assert not w.batch_box.isVisible()
    # Hidden until real inputs exist, but project QA and stage navigation remain.
    assert w.stage_tabs.count() == 4
    assert w.copy_project_qa_btn.isVisible()
    assert w.sources_btn.isVisible()


def test_advanced_photometry_preserves_multiband_controls_without_fits(qtbot):
    from lsbpipe.gui.photometry_workspace import PhotometryWorkspace
    w = PhotometryWorkspace()
    qtbot.addWidget(w)
    w.resize(1366, 768)
    w.show()
    assert not w.batch_box.isVisible()
    w.view_mode.setCurrentIndex(w.view_mode.findData("research"))
    assert w.batch_box.isVisible()
    assert w.save_full_qa_btn.isVisible()
    assert w.copy_full_qa_btn.isVisible()
    w.view_mode.setCurrentIndex(w.view_mode.findData("guided"))
    assert not w.batch_box.isVisible()



def test_photometry_stage_shortcuts_are_scoped_and_side_effect_free(qtbot):
    from PySide6.QtCore import Qt
    from PySide6.QtGui import QShortcut
    from lsbpipe.gui.photometry_workspace import PhotometryWorkspace

    w=PhotometryWorkspace()
    qtbot.addWidget(w)
    w.resize(1000,700)
    w.show()
    assert w.stage_tabs.count()==4
    assert len(w._stage_shortcuts)==4
    assert "Ctrl+2" in w.stage_tabs.tabToolTip(1)
    assert "Ctrl+4" in w.stage_tabs.tabToolTip(3)
    assert w.stage_tabs.tabBar().accessibleName().startswith(
        "Photometry stages; Ctrl+1")
    for key,expected in enumerate(range(4),start=1):
        shortcut=w._stage_shortcuts[key]
        assert isinstance(shortcut,QShortcut)
        assert shortcut.key().toString()==f"Ctrl+{key}"
        assert shortcut.context()==Qt.WidgetWithChildrenShortcut
        shortcut.activated.emit()
        assert w.stage_tabs.currentIndex()==expected
        assert w.scene_editor._loaded is None
        assert w.scene_editor._future is None
        assert w.sky_stage._future is None
        assert w.aperture_stage._future is None
    assert w.current_navigation()==("aperture","g")
    w._stage_shortcuts[1].activated.emit()
    assert w.current_navigation()==("psf0","g")


def test_psf_qa_export_suggests_bound_galaxy_not_leo_p(
        qtbot, tmp_path, monkeypatch):
    """A different galaxy must never inherit Leo P's debug export name."""
    from pathlib import Path
    from lsbpipe.gui import photometry_workspace as module

    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path / "session"))
    project = create_object_project(
        "UGC 2162 (candidate)", tmp_path / "UGC_2162",
        field=SkyField(40.096208, 1.229250, 5.))
    widget = module.PhotometryWorkspace()
    qtbot.addWidget(widget)
    widget.bind_project(project)
    widget.show()
    dialogs = []

    def cancel_dialog(*args):
        dialogs.append(args)
        return "", ""

    monkeypatch.setattr(module.QFileDialog, "getSaveFileName", cancel_dialog)
    widget.save_full_qa()
    assert len(dialogs) == 1
    proposed = Path(dialogs[0][2])
    assert proposed.parent == project.root
    assert proposed.name.startswith("Desidera_UGC_2162_candidate_PSF0_griz_debug_")
    assert proposed.suffix == ".json"
    assert "LeoP" not in proposed.name
    assert not list(project.root.glob("Desidera_*.json"))


def test_unbound_psf_feedback_is_galaxy_neutral(qtbot):
    from lsbpipe.gui.photometry_workspace import PhotometryWorkspace
    widget = PhotometryWorkspace()
    qtbot.addWidget(widget)
    widget._request_all_psfs()
    assert "galaxy project" in widget.batch_status.text()
    assert "Leo P" not in widget.batch_status.text()
    widget.save_full_qa()
    assert "galaxy project" in widget.batch_status.text()
    assert "Leo P" not in widget.batch_status.text()


def test_psf_dispatch_rejects_foreign_project_without_leo_p_message(
        qtbot, tmp_path, monkeypatch):
    """A foreign project request fails closed without misleading target name."""
    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path / "session"))
    project = create_object_project(
        "UGC 2162", tmp_path / "UGC_2162",
        field=SkyField(40.096208, 1.229250, 5.))
    with offline_legacy_viewer():
        window = SurveyExplorer(restore_session=False, launch_remote=False)
    qtbot.addWidget(window)
    window.attach_object_project(project)
    panel = window.photometry_workspace.psf0_panels["r"]
    window._run_project_psf0({
        "project_root": str(tmp_path / "foreign"), "band": "r"})
    assert "Galaxy project identity changed" in panel.qa.toPlainText()
    assert "Leo P" not in panel.qa.toPlainText()
    assert window.thread is None
    assert "current galaxy project" in panel.get_tractor_btn.toolTip()
    assert "selected galaxy" in panel.protection.toolTip()
    window._returning_to_hub = True
    window.close()


def test_photometry_stage_icons_and_theme_do_not_start_work(qtbot):
    from lsbpipe.gui.photometry_workspace import PhotometryWorkspace
    from lsbpipe.gui.visual_theme import theme_colors
    widget=PhotometryWorkspace()
    qtbot.addWidget(widget)
    assert all(not widget.stage_tabs.tabIcon(i).isNull()
               for i in range(4))
    assert all(not button.icon().isNull() for button in (
        widget.back_btn,widget.sources_btn,
        widget.refresh_btn,widget.copy_project_qa_btn))
    assert theme_colors(widget)["ACCENT"] in widget.styleSheet()
    assert widget.stage_tabs.tabBar().accessibleName().startswith(
        "Photometry stages")
    assert widget.scene_editor._loaded is None
    assert widget.scene_editor._future is None
    assert widget.sky_stage._future is None
    assert widget.aperture_stage._future is None


def test_observatory_vector_icons_render_full_shapes_not_cropped(qtbot):
    """HiDPI QPixmap rendering must not double-scale and clip icons."""
    from PySide6.QtCore import QSize
    from lsbpipe.gui.visual_theme import icon_for
    for kind in ("galaxy","folder","plus","back","database","chart",
                 "download","search","refresh","copy","mask","sky",
                 "aperture","open","save"):
        icon=icon_for(kind,color="#215CAA",size=24)
        assert not icon.isNull(),kind
        pix=icon.pixmap(QSize(24,24))
        img=pix.toImage()
        points=[(x,y) for y in range(img.height())
                for x in range(img.width())
                if img.pixelColor(x,y).alpha()>24]
        assert points,kind
        xs=[x for x,_ in points]
        ys=[y for _,y in points]
        assert max(xs)-min(xs)>img.width()*.40,kind
        assert max(ys)-min(ys)>img.height()*.40,kind
