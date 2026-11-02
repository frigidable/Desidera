"""UX sequence and disclosure regressions (no network or original FITS writes)."""
from __future__ import annotations
import pytest

pytest.importorskip("pytestqt")
pytest.importorskip("PySide6.QtWebEngineWidgets")

from scripts.offline_gui_fixture import offline_legacy_viewer, apply_offline_state


def test_step_labels_follow_conservative_workflow_contract():
    from lsbpipe.gui.workflow_ux import workflow_cue
    def stage(**kwargs):
        args = dict(has_target=True, has_catalog=False,
                    selected_brick=False, has_local_fits=False,
                    has_interactive_image=False)
        args.update(kwargs)
        return workflow_cue(**args)
    assert stage(has_target=False).step == 1
    assert stage().step == 1
    assert stage(has_catalog=True).step == 2
    assert stage(has_catalog=True, selected_brick=True).step == 2
    assert stage(has_catalog=True, selected_brick=True,
                 has_local_fits=True).step == 3
    assert "Open" in stage(
        has_catalog=True, selected_brick=True,
        has_local_fits=True).instruction
    assert stage(has_catalog=True, selected_brick=True,
                 has_local_fits=True,
                 has_interactive_image=True).step == 3


def test_clean_initial_survey_hides_unneeded_elements(qtbot):
    from lsbpipe.gui.survey_explorer import SurveyExplorer
    with offline_legacy_viewer():
        w = SurveyExplorer(restore_session=False, launch_remote=False)
    qtbot.addWidget(w)
    w.show()
    assert w.workflow_stage.text().startswith("01 / 03")
    assert w.search_btn.isVisible()
    assert w.center_btn.isVisible()
    assert not w.advanced_target.isVisible()
    for dr in ("dr9", "dr10", "dr11"):
        assert w.release_boxes[dr].isVisible()
    assert not w.project_status.isVisible()
    assert not w.bricks.isVisible()
    assert not w.data_drawer.isVisible()
    assert not w.progress_widget.isVisible()
    assert not w.comparison_btn.isVisible()
    w.target_options.click()
    assert w.advanced_target.isVisible()
    assert w.paste_target_btn.isVisible()
    assert w.release_boxes["dr10"].isVisible()
    assert w.view_mode.isVisible()
    w.target_options.click()
    assert not w.advanced_target.isVisible()


def test_discovery_moves_to_data_step_without_auto_downloading(qtbot):
    from lsbpipe.gui.survey_explorer import SurveyExplorer
    with offline_legacy_viewer():
        w = SurveyExplorer(restore_session=False, launch_remote=False)
    qtbot.addWidget(w)
    w.show()
    apply_offline_state(w, "S1")
    assert w.workflow_stage.text().startswith("02 / 03")
    assert w.bricks.isVisible()
    assert w.files.isVisible()
    assert w.files.rowCount() == len(w.manifest["products"])
    assert w.data_drawer.isVisible()
    assert not w.advanced_files.isVisible()
    assert len(w._selected()) == 0
    assert not w.download_btn.isEnabled()
    from PySide6.QtCore import QPoint
    w.resize(1366, 768)
    qtbot.wait(30)
    y = w.selected_label.mapTo(w.data_scroll.viewport(), QPoint()).y()
    assert 0 <= y < 210, "Core archive actions must be top-aligned"
    w.files_options.click()
    assert w.advanced_files.isVisible()
    assert w.files.isVisible()
    assert w.profile.isVisible()
    assert w.verify_btn.isVisible()
    assert w.qa_toggle.isVisible()
    w.files_options.click()
    assert not w.advanced_files.isVisible()


def test_science_panel_has_single_visible_primary_action(qtbot):
    from lsbpipe.gui.qa_workbench import ScienceQAControls
    qa = ScienceQAControls()
    qtbot.addWidget(qa)
    qa.resize(380, 650)
    qa.show()
    assert qa.dataset_audit.isVisible()
    assert qa.buttons["explore"].isVisible()
    assert not qa.steps.isVisible()
    assert not qa.report_panel.isVisible()
    qa.advanced_toggle.click()
    assert qa.steps.isVisible()
    assert qa.buttons["square"].isVisible()
    qa.advanced_toggle.click()
    assert not qa.steps.isVisible()
    qa.result_toggle.click()
    assert qa.report_panel.isVisible()
    qa.result_toggle.click()
    assert not qa.report_panel.isVisible()


def test_busy_indicator_only_visible_while_worker_runs(qtbot):
    from lsbpipe.gui.survey_explorer import SurveyExplorer
    with offline_legacy_viewer():
        w = SurveyExplorer(restore_session=False, launch_remote=False)
    qtbot.addWidget(w)
    w.show()
    result = []
    assert not w.progress_widget.isVisible()
    w._start(lambda progress: "safe", result.append)
    assert w.progress_widget.isVisible()
    qtbot.waitUntil(lambda: w.thread is None and result == ["safe"],
                    timeout=15000)
    assert not w.progress_widget.isVisible()
    assert not w.bar.isVisible()

def test_release_toggle_invalidates_catalog_without_touching_project_files(qtbot):
    from lsbpipe.gui.survey_explorer import SurveyExplorer
    with offline_legacy_viewer():
        window = SurveyExplorer(restore_session=False, launch_remote=False)
    qtbot.addWidget(window)
    window.show()
    apply_offline_state(window, "S3")
    assert window.manifest is not None
    assert len(window._selected()) == 4
    assert window.files.isVisible()
    window.release_boxes["dr10"].setChecked(True)
    assert window.manifest is None
    assert window.hits == []
    assert window.bricks.rowCount() == 0
    assert window.files.rowCount() == 0
    assert not window.files.isVisible()
    assert not window.download_btn.isEnabled() or window._selected() == set()
    assert "Find data" in window.message.text()


def test_desidera_launch_compatibility_aliases():
    from desidera.cli import main as new_cli
    from desidera.gui.__main__ import main as new_hub
    from lsbpipe.cli import main as old_cli
    from lsbpipe.gui.project_hub import main as old_hub
    assert new_cli is old_cli
    assert new_hub is old_hub

def test_same_brick_across_releases_has_separate_candidate_rows(qtbot):
    from dataclasses import replace
    from lsbpipe.gui.survey_explorer import SurveyExplorer
    from lsbpipe.surveys.legacy import SkyField, build_manifest
    from scripts.offline_gui_fixture import make_offline_catalog
    with offline_legacy_viewer():
        w = SurveyExplorer(restore_session=False, launch_remote=False)
    qtbot.addWidget(w)
    w.show()
    hits, _ = make_offline_catalog()
    dr9 = hits[0]
    dr11 = replace(dr9, release="dr11")
    all_hits = [dr9, dr11]
    field = SkyField(w.ra.value(), w.dec.value(), w.field_width.value())
    manifest = build_manifest(field, all_hits)
    w._receive_search((all_hits, manifest))
    assert w.files.isVisible()
    assert len(manifest["products"]) > 4
    w.bricks.selectRow(0)
    assert tuple(w.brick_filter.currentData()) == (
        dr9.release, dr9.region, dr9.name)
    shown = [manifest["products"][i] for i in range(w.files.rowCount())
             if not w.files.isRowHidden(i)]
    assert shown and {e["release"] for e in shown} == {"dr9"}
    w.bricks.selectRow(1)
    assert tuple(w.brick_filter.currentData()) == (
        dr11.release, dr11.region, dr11.name)
    shown = [manifest["products"][i] for i in range(w.files.rowCount())
             if not w.files.isRowHidden(i)]
    assert shown and {e["release"] for e in shown} == {"dr11"}
    assert not w._selected(), "DR selection must not trigger a download"

def test_remote_webengine_isolation_preserves_catalog_discovery(qtbot, monkeypatch):
    """No Chromium widget is created, yet science FITS selection remains."""
    from lsbpipe.gui.survey_explorer import SurveyExplorer, _NoWebEngineSky
    monkeypatch.setenv("DESIDERA_DISABLE_WEBENGINE", "1")
    with offline_legacy_viewer():
        w = SurveyExplorer(restore_session=False, launch_remote=False)
    qtbot.addWidget(w)
    w.show()
    assert isinstance(w.sky, _NoWebEngineSky)
    assert w.sky.url().isEmpty()
    apply_offline_state(w, "S1")
    assert w.manifest is not None
    assert len(w.manifest["products"]) == w.files.rowCount()
    assert w.files.isVisible()
    assert w.sky.url().isEmpty()
    w.bricks.selectRow(0)
    assert w.science_btn.isEnabled()
    assert not w.download_btn.isEnabled()

def test_optional_science_pages_do_not_paint_over_navigation(qtbot):
    """Regression for the reported top-left ghost FITS/Data Manager controls.

    Before their first addTab, the three constructed QWidget pages were
    visible, sized 100x30 at workspace (0,0) because their parent was the
    QMainWindow itself. Hidden pages must remain invisible until activated,
    including across removal/reopening.
    """
    from lsbpipe.gui.survey_explorer import SurveyExplorer
    from PySide6.QtCore import QPoint
    from PySide6.QtWidgets import QApplication
    with offline_legacy_viewer():
        w = SurveyExplorer(restore_session=False, launch_remote=False)
    qtbot.addWidget(w)
    w.show()
    w.back_to_projects.show()
    qtbot.wait(50)
    pages = (w.four_layer_comparison, w.multiband_comparison,
             w.dataset_manager)
    for page in pages:
        assert w.preview_tabs.indexOf(page) == -1
        assert page.parent() is w.preview_tabs
        assert not page.isVisible()
        assert not any(child.isVisible() for child
                       in page.findChildren(type(w.back_to_projects)))
    # The navigation button owns its own top-left hit-testing surface.
    hit = QApplication.widgetAt(w.mapToGlobal(QPoint(65, 35)))
    assert hit is w.back_to_projects

    tab = w.preview_tabs.addTab(w.multiband_comparison, "Compare optical bands")
    w.preview_tabs.setCurrentIndex(tab)
    assert w.multiband_comparison.isVisible()
    w._close_multiband()
    assert w.preview_tabs.indexOf(w.multiband_comparison) == -1
    assert not w.multiband_comparison.isVisible()

    tab = w.preview_tabs.addTab(w.four_layer_comparison, "Data quality")
    w.preview_tabs.setCurrentIndex(tab)
    assert w.four_layer_comparison.isVisible()
    w._close_four_comparison()
    assert w.preview_tabs.indexOf(w.four_layer_comparison) == -1
    assert not w.four_layer_comparison.isVisible()

    tab = w.preview_tabs.addTab(w.dataset_manager, "Data Manager")
    w.preview_tabs.setCurrentIndex(tab)
    assert w.dataset_manager.isVisible()
    w.preview_tabs.setCurrentWidget(w.sky)
    assert not w.dataset_manager.isVisible()
