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
    assert not w.files.isVisible()
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
