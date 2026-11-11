"""Stage 03 Qt smoke tests: project scope, explicit execution, stale QA."""
from pathlib import Path
from types import SimpleNamespace

import pytest

pytest.importorskip("pytestqt")
pytest.importorskip("PySide6")

from PySide6.QtWidgets import QApplication
from lsbpipe.gui.photometry_workspace import PhotometryWorkspace
from lsbpipe.gui.sky_workspace import SkyWorkspace
import lsbpipe.gui.sky_workspace as sky_module


def test_stage03_is_real_and_resumable_without_loading_or_running(qtbot):
    w = PhotometryWorkspace()
    qtbot.addWidget(w)
    assert w.stage_tabs.count() == 3
    assert w.stage_tabs.tabText(2) == "3 · Background"
    assert w.sky_stage._future is None
    assert w.sky_stage._result is None
    assert not w.sky_stage.run_btn.isEnabled()
    w.restore_navigation(stage="sky", band="r")
    assert w.stage_tabs.currentWidget() is w.sky_stage
    assert w.current_navigation() == ("sky", "r")
    assert w.scene_editor._loaded is None
    assert w.sky_stage._future is None
    assert "current" in w.sky_stage.context_label.text().lower()
    w.restore_navigation(stage="scene",band="i")
    assert w.current_navigation() == ("scene", "i")
    w.restore_navigation(stage="sky",band="i")
    assert w.current_navigation() == ("sky", "i")
    # Returning from Stage 03 must not overwrite the native Scene band
    # from a stale PSF tab selector.
    w.band_tabs.setCurrentIndex(0)
    w.stage_tabs.setCurrentWidget(w.scene_editor)
    assert w.current_navigation() == ("scene", "i")


def test_stage03_requires_current_scene_and_discards_stale_result(qtbot,monkeypatch):
    project="test-sky-project"
    source="dr11/south/example/r"
    scene=SimpleNamespace(
        _band="r", _future=None, _project_root=Path("/not/read"),
        _project_id=project,
        _loaded={"project_id":project,"band":"r","source":source},
        _draft=SimpleNamespace(project_id=project,digest="a"*64),
        _contamination={
            "status":"CATALOGUE_CONTAMINATION_MASK_NOT_ADOPTED",
            "project_id":project,"band":"r","dataset_source":source,
            "draft_digest":"a"*64,"catalog_sha256":"b"*64},
    )
    def fake(dataset,draft,contamination):
        return {
            "status":"INSUFFICIENT_AUTOMATIC_SCENE_SKY_GEOMETRY",
            "science_status":"UNADOPTED_DIAGNOSTIC_ONLY",
            "project_id":dataset["project_id"],"dataset_source":dataset["source"],
            "band":dataset["band"],"draft_sha256":draft.digest,
            "catalog_sha256":contamination["catalog_sha256"],
            "sky_train_effective_pixels":0,
            "sky_holdout_effective_pixels":0,
            "methods":{},"method_comparison":None,
            "automatic_sky_selection_qa":{
                "strategy":"SCENE_MASKS_FIXED_SPATIAL_TILES_V1",
                "training_tiles":0,"holdout_tiles":0,
                "sky_regions_drawn_by_user_required":False,
            },
            "sky_model_adopted":False,
            "original_fits_or_project_draft_modified":False,
        }
    monkeypatch.setattr(sky_module,"diagnose_current_project_residual_sky",fake)
    w=SkyWorkspace(scene)
    qtbot.addWidget(w)
    assert w.run_btn.isEnabled()
    assert w._result is None
    w.run_diagnostic()
    qtbot.waitUntil(lambda: w._future is None, timeout=4000)
    assert w._result is not None
    assert w._result["methods"] == {}
    assert w.copy_btn.isEnabled()
    w.copy_qa()
    assert "INSUFFICIENT_AUTOMATIC_SCENE_SKY_GEOMETRY" in (
        QApplication.clipboard().text())
    assert "sky_regions_drawn_by_user_required" in (
        QApplication.clipboard().text())
    assert "No manual sky regions required" in w.readout.toPlainText() or (
        "no manual sky regions required" in w.readout.toPlainText().lower())
    assert "No numerical sky fit performed" in w.readout.toPlainText()
    # Current background model cannot survive a new draft revision.
    scene._draft=SimpleNamespace(project_id=project,digest="c"*64)
    scene._contamination=dict(scene._contamination,draft_digest="c"*64)
    w.refresh_context()
    assert w._result is None
    assert not w.copy_btn.isEnabled()
    assert "invalidated" in w.readout.toPlainText()
    # No stale former result may be copied after the changed band.
    scene._band="g"
    w.copy_qa()
    assert w._result is None
