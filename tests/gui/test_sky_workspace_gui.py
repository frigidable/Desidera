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
            "status":"EXPLICIT_DISJOINT_SKY_REGIONS_REQUIRED",
            "science_status":"UNADOPTED_DIAGNOSTIC_ONLY",
            "project_id":dataset["project_id"],"dataset_source":dataset["source"],
            "band":dataset["band"],"draft_sha256":draft.digest,
            "catalog_sha256":contamination["catalog_sha256"],
            "sky_train_effective_pixels":0,
            "sky_holdout_effective_pixels":0,
            "methods":{},"method_comparison":None,
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
    assert "EXPLICIT_DISJOINT_SKY_REGIONS_REQUIRED" in (
        QApplication.clipboard().text())
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
