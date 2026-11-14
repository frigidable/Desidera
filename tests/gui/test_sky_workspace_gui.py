"""Stage 03 Qt smoke tests: project scope, explicit execution, stale QA."""
from pathlib import Path
from types import SimpleNamespace

import pytest

pytest.importorskip("pytestqt")
pytest.importorskip("PySide6")

from PySide6.QtWidgets import QApplication, QGroupBox
from lsbpipe.gui.photometry_workspace import PhotometryWorkspace
from lsbpipe.gui.sky_workspace import SkyWorkspace
import lsbpipe.gui.sky_workspace as sky_module


def test_stage03_is_real_and_resumable_without_loading_or_running(qtbot):
    w = PhotometryWorkspace()
    qtbot.addWidget(w)
    assert w.stage_tabs.count() == 4
    assert w.stage_tabs.tabText(2) == "3 · Background"
    assert w.stage_tabs.tabText(3) == "4 · Aperture Photometry"
    assert w.sky_stage._future is None
    assert w.sky_stage._result is None
    assert not w.sky_stage.run_btn.isEnabled()
    w.restore_navigation(stage="sky", band="r")
    assert w.stage_tabs.currentWidget() is w.sky_stage
    assert w.current_navigation() == ("sky", "r")
    assert w.scene_editor._loaded is None
    assert w.sky_stage._future is None
    assert "current" in w.sky_stage.context_label.text().lower()
    assert w.sky_stage.findChild(QGroupBox).title() == (
        "Use current Scene and Masks evidence")
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



def test_stage03_snapshot_refuses_moved_target_or_foreign_brick(qtbot):
    """Do not pair old original pixels with new centre or Tractor masks."""
    identity=("dr11","south","1554p180","r")
    pid="p"
    ds={"project_id":pid,"band":"r","source":"DR11 south · 1554p180 · r",
        "source_quartet_key":list(identity),
        "loaded_target_icrs":[155.4,18.]}
    contam={
        "status":"CATALOGUE_CONTAMINATION_MASK_NOT_ADOPTED",
        "project_id":pid,"band":"r","dataset_source":ds["source"],
        "source_quartet_key":list(identity),
        "draft_digest":"a"*64,"catalog_sha256":"b"*64}
    scene=SimpleNamespace(
        _band="r",_future=None,_project_root=Path("/not/read"),
        _project_id=pid,_project_target=(155.4,18.),
        _source_key=("root",pid,identity,"root",None),
        _loaded=ds,
        _draft=SimpleNamespace(project_id=pid,digest="a"*64),
        _contamination=contam)
    w=SkyWorkspace(scene)
    qtbot.addWidget(w)
    assert w._snapshot() is not None
    assert w.run_btn.isEnabled()
    scene._project_target=(155.401,18.)
    assert w._snapshot() is None
    scene._project_target=(155.4,18.)
    contam["source_quartet_key"]=["dr11","south","1555p180","r"]
    assert w._snapshot() is None
    contam["source_quartet_key"]=None
    assert w._snapshot() is None
    contam["source_quartet_key"]=list(identity)
    scene._source_key=("root",pid,("dr11","south","1555p180","r"),
                       "root",None)
    assert w._snapshot() is None
    scene._source_key=("root",pid,identity,"root",None)
    assert w._snapshot() is not None


def test_stage03_readout_exposes_synthetic_sky_leakage_not_only_recovery():
    result = {
        "dataset_source": "DR11 south sample r",
        "draft_sha256": "a"*64,
        "sky_train_effective_pixels": 120,
        "sky_holdout_effective_pixels": 64,
        "methods": {},
        "method_comparison": None,
        "sky_injection_recovery_qa": {
            "status": "SAME_COADD_FROZEN_MASK_INJECTION_DIAGNOSTIC_ONLY",
            "scenarios": [{
                "scale_length_as_protected_area_radius": 1.0,
                "scale_length_arcsec": 60.,
                "unmasked_aperture_pixels": 700,
                "geometric_aperture_pixels": 1000,
                "unmasked_aperture_fraction": .7,
                "models": {
                    "constant": {
                        "injection_flux_recovery_fraction_fixed_aperture": .8,
                    }
                },
                "injected_signal_on_frozen_sky_samples": {
                    "training_median_nMgy_per_pixel": .0004,
                    "holdout_median_nMgy_per_pixel": .0005,
                },
            }],
        },
    }
    displayed = SkyWorkspace._format_readout(result)
    assert "h=60.00 arcsec" in displayed
    assert "700/1000 pixels (70.0%)" in displayed
    assert "training median=0.0004" in displayed
    assert "holdout median=0.0005" in displayed
    assert "NOT total galaxy light" in displayed
    assert "signal leakage, NOT a sky pedestal" in displayed
    assert "SKY SCIENCE HOLD" in displayed


def test_sky_summary_is_compact_and_full_science_report_on_demand(qtbot):
    from types import SimpleNamespace
    scene=SimpleNamespace(
        _band="r",_future=None,_project_root=Path("/not/read"),
        _project_id="p",_loaded=None,_draft=None,_contamination=None)
    w=SkyWorkspace(scene)
    qtbot.addWidget(w)
    w.resize(1100,700)
    w.show()
    assert not w.run_btn.isEnabled()
    assert not w.report_toggle.isVisible()
    assert not w.readout.isVisible()
    assert w.status.isVisible()
    assert w.status.text() == "Not evaluated."
    from PySide6.QtWidgets import QLabel
    assert any("science HOLD" in label.text()
               for label in w.findChildren(QLabel))
    assert "Scene" in w.context_label.text()
    # A numeric readout remains stored (and can be exported) but is not
    # a permanently open wall of text in the initial workspace.
    assert "No background result" in w.readout.toPlainText()
    w._result={"status":"INSUFFICIENT_AUTOMATIC_SCENE_SKY_GEOMETRY"}
    w.report_toggle.show()
    w.report_toggle.click()
    assert w.readout.isVisible()
    w._reset_report_visibility()
    assert not w.report_toggle.isVisible()
    assert not w.readout.isVisible()


def test_sky_empty_state_has_one_active_navigation_action(qtbot):
    w=PhotometryWorkspace()
    qtbot.addWidget(w)
    w.resize(1366,768)
    w.show()
    w.restore_navigation(stage="sky",band="g")
    assert w.sky_stage.open_scene_btn.isVisible()
    assert w.sky_stage.open_scene_btn.isEnabled()
    assert not w.sky_stage.run_btn.isVisible()
    assert not w.sky_stage.copy_btn.isVisible()
    assert not w.sky_stage.report_toggle.isVisible()
    w.sky_stage.open_scene_btn.click()
    assert w.stage_tabs.currentWidget() is w.scene_editor
    assert w.scene_editor._loaded is None
