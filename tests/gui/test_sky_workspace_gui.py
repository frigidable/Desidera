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
    assert w.stage_tabs.tabText(2) == "3 · Sky"
    assert w.stage_tabs.tabText(3) == "4 · Apertures"
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


def test_stage03_user_discard_during_running_worker_is_safe_and_retryable(
        qtbot, monkeypatch):
    """A running Future finishes, but its result never reaches QA or cache."""
    from concurrent.futures import Future
    from types import SimpleNamespace

    scene=SimpleNamespace(
        _band="r",_loaded=None,_draft=None,_contamination=None)
    stage=SkyWorkspace(scene)
    qtbot.addWidget(stage)
    stage.show()
    dataset={"project_id":"p","source":"r-native","band":"r"}
    draft=SimpleNamespace(digest="a"*64)
    contamination={"catalog_sha256":"b"*64}
    identity=("p","r-native","r",draft.digest,"b"*64,
              id(dataset),id(contamination))
    snap=(identity,dataset,draft,contamination)
    monkeypatch.setattr(stage,"_snapshot",lambda:snap)
    stage.refresh_context()
    future=Future()
    assert future.set_running_or_notify_cancel()
    monkeypatch.setattr(stage._executor,"submit",lambda *a,**kw:future)
    stage.run_diagnostic()
    assert stage._future is future
    assert stage.discard_btn.isVisible()
    assert not stage.run_btn.isEnabled()
    stage.discard_btn.click()
    assert not stage.discard_btn.isEnabled()
    assert "will finish" in stage.status.text()
    assert not stage.copy_btn.isEnabled()
    assert not stage.run_btn.isEnabled()
    future.set_result({"status":"SHOULD_NEVER_ENTER_QA"})
    stage._poll_worker()
    assert stage._result is None
    assert stage._result_identity is None
    assert stage._last_failure is None
    assert "discarded by user" in stage.status.text()
    assert stage.run_btn.isEnabled()
    assert not stage.discard_btn.isVisible()
    assert not stage.copy_btn.isVisible()
    assert "SHOULD_NEVER_ENTER_QA" not in stage.readout.toPlainText()
    assert not stage._timer.isActive()


def test_stage03_queued_diagnostic_is_cancelled_before_execution(qtbot,monkeypatch):
    from concurrent.futures import Future
    from types import SimpleNamespace
    scene=SimpleNamespace(
        _band="r",_loaded=None,_draft=None,_contamination=None)
    stage=SkyWorkspace(scene)
    qtbot.addWidget(stage)
    stage.show()
    dataset={"project_id":"p","source":"r-native","band":"r"}
    draft=SimpleNamespace(digest="a"*64)
    mask={"catalog_sha256":"b"*64}
    snap=(("p","r-native","r",draft.digest,"b"*64,
           id(dataset),id(mask)),dataset,draft,mask)
    monkeypatch.setattr(stage,"_snapshot",lambda:snap)
    stage.refresh_context()
    future=Future()
    monkeypatch.setattr(stage._executor,"submit",lambda *a,**kw:future)
    stage.run_diagnostic()
    stage.discard_btn.click()
    assert future.cancelled()
    assert stage._future is None
    assert stage.run_btn.isEnabled()
    assert "discarded by user" in stage.status.text()


def test_sky_worker_failure_is_actionable_and_can_be_retried_without_new_scene(
        qtbot,monkeypatch):
    """Failed diagnostic never consumes original Scene/mask evidence."""
    from concurrent.futures import Future
    from types import SimpleNamespace
    scene=SimpleNamespace(
        _band="r",_loaded=None,_draft=None,_contamination=None)
    stage=SkyWorkspace(scene)
    qtbot.addWidget(stage)
    stage.show()
    dataset={"project_id":"p","source":"native-r","band":"r"}
    draft=SimpleNamespace(digest="d"*64)
    mask={"catalog_sha256":"c"*64}
    identity=("p","native-r","r",draft.digest,mask["catalog_sha256"],
              id(dataset),id(mask))
    snap=(identity,dataset,draft,mask)
    monkeypatch.setattr(stage,"_snapshot",lambda:snap)
    stage.refresh_context()
    failed=Future()
    failed.set_exception(ValueError("test-only insufficient tiles"))
    retried=Future()
    tasks=[failed,retried]
    monkeypatch.setattr(stage._executor,"submit",lambda *args,**kw:tasks.pop(0))
    stage.run_diagnostic()
    stage._poll_worker()
    assert stage._result is None
    assert stage._last_failure["status"]=="WORKER_FAILED_CLOSED"
    assert "test-only insufficient tiles" in stage.status.text()
    assert "Evaluate residual sky to retry" in stage.status.text()
    assert stage.run_btn.isEnabled()
    assert not stage.copy_btn.isEnabled()
    assert not stage.report_toggle.isVisible()
    assert stage._snapshot() is snap
    stage.run_diagnostic()
    assert stage._future is retried
    assert stage._last_failure is None
    assert not stage.run_btn.isEnabled()
    stage.discard_pending_result()
    assert retried.cancelled()
    assert stage._future is None
    assert stage.run_btn.isEnabled()


def test_sky_executor_start_failure_is_reported_without_qt_exception(
        qtbot,monkeypatch):
    from types import SimpleNamespace
    scene=SimpleNamespace(_band="r",_loaded=None,_draft=None,
                          _contamination=None)
    stage=SkyWorkspace(scene)
    qtbot.addWidget(stage)
    ds={"project_id":"p","source":"native-r","band":"r"}
    draft=SimpleNamespace(digest="d"*64)
    mask={"catalog_sha256":"c"*64}
    snap=(("p","native-r","r",draft.digest,"c"*64,
           id(ds),id(mask)),ds,draft,mask)
    monkeypatch.setattr(stage,"_snapshot",lambda:snap)
    monkeypatch.setattr(
        stage._executor,"submit",
        lambda *args,**kw:(_ for _ in ()).throw(
            RuntimeError("test-only pool unavailable")))
    stage.refresh_context()
    stage.run_diagnostic()
    assert stage._future is None
    assert stage._active_identity is None
    assert stage.run_btn.isEnabled()
    assert stage._last_failure["status"]=="WORKER_START_FAILED_CLOSED"
    assert "test-only pool unavailable" in stage.status.text()
    assert "to retry" in stage.status.text()
    assert not stage.copy_btn.isEnabled()


def test_sky_copy_failure_report_redacts_exception_and_invalidates_on_rebind(
        qtbot,monkeypatch):
    import json
    from concurrent.futures import Future
    from types import SimpleNamespace
    from PySide6.QtWidgets import QApplication

    scene=SimpleNamespace(_band="r",_loaded=None,_draft=None,
                          _contamination=None)
    stage=SkyWorkspace(scene)
    qtbot.addWidget(stage)
    stage.show()
    ds={"project_id":"private-galaxy","source":"native-r","band":"r"}
    draft=SimpleNamespace(digest="d"*64)
    mask={"catalog_sha256":"c"*64}
    identity=(ds["project_id"],ds["source"],"r",draft.digest,
              mask["catalog_sha256"],id(ds),id(mask))
    snapshot=(identity,ds,draft,mask)
    monkeypatch.setattr(stage,"_snapshot",lambda:snapshot)
    stage.refresh_context()
    future=Future()
    future.set_exception(RuntimeError(
        "sensitive /home/researcher/secret/FITS and token=private"))
    monkeypatch.setattr(stage._executor,"submit",lambda *a,**kw:future)
    stage.run_diagnostic()
    stage._poll_worker()
    assert stage._last_failure is not None
    assert stage.copy_error_btn.isVisible()
    assert not stage.copy_btn.isVisible()
    stage.copy_error_btn.click()
    report=QApplication.clipboard().text()
    payload=json.loads(report)
    assert payload["operation"]=="evaluate_residual_sky"
    assert payload["error_type"]=="RuntimeError"
    assert payload["status"]=="WORKER_FAILED_CLOSED"
    assert payload["band"]=="r"
    assert "/home/researcher" not in report
    assert "token=private" not in report
    assert "private-galaxy" not in report
    assert "traceback" not in payload
    assert payload["science_status"]=="UNADOPTED_DIAGNOSTIC_ONLY"
    monkeypatch.setattr(stage,"_snapshot",lambda:None)
    stage.refresh_context()
    assert stage._last_failure is None
    assert not stage.copy_error_btn.isVisible()
    previous=QApplication.clipboard().text()
    stage.copy_failure_report()
    assert QApplication.clipboard().text()==previous


def test_failure_report_formatter_is_strict_allowlist():
    import json
    import pytest
    from lsbpipe.gui.workflow_ux import format_gui_failure_report
    report=json.loads(format_gui_failure_report(
        stage="background", error_type="Access /home/secret token",
        failure_status="WORKER_START_FAILED_CLOSED",
        band="../private/dir"))
    assert report["error_type"]=="Exception"
    assert report["band"] is None
    assert "secret" not in json.dumps(report)
    assert "background worker" in report["possible_cause"].lower()
    with pytest.raises(ValueError):
        format_gui_failure_report(
            stage="/home/private",error_type="ValueError",
            failure_status="WORKER_FAILED_CLOSED")
    with pytest.raises(ValueError):
        format_gui_failure_report(
            stage="scene",error_type="ValueError",
            failure_status="RAW_PRIVATE_TRACEBACK")


def test_sky_failure_keyboard_copy_shortcut_scoped_and_safe(qtbot,monkeypatch):
    from PySide6.QtCore import Qt
    from PySide6.QtGui import QFont
    from PySide6.QtWidgets import QApplication
    from concurrent.futures import Future
    from types import SimpleNamespace
    import json
    stage=SkyWorkspace(SimpleNamespace(
        _band="r",_loaded=None,_draft=None,_contamination=None))
    qtbot.addWidget(stage)
    stage.resize(760,580)
    stage.setFont(QFont("Arial",16))
    stage.show()
    data={"project_id":"p","source":"r","band":"r"}
    draft=SimpleNamespace(digest="d"*64)
    mask={"catalog_sha256":"c"*64}
    identity=("p","r","r",draft.digest,"c"*64,id(data),id(mask))
    monkeypatch.setattr(stage,"_snapshot",lambda:(identity,data,draft,mask))
    stage.refresh_context()
    failed=Future()
    failed.set_exception(ValueError("test-only failure"))
    monkeypatch.setattr(stage._executor,"submit",lambda *a,**kw:failed)
    stage.run_diagnostic()
    stage._poll_worker()
    qtbot.wait(20)
    assert stage.copy_error_btn.isVisible()
    assert stage.run_btn.isVisible()
    assert stage._failure_shortcut.context()==Qt.WidgetWithChildrenShortcut
    assert stage._failure_shortcut.key().toString()=="Ctrl+Shift+E"
    for button in (stage.run_btn,stage.copy_error_btn):
        point=button.mapTo(stage,button.rect().center())
        assert stage.rect().contains(point)
        assert button.width()>70
    QApplication.clipboard().setText("sentinel")
    stage._failure_shortcut.activated.emit()
    assert json.loads(QApplication.clipboard().text())["operation"]=="evaluate_residual_sky"
    monkeypatch.setattr(stage,"_snapshot",lambda:None)
    stage.refresh_context()
    QApplication.clipboard().setText("sentinel")
    stage._failure_shortcut.activated.emit()
    assert QApplication.clipboard().text()=="sentinel"


def test_sky_visual_error_tone_is_not_science_acceptance(qtbot,monkeypatch):
    from concurrent.futures import Future
    from types import SimpleNamespace
    stage=SkyWorkspace(SimpleNamespace(
        _band="r",_loaded=None,_draft=None,_contamination=None))
    qtbot.addWidget(stage)
    ds={"project_id":"p","source":"r","band":"r"}
    draft=SimpleNamespace(digest="d"*64)
    mask={"catalog_sha256":"c"*64}
    identity=("p","r","r",draft.digest,"c"*64,id(ds),id(mask))
    monkeypatch.setattr(stage,"_snapshot",lambda:(identity,ds,draft,mask))
    failed=Future()
    failed.set_exception(ValueError("test-only sky failure"))
    monkeypatch.setattr(stage._executor,"submit",lambda *a,**k:failed)
    stage.run_diagnostic()
    stage._poll_worker()
    assert stage.status.property("visualState")=="error"
    assert stage.run_btn.property("visualRole")=="primary"
    assert not stage.run_btn.icon().isNull()
    assert stage._result is None
    assert not stage.copy_btn.isEnabled()
    stage._snapshot=lambda:None
    stage.refresh_context()
    assert stage.status.property("visualState")=="normal"


def test_expanded_sky_source_identity_and_stale_probe(qtbot,monkeypatch):
    scene=SimpleNamespace(_band="r",_future=None,_project_root=Path("/not/read"),
        _project_id="p",_loaded={"project_id":"p","band":"r","source":"s",
            "grid_fingerprint":"grid","original_fits_sha256_local":{"image":"hash"}},
        _draft=SimpleNamespace(project_id="p",digest="d"),
        _contamination={"status":"CATALOGUE_CONTAMINATION_MASK_NOT_ADOPTED",
            "project_id":"p","band":"r","dataset_source":"s","draft_digest":"d",
            "catalog_sha256":"cat","effective_source_mask_sha256":"m"})
    w=SkyWorkspace(scene); qtbot.addWidget(w)
    assert not w.expanded_btn.isEnabled()
    scene._foreground_psf_result=dict(project_id="p",dataset_source="s",band="r",draft_digest="d",
        tractor_catalog_sha256="cat",original_scene_grid_fingerprint="grid",
        original_fits_sha256_local={"image":"hash"},sh1c_preview=dict(
            status="SH1C_CONDITIONAL_PREVIEW_NOT_ADOPTED",grid_fingerprint="grid",draft_digest="d",
            effective_source_mask_sha256="m",sky1_inputs={"image":"placeholder"}))
    w.refresh_context()
    assert w.expanded_btn.isEnabled()
    ident=w._snapshot()[0]+(id(scene._foreground_psf_result),)
    assert w._current_identity(ident,w._snapshot())==ident
    import numpy as np
    import json
    a=np.ones((64,64))
    def fake(ds,draft,cont,probe,other):
        return dict(status="SKY1_NOT_ADOPTED",project_id="p",dataset_source="s",band="r",
            draft_sha256="d",catalog_sha256="cat",foreground_source_objid=360,
            expanded_sky=dict(expanded_shape=[64,64],expanded_target_xy=[32.,32.],
                pixel_scale_arcsec=.262,warning="conditional",expanded_layers={"original":a,"SH1C":a},
                expanded_safe_mask=np.ones_like(a,dtype=bool),layers={},
                source_entries=[{"objid":5,"x":12.,"y":12.}],
                trials=[dict(halo_scenario="SH1C",target_guard_arcsec=60.,tile_size_pix=32,
                    candidates=[],blocks=[dict(tile_x_pix=0,tile_y_pix=0,x_pix=15.,y_pix=15.,
                        mean_nMgy_pix=.001,median_nMgy_pix=.0009,n=1024)])]))
    monkeypatch.setattr(sky_module,"run_expanded_sky",fake)
    w.run_expanded_diagnostic()
    qtbot.waitUntil(lambda:w._future is None,timeout=4000)
    assert w.current_expanded_result() is not None
    assert w.sky_graphics.items()
    w._inspect_block(12.5,51.5)
    assert "OBJID 5" in w.block_caption.text()
    assert "N=1024" in w.block_caption.text()
    w.copy_qa()
    qa=json.loads(QApplication.clipboard().text())
    assert qa["expanded_sky"]["expanded_layers"]["original"]["pixel_array_omitted"]
    scene._foreground_psf_result=dict(scene._foreground_psf_result)
    assert w._current_identity(ident,w._snapshot())!=ident
    scene._foreground_psf_result["draft_digest"]="stale"
    w.refresh_context()
    assert not w.expanded_btn.isEnabled()

    assert w.current_expanded_result() is None
    assert not w.copy_btn.isEnabled()
