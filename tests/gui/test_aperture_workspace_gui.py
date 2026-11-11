"""Stage 04 Qt smoke: no work on entry, source-bound worker and stale QA."""
from types import SimpleNamespace

import numpy as np
import pytest

pytest.importorskip("pytestqt")
pytest.importorskip("PySide6")
pytest.importorskip("photutils")
from PySide6.QtWidgets import QApplication

from lsbpipe.gui.photometry_workspace import PhotometryWorkspace
from lsbpipe.gui.aperture_workspace import (
    ApertureWorkspace, build_project_aperture_growth)
from lsbpipe.photometry.sky_masks import SkyMaskDraft
from test_aperture_flux import frame


def _inputs():
    f=frame(value=.02)
    uuid="stage04-synthetic-project"
    draft=SkyMaskDraft(uuid)
    dataset={
        "status":"NATIVE_SCENE_CUTOUT_NOT_SCIENCE_ADOPTED",
        "project_id":uuid,
        "band":"r","source":"DR11 south · 1554p180 · r",
        "pixel_shape":f.shape,"wcs":f.wcs,
        "native_hard":~f.valid,
        "data":{
            "image":f.image_nMgy,"invvar":f.invvar_nMgy2,
            "nexp":f.nexp,"maskbits":f.maskbits,
        },
    }
    yy,xx=np.indices(f.shape)
    protected=(xx-60)**2+(yy-60)**2 <= 12**2
    contamination={
        "status":"CATALOGUE_CONTAMINATION_MASK_NOT_ADOPTED",
        "project_id":uuid,"band":"r",
        "dataset_source":dataset["source"],
        "draft_digest":draft.digest,"catalog_sha256":"a"*64,
        "auto_mask":np.zeros(f.shape,bool),
        "protected":protected,
    }
    return dataset,draft,contamination


def test_aperture_fourth_tab_restore_no_automatic_fits_work(qtbot):
    widget=PhotometryWorkspace()
    qtbot.addWidget(widget)
    assert widget.stage_tabs.count()==4
    assert widget.stage_tabs.tabText(3)=="4 · Aperture Photometry"
    assert widget.aperture_stage._future is None
    assert widget.aperture_stage._result is None
    assert not widget.aperture_stage.run_btn.isEnabled()
    widget.restore_navigation(stage="aperture",band="r")
    assert widget.stage_tabs.currentWidget() is widget.aperture_stage
    assert widget.current_navigation()==("aperture","r")
    assert widget.scene_editor._loaded is None
    assert widget.aperture_stage._future is None
    widget.restore_navigation(stage="scene",band="r")
    assert widget.current_navigation()==("scene","r")


def test_aperture_native_worker_retains_masks_and_zero_sky():
    ds,draft,contam=_inputs()
    original=ds["data"]["image"].copy()
    good=build_project_aperture_growth(
        ds,draft,contam,source_triplet=("dr11","south","1554p180"),
        target_icrs=(155.4,18.))
    assert good["status"]=="DIRECT_FINITE_APERTURE_DIAGNOSTICS_NOT_ADOPTED"
    assert good["band"]=="r"
    assert len(good["curves"])==3
    assert good["sky_method"]=="zero_additional_correction"
    assert not good["adopted_for_science"]
    assert not good["masked_light_reconstructed"]
    assert not good["total_galaxy_magnitude_inferred"]
    assert good["source_fits_hashes_verified"] is False
    assert good["systematics_04_3"]["status"] == (
        "CONDITIONAL_SYSTEMATICS_DIAGNOSTIC_ONLY")
    assert good["systematics_04_3"]["blank_apertures"][
        "covariance_calibrated_for_science"] is False
    assert good["systematics_04_3"]["sky_response"][
        "sky_level_uncertainty_assumed"] is None
    assert good["n_sky_target_protection_pixels_not_masked_for_photometry"]>200
    assert all(len(c["cumulative"])==10 for c in good["curves"])
    assert np.array_equal(ds["data"]["image"],original)
    with pytest.raises(ValueError,match="stale"):
        build_project_aperture_growth(
            ds,draft,dict(contam,draft_digest="b"*64),
            source_triplet=("dr11","south","1554p180"),
            target_icrs=(155.4,18.))


def test_aperture_gui_worker_copy_and_invalidate(qtbot):
    ds,draft,contam=_inputs()
    scene=SimpleNamespace(
        _project_target=(155.4,18.),
        _source_key=("root",ds["project_id"],
                     ("dr11","south","1554p180","r"),"source"),
        _sga_candidate=None,
        sga_use_prior=SimpleNamespace(isChecked=lambda:False),
        sga_use_seed=SimpleNamespace(isChecked=lambda:False),
    )
    def sky_snapshot():
        return ((
            ds["project_id"],ds["source"],ds["band"],draft.digest,
            contam["catalog_sha256"],id(ds),id(contam)),
            ds,draft,contam)
    sky=SimpleNamespace(_snapshot=sky_snapshot)
    w=ApertureWorkspace(scene,sky)
    qtbot.addWidget(w)
    assert w.run_btn.isEnabled()
    assert not w.copy_btn.isEnabled()
    assert w._future is None
    w.run_diagnostic()
    qtbot.waitUntil(lambda:w._future is None,timeout=15000)
    assert w._result is not None
    assert w.copy_btn.isEnabled()
    assert w._result["status"]=="DIRECT_FINITE_APERTURE_DIAGNOSTICS_NOT_ADOPTED"
    assert "No masked galaxy light" not in w.readout.toPlainText()
    assert "UNADOPTED" in w.readout.toPlainText()
    assert "04.3 · Conditional uncertainty diagnostics" in w.readout.toPlainText()
    assert len(w.plot._curves)==3
    w.copy_qa()
    copied=QApplication.clipboard().text()
    assert "DIRECT_FINITE_APERTURE_DIAGNOSTICS_NOT_ADOPTED" in copied
    assert '"adopted_for_science": false' in copied
    # Stage 02 re-prepared a different draft: old QA can no longer be copied.
    draft=SkyMaskDraft(ds["project_id"]).with_edit(
        __import__("lsbpipe.photometry.sky_masks",fromlist=["SkyMaskEdit","SkyMaskShape"]).SkyMaskEdit(
            "new_manual","exclude",
            __import__("lsbpipe.photometry.sky_masks",fromlist=["SkyMaskShape"]).SkyMaskShape(
                "circle",(155.4,18.),radius_arcsec=1.),band="r"))
    contam=dict(contam,draft_digest=draft.digest)
    w.refresh_context()
    assert w._result is None
    assert not w.copy_btn.isEnabled()
    assert not w.plot._curves
