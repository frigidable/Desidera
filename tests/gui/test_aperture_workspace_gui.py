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
    assert good["source_quartet_key"] == ["dr11","south","1554p180","r"]
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


def test_aperture_refuses_foreign_source_quartet_even_if_display_label_matches():
    ds,draft,contam=_inputs()
    ds["source_quartet_key"]=["dr11","south","1555p180","r"]
    with pytest.raises(ValueError,match="source-triplet mismatch"):
        build_project_aperture_growth(
            ds,draft,contam,source_triplet=("dr11","south","1554p180"),
            target_icrs=(155.4,18.))
    ds["source_quartet_key"]=["dr11","south","1554p180","r"]
    contam["source_quartet_key"]=["dr11","south","1555p180","r"]
    with pytest.raises(ValueError,match="source-triplet mismatch"):
        build_project_aperture_growth(
            ds,draft,contam,source_triplet=("dr11","south","1554p180"),
            target_icrs=(155.4,18.))
    contam["source_quartet_key"]=list(ds["source_quartet_key"])
    ds["loaded_target_icrs"]=[155.399,18.]
    with pytest.raises(ValueError,match="target ICRS differs"):
        build_project_aperture_growth(
            ds,draft,contam,source_triplet=("dr11","south","1554p180"),
            target_icrs=(155.4,18.))
    ds["loaded_target_icrs"]=[155.4,18.]
    good=build_project_aperture_growth(
        ds,draft,contam,source_triplet=("dr11","south","1554p180"),
        target_icrs=(155.4,18.))
    assert good["source_quartet_key"] == ds["source_quartet_key"]


def test_aperture_recipe_exposes_changed_scene_mask_geometry_policy():
    """A fixed catalogue/draft with a new raster recipe is a new product."""
    ds,draft,contamination=_inputs()
    original=build_project_aperture_growth(
        ds,draft,contamination,source_triplet=("dr11","south","1554p180"),
        target_icrs=(155.4,18.))
    contamination["source_mask_policy_id"]=(
        "PROVISIONAL_SCENE_EFFECTIVE_AUTO_PLUS_EXCLUDE_V2_TAN_WCS")
    updated=build_project_aperture_growth(
        ds,draft,contamination,source_triplet=("dr11","south","1554p180"),
        target_icrs=(155.4,18.))
    assert updated["source_policy_id"] == contamination["source_mask_policy_id"]
    assert updated["source_policy_id"] != original["source_policy_id"]
    assert updated["curves"][0]["cumulative"][0]["recipe_sha256"] != (
        original["curves"][0]["cumulative"][0]["recipe_sha256"])
    assert not updated["adopted_for_science"]


def test_stage04_exact_mask_receipts_reject_changed_source_or_manual_pixels():
    from lsbpipe.photometry.sky_masks import native_boolean_mask_sha256
    ds,draft,contamination=_inputs()
    source_hash=native_boolean_mask_sha256(contamination["auto_mask"])
    manual_hash=native_boolean_mask_sha256(
        np.zeros(ds["pixel_shape"],dtype=bool))
    contamination["effective_source_mask_sha256"]=source_hash
    contamination["manual_exclude_mask_sha256"]=manual_hash
    good=build_project_aperture_growth(
        ds,draft,contamination,source_triplet=("dr11","south","1554p180"),
        target_icrs=(155.4,18.))
    assert good["effective_source_mask_sha256"]==source_hash
    assert good["manual_exclude_mask_sha256"]==manual_hash
    assert good["adopted_for_science"] is False
    contamination["effective_source_mask_sha256"]="f"*64
    with pytest.raises(ValueError,match="mask geometry fingerprint mismatch"):
        build_project_aperture_growth(
            ds,draft,contamination,
            source_triplet=("dr11","south","1554p180"),
            target_icrs=(155.4,18.))
    contamination["effective_source_mask_sha256"]=source_hash
    contamination["manual_exclude_mask_sha256"]="e"*64
    with pytest.raises(ValueError,match="mask geometry fingerprint mismatch"):
        build_project_aperture_growth(
            ds,draft,contamination,
            source_triplet=("dr11","south","1554p180"),
            target_icrs=(155.4,18.))


def test_aperture_source_sha_propagates_into_each_exact_flux_recipe():
    """The full original-FITS SHA identity follows Scene -> Stage04."""
    ds,draft,contam=_inputs()
    sha={k:(str(i)*64) for i,k in enumerate(
        ("image","invvar","nexp","maskbits"),1)}
    ds["original_fits_sha256_local"]=sha
    ds["source_quartet_key"]=["dr11","south","1554p180","r"]
    contam["source_quartet_key"]=list(ds["source_quartet_key"])
    ds["image_unit_evidence"]="FITS_IMAGE_BUNIT_NANOMAGGIES_CONFIRMED"
    ds["invvar_unit_evidence"]="FITS_INVVAR_BUNIT_CONFIRMED"
    ds["original_fits_provenance_status"]=(
        "FULL_ORIGINAL_FITS_SHA256_VERIFIED_AT_SCENE_LOAD")
    result=build_project_aperture_growth(
        ds,draft,contam,source_triplet=("dr11","south","1554p180"),
        target_icrs=(155.4,18.))
    assert result["source_fits_hashes_verified"] is True
    assert result["original_fits_sha256_local"] == sha
    assert result["adopted_for_science"] is False
    first=result["curves"][0]["cumulative"][0]
    assert "FULL_SOURCE_FITS_HASHES_NOT_VERIFIED" not in first["warnings"]
    before_recipe=first["recipe_sha256"]
    second=dict(sha,image="a"*64)
    ds["original_fits_sha256_local"]=second
    changed=build_project_aperture_growth(
        ds,draft,contam,source_triplet=("dr11","south","1554p180"),
        target_icrs=(155.4,18.))
    assert changed["curves"][0]["cumulative"][0]["recipe_sha256"] != (
        before_recipe)
    for bad in (
            {"image":"a"*64}, {"image":"not-a-hash",
                               "invvar":"2"*64,"nexp":"3"*64,
                               "maskbits":"4"*64}):
        ds["original_fits_sha256_local"]=bad
        with pytest.raises(ValueError,match="Malformed or inconsistent"):
            build_project_aperture_growth(
                ds,draft,contam,
                source_triplet=("dr11","south","1554p180"),
                target_icrs=(155.4,18.))
    ds["original_fits_sha256_local"]=sha
    ds["original_fits_provenance_status"]="NOT_HASHED"
    with pytest.raises(ValueError,match="Malformed or inconsistent"):
        build_project_aperture_growth(
            ds,draft,contam,source_triplet=("dr11","south","1554p180"),
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


def test_sga_aperture_geometry_uses_seed_not_protection(qtbot):
    """Stage02 protection and shape-seed toggles must stay independent."""
    ds, draft, contamination = _inputs()
    scene = SimpleNamespace(
        _project_target=(155.4, 18.),
        _source_key=("root", ds["project_id"],
                     ("dr11", "south", "1554p180", "r"), "source"),
        _sga_candidate=SimpleNamespace(ba=.47, pa_deg=73.),
        sga_use_prior=SimpleNamespace(isChecked=lambda: False),
        sga_use_seed=SimpleNamespace(isChecked=lambda: True),
    )
    def sky_snapshot():
        return ((ds["project_id"], ds["source"], ds["band"],
                 draft.digest, contamination["catalog_sha256"],
                 id(ds), id(contamination)),
                ds, draft, contamination)
    sky = SimpleNamespace(_snapshot=sky_snapshot)
    w = ApertureWorkspace(scene, sky)
    qtbot.addWidget(w)
    snapshot = w._snapshot()
    assert snapshot is not None
    assert snapshot[-2:] == (.47, 73.)
    scene.sga_use_seed = SimpleNamespace(isChecked=lambda: False)
    scene.sga_use_prior = SimpleNamespace(isChecked=lambda: True)
    snapshot2 = w._snapshot()
    assert snapshot2[-2:] == (None, None)
    assert snapshot2[0] != snapshot[0]
