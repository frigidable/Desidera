"""Offline E2E: actual FITS bytes -> native Scene -> Tractor -> Sky -> Aperture.

Every FITS pixel and Tractor row here is SYNTHETIC. Unlike GUI metadata
mocks, the actual disk readers, WCS, source identities, mask preparation,
sky worker, aperture worker and Qt transitions execute without substitutes.
This is not a real DESI observation, archival validation or science PASS.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("pytestqt")
pytest.importorskip("PySide6")
pytest.importorskip("astropy")
pytest.importorskip("photutils")

from astropy.io import fits
from astropy.wcs import WCS
from PySide6.QtWidgets import QApplication

from lsbpipe.gui.photometry_workspace import PhotometryWorkspace
from lsbpipe.photometry.mask_drafts import read_active_mask_draft
from lsbpipe.surveys.dataset_integrity import _expected_paths
from lsbpipe.surveys.legacy import SkyField
from lsbpipe.surveys.object_projects import (
    create_object_project, save_object_project,
)


def _write_original_synthetic_quartet(project, root):
    """Generate actual on-disk co-registered native-style FITS (not fake arrays)."""
    size = 384
    yy, xx = np.indices((size, size), dtype=float)
    center = (size-1)/2
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [center+1., center+1.]
    wcs.wcs.crval = [project.field.ra_deg, project.field.dec_deg]
    wcs.wcs.cd = np.diag([-.262, .262])/3600.
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    image = (0.00035 +
             0.024*np.exp(-.5*((xx-center)**2+(yy-center)**2)/22.**2) +
             0.42*np.exp(-.5*((xx-305.)**2+(yy-245.)**2)/2.4**2))
    noise = np.random.default_rng(20261114).normal(0,.002,(size,size))
    image = (image+noise).astype("f4")
    invvar = np.full((size,size),1/.002**2,dtype="f4")
    nexp = np.full((size,size),2,dtype="i2")
    maskbits = np.zeros((size,size),dtype="i2")
    invvar[14:18, 11:18] = 0
    maskbits[90:96, 102:108] = 1 << 1  # caution only; not a hard science mask
    paths = _expected_paths(root,project.selected_dataset)
    arrays = {"image":image,"invvar":invvar,"nexp":nexp,
              "maskbits":maskbits}
    for kind,path in paths.items():
        path.parent.mkdir(parents=True,exist_ok=True)
        hdr = wcs.to_header()
        if kind == "image":
            hdr["BUNIT"]="nanomaggies"
        if kind == "invvar":
            hdr["BUNIT"]="1/nanomaggies^2"
        if kind == "maskbits":
            # The official optical MASKBITS is HDU1. A generic primary
            # raster is invalid even when its FITS bytes and WCS are sound.
            fits.HDUList([
                fits.PrimaryHDU(),
                fits.ImageHDU(data=arrays[kind],header=hdr,name="MASKBITS"),
            ]).writeto(path)
        else:
            fits.PrimaryHDU(data=arrays[kind],header=hdr).writeto(path)
    return wcs,paths


def _write_synthetic_same_brick_tractor(root,wcs):
    """Real FITS binary table with required native Tractor/Gaia schema."""
    coords = np.array([[191.5,191.5],[305.,245.]],dtype=float)
    ra,dec=wcs.all_pix2world(coords[:,0],coords[:,1],0)
    cols=[
        fits.Column(name="BRICKNAME",format="8A",array=["1554p180"]*2),
        fits.Column(name="BRICK_PRIMARY",format="L",array=[True,True]),
        fits.Column(name="RA",format="D",array=ra),
        fits.Column(name="DEC",format="D",array=dec),
        fits.Column(name="REF_CAT",format="2A",array=["","G3"]),
        fits.Column(name="REF_ID",format="K",array=[0,123456789]),
        fits.Column(name="TYPE",format="4A",array=["PSF","PSF"]),
        fits.Column(name="OBJID",format="J",array=[101,102]),
        fits.Column(name="FITBITS",format="J",array=[0,0]),
        fits.Column(name="MASKBITS",format="J",array=[0,0]),
        fits.Column(name="GAIA_DUPLICATED_SOURCE",format="L",
                    array=[False,False]),
        fits.Column(name="FLUX_R",format="D",array=[80.,120.]),
        fits.Column(name="FLUX_IVAR_R",format="D",array=[4.,4.]),
    ]
    path=root/"dr11"/"south"/"tractor"/"155"/"tractor-1554p180.fits"
    path.parent.mkdir(parents=True,exist_ok=True)
    fits.HDUList([fits.PrimaryHDU(),
                  fits.BinTableHDU.from_columns(cols)]).writeto(path)
    return path


def _draw_native_circle(scene,x,y,radius_pix=7.):
    """Actual two-click Scene geometry path (FITS rows inverted for Qt)."""
    scene.scene_task_buttons["draw"].click()
    scene.mode.setCurrentIndex(scene.mode.findData("circle"))
    scene.purpose.setCurrentIndex(scene.purpose.findData("exclude"))
    scene.scope.setCurrentIndex(scene.scope.findData("band"))
    ny=scene._loaded["pixel_shape"][0]
    scene._on_click(x+.5,ny-y-.5)
    scene._on_click(x+radius_pix+.5,ny-y-.5)


def test_synthetic_original_fits_full_photometry_gui_route(
        qtbot,tmp_path,monkeypatch):
    monkeypatch.setenv("LSBPIPE_SESSION_DIR",str(tmp_path/"session"))
    project=create_object_project(
        "Synthetic native E2E",tmp_path/"galaxy",
        field=SkyField(155.438012,18.088025,1.3))
    project.selected_dataset=("dr11","south","1554p180","r")
    project.data_root=tmp_path/"originals"
    project.data_root.mkdir()
    project.target_radius_arcsec=12.
    save_object_project(project)
    wcs,paths=_write_original_synthetic_quartet(project,project.data_root)
    tractor=_write_synthetic_same_brick_tractor(project.data_root,wcs)
    source_files={**paths,"tractor":tractor}
    before={key:hashlib.sha256(path.read_bytes()).hexdigest()
            for key,path in source_files.items()}
    before_passport=(project.root/"lsbpipe-object.json").read_bytes()

    widget=PhotometryWorkspace()
    qtbot.addWidget(widget)
    widget.resize(1366,768)
    widget.show()
    widget.bind_project(project)
    scene=widget.scene_editor
    sky=widget.sky_stage
    aperture=widget.aperture_stage
    assert scene._loaded is None
    assert sky._result is None and aperture._result is None
    assert scene._future is None
    widget.restore_navigation(stage="scene",band="r")
    assert scene._band=="r"
    assert scene.load_btn.isEnabled()

    # No mock reader: read four actual FITS files and verify exact hashes.
    scene.load_btn.click()
    qtbot.waitUntil(lambda:scene._future is None,timeout=30000)
    assert scene._loaded is not None,scene.status.text()
    native=scene._loaded
    assert native["original_fits_provenance_status"] == (
        "FULL_ORIGINAL_FITS_SHA256_VERIFIED_AT_SCENE_LOAD")
    assert native["source_quartet_key"]==["dr11","south","1554p180","r"]
    assert native["original_fits_sha256_local"] == {
        key:before[key] for key in paths}
    assert all(not native["data"][kind].flags.writeable for kind in paths)
    native_image=native["data"]["image"].copy()

    # Draw and save an actual ICRS region through the Scene GUI callbacks.
    _draw_native_circle(scene,90.,180.)
    assert len(scene.draft.edits)==1
    assert scene.draft.edits[0].purpose=="exclude"
    scene.save_btn.click()
    saved,digest=read_active_mask_draft(project.root)
    assert saved is not None and len(saved.edits)==1
    assert digest==scene.draft.digest

    # Prepare masks from the actual local binary Tractor FITS, not a
    # monkeypatched catalogue object, synthetic JSON or a network lookup.
    scene.scene_task_buttons["isophote"].click()
    scene.iso_prepare.click()
    qtbot.waitUntil(lambda:scene._future is None,timeout=30000)
    assert scene._contamination is not None,scene.iso_mask_status.text()
    mask=scene._contamination
    assert mask["source_quartet_key"]==native["source_quartet_key"]
    assert mask["catalog_sha256"]==before["tractor"]
    assert mask["draft_digest"]==digest
    assert mask["stats"]["catalog_rows"]==2
    assert mask["stats"]["masked_catalog_sources"]>=1
    assert not np.array_equal(mask["auto_mask"],mask["protected"])
    assert sky._snapshot() is not None
    assert aperture._snapshot() is not None

    widget.restore_navigation(stage="sky",band="r")
    assert sky.run_btn.isEnabled()
    assert not sky.open_scene_btn.isVisible()
    sky.run_btn.click()
    qtbot.waitUntil(lambda:sky._future is None,timeout=90000)
    assert sky._result is not None,sky.status.text()
    assert sky._result["catalog_sha256"]==before["tractor"]
    assert sky._result["science_status"]=="UNADOPTED_DIAGNOSTIC_ONLY"
    assert not sky._result["sky_model_adopted"]
    assert not sky.readout.isVisible()
    sky.report_toggle.click()
    assert sky.readout.isVisible()
    sky.copy_btn.click()
    sky_json=json.loads(QApplication.clipboard().text())
    assert sky_json["catalog_sha256"]==before["tractor"]
    assert not sky_json["sky_model_adopted"]

    widget.restore_navigation(stage="aperture",band="r")
    assert aperture.run_btn.isEnabled()
    aperture.run_btn.click()
    qtbot.waitUntil(lambda:aperture._future is None,timeout=90000)
    assert aperture._result is not None,aperture.status.text()
    result=aperture._result
    assert result["status"]=="DIRECT_FINITE_APERTURE_DIAGNOSTICS_NOT_ADOPTED"
    assert result["source_quartet_key"]==native["source_quartet_key"]
    assert result["source_fits_hashes_verified"]
    assert result["original_fits_sha256_local"]=={
        key:before[key] for key in paths}
    assert result["n_measurable_geometry_hypotheses"]>=1
    assert result["sky_method"]=="zero_additional_correction"
    assert not result["adopted_for_science"]
    assert not result["masked_light_reconstructed"]
    assert aperture.plots_panel.isVisible()
    aperture.copy_btn.click()
    aperture_json=json.loads(QApplication.clipboard().text())
    assert not aperture_json["adopted_for_science"]

    # Capture the *actual successful synthetic data-bearing stages*
    # before invalidating their results; these are not archival observations.
    capture_dir=os.environ.get("DESIDERA_E2E_SCREENSHOT_DIR")
    if capture_dir:
        screenshot_root=Path(capture_dir)
        screenshot_root.mkdir(parents=True,exist_ok=True)
        for stage in ("scene","sky","aperture"):
            widget.restore_navigation(stage=stage,band="r")
            QApplication.processEvents()
            pix=widget.grab()
            if pix.isNull() or not pix.save(str(
                    screenshot_root/f"synthetic_native_{stage}_1366x768.png"),
                    "PNG"):
                raise AssertionError(f"Cannot capture {stage} fixture")

    # A second manual edit immediately invalidates BOTH derived QA views.
    widget.restore_navigation(stage="scene",band="r")
    _draw_native_circle(scene,100.,120.)
    assert len(scene.draft.edits)==2
    assert scene._contamination is None
    sky.refresh_context()
    aperture.refresh_context()
    assert sky._result is None and aperture._result is None
    assert not sky.copy_btn.isEnabled()
    assert not aperture.copy_btn.isEnabled()
    assert not aperture.plots_panel.isVisible()
    scene.save_btn.click()
    updated,updated_digest=read_active_mask_draft(project.root)
    assert len(updated.edits)==2 and updated_digest!=digest
    np.testing.assert_array_equal(native["data"]["image"],native_image)

    # The only on-disk edits are the explicit mask draft and GUI resume
    # receipts; the original four FITS, Tractor, and project passport
    # are byte-identical after every stage.
    assert (project.root/"lsbpipe-object.json").read_bytes()==before_passport
    assert {key:hashlib.sha256(path.read_bytes()).hexdigest()
            for key,path in source_files.items()}==before

    if capture_dir:
        widget.restore_navigation(stage="aperture",band="r")
        QApplication.processEvents()
        stale=widget.grab()
        if stale.isNull() or not stale.save(
                str(screenshot_root/"synthetic_native_aperture_stale_1366x768.png"),
                "PNG"):
            raise AssertionError("Cannot capture invalidated aperture")
        (screenshot_root/"synthetic_native_receipt.json").write_text(
            json.dumps({
                "fixture":"GENERATED_SYNTHETIC_FITS_NOT_DESI_OBSERVATION",
                "science_adopted":False,
                "original_fits_unchanged":True,
                "original_fits_sha256":{key:before[key] for key in paths},
                "tractor_sha256":before["tractor"],
                "sky_status":sky_json["status"],
                "sky_training_pixels":sky_json["sky_train_effective_pixels"],
                "sky_holdout_pixels":sky_json["sky_holdout_effective_pixels"],
                "aperture_status":aperture_json["status"],
                "measurable_geometry_count":result[
                    "n_measurable_geometry_hypotheses"],
                "sky_model_adopted":False,
                "aperture_adopted_for_science":False,
                "real_monitor_dpi_tested":False,
            },indent=2,sort_keys=True),encoding="utf-8")
