"""Offline four-filter PSF GUI: independent results, batch seeing, QA and errors."""
import json
import pytest
pytest.importorskip("pytestqt")
pytest.importorskip("PySide6.QtWebEngineWidgets")

from PySide6.QtWidgets import QApplication
from lsbpipe.gui.survey_explorer import SurveyExplorer
from lsbpipe.surveys.object_projects import (
    create_object_project, save_object_project)
from lsbpipe.surveys.dataset_integrity import _expected_paths
from lsbpipe.surveys.legacy import SkyField
from scripts.offline_gui_fixture import offline_legacy_viewer


def _setup(qtbot, tmp_path, monkeypatch):
    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path/"session"))
    project = create_object_project(
        "Leo P", tmp_path/"Leo_P",
        field=SkyField(155.438012, 18.088025, 4.))
    with offline_legacy_viewer():
        window = SurveyExplorer(restore_session=False, launch_remote=False)
    qtbot.addWidget(window)
    window.attach_object_project(project)
    root = project.root/"originals"
    root.mkdir()
    project.data_root = root
    project.selected_dataset = ("dr11", "south", "1554p180", "g")
    save_object_project(project)
    for band in "griz":
        for src in _expected_paths(root, ("dr11","south","1554p180",band)).values():
            src.parent.mkdir(parents=True, exist_ok=True)
            src.write_bytes(b"offline placeholder -- NEVER accepted as science")
    window.photometry_btn.click()
    return window, project


def _rejected(root, band):
    return {
        "project_root":str(root), "review_band":band,
        "status":"PSF0_QUALITY_REJECTED_NOT_ADOPTED",
        "reason":"Synthetic signed ringing diagnostic -- NOT ADOPTED",
        "diagnostics":{
            "negative_absolute_mass_fraction":.05,
            "blocking_threshold_fraction":.025,
            "n_train":12,"n_holdout":4,
            "training":[],"holdout":[],"rejected_reasons":{},
            "signed_holdout_summary":{},"signed_kernel_qa":{},
            "solver_diagnostics":{},
            "notes":["offline fixture, no scientific FITS"],
            "catalog_selection_qa":{},
        },
    }


def test_four_band_batch_preserves_independent_results_and_reports(
        qtbot, tmp_path, monkeypatch):
    w, project = _setup(qtbot, tmp_path, monkeypatch)
    workspace = w.photometry_workspace
    assert list(workspace.psf0_panels) == list("griz")
    assert workspace.psf0_panel is workspace.psf0_panels["g"]
    assert all(p.experiment.count() == 3
               for p in workspace.psf0_panels.values())
    for band, panel in workspace.psf0_panels.items():
        assert panel.band.currentData() == band
        assert not panel.band.isEnabled()
        assert panel.run_btn.isEnabled()
    from lsbpipe.gui import psf0_batch_controller as batch
    seen = []
    fwhms = {"g":6.57, "r":5.41, "i":4.33, "z":3.75}

    def fake_seeing(root, *, band):
        return {
            "project_root":str(project.root),"band":band,
            "seeing":{
                "available":True,"source":"LOCAL_NATIVE_PSF_SIZE_MAP",
                "psfsize_arcsec":fwhms[band]*.262,
                "image_pixel_scale_arcsec":.262,
                "suggested_detection_fwhm_pix":fwhms[band],
                "path":f"/original/psfsize-{band}.fits.fz",
                "sample":"5x5 local median",
            },
        }

    def fake_run(root, *, band, config, protection_arcsec):
        seen.append((band, config.guessed_fwhm_pix, config.cutout_size_pix))
        if band == "r":
            raise ValueError("synthetic bad r image WCS")
        return _rejected(project.root, band)

    monkeypatch.setattr(batch, "suggest_gui_psf0_fwhm", fake_seeing)
    monkeypatch.setattr(batch, "run_gui_psf0", fake_run)
    workspace.psf0_panels["i"].experiment.setCurrentIndex(
        workspace.psf0_panels["i"].experiment.findData("outer_sky"))
    workspace.build_all_psfs_btn.click()
    qtbot.waitUntil(lambda: w.thread is None, timeout=30000)
    assert [x[0] for x in seen] == ["g","r","i","z"]
    assert [x[1] for x in seen] == pytest.approx(
        [6.57,5.41,4.33,3.75])
    assert [x[2] for x in seen] == [31,31,55,31]
    for band in ("g","i","z"):
        assert workspace.psf0_panels[band]._qa_payload is not None
        assert "REJECTED" in workspace.psf0_panels[band].status.text()
    assert "synthetic bad r image WCS" in (
        workspace.psf0_panels["r"].qa.toPlainText())
    assert workspace.psf0_panels["r"]._qa_payload is None
    assert "r: ERROR" in workspace.batch_status.text()
    assert "z: QUALITY REJECTED" in workspace.batch_status.text()
    assert workspace.psf0_panels["g"].fwhm.value() == pytest.approx(6.57)
    assert workspace.psf0_panels["z"].fwhm.value() == pytest.approx(3.75)
    workspace.copy_full_qa_btn.click()
    qa = json.loads(QApplication.clipboard().text())
    assert qa["science_adopted"] is False
    assert sorted(qa["bands"]) == list("girz")  # sorted alphabetically
    assert qa["bands"]["g"]["science_report_or_rejection"] is not None
    assert "synthetic bad r image WCS" in qa["bands"]["r"]["full_visible_qa"]
    assert qa["bands"]["z"]["detection_fwhm_guess_pix_not_measured"] == 3.75
    # Re-enter the same existing project, leaving all four independent
    # preview/rejection results untouched; no new filesystem/project.
    original = workspace.psf0_panels["g"].qa.toPlainText()
    w.photometry_btn.click()
    assert workspace.psf0_panels["g"].qa.toPlainText() == original
    assert not (project.root/"analyses"/"psf-pilots").exists()
    w._returning_to_hub = True
    w.close()


def test_manual_band_fwhm_override_is_preserved_in_batch(
        qtbot, tmp_path, monkeypatch):
    window, project = _setup(qtbot, tmp_path, monkeypatch)
    from lsbpipe.gui import psf0_batch_controller as controller
    workspace = window.photometry_workspace
    workspace.psf0_panels["r"].fwhm.setValue(5.25)
    workspace.psf0_panels["z"].experiment.setCurrentIndex(
        workspace.psf0_panels["z"].experiment.findData("outer_sky_robust"))
    measurements = []
    seeing_calls = []

    def fake_seeing(root, *, band):
        seeing_calls.append(band)
        return {
            "project_root":str(project.root), "band":band,
            "seeing":{
                "available":True,"source":"LOCAL_NATIVE_PSF_SIZE_MAP",
                "psfsize_arcsec":1.31,
                "image_pixel_scale_arcsec":.262,
                "suggested_detection_fwhm_pix":5.0,
                "path":f"/synthetic/psfsize-{band}.fits.fz",
                "sample":"synthetic fixture",
            },
        }

    def fake_run(root, *, band, config, protection_arcsec):
        measurements.append((band,config.guessed_fwhm_pix,
                             config.training_outlier_passes))
        return _rejected(project.root,band)

    monkeypatch.setattr(controller,"suggest_gui_psf0_fwhm",fake_seeing)
    monkeypatch.setattr(controller,"run_gui_psf0",fake_run)
    workspace.build_all_psfs_btn.click()
    qtbot.waitUntil(lambda:window.thread is None, timeout=30000)
    assert seeing_calls == ["g","i","z"]
    assert measurements == [
        ("g",5.0,0),("r",5.25,0),("i",5.0,0),("z",5.0,2)]
    qa=workspace.full_qa_report()
    assert qa["bands"]["r"]["fwhm_strategy"] == "manual"
    assert qa["bands"]["g"]["fwhm_strategy"] == "local_psfsize"
    assert qa["bands"]["z"]["experiment"] == "outer_sky_robust"
    window._returning_to_hub=True
    window.close()


def test_full_griz_report_saves_new_atomic_json_without_overwriting(
        qtbot, tmp_path, monkeypatch):
    window, project = _setup(qtbot, tmp_path, monkeypatch)
    workspace = window.photometry_workspace
    target = tmp_path / "PSF0_all_filters_for_ChatGPT.json"
    from lsbpipe.gui import photometry_workspace as module
    monkeypatch.setattr(
        module.QFileDialog, "getSaveFileName",
        lambda *args, **kwargs: (str(target), "JSON (*.json)"))
    workspace.save_full_qa_btn.click()
    assert target.is_file()
    raw = target.read_bytes()
    report = json.loads(raw.decode("utf-8"))
    assert report["science_adopted"] is False
    assert sorted(report["bands"]) == list("girz")
    assert report["no_raw_fits_or_stamp_pixels_in_clipboard"] is True
    workspace.save_full_qa_btn.click()
    assert target.read_bytes() == raw
    assert "Existing files" in workspace.batch_status.text()
    assert not (project.root/"analyses"/"psf-pilots").exists()
    window._returning_to_hub=True
    window.close()
