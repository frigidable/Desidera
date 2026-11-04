"""Offline GUI PSF-0 navigation and scientifically rejected Leo P pilot."""
from __future__ import annotations

import pytest
import numpy as np
pytest.importorskip("pytestqt")
pytest.importorskip("PySide6.QtWebEngineWidgets")

from lsbpipe.gui.survey_explorer import SurveyExplorer
from lsbpipe.gui.psf0_panel import PSF0Panel
from PySide6.QtWidgets import QApplication
from lsbpipe.photometry.psf_signed_qa import signed_psf_qa
from lsbpipe.surveys.object_projects import create_object_project, save_object_project
from lsbpipe.surveys.dataset_integrity import _expected_paths
from lsbpipe.surveys.legacy import SkyField
from scripts.offline_gui_fixture import offline_legacy_viewer


def test_psf_gui_uses_existing_project_and_chooses_present_r_instead_of_saved_g(
        qtbot, tmp_path, monkeypatch):
    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path / "session"))
    project = create_object_project(
        "Leo P", tmp_path / "Leo_P",
        field=SkyField(155.438012, 18.088025, 5.))
    with offline_legacy_viewer():
        window = SurveyExplorer(restore_session=False, launch_remote=False)
    qtbot.addWidget(window)
    window.attach_object_project(project)
    window.show()

    # Attach local inventory AFTER opening the window, avoiding any
    # legacy auto-index attempt on intentionally non-FITS test placeholders.
    project.data_root = tmp_path / "Leo_P" / "originals"
    project.data_root.mkdir()
    project.selected_dataset = ("dr11", "south", "1554p180", "g")
    save_object_project(project)
    paths = _expected_paths(project.data_root, ("dr11", "south", "1554p180", "r"))
    for p in paths.values():
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(b"not a FITS; UI existence-only fixture")

    configurations = []

    def fake_backend(root, *, band, config, protection_arcsec, progress=None):
        configurations.append(config)
        if progress is not None:
            progress((0, 1, "Synthetic QA-rejected star sample"))
        assert root == str(project.root)
        assert band == "r"
        assert config.guessed_fwhm_pix == pytest.approx(4.8)
        assert protection_arcsec == pytest.approx(90.)
        signed = np.zeros((21,21), float)
        signed[10,10] = 1.
        signed[10,15] = -(.075/.925)
        signed_qa = signed_psf_qa(signed, 4.8)
        return {
            "project_root": root, "project_id": project.project_id,
            "dataset": project.selected_dataset,
            "status": "PSF0_QUALITY_REJECTED_NOT_ADOPTED",
            "reason": "Strong negative EPSF ringing: 7.50% exceeds 2.50%",
            "diagnostics": {
                "negative_absolute_mass_fraction": signed_qa["negative_absolute_mass_fraction"],
                "blocking_threshold_fraction": .025,
                "signed_kernel_qa": signed_qa,
                "signed_kernel_values_for_qa_only": signed.tolist(),
                "n_train": 10, "n_holdout": 3,
                "training": [{"x_pix": 20., "y_pix": 30.,
                              "cutout_sky_nMgy": .2,
                              "gradient_sigma_span": .7}],
                "holdout": [],
                "rejected_reasons": {"NEIGHBOR_WITHIN_ISOLATION_RADIUS": 2},
                "notes": ["Not a science PSF"],
            },
        }

    import lsbpipe.photometry.psf_gui_backend as backend
    monkeypatch.setattr(backend, "run_gui_psf0", fake_backend)
    window.photometry_btn.click()
    panel = window.photometry_workspace.psf0_panel
    assert window.preview_tabs.currentWidget() is window.photometry_workspace
    assert panel.run_btn.isEnabled()
    assert panel.experiment.count() == 3
    assert tuple(window.photometry_workspace.psf0_panels) == ("g", "r", "i", "z")
    assert window.photometry_workspace.psf0_panels["r"] is panel
    assert panel.experiment.currentData() == "baseline"
    assert panel.band.currentData() == "r"
    assert panel.protection.value() == pytest.approx(90.)
    assert not panel.save_btn.isEnabled()
    before = {key:path.read_bytes() for key,path in paths.items()}
    panel.run_btn.click()
    qtbot.waitUntil(lambda: window.thread is None, timeout=15000)
    assert "REJECTED" in panel.status.text()
    assert "7.50%" in panel.summary.text()
    assert "Strong negative EPSF ringing" in panel.qa.toPlainText()
    assert panel.stars.rowCount() == 1
    assert panel.stars.columnCount() == 6
    assert panel.stars.horizontalHeaderItem(0).text() == "Star ID"
    assert panel.kernel_view.pixmap() is not None
    assert not panel.kernel_view.pixmap().isNull()
    assert "Intermediate" in panel.qa.toPlainText()
    assert panel.copy_qa_btn.isEnabled()
    panel.copy_qa_btn.click()
    copied = QApplication.clipboard().text()
    assert "PSF0_QUALITY_REJECTED_NOT_ADOPTED" in copied
    assert "signed_kernel_qa" in copied
    assert "signed_kernel_values_for_qa_only" not in copied
    assert '"comparison_trials_not_adopted"' in copied
    assert "baseline: neg 7.50%" in panel.compare.text()
    assert not panel.save_btn.isEnabled()
    # A second reconstruction can change the outer-sky annulus while the
    # historical baseline remains visible and traceable as non-adopted.
    panel.experiment.setCurrentIndex(panel.experiment.findData("outer_sky"))
    assert not panel.copy_qa_btn.isEnabled()
    assert "baseline: neg 7.50%" in panel.compare.text()
    panel.run_btn.click()
    qtbot.waitUntil(lambda: window.thread is None, timeout=15000)
    assert len(configurations) == 2
    assert configurations[-1].cutout_size_pix == 51
    assert configurations[-1].sky_annulus_inner_pix == 19.
    assert "outer_sky: neg 7.50%" in panel.compare.text()
    assert "baseline: neg 7.50%" in panel.compare.text()
    assert panel.copy_qa_btn.isEnabled()
    panel.copy_qa_btn.click()
    copied_comparison = QApplication.clipboard().text()
    assert '"experiment": "baseline"' in copied_comparison
    assert '"experiment": "outer_sky"' in copied_comparison
    # Changing the recipe must invalidate rejection evidence, not just a
    # successful result, to prevent accidental cross-recipe reporting.
    panel.fwhm.setValue(5.5)
    assert "baseline: neg" not in panel.compare.text()
    assert "outer_sky: neg" not in panel.compare.text()
    assert not panel.copy_qa_btn.isEnabled()
    assert not (project.root / "analyses" / "psf-pilots").exists()
    assert all(paths[k].read_bytes() == data for k,data in before.items())
    assert project.selected_dataset == ("dr11", "south", "1554p180", "g")
    window._returning_to_hub = True
    window.close()


def test_standalone_panel_has_no_results_or_actions_without_project(qtbot):
    panel = PSF0Panel()
    qtbot.addWidget(panel)
    panel.show()
    assert not panel.run_btn.isEnabled()
    assert not panel.save_btn.isEnabled()
    assert panel.active_preview is None
    assert not panel.copy_qa_btn.isEnabled()


def test_star_score_history_is_not_dumped_into_on_screen_catalog_summary():
    from lsbpipe.gui.psf0_panel import PSF0Panel
    compact=PSF0Panel._catalog_screening_summary({
        "policy":"gaia_strict",
        "catalog_source_sha256":"a"*64,
        "training_residual_screening_history":[{
            "screened_stars":[{"star_id":"star-0001",
                              "score":20.,"role":"train"}]}],
        "automatic_training_rejections":["star-0001"],
        "brightness_score_coupling_blocked_auto_trim":True,
    })
    assert "gaia_strict" in compact
    assert "brightness_score_coupling_blocked_auto_trim" in compact
    assert "screened_stars" not in compact
    assert "star-0001" not in compact


def test_gui_psfsize_suggestion_uses_real_provenance_and_updates_spinbox(qtbot):
    from pathlib import Path
    panel = PSF0Panel()
    qtbot.addWidget(panel)
    panel._project_root = Path("C:/Leo_P")
    panel.band.addItem("g · native", "g")
    result = {
        "project_root":str(panel._project_root), "band":"g",
        "seeing":{
            "available":True,
            "source":"LOCAL_DR11_TRACTOR_NEARBY_PSF_SIZE_NOT_PIXEL_MAP",
            "path":"C:/Leo_P/dr11/south/tractor/155/tractor-1554p180.fits",
            "sample":"median of 40 local BRICK_PRIMARY PSFSIZE_G entries",
            "psfsize_arcsec":1.2, "image_pixel_scale_arcsec":.262,
            "suggested_detection_fwhm_pix":1.2/.262,
            "native_map_status":"Native PSFSIZE map rejected: WCS mismatch",
        },
    }
    panel.show_fwhm_suggestion(result)
    assert panel.fwhm.value() == pytest.approx(4.58)
    assert "Tractor" in panel.prior_source.text()
    assert "WCS mismatch" in panel.prior_source.text()
    assert "4.58 px" in panel.status.text()
    assert "Not an empirical ePSF" in panel.prior_source.text()


def test_failed_psfsize_preserves_manual_guess_with_visible_reason(qtbot):
    from pathlib import Path
    panel=PSF0Panel()
    qtbot.addWidget(panel)
    panel._project_root=Path("C:/Leo_P")
    panel.band.addItem("g", "g")
    panel.fwhm.setValue(5.25)
    panel.show_fwhm_suggestion({
        "project_root":str(panel._project_root),"band":"g",
        "seeing":{"available":False,
                  "source":"NO_VALID_LOCAL_DR11_PSF_SIZE_PRIOR",
                  "message":"Expected native PSFSIZE file not found"},
    })
    assert panel.fwhm.value()==pytest.approx(5.25)
    assert "Expected native PSFSIZE file not found" in (
        panel.prior_source.text())
