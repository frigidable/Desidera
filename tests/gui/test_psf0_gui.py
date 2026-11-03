"""Offline GUI PSF-0 navigation and scientifically rejected Leo P pilot."""
from __future__ import annotations

import pytest
pytest.importorskip("pytestqt")
pytest.importorskip("PySide6.QtWebEngineWidgets")

from lsbpipe.gui.survey_explorer import SurveyExplorer
from lsbpipe.gui.psf0_panel import PSF0Panel
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

    def fake_backend(root, *, band, config, protection_arcsec, progress=None):
        if progress is not None:
            progress((0, 1, "Synthetic QA-rejected star sample"))
        assert root == str(project.root)
        assert band == "r"
        assert config.guessed_fwhm_pix == pytest.approx(4.8)
        assert protection_arcsec == pytest.approx(90.)
        return {
            "project_root": root, "project_id": project.project_id,
            "dataset": project.selected_dataset,
            "status": "PSF0_QUALITY_REJECTED_NOT_ADOPTED",
            "reason": "Strong negative EPSF ringing: 7.50% exceeds 2.50%",
            "diagnostics": {
                "negative_absolute_mass_fraction": .075,
                "blocking_threshold_fraction": .025,
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
    assert not panel.save_btn.isEnabled()
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
