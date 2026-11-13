"""Real Qt checks for the offline object-centric launcher."""
import pytest

pytest.importorskip("pytestqt")
pytest.importorskip("PySide6.QtWebEngineWidgets")

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QDialog
from scripts.offline_gui_fixture import offline_legacy_viewer
from lsbpipe.gui.project_hub import ProjectHub, CreateProjectDialog
from lsbpipe.surveys.object_projects import (
    create_object_project, read_object_project, recent_projects,
)
from lsbpipe.surveys.legacy import SkyField


@pytest.fixture
def hub(qtbot):
    with offline_legacy_viewer():
        window = ProjectHub()
    qtbot.addWidget(window)
    window.show()
    return window


def test_launcher_is_empty_and_never_opens_viewer_on_first_launch(hub):
    assert hub.objectName() == "lsbpipe_project_hub"
    assert hub.projects.count() == 0
    assert hub.empty_panel.isVisible()
    assert not hub.projects.isVisible()
    assert hub.workspace is None
    assert not hub.continue_btn.isEnabled()
    assert hub.new_btn.isEnabled()
    assert hub.import_btn.isEnabled()


def test_single_page_creation_preserves_unidentified_object(qtbot, tmp_path):
    dialog = CreateProjectDialog(default_root=tmp_path)
    qtbot.addWidget(dialog)
    assert not dialog.create_btn.isEnabled()
    dialog.name.setText("Unknown LSB dwarf")
    assert dialog.create_btn.isEnabled()
    assert "Unknown_LSB_dwarf" in dialog.directory.text()
    assert "not specified" in dialog.description.text()
    dialog._create()
    assert dialog.result() == QDialog.Accepted
    project = dialog.project
    assert project is not None and project.field is None
    assert not (project.root / "data").exists()
    assert read_object_project(project.root).project_id == project.project_id


def test_create_wizard_parses_icrs_and_rejects_partial_or_invalid_input(
        qtbot, monkeypatch, tmp_path):
    from PySide6.QtWidgets import QMessageBox
    warnings = []
    monkeypatch.setattr(
        QMessageBox, "warning",
        lambda *args, **kwargs: warnings.append(args))
    dialog = CreateProjectDialog(default_root=tmp_path)
    qtbot.addWidget(dialog)
    dialog.name.setText("Coordinate target")
    dialog.ra.setText("02:40:23.09")
    dialog._create()
    assert dialog.project is None
    assert warnings and "both RA and Dec" in dialog.description.text()
    assert not (tmp_path / "Coordinate_target").exists()
    dialog.dec.setText("+01:13:45.3")
    dialog._create()
    assert dialog.result() == QDialog.Accepted
    assert dialog.project.field.ra_deg == pytest.approx(40.0962083, abs=1e-5)
    assert dialog.project.field.dec_deg == pytest.approx(1.22925, abs=1e-5)


def test_create_wizard_refuses_existing_science_folder_without_mutation(
        qtbot, monkeypatch, tmp_path):
    from PySide6.QtWidgets import QMessageBox
    existing = tmp_path / "existing"
    existing.mkdir()
    original = existing / "image.fits.fz"
    original.write_bytes(b"immutable science bytes")
    dialog = CreateProjectDialog(default_root=tmp_path)
    qtbot.addWidget(dialog)
    dialog.name.setText("LSB target")
    dialog.directory.setText(str(existing))
    warnings = []
    monkeypatch.setattr(
        QMessageBox, "warning",
        lambda *args, **kwargs: warnings.append(args))
    dialog._create()
    assert dialog.project is None
    assert warnings
    assert original.read_bytes() == b"immutable science bytes"
    assert not (existing / "lsbpipe-object.json").exists()


def test_hub_open_create_and_return_preserves_object_identity(
        hub, qtbot, tmp_path):
    project = create_object_project(
        "UGC 2162", tmp_path / "UGC_2162",
        field=SkyField(40.096208, 1.22925, 5.))
    with offline_legacy_viewer():
        hub.refresh()
        assert hub.projects.count() == 1
        hub._open_selected()
    window = hub.workspace
    assert window is not None
    assert window._object_project.project_id == project.project_id
    assert window._local_project is None
    assert window.ra.value() == pytest.approx(40.096208)
    assert window.dec.value() == pytest.approx(1.22925)
    assert window.bricks.rowCount() == 0
    assert window.back_to_projects.isVisible()
    assert not hub.isVisible()
    window._go_back_to_projects()
    assert hub.workspace is None
    assert hub.isVisible()
    assert read_object_project(project.root).project_id == project.project_id
    hub.refresh()
    assert hub.projects.count() == 1


def test_unresolved_draft_has_no_fake_science_target(hub, tmp_path):
    project = create_object_project("NewTarget", tmp_path / "draft")
    with offline_legacy_viewer():
        hub._enter_project(project)
    w = hub.workspace
    assert w.ra.text() == ""
    assert w.dec.text() == ""
    assert w.bricks.rowCount() == 0
    assert w.manifest is None
    assert "No ICRS position" in w.target_hint.text()
    w._go_back_to_projects()
    assert read_object_project(project.root).field is None


def test_hub_details_show_saved_dataset_without_claiming_integrity(
        hub, tmp_path):
    from lsbpipe.surveys.object_projects import save_object_project
    p = create_object_project(
        "UGC 2162", tmp_path / "UGC_2162",
        field=SkyField(40.096208, 1.22925, 5.))
    data = p.root / "data"
    data.mkdir()
    p.data_root = data
    p.selected_dataset = ("dr9", "north", "0401p012", "r")
    save_object_project(p)
    hub.refresh()
    assert hub.detail_title.text() == "UGC 2162"
    assert "integrity not checked" in hub.detail_status.text()
    assert "DR9 / north / 0401p012 / r" in hub.detail_meta.text()
    assert str(p.root) in hub.detail_path.text()
    assert "not automatically adopted" in hub.hint.text()


def test_hub_filter_selects_first_visible_project_and_shows_empty_state(
        hub, tmp_path):
    create_object_project("Alpha", tmp_path / "alpha")
    create_object_project("Beta", tmp_path / "beta")
    hub.refresh()
    hub.filter.setText("alpha")
    assert hub.continue_btn.isEnabled()
    assert "Alpha" in hub.detail_title.text()
    assert hub.projects.isVisible()
    assert not hub.empty_panel.isVisible()
    hub.filter.setText("does-not-exist")
    assert not hub.continue_btn.isEnabled()
    assert not hub.remove_btn.isEnabled()
    assert hub.empty_label.isVisible()
    assert "No matching" in hub.empty_label.text()
    hub.filter.clear()
    assert hub.continue_btn.isEnabled()


def test_filters_and_missing_project_stays_in_library(hub, tmp_path):
    create_object_project("NGC 2146 satellite", tmp_path / "LSB1")
    create_object_project("UGC 2162", tmp_path / "LSB2")
    hub.refresh()
    assert hub.projects.count() == 2
    hub.filter.setText("2162")
    assert sum(not hub.projects.item(i).isHidden()
               for i in range(hub.projects.count())) == 1
    hub.filter.clear()
    root = tmp_path / "LSB2"
    root.rename(tmp_path / "moved")
    hub.refresh()
    assert any("Missing" in hub.projects.item(i).text()
               for i in range(hub.projects.count()))


def test_import_legacy_folder_into_hub_without_network(hub, tmp_path, monkeypatch):
    from lsbpipe.surveys.legacy import BrickHit, product_for
    import json
    root = tmp_path / "UGC_2162"
    data = root / "data"
    hit = BrickHit("dr9", "south", "0401p012", 0, 0, 0, 0, {"r": 1})
    for kind in ("image", "invvar", "nexp", "maskbits"):
        p = product_for(hit, kind, None if kind == "maskbits" else "r")
        path = data / p.relative_path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"existing FITS placeholder")
    (data / "download-report.json").write_text(json.dumps({
        "schema": "lsbpipe.legacy-download-report/0.1",
        "field": {"ra_deg": 40.096208, "dec_deg": 1.22925,
                  "width_arcmin": 5.}, "results": []}))
    def forbid(*args, **kwargs):
        raise AssertionError("No network request permitted")
    monkeypatch.setattr(
        "lsbpipe.gui.survey_explorer.discover_bricks", forbid)
    with offline_legacy_viewer():
        assert hub.import_folder(root)
    assert hub.workspace is not None
    assert hub.workspace._local_project is not None
    assert hub.workspace.bricks.rowCount() == 1
    assert hub.workspace._object_project.display_name == "UGC_2162"
    hub.workspace._go_back_to_projects()
    assert len(recent_projects()) == 1


def test_object_passport_wins_over_stale_legacy_on_gui_reopen(
        hub, qtbot, tmp_path):
    """Do not replace a revised galaxy center with old download-report ICRS."""
    import json
    from lsbpipe.surveys.legacy import BrickHit, product_for
    from lsbpipe.surveys.object_projects import (
        import_object_project, save_object_project)
    root = tmp_path / "UGC_2162"
    data = root / "data"
    for region in ("north", "south"):
        hit = BrickHit("dr9", region, "0401p012", 0, 0, 0, 0, {"r": 1})
        for kind in ("image", "invvar", "nexp", "maskbits"):
            p = product_for(
                hit, kind, None if kind == "maskbits" else "r")
            filename = data / p.relative_path
            filename.parent.mkdir(parents=True, exist_ok=True)
            filename.write_bytes(b"unaltered original FITS placeholder")
    legacy = root / "lsbpipe-project.json"
    legacy.write_text(json.dumps({
        "schema": "lsbpipe.local-project/0.1",
        "field": {"ra_deg": 40.096208, "dec_deg": 1.22925,
                  "width_arcmin": 5.},
        "selected": ["dr9", "north", "0401p012", "r"],
    }), encoding="utf-8")
    old_bytes = legacy.read_bytes()
    project = import_object_project(root)
    revised = SkyField(40.101234, 1.225678, 6.)
    project.field = revised
    project.selected_dataset = ("dr9", "south", "0401p012", "r")
    project.target_radius_arcsec = 35.
    project.sky_inner_arcsec = 65.
    project.sky_outer_arcsec = 110.
    save_object_project(project)
    with offline_legacy_viewer():
        hub._enter_project(project)
    w = hub.workspace
    assert w is not None and w._local_project is not None
    assert w.ra.value() == pytest.approx(revised.ra_deg)
    assert w.dec.value() == pytest.approx(revised.dec_deg)
    assert w.field_width.value() == pytest.approx(6.)
    assert w._local_project.field == revised
    assert w._local_project.selected == project.selected_dataset
    assert w.hits[w.bricks.currentRow()].region == "south"
    assert (w.target_radius.value(), w.sky_inner.value(),
            w.sky_outer.value()) == (35., 65., 110.)
    assert legacy.read_bytes() == old_bytes

    # Editing the ICRS center never updates the original/legacy files.
    w.ra.setText("40.111111")
    w.target_radius.setValue(40.)
    w._go_back_to_projects()
    assert hub.workspace is None
    current = read_object_project(root)
    assert current.field.ra_deg == pytest.approx(40.111111)
    assert current.field.dec_deg == pytest.approx(revised.dec_deg)
    assert current.selected_dataset == ("dr9", "south", "0401p012", "r")
    assert current.target_radius_arcsec == pytest.approx(40.)
    assert legacy.read_bytes() == old_bytes

    # Closing/reopening must not re-import the original download position.
    with offline_legacy_viewer():
        hub._enter_project(current)
    reopened = hub.workspace
    assert reopened.ra.value() == pytest.approx(40.111111)
    assert reopened.dec.value() == pytest.approx(revised.dec_deg)
    assert reopened._local_project.field.ra_deg == pytest.approx(40.111111)
    assert reopened.hits[reopened.bricks.currentRow()].region == "south"
    assert legacy.read_bytes() == old_bytes
    reopened._go_back_to_projects()


def test_empty_galaxy_qa_radii_are_persisted_without_fits(hub, tmp_path):
    project = create_object_project("Unobserved dwarf", tmp_path / "draft")
    with offline_legacy_viewer():
        hub._enter_project(project)
    w = hub.workspace
    w.target_radius.setValue(35.)
    w.sky_inner.setValue(65.)
    w.sky_outer.setValue(110.)
    w._go_back_to_projects()
    again = read_object_project(project.root)
    assert (again.target_radius_arcsec, again.sky_inner_arcsec,
            again.sky_outer_arcsec) == (35., 65., 110.)
    with offline_legacy_viewer():
        hub._enter_project(again)
    reopened = hub.workspace
    assert (reopened.target_radius.value(), reopened.sky_inner.value(),
            reopened.sky_outer.value()) == (35., 65., 110.)
    reopened._go_back_to_projects()


def test_missing_saved_dataset_requires_explicit_new_selection(
        hub, tmp_path):
    import json
    from lsbpipe.surveys.legacy import BrickHit, product_for
    from lsbpipe.surveys.object_projects import (
        import_object_project, save_object_project)
    root = tmp_path / "UGC_2162"
    data = root / "data"
    hit = BrickHit("dr9", "south", "0401p012", 0, 0, 0, 0, {"r": 1})
    for kind in ("image", "invvar", "nexp", "maskbits"):
        product = product_for(
            hit, kind, None if kind == "maskbits" else "r")
        path = data / product.relative_path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"read only")
    (data / "download-report.json").write_text(json.dumps({
        "schema": "lsbpipe.legacy-download-report/0.1",
        "field": {"ra_deg": 40.096208, "dec_deg": 1.22925,
                  "width_arcmin": 5.0},
    }), encoding="utf-8")
    project = import_object_project(root)
    missing = ("dr9", "north", "0401p012", "r")
    project.selected_dataset = missing
    save_object_project(project)
    with offline_legacy_viewer():
        hub._enter_project(project)
    w = hub.workspace
    assert w is not None
    assert w.bricks.rowCount() == 1
    assert not w.bricks.selectedItems()
    assert "Saved dataset missing" in w.project_status.text()
    w._go_back_to_projects()
    assert read_object_project(root).selected_dataset == missing

    # Only deliberate row selection changes the object passport.
    with offline_legacy_viewer():
        hub._enter_project(read_object_project(root))
    w = hub.workspace
    w.bricks.selectRow(0)
    w._go_back_to_projects()
    assert read_object_project(root).selected_dataset == (
        "dr9", "south", "0401p012", "r")


def test_project_wide_dataset_audit_runs_without_brick_selection(
        hub, qtbot, tmp_path):
    """Real worker: audit original FITS and save receipt inside project root."""
    import json
    import numpy as np
    from astropy.io import fits
    from astropy.wcs import WCS
    from lsbpipe.surveys.legacy import BrickHit, product_for
    field = SkyField(40.096208, 1.22925, 1.)
    p = create_object_project("UGC 2162 deep audit", tmp_path / "UGC_2162",
                              field=field)
    data = p.root / "data"
    hit = BrickHit("dr9", "south", "0401p012", 0, 0, 0, 0, {"r": 2})
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [32.5, 32.5]
    wcs.wcs.cdelt = [-1/3600, 1/3600]
    wcs.wcs.crval = [field.ra_deg, field.dec_deg]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    arrays = {
        "image": np.ones((64, 64), dtype=np.float32),
        "invvar": np.ones((64, 64), dtype=np.float32),
        "nexp": np.ones((64, 64), dtype=np.int16),
        "maskbits": np.zeros((64, 64), dtype=np.int16),
    }
    for kind, arr in arrays.items():
        item = product_for(hit, kind, None if kind == "maskbits" else "r")
        file = data / item.relative_path
        file.parent.mkdir(parents=True, exist_ok=True)
        if kind == "maskbits":
            fits.HDUList([fits.PrimaryHDU(),
                fits.ImageHDU(data=arr, header=wcs.to_header())]).writeto(file)
        else:
            fits.PrimaryHDU(data=arr, header=wcs.to_header()).writeto(file)
    with offline_legacy_viewer():
        hub._enter_project(p)
    w = hub.workspace
    assert w is not None
    assert w.bricks.rowCount() == 0
    assert w.science_qa.dataset_audit.isEnabled()
    w.science_qa.dataset_audit.click()
    qtbot.waitUntil(
        lambda: w.thread is None and
        "Dataset integrity audit saved" in w.message.text(), timeout=25000)
    text = w.science_qa.result.toPlainText()
    assert "Local datasets: 1; FITS/WCS validated: 1" in text
    assert "DIAGNOSTIC_NOT_ADOPTED" in text
    receipts = list(p.root.glob("dataset-audits/*/audit.json"))
    assert len(receipts) == 1
    payload = json.loads(receipts[0].read_text("utf-8"))
    assert payload["project_id"] == p.project_id
    assert payload["counts"]["validated"] == 1
    assert payload["multi_brick_groups"][0]["region"] == "south"
    assert not list(data.rglob("audit.json"))
    w._go_back_to_projects()
    assert hub.workspace is None



def test_return_to_project_hub_explicitly_shuts_down_child_fits_renderer(
        hub, tmp_path):
    project = create_object_project(
        "UI worker teardown", tmp_path / "ui-worker-teardown")
    with offline_legacy_viewer():
        hub._enter_project(project)
    window = hub.workspace
    assert window is not None
    viewer = window.fits_explorer
    assert not viewer._renderer_closed
    window._go_back_to_projects()
    assert viewer._renderer_closed
    assert not viewer._render_timer.isActive()
    assert hub.workspace is None
    assert hub.isVisible()
    with offline_legacy_viewer():
        hub._enter_project(read_object_project(project.root))
    assert hub.workspace is not None
    assert not hub.workspace.fits_explorer._renderer_closed
    hub.workspace._go_back_to_projects()



def test_hub_removes_selected_stale_shortcut_after_confirmation_without_files(
        hub, monkeypatch, tmp_path):
    from PySide6.QtWidgets import QMessageBox
    a = create_object_project("UGC 2162 usable", tmp_path / "current")
    b = create_object_project("UGC_2162 stale", tmp_path / "stale")
    moved = tmp_path / "stale_moved"
    b.root.rename(moved)
    hub.refresh()
    assert hub.projects.count() == 2
    stale_index = next(i for i in range(hub.projects.count()) if
                       hub.projects.item(i).data(Qt.UserRole)["project_id"] ==
                       b.project_id)
    hub.projects.setCurrentRow(stale_index)
    assert hub.remove_btn.isEnabled()
    assert "Missing" in hub.projects.currentItem().text()
    choices = []

    def confirm(*args, **kwargs):
        choices.append(args)
        return QMessageBox.Cancel if len(choices) == 1 else QMessageBox.Yes

    monkeypatch.setattr(QMessageBox, "question", confirm)
    hub.remove_btn.click()
    assert hub.projects.count() == 2
    assert read_object_project(moved).project_id == b.project_id
    assert len(choices) == 1
    hub.remove_btn.click()
    assert hub.projects.count() == 1
    assert hub.projects.item(0).data(Qt.UserRole)["project_id"] == a.project_id
    assert read_object_project(moved).project_id == b.project_id
    assert read_object_project(a.root).project_id == a.project_id
    assert "No files were deleted" in hub.hint.text()
    assert all("original FITS" not in str(x) or "NOT be deleted" in str(x)
               for x in choices)


def test_hub_remove_available_project_preserves_data_and_can_reimport(
        hub, monkeypatch, tmp_path):
    from PySide6.QtWidgets import QMessageBox
    p = create_object_project("LSB-1", tmp_path / "LSB-1")
    image = p.root / "original.fits.fz"
    image.write_bytes(b"read only original")
    hub.refresh()
    monkeypatch.setattr(
        QMessageBox, "question", lambda *args, **kwargs: QMessageBox.Yes)
    hub.remove_btn.click()
    assert hub.projects.count() == 0
    assert image.read_bytes() == b"read only original"
    assert read_object_project(p.root).project_id == p.project_id
    with offline_legacy_viewer():
        assert hub.import_folder(p.root)
    assert hub.workspace._object_project.project_id == p.project_id
    hub.workspace._go_back_to_projects()
    assert hub.projects.count() == 1


def test_hub_remove_disabled_when_no_visible_selection_and_with_active_workspace(
        hub, monkeypatch, tmp_path):
    from PySide6.QtWidgets import QMessageBox
    p = create_object_project("Valid", tmp_path / "Valid")
    hub.refresh()
    hub.filter.setText("not-the-object")
    assert not hub.remove_btn.isEnabled()
    hub.filter.clear()
    with offline_legacy_viewer():
        hub._enter_project(p)
    assert hub.workspace is not None
    hub._selection_changed()
    assert not hub.remove_btn.isEnabled()
    hub._remove_selected()
    assert recent_projects()[0]["project_id"] == p.project_id
    hub.workspace._go_back_to_projects()
    assert hub.remove_btn.isEnabled()



def test_new_project_wizard_opens_offline_sky_without_network(
        hub, monkeypatch, tmp_path):
    from lsbpipe.gui import project_hub as hub_module
    project = create_object_project(
        "UGC 2162 new", tmp_path / "ugc_new",
        field=SkyField(40.096208, 1.22925, 5.))

    class AcceptedWizard:
        def __init__(self, parent):
            self.project = project

        def exec(self):
            return QDialog.Accepted

    monkeypatch.setattr(hub_module, "CreateProjectDialog", AcceptedWizard)
    with offline_legacy_viewer():
        hub._new_project()
    w = hub.workspace
    assert w is not None
    assert w.preview_tabs.currentWidget() is w.sky
    assert w.inspector_tabs.currentIndex() == 0
    assert not w.sky.url().toString().startswith("https://")
    assert not w.fits_explorer._dataset
    assert w.manifest is None and w.bricks.rowCount() == 0
    assert w.back_to_projects.isVisible()
    target_box = w.ra.parentWidget()
    assert w.back_to_projects.parentWidget() is target_box.parentWidget()
    assert w.back_to_projects.parentWidget() is not target_box
    w._go_back_to_projects()


def test_existing_project_without_fits_opens_survey_tab_offline(
        hub, tmp_path):
    project = create_object_project(
        "Blank target", tmp_path / "blank_target",
        field=SkyField(40.096208, 1.22925, 5.))
    with offline_legacy_viewer():
        hub._enter_project(project)
    w = hub.workspace
    assert w.preview_tabs.currentWidget() is w.sky
    assert w.inspector_tabs.currentIndex() == 0
    assert w.manifest is None
    assert not w.sky.url().toString().startswith("https://")
    w._go_back_to_projects()


def test_project_root_data_directory_allows_background_job_journaling(
        hub, qtbot, tmp_path):
    from lsbpipe.surveys.object_projects import save_object_project
    project = create_object_project(
        "Root data project", tmp_path / "root_data",
        field=SkyField(40.096208, 1.22925, 5.))
    project.data_root = project.root
    save_object_project(project)
    with offline_legacy_viewer():
        hub._enter_project(read_object_project(project.root))
    w = hub.workspace
    finished = []
    w._start(lambda progress: "safe", finished.append,
             operation_name="survey_search")
    qtbot.waitUntil(
        lambda: w.thread is None and finished == ["safe"],
        timeout=12000)
    assert "Cannot journal operation" not in w.project_status.text()
    assert (project.root / "operations").is_dir()
    w._go_back_to_projects()


def test_opt_in_startup_resumes_recent_galaxy_at_saved_scene(
        qtbot, tmp_path, monkeypatch):
    from PySide6.QtCore import QSettings
    import lsbpipe.gui.project_hub as hub_module
    from lsbpipe.photometry.project_resume import save_project_resume

    # Do not write the real OS application preferences in CI.
    settings_path = str(tmp_path / "test-galaxy-startup.ini")
    monkeypatch.setattr(
        hub_module, "QSettings",
        lambda *_: QSettings(settings_path, QSettings.IniFormat))
    project = create_object_project(
        "Leo P", tmp_path / "startup_leo",
        field=SkyField(155.438012, 18.088025, 5.))
    save_project_resume(
        project.root, stage="scene", band="r", location="photometry")
    with offline_legacy_viewer():
        h = ProjectHub()
    qtbot.addWidget(h)
    h.show()
    assert not h.auto_resume.isChecked()
    h._auto_resume_last_project()
    assert h.workspace is None  # preference explicitly opt-in
    h.auto_resume.setChecked(True)
    assert QSettings(settings_path, QSettings.IniFormat).value(
        "startup/auto_resume_last", False, type=bool) is True
    with offline_legacy_viewer():
        h._auto_resume_last_project()
    w = h.workspace
    assert w is not None
    assert w._object_project.project_id == project.project_id
    assert w.preview_tabs.currentWidget() is w.photometry_workspace
    assert w.photometry_workspace.current_navigation() == ("scene", "r")
    w._go_back_to_projects()
    assert h.workspace is None
    assert h.isVisible()
    h.auto_resume.setChecked(False)
