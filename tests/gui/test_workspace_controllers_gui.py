"""P0.4 controller extraction: real Qt widget behavior, no network."""
from __future__ import annotations

import pytest

pytest.importorskip("pytestqt")
pytest.importorskip("PySide6.QtWebEngineWidgets")

from PySide6.QtWidgets import QCheckBox, QInputDialog
from scripts.offline_gui_fixture import offline_legacy_viewer
from lsbpipe.gui.survey_explorer import SurveyExplorer
from lsbpipe.gui.target_controller import TargetController, TargetTextField
from lsbpipe.gui.survey_data_controller import SurveyDataController
from lsbpipe.surveys.legacy import BrickHit, SkyField, build_manifest


@pytest.fixture
def workspace(qtbot):
    with offline_legacy_viewer():
        w = SurveyExplorer(restore_session=False, launch_remote=False)
    qtbot.addWidget(w)
    w.show()
    yield w
    if w.thread is not None:
        w._request_stop()
        qtbot.waitUntil(lambda: w.thread is None, timeout=12000)


def test_target_controller_is_real_owner_and_wrappers_preserve_strict_icrs(
        workspace):
    w = workspace
    assert isinstance(w._target_controller, TargetController)
    assert isinstance(w.ra, TargetTextField)
    assert w._field() == w._target_controller.field()
    w.ra.setText("not-a-position")
    with pytest.raises(ValueError):
        w._field()
    w.ra.setValue(40.096208)
    w.dec.setValue(1.229250)
    assert w._field().ra_deg == pytest.approx(40.096208)
    assert w._preferred_viewer_region("dr9") == "south"
    w.dec.setValue(35.)
    assert w._preferred_viewer_region("dr9") == "north"
    w.dec.setValue(32.375)
    assert w._preferred_viewer_region("dr9") is None


def test_data_controller_preserves_hidden_selected_products_and_quartet(
        workspace, monkeypatch):
    w = workspace
    assert isinstance(w._data_controller, SurveyDataController)
    hit = BrickHit("dr9", "south", "0401p012",
                   40., 40.2, 1., 1.4, {"g": 1, "r": 1})
    field = SkyField(40.096208, 1.22925, 5.)
    w.hits = [hit]
    w.manifest = build_manifest(
        field, [hit], kinds=("image", "invvar", "nexp", "maskbits"),
        bands=("g", "r"))
    w.bricks.setRowCount(1)
    w.bricks.selectRow(0)
    w._fill_brick_filters()
    w._fill_files()
    assert len(w._selected()) == len(w.manifest["products"])
    all_paths = w._selected()
    w.band_filter.setCurrentText("r")
    w._filter_products()
    assert w._selected() == all_paths, "visual filter must not alter selection"
    w._set_visible_selection(False)
    assert len(w._selected()) < len(all_paths)
    monkeypatch.setattr(QInputDialog, "getItem",
                        lambda *args, **kwargs: ("r", True))
    w.select_science_quartet()
    selected = [
        p for p in w.manifest["products"]
        if p["relative_path"] in w._selected()]
    assert len(selected) == 4
    assert {p["kind"] for p in selected} == {
        "image", "invvar", "nexp", "maskbits"}
    assert all(p["band"] in ("r", None) for p in selected)
    w._clear_all_selection()
    assert w._selected() == set()


def test_project_hub_can_create_three_sequential_independent_workspaces(
        qtbot, tmp_path):
    from lsbpipe.gui.project_hub import ProjectHub
    from lsbpipe.surveys.object_projects import (
        create_object_project, read_object_project)
    projects = [
        create_object_project("Candidate " + str(i), tmp_path / ("candidate" + str(i)))
        for i in range(3)
    ]
    with offline_legacy_viewer():
        hub = ProjectHub()
    qtbot.addWidget(hub)
    hub.show()
    seen = set()
    for project in projects:
        with offline_legacy_viewer():
            hub._enter_project(read_object_project(project.root))
        viewer = hub.workspace.fits_explorer
        assert not viewer._renderer_closed
        seen.add(hub.workspace._object_project.project_id)
        hub.workspace._go_back_to_projects()
        assert viewer._renderer_closed
        assert hub.workspace is None
        assert hub.isVisible()
    assert len(seen) == 3



def test_science_context_controller_freezes_project_and_rejects_old_selection(
        workspace, tmp_path, monkeypatch):
    from PySide6.QtWidgets import QTableWidgetItem
    from lsbpipe.gui.science_context_controller import ScienceContextController
    from lsbpipe.surveys.object_projects import create_object_project
    w = workspace
    project = create_object_project(
        "Science candidate", tmp_path / "science_candidate",
        field=SkyField(40.096208, 1.22925, 5.))
    w._object_project = project
    assert isinstance(w._science_context_controller, ScienceContextController)
    w.ra.setValue(project.field.ra_deg)
    w.dec.setValue(project.field.dec_deg)
    w.field_width.setValue(project.field.width_arcmin)
    hit = BrickHit("dr9", "south", "0401p012",
                   40., 40.2, 1., 1.4, {"r": 1})
    w.hits = [hit]
    w.bricks.setRowCount(1)
    w.bricks.setItem(0, 0, QTableWidgetItem("DR9"))
    w.bricks.selectRow(0)
    w.science_qa.band.setCurrentText("r")
    w.science_qa.directory.setText(str(project.root / "data"))
    monkeypatch.setattr(w, "_brick_result_matches_target", lambda: True)
    context = w._current_science_context("census")
    assert context is not None
    assert context.project_id == project.project_id
    ticket = w._qa_generation
    assert w._science_ticket_current(context=context, ticket=ticket)
    assert w._science_context_controller.ticket_current(
        context=context, ticket=ticket)
    assert w._qa_context_current({
        "field": project.field,
        "source": f"DR9 south · 0401p012 · r"})
    w.target_radius.setValue(w.target_radius.value() + 1)
    assert not w._science_ticket_current(context=context, ticket=ticket)
    w.target_radius.setValue(context.target_radius_arcsec)
    w.ra.setValue(40.111111)
    assert not w._qa_context_current({
        "field": project.field,
        "source": f"DR9 south · 0401p012 · r"})
    assert not w._science_ticket_current(context=context, ticket=ticket)
