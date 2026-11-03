"""Production Project Hub integration of the compact right-hand passport.

Uses a separate per-test project registry, no network, no science jobs.
"""
import pytest

pytest.importorskip("pytestqt")
pytest.importorskip("PySide6.QtWebEngineWidgets")

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QSplitter

from scripts.offline_gui_fixture import offline_legacy_viewer
from lsbpipe.gui.project_hub import ProjectHub
from lsbpipe.gui.galaxy_cards import GalaxyDetailCard
from lsbpipe.surveys.object_projects import (
    create_object_project, read_object_project, save_object_project,
)
from lsbpipe.surveys.legacy import SkyField


@pytest.fixture
def hub(qtbot):
    with offline_legacy_viewer():
        window = ProjectHub()
    qtbot.addWidget(window)
    window.show()
    return window


def test_passport_replaces_right_labels_without_replacing_left_list(hub, tmp_path):
    p = create_object_project(
        "UGC 2162", tmp_path / "UGC_2162",
        field=SkyField(40.096208, 1.22925, 5.),
    )
    original = p.root / "test_image.fits"
    original.write_bytes(b"original observations must remain untouched")
    hub.refresh()

    assert hub.projects.objectName() == "hub_project_list"
    assert hub.projects.count() == 1
    assert isinstance(hub.splitter, QSplitter)
    assert len(hub.findChildren(GalaxyDetailCard)) == 1
    assert hub.galaxy_card.snapshot.project_id == p.project_id
    assert hub.galaxy_card.name.text() == "UGC 2162"
    assert "40.096208" in hub.galaxy_card.coordinates.text()
    assert hub.galaxy_card.preview.text() == "No image yet"
    assert hub.galaxy_card.open_button.isVisible()
    assert not hub.continue_btn.isVisible()
    assert not hub.remove_btn.isVisible()
    assert hub.tools_menu_btn.isVisible()
    assert hub.remove_action.isEnabled()
    assert original.read_bytes() == b"original observations must remain untouched"


def test_passport_selection_updates_with_existing_list_and_search(hub, tmp_path):
    a = create_object_project("Leo P", tmp_path / "leo")
    b = create_object_project("KKR 25", tmp_path / "kkr")
    hub.refresh()
    hub.filter.setText("leo")
    assert hub.galaxy_card.snapshot.project_id == a.project_id
    assert hub.galaxy_card.name.text() == "Leo P"
    hub.filter.setText("no-such-object")
    assert hub.galaxy_card.snapshot is None
    assert not hub.galaxy_card.open_button.isEnabled()
    hub.filter.clear()
    hub.projects.setCurrentRow(next(
        i for i in range(hub.projects.count())
        if hub.projects.item(i).data(Qt.UserRole)["project_id"] == b.project_id
    ))
    assert hub.galaxy_card.snapshot.project_id == b.project_id


def test_passport_shows_dataset_qa_caveat_and_survives_legacy_widgets(
        hub, tmp_path):
    p = create_object_project(
        "Leo P", tmp_path / "leo",
        field=SkyField(155.438012, 18.088025, 4.),
    )
    data = p.root / "data"
    data.mkdir()
    p.data_root = data
    p.selected_dataset = ("dr11", "south", "1554p180", "r")
    save_object_project(p)
    hub.refresh()
    assert "dr11 / south / 1554p180 / r" in hub.galaxy_card.dataset.text()
    assert "No original FITS found" in hub.galaxy_card.fits.text()
    assert "get original FITS" in hub.galaxy_card.open_button.text()
    assert "integrity not checked" in hub.detail_status.text()
    assert read_object_project(p.root).selected_dataset == p.selected_dataset


def test_card_open_uses_same_project_action_and_stable_selection(
        hub, tmp_path, monkeypatch, qtbot):
    p = create_object_project("Leo P", tmp_path / "leo")
    hub.refresh()
    hits = []
    monkeypatch.setattr(hub, "_open_selected", lambda: hits.append(
        hub.projects.currentItem().data(Qt.UserRole)["project_id"]
    ))
    qtbot.mouseClick(hub.galaxy_card.open_button, Qt.LeftButton)
    assert hits == [p.project_id]
    hub._open_card_project("different-project-id")
    assert hits == [p.project_id]
    hub._open_card_project(p.project_id)
    assert hits == [p.project_id, p.project_id]


def test_missing_project_relinks_instead_of_fabricating_fits_status(
        hub, tmp_path):
    p = create_object_project("Missing", tmp_path / "moved")
    moved = tmp_path / "elsewhere"
    p.root.rename(moved)
    hub.refresh()
    assert hub.galaxy_card.snapshot.project_id == p.project_id
    assert "Relink" in hub.galaxy_card.open_button.text()
    assert "Location missing" in hub.galaxy_card.fits.text()
    assert read_object_project(moved).project_id == p.project_id
