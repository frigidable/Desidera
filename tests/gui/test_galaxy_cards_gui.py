"""Qt regression tests for a single passport in the right-hand inspector."""
from pathlib import Path

import pytest

pytest.importorskip("pytestqt")
pytest.importorskip("PySide6")

from PySide6.QtCore import Qt, QPoint
from PySide6.QtWidgets import QSplitter, QScrollArea

from lsbpipe.gui.galaxy_card_model import GalaxyCardSnapshot
from lsbpipe.gui.galaxy_cards import (
    GalaxyCardsPreview, GalaxyDetailCard, INITIAL_SPLITTER_SIZES,
)


def sample(name, key, *, available=True, aliases=()):
    return GalaxyCardSnapshot(
        project_id=key,
        display_name=name,
        root=Path("/sample") / name,
        available=available,
        legacy=False,
        coordinates="RA 40.096208° · Dec +1.229250° (ICRS)",
        dataset="dr11 / south / 0400p012 / r",
        fits_status="No local FITS link",
        last_opened="2026-11-02",
        aliases=aliases,
    )


def test_existing_left_list_and_single_right_passport(qtbot):
    view = GalaxyCardsPreview()
    qtbot.addWidget(view)
    view.resize(1100, 720)
    view.show()
    view.set_snapshots([sample("Leo P", "leo"), sample("UGC 2162", "ugc")])
    assert isinstance(view.splitter, QSplitter)
    assert tuple(INITIAL_SPLITTER_SIZES) == (480, 540)
    assert view.projects.count() == 2
    assert view.detail.name.text() == "Leo P"
    assert view.findChildren(GalaxyDetailCard) == [view.detail]
    assert view.detail.minimumWidth() <= view.splitter.sizes()[1]
    assert view.detail.preview.width() == view.detail.preview.height()
    assert view.detail.preview.width() in (124, 184)
    assert (view.detail.preview.mapToGlobal(QPoint(0, 0)).x()
            < view.detail.name.mapToGlobal(QPoint(0, 0)).x())
    assert view.detail.preview.size().width() <= 184


def test_right_passport_tracks_selected_left_row(qtbot):
    view = GalaxyCardsPreview()
    qtbot.addWidget(view)
    view.show()
    view.set_snapshots([sample("Leo P", "leo"), sample("UGC 2162", "ugc")])
    with qtbot.waitSignal(view.selectionChanged, timeout=1500) as event:
        view.projects.setCurrentRow(1)
    assert event.args == ["ugc"]
    assert view.detail.snapshot.project_id == "ugc"
    assert "0400p012" in view.detail.dataset.text()


def test_search_name_alias_and_empty_state(qtbot):
    view = GalaxyCardsPreview()
    qtbot.addWidget(view)
    view.show()
    view.set_snapshots([
        sample("Leo P", "leo", aliases=("AGC 208583",)),
        sample("KKR 25", "kkr"),
    ])
    view.search.setText("208583")
    assert not view.projects.item(0).isHidden()
    assert view.projects.item(1).isHidden()
    assert view.selected_id == "leo"
    view.search.setText("nothing like this")
    assert view.selected_id is None
    assert view.detail.snapshot is None
    assert not view.detail.open_button.isEnabled()


def test_open_action_emits_immutable_id_without_side_effects(qtbot):
    view = GalaxyCardsPreview()
    qtbot.addWidget(view)
    view.show()
    view.set_snapshots([sample("Leo P", "stable-uuid")])
    with qtbot.waitSignal(view.openRequested, timeout=1500) as event:
        qtbot.mouseClick(view.detail.open_button, Qt.LeftButton)
    assert event.args == ["stable-uuid"]


def test_missing_project_shows_relink_without_faking_data(qtbot):
    card = GalaxyDetailCard()
    qtbot.addWidget(card)
    card.resize(420, 640)
    card.show()
    card.set_snapshot(sample("KKR 25", "missing", available=False))
    assert "Relink" in card.open_button.text()
    assert "Missing" in card.fits.text()
    assert "No image yet" in card.preview.text()
    assert "No local FITS link" in card.fits.text()
    card.set_snapshot(None)
    assert not card.open_button.isEnabled()
    assert card.snapshot is None


def test_right_inspector_remains_resizable_and_scrollable(qtbot):
    view = GalaxyCardsPreview()
    qtbot.addWidget(view)
    view.resize(760, 560)
    view.show()
    view.set_snapshots([sample("Leo P", "leo")])
    view.splitter.setSizes([300, 420])
    assert view.detail.width() >= 300
    assert view.detail.findChild(QScrollArea, "galaxy_card_scroll") is not None
    assert not view.detail.safety.isVisible()
    assert "scientific QA" in view.detail.open_button.toolTip()
    assert not view.detail.details_panel.isVisible()



def test_empty_galaxy_passport_removes_redundant_explanations(qtbot):
    from types import SimpleNamespace
    card=GalaxyDetailCard()
    qtbot.addWidget(card)
    card.resize(600,650)
    card.show()
    snapshot=sample("UGC 2162", "ugc")
    # A selected object with no confirmed catalogue has one short status,
    # with complete research help in the accessible tooltip.
    card.set_snapshot(snapshot)
    assert not card.catalogue_help.isVisible()
    assert "Why identify this galaxy" in card.catalogue.toolTip() or (
        "Why identify this galaxy" in card.catalogue.accessibleDescription())
    assert card.preview.width()==124
    assert card.preview.height()==124
    # A source-bearing snapshot remains fully visible and keeps the facts.
    assert "0400p012" in card.dataset.text()
    assert "No local FITS link" in card.fits.text()
    assert card.fits.isVisible()
