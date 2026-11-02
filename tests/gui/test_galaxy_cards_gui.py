"""GUI smoke checks for the isolated Galaxy Cards module."""
import pytest

pytest.importorskip("pytestqt")
pytest.importorskip("PySide6")

from pathlib import Path
from PySide6.QtCore import Qt
from lsbpipe.gui.galaxy_card_model import GalaxyCardSnapshot
from lsbpipe.gui.galaxy_cards import GalaxyCardGallery


def sample(name, key, *, available=True):
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
    )


def test_gallery_filters_and_preserves_order(qtbot):
    gallery = GalaxyCardGallery()
    qtbot.addWidget(gallery)
    gallery.resize(940, 700)
    gallery.show()
    gallery.set_snapshots([sample("Leo P", "leo"), sample("UGC 2162", "ugc")])
    assert [s.project_id for s in gallery.filtered_snapshots()] == ["leo", "ugc"]
    gallery.search.setText("ugc")
    assert [s.project_id for s in gallery.filtered_snapshots()] == ["ugc"]
    assert gallery.counter.text() == "1 / 2 projects"
    gallery.search.setText("unknown")
    assert not gallery.filtered_snapshots()
    assert gallery.empty.isVisible()
    gallery.search.clear()
    assert len(gallery.filtered_snapshots()) == 2


def test_card_click_selects_and_enter_requests_open(qtbot):
    gallery = GalaxyCardGallery()
    qtbot.addWidget(gallery)
    gallery.resize(800, 640)
    gallery.show()
    gallery.set_snapshots([sample("Leo P", "leo")])
    card = gallery._cards["leo"]
    with qtbot.waitSignal(gallery.selectionChanged, timeout=1500) as selected:
        qtbot.mouseClick(card, Qt.LeftButton)
    assert selected.args == ["leo"]
    assert gallery.selected_id == "leo"
    with qtbot.waitSignal(gallery.openRequested, timeout=1500) as opened:
        qtbot.keyClick(card, Qt.Key_Return)
    assert opened.args == ["leo"]


def test_missing_object_is_readable_without_thumbnail_or_network(qtbot):
    gallery = GalaxyCardGallery()
    qtbot.addWidget(gallery)
    gallery.show()
    gallery.set_snapshots([sample("KKR 25", "missing", available=False)])
    card = gallery._cards["missing"]
    assert "Missing" in card.availability.text()
    assert card.preview.text() == "No cached image"
    assert "KKR 25" in card.accessibleName()
