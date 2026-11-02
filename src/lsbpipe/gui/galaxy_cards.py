"""Standalone, read-only PySide6 gallery for Galaxy Projects.

The existing ProjectHub remains the production launcher. This widget is
intentionally developed separately so it can be integrated after the active
GUI / data-manager branches have settled.
"""
from __future__ import annotations

import sys

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import (
    QApplication, QFrame, QGridLayout, QHBoxLayout, QLabel, QLineEdit,
    QMainWindow, QScrollArea, QSizePolicy, QVBoxLayout, QWidget,
)

from .galaxy_card_model import GalaxyCardSnapshot, snapshots_from_records
from lsbpipe.surveys.object_projects import recent_projects


CARD_MIN_WIDTH = 256

_CARD_STYLE = """
QFrame#galaxy_card {
    border: 1px solid palette(mid);
    border-radius: 8px;
    background: palette(base);
}
QFrame#galaxy_card[selectedCard="true"] {
    border: 2px solid palette(highlight);
}
QLabel#card_title { font-size: 17px; font-weight: 600; }
QLabel#card_meta { color: palette(text); }
QLabel#card_preview {
    background: palette(alternate-base);
    border: 0px;
    border-radius: 5px;
}
"""


class GalaxyCard(QFrame):
    """Keyboard-accessible visual representation of one saved project."""

    selected = Signal(str)
    activated = Signal(str)

    def __init__(self, snapshot: GalaxyCardSnapshot, parent=None):
        super().__init__(parent)
        self.snapshot = snapshot
        self.setObjectName("galaxy_card")
        self.setProperty("selectedCard", False)
        self.setFocusPolicy(Qt.StrongFocus)
        self.setCursor(Qt.PointingHandCursor)
        self.setAccessibleName(f"Galaxy project: {snapshot.display_name}")
        self.setToolTip(str(snapshot.root))
        self.setMinimumWidth(CARD_MIN_WIDTH)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.setStyleSheet(_CARD_STYLE)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(8)

        self.preview = QLabel("No cached image")
        self.preview.setObjectName("card_preview")
        self.preview.setAlignment(Qt.AlignCenter)
        self.preview.setFixedHeight(124)
        self.preview.setWordWrap(True)
        layout.addWidget(self.preview)
        self._pixmap = None
        if snapshot.thumbnail is not None:
            candidate = QPixmap(str(snapshot.thumbnail))
            if not candidate.isNull():
                self._pixmap = candidate
                self._update_preview()

        self.title = QLabel(snapshot.display_name)
        self.title.setObjectName("card_title")
        self.title.setWordWrap(True)
        layout.addWidget(self.title)

        self.coordinates = QLabel(snapshot.coordinates)
        self.coordinates.setObjectName("card_meta")
        self.coordinates.setWordWrap(True)
        self.coordinates.setTextInteractionFlags(Qt.TextSelectableByMouse)
        layout.addWidget(self.coordinates)

        self.dataset = QLabel(f"Dataset: {snapshot.dataset}")
        self.dataset.setObjectName("card_dataset")
        self.dataset.setWordWrap(True)
        layout.addWidget(self.dataset)

        status = "On disk" if snapshot.available else "Missing · relink required"
        self.availability = QLabel(f"{status}  ·  {snapshot.fits_status}")
        self.availability.setObjectName("card_availability")
        self.availability.setWordWrap(True)
        layout.addWidget(self.availability)
        if snapshot.note:
            self.note = QLabel(snapshot.note)
            self.note.setObjectName("card_note")
            self.note.setWordWrap(True)
            layout.addWidget(self.note)

    def _update_preview(self) -> None:
        if self._pixmap is not None:
            self.preview.setPixmap(self._pixmap.scaled(
                self.preview.size(), Qt.KeepAspectRatio,
                Qt.SmoothTransformation,
            ))

    def resizeEvent(self, event) -> None:
        super().resizeEvent(event)
        self._update_preview()

    def set_selected(self, value: bool) -> None:
        self.setProperty("selectedCard", bool(value))
        self.style().unpolish(self)
        self.style().polish(self)
        self.update()

    def mousePressEvent(self, event) -> None:
        if event.button() == Qt.LeftButton:
            self.selected.emit(self.snapshot.project_id)
            self.setFocus(Qt.MouseFocusReason)
        super().mousePressEvent(event)

    def mouseDoubleClickEvent(self, event) -> None:
        if event.button() == Qt.LeftButton:
            self.activated.emit(self.snapshot.project_id)
        super().mouseDoubleClickEvent(event)

    def keyPressEvent(self, event) -> None:
        if event.key() in (Qt.Key_Return, Qt.Key_Enter):
            self.activated.emit(self.snapshot.project_id)
            event.accept()
            return
        if event.key() == Qt.Key_Space:
            self.selected.emit(self.snapshot.project_id)
            event.accept()
            return
        super().keyPressEvent(event)


class GalaxyCardGallery(QWidget):
    """Responsive, searchable card grid; no launcher or project mutations."""

    selectionChanged = Signal(str)
    openRequested = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("galaxy_card_gallery")
        self._snapshots: list[GalaxyCardSnapshot] = []
        self._cards: dict[str, GalaxyCard] = {}
        self.selected_id: str | None = None
        self._columns = 0

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(10)

        header = QHBoxLayout()
        label = QLabel("Galaxy library")
        label.setStyleSheet("font-size: 17px; font-weight: 600")
        header.addWidget(label)
        header.addStretch()
        self.counter = QLabel("0 projects")
        header.addWidget(self.counter)
        outer.addLayout(header)

        self.search = QLineEdit()
        self.search.setObjectName("galaxy_cards_search")
        self.search.setClearButtonEnabled(True)
        self.search.setPlaceholderText("Search name, aliases, coordinates or folder…")
        self.search.textChanged.connect(self._reflow)
        outer.addWidget(self.search)

        self.scroll = QScrollArea()
        self.scroll.setObjectName("galaxy_cards_scroll")
        self.scroll.setWidgetResizable(True)
        self.scroll.setFrameShape(QFrame.NoFrame)
        self.surface = QWidget()
        self.grid = QGridLayout(self.surface)
        self.grid.setContentsMargins(0, 0, 0, 0)
        self.grid.setHorizontalSpacing(12)
        self.grid.setVerticalSpacing(12)
        self.grid.setAlignment(Qt.AlignTop)
        self.empty = QLabel("No matching galaxy projects")
        self.empty.setObjectName("galaxy_cards_empty")
        self.empty.setAlignment(Qt.AlignCenter)
        self.empty.setWordWrap(True)
        self.scroll.setWidget(self.surface)
        outer.addWidget(self.scroll, 1)

    def set_snapshots(self, snapshots: list[GalaxyCardSnapshot]) -> None:
        for card in self._cards.values():
            self.grid.removeWidget(card)
            card.deleteLater()
        self._cards.clear()
        self._snapshots = list(snapshots)
        for snapshot in self._snapshots:
            card = GalaxyCard(snapshot, self.surface)
            card.selected.connect(self.select_project)
            card.activated.connect(self.openRequested.emit)
            self._cards[snapshot.project_id] = card
        if self.selected_id not in self._cards:
            self.selected_id = None
        self._reflow()

    def select_project(self, project_id: str) -> None:
        if project_id not in self._cards:
            return
        self.selected_id = project_id
        for key, card in self._cards.items():
            card.set_selected(key == project_id)
        self.selectionChanged.emit(project_id)

    def filtered_snapshots(self) -> list[GalaxyCardSnapshot]:
        return [s for s in self._snapshots if s.matches(self.search.text())]

    def _reflow(self, *_args) -> None:
        while self.grid.count():
            item = self.grid.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.hide()
        items = self.filtered_snapshots()
        self.counter.setText(f"{len(items)} / {len(self._snapshots)} projects")
        width = max(self.scroll.viewport().width() - 16, CARD_MIN_WIDTH)
        self._columns = max(1, width // (CARD_MIN_WIDTH + 12))
        if not items:
            self.grid.addWidget(self.empty, 0, 0)
            self.empty.show()
            return
        self.empty.hide()
        for i, item in enumerate(items):
            card = self._cards[item.project_id]
            self.grid.addWidget(card, i // self._columns, i % self._columns)
            card.show()
            card.set_selected(item.project_id == self.selected_id)

    def resizeEvent(self, event) -> None:
        super().resizeEvent(event)
        columns = max(
            1, (self.scroll.viewport().width() - 16) // (CARD_MIN_WIDTH + 12)
        )
        if columns != self._columns:
            self._reflow()


def main() -> None:
    """Opt-in standalone visual smoke preview, leaving ProjectHub unchanged."""
    app = QApplication.instance() or QApplication(sys.argv)
    window = QMainWindow()
    window.setWindowTitle("DESIDERA · Galaxy Cards (isolated preview)")
    window.resize(1040, 760)
    frame = QWidget()
    layout = QVBoxLayout(frame)
    gallery = GalaxyCardGallery()
    gallery.set_snapshots(snapshots_from_records(recent_projects()))
    layout.addWidget(gallery)
    hint = QLabel(
        "Read-only preview. Open projects from the regular DESIDERA launcher; "
        "this view does not launch analyses or download data."
    )
    hint.setWordWrap(True)
    layout.addWidget(hint)
    window.setCentralWidget(frame)
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
