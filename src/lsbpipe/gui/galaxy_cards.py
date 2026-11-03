"""Compact, read-only galaxy passport for the existing Project Hub inspector.

Production ProjectHub is intentionally unchanged. This module provides a
right-hand detail widget and a standalone two-pane integration preview.
No catalogue requests, science runs, FITS decoding or manifest writes occur.
"""
from __future__ import annotations

import sys

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import (
    QApplication, QFrame, QHBoxLayout, QLabel, QLineEdit, QListWidget,
    QListWidgetItem, QMainWindow, QMessageBox, QPushButton, QScrollArea,
    QSizePolicy, QSplitter, QVBoxLayout, QWidget,
)

from .galaxy_card_model import GalaxyCardSnapshot, snapshots_from_records
from lsbpipe.surveys.object_projects import recent_projects


# Widening the previous ~40% right-hand inspector by about 30% yields ~52%.
INITIAL_SPLITTER_SIZES = (480, 540)

_INSPECTOR_STYLE = """
QFrame#galaxy_detail_card {
    border: 1px solid palette(mid);
    border-radius: 8px;
}
QLabel#galaxy_card_preview {
    background: palette(alternate-base);
    border: 1px solid palette(mid);
    border-radius: 6px;
}
QLabel#galaxy_card_name { font-size: 20px; font-weight: 600; }
QLabel#galaxy_card_section { font-size: 11px; font-weight: 600; }
"""


def _metadata(label: str, *, selectable: bool = False) -> QLabel:
    widget = QLabel(label)
    widget.setWordWrap(True)
    if selectable:
        widget.setTextInteractionFlags(Qt.TextSelectableByMouse)
    return widget


class GalaxyDetailCard(QFrame):
    """One selected object, sized for the CURRENT right-hand inspector.

    Host integrations can connect openRequested to the existing ProjectHub
    open/relink behavior by immutable project ID. Rendering never mutates a
    project, and a linked FITS path is explicitly not a scientific QA pass.
    """

    openRequested = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("galaxy_detail_card")
        self.setFrameShape(QFrame.StyledPanel)
        self.setMinimumWidth(300)
        self.setStyleSheet(_INSPECTOR_STYLE)
        self.snapshot: GalaxyCardSnapshot | None = None
        self._pixmap: QPixmap | None = None

        container = QVBoxLayout(self)
        container.setContentsMargins(12, 12, 12, 12)
        container.setSpacing(8)

        kicker = QLabel("SELECTED GALAXY")
        kicker.setObjectName("galaxy_card_section")
        container.addWidget(kicker)

        # Preview and passport share one horizontally arranged, scrollable
        # content surface. Keep the image square at every inspector width.
        scroll = QScrollArea(self)
        scroll.setObjectName("galaxy_card_scroll")
        scroll.setFrameShape(QFrame.NoFrame)
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        inner = QWidget()
        detail = QHBoxLayout(inner)
        detail.setContentsMargins(0, 0, 0, 0)
        detail.setSpacing(12)

        self.preview_holder = QWidget()
        self.preview_holder.setObjectName("galaxy_card_preview_column")
        preview_col = QVBoxLayout(self.preview_holder)
        preview_col.setContentsMargins(0, 0, 0, 0)
        preview_col.setSpacing(5)

        self.preview = QLabel("No preview")
        self.preview.setObjectName("galaxy_card_preview")
        self.preview.setAlignment(Qt.AlignCenter)
        self.preview.setWordWrap(True)
        self.preview.setTextFormat(Qt.PlainText)
        self.preview.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        preview_col.addWidget(self.preview)

        self.preview_note = _metadata("")
        self.preview_note.setObjectName("galaxy_card_preview_note")
        self.preview_note.setTextFormat(Qt.PlainText)
        preview_col.addWidget(self.preview_note)
        preview_col.addStretch(1)
        detail.addWidget(self.preview_holder, 0, Qt.AlignTop)

        metadata = QWidget()
        metadata.setObjectName("galaxy_card_metadata_column")
        passport = QVBoxLayout(metadata)
        passport.setContentsMargins(0, 0, 0, 0)
        passport.setSpacing(8)

        self.name = QLabel("Choose a galaxy")
        self.name.setObjectName("galaxy_card_name")
        self.name.setWordWrap(True)
        passport.addWidget(self.name)

        self.coordinates = _metadata("Select a project from the list.")
        self.coordinates.setObjectName("galaxy_card_coordinates")
        self.coordinates.setTextInteractionFlags(Qt.TextSelectableByMouse)
        passport.addWidget(self.coordinates)

        knowledge_header = QLabel("CATALOGUE IDENTITY")
        knowledge_header.setObjectName("galaxy_card_section")
        passport.addWidget(knowledge_header)
        self.catalogue = _metadata("No user-confirmed catalogue identity")
        self.catalogue.setObjectName("galaxy_card_catalogue")
        self.catalogue.setTextFormat(Qt.PlainText)
        passport.addWidget(self.catalogue)
        self.catalogue_provenance = _metadata("", selectable=True)
        self.catalogue_provenance.setObjectName(
            "galaxy_card_catalogue_provenance")
        self.catalogue_provenance.setTextFormat(Qt.PlainText)
        passport.addWidget(self.catalogue_provenance)

        observation = QLabel("OBSERVATIONS")
        observation.setObjectName("galaxy_card_section")
        passport.addWidget(observation)
        self.dataset = _metadata("No dataset selected", selectable=True)
        self.dataset.setObjectName("galaxy_card_dataset")
        passport.addWidget(self.dataset)
        self.fits = _metadata("FITS: not inspected")
        self.fits.setObjectName("galaxy_card_fits")
        passport.addWidget(self.fits)

        research = QLabel("PROJECT INFORMATION")
        research.setObjectName("galaxy_card_section")
        passport.addWidget(research)
        self.aliases = _metadata("")
        self.aliases.setObjectName("galaxy_card_aliases")
        passport.addWidget(self.aliases)
        self.last_opened = _metadata("")
        self.last_opened.setObjectName("galaxy_card_last_opened")
        passport.addWidget(self.last_opened)
        self.path = _metadata("", selectable=True)
        self.path.setObjectName("galaxy_card_path")
        passport.addWidget(self.path)
        self.note = _metadata("")
        self.note.setObjectName("galaxy_card_note")
        passport.addWidget(self.note)
        passport.addStretch(1)
        detail.addWidget(metadata, 1)

        scroll.setWidget(inner)
        container.addWidget(scroll, 1)

        self.open_button = QPushButton("Open project")
        self.open_button.setObjectName("galaxy_card_open")
        self.open_button.setEnabled(False)
        self.open_button.clicked.connect(self._request_open)
        container.addWidget(self.open_button)
        self.safety = _metadata(
            "Opening does not download FITS or approve previous scientific QA.")
        self.safety.setObjectName("galaxy_card_safety")
        container.addWidget(self.safety)
        self._set_preview_side(184)
        self.set_snapshot(None)

    def _set_preview_side(self, side: int):
        """Square preview, responsive down to the 300 px narrow inspector."""
        if self.preview.width() == side and self.preview.height() == side:
            return
        self.preview.setFixedSize(side, side)
        self.preview_holder.setFixedWidth(side)
        if self._pixmap is not None:
            self._update_preview()

    def set_snapshot(self, snapshot: GalaxyCardSnapshot | None) -> None:
        self.snapshot = snapshot
        self._pixmap = None
        self.open_button.setEnabled(snapshot is not None)
        if snapshot is None:
            self.name.setText("Choose a galaxy")
            self.coordinates.setText("Select a project from the list.")
            self.preview.setPixmap(QPixmap())
            self.preview.setText("No galaxy selected")
            self.preview_note.hide()
            self.dataset.clear()
            self.fits.clear()
            self.aliases.hide()
            self.catalogue.setText("No galaxy selected")
            self.catalogue_provenance.hide()
            self.last_opened.hide()
            self.path.clear()
            self.note.hide()
            return

        self.name.setText(snapshot.display_name)
        self.coordinates.setText(snapshot.coordinates)
        self.dataset.setText(f"Dataset: {snapshot.dataset}")
        catalogues = snapshot.catalogue_identities or (
            (snapshot.catalogue_identity,) if snapshot.catalogue_identity else ())
        self.catalogue_provenance.hide()
        if snapshot.knowledge_status == "EVIDENCE_ERROR":
            self.catalogue.setText("Catalogue evidence requires inspection")
            self.catalogue_provenance.setText(
                "Stored catalogue evidence is invalid or unreadable: "
                + snapshot.knowledge_error)
            self.catalogue_provenance.show()
        elif catalogues:
            self.catalogue.setText("\n".join(
                f"{entry.provider}: {entry.main_id} · "
                "user-confirmed association"
                + (f" · type: {entry.object_type}" if entry.object_type else "")
                for entry in catalogues))
            details = []
            for entry in catalogues:
                difference = (
                    f"Separation from project centre: "
                    f"{entry.project_offset_arcsec:.3f} arcsec. "
                    if entry.project_offset_arcsec is not None else
                    "Project position unavailable for comparison. ")
                note = (
                    f"{entry.provider} ICRS: {entry.ra_deg:.8f}°, "
                    f"{entry.dec_deg:+.8f}°. "
                    + difference
                    + f"Evidence: {entry.confirmation_sha256[:12]}. ")
                if entry.provider == "NED" and entry.redshift is not None:
                    note += (f"NED catalogue heliocentric z: "
                             f"{entry.redshift:.7g} "
                             f"(refcode: {entry.redshift_refcode or 'unspecified'}; "
                             "unadopted). ")
                details.append(note)
            if snapshot.provider_offset_arcsec is not None:
                details.append(
                    f"SIMBAD–NED positions differ by "
                    f"{snapshot.provider_offset_arcsec:.3f} arcsec. "
                    "This is a descriptive comparison, not an automatic "
                    "source reconciliation. ")
            details.append(
                "Project coordinates unchanged. Catalogue associations and "
                "literature values are not adopted scientific measurements.")
            self.catalogue_provenance.setText("\n".join(details))
            self.catalogue_provenance.show()
        elif snapshot.knowledge_status == "UNVERIFIED_CANDIDATES_ONLY":
            self.catalogue.setText(
                "SIMBAD candidates stored · no identity confirmed"
                + (f" · NED: {snapshot.ned_status}"
                   if snapshot.ned_status != "NOT_QUERIED" else ""))
        else:
            self.catalogue.setText(
                "No user-confirmed catalogue identity"
                + (f" · NED: {snapshot.ned_status}"
                   if snapshot.ned_status != "NOT_QUERIED" else ""))
        project_state = (
            "On disk" if snapshot.available else "Missing · relink required"
        )
        self.fits.setText(f"{project_state} · {snapshot.fits_status}")
        self.aliases.setVisible(bool(snapshot.aliases))
        self.aliases.setText("Also known as: " + ", ".join(snapshot.aliases))
        self.last_opened.setVisible(bool(snapshot.last_opened))
        self.last_opened.setText("Last opened: " + snapshot.last_opened)
        self.path.setText("Project folder: " + str(snapshot.root))
        self.note.setVisible(bool(snapshot.note))
        self.note.setText(snapshot.note)
        self.open_button.setText(
            "Relink project…" if not snapshot.available else "Open project"
        )
        self.preview.setPixmap(QPixmap())
        self.preview.setText("No cached local preview")
        self.preview_note.setVisible(snapshot.thumbnail is not None)
        if snapshot.preview_source == "REMOTE_LEGACY_RGB":
            self.preview_note.setText(
                "Legacy Surveys RGB · display only · cached")
        elif snapshot.preview_status == "READY_RECENT_STAT_NOT_FITS_REHASHED":
            self.preview_note.setText("Local FITS · display only")
        else:
            self.preview_note.setText(
                "Local preview · unverified provenance")
        if snapshot.thumbnail is not None:
            image = QPixmap(str(snapshot.thumbnail))
            if (image.isNull()
                    and snapshot.preview_source == "REMOTE_LEGACY_RGB"
                    and snapshot.local_thumbnail is not None):
                image = QPixmap(str(snapshot.local_thumbnail))
                self.preview_note.setText(
                    "Remote RGB unreadable · local FITS fallback · display only")
            if not image.isNull():
                self._pixmap = image
                self._update_preview()

    def _update_preview(self) -> None:
        if self._pixmap is not None:
            self.preview.setPixmap(self._pixmap.scaled(
                self.preview.size(), Qt.KeepAspectRatio,
                Qt.SmoothTransformation,
            ))

    def resizeEvent(self, event) -> None:
        super().resizeEvent(event)
        self._set_preview_side(124 if self.width() < 540 else 184)
        self._update_preview()

    def _request_open(self) -> None:
        if self.snapshot is not None:
            self.openRequested.emit(self.snapshot.project_id)


class GalaxyCardsPreview(QWidget):
    """Standalone two-pane preview mirroring the CURRENT Galaxy Projects UI.

    Left: an ordinary searchable project list. Right: a single compact card.
    The splitter defaults near 48/52, remains user-resizable and never
    creates a second card grid or changes the production launcher.
    """

    selectionChanged = Signal(str)
    openRequested = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("galaxy_cards_preview")
        self._snapshots: dict[str, GalaxyCardSnapshot] = {}
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(8)

        self.splitter = QSplitter(Qt.Horizontal)
        self.splitter.setObjectName("galaxy_cards_splitter")
        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(0, 0, 0, 0)
        left_layout.addWidget(QLabel("Recent galaxies"))
        self.search = QLineEdit()
        self.search.setObjectName("galaxy_cards_search")
        self.search.setClearButtonEnabled(True)
        self.search.setPlaceholderText("Search galaxies or folder paths…")
        left_layout.addWidget(self.search)
        self.projects = QListWidget()
        self.projects.setObjectName("galaxy_cards_project_list")
        self.projects.setWordWrap(True)
        self.projects.setSpacing(4)
        left_layout.addWidget(self.projects, 1)
        left.setMinimumWidth(280)
        self.splitter.addWidget(left)

        self.detail = GalaxyDetailCard()
        self.splitter.addWidget(self.detail)
        self.splitter.setChildrenCollapsible(False)
        self.splitter.setStretchFactor(0, 48)
        self.splitter.setStretchFactor(1, 52)
        self.splitter.setSizes(list(INITIAL_SPLITTER_SIZES))
        outer.addWidget(self.splitter, 1)

        self.search.textChanged.connect(self._filter_rows)
        self.projects.itemSelectionChanged.connect(self._selection_changed)
        self.projects.itemDoubleClicked.connect(lambda _: self._activate())
        self.detail.openRequested.connect(self.openRequested.emit)

    def set_snapshots(self, snapshots: list[GalaxyCardSnapshot]) -> None:
        previous = self.selected_id
        self._snapshots = {s.project_id: s for s in snapshots}
        self.projects.clear()
        for snapshot in snapshots:
            state = "On disk" if snapshot.available else "Missing"
            item = QListWidgetItem(
                f"{snapshot.display_name} · {state}\n{snapshot.coordinates}"
            )
            item.setData(Qt.UserRole, snapshot.project_id)
            item.setToolTip(str(snapshot.root))
            self.projects.addItem(item)
        self._filter_rows()
        if previous in self._snapshots:
            for i in range(self.projects.count()):
                item = self.projects.item(i)
                if item.data(Qt.UserRole) == previous and not item.isHidden():
                    self.projects.setCurrentRow(i)
                    break

    @property
    def selected_id(self) -> str | None:
        item = self.projects.currentItem()
        if item is None or item.isHidden():
            return None
        return item.data(Qt.UserRole)

    def _filter_rows(self, *_args) -> None:
        query = self.search.text()
        first_visible = -1
        for i in range(self.projects.count()):
            item = self.projects.item(i)
            snapshot = self._snapshots[item.data(Qt.UserRole)]
            matches = snapshot.matches(query)
            item.setHidden(not matches)
            if matches and first_visible == -1:
                first_visible = i
        if self.selected_id is None:
            self.projects.setCurrentRow(first_visible)
        self._selection_changed()

    def _selection_changed(self) -> None:
        selected = self.selected_id
        self.detail.set_snapshot(
            self._snapshots.get(selected) if selected is not None else None
        )
        if selected is not None:
            self.selectionChanged.emit(selected)

    def _activate(self) -> None:
        if self.selected_id is not None:
            self.openRequested.emit(self.selected_id)


def main() -> None:
    """Opt-in preview. Does not change or replace the real ProjectHub."""
    app = QApplication.instance() or QApplication(sys.argv)
    window = QMainWindow()
    window.setWindowTitle("DESIDERA · Compact Galaxy Passport (preview)")
    window.resize(1100, 720)
    host = QWidget()
    layout = QVBoxLayout(host)
    preview = GalaxyCardsPreview()
    preview.set_snapshots(snapshots_from_records(recent_projects()))
    preview.openRequested.connect(lambda _id: QMessageBox.information(
        window, "Preview only",
        "To open or relink a project, use the existing DESIDERA launcher."
    ))
    layout.addWidget(preview)
    hint = QLabel(
        "Read-only UI preview · no downloads, FITS changes or scientific QA."
    )
    hint.setWordWrap(True)
    layout.addWidget(hint)
    window.setCentralWidget(host)
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
