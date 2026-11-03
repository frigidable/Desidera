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
from .catalogue_links import confirmed_name_link
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
        container.setContentsMargins(14, 14, 14, 12)
        container.setSpacing(8)

        scroll = QScrollArea(self)
        scroll.setObjectName("galaxy_card_scroll")
        scroll.setFrameShape(QFrame.NoFrame)
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        inner = QWidget()
        body = QVBoxLayout(inner)
        body.setContentsMargins(0, 0, 0, 0)
        body.setSpacing(12)

        # Hero: concise identity flows alongside the square image. Subsequent
        # sections span the full card width directly below the image, instead
        # of leaving a permanently empty narrow column on the left.
        hero = QWidget()
        hero.setObjectName("galaxy_card_hero")
        hero_row = QHBoxLayout(hero)
        hero_row.setContentsMargins(0, 0, 0, 0)
        hero_row.setSpacing(14)
        self.preview_holder = QWidget()
        self.preview_holder.setObjectName("galaxy_card_preview_column")
        image_layout = QVBoxLayout(self.preview_holder)
        image_layout.setContentsMargins(0, 0, 0, 0)
        self.preview = QLabel("No preview")
        self.preview.setObjectName("galaxy_card_preview")
        self.preview.setAlignment(Qt.AlignCenter)
        self.preview.setTextFormat(Qt.PlainText)
        self.preview.setWordWrap(True)
        self.preview.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        image_layout.addWidget(self.preview)
        hero_row.addWidget(self.preview_holder, 0, Qt.AlignTop)

        # Kept as a nonvisual compatibility/provenance hook. A thumbnail has
        # no caption; detailed source information belongs in its tooltip.
        self.preview_note = _metadata("")
        self.preview_note.setObjectName("galaxy_card_preview_note")
        self.preview_note.hide()

        identity = QWidget()
        identity.setObjectName("galaxy_card_metadata_column")
        identity_col = QVBoxLayout(identity)
        identity_col.setContentsMargins(0, 0, 0, 0)
        identity_col.setSpacing(7)
        kicker = QLabel("SELECTED GALAXY")
        kicker.setObjectName("galaxy_card_section")
        identity_col.addWidget(kicker)
        self.name = QLabel("Choose a galaxy")
        self.name.setObjectName("galaxy_card_name")
        self.name.setWordWrap(True)
        identity_col.addWidget(self.name)
        self.coordinates = _metadata("Select a project from the list.")
        self.coordinates.setObjectName("galaxy_card_coordinates")
        self.coordinates.setTextInteractionFlags(Qt.TextSelectableByMouse)
        identity_col.addWidget(self.coordinates)
        self.catalogue = _metadata("No confirmed catalogue identity")
        self.catalogue.setObjectName("galaxy_card_catalogue")
        self.catalogue.setTextFormat(Qt.RichText)
        self.catalogue.setTextInteractionFlags(
            Qt.TextBrowserInteraction | Qt.TextSelectableByMouse)
        self.catalogue.setOpenExternalLinks(True)
        self.catalogue.setToolTip(
            "Only user-confirmed catalogue identities. External evidence "
            "does not change project coordinates or adopt measurements.")
        identity_col.addWidget(self.catalogue)
        identity_col.addStretch(1)
        hero_row.addWidget(identity, 1)
        body.addWidget(hero)

        # Explain what catalogue confirmation unlocks before any link exists.
        # This is informational; it never sends network requests on display.
        self.catalogue_help = _metadata("")
        self.catalogue_help.setObjectName("galaxy_card_catalogue_help")
        body.addWidget(self.catalogue_help)
        self.catalogue_help.hide()

        self.literature_heading = QLabel("LITERATURE · EXTERNAL VALUES")
        self.literature_heading.setObjectName("galaxy_card_section")
        body.addWidget(self.literature_heading)
        self.literature = _metadata("")
        self.literature.setObjectName("galaxy_card_literature")
        self.literature.setTextInteractionFlags(Qt.TextSelectableByMouse)
        body.addWidget(self.literature)
        self.cross_catalogue = _metadata("")
        self.cross_catalogue.setObjectName("galaxy_card_catalogue_offset")
        body.addWidget(self.cross_catalogue)

        observations = QLabel("OBSERVATIONS")
        observations.setObjectName("galaxy_card_section")
        body.addWidget(observations)
        self.dataset = _metadata("No dataset selected", selectable=True)
        self.dataset.setObjectName("galaxy_card_dataset")
        body.addWidget(self.dataset)
        self.fits = _metadata("FITS: not inspected")
        self.fits.setObjectName("galaxy_card_fits")
        body.addWidget(self.fits)
        self.observation_explanation = _metadata("")
        self.observation_explanation.setObjectName(
            "galaxy_card_observation_explanation")
        body.addWidget(self.observation_explanation)

        next_title = QLabel("NEXT STEP")
        next_title.setObjectName("galaxy_card_section")
        body.addWidget(next_title)
        self.next_step_title = next_title
        self.next_step = _metadata("")
        self.next_step.setObjectName("galaxy_card_next_step")
        body.addWidget(self.next_step)

        self.details_toggle = QPushButton("Catalogue & project details  ▾")
        self.details_toggle.setObjectName("galaxy_card_details_toggle")
        self.details_toggle.setCheckable(True)
        self.details_toggle.setChecked(False)
        body.addWidget(self.details_toggle)

        self.details_panel = QWidget()
        self.details_panel.setObjectName("galaxy_card_details_panel")
        extended = QVBoxLayout(self.details_panel)
        extended.setContentsMargins(6, 2, 6, 2)
        extended.setSpacing(9)

        provenance_title = QLabel("CATALOGUE SOURCES")
        provenance_title.setObjectName("galaxy_card_section")
        extended.addWidget(provenance_title)
        self.catalogue_provenance = _metadata("", selectable=True)
        self.catalogue_provenance.setObjectName(
            "galaxy_card_catalogue_provenance")
        self.catalogue_provenance.setTextFormat(Qt.PlainText)
        extended.addWidget(self.catalogue_provenance)

        project_title = QLabel("PROJECT INFORMATION")
        project_title.setObjectName("galaxy_card_section")
        extended.addWidget(project_title)
        self.aliases = _metadata("")
        self.aliases.setObjectName("galaxy_card_aliases")
        extended.addWidget(self.aliases)
        self.last_opened = _metadata("")
        self.last_opened.setObjectName("galaxy_card_last_opened")
        extended.addWidget(self.last_opened)
        self.path = _metadata("", selectable=True)
        self.path.setObjectName("galaxy_card_path")
        extended.addWidget(self.path)
        self.note = _metadata("")
        self.note.setObjectName("galaxy_card_note")
        extended.addWidget(self.note)
        body.addWidget(self.details_panel)
        self.details_panel.hide()
        self.details_toggle.toggled.connect(self._details_toggled)
        body.addStretch(1)

        scroll.setWidget(inner)
        container.addWidget(scroll, 1)
        self.open_button = QPushButton("Open project")
        self.open_button.setObjectName("galaxy_card_open")
        self.open_button.setEnabled(False)
        self.open_button.setToolTip(
            "Opens the project offline. Original FITS and previous scientific "
            "QA are never automatically adopted.")
        self.open_button.clicked.connect(self._request_open)
        container.addWidget(self.open_button)

        # Kept for older integrations but deliberately absent from the
        # displayed card: the repetitive caveat is now a button tooltip.
        self.safety = _metadata("")
        self.safety.setObjectName("galaxy_card_safety")
        self.safety.hide()
        self._set_preview_side(184)
        self.set_snapshot(None)

    def _details_toggled(self, checked: bool):
        self.details_panel.setVisible(checked)
        self.details_toggle.setText(
            "Hide catalogue & project details  ▴" if checked
            else "Catalogue & project details  ▾")

    def _set_preview_side(self, side: int):
        """Square preview, responsive down to the 300 px narrow inspector."""
        if self.preview.width() == side and self.preview.height() == side:
            return
        self.preview.setFixedSize(side, side)
        self.preview_holder.setFixedWidth(side)
        if self._pixmap is not None:
            self._update_preview()

    def set_snapshot(self, snapshot: GalaxyCardSnapshot | None) -> None:
        previous_id = self.snapshot.project_id if self.snapshot else None
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
            self.catalogue_help.hide()
            self.literature_heading.hide()
            self.literature.hide()
            self.cross_catalogue.hide()
            self.observation_explanation.hide()
            self.next_step_title.hide()
            self.next_step.hide()
            self.last_opened.hide()
            self.path.clear()
            self.note.hide()
            return

        self.name.setText(snapshot.display_name)
        self.coordinates.setText(snapshot.coordinates)
        self.next_step.show()
        if (previous_id != snapshot.project_id
                and snapshot.knowledge_status != "EVIDENCE_ERROR"):
            self.details_toggle.setChecked(False)
        self.dataset.setText(f"Dataset: {snapshot.dataset}")
        self.literature_heading.hide()
        self.literature.hide()
        self.cross_catalogue.hide()
        self.catalogue_help.hide()
        self.next_step_title.show()
        self.observation_explanation.hide()
        catalogues = snapshot.catalogue_identities or (
            (snapshot.catalogue_identity,) if snapshot.catalogue_identity else ())
        self.catalogue_provenance.hide()
        if snapshot.knowledge_status == "EVIDENCE_ERROR":
            self.catalogue.setText("Catalogue evidence requires inspection")
            self.catalogue_provenance.setText(
                "Stored catalogue evidence is invalid or unreadable: "
                + snapshot.knowledge_error)
            self.catalogue_provenance.show()
            self.details_toggle.setChecked(True)
        elif catalogues:
            # CatalogueIdentity instances are produced exclusively by the
            # receipt-validated knowledge read model. Neither unidentified
            # cone candidates nor corrupt receipts receive browser links.
            self.catalogue.setText("<br>".join(
                f"{entry.provider}: "
                + confirmed_name_link(
                    entry.provider, entry.main_id,
                    endpoint_host=entry.endpoint_host)
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
                "Catalogue candidates found · no identity confirmed")
        else:
            self.catalogue.setText("Not yet linked to SIMBAD or NED")

        if snapshot.knowledge_status != "EVIDENCE_ERROR":
            providers = {entry.provider for entry in catalogues}
            if not providers:
                self.catalogue_help.setText(
                    "Why identify this galaxy? Confirm a SIMBAD/NED match "
                    "to connect published names, source references and "
                    "comparable catalogue values to this project. "
                    "Use Catalogues ▾ to review candidates; no scientific "
                    "results or ICRS coordinates are adopted automatically.")
                self.catalogue_help.show()
            elif "SIMBAD" not in providers:
                self.catalogue_help.setText(
                    "SIMBAD is not linked yet. Confirm a match to connect "
                    "literature names and bibliography to this project.")
                self.catalogue_help.show()
            elif "NED" not in providers:
                self.catalogue_help.setText(
                    "NED is not linked yet. Confirm a match to review "
                    "published extragalactic properties with source references.")
                self.catalogue_help.show()
            ned = next((entry for entry in catalogues
                        if entry.provider == "NED"), None)
            if ned is not None and (ned.object_type or ned.redshift is not None):
                literature = []
                if ned.object_type:
                    meaning = " (galaxy)" if ned.object_type == "G" else ""
                    literature.append(
                        f"NED classification: {ned.object_type}{meaning}")
                if ned.redshift is not None:
                    literature.append(
                        f"NED heliocentric z = {ned.redshift:.7g} "
                        f"· refcode {ned.redshift_refcode or 'unspecified'}")
                literature.append(
                    "External catalogue statements · not adopted as "
                    "DESIDERA measurements or distance estimates.")
                self.literature.setText("\n".join(literature))
                self.literature_heading.show()
                self.literature.show()
            if snapshot.provider_offset_arcsec is not None:
                self.cross_catalogue.setText(
                    f"SIMBAD–NED positions differ by "
                    f"{snapshot.provider_offset_arcsec:.2f}″. "
                    "Compare sources in details; this does not establish "
                    "a wrong identification.")
                self.cross_catalogue.show()

        project_state = (
            "On disk" if snapshot.available else "Missing · relink required"
        )
        self.fits.setText(f"{project_state} · {snapshot.fits_status}")
        self.next_step.setText("Open the project to review its current data.")
        self.open_button.setText(
            "Relink project…" if not snapshot.available else "Open project")
        if snapshot.readiness is not None:
            self.fits.setText(snapshot.readiness.observation_label)
            self.observation_explanation.setText(
                snapshot.readiness.explanation)
            self.observation_explanation.show()
            self.next_step.setText(snapshot.readiness.next_action_label)
            self.open_button.setText(snapshot.readiness.next_action_label)
        elif not snapshot.available:
            self.next_step.setText("Relink the project folder to continue.")
        self.aliases.setVisible(bool(snapshot.aliases))
        self.aliases.setText("Also known as: " + ", ".join(snapshot.aliases))
        self.last_opened.setVisible(bool(snapshot.last_opened))
        self.last_opened.setText("Last opened: " + snapshot.last_opened)
        self.path.setText("Project folder: " + str(snapshot.root))
        self.note.setVisible(bool(snapshot.note))
        self.note.setText(snapshot.note)
        self.preview.setPixmap(QPixmap())
        self.preview.setText("No image yet")
        self.preview_note.hide()
        if snapshot.preview_source == "REMOTE_LEGACY_RGB":
            tooltip = ("Legacy Surveys Viewer RGB cutout; display only. "
                       "Source and SHA-256 in the local preview receipt.")
        elif snapshot.preview_status == "READY_RECENT_STAT_NOT_FITS_REHASHED":
            tooltip = ("Local original-FITS image; display only. "
                       "Source FITS integrity not rehashed at startup.")
        else:
            tooltip = "Preview only. No photometric measurement or QA adoption."
        if snapshot.thumbnail is not None:
            image = QPixmap(str(snapshot.thumbnail))
            if (image.isNull()
                    and snapshot.preview_source == "REMOTE_LEGACY_RGB"
                    and snapshot.local_thumbnail is not None):
                image = QPixmap(str(snapshot.local_thumbnail))
                tooltip = "Remote RGB unreadable; local FITS fallback (display only)"
            if not image.isNull():
                self._pixmap = image
                self._update_preview()
        if snapshot.thumbnail is None:
            tooltip += " Link local FITS or select a dataset to obtain a preview."
        self.preview.setToolTip(tooltip)

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
