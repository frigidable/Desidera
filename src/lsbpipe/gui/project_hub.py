"""Offline-first startup Project Hub and object creation wizard."""
from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from PySide6.QtCore import Qt, Signal, QSize, QTimer
from PySide6.QtWidgets import (
    QDialog, QDialogButtonBox, QFileDialog, QFormLayout, QFrame,
    QHBoxLayout, QLabel, QLineEdit, QListWidget, QListWidgetItem,
    QMainWindow, QMessageBox, QPushButton, QStackedWidget, QVBoxLayout,
    QInputDialog, QSplitter,
    QWidget, QDoubleSpinBox, QApplication,
)

from lsbpipe.surveys.legacy import SkyField
from lsbpipe.surveys.object_projects import (
    create_object_project, import_object_project, read_object_project,
    recent_projects, relink_object_project, touch_object_project,
    forget_recent_project,
)
from lsbpipe.surveys.local_project import open_local_project
from lsbpipe.gui.target_input import parse_coordinate_pair
from .survey_explorer import SurveyExplorer
from .galaxy_cards import GalaxyDetailCard
from .galaxy_card_model import snapshot_from_record
from .card_previews import build_local_preview, has_local_image_source
from .galaxy_knowledge_dialog import GalaxyKnowledgeDialog
from .ned_knowledge_dialog import NedKnowledgeDialog
from .analysis_configurations_dialog import AnalysisConfigurationsDialog
from lsbpipe.surveys.galaxy_knowledge import latest_simbad_confirmation
from lsbpipe.surveys.knowledge_summary import load_galaxy_knowledge


_PREVIEW_POOL = ThreadPoolExecutor(
    max_workers=1, thread_name_prefix="desidera-card-preview")


class CreateProjectDialog(QDialog):
    """One-page object creation. No network or FITS access."""

    def __init__(self, parent=None, *, default_root: Path | None = None):
        super().__init__(parent)
        self.setWindowTitle("New galaxy project · DESIDERA")
        self.setMinimumWidth(540)
        self.resize(620, 510)
        self.project = None
        self.parent_root = default_root or Path.home() / "DESIDERAProjects"
        self._directory_user_edited = False

        self.name = QLineEdit()
        self.name.setObjectName("project_new_name")
        self.name.setPlaceholderText("e.g. UGC 2162")
        self.ra = QLineEdit()
        self.ra.setObjectName("project_new_ra")
        self.ra.setPlaceholderText("Optional · degrees or hh:mm:ss")
        self.dec = QLineEdit()
        self.dec.setObjectName("project_new_dec")
        self.dec.setPlaceholderText("Optional · degrees or ±dd:mm:ss")
        self.width_arcmin = QDoubleSpinBox()
        self.width_arcmin.setObjectName("project_new_width")
        self.width_arcmin.setRange(0.1, 120.)
        self.width_arcmin.setValue(5.)
        self.width_arcmin.setSuffix(" arcmin")
        self.directory = QLineEdit()
        self.directory.setObjectName("project_new_root")
        self.directory.setPlaceholderText("Local folder for this galaxy")
        self.description = QLabel()
        self.description.setObjectName("project_new_summary")
        self.description.setWordWrap(True)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(24, 20, 24, 20)
        layout.setSpacing(12)
        title = QLabel("Create galaxy project")
        title.setStyleSheet("font-size: 20px; font-weight: 600")
        layout.addWidget(title)
        subtitle = QLabel(
            "Each galaxy has its own local workspace. Original FITS files "
            "can be connected later without downloading them again.")
        subtitle.setWordWrap(True)
        layout.addWidget(subtitle)

        identity = QFormLayout()
        identity.setVerticalSpacing(10)
        identity.addRow("Object / project name", self.name)
        identity.addRow("RA (ICRS)", self.ra)
        identity.addRow("Dec (ICRS)", self.dec)
        identity.addRow("Search field", self.width_arcmin)
        layout.addLayout(identity)
        note = QLabel(
            "Coordinates are optional for an unidentified object. "
            "Provide both RA and Dec to define its sky position.")
        note.setWordWrap(True)
        layout.addWidget(note)

        layout.addWidget(QLabel("Local project storage"))
        storage = QHBoxLayout()
        storage.addWidget(self.directory, 1)
        browse = QPushButton("Browse…")
        browse.setObjectName("project_new_browse")
        browse.clicked.connect(self._browse)
        storage.addWidget(browse)
        layout.addLayout(storage)
        self.description.setTextInteractionFlags(Qt.TextSelectableByMouse)
        layout.addWidget(self.description)

        layout.addStretch(1)
        buttons = QHBoxLayout()
        self.cancel = QPushButton("Cancel")
        self.cancel.clicked.connect(self.reject)
        self.create_btn = QPushButton("Create project")
        self.create_btn.setObjectName("project_wizard_create")
        self.create_btn.setDefault(True)
        self.create_btn.clicked.connect(self._create)
        buttons.addStretch(1)
        buttons.addWidget(self.cancel)
        buttons.addWidget(self.create_btn)
        layout.addLayout(buttons)

        self.name.textChanged.connect(self._suggest_directory)
        self.name.textChanged.connect(self._update_preview)
        self.ra.textChanged.connect(self._update_preview)
        self.dec.textChanged.connect(self._update_preview)
        self.width_arcmin.valueChanged.connect(self._update_preview)
        self.directory.textEdited.connect(self._note_directory_edit)
        self.directory.textChanged.connect(self._update_preview)
        self._update_preview()
        self.name.setFocus()

    def _suggest_directory(self, _=None):
        if self._directory_user_edited:
            return
        value = self.name.text().strip()
        safe = "".join(
            ch if ch.isalnum() or ch in ("_", "-", " ") else "_"
            for ch in value).strip().replace(" ", "_")
        self.directory.setText(str(self.parent_root / safe) if safe else "")

    def _note_directory_edit(self, *_):
        self._directory_user_edited = True

    def _browse(self):
        folder = QFileDialog.getExistingDirectory(
            self, "Choose parent folder for the new project",
            str(self.parent_root))
        if folder:
            self._directory_user_edited = True
            leaf = self.name.text().strip().replace(" ", "_") or "NewGalaxy"
            self.directory.setText(str(Path(folder) / leaf))

    def _field(self) -> SkyField | None:
        ra, dec = self.ra.text().strip(), self.dec.text().strip()
        if not ra and not dec:
            return None
        if not ra or not dec:
            raise ValueError("Provide both RA and Dec or leave both blank")
        from lsbpipe.gui.target_input import parse_ra_deg, parse_dec_deg
        return SkyField(parse_ra_deg(ra), parse_dec_deg(dec),
                        self.width_arcmin.value())

    def _update_preview(self, *_):
        name = self.name.text().strip()
        directory = self.directory.text().strip()
        try:
            field = self._field()
            position = ("Position: not specified (can be added later)"
                        if field is None else
                        f"ICRS: {field.ra_deg:.8f}°, {field.dec_deg:+.8f}°")
        except ValueError as exc:
            position = f"Check coordinates: {exc}"
        self.description.setText(
            f"{position}\n"
            f"Project directory: {directory or 'not selected'}\n"
            "Creates local metadata only · no network, download or FITS changes.")
        self.create_btn.setEnabled(bool(name and directory))

    def _create(self):
        try:
            name = self.name.text().strip()
            if not name:
                raise ValueError("Enter an object or project name")
            directory = self.directory.text().strip()
            if not directory:
                raise ValueError("Choose a project directory")
            folder = Path(directory).expanduser()
            if folder.exists() and any(folder.iterdir()):
                raise FileExistsError(
                    "Folder contains existing files. Use Import existing "
                    "folder in Project Hub to preserve those files.")
            self.project = create_object_project(
                name, folder, field=self._field())
        except (ValueError, OSError) as exc:
            self.description.setText(f"Cannot create project: {exc}")
            QMessageBox.warning(self, "Cannot create project", str(exc))
            return
        self.accept()


class ProjectHub(QMainWindow):
    """Launch surface. Explorer is instantiated only for a selected project."""
    projectOpened = Signal(str)

    def __init__(self):
        super().__init__()
        self.setObjectName("lsbpipe_project_hub")
        self.setWindowTitle("DESIDERA · Galaxy Projects")
        self.resize(1100, 720)
        self.workspace = None
        self._preview_future = None
        self._preview_project_id = None
        self._preview_errors = {}
        self._preview_timer = QTimer(self)
        self._preview_timer.setInterval(100)
        self._preview_timer.timeout.connect(self._finish_preview)

        content = QWidget()
        outer = QVBoxLayout(content)
        outer.setContentsMargins(24, 20, 24, 18)
        outer.setSpacing(14)

        heading = QLabel("DESIDERA")
        heading.setStyleSheet("font-size: 26px; font-weight: 700")
        subheading = QLabel(
            "Galaxy projects · local-first scientific workspace")
        outer.addWidget(heading)
        outer.addWidget(subheading)

        actions = QHBoxLayout()
        self.new_btn = QPushButton("New galaxy project")
        self.new_btn.setObjectName("hub_new_project")
        self.new_btn.clicked.connect(self._new_project)
        self.import_btn = QPushButton("Import existing folder…")
        self.import_btn.setObjectName("hub_import_project")
        self.import_btn.clicked.connect(self._import_dialog)
        actions.addWidget(self.new_btn)
        actions.addWidget(self.import_btn)
        actions.addStretch(1)
        outer.addLayout(actions)

        splitter = QSplitter(Qt.Horizontal)
        splitter.setObjectName("hub_project_splitter")
        self.splitter = splitter
        library = QWidget()
        left = QVBoxLayout(library)
        left.setContentsMargins(0, 0, 0, 0)
        left.setSpacing(10)
        left.addWidget(QLabel("Recent galaxies"))
        self.filter = QLineEdit()
        self.filter.setObjectName("hub_project_filter")
        self.filter.setClearButtonEnabled(True)
        self.filter.setPlaceholderText("Search galaxies or folder paths…")
        self.filter.textChanged.connect(self._filter_rows)
        left.addWidget(self.filter)
        self.projects = QListWidget()
        self.projects.setObjectName("hub_project_list")
        self.projects.setWordWrap(True)
        self.projects.setSpacing(5)
        self.projects.setStyleSheet(
            "QListWidget { border: 1px solid palette(mid); "
            "border-radius: 8px; padding: 6px; } "
            "QListWidget::item { padding: 6px; border-radius: 5px; } "
            "QListWidget::item:selected { background: palette(highlight); "
            "color: palette(highlighted-text); }")
        self.projects.itemDoubleClicked.connect(
            lambda _: self._open_selected())
        self.projects.itemSelectionChanged.connect(self._selection_changed)
        left.addWidget(self.projects, 1)
        self.empty_panel = QWidget()
        self.empty_panel.setObjectName("hub_empty_panel")
        empty_layout = QVBoxLayout(self.empty_panel)
        empty_layout.addStretch(1)
        self.empty_label = QLabel(
            "No projects yet.\nCreate a galaxy project or import an "
            "existing folder to begin.")
        self.empty_label.setObjectName("hub_empty_state")
        self.empty_label.setWordWrap(True)
        self.empty_label.setAlignment(Qt.AlignCenter)
        self.empty_label.setStyleSheet("font-size: 16px")
        empty_layout.addWidget(self.empty_label)
        empty_layout.addStretch(1)
        left.addWidget(self.empty_panel, 1)
        library.setMinimumWidth(310)
        splitter.addWidget(library)

        details = QFrame()
        details.setObjectName("hub_project_details")
        details.setFrameShape(QFrame.NoFrame)
        right = QVBoxLayout(details)
        right.setContentsMargins(0, 0, 0, 0)
        right.setSpacing(12)
        self.detail_heading = QLabel("SELECTED GALAXY")
        right.addWidget(self.detail_heading)
        self.detail_title = QLabel("Choose a project")
        self.detail_title.setObjectName("hub_detail_title")
        self.detail_title.setStyleSheet(
            "font-size: 19px; font-weight: 600")
        self.detail_title.setWordWrap(True)
        right.addWidget(self.detail_title)
        self.detail_status = QLabel()
        self.detail_status.setObjectName("hub_detail_status")
        self.detail_status.setWordWrap(True)
        right.addWidget(self.detail_status)
        self.detail_meta = QLabel()
        self.detail_meta.setObjectName("hub_detail_meta")
        self.detail_meta.setWordWrap(True)
        self.detail_meta.setTextInteractionFlags(Qt.TextSelectableByMouse)
        right.addWidget(self.detail_meta)
        self.detail_path = QLabel()
        self.detail_path.setObjectName("hub_detail_path")
        self.detail_path.setWordWrap(True)
        self.detail_path.setTextInteractionFlags(Qt.TextSelectableByMouse)
        right.addWidget(self.detail_path)
        self.detail_note = QLabel()
        self.detail_note.setObjectName("hub_detail_note")
        self.detail_note.setWordWrap(True)
        right.addWidget(self.detail_note)

        self.continue_btn = QPushButton("Open project")
        self.continue_btn.setObjectName("hub_continue")
        self.continue_btn.clicked.connect(self._open_selected)
        right.addWidget(self.continue_btn)
        self.preview_btn = QPushButton("Regenerate local FITS preview")
        self.preview_btn.setObjectName("hub_generate_card_preview")
        self.preview_btn.setToolTip(
            "Generate a cached monochrome IMAGE thumbnail from original "
            "local Legacy FITS; never alter originals or run photometry.")
        self.preview_btn.clicked.connect(self._regenerate_preview)
        right.addWidget(self.preview_btn)
        self.knowledge_btn = QPushButton("Galaxy Knowledge · SIMBAD…")
        self.knowledge_btn.setObjectName("hub_knowledge_simbad")
        self.knowledge_btn.setToolTip(
            "Open cached identity information offline; remote SIMBAD query "
            "requires explicit Search and cannot change target coordinates.")
        self.knowledge_btn.clicked.connect(self._open_knowledge)
        right.addWidget(self.knowledge_btn)
        self.ned_btn = QPushButton("Galaxy Knowledge · NED…")
        self.ned_btn.setObjectName("hub_knowledge_ned")
        self.ned_btn.setToolTip(
            "Open saved NED evidence offline; remote ObjectLookup is "
            "opt-in and never modifies project coordinates or science.")
        self.ned_btn.clicked.connect(self._open_ned)
        right.addWidget(self.ned_btn)
        self.config_btn = QPushButton("Analysis configurations…")
        self.config_btn.setObjectName("hub_analysis_configurations")
        self.config_btn.setToolTip(
            "Review immutable project setup snapshots. Freezing requires "
            "explicit action and never runs or adopts scientific results.")
        self.config_btn.clicked.connect(self._open_configurations)
        right.addWidget(self.config_btn)
        self.remove_btn = QPushButton("Remove from recent…")
        self.remove_btn.setObjectName("hub_remove_recent")
        self.remove_btn.setToolTip(
            "Forget only the Project Hub entry. Original FITS, directories, "
            "project metadata and science reports are never deleted.")
        self.remove_btn.clicked.connect(self._remove_selected)
        right.addWidget(self.remove_btn)
        # Integrate the reviewed compact passport without removing the
        # independently reviewed SIMBAD/configuration actions.
        self.galaxy_card = GalaxyDetailCard()
        right.insertWidget(0, self.galaxy_card, 1)
        self.galaxy_card.openRequested.connect(self._open_card_project)
        for old_widget in (
            self.detail_heading, self.detail_title, self.detail_status,
            self.detail_meta, self.detail_path, self.detail_note,
            self.continue_btn,
        ):
            old_widget.hide()
        details.setMinimumWidth(300)
        splitter.addWidget(details)
        splitter.setStretchFactor(0, 48)
        splitter.setStretchFactor(1, 52)
        splitter.setSizes([480, 540])
        outer.addWidget(splitter, 1)

        self.hint = QLabel(
            "Projects open offline. Local FITS are read-only; earlier "
            "diagnostics are not automatically adopted as scientific results.")
        self.hint.setWordWrap(True)
        outer.addWidget(self.hint)
        self.setCentralWidget(content)
        self.refresh()

    def refresh(self):
        prior = self.projects.currentItem()
        prior_id = (prior.data(Qt.UserRole).get("project_id")
                    if prior is not None else None)
        self.projects.clear()
        try:
            rows = recent_projects()
        except (ValueError, OSError) as exc:
            self.hint.setText(
                f"Project library index could not be read: {exc}. "
                "Existing project folders are untouched; import them again.")
            rows = []
        for record in rows:
            exists = bool(record["available"])
            name = record.get("display_name") or "Unnamed"
            status = "On disk" if exists else "Missing · relink needed"
            detail = ""
            search_identifiers = []
            if exists:
                try:
                    if record.get("legacy"):
                        prior_project = open_local_project(record["root"])
                        detail = (
                            f"Older project · "
                            f"{prior_project.field.ra_deg:.6f}°, "
                            f"{prior_project.field.dec_deg:+.6f}°")
                    else:
                        entry = read_object_project(record["root"])
                        search_identifiers.extend(entry.aliases)
                        try:
                            knowledge = load_galaxy_knowledge(
                                entry.root, project=entry)
                            # Only explicitly confirmed names are searchable.
                            # Unverified nearby stars must never be indexed as
                            # aliases of the galaxy.
                            search_identifiers.extend(knowledge.confirmed_ids)
                        except (OSError, ValueError, KeyError, TypeError):
                            pass  # The selected card shows an evidence error.
                        detail = (
                            f"{entry.field.ra_deg:.6f}°, "
                            f"{entry.field.dec_deg:+.6f}°"
                            if entry.field is not None else
                            "Coordinates not specified")
                except (OSError, ValueError, KeyError, TypeError):
                    detail = "Project metadata needs attention"
            item = QListWidgetItem(
                f"{name}  ·  {status}\n{detail}")
            item.setSizeHint(QSize(100, 64))
            item.setData(Qt.UserRole, record)
            item.setData(Qt.UserRole + 1, tuple(search_identifiers))
            item.setToolTip(record["root"])
            self.projects.addItem(item)
        selection = next(
            (i for i in range(self.projects.count())
             if self.projects.item(i).data(Qt.UserRole).get("project_id")
             == prior_id), 0)
        if self.projects.count():
            self.projects.setCurrentRow(selection)
        self._filter_rows()

    def _filter_rows(self):
        query = self.filter.text().strip().casefold()
        visible = []
        for index in range(self.projects.count()):
            item = self.projects.item(index)
            record = item.data(Qt.UserRole)
            matches = (query in item.text().casefold() or
                       query in str(record.get("root", "")).casefold() or
                       query in str(record.get("project_id", "")).casefold() or
                       any(query in value.casefold() for value in
                           (item.data(Qt.UserRole + 1) or ())))
            item.setHidden(not matches)
            if matches:
                visible.append(index)
        current = self.projects.currentItem()
        if current is None or current.isHidden():
            self.projects.setCurrentRow(visible[0] if visible else -1)
        self.projects.setVisible(bool(visible))
        self.empty_panel.setVisible(not visible)
        self.empty_label.setText(
            "No matching galaxies.\nClear the search to see all projects."
            if query and not visible else
            "No projects yet.\nCreate a galaxy project or import an "
            "existing folder to begin.")
        self._selection_changed()

    def _selection_changed(self):
        item = self.projects.currentItem()
        record = (item.data(Qt.UserRole)
                  if item is not None and not item.isHidden() else None)
        selected = isinstance(record, dict)
        card = snapshot_from_record(record) if selected else None
        self.galaxy_card.set_snapshot(card)
        self.galaxy_card.open_button.setEnabled(
            selected and self.workspace is None)
        self.continue_btn.setEnabled(selected)
        can_preview = (
            selected and self.workspace is None
            and bool(record.get("available")) and not bool(record.get("legacy")))
        self.preview_btn.setEnabled(
            can_preview and self._preview_future is None)
        self.remove_btn.setEnabled(selected and self.workspace is None)
        if can_preview and card is not None:
            self._auto_preview_selected(record, card)
        self.knowledge_btn.setEnabled(
            selected and self.workspace is None and
            bool(record.get("available")) and not bool(record.get("legacy")))
        self.config_btn.setEnabled(
            selected and self.workspace is None and
            bool(record.get("available")) and not bool(record.get("legacy")))
        self.ned_btn.setEnabled(
            selected and self.workspace is None and
            bool(record.get("available")) and not bool(record.get("legacy")))
        self.continue_btn.setText(
            "Relink project…" if selected and not record["available"]
            else "Open project")

        if not selected:
            self.detail_title.setText("Choose a galaxy")
            self.detail_status.setText(
                "Select an existing project or create a new one.")
            self.detail_meta.clear()
            self.detail_path.clear()
            self.detail_note.clear()
            return
        self.detail_title.setText(record.get("display_name") or "Unnamed")
        self.detail_path.setText(
            "Project folder:\n" + str(record.get("root", "")))
        if not record["available"]:
            self.detail_status.setText("Project folder is missing.")
            self.detail_meta.setText("The saved project identity is preserved.")
            self.detail_note.setText(
                "Relink the moved project folder. No local data are deleted.")
            return
        if record.get("legacy"):
            self.detail_status.setText("Older local project")
            self.detail_meta.setText("Existing survey data can be imported.")
            self.detail_note.setText(
                "Opening creates a project identity without changing "
                "original survey FITS.")
            return
        try:
            project = read_object_project(record["root"])
        except (OSError, ValueError, KeyError, TypeError) as exc:
            self.detail_status.setText("Project metadata needs attention")
            self.detail_meta.setText(str(exc))
            self.detail_note.setText(
                "Opening will report the error; original files are untouched.")
            return

        position = (
            f"ICRS: {project.field.ra_deg:.8f}°, "
            f"{project.field.dec_deg:+.8f}°\n"
            f"Search field: {project.field.width_arcmin:g} arcmin"
            if project.field is not None else
            "ICRS: not specified yet")
        if project.selected_dataset:
            release, region, brick, band = project.selected_dataset
            dataset = f"{release.upper()} / {region} / {brick} / {band}"
        else:
            dataset = "None selected"
        if project.data_root is None:
            state = "No saved local FITS link"
            location = "FITS location: not linked"
        elif not project.data_root.is_dir():
            state = "Saved FITS location is unavailable"
            location = f"FITS location: {project.data_root}"
        else:
            state = "Local FITS location linked · integrity not checked"
            location = f"FITS location: {project.data_root}"
        self.detail_status.setText(state)
        self.detail_meta.setText(
            f"{position}\nSelected dataset: {dataset}\n{location}")
        note = ("Opening does not download files or approve previous QA. "
                "Inspect original FITS and run diagnostics when needed.")
        try:
            confirmed = latest_simbad_confirmation(project.root)
            if confirmed:
                note += ("\nSIMBAD (user-confirmed catalogue association): "
                         + confirmed["candidate"]["main_id"]
                         + " · ICRS position unchanged")
        except (OSError, ValueError, KeyError, TypeError) as exc:
            note += f"\nKnowledge evidence unavailable or corrupt: {exc}"
        self.detail_note.setText(note)

    def _selected_preview_record(self):
        item = self.projects.currentItem()
        if item is None or item.isHidden() or self.workspace is not None:
            return None
        record = item.data(Qt.UserRole)
        if (not isinstance(record, dict) or not record.get("available")
                or record.get("legacy")):
            return None
        return record

    def _auto_preview_selected(self, record, card):
        if card.thumbnail is not None or self._preview_future is not None:
            return
        project_id = record.get("project_id")
        if project_id in self._preview_errors:
            self.galaxy_card.preview.setText(
                "Local preview failed: " + self._preview_errors[project_id]
                + ". Use Regenerate local FITS preview to retry.")
            return
        try:
            project = read_object_project(record["root"])
            if not has_local_image_source(project):
                return
        except (OSError, ValueError, KeyError, TypeError, AttributeError):
            return
        self._start_preview(record)

    def _start_preview(self, record):
        if self._preview_future is not None:
            return
        self._preview_project_id = record["project_id"]
        self._preview_future = _PREVIEW_POOL.submit(
            build_local_preview, record["root"])
        self.preview_btn.setEnabled(False)
        self.galaxy_card.preview.setText(
            "Building local FITS thumbnail (display only)…")
        self._preview_timer.start()

    def _finish_preview(self):
        if self._preview_future is None or not self._preview_future.done():
            return
        task, completed_id = self._preview_future, self._preview_project_id
        self._preview_future = None
        self._preview_project_id = None
        self._preview_timer.stop()
        try:
            task.result()
            self._preview_errors.pop(completed_id, None)
        except Exception as exc:
            self._preview_errors[completed_id] = (
                f"{type(exc).__name__}: {exc}")
        if self._selected_preview_record() is not None:
            self._selection_changed()

    def _regenerate_preview(self):
        record = self._selected_preview_record()
        if record is None or self._preview_future is not None:
            return
        self._preview_errors.pop(record["project_id"], None)
        self._start_preview(record)

    def _open_card_project(self, project_id: str):
        """Activate the currently selected original project by immutable ID."""
        item = self.projects.currentItem()
        record = item.data(Qt.UserRole) if item is not None else None
        if (self.workspace is None and item is not None and not item.isHidden()
                and isinstance(record, dict)
                and record.get("project_id") == project_id):
            self._open_selected()

    def _remove_selected(self):
        """Forget a recent entry, never delete a research project on disk."""
        if self.workspace is not None:
            return
        item = self.projects.currentItem()
        if item is None or item.isHidden():
            return
        record = item.data(Qt.UserRole)
        if not isinstance(record, dict):
            return
        name = record.get("display_name") or "Unnamed"
        root = record["root"]
        response = QMessageBox.question(
            self, "Remove from recent projects?",
            f"Remove this entry from the DESIDERA recent-project list?\n\n"
            f"Project: {name}\n"
            f"Location: {root}\n\n"
            "Only the list entry will be removed. The directory, original "
            "FITS files, project metadata and science reports will NOT "
            "be deleted. You can import the project again at any time.",
            QMessageBox.Yes | QMessageBox.Cancel, QMessageBox.Cancel)
        if response != QMessageBox.Yes:
            return
        try:
            forget_recent_project(record["project_id"], root)
        except (ValueError, OSError, KeyError, TypeError) as exc:
            QMessageBox.warning(
                self, "Could not remove recent project", str(exc))
            return
        self.refresh()
        self.hint.setText(
            f"Removed {name} from Recent projects. No files were deleted.")

    def _new_project(self):
        dialog = CreateProjectDialog(self)
        if dialog.exec() == QDialog.Accepted and dialog.project is not None:
            # A new object starts offline. The hosted sky preview is an
            # explicit researcher action inside the workspace.
            self._enter_project(dialog.project, show_remote=False)

    def _import_dialog(self):
        directory = QFileDialog.getExistingDirectory(
            self, "Import existing DESIDERA or Legacy FITS folder",
            str(Path.home()))
        if directory:
            try:
                root = Path(directory).expanduser().resolve()
                manifest = root / "lsbpipe-object.json"
                if manifest.is_file():
                    project = read_object_project(root)
                    description = (
                        f"Existing project: {project.display_name}\n"
                        f"Identifier: {project.project_id}\n"
                        f"Directory: {root}\n\n"
                        "Open this project without downloading FITS?")
                else:
                    preview = open_local_project(root)
                    complete = sum(preview.complete.values())
                    description = (
                        f"Existing local data: {preview.root}\n"
                        f"Target RA {preview.field.ra_deg:.6f}°, "
                        f"Dec {preview.field.dec_deg:+.6f}°\n"
                        f"Bricks: {len(preview.hits)}; "
                        f"complete band quartets: "
                        f"{complete}/{len(preview.complete)}\n\n"
                        "Create an object-project identity beside the "
                        "existing FITS? No originals will be changed.")
                if QMessageBox.question(
                        self, "Review local project import", description,
                        QMessageBox.Yes | QMessageBox.No,
                        QMessageBox.No) == QMessageBox.Yes:
                    self.import_folder(root)
            except (OSError, ValueError, KeyError) as exc:
                # A metadata-free legacy folder cannot imply an ICRS center.
                # Ask explicitly instead of inventing a coordinate from brick.
                if "No saved target coordinates" not in str(exc):
                    QMessageBox.warning(
                        self, "Cannot inspect local folder", str(exc))
                    return
                pair, ok = QInputDialog.getText(
                    self, "Missing target position",
                    "Enter ICRS RA and Dec (degrees or sexagesimal):")
                if not ok:
                    return
                try:
                    ra, dec = parse_coordinate_pair(pair)
                    fallback = SkyField(ra, dec, 5.)
                except ValueError as invalid:
                    QMessageBox.warning(self, "Invalid ICRS position", str(invalid))
                    return
                self.import_folder(root, fallback_field=fallback)

    def import_folder(self, directory: str | Path, *,
                      fallback_field: SkyField | None = None) -> bool:
        try:
            project = import_object_project(
                directory, fallback_field=fallback_field)
        except (OSError, ValueError, KeyError) as exc:
            self.hint.setText(f"Import failed: {exc}")
            QMessageBox.warning(self, "Import was not completed", str(exc))
            return False
        self._enter_project(project)
        return True

    def _open_configurations(self):
        item = self.projects.currentItem()
        if item is None or item.isHidden() or self.workspace is not None:
            return
        record = item.data(Qt.UserRole)
        if (not isinstance(record, dict) or not record.get("available")
                or record.get("legacy")):
            return
        try:
            dialog = AnalysisConfigurationsDialog(record["root"], self)
            dialog.exec()
        except (OSError, ValueError, TypeError, KeyError) as exc:
            QMessageBox.warning(
                self, "Cannot open analysis configurations", str(exc))

    def _open_ned(self):
        item = self.projects.currentItem()
        if item is None or item.isHidden() or self.workspace is not None:
            return
        record = item.data(Qt.UserRole)
        if (not isinstance(record, dict) or not record.get("available")
                or record.get("legacy")):
            return
        try:
            dialog = NedKnowledgeDialog(record["root"], self)
            dialog.exec()
            self.refresh()
        except (OSError, ValueError, TypeError, KeyError) as exc:
            QMessageBox.warning(self, "Cannot open NED Knowledge", str(exc))

    def _open_knowledge(self):
        item = self.projects.currentItem()
        if item is None or item.isHidden() or self.workspace is not None:
            return
        record = item.data(Qt.UserRole)
        if (not isinstance(record, dict) or not record.get("available")
                or record.get("legacy")):
            return
        try:
            dialog = GalaxyKnowledgeDialog(record["root"], self)
            dialog.exec()
            # Rebuild aliases/confirmed-ID search and the selected card from
            # immutable local evidence. No network, FITS or science jobs.
            self.refresh()
        except (OSError, ValueError, TypeError, KeyError) as exc:
            QMessageBox.warning(self, "Cannot open Galaxy Knowledge", str(exc))

    def _open_selected(self):
        item = self.projects.currentItem()
        if not item or item.isHidden():
            return
        row = item.data(Qt.UserRole)
        if not row["available"]:
            folder = QFileDialog.getExistingDirectory(
                self, "Relink moved project", str(Path.home()))
            if folder:
                try:
                    project = relink_object_project(row["project_id"], folder)
                    self._enter_project(project)
                except (OSError, ValueError) as exc:
                    QMessageBox.warning(self, "Cannot relink project", str(exc))
            return
        if row.get("legacy"):
            self.import_folder(row["root"])
            return
        try:
            self._enter_project(read_object_project(row["root"]))
        except (OSError, ValueError, KeyError) as exc:
            QMessageBox.warning(self, "Cannot open project", str(exc))

    def _enter_project(self, project, *, show_remote=False):
        if self.workspace is not None:
            return
        try:
            touch_object_project(project)
            workspace = SurveyExplorer(
                restore_session=False, launch_remote=False)
            workspace.attach_object_project(project)
            if show_remote and project.field is not None:
                # The user explicitly created a target; display the original
                # remote Sky Explorer, not a vacant FITS panel. Reopening an
                # existing project remains offline unless requested.
                workspace.preview_tabs.setCurrentWidget(workspace.sky)
                workspace.inspector_tabs.setCurrentWidget(
                    workspace.inspector_tabs.widget(0))
                workspace.show_viewer()
            from lsbpipe.surveys.operation_journal import interrupted_operations
            interrupted = interrupted_operations(project.root)
            if interrupted:
                workspace.project_status.setText(
                    f"{len(interrupted)} uncertain prior operation(s) · "
                    "no automatic replay or science adoption")
                workspace.message.setText(
                    "The previous session left unfinished operation records. "
                    "Inspect the operations/ directory and immutable QA "
                    "history before using their outputs. Original FITS "
                    "were not modified by recovery.")
        except (OSError, ValueError, KeyError, RuntimeError) as exc:
            if "workspace" in locals():
                workspace.close()
                workspace.deleteLater()
            self.hint.setText(f"Could not open workspace: {exc}")
            QMessageBox.warning(
                self, "Cannot open galaxy workspace", str(exc))
            return
        self.workspace = workspace
        workspace.returnToProjects.connect(self._return_from_workspace)
        # Windowed maximized mode retains OS taskbar and title bar.
        # Do not use frameless/full-screen presentation.
        workspace.showMaximized()
        if os.environ.get("DESIDERA_UI_DIAG") == "1":
            from .ui_diagnostics import capture_workspace_corner
            QTimer.singleShot(
                1100,
                lambda w=workspace: (
                    capture_workspace_corner(w)
                    if self.workspace is w and w.isVisible() else None))
        self.hide()
        self.projectOpened.emit(project.project_id)

    def _return_from_workspace(self):
        window = self.workspace
        if window is None:
            return
        if window.thread is not None:
            QMessageBox.information(
                window, "Operation in progress",
                "Finish the running task before switching projects.")
            return
        window._returning_to_hub = True
        window._save_project_state()
        self.show()
        window.close()
        window.deleteLater()
        self.workspace = None
        self.refresh()


def main():
    app = QApplication.instance() or QApplication([])
    hub = ProjectHub()
    hub.showMaximized()
    return app.exec()
