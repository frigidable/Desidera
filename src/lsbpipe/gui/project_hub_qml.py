"""Opt-in Qt Quick Project Hub comparison, backed by the LIVE Widgets library.

This is an experimental presentation only: a second registry, project
serializer or scientific-state estimator must never be introduced here.
Qt Quick is loaded only on explicit user request; normal startup and all
existing FITS/QtWebEngine workspaces remain Qt Widgets.
"""
from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import (
    QAbstractListModel, QModelIndex, QObject, Property, Qt, QUrl, Signal, Slot,
)
from PySide6.QtQuick import QQuickView
from PySide6.QtQuickControls2 import QQuickStyle

from .visual_theme import theme_colors


_QML_STYLE_READY = False


class HubProjectListModel(QAbstractListModel):
    """Read-only QML projection of the *existing* Project Hub list roles."""

    GalaxyName = Qt.UserRole + 10
    StatusText = Qt.UserRole + 11
    Coordinates = Qt.UserRole + 12
    ProjectId = Qt.UserRole + 13
    Available = Qt.UserRole + 14

    _NAMES = {
        GalaxyName: b"galaxyName",
        StatusText: b"statusText",
        Coordinates: b"coordinates",
        ProjectId: b"projectId",
        Available: b"isAvailable",
    }

    def __init__(self, hub, parent=None):
        super().__init__(parent)
        self._hub = hub
        self._rows = []
        self.sync()

    def roleNames(self):
        return self._NAMES.copy()

    def rowCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self._rows)

    def data(self, index, role=Qt.DisplayRole):
        if not index.isValid() or not 0 <= index.row() < len(self._rows):
            return None
        record = self._rows[index.row()]
        if role == Qt.DisplayRole:
            return record["galaxyName"]
        name = self._NAMES.get(role)
        return record.get(name.decode("ascii")) if name is not None else None

    def sync(self, *_):
        """No FITS, catalogue, files or network: reuse visible Widgets rows."""
        rows = []
        for index in range(self._hub.projects.count()):
            item = self._hub.projects.item(index)
            if item.isHidden():
                continue
            name, status, coords = item.data(Qt.UserRole + 2)
            record = item.data(Qt.UserRole)
            rows.append({
                "galaxyName": name,
                "statusText": status,
                "coordinates": coords,
                "projectId": record["project_id"],
                "isAvailable": bool(record.get("available")),
            })
        self.beginResetModel()
        self._rows = rows
        self.endResetModel()


class HubQmlBridge(QObject):
    """QML commands delegate to the SAME tested ProjectHub functions."""

    changed = Signal()

    def __init__(self, hub, model, window, parent=None):
        super().__init__(parent)
        self._hub = hub
        self._model = model
        self._window = window

    @Property(bool, notify=changed)
    def darkTheme(self):
        return theme_colors(self._hub)["BG"] == "#0C1424"

    @Property(str, notify=changed)
    def selectedId(self):
        item = self._hub.projects.currentItem()
        if item is None or item.isHidden():
            return ""
        record = item.data(Qt.UserRole)
        return str(record.get("project_id", ""))

    @Property(str, notify=changed)
    def selectedName(self):
        snap = self._hub.galaxy_card.snapshot
        return snap.display_name if snap is not None else "Choose a galaxy"

    @Property(str, notify=changed)
    def selectedCoords(self):
        snap = self._hub.galaxy_card.snapshot
        return snap.coordinates if snap is not None else "Select a saved project"

    @Property(str, notify=changed)
    def selectedStatus(self):
        item = self._hub.projects.currentItem()
        if item is None or item.isHidden():
            return "No project selected"
        fields = item.data(Qt.UserRole + 2)
        return str(fields[1]) if fields else "Status unavailable"

    @Property(str, notify=changed)
    def datasetLabel(self):
        snap = self._hub.galaxy_card.snapshot
        return snap.dataset if snap is not None else "Not selected"

    @Property(str, notify=changed)
    def fitsLabel(self):
        snap = self._hub.galaxy_card.snapshot
        return snap.fits_status if snap is not None else "No local FITS evidence"

    @Property(str, notify=changed)
    def nextAction(self):
        snap = self._hub.galaxy_card.snapshot
        if snap is None:
            return "Select a galaxy to inspect saved data and preparation."
        if not snap.available:
            return "Relink the original project folder; do not substitute another galaxy."
        if snap.readiness is not None:
            return snap.readiness.next_action_label
        return "Open the project to inspect its data."

    @Property(str, notify=changed)
    def openLabel(self):
        return ("Select a galaxy" if not self.selectedId else
                self._hub.galaxy_card.open_button.text())

    @Property(str, notify=changed)
    def catalogueLabel(self):
        """Only receipt-confirmed identities are named as associations."""
        snap = self._hub.galaxy_card.snapshot
        if snap is None:
            return "No galaxy selected"
        if snap.knowledge_status == "EVIDENCE_ERROR":
            return "Catalogue evidence requires inspection"
        identities = snap.catalogue_identities or (
            (snap.catalogue_identity,) if snap.catalogue_identity else ())
        if identities:
            return "  ·  ".join(
                f"{entry.provider}: {entry.main_id}" for entry in identities)
        if snap.knowledge_status == "UNVERIFIED_CANDIDATES_ONLY":
            return "Unverified candidates · no confirmed identity"
        return "No confirmed SIMBAD / NED association"

    @Property(str, notify=changed)
    def catalogueEvidence(self):
        card = self._hub.galaxy_card
        snap = card.snapshot
        if snap is None:
            return "Select a project to inspect catalogue sources."
        if snap.knowledge_status == "EVIDENCE_ERROR":
            return "Stored catalogue evidence unreadable: " + snap.knowledge_error
        if not card.catalogue_provenance.isHidden():
            return card.catalogue_provenance.text()
        return ("No confirmed catalogue identity. Remote queries are opt-in; "
                "project coordinates remain unchanged.")

    @Property(str, notify=changed)
    def literatureLabel(self):
        label = self._hub.galaxy_card.literature
        return label.text() if not label.isHidden() else ""

    @Property(str, notify=changed)
    def qaSummary(self):
        label = self._hub.galaxy_card.qa_status
        return label.text() if not label.isHidden() else (
            "No saved FITS/WCS audit for this selection")

    @Property(str, notify=changed)
    def qaEvidence(self):
        """Expanded source evidence, independent of the Widgets disclosure."""
        card = self._hub.galaxy_card
        snap = card.snapshot
        if snap is None:
            return ""
        lines = []
        if snap.audit is not None:
            audit = snap.audit
            lines.extend((card.qa_status.text() or audit.label,
                          audit.explanation))
            if audit.receipt_relative_path:
                lines.append(
                    "Saved audit receipt: " + audit.receipt_relative_path
                    + (" · " + audit.created_utc if audit.created_utc else ""))
        if snap.readiness is not None:
            lines.append(snap.readiness.explanation)
        return "\n".join(line for line in lines if line) or (
            "No current FITS/WCS diagnostic is documented for this selection.")

    @Property(str, notify=changed)
    def projectEvidence(self):
        card = self._hub.galaxy_card
        if card.snapshot is None:
            return ""
        lines = []
        for label in (card.aliases, card.last_opened, card.path, card.note):
            if not label.isHidden() and label.text():
                lines.append(label.text())
        return "\n".join(lines)

    @Property(str, notify=changed)
    def previewUrl(self):
        """A cached display-only preview, never an original FITS path."""
        snap = self._hub.galaxy_card.snapshot
        if snap is None or snap.thumbnail is None:
            return ""
        return QUrl.fromLocalFile(str(snap.thumbnail)).toString()

    @Property(str, notify=changed)
    def previewLabel(self):
        snap = self._hub.galaxy_card.snapshot
        if snap is None:
            return "No image selected"
        if snap.preview_source == "REMOTE_LEGACY_RGB":
            return "Cached Legacy RGB · display only"
        if snap.preview_source == "LOCAL_FITS_GRAYSCALE":
            return "Cached local FITS thumbnail · display only"
        return "No cached preview · no science result"

    @Property(bool, notify=changed)
    def canCatalogue(self):
        return self._hub.catalogues_menu_btn.isEnabled()

    @Property(bool, notify=changed)
    def canConfigure(self):
        return self._hub.config_action.isEnabled()

    @Property(bool, notify=changed)
    def canLocalPreview(self):
        return self._hub.local_preview_action.isEnabled()

    @Property(bool, notify=changed)
    def canRemotePreview(self):
        return self._hub.rgb_action.isEnabled()

    @Property(bool, notify=changed)
    def canRemove(self):
        return self._hub.remove_action.isEnabled()

    @Property(int, notify=changed)
    def selectedIndex(self):
        """Current visible row, without changing the saved galaxy selection."""
        selected = self.selectedId
        return next((i for i, row in enumerate(self._model._rows)
                     if row["projectId"] == selected), -1)

    @Slot(int, result=str)
    def idAt(self, row):
        if not 0 <= row < self._model.rowCount():
            return ""
        return self._model._rows[row]["projectId"]

    @Property(int, notify=changed)
    def projectCount(self):
        return self._hub.projects.count()

    @Property(int, notify=changed)
    def visibleCount(self):
        return self._model.rowCount()

    @Slot(str)
    def setFilter(self, query):
        if self._hub.filter.text() != query:
            self._hub.filter.setText(query)

    @Slot(str)
    def choose(self, project_id):
        for row in range(self._hub.projects.count()):
            item = self._hub.projects.item(row)
            record = item.data(Qt.UserRole)
            if (not item.isHidden()
                    and str(record.get("project_id")) == project_id):
                self._hub.projects.setCurrentRow(row)
                self.changed.emit()
                return

    @Slot()
    def openSelected(self):
        project_id = self.selectedId
        if not project_id:
            return
        # Project ID must still match the native current selection.
        self._hub._open_card_project(project_id)

    @Slot(str)
    def runTool(self, command):
        """Route an explicitly requested action through the native guards.

        The experiment must not invent a second network/download flow or
        bypass UUID, readiness or confirmation checks of Project Hub.
        """
        if not self.selectedId or self._hub.workspace is not None:
            return
        allowed = {
            "simbad": (self._hub.simbad_action, self._hub._open_knowledge),
            "ned": (self._hub.ned_action, self._hub._open_ned),
            "config": (self._hub.config_action, self._hub._open_configurations),
            "localPreview": (self._hub.local_preview_action,
                             self._hub._regenerate_preview),
            "remotePreview": (self._hub.rgb_action,
                              self._hub._fetch_rgb_preview),
            "remove": (self._hub.remove_action, self._hub._remove_selected),
        }
        route = allowed.get(command)
        if route is None or not route[0].isEnabled():
            return
        self._window.hide()
        self._hub.show()
        self._hub.raise_()
        try:
            route[1]()
        finally:
            if self._hub.workspace is None:
                self._window.sync()
                self._window.show()
                self._window.raise_()
                self._window.requestActivate()

    @Slot()
    def createProject(self):
        self._window.close()
        self._hub.show()
        self._hub.raise_()
        self._hub._new_project()

    @Slot()
    def importFolder(self):
        self._window.close()
        self._hub.show()
        self._hub.raise_()
        self._hub._import_dialog()

    @Slot()
    def closePreview(self):
        self._window.close()


class QmlProjectHubPreview(QQuickView):
    """Isolated Qt Quick comparison window; never replaces the real hub."""

    def __init__(self, hub):
        # Qt's native Windows controls ignore customized background/contentItem.
        # This experimental window owns its QML scene; it never changes
        # Qt Widgets/FITS application palettes or scientific rendering.
        global _QML_STYLE_READY
        if not _QML_STYLE_READY:
            QQuickStyle.setStyle("Basic")
            _QML_STYLE_READY = True
        super().__init__()
        self.setObjectName("desidera_qml_hub_preview")
        self.setTitle("DESIDERA · Observatory Studio · experimental Qt Quick")
        self.setResizeMode(QQuickView.SizeRootObjectToView)
        self.model = HubProjectListModel(hub, self)
        self.bridge = HubQmlBridge(hub, self.model, self, self)
        self.rootContext().setContextProperty("projectModel", self.model)
        self.rootContext().setContextProperty("hubBridge", self.bridge)
        source = Path(__file__).with_name("project_hub_preview.qml")
        self.setSource(QUrl.fromLocalFile(str(source)))
        if self.status() != QQuickView.Ready:
            errors = "; ".join(str(error) for error in self.errors())
            raise RuntimeError("QML Project Hub failed to load: " + errors)
        self.resize(1180, 740)
        hub.filter.textChanged.connect(self.sync)
        hub.projects.itemSelectionChanged.connect(self.bridge.changed.emit)
        hub.projectOpened.connect(self.close)
        hub.destroyed.connect(self.close)
        self.sync()

    def sync(self, *_):
        self.model.sync()
        self.bridge.changed.emit()
