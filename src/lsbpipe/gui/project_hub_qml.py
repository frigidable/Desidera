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

from .visual_theme import theme_colors


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
        if not self.selectedId:
            return
        self._hub._open_selected()  # Same availability/UUID/relink guards.

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
