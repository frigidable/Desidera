"""Project-bound, offline Dataset Manager: inspect real files, not fake releases.

Registry reading is cheap. Heavy FITS/metadata operations run via the
existing workspace background-task controller; no downloads or coaddition.
"""
from __future__ import annotations

from collections import defaultdict

from PySide6.QtCore import Qt, Signal, QPointF
from PySide6.QtGui import QBrush, QColor, QPen, QPolygonF, QPainter
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QSplitter, QLabel, QPushButton,
    QTreeWidget, QTreeWidgetItem, QGraphicsScene, QGraphicsView,
    QGraphicsPolygonItem, QGraphicsRectItem, QGraphicsEllipseItem,
    QHeaderView, QFrame,
)


SOURCE_COLORS = (
    "#64A5E2", "#E5A75C", "#8BC6A2", "#CC91DB",
    "#E07983", "#91B8DA", "#D2C877", "#77BFD1",
)


class _SkyFootprintView(QGraphicsView):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setDragMode(QGraphicsView.ScrollHandDrag)
        self.setRenderHint(QPainter.Antialiasing, True)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.setObjectName("dataset_coverage_map_view")

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if not self.scene().sceneRect().isEmpty():
            self.fitInView(self.scene().sceneRect(), Qt.KeepAspectRatio)


class DatasetManager(QWidget):
    """A view of *registered* source provenance with explicit P0 QA limits."""

    registerRequested = Signal(str)
    coverageRequested = Signal()
    checksumRequested = Signal(str)
    backRequested = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("dataset_manager")
        self.project = None
        self.index = None
        self.coverage = None
        self._selected_source_id = None
        self._colors = {}
        root = QVBoxLayout(self)
        root.setContentsMargins(12, 12, 12, 12)
        root.setSpacing(9)

        heading = QHBoxLayout()
        title = QLabel("Data Manager")
        title.setStyleSheet("font-size: 20px; font-weight: 600")
        heading.addWidget(title)
        heading.addStretch(1)
        self.back_btn = QPushButton("Back to FITS Explorer")
        self.back_btn.setObjectName("dataset_manager_back")
        self.back_btn.clicked.connect(self.backRequested.emit)
        heading.addWidget(self.back_btn)
        root.addLayout(heading)
        subtitle = QLabel(
            "Original observational sources for this galaxy · "
            "separate releases, bricks and filters · no scientific adoption")
        subtitle.setWordWrap(True)
        root.addWidget(subtitle)

        actions = QHBoxLayout()
        self.import_btn = QPushButton("Index local FITS folder…")
        self.import_btn.setObjectName("dataset_manager_import")
        self.import_btn.setToolTip(
            "Offline metadata-only registration. Existing FITS stay where "
            "they are; no download and no automatic scientific acceptance.")
        self.import_btn.clicked.connect(self._request_register)
        actions.addWidget(self.import_btn)
        self.map_btn = QPushButton("Inspect WCS footprints")
        self.map_btn.setObjectName("dataset_manager_coverage")
        self.map_btn.setToolTip(
            "Read original IMAGE FITS headers only; this is geometric "
            "coverage, not valid data fraction or photometric quality.")
        self.map_btn.clicked.connect(self.coverageRequested.emit)
        actions.addWidget(self.map_btn)
        self.hash_btn = QPushButton("Verify selected file SHA-256")
        self.hash_btn.setObjectName("dataset_manager_hash")
        self.hash_btn.setToolTip(
            "Compute SHA-256 without modifying the source file. Does not "
            "establish FITS validity or scientific quality.")
        self.hash_btn.clicked.connect(self._request_hash)
        actions.addWidget(self.hash_btn)
        actions.addStretch(1)
        root.addLayout(actions)

        split = QSplitter(Qt.Horizontal)
        sources_widget = QWidget()
        source_layout = QVBoxLayout(sources_widget)
        source_layout.setContentsMargins(0, 0, 0, 0)
        source_layout.setSpacing(5)
        source_layout.addWidget(QLabel("REGISTERED SOURCES · LOCAL FILES"))
        self.tree = QTreeWidget()
        self.tree.setObjectName("dataset_manager_sources")
        self.tree.setHeaderLabels(["Source / original product", "Band", "Evidence"])
        self.tree.header().setSectionResizeMode(0, QHeaderView.Stretch)
        self.tree.header().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        self.tree.header().setSectionResizeMode(2, QHeaderView.ResizeToContents)
        self.tree.setRootIsDecorated(True)
        self.tree.setAlternatingRowColors(True)
        self.tree.currentItemChanged.connect(self._selected)
        source_layout.addWidget(self.tree, 1)
        self.source_summary = QLabel(
            "No registered sources. Existing original FITS can be indexed "
            "in place without downloading them again.")
        self.source_summary.setObjectName("dataset_manager_source_summary")
        self.source_summary.setWordWrap(True)
        source_layout.addWidget(self.source_summary)
        split.addWidget(sources_widget)

        map_panel = QFrame()
        map_panel.setFrameShape(QFrame.StyledPanel)
        map_layout = QVBoxLayout(map_panel)
        map_layout.setContentsMargins(8, 8, 8, 8)
        map_layout.addWidget(QLabel("SKY FOOTPRINTS · ICRS · north up, east left"))
        self.scene = QGraphicsScene(self)
        self.map_view = _SkyFootprintView()
        self.map_view.setScene(self.scene)
        map_layout.addWidget(self.map_view, 1)
        self.map_status = QLabel(
            "Footprints not checked. Select Inspect WCS footprints to "
            "read actual local FITS headers.")
        self.map_status.setObjectName("dataset_manager_map_status")
        self.map_status.setWordWrap(True)
        map_layout.addWidget(self.map_status)
        split.addWidget(map_panel)
        split.setStretchFactor(0, 3)
        split.setStretchFactor(1, 4)
        root.addWidget(split, 1)

        self.details = QLabel(
            "Dataset evidence: file present ≠ valid FITS ≠ photometry-ready. "
            "Joint coaddition is disabled without exposure lineage and "
            "calibration validation.")
        self.details.setObjectName("dataset_manager_details")
        self.details.setWordWrap(True)
        root.addWidget(self.details)
        self._draw_coverage()

    def bind_project(self, project, index: dict):
        self.project = project
        self.set_index(index)

    def set_index(self, index: dict):
        self.index = index
        self.coverage = None
        self.tree.blockSignals(True)
        self.tree.clear()
        self._colors.clear()
        sources = index.get("sources", [])
        assets = index.get("assets", [])
        by_source = defaultdict(list)
        for a in assets:
            by_source[a["source_id"]].append(a)
        selected_id = self._selected_source_id
        for n, source in enumerate(sources):
            sid = source["source_id"]
            self._colors[sid] = SOURCE_COLORS[n % len(SOURCE_COLORS)]
            label = " · ".join(str(x) for x in (
                source.get("provider", "External"), source.get("release"),
                source.get("region")) if x)
            children = by_source.get(sid, [])
            bands = ", ".join(sorted({
                a["band"] for a in children if a.get("band")}))
            item = QTreeWidgetItem([
                f"{label} · {len(children)} files",
                bands or "—",
                source.get("calibration_status", "UNASSESSED"),
            ])
            item.setData(0, Qt.UserRole, ("source", sid))
            self.tree.addTopLevelItem(item)
            for asset in sorted(children, key=lambda a: (
                    str(a.get("brick") or ""),
                    str(a.get("band") or ""),
                    str(a.get("product_kind") or ""))):
                name = " / ".join(str(x) for x in (
                    asset.get("brick"), asset.get("product_kind")) if x)
                label = (name or "External FITS")
                status = asset.get("checksum_status", "NOT_HASHED")
                if asset.get("fits_status") != "FITS_WCS_VALIDATED_DIAGNOSTIC":
                    status += " · FITS unverified"
                child = QTreeWidgetItem([
                    label, str(asset.get("band") or "shared"),
                    status])
                child.setData(0, Qt.UserRole, ("asset", asset["asset_id"], sid))
                child.setToolTip(0, str(asset.get("path_ref", {})))
                item.addChild(child)
            if selected_id == sid:
                self.tree.setCurrentItem(item)
        self.tree.blockSignals(False)
        if sources and not self.tree.currentItem():
            self.tree.setCurrentItem(self.tree.topLevelItem(0))
        self.source_summary.setText(
            f"{len(sources)} source groups · {len(assets)} original files "
            "· stored evidence only · not automatically adopted"
            if sources else
            "No indexed sources. Use Index local FITS folder to register "
            "existing Legacy files without copying or downloading.")
        self.hash_btn.setEnabled(False)
        self.map_btn.setEnabled(bool(self.project and self.project.field and assets))
        self._draw_coverage()

    def _selected(self, item, _previous):
        meta = item.data(0, Qt.UserRole) if item is not None else None
        if not meta:
            self.hash_btn.setEnabled(False)
            return
        self._selected_source_id = meta[-1]
        self.hash_btn.setEnabled(meta[0] == "asset")
        if self.index and meta[0] == "asset":
            record = next(
                (x for x in self.index["assets"] if x["asset_id"] == meta[1]), None)
            if record:
                self.details.setText(
                    f"Original product: {record.get('product_kind')} · "
                    f"SHA-256: {record.get('sha256') or 'not yet verified'} · "
                    f"FITS/WCS: {record.get('fits_status')} · "
                    "no PSF/sky/LSB suitability certification")
        else:
            self.details.setText(
                "Selected source: independent lineage UNKNOWN; "
                "photometric calibration UNASSESSED. Primary/compare "
                "analysis configurations will be created separately.")
        if self.coverage is not None:
            self._draw_coverage()

    def _request_register(self):
        from PySide6.QtWidgets import QFileDialog
        base = str(self.project.data_root or self.project.root) if self.project else ""
        root = QFileDialog.getExistingDirectory(
            self, "Index existing original Legacy FITS directory", base)
        if root:
            self.registerRequested.emit(root)

    def _request_hash(self):
        item = self.tree.currentItem()
        meta = item.data(0, Qt.UserRole) if item else None
        if meta and meta[0] == "asset":
            self.checksumRequested.emit(meta[1])

    def set_busy(self, busy: bool):
        for button in (self.import_btn, self.map_btn, self.hash_btn):
            button.setEnabled(not busy)
        if not busy:
            self.hash_btn.setEnabled(
                bool(self.tree.currentItem() and
                     self.tree.currentItem().data(0, Qt.UserRole) and
                     self.tree.currentItem().data(0, Qt.UserRole)[0] == "asset"))

    def set_coverage(self, result: dict):
        self.coverage = result
        self._draw_coverage()

    def _draw_coverage(self):
        self.scene.clear()
        if self.project is None or self.project.field is None:
            self.map_status.setText(
                "ICRS coordinates required for FITS coverage geometry.")
            return
        width = self.project.field.width_arcmin * 60
        field_rect = (-width/2, -width/2, width, width)
        pen = QPen(QColor("#8794A8"), 1.2, Qt.DashLine)
        pen.setCosmetic(True)
        frame = self.scene.addRect(*field_rect, pen)
        frame.setZValue(9)
        radii = [
            (self.project.sky_outer_arcsec, "#A1A9B7"),
            (self.project.sky_inner_arcsec, "#A1A9B7"),
            (self.project.target_radius_arcsec, "#DBB867"),
        ]
        for rad, color in radii:
            circle = self.scene.addEllipse(
                -rad, -rad, 2*rad, 2*rad,
                QPen(QColor(color), 1.1, Qt.DashLine))
            circle.setZValue(10)
        cross = self.scene.addEllipse(
            -2, -2, 4, 4, QPen(QColor("#FFFFFF"), 1),
            QBrush(QColor("#FFFFFF")))
        cross.setZValue(11)
        if self.coverage is not None:
            for record in self.coverage.get("footprints", []):
                points = record["points"]
                if len(points) < 3:
                    continue
                selected = self._selected_source_id
                highlight = selected is None or record["source_id"] == selected
                color = QColor(self._colors.get(record["source_id"], "#81A8CE"))
                fill = QColor(color)
                fill.setAlpha(35 if highlight else 7)
                polygon = QGraphicsPolygonItem(
                    QPolygonF([QPointF(*xy) for xy in points]))
                polygon.setPen(QPen(color if highlight else
                                    QColor("#596474"),
                                    1.4 if highlight else 0.5))
                polygon.setBrush(QBrush(fill))
                polygon.setZValue(1 if highlight else 0)
                label = (record.get("release") or
                         record.get("provider") or "External")
                polygon.setToolTip(
                    f"{label} · {record.get('region') or '—'} · "
                    f"{record.get('brick') or 'external'} · "
                    f"{record.get('band') or '—'} · WCS header only")
                self.scene.addItem(polygon)
            good = len(self.coverage.get("footprints", []))
            rejected = len(self.coverage.get("rejected", []))
            omitted = self.coverage.get("omitted_count", 0)
            self.map_status.setText(
                f"Original IMAGE WCS headers: {good} footprints, "
                f"{rejected} invalid/missing, {omitted} omitted by "
                "display cap. Geometry only; no valid-pixel or sky evidence.")
        else:
            self.map_status.setText(
                "ICRS target and sky rings shown. FITS sky footprints "
                "not inspected; no geometric coverage claim.")
        bounds = self.scene.itemsBoundingRect()
        bounds.adjust(-width*.08, -width*.08, width*.08, width*.08)
        self.scene.setSceneRect(bounds)
        if self.map_view.viewport().width() > 0:
            self.map_view.fitInView(bounds, Qt.KeepAspectRatio)
