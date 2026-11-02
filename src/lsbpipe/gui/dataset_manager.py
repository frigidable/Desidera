"""Project-bound, offline Dataset Manager: inspect real files, not fake releases.

Registry reading is cheap. Heavy FITS/metadata operations run via the
existing workspace background-task controller; no downloads or coaddition.
"""
from __future__ import annotations

from collections import defaultdict

from PySide6.QtCore import Qt, Signal, QPointF, QRectF
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


def _legacy_quartets(assets: list[dict]) -> tuple[int, int, list[tuple[str, str, str]]]:
    """Index-presence check, NOT FITS/WCS or scientific completeness.

    IMAGE/INVVAR/NEXP belong to each brick-band, while one MASKBITS
    may be shared between all bands in the same brick.
    """
    bricks = defaultdict(lambda: {"masks": False, "bands": defaultdict(set)})
    for asset in assets:
        brick = asset.get("brick")
        if not brick:
            continue
        kind = asset.get("product_kind")
        if kind == "maskbits":
            bricks[brick]["masks"] = True
        elif asset.get("band"):
            bricks[brick]["bands"][asset["band"]].add(kind)
    total = complete = 0
    missing = []
    for brick, data in sorted(bricks.items()):
        for band, products in sorted(data["bands"].items()):
            total += 1
            required = {"image", "invvar", "nexp"}
            if not data["masks"]:
                required.add("maskbits")
            absent = sorted(required - products)
            if not absent:
                complete += 1
            else:
                missing.append((brick, band, ", ".join(absent)))
    return complete, total, missing


class _SkyFootprintView(QGraphicsView):
    resized = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setDragMode(QGraphicsView.ScrollHandDrag)
        self.setRenderHint(QPainter.Antialiasing, True)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.setObjectName("dataset_coverage_map_view")

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self.resized.emit()


class DatasetManager(QWidget):
    """A view of *registered* source provenance with explicit P0 QA limits."""

    registerRequested = Signal(str)
    externalRequested = Signal(str, str, str)
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
        self._selected_asset_id = None
        self._busy = False
        self._map_fit_mode = "target"
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
        self.external_btn = QPushButton("Attach external FITS…")
        self.external_btn.setObjectName("dataset_manager_external")
        self.external_btn.setToolTip(
            "Register one existing FITS with an explicit provider and band. "
            "No copying, calibration assumption, or network access.")
        self.external_btn.clicked.connect(self._request_external)
        actions.addWidget(self.external_btn)
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
        map_tools = QHBoxLayout()
        map_tools.addWidget(QLabel("SKY FOOTPRINTS · ICRS · N up, E left"), 1)
        self.fit_target_btn = QPushButton("Target + sky")
        self.fit_target_btn.setObjectName("dataset_map_fit_target")
        self.fit_target_btn.setToolTip(
            "Focus on galaxy field and diagnostic sky annulus, "
            "which may crop outer parts of a large original FITS frame.")
        self.fit_target_btn.clicked.connect(
            lambda: self._set_map_fit_mode("target"))
        map_tools.addWidget(self.fit_target_btn)
        self.fit_all_btn = QPushButton("All footprints")
        self.fit_all_btn.setObjectName("dataset_map_fit_all")
        self.fit_all_btn.setToolTip(
            "Fit every inspected original WCS footprint and the target region.")
        self.fit_all_btn.clicked.connect(
            lambda: self._set_map_fit_mode("all"))
        map_tools.addWidget(self.fit_all_btn)
        map_layout.addLayout(map_tools)
        self.scene = QGraphicsScene(self)
        self.map_view = _SkyFootprintView()
        self.map_view.setScene(self.scene)
        self.map_view.resized.connect(self._fit_map)
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
        if (self.project is None or
                self.project.project_id != project.project_id):
            self._selected_source_id = None
            self._selected_asset_id = None
            self._map_fit_mode = "target"
        self.project = project
        self.set_index(index)

    def set_index(self, index: dict):
        self.index = index
        self.coverage = None
        selected_item = self.tree.currentItem()
        selected_meta = (
            selected_item.data(0, Qt.UserRole) if selected_item else None)
        if selected_meta and selected_meta[0] == "asset":
            self._selected_asset_id = selected_meta[1]
        self.tree.blockSignals(True)
        self.tree.clear()
        self._colors.clear()
        sources = index.get("sources", [])
        assets = index.get("assets", [])
        by_source = defaultdict(list)
        for a in assets:
            by_source[a["source_id"]].append(a)
        selected_id = self._selected_source_id
        restore_item = None
        incomplete_count = 0
        for n, source in enumerate(sources):
            sid = source["source_id"]
            self._colors[sid] = SOURCE_COLORS[n % len(SOURCE_COLORS)]
            label = " · ".join(str(x) for x in (
                source.get("provider", "External"), source.get("release"),
                source.get("region")) if x)
            children = by_source.get(sid, [])
            bands = ", ".join(sorted({
                a["band"] for a in children if a.get("band")}))
            legacy = source.get("provider") == "legacy_surveys"
            complete, total, missing = (
                _legacy_quartets(children) if legacy else (0, 0, []))
            incomplete_count += len(missing)
            evidence = (
                f"{complete}/{total} quartets indexed"
                if legacy and total else
                "Source unassessed")
            item = QTreeWidgetItem([
                f"{label} · {len(children)} files",
                bands or "—", evidence,
            ])
            item.setToolTip(
                2, "Index presence only, not FITS/WCS or photometric QA.")
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
                if asset["asset_id"] == self._selected_asset_id:
                    restore_item = child
            for brick, band, absent in missing:
                warning = QTreeWidgetItem([
                    f"{brick} / missing {absent}",
                    band, "NOT INDEXED · no FITS file claim"])
                warning.setData(0, Qt.UserRole, ("missing", sid))
                warning.setForeground(0, QBrush(QColor("#E0B66C")))
                item.addChild(warning)
                item.setExpanded(True)
            if selected_id == sid and restore_item is None:
                restore_item = item
        if restore_item is None and sources:
            restore_item = self.tree.topLevelItem(0)
        if restore_item is not None:
            self.tree.setCurrentItem(restore_item)
        self.tree.blockSignals(False)
        self._selected(self.tree.currentItem(), None)
        self.source_summary.setText(
            f"{len(sources)} source groups · {len(assets)} registered FITS · "
            f"{incomplete_count} incomplete indexed brick/band sets · "
            "FITS/WCS & LSB suitability NOT certified"
            if sources else
            "No indexed sources. Use Index local FITS folder to register "
            "existing Legacy files without copying or downloading.")
        self.hash_btn.setEnabled(
            bool(not self._busy and self.tree.currentItem() and
                 self.tree.currentItem().data(0, Qt.UserRole) and
                 self.tree.currentItem().data(0, Qt.UserRole)[0] == "asset"))
        self.map_btn.setEnabled(
            bool(not self._busy and self.project and
                 self.project.field and assets))
        self._draw_coverage()

    def _selected(self, item, _previous):
        meta = item.data(0, Qt.UserRole) if item is not None else None
        if not meta:
            self.hash_btn.setEnabled(False)
            return
        self._selected_source_id = meta[-1]
        self._selected_asset_id = (
            meta[1] if meta[0] == "asset" else None)
        self.hash_btn.setEnabled(not self._busy and meta[0] == "asset")
        if meta[0] == "missing":
            self.details.setText(
                "Incomplete indexed Legacy quartet: this original product "
                "is not registered locally. Do not run four-layer "
                "photometric QA until the missing input is supplied.")
        elif self.index and meta[0] == "asset":
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

    def _request_external(self):
        from PySide6.QtWidgets import QFileDialog, QInputDialog
        filename, _ = QFileDialog.getOpenFileName(
            self, "Select existing original FITS", "",
            "FITS images (*.fits *.fit *.fts *.fits.fz)")
        if not filename:
            return
        provider, accepted = QInputDialog.getText(
            self, "Observation source", "Survey / facility / instrument:")
        if not accepted or not provider.strip():
            return
        band, accepted = QInputDialog.getText(
            self, "Passband", "Native filter or passband name:")
        if accepted and band.strip():
            self.externalRequested.emit(filename, provider.strip(), band.strip())

    def _request_hash(self):
        item = self.tree.currentItem()
        meta = item.data(0, Qt.UserRole) if item else None
        if meta and meta[0] == "asset":
            self.checksumRequested.emit(meta[1])

    def set_busy(self, busy: bool):
        self._busy = bool(busy)
        for button in (
                self.import_btn, self.external_btn, self.map_btn, self.hash_btn):
            button.setEnabled(not busy)
        if not busy:
            self.map_btn.setEnabled(bool(
                self.project and self.project.field and self.index and
                self.index.get("assets")))
            self.hash_btn.setEnabled(
                bool(self.tree.currentItem() and
                     self.tree.currentItem().data(0, Qt.UserRole) and
                     self.tree.currentItem().data(0, Qt.UserRole)[0] == "asset"))

    def set_coverage(self, result: dict):
        self.coverage = result
        self._draw_coverage()

    def _set_map_fit_mode(self, mode: str):
        self._map_fit_mode = mode
        self._fit_map()

    def _fit_map(self):
        if self.project is None or self.project.field is None:
            return
        if min(self.map_view.viewport().width(),
               self.map_view.viewport().height()) < 30:
            return
        width = float(self.project.field.width_arcmin * 60)
        rad = max(width/2., float(self.project.sky_outer_arcsec),
                  float(self.project.target_radius_arcsec))
        # Target mode keeps the whole sky ring visible, without reducing the
        # galaxy to a tiny dot in a much larger FITS footprint.
        bounds = (
            QRectF(-rad * 1.16, -rad * 1.16,
                   2 * rad * 1.16, 2 * rad * 1.16)
            if self._map_fit_mode == "target"
            else self.scene.sceneRect())
        if not bounds.isEmpty():
            self.map_view.fitInView(bounds, Qt.KeepAspectRatio)

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
            # Distinct IMAGE assets may share exactly the same cutout WCS.
            # Merge *identical* sky polygons visually within a source, but
            # retain every original image record and its own band in the
            # tooltip and header count. Never conflate North/South sources.
            unique = {}
            for record in self.coverage.get("footprints", []):
                points = record.get("points", [])
                if len(points) < 3:
                    continue
                signature = (record["source_id"], tuple(
                    (round(x, 3), round(y, 3)) for x, y in points))
                if signature not in unique:
                    unique[signature] = [record]
                else:
                    unique[signature].append(record)
            for group in unique.values():
                record = group[0]
                points = record["points"]
                selected = self._selected_source_id
                highlight = selected is None or record["source_id"] == selected
                color = QColor(self._colors.get(record["source_id"], "#81A8CE"))
                pen = QPen(
                    color if highlight else QColor("#637183"),
                    1.9 if highlight else 1.0)
                pen.setCosmetic(True)
                polygon = QGraphicsPolygonItem(
                    QPolygonF([QPointF(*xy) for xy in points]))
                polygon.setPen(pen)
                polygon.setBrush(Qt.NoBrush)  # No opaque accumulated fills
                polygon.setZValue(2 if highlight else 1)
                label = (record.get("release") or
                         record.get("provider") or "External")
                bands = ", ".join(sorted({
                    str(x.get("band") or "unknown") for x in group}))
                polygon.setToolTip(
                    f"{label} · {record.get('region') or '—'} · "
                    f"{record.get('brick') or 'external'} · "
                    f"bands {bands} · {len(group)} original IMAGE files · "
                    "WCS header only")
                self.scene.addItem(polygon)
            good = len(self.coverage.get("footprints", []))
            rejected = len(self.coverage.get("rejected", []))
            omitted = self.coverage.get("omitted_count", 0)
            self.map_status.setText(
                f"{good} original IMAGE WCS headers · "
                f"{len(unique)} distinct source outlines · "
                f"{rejected} invalid/missing · {omitted} omitted · "
                "geometry only, not usable-pixel coverage. "
                "Target + sky may crop outer FITS boundaries.")
        else:
            self.map_status.setText(
                "ICRS target and sky rings shown. FITS sky footprints "
                "not inspected; no geometric coverage claim.")
        bounds = self.scene.itemsBoundingRect()
        bounds.adjust(-width*.08, -width*.08, width*.08, width*.08)
        self.scene.setSceneRect(bounds)
        self._fit_map()
