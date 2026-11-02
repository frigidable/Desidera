"""Interactive four-layer inspection of a single read-only Legacy FITS quartet.

All four native arrays share one WCS/grid. The comparison does not make a
mosaic, resample, photometrically adopt a mask or load any saved PNG.
NumPy RGB conversion is offloaded; QPixmaps are created on the GUI thread.
"""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor

import numpy as np
from PySide6.QtCore import Qt, QTimer, Signal, QPointF
from PySide6.QtGui import QColor, QPainterPath, QPen, QPixmap
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QCheckBox,
    QLabel, QPushButton, QGraphicsScene, QGraphicsPixmapItem,
    QGraphicsPathItem, QFrame,
)

from .fits_explorer import FITSGraphicsView, FITSExplorer
from .fits_raster import fits_pixel_from_scene, rgb_for_layer
from lsbpipe.survey import maskbits_to_uint


class _LinkedView(FITSGraphicsView):
    zoomChanged = Signal(object)

    def wheelEvent(self, event):
        old = self.transform()
        super().wheelEvent(event)
        if self.transform() != old:
            self.zoomChanged.emit(self)

    def follow(self, source):
        self.setTransform(source.transform())
        self._relative_zoom = source._relative_zoom
        center = source.mapToScene(source.viewport().rect().center())
        self.centerOn(center)


class FourLayerComparison(QWidget):
    """One background conversion, four independently pannable/zoomable views."""

    LAYERS = (
        ("image", "IMAGE · nMgy"),
        ("invvar", "INVVAR · inverse variance"),
        ("nexp", "NEXP · exposures"),
        ("maskbits", "MASKBITS · flags"),
    )

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("fits_four_layer_comparison")
        self._dataset = None
        self._generation = 0
        self._future = None
        self._closing = False
        self._executor = ThreadPoolExecutor(
            max_workers=1, thread_name_prefix="lsbpipe-four-layer-qa")
        self._timer = QTimer(self)
        self._timer.setInterval(30)
        self._timer.timeout.connect(self._poll)
        self.views = {}
        self.scenes = {}
        self.images = {}
        self.outline_items = {}

        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(5)
        tools = QHBoxLayout()
        heading = QLabel("Data quality · one optical band")
        heading.setStyleSheet("font-weight: 600;")
        tools.addWidget(heading)
        tools.addStretch(1)
        self.link_zoom = QCheckBox("Link zoom")
        self.link_zoom.setObjectName("fits_four_link_zoom")
        self.link_zoom.setChecked(True)
        tools.addWidget(self.link_zoom)
        self.fit_all = QPushButton("Fit all")
        self.fit_all.setObjectName("fits_four_fit_all")
        self.fit_all.clicked.connect(self.fit_views)
        tools.addWidget(self.fit_all)
        layout.addLayout(tools)
        guide = QLabel(
            "IMAGE: recorded flux (nanomaggy) · "
            "INVVAR: inverse-variance weight · "
            "NEXP: contributing exposures · "
            "MASKBITS: survey quality flags. "
            "These are four maps of ONE filter, not four optical filters.")
        guide.setWordWrap(True)
        guide.setObjectName("fits_four_layer_explanation")
        layout.addWidget(guide)

        grid = QGridLayout()
        grid.setSpacing(5)
        for i, (kind, title) in enumerate(self.LAYERS):
            cell = QFrame()
            cell.setObjectName(f"fits_four_panel_{kind}")
            cell.setFrameShape(QFrame.StyledPanel)
            box = QVBoxLayout(cell)
            box.setContentsMargins(2, 2, 2, 2)
            box.setSpacing(3)
            caption = QLabel(title)
            caption.setStyleSheet("font-weight: 600;")
            box.addWidget(caption)
            view = _LinkedView()
            view.setObjectName(f"fits_four_view_{kind}")
            scene = QGraphicsScene(view)
            view.setScene(scene)
            image = QGraphicsPixmapItem()
            scene.addItem(image)
            view.pointHovered.connect(
                lambda x, y, k=kind: self._inspect(k, x, y))
            view.zoomChanged.connect(self._sync_zoom)
            box.addWidget(view, 1)
            grid.addWidget(cell, i // 2, i % 2)
            self.views[kind] = view
            self.scenes[kind] = scene
            self.images[kind] = image
            self.outline_items[kind] = []
        grid.setRowStretch(0, 1)
        grid.setRowStretch(1, 1)
        grid.setColumnStretch(0, 1)
        grid.setColumnStretch(1, 1)
        layout.addLayout(grid, 1)

        self.readout = QLabel(
            "Zoom with the wheel · drag to pan · hover for original pixel and ICRS")
        self.readout.setObjectName("fits_four_pixel_readout")
        self.readout.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.readout.setWordWrap(True)
        layout.addWidget(self.readout)
        self.status = QLabel("No native quartet loaded")
        self.status.setObjectName("fits_four_status")
        self.status.setStyleSheet("font-size: 11px;")
        layout.addWidget(self.status)

    @staticmethod
    def _rgb_layers(dataset):
        previews = dataset.get("previews", {})
        for kind, _ in FourLayerComparison.LAYERS:
            if kind in previews:
                rgb = previews[kind]
            else:
                rgb = rgb_for_layer(
                    dataset["data"][kind], kind, dataset["band"])
            yield kind, rgb

    def set_dataset(self, dataset: dict):
        if self._closing:
            raise RuntimeError("Four-layer renderer has been shut down")
        if dataset.get("status") != "DIAGNOSTIC_NOT_ADOPTED":
            raise ValueError("Only diagnostic original FITS accepted")
        shape = dataset["pixel_shape"]
        if (len(shape) != 2 or int(shape[0]) * int(shape[1]) > 6_000_000
                or any(int(d) <= 0 for d in shape)):
            raise ValueError("Native quartet exceeds safe viewer pixel bounds")
        if any(np.shape(dataset["data"][key]) != tuple(shape)
               for key, _ in self.LAYERS):
            raise ValueError("Four original FITS layers must share a grid")
        self.clear_data()
        self._dataset = dataset
        ny, nx = shape
        for kind, _ in self.LAYERS:
            self.scenes[kind].setSceneRect(0, 0, nx, ny)
            self._make_outlines(kind)
        self.status.setText(
            f"{dataset['source']} · {nx} × {ny} original pixels · "
            "diagnostic, not adopted")
        generation = self._generation
        # At most one queued conversion (no per-layer independent workers).
        self._future = self._executor.submit(
            lambda: list(self._rgb_layers(dataset)))
        self._timer.start()

    def _make_outlines(self, kind):
        if self._dataset is None:
            return
        ny = self._dataset["pixel_shape"][0]
        for name, coords in self._dataset["outlines"].items():
            if not coords:
                continue
            path = QPainterPath(
                QPointF(coords[0][0]+.5, ny-1-coords[0][1]+.5))
            for x, y in coords[1:]:
                path.lineTo(x+.5, ny-1-y+.5)
            path.closeSubpath()
            color = QColor("#E45EA8" if name == "square"
                           else "#D6F38E" if name == "target"
                           else "#F6D06F")
            pen = QPen(color, 1.1, Qt.DashLine)
            pen.setCosmetic(True)
            item = QGraphicsPathItem(path)
            item.setPen(pen)
            item.setZValue(2)
            self.scenes[kind].addItem(item)
            self.outline_items[kind].append(item)

    def _poll(self):
        if self._closing:
            self._timer.stop()
            return
        future = self._future
        if future is None or not future.done():
            return
        self._future = None
        self._timer.stop()
        try:
            layers = future.result()
        except Exception as exc:
            self.status.setText(
                f"FITS display conversion failed: {type(exc).__name__}: {exc}")
            return
        if self._dataset is None:
            return
        for kind, rgb in layers:
            self.images[kind].setPixmap(FITSExplorer._pixmap(rgb))
        self.fit_views()
        self.status.setText(
            f"{self._dataset['source']} · interactive native pixels · "
            "display only · not adopted")

    def _sync_zoom(self, source):
        if self.link_zoom.isChecked():
            for view in self.views.values():
                if view is not source:
                    view.follow(source)

    def fit_views(self):
        for view in self.views.values():
            view.fit_image()

    def _inspect(self, kind, x_scene, y_scene):
        dataset = self._dataset
        if dataset is None:
            return
        pixel = fits_pixel_from_scene(
            x_scene, y_scene, dataset["pixel_shape"])
        if pixel is None:
            self.readout.setText("Outside original FITS cutout")
            return
        x, y = pixel
        ra, dec = dataset["wcs"].all_pix2world(float(x), float(y), 0)
        raw = dataset["data"][kind][y, x]
        value = (int(maskbits_to_uint(
            dataset["data"]["maskbits"][y:y+1, x:x+1])[0, 0])
            if kind == "maskbits" else float(raw))
        self.readout.setText(
            f"{kind.upper()} ({x}, {y}) · "
            f"ICRS {ra:.7f}°, {dec:+.7f}° · native value {value:.7g}")

    def clear_data(self):
        self._generation += 1
        self._dataset = None
        self._timer.stop()
        # A still-running conversion may finish, but never publish a result
        # for the old dataset. Do not manipulate it from the GUI thread.
        self._future = None
        for kind, _ in self.LAYERS:
            self.images[kind].setPixmap(QPixmap())
            for item in self.outline_items[kind]:
                self.scenes[kind].removeItem(item)
            self.outline_items[kind].clear()
            self.scenes[kind].setSceneRect(0, 0, 0, 0)
        self.readout.setText("No original FITS data loaded")
        self.status.setText("No native quartet loaded")

    def shutdown_renderer(self):
        if self._closing:
            return
        self._closing = True
        self.clear_data()
        self._executor.shutdown(wait=False, cancel_futures=True)

    def closeEvent(self, event):
        self.shutdown_renderer()
        super().closeEvent(event)
