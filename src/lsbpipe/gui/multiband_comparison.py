"""Read-only, WCS-aware comparison of original IMAGE data across optical bands.

No pixel coaddition/reprojection or synthetic RGB. Native arrays remain
unchanged; individual display stretches are not calibrated colors.
"""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import math
import numpy as np

from PySide6.QtCore import Qt, QTimer, Signal, QPointF
from PySide6.QtGui import QPixmap, QColor, QPen, QPainterPath
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, QPushButton,
    QComboBox, QCheckBox, QFrame, QGraphicsScene, QGraphicsPixmapItem,
    QGraphicsPathItem, QGraphicsScene,
)

from .fits_explorer import FITSExplorer, FITSGraphicsView
from .fits_raster import rgb_for_layer, fits_pixel_from_scene


BANDS = ("g", "r", "i", "z")
MAX_TOTAL_PIXELS = 8_000_000


def load_multiband_originals(*, directory, field, hit, radii, checkpoint=None):
    """Prepare only complete quartets from the *same* release/region/brick.

    Missing bands remain explicitly missing, never substituted with another
    release/region. The full read happens in a GUI-owned background worker.
    """
    from .qa_records import original_quartet
    from .fits_raster import prepare_fits_explorer
    from pathlib import Path

    root = Path(directory).expanduser()
    if not str(directory).strip():
        raise ValueError("Choose a local original-FITS directory first")
    found = {}
    missing = {}
    total = 0
    for band in BANDS:
        if checkpoint is not None:
            checkpoint()
        # Check the exact four-product set, without opening any alternate brick
        # or treating catalog availability as local FITS availability.
        candidate = None
        for path in (root, root / "data"):
            files = original_quartet(path, hit, band)
            if all(f.is_file() for f in files.values()):
                candidate = path
                break
        if candidate is None:
            missing[band] = "Complete IMAGE / INVVAR / NEXP / MASKBITS set not found"
            continue
        # Delegate actual scientific WCS/shape/validity checks to the existing
        # FITS loader, not to this display-only comparison.
        dataset = prepare_fits_explorer(
            directory=candidate, field=field, hit=hit, band=band,
            target_radius_arcsec=radii[0], sky_inner_arcsec=radii[1],
            sky_outer_arcsec=radii[2])
        pixels = int(np.prod(dataset["pixel_shape"]))
        total += pixels
        if total > MAX_TOTAL_PIXELS:
            raise ValueError(
                "Multi-band viewer exceeds 8 million original pixels in total. "
                "Choose a narrower search field; no silent downsampling.")
        found[band] = dataset
    if not found:
        raise FileNotFoundError(
            "No complete original-FITS quartet for this release/region/brick "
            "in the selected local folder. No remote search was performed.")
    return {"datasets": found, "missing": missing, "status": "DIAGNOSTIC_NOT_ADOPTED"}


class _BandView(FITSGraphicsView):
    navigationChanged = Signal(object)

    def wheelEvent(self, event):
        before = self.transform()
        super().wheelEvent(event)
        if self.transform() != before:
            self.navigationChanged.emit(self)

    def mouseReleaseEvent(self, event):
        super().mouseReleaseEvent(event)
        if event.button() == Qt.LeftButton:
            self.navigationChanged.emit(self)


class MultibandComparison(QWidget):
    """Side-by-side original g/r/i/z, linked by sky coordinates, not pixel IDs."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("fits_multiband_comparison")
        self.datasets = {}
        self.missing = {}
        self.views = {}
        self.scenes = {}
        self.images = {}
        self.outlines = {}
        self.placeholders = {}
        self._executor = ThreadPoolExecutor(
            max_workers=1, thread_name_prefix="lsbpipe-band-display")
        self._future = None
        self._generation = 0
        self._closing = False
        self._timer = QTimer(self)
        self._timer.setInterval(30)
        self._timer.timeout.connect(self._poll)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(6)
        tools = QHBoxLayout()
        title = QLabel("Compare original optical bands")
        title.setStyleSheet("font-weight: 600;")
        tools.addWidget(title)
        tools.addStretch(1)
        self.stretch = QComboBox()
        self.stretch.setObjectName("fits_multiband_stretch")
        self.stretch.addItem("Standard", "standard")
        self.stretch.addItem("Faint structures", "faint")
        self.stretch.currentIndexChanged.connect(self._render)
        tools.addWidget(self.stretch)
        self.link = QCheckBox("Link sky position and scale")
        self.link.setObjectName("fits_multiband_link")
        self.link.setChecked(True)
        tools.addWidget(self.link)
        self.fit = QPushButton("Fit all")
        self.fit.setObjectName("fits_multiband_fit")
        self.fit.clicked.connect(self.fit_views)
        tools.addWidget(self.fit)
        layout.addLayout(tools)

        grid = QGridLayout()
        grid.setSpacing(5)
        for n, band in enumerate(BANDS):
            cell = QFrame()
            cell.setObjectName("fits_multiband_panel_" + band)
            cell.setFrameShape(QFrame.StyledPanel)
            column = QVBoxLayout(cell)
            column.setContentsMargins(4, 4, 4, 4)
            column.setSpacing(3)
            column.addWidget(QLabel(band.upper() + " · original IMAGE (nanomaggy)"))
            view = _BandView()
            view.setObjectName("fits_multiband_view_" + band)
            scene = QGraphicsScene(view)
            image = QGraphicsPixmapItem()
            scene.addItem(image)
            view.setScene(scene)
            view.pointHovered.connect(
                lambda x, y, b=band: self._inspect(b, x, y))
            view.navigationChanged.connect(self._sync)
            column.addWidget(view, 1)
            placeholder = QLabel("Not loaded")
            placeholder.setAlignment(Qt.AlignCenter)
            placeholder.setWordWrap(True)
            column.addWidget(placeholder)
            grid.addWidget(cell, n // 2, n % 2)
            self.views[band] = view
            self.scenes[band] = scene
            self.images[band] = image
            self.placeholders[band] = placeholder
            self.outlines[band] = []
        grid.setRowStretch(0, 1)
        grid.setRowStretch(1, 1)
        grid.setColumnStretch(0, 1)
        grid.setColumnStretch(1, 1)
        layout.addLayout(grid, 1)
        self.readout = QLabel(
            "Hover to compare native flux at the same ICRS position. "
            "Different band stretches are independent.")
        self.readout.setObjectName("fits_multiband_readout")
        self.readout.setWordWrap(True)
        self.readout.setTextInteractionFlags(Qt.TextSelectableByMouse)
        layout.addWidget(self.readout)
        self.status = QLabel(
            "Display only · no PSF matching, reprojection, or science adoption")
        self.status.setObjectName("fits_multiband_status")
        self.status.setWordWrap(True)
        layout.addWidget(self.status)

    def set_data(self, result):
        if self._closing:
            raise RuntimeError("Multiband display already shut down")
        if result.get("status") != "DIAGNOSTIC_NOT_ADOPTED":
            raise ValueError("Only read-only diagnostic FITS are permitted")
        datasets = result.get("datasets", {})
        if not datasets or set(datasets) - set(BANDS):
            raise ValueError("No supported original optical bands")
        if sum(int(np.prod(d["pixel_shape"])) for d in datasets.values()
               ) > MAX_TOTAL_PIXELS:
            raise ValueError("Multi-band viewer memory limit exceeded")
        identity = None
        for band, d in datasets.items():
            if d.get("status") != "DIAGNOSTIC_NOT_ADOPTED" or d["band"] != band:
                raise ValueError("Incorrect band or diagnostic status")
            if any(np.shape(v) != tuple(d["pixel_shape"])
                   for v in d["data"].values()):
                raise ValueError("Unregistered original FITS layers")
            key = (d["field"].ra_deg, d["field"].dec_deg,
                   d["field"].width_arcmin, d["source"].rsplit(" · ", 1)[0])
            if identity is not None and identity != key:
                raise ValueError("Cannot compare different objects or bricks")
            identity = key
        self.clear_data()
        self.datasets = dict(datasets)
        self.missing = dict(result.get("missing", {}))
        for band in BANDS:
            d = self.datasets.get(band)
            self.views[band].setVisible(d is not None)
            self.placeholders[band].setVisible(d is None)
            if d is None:
                self.placeholders[band].setText(
                    band.upper() + " unavailable locally · no substitution")
                continue
            ny, nx = d["pixel_shape"]
            self.scenes[band].setSceneRect(0, 0, nx, ny)
            self._make_outlines(band)
        self._render()

    def _make_outlines(self, band):
        d = self.datasets[band]
        ny = d["pixel_shape"][0]
        for name in ("square", "target"):
            points = d["outlines"].get(name, ())
            if not points:
                continue
            path = QPainterPath(
                QPointF(points[0][0] + .5, ny - 1 - points[0][1] + .5))
            for x, y in points[1:]:
                path.lineTo(x + .5, ny - 1 - y + .5)
            item = QGraphicsPathItem(path)
            pen = QPen(QColor("#E45EA8" if name == "square" else "#D6F38E"),
                       1.0, Qt.DashLine)
            pen.setCosmetic(True)
            item.setPen(pen)
            item.setZValue(2)
            self.scenes[band].addItem(item)
            self.outlines[band].append(item)

    def _render(self, *_):
        if self._closing or not self.datasets:
            return
        self._generation += 1
        generation = self._generation
        # A single bounded worker converts grayscale previews off the UI
        # thread. Source FITS arrays are read-only and are never normalized
        # in place or reprojected.
        data = dict(self.datasets)
        contrast = self.stretch.currentData()
        self._future = self._executor.submit(
            lambda: {band: rgb_for_layer(
                d["data"]["image"], "image", band, contrast=contrast)
                for band, d in data.items()})
        self._pending_generation = generation
        self.status.setText("Preparing original IMAGE previews · display only")
        self._timer.start()

    def _poll(self):
        future = self._future
        if self._closing or future is None or not future.done():
            return
        self._future = None
        self._timer.stop()
        try:
            images = future.result()
        except Exception as exc:
            self.status.setText(f"Display conversion failed: {exc}")
            return
        if not self.datasets:
            return
        for band, arr in images.items():
            if band in self.datasets:
                self.images[band].setPixmap(FITSExplorer._pixmap(arr))
        self.fit_views()
        present = ", ".join(b.upper() for b in BANDS if b in self.datasets)
        missing = ", ".join(b.upper() for b in BANDS if b not in self.datasets)
        self.status.setText(
            f"Loaded: {present} · unavailable: {missing or 'none'} · "
            "individual display stretches · original FITS only · not adopted")

    @staticmethod
    def _pixel_scale(d):
        # WCS pixel scale in degrees per pixel, used for equal angular
        # viewport width; pixel positions are never assumed to coincide.
        from astropy.wcs.utils import proj_plane_pixel_scales
        scales = proj_plane_pixel_scales(d["wcs"].celestial)
        scale = float(np.mean(np.abs(scales[:2])))
        if not math.isfinite(scale) or scale <= 0:
            raise ValueError("Cannot establish WCS angular pixel scale")
        return scale

    def _sync(self, source):
        if not self.link.isChecked() or len(self.datasets) < 2:
            return
        band = next((b for b, v in self.views.items() if v is source), None)
        if band not in self.datasets:
            return
        d = self.datasets[band]
        center = source.mapToScene(source.viewport().rect().center())
        ny = d["pixel_shape"][0]
        x, y = center.x(), ny - 1 - center.y()
        try:
            ra, dec = d["wcs"].all_pix2world(x, y, 0)
            if not all(math.isfinite(z) for z in (ra, dec)):
                return
            angular_width = (
                source.viewport().width() / abs(source.transform().m11())
                * self._pixel_scale(d))
            for other_band, other_d in self.datasets.items():
                if other_band == band:
                    continue
                target = self.views[other_band]
                tx, ty = other_d["wcs"].all_world2pix(ra, dec, 0)
                if not all(math.isfinite(z) for z in (tx, ty)):
                    continue
                scale = (target.viewport().width()
                         * self._pixel_scale(other_d) / angular_width)
                if math.isfinite(scale) and scale > 0:
                    target.resetTransform()
                    target.scale(scale, scale)
                    target._relative_zoom = 2.  # Preserve focus on resize
                    target.centerOn(tx, other_d["pixel_shape"][0] - 1 - ty)
        except (ValueError, ArithmeticError):
            self.status.setText(
                "WCS-linked navigation unavailable; no pixel-coordinate "
                "alignment has been assumed")

    def fit_views(self):
        for band in self.datasets:
            self.views[band].fit_image()
        if self.datasets and self.link.isChecked():
            self._sync(self.views[next(iter(self.datasets))])

    def _inspect(self, band, x_scene, y_scene):
        d = self.datasets.get(band)
        if d is None:
            return
        pixel = fits_pixel_from_scene(x_scene, y_scene, d["pixel_shape"])
        if pixel is None:
            self.readout.setText("Outside the original " + band + "-band cutout")
            return
        ra, dec = d["wcs"].all_pix2world(*map(float, pixel), 0)
        readings = []
        for target in BANDS:
            entry = self.datasets.get(target)
            if entry is None:
                readings.append(target + ": unavailable")
                continue
            xx, yy = entry["wcs"].all_world2pix(ra, dec, 0)
            xi, yi = int(round(xx)), int(round(yy))
            ny, nx = entry["pixel_shape"]
            if not (0 <= xi < nx and 0 <= yi < ny):
                readings.append(target + ": outside")
                continue
            flux = float(entry["data"]["image"][yi, xi])
            invvar = float(entry["data"]["invvar"][yi, xi])
            nexp = float(entry["data"]["nexp"][yi, xi])
            if not math.isfinite(flux) or not math.isfinite(invvar) or invvar <= 0 or nexp <= 0:
                readings.append(target + ": invalid sample")
            else:
                readings.append(f"{target}: {flux:.4g} nMgy")
        self.readout.setText(
            f"ICRS {ra:.7f}°, {dec:+.7f}°  ·  "
            + "  |  ".join(readings)
            + "  ·  nearest native pixels, not PSF-matched photometry")

    def clear_data(self):
        self._generation += 1
        self._future = None
        self._timer.stop()
        self.datasets = {}
        self.missing = {}
        for band in BANDS:
            self.images[band].setPixmap(QPixmap())
            for item in self.outlines[band]:
                self.scenes[band].removeItem(item)
            self.outlines[band].clear()
            self.scenes[band].setSceneRect(0, 0, 0, 0)
            self.views[band].setVisible(True)
            self.placeholders[band].setVisible(False)
        self.status.setText("No original bands loaded")

    def shutdown_renderer(self):
        if self._closing:
            return
        self._closing = True
        self.clear_data()
        self._executor.shutdown(wait=False, cancel_futures=True)

    def closeEvent(self, event):
        self.shutdown_renderer()
        super().closeEvent(event)
