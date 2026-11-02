"""Read-only, WCS-aware comparison of original IMAGE data across optical bands.

No pixel coaddition/reprojection or synthetic RGB. Native arrays remain
unchanged; individual display stretches are not calibrated colors.
"""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import math
import numpy as np

from PySide6.QtCore import Qt, QTimer, Signal, QPointF, QRectF
from PySide6.QtGui import QPixmap, QColor, QPen, QPainterPath, QTransform
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, QPushButton,
    QComboBox, QCheckBox, QFrame, QGraphicsScene, QGraphicsPixmapItem,
    QGraphicsPathItem,
)

from .fits_explorer import FITSExplorer, FITSGraphicsView
from .fits_raster import rgb_for_layer, fits_pixel_from_scene
from .display_transfer import display_range
from .multiband_preview import (
    MAX_PREVIEW_SIDE, MAX_NATIVE_BYTES, native_dataset_bytes, display_preview,
)

BANDS = ("g", "r", "i", "z")


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
    used_bytes = 0
    for band in BANDS:
        if checkpoint is not None:
            checkpoint()
        # Release/region restrictions are checked BEFORE building paths:
        # DR9 North, for example, cannot provide an i-band quartet.
        if band not in hit.median_exposures:
            missing[band] = "Band unavailable for selected release / region"
            continue
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
        try:
            dataset = prepare_fits_explorer(
                directory=candidate, field=field, hit=hit, band=band,
                target_radius_arcsec=radii[0], sky_inner_arcsec=radii[1],
                sky_outer_arcsec=radii[2])
        except (OSError, ValueError) as exc:
            # One corrupt/incompatible band should be identified and
            # excluded, not presented as valid or block other valid bands.
            missing[band] = (
                "Local quartet failed FITS/WCS validation: "
                f"{type(exc).__name__}: {exc}")
            continue
        required = native_dataset_bytes(dataset)
        if used_bytes + required > MAX_NATIVE_BYTES:
            # Display previews are bounded separately, but native arrays are
            # retained for exact original-pixel hover. Do not overcommit RAM.
            missing[band] = (
                f"Native FITS memory budget exceeded: "
                f"{(used_bytes + required) / (1024**2):.0f} MiB requested, "
                f"{MAX_NATIVE_BYTES // (1024**2)} MiB supported. "
                "Reduce field size or inspect this band individually.")
            continue
        used_bytes += required
        found[band] = dataset
    if not found:
        raise FileNotFoundError(
            "No complete original-FITS quartet for this release/region/brick "
            "in the selected local folder. No remote search was performed.")
    return {"datasets": found, "missing": missing, "status": "DIAGNOSTIC_NOT_ADOPTED"}


class _BandView(FITSGraphicsView):
    navigationChanged = Signal(object)
    viewportResized = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        # Scrollbars obscure the image in narrow three-band layouts.
        # ScrollHandDrag still moves the viewport at user zoom.
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self.viewportResized.emit()

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
        self.cells = {}
        self.preview_info = {}
        self.preview_factors = {}
        self._executor = ThreadPoolExecutor(
            max_workers=1, thread_name_prefix="lsbpipe-band-display")
        self._future = None
        self._generation = 0
        self._closing = False
        self._fit_queued = False
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
        self.stretch.setCurrentIndex(1)
        self.stretch.currentIndexChanged.connect(self._render)
        self.stretch.setToolTip(
            "Display-only percentile limits, independently determined "
            "for each original optical band.")
        tools.addWidget(self.stretch)
        self.scaling = QComboBox()
        self.scaling.setObjectName("fits_multiband_scale")
        self.scaling.addItem("Log10 · positive", "log")
        self.scaling.addItem("Asinh · signed", "asinh")
        self.scaling.addItem("Linear", "linear")
        self.scaling.setToolTip(
            "Log10 shows only strictly positive sky-subtracted IMAGE "
            "pixels; zero and negative samples are black, not missing FITS. "
            "Asinh and linear retain finite signed pixels. Display only.")
        self.scaling.currentIndexChanged.connect(self._render)
        tools.addWidget(self.scaling)
        self.show_missing = QCheckBox("Show missing bands")
        self.show_missing.setObjectName("fits_multiband_show_missing")
        self.show_missing.setChecked(False)
        self.show_missing.toggled.connect(self._layout_panels)
        tools.addWidget(self.show_missing)
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
        self.band_grid = grid
        for n, band in enumerate(BANDS):
            cell = QFrame()
            cell.setObjectName("fits_multiband_panel_" + band)
            cell.setFrameShape(QFrame.StyledPanel)
            column = QVBoxLayout(cell)
            column.setContentsMargins(4, 4, 4, 4)
            column.setSpacing(3)
            column.addWidget(QLabel(band.upper() + " · original IMAGE (nanomaggy)"))
            info = QLabel("Screen overview · preparing")
            info.setObjectName("fits_multiband_overview_" + band)
            info.setStyleSheet("font-size: 11px; color: palette(text)")
            column.addWidget(info)
            view = _BandView()
            view.setObjectName("fits_multiband_view_" + band)
            scene = QGraphicsScene(view)
            image = QGraphicsPixmapItem()
            scene.addItem(image)
            view.setScene(scene)
            view.pointHovered.connect(
                lambda x, y, b=band: self._inspect(b, x, y))
            view.navigationChanged.connect(self._sync)
            view.viewportResized.connect(self._queue_fit_for_resize)
            column.addWidget(view, 1)
            placeholder = QLabel("Not loaded")
            placeholder.setAlignment(Qt.AlignCenter)
            placeholder.setWordWrap(True)
            column.addWidget(placeholder, 1)
            grid.addWidget(cell, n // 2, n % 2)
            self.cells[band] = cell
            self.views[band] = view
            self.scenes[band] = scene
            self.images[band] = image
            self.placeholders[band] = placeholder
            self.preview_info[band] = info
            self.outlines[band] = []
        grid.setRowStretch(0, 1)
        grid.setRowStretch(1, 1)
        grid.setColumnStretch(0, 1)
        grid.setColumnStretch(1, 1)
        layout.addLayout(grid, 1)
        self.readout = QLabel(
            "Hover to compare native flux at the same ICRS position. "
            "Each band's display limits are independently computed; "
            "pixel readouts remain original nanomaggy.")
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
        total_bytes = sum(native_dataset_bytes(d) for d in datasets.values())
        if total_bytes > MAX_NATIVE_BYTES:
            raise ValueError(
                f"Original FITS arrays need {total_bytes / 2**20:.0f} MiB, "
                f"exceeding {MAX_NATIVE_BYTES // 2**20} MiB viewer budget")
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
            self.preview_info[band].setVisible(d is not None)
            if d is None:
                self.placeholders[band].setText(
                    band.upper() + ": " +
                    self.missing.get(band, "Unavailable locally") +
                    "\nNo substitution from another dataset")
                continue
            ny, nx = d["pixel_shape"]
            # The QGraphicsScene retains the exact native FITS cutout.
            # QGraphicsView needs a scrollable *virtual* margin around it:
            # otherwise Qt silently ignores centerOn() for images smaller
            # than the viewport and WCS centers drift between filters.
            self.scenes[band].setSceneRect(0, 0, nx, ny)
            self.views[band].setSceneRect(
                QRectF(-nx / 2., -ny / 2., 2. * nx, 2. * ny))
            self._make_outlines(band)
        self._layout_panels()
        self._render()

    def _layout_panels(self, *_):
        """Use all available width: grz becomes three adjacent FITS panels.

        A missing i band is reported in status instead of permanently
        occupying one quarter of the science workspace. Users can opt to
        show missing tiles explicitly.
        """
        present = [b for b in BANDS if b in self.datasets]
        visible = (list(BANDS) if self.show_missing.isChecked()
                   else present)
        if not visible:
            visible = list(BANDS)
        columns = (3 if len(visible) == 3 else
                   2 if len(visible) in (2, 4) else 1)
        for band in BANDS:
            self.band_grid.removeWidget(self.cells[band])
            self.cells[band].setVisible(band in visible)
        for col in range(4):
            self.band_grid.setColumnStretch(col, 1 if col < columns else 0)
        for row in range(4):
            self.band_grid.setRowStretch(row, 1 if row < math.ceil(len(visible) / columns) else 0)
        for i, band in enumerate(visible):
            self.band_grid.addWidget(
                self.cells[band], i // columns, i % columns)
        self.band_grid.invalidate()
        self._queue_fit_for_resize()

    def _queue_fit_for_resize(self):
        """Defer fitting until Qt has allocated final three-column geometry."""
        if self._closing or not self.datasets or self._fit_queued:
            return
        # An explicit user zoom is retained when resizing the window.
        if any(abs(self.views[b]._relative_zoom - 1.) > 1e-6
               for b in self.datasets):
            return
        self._fit_queued = True
        QTimer.singleShot(0, self._fit_after_layout)

    def _fit_after_layout(self):
        self._fit_queued = False
        if self._closing or not self.datasets:
            return
        if any(abs(self.views[b]._relative_zoom - 1.) > 1e-6
               for b in self.datasets):
            return
        self.fit_views()

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
        scaling = self.scaling.currentData()
        self._future = self._executor.submit(
            lambda: {band: self._prepare_overview(
                d, band, contrast, scaling) for band, d in data.items()})
        self._pending_generation = generation
        self.status.setText("Preparing original IMAGE previews · display only")
        self._timer.start()

    @staticmethod
    def _prepare_overview(dataset, band, contrast, scaling="log"):
        overview, factor = display_preview(dataset["data"]["image"])
        # Quantitative display limits belong to the SCREEN overview, not to
        # science photometry. They are independent for each optical band.
        lo, hi, mode = display_range(
            overview, "image", contrast=contrast, scaling=scaling)
        finite_limits = (math.isfinite(lo) and math.isfinite(hi))
        rgb = rgb_for_layer(
            overview, "image", band, contrast=contrast, scaling=scaling,
            limits=(lo, hi) if finite_limits else None)
        return rgb, factor, lo, hi, mode

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
        self.preview_factors.clear()
        for band, (rgb, factor, lo, hi, mode) in images.items():
            if band not in self.datasets:
                continue
            ny, nx = self.datasets[band]["pixel_shape"]
            preview_y, preview_x = rgb.shape[:2]
            item = self.images[band]
            item.setPixmap(FITSExplorer._pixmap(rgb))
            # The graphics SCENE and WCS remain in ORIGINAL pixel units.
            # Only the display raster is averaged and expanded to cover
            # the native field of view, with no science-array mutation.
            item.setTransform(QTransform.fromScale(
                nx / preview_x, ny / preview_y))
            self.preview_factors[band] = factor
            scale_note = (
                f"{mode} {lo:.3g}…{hi:.3g} nMgy"
                if math.isfinite(lo) and math.isfinite(hi)
                else f"{mode}: no valid positive pixels")
            extent = (
                f"{preview_x}×{preview_y} screen / {nx}×{ny} native "
                f"({factor}×{factor} bin, area mean)"
                if factor > 1 else
                f"{nx}×{ny} native screen")
            self.preview_info[band].setText(
                f"{scale_note} · {extent}")
        self._queue_fit_for_resize()
        present = ", ".join(b.upper() for b in BANDS if b in self.datasets)
        missing = ", ".join(b.upper() for b in BANDS if b not in self.datasets)
        self.status.setText(
            f"Loaded: {present} · unavailable: {missing or 'none'} · "
            f"{self.scaling.currentText()} display (independent limits); "
            "screen area averages where needed; "
            "native FITS readout by WCS · display only · not adopted")

    @staticmethod
    def _scalar(value):
        """Astropy WCS may return scalar-shaped ndarray across versions."""
        array = np.asarray(value)
        if array.size != 1:
            raise ValueError("WCS did not return a scalar coordinate")
        return float(array.reshape(-1)[0])

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
            ra, dec = self._scalar(ra), self._scalar(dec)
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
                tx, ty = self._scalar(tx), self._scalar(ty)
                if not all(math.isfinite(z) for z in (tx, ty)):
                    continue
                scale = (target.viewport().width()
                         * self._pixel_scale(other_d) / angular_width)
                if math.isfinite(scale) and scale > 0:
                    target.resetTransform()
                    target.scale(scale, scale)
                    target._relative_zoom = source._relative_zoom
                    target.centerOn(tx, other_d["pixel_shape"][0] - 1 - ty)
        except (ValueError, ArithmeticError):
            self.status.setText(
                "WCS-linked navigation unavailable; no pixel-coordinate "
                "alignment has been assumed")

    def fit_views(self):
        """Fit complete native cutouts at one angular scale and WCS center.

        The former fit-then-_sync algorithm could magnify target panels
        after fitting, cropping most of a native frame in narrow 3-column
        windows. Compute a common screen angular scale that guarantees all
        rectangular scenes fit, before applying any WCS centering.
        """
        if not self.datasets:
            return
        # Widget geometry may not be stable until the event loop lays out
        # dynamic grz columns. Never force a fit while viewports are zero.
        if any(min(self.views[b].viewport().width(),
                   self.views[b].viewport().height()) < 32
               for b in self.datasets):
            # Await a real viewport resize; never reschedule an idle loop
            # while a hidden workspace has zero-size tiles.
            return
        if not self.link.isChecked() or len(self.datasets) == 1:
            for band in self.datasets:
                self.views[band].fit_image()
            return
        try:
            ref = self.datasets[next(iter(self.datasets))]
            ref_ny, ref_nx = ref["pixel_shape"]
            ra, dec = ref["wcs"].all_pix2world(
                (ref_nx - 1.) / 2., (ref_ny - 1.) / 2., 0)
            ra, dec = self._scalar(ra), self._scalar(dec)
            centers = {}
            min_deg_per_display_pixel = 0.
            for band, d in self.datasets.items():
                ny, nx = d["pixel_shape"]
                xx, yy = d["wcs"].all_world2pix(ra, dec, 0)
                xx, yy = self._scalar(xx), self._scalar(yy)
                if not all(math.isfinite(v) for v in (xx, yy)):
                    raise ValueError("Nonfinite WCS view center")
                scene_y = ny - 1. - yy
                centers[band] = (xx, scene_y)
                view = self.views[band]
                # Account for an off-center WCS match and retain a small
                # safety margin for frame edges / pixel item bounds.
                scene_width = max(nx, 2 * max(xx, nx - xx))
                scene_height = max(ny, 2 * max(scene_y, ny - scene_y))
                deg_per_pixel = self._pixel_scale(d)
                needed = max(
                    scene_width * deg_per_pixel /
                    max(1, view.viewport().width() - 14),
                    scene_height * deg_per_pixel /
                    max(1, view.viewport().height() - 14))
                min_deg_per_display_pixel = max(
                    min_deg_per_display_pixel, needed)
            if not (math.isfinite(min_deg_per_display_pixel)
                    and min_deg_per_display_pixel > 0):
                raise ValueError("Invalid WCS screen scale")
            for band, d in self.datasets.items():
                view = self.views[band]
                scale = self._pixel_scale(d) / min_deg_per_display_pixel
                view.resetTransform()
                view.scale(scale, scale)
                view._relative_zoom = 1.
                view.centerOn(*centers[band])
        except (ValueError, ArithmeticError) as exc:
            # Fail visibly to independent *complete* per-band fitting.
            for band in self.datasets:
                self.views[band].fit_image()
            self.status.setText(
                f"Common WCS fit unavailable ({exc}); "
                "individual complete frames fitted instead.")

    def _inspect(self, band, x_scene, y_scene):
        d = self.datasets.get(band)
        if d is None:
            return
        pixel = fits_pixel_from_scene(x_scene, y_scene, d["pixel_shape"])
        if pixel is None:
            self.readout.setText("Outside the original " + band + "-band cutout")
            return
        ra, dec = d["wcs"].all_pix2world(*map(float, pixel), 0)
        ra, dec = self._scalar(ra), self._scalar(dec)
        readings = []
        for target in BANDS:
            entry = self.datasets.get(target)
            if entry is None:
                readings.append(target + ": unavailable")
                continue
            xx, yy = entry["wcs"].all_world2pix(ra, dec, 0)
            xx, yy = self._scalar(xx), self._scalar(yy)
            if not (math.isfinite(xx) and math.isfinite(yy)):
                readings.append(target + ": outside")
                continue
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
        self.preview_factors.clear()
        for band in BANDS:
            self.images[band].setPixmap(QPixmap())
            self.images[band].setTransform(QTransform())
            self.preview_info[band].setText("Screen overview · preparing")
            for item in self.outlines[band]:
                self.scenes[band].removeItem(item)
            self.outlines[band].clear()
            self.scenes[band].setSceneRect(0, 0, 0, 0)
            self.views[band].setSceneRect(QRectF())
            self.views[band].setVisible(True)
            self.placeholders[band].setVisible(False)
        self._layout_panels()
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
