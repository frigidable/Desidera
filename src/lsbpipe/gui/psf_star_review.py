"""Interactive in-memory PSF training-star review for the existing GUI.

Cards show the exact local-sky-subtracted native cutout, the signed
candidate ePSF evaluated on the same pixel grid, and a *square* radial
diagnostic plot (observed vs model). A reviewer can exclude/restore only
training stars; independent holdouts are not mutable through this UI.

Images are diagnostic previews built from arrays already in worker memory.
No source FITS are modified, no manual decision is silently applied,
and no rejected model is adopted as a scientific PSF.
"""
from __future__ import annotations

import math
import numpy as np

from PySide6.QtCore import Qt, QRectF, QPointF
from PySide6.QtGui import (
    QColor, QImage, QPainter, QPen, QPixmap,
)
from PySide6.QtWidgets import (
    QCheckBox, QDialog, QDialogButtonBox, QFrame, QGridLayout,
    QHBoxLayout, QLabel, QPushButton, QScrollArea, QVBoxLayout, QWidget,
)


IMAGE_SIDE = 126


def _square_image(values, *, relative_peak, valid=None):
    """Signed asinh contrast; equal angular/pixel scale for image and model."""
    data = np.asarray(values, dtype=float)
    if data.ndim != 2 or min(data.shape) < 9 or not np.isfinite(data).all():
        raise ValueError("Star review image requires a finite 2D cutout")
    peak = max(float(relative_peak), 1e-10)
    abs_value = np.clip(np.abs(data)/peak, 0, 100.)
    transfer = np.arcsinh(15. * abs_value) / math.asinh(15.)
    scaled = np.clip(np.rint(205*transfer), 0, 205).astype(np.int16)
    rgb = np.zeros((*data.shape, 3), dtype=np.uint8)
    rgb[:] = (30, 32, 38)
    pos = data >= 0
    neg = ~pos
    rgb[pos, 0] = 34
    rgb[pos, 1] = 50 + scaled[pos] * 3//4
    rgb[pos, 2] = 50 + scaled[pos]
    rgb[neg, 0] = 50 + scaled[neg]
    rgb[neg, 1] = 34
    rgb[neg, 2] = 40
    if valid is not None:
        mask = np.asarray(valid, dtype=bool)
        if mask.shape != data.shape:
            raise ValueError("Review mask and original cutout sizes disagree")
        rgb[~mask] = (95, 65, 12)
    ny, nx = data.shape
    return QPixmap.fromImage(
        QImage(rgb.tobytes(), nx, ny, nx*3, QImage.Format_RGB888).copy()
    ).scaled(IMAGE_SIDE, IMAGE_SIDE, Qt.KeepAspectRatio, Qt.FastTransformation)


def _radial_profile(values, center, *, peak):
    data = np.asarray(values, dtype=float)
    cy, cx = float(center[1]), float(center[0])
    y, x = np.indices(data.shape)
    rr = np.hypot(x-cx, y-cy)
    max_r = min(data.shape)//2.
    radii = np.arange(0., max_r, 1.)
    med = np.full(len(radii), np.nan)
    for i, r in enumerate(radii):
        sample = (rr >= r) & (rr < r+1) & np.isfinite(data)
        if np.any(sample):
            med[i] = np.median(data[sample])/max(float(peak), 1e-10)
    return radii, med


def _radial_plot(observed, modeled, center, *, peak):
    """Square QPainter profile; same cutout radius and both flux scales."""
    from PySide6.QtCore import QPointF
    image = QImage(IMAGE_SIDE, IMAGE_SIDE, QImage.Format_RGB32)
    image.fill(QColor(25, 28, 33))
    painter = QPainter(image)
    painter.setRenderHint(QPainter.Antialiasing, True)
    left, right, top, bottom = 17., 5., 8., 16.
    axis_w = IMAGE_SIDE-left-right
    axis_h = IMAGE_SIDE-top-bottom
    def to_screen(x, y, xmax):
        # Signed flux normalized by measured star peak, not total PSF flux.
        return QPointF(
            left + (x/max(xmax, 1.))*axis_w,
            top + ((1.08 - np.clip(y, -.16, 1.08))/1.24)*axis_h)
    painter.setPen(QPen(QColor(95, 99, 106), 1))
    painter.drawLine(QPointF(left, top), QPointF(left, top+axis_h))
    painter.drawLine(QPointF(left, top+axis_h), QPointF(left+axis_w, top+axis_h))
    painter.setPen(QPen(QColor(85, 85, 85), 1, Qt.DashLine))
    y0 = to_screen(0., 0., 1).y()
    painter.drawLine(QPointF(left, y0), QPointF(left+axis_w, y0))
    series = (
        (observed, QColor(88, 210, 233)),
        (modeled, QColor(244, 155, 93)),
    )
    for values, color in series:
        xx, yy = _radial_profile(values, center, peak=peak)
        painter.setPen(QPen(color, 1.7))
        for index in range(1, len(xx)):
            if np.isfinite(yy[index-1]) and np.isfinite(yy[index]):
                painter.drawLine(
                    to_screen(xx[index-1], yy[index-1], xx[-1]+1),
                    to_screen(xx[index], yy[index], xx[-1]+1))
    painter.setPen(QColor(174, 178, 187))
    painter.drawText(QPointF(2, 12), "1")
    painter.drawText(QPointF(3, y0+3), "0")
    painter.drawText(QPointF(IMAGE_SIDE-25, IMAGE_SIDE-3), "r/pix")
    painter.end()
    return QPixmap.fromImage(image)


class PSFStarReviewDialog(QDialog):
    """All training cutouts, side-by-side native PSF and radial QA."""
    def __init__(self, snapshot: dict, *, parent=None):
        super().__init__(parent)
        if snapshot.get("schema") != "desidera.psf0-star-review-in-memory/1":
            raise ValueError("Unsupported or missing PSF training-star review data")
        self.snapshot = snapshot
        self._controls = {}
        self.setWindowTitle("Desidera · PSF training-star mosaic · manual review")
        self.setObjectName("psf0_training_star_mosaic")
        self.resize(1160, 770)
        root = QVBoxLayout(self)
        heading = QLabel(
            "Training-star review · cyan positive, red negative, ochre masked. "
            "Each card: original sky-subtracted cutout / aligned candidate PSF / "
            "observed and model radial profiles (all square).")
        heading.setWordWrap(True)
        root.addWidget(heading)
        notice = QLabel(
            "Exclusions affect ONLY the ePSF training sample. Independent "
            "holdout stars stay untouched. Changes apply only after clicking "
            "Apply exclusions & rebuild; current model remains non-adopted.")
        notice.setWordWrap(True)
        root.addWidget(notice)
        self.count_label = QLabel()
        self.count_label.setObjectName("psf0_review_counts")
        root.addWidget(self.count_label)
        scroll = QScrollArea()
        scroll.setObjectName("psf0_star_mosaic_scroll")
        scroll.setWidgetResizable(True)
        content = QWidget()
        grid = QGridLayout(content)
        grid.setContentsMargins(10, 10, 10, 10)
        grid.setHorizontalSpacing(9)
        grid.setVerticalSpacing(11)
        included = 0
        for index, star in enumerate(snapshot["stars"]):
            role = star["role"]
            if role not in ("train", "excluded"):
                continue  # holdout stars are never a selectable review control
            included += 1
            sid = star["star_id"]
            card = QFrame(content)
            card.setFrameShape(QFrame.StyledPanel)
            card.setObjectName("psf0_star_card_" + sid)
            box = QVBoxLayout(card)
            x, y = star["x_pix"], star["y_pix"]
            label = QLabel(f"{sid} · ({x:.1f}, {y:.1f}) pix")
            label.setToolTip(
                f"ICRS {star['ra_deg']:.7f}°, {star['dec_deg']:.7f}°"
                " · source pixels are not changed")
            box.addWidget(label)
            plots = QHBoxLayout()
            observed = np.asarray(star["cutout"], dtype=float)
            model = np.asarray(star["model"], dtype=float)
            if observed.shape != model.shape:
                raise ValueError("Star cutout and ePSF model grid differ")
            peak = star["peak_nMgy"]
            image_pixmaps = (
                ("Cutout", _square_image(
                    observed, relative_peak=peak, valid=star["valid"])),
                ("PSF model", _square_image(
                    model, relative_peak=peak, valid=star["valid"])),
                ("Radial PSF", _radial_plot(
                    observed, model, star["center_in_cutout"], peak=peak)),
            )
            for title, pixmap in image_pixmaps:
                panel = QVBoxLayout()
                legend = QLabel(title)
                legend.setAlignment(Qt.AlignHCenter)
                panel.addWidget(legend)
                graphic = QLabel()
                graphic.setObjectName("psf0_review_" + title.lower().replace(" ", "_") + "_" + sid)
                graphic.setPixmap(pixmap)
                graphic.setFixedSize(IMAGE_SIDE, IMAGE_SIDE)
                graphic.setToolTip(
                    "Original pixel grid, no interpolation of observed cutout. "
                    "PSF alignment/model only for review, not science-adopted.")
                panel.addWidget(graphic)
                plots.addLayout(panel)
            box.addLayout(plots)
            check = QCheckBox("Exclude from PSF training")
            check.setObjectName("psf0_exclude_" + sid)
            check.setChecked(role == "excluded")
            check.toggled.connect(self._update_counts)
            self._controls[sid] = check
            box.addWidget(check)
            footer = QLabel(
                f"Sky {star['sky_nMgy']:.3g} nMgy/pix  · "
                f"gradient/RMS {star['gradient_sigma_span']:.2f}")
            footer.setWordWrap(True)
            box.addWidget(footer)
            grid.addWidget(card, index // 2, index % 2)
        if included < 1:
            raise ValueError("No original training stars available for review")
        scroll.setWidget(content)
        root.addWidget(scroll, stretch=1)
        footer = QHBoxLayout()
        reset = QPushButton("Restore all training stars")
        reset.setObjectName("psf0_restore_training")
        reset.clicked.connect(self._restore_all)
        footer.addWidget(reset)
        footer.addStretch(1)
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.clicked.connect(self.reject)
        footer.addWidget(self.cancel_button)
        self.apply_button = QPushButton("Apply exclusions & rebuild")
        self.apply_button.setObjectName("psf0_apply_exclusions")
        self.apply_button.clicked.connect(self.accept)
        footer.addWidget(self.apply_button)
        root.addLayout(footer)
        self._update_counts()

    def _restore_all(self):
        for check in self._controls.values():
            check.setChecked(False)

    def _update_counts(self, *_):
        count = len(self._controls)
        excluded = len(self.excluded_ids())
        self.count_label.setText(
            f"Training candidate cards: {count} · excluded manually: "
            f"{excluded} · to rebuild: {count-excluded}. "
            "No holdout object is selectable.")

    def excluded_ids(self):
        return tuple(sorted(k for k, ctrl in self._controls.items()
                            if ctrl.isChecked()))

    def changed(self):
        return self.excluded_ids() != tuple(
            self.snapshot.get("excluded_training_star_ids", ()))
