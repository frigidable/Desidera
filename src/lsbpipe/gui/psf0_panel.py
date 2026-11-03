"""Interactive PSF-0 control embedded in the existing Galaxy Project GUI.

The researcher selects an already-present band, performs the computation
in the existing Qt worker, inspects its non-adopted diagnostics, and may
save explicit immutable derivatives. No filesystem, FITS or network side
effects occur in this widget itself.
"""
from __future__ import annotations

import numpy as np

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QImage, QPixmap
from PySide6.QtWidgets import (
    QAbstractItemView, QComboBox, QDoubleSpinBox, QFormLayout, QGroupBox,
    QHBoxLayout, QHeaderView, QLabel, QPlainTextEdit, QPushButton,
    QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget,
)

from lsbpipe.photometry.psf_gui_backend import inspect_psf0_bands


class PSF0Panel(QGroupBox):
    """Project-bound, band-specific **diagnostic** ePSF workbench."""

    runRequested = Signal(object)
    saveRequested = Signal(object)

    def __init__(self, parent=None):
        super().__init__("04 · Empirical PSF · first pilot Leo P", parent)
        self.setObjectName("photometry_psf0_panel")
        self._project_root = None
        self._snapshot = None
        self._preview = None
        self._busy = False

        layout = QVBoxLayout(self)
        layout.setSpacing(8)
        explanation = QLabel(
            "Uses the selected galaxy project's original FITS, not a web preview. "
            "Run only the available band; preview is read-only. "
            "A model is experimental until stellar, wing and spatial QA pass.")
        explanation.setWordWrap(True)
        layout.addWidget(explanation)

        form = QFormLayout()
        self.band = QComboBox()
        self.band.setObjectName("psf0_band")
        self.band.setAccessibleName("Native optical band for provisional PSF")
        self.band.currentIndexChanged.connect(self._band_changed)
        form.addRow("Native band", self.band)
        self.fwhm = QDoubleSpinBox()
        self.fwhm.setObjectName("psf0_fwhm_guess")
        self.fwhm.setRange(2.5, 20.)
        self.fwhm.setDecimals(2)
        self.fwhm.setSingleStep(.25)
        self.fwhm.setValue(4.8)
        self.fwhm.setSuffix(" px")
        self.fwhm.setToolTip(
            "DAOStarFinder detection-width GUESS, not a measured FWHM. "
            "Changing it does not repair a scientifically rejected ePSF.")
        form.addRow("Detection FWHM guess", self.fwhm)
        self.protection = QDoubleSpinBox()
        self.protection.setObjectName("psf0_target_protection")
        self.protection.setRange(1., 1800.)
        self.protection.setDecimals(1)
        self.protection.setValue(90.)
        self.fwhm.valueChanged.connect(self._parameters_changed)
        self.protection.valueChanged.connect(self._parameters_changed)
        self.protection.setSuffix(" arcsec")
        self.protection.setToolTip(
            "Preliminary exclusion circle around Leo P. Must be reviewed "
            "against diffuse galaxy light before science-quality PSF.")
        form.addRow("Protect target radius", self.protection)
        layout.addLayout(form)

        actions = QHBoxLayout()
        self.run_btn = QPushButton("Build PSF preview")
        self.run_btn.setObjectName("psf0_run_button")
        self.run_btn.setToolTip(
            "Runs real FITS validation, star selection and ePSF calculation "
            "in the existing background worker; does not write a PSF.")
        self.run_btn.clicked.connect(self._request_run)
        actions.addWidget(self.run_btn)
        self.save_btn = QPushButton("Save diagnostic FITS + report")
        self.save_btn.setObjectName("psf0_save_button")
        self.save_btn.setEnabled(False)
        self.save_btn.setToolTip(
            "Only after a successful preview: verify original FITS hashes "
            "again and save an immutable non-adopted PSF derivative.")
        self.save_btn.clicked.connect(self._request_save)
        actions.addWidget(self.save_btn)
        actions.addStretch(1)
        layout.addLayout(actions)

        self.status = QLabel(
            "Open saved Leo P project; source bands have not been checked.")
        self.status.setObjectName("psf0_status")
        self.status.setWordWrap(True)
        layout.addWidget(self.status)

        self.summary = QLabel("No ePSF result. No stellar sample has been verified.")
        self.summary.setObjectName("psf0_summary")
        self.summary.setWordWrap(True)
        self.summary.setTextInteractionFlags(Qt.TextSelectableByMouse)
        layout.addWidget(self.summary)

        details = QHBoxLayout()
        self.kernel_view = QLabel("No empirical PSF image")
        self.kernel_view.setObjectName("psf0_kernel_preview")
        self.kernel_view.setMinimumSize(160, 160)
        self.kernel_view.setMaximumSize(245, 245)
        self.kernel_view.setAlignment(Qt.AlignCenter)
        self.kernel_view.setWordWrap(True)
        self.kernel_view.setToolTip(
            "Contrast-stretched diagnostic discrete kernel, NOT measured "
            "seeing or a calibrated wide-angle PSF.")
        details.addWidget(self.kernel_view, 0)
        self.stars = QTableWidget(0, 5)
        self.stars.setObjectName("psf0_star_audit")
        self.stars.setHorizontalHeaderLabels(
            ["Role", "X (pix)", "Y (pix)", "Cutout sky", "Gradient/RMS"])
        self.stars.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.stars.verticalHeader().setVisible(False)
        self.stars.setSelectionMode(QAbstractItemView.NoSelection)
        self.stars.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.stars.setMinimumHeight(180)
        self.stars.setToolTip(
            "Training and held-out PSF stars. Coordinates are native FITS "
            "pixel indices (0-based); fields are not Gaia-confirmed.")
        details.addWidget(self.stars, 1)
        layout.addLayout(details)

        self.qa = QPlainTextEdit()
        self.qa.setObjectName("psf0_qa_details")
        self.qa.setReadOnly(True)
        self.qa.setMinimumHeight(110)
        self.qa.setMaximumHeight(210)
        self.qa.setPlaceholderText(
            "PSF quality flags and source limitations appear here.")
        layout.addWidget(self.qa)

    @property
    def active_preview(self):
        return self._preview

    def bind_project(self, project_root):
        self._project_root = project_root
        self.refresh_sources()

    def refresh_sources(self):
        self._preview = None
        self.save_btn.setEnabled(False)
        self.stars.setRowCount(0)
        self.kernel_view.setPixmap(QPixmap())
        self.kernel_view.setText("No empirical PSF image")
        self.qa.clear()
        self.summary.setText("No ePSF measurement has been made or adopted.")
        current_band = self.band.currentData()
        self.band.blockSignals(True)
        self.band.clear()
        self._snapshot = None
        try:
            if self._project_root is None:
                raise ValueError("No Galaxy Project is open")
            snapshot = inspect_psf0_bands(self._project_root)
        except (OSError, ValueError, TypeError, KeyError) as exc:
            self.band.addItem("Unavailable", None)
            self.status.setText(f"PSF-0 unavailable: {exc}")
        else:
            self._snapshot = snapshot
            available = [e for e in snapshot["bands"]
                         if e["paths_present_unverified"]]
            for entry in snapshot["bands"]:
                label = (f"{entry['band']} · 4/4 files found (unverified)"
                         if entry["paths_present_unverified"]
                         else f"{entry['band']} · {entry['present_count']}/4 files")
                self.band.addItem(label, entry["band"])
                model = self.band.model()
                item = model.item(self.band.count()-1)
                if item is not None:
                    item.setEnabled(entry["paths_present_unverified"])
            choices = [e["band"] for e in available]
            chosen = (current_band if current_band in choices else
                      snapshot["dataset"][3] if snapshot["dataset"][3] in choices
                      else choices[0] if choices else None)
            if chosen is not None:
                self.band.setCurrentIndex(self.band.findData(chosen))
            self.protection.setValue(
                min(1800., snapshot["target_protection_default_arcsec"]))
            self.status.setText(
                "Original FITS paths found, but not yet validated. "
                "Click Build PSF preview to verify FITS/WCS/units and run.")
            if not available:
                self.status.setText(
                    "No complete local native FITS quartet for PSF-0. "
                    "Use Manage sources; no automatic download.")
        self.band.blockSignals(False)
        self._update_buttons()

    def _parameters_changed(self, *_):
        if self._preview is not None:
            self._band_changed()

    def _band_changed(self, *_):
        # Changing filter invalidates any earlier successful kernel.
        self._preview = None
        self.save_btn.setEnabled(False)
        self.stars.setRowCount(0)
        self.kernel_view.setPixmap(QPixmap())
        self.kernel_view.setText("No empirical PSF image")
        self.summary.setText("No PSF preview for this band.")
        self.qa.clear()
        self._update_buttons()

    def _update_buttons(self):
        selected = self.band.currentData()
        eligible = bool(self._snapshot) and any(
            e["band"] == selected and e["paths_present_unverified"]
            for e in self._snapshot["bands"])
        self.run_btn.setEnabled(eligible and not self._busy)
        self.save_btn.setEnabled(self._preview is not None and not self._busy
                                 and self._preview["result"].lineage_id is not None)
        self.band.setEnabled(not self._busy)
        self.fwhm.setEnabled(not self._busy)
        self.protection.setEnabled(not self._busy)

    def set_busy(self, state: bool, message: str = ""):
        self._busy = bool(state)
        if message:
            self.status.setText(message)
        self._update_buttons()

    def _request_run(self):
        if self._project_root is None or self._busy or self.band.currentData() is None:
            return
        self._preview = None
        self.save_btn.setEnabled(False)
        self.stars.setRowCount(0)
        self.kernel_view.setPixmap(QPixmap())
        self.kernel_view.setText("PSF calculation running…")
        self.summary.setText("No accepted PSF; calculating provisional model.")
        self.qa.clear()
        self.runRequested.emit({
            "project_root": str(self._project_root),
            "band": self.band.currentData(),
            "fwhm_pix": self.fwhm.value(),
            "target_protection_arcsec": self.protection.value(),
        })

    def _request_save(self):
        if self._busy or self._preview is None:
            return
        self.saveRequested.emit(self._preview)

    def show_failure(self, error: str):
        self._preview = None
        self.save_btn.setEnabled(False)
        self.kernel_view.setPixmap(QPixmap())
        self.kernel_view.setText("Rejected ePSF · no kernel adopted")
        self.stars.setRowCount(0)
        self.status.setText("PSF-0 failed scientific preparation · NOT ADOPTED")
        if "negative EPSF ringing" in error:
            explanation = (
                "This is a scientifically blocking ePSF-quality failure. "
                "Inspect isolated-star selection, background annuli, "
                "subpixel alignment and the distribution of negative model values. "
                "Do NOT solve it by clipping/relaxing thresholds alone.\n\n")
        else:
            explanation = ""
        self.qa.setPlainText(explanation + error)
        self.summary.setText("No PSF saved or adopted. Original FITS unchanged.")
        self._update_buttons()

    def show_quality_rejection(self, payload: dict):
        """Rejected ePSF: inspect stars and failure metrics without a kernel."""
        self._preview = None
        self.kernel_view.setPixmap(QPixmap())
        self.kernel_view.setText("REJECTED ePSF · no kernel accepted")
        self.status.setText(
            "PSF-0 QUALITY REJECTED · no PSF model adopted or saved")
        diagnostic = payload.get("diagnostics", {})
        negative = diagnostic.get("negative_absolute_mass_fraction")
        if negative is not None:
            self.summary.setText(
                f"Negative absolute ePSF mass: {negative:.2%} "
                f"(experimental limit {diagnostic.get('blocking_threshold_fraction', .025):.2%}); "
                f"training stars {diagnostic.get('n_train', 0)}, "
                f"holdout stars {diagnostic.get('n_holdout', 0)}. "
                "Review the source selection and local backgrounds.")
        else:
            self.summary.setText("Provisional PSF rejected before adoption.")
        rows = [(item, "train") for item in diagnostic.get("training", [])]
        rows += [(item, "holdout") for item in diagnostic.get("holdout", [])]
        self.stars.setRowCount(len(rows))
        for index, (star, role) in enumerate(rows):
            for col, label in enumerate((
                    role, f"{star['x_pix']:.2f}", f"{star['y_pix']:.2f}",
                    f"{star['cutout_sky_nMgy']:.4g}",
                    f"{star['gradient_sigma_span']:.2f}")):
                self.stars.setItem(index, col, QTableWidgetItem(label))
        self.qa.setPlainText(
            payload.get("reason", "PSF reconstruction quality failure")
            + "\n\nRejected candidates: "
            + str(diagnostic.get("rejected_reasons", {}))
            + "\n\n"
            + "\n".join(diagnostic.get("notes", []))
            + "\n\nNo invalid PSF kernel or photometric measurement was saved.")
        self._update_buttons()

    def show_result(self, preview: dict):
        if (self._project_root is None or
                str(self._project_root) != str(preview["project_root"])):
            self.show_failure("Result project is different from the open Galaxy Project")
            return
        result = preview["result"]
        self._preview = preview
        report = result.report()
        self.status.setText(
            f"{result.dataset_id} · EXPERIMENTAL_PSF_NOT_ADOPTED")
        self.summary.setText(
            f"Training stars: {len(result.training)} · Independent holdout: "
            f"{len(result.holdout)} · Detected candidates: {result.detection_count} · "
            f"Median heldout core |residual|/peak: "
            f"{result.metrics['holdout_median_abs_core_residual_over_peak']:.3g}. "
            "This is NOT a calibrated science PSF or FWHM.")
        rows = result.training + result.holdout
        self.stars.setRowCount(len(rows))
        for i, star in enumerate(rows):
            for j, text in enumerate((
                    star.role, f"{star.x_pix:.2f}", f"{star.y_pix:.2f}",
                    f"{star.cutout_sky_nMgy:.4g}",
                    f"{star.gradient_sigma_span:.2f}")):
                self.stars.setItem(i, j, QTableWidgetItem(text))
        self.qa.setPlainText(
            "Scientific limits:\n" + "\n".join(f"• {f}" for f in result.flags)
            + "\n\nRejected star counts: " + str(result.rejected_reasons)
            + "\nSource SHA-256 lineage: "
            + (result.lineage_id or "not verified"))
        kernel = result.kernel
        mapped = np.log1p(np.maximum(kernel, 0.) * 1000.)
        peak = float(np.max(mapped))
        if peak > 0:
            mapped /= peak
        pixels = np.clip(np.rint(mapped*255), 0, 255).astype(np.uint8)
        h, w = pixels.shape
        image = QImage(pixels.tobytes(), w, h, w, QImage.Format_Grayscale8).copy()
        self.kernel_view.setPixmap(
            QPixmap.fromImage(image).scaled(
                230, 230, Qt.KeepAspectRatio, Qt.FastTransformation))
        self._update_buttons()

    def show_saved(self, output_path: str):
        self.status.setText(
            "Diagnostic ePSF + QA saved (not adopted). "
            "Original FITS were rehashed and remain unchanged.")
        self.qa.appendPlainText("\nSaved diagnostic path: " + output_path)
        # Saving a duplicate is intentionally disallowed by the immutable
        # exporter, so discourage a repeated save click in the same session.
        self._preview = None
        self._update_buttons()
