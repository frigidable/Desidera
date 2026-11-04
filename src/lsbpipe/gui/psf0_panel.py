"""Interactive PSF-0 control embedded in the existing Galaxy Project GUI.

The researcher selects an already-present band, performs the computation
in the existing Qt worker, inspects its non-adopted diagnostics, and may
save explicit immutable derivatives. No filesystem, FITS or network side
effects occur in this widget itself.
"""
from __future__ import annotations

import json
import numpy as np

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QImage, QPixmap
from PySide6.QtWidgets import (
    QAbstractItemView, QApplication, QComboBox, QDoubleSpinBox,
    QFormLayout, QGroupBox, QHBoxLayout, QHeaderView, QLabel,
    QPlainTextEdit, QPushButton,
    QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget,
)

from lsbpipe.photometry.psf_gui_backend import (
    inspect_psf0_bands, PSF0_DIAGNOSTIC_EXPERIMENTS,
)
from lsbpipe.gui.psf_star_review import PSFStarReviewDialog


class PSF0Panel(QGroupBox):
    """Project-bound, band-specific **diagnostic** ePSF workbench."""

    runRequested = Signal(object)
    saveRequested = Signal(object)
    autoFwhmRequested = Signal(object)
    tractorRequested = Signal(object)

    def __init__(self, parent=None):
        super().__init__("04 · Empirical PSF · first pilot Leo P", parent)
        self.setObjectName("photometry_psf0_panel")
        self._project_root = None
        self._snapshot = None
        self._preview = None
        self._qa_payload = None
        self._busy = False
        self._experiment_history = []
        self._history_band = None
        self._history_policy = None
        self._star_review_data = None
        self._active_exclusions = ()

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
        catalogue_row = QHBoxLayout()
        self.catalog_policy = QComboBox()
        self.catalog_policy.setObjectName("psf0_star_catalog_policy")
        self.catalog_policy.addItem(
            "Gaia DR3 + Tractor PSF + clean-image QA (strict)", "gaia_strict")
        self.catalog_policy.addItem(
            "DAO detections only · unverified diagnostic", "dao_diagnostic")
        self.catalog_policy.currentIndexChanged.connect(self._band_changed)
        self.catalog_policy.setToolTip(
            "Gaia DR3 REF_CAT=G3 and TYPE=PSF plus native-image quality "
            "checks. Forced/frozen models, duplicates, masks, saturation, "
            "ambiguous matches and abnormal shapes are rejected. "
            "Missing Tractor catalogue blocks strict runs.")
        catalogue_row.addWidget(self.catalog_policy, 2)
        self.get_tractor_btn = QPushButton("Get DR11 Tractor + Gaia…")
        self.get_tractor_btn.setObjectName("psf0_get_tractor_gaia")
        self.get_tractor_btn.clicked.connect(self._request_tractor_catalog)
        self.get_tractor_btn.setToolTip(
            "Explicitly download the ORIGINAL DR11 south Tractor brick "
            "catalogue from the official Legacy Surveys archive into "
            "the current Leo P project's data root. No new project.")
        catalogue_row.addWidget(self.get_tractor_btn, 1)
        form.addRow("PSF star selection", catalogue_row)
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
        fwhm_row = QHBoxLayout()
        fwhm_row.addWidget(self.fwhm, 2)
        self.auto_fwhm_btn = QPushButton("From DESI PSFSIZE…")
        self.auto_fwhm_btn.setObjectName("psf0_auto_fwhm")
        self.auto_fwhm_btn.setToolTip(
            "Check the exact native DR11 coadd PSFSIZE map first. If missing "
            "or invalid, try a separately identified nearby sample of "
            "PSFSIZE values in the already-downloaded original DR11 Tractor "
            "catalogue. Every source and any map rejection is shown. "
            "Never downloads or adopts a scientific PSF.")
        self.auto_fwhm_btn.clicked.connect(self._request_auto_fwhm)
        fwhm_row.addWidget(self.auto_fwhm_btn, 1)
        form.addRow("Detection FWHM guess", fwhm_row)
        self.experiment = QComboBox()
        self.experiment.setObjectName("psf0_diagnostic_experiment")
        self.experiment.setAccessibleName("Controlled PSF reconstruction diagnostic")
        for key, info in PSF0_DIAGNOSTIC_EXPERIMENTS.items():
            self.experiment.addItem(info["label"], key)
        self.experiment.currentIndexChanged.connect(self._band_changed)
        self.experiment.setToolTip(
            "Controlled one-factor comparisons to investigate the actual "
            "Leo P outer-wing failure. These tests are not scientific PSF "
            "selection and may change the selected star sample.")
        form.addRow("Diagnostic experiment", self.experiment)
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
        self.prior_source = QLabel(
            "FWHM: 4.80 px is a USER-OVERRIDABLE DEFAULT, not measured seeing. "
            "Target protection: project radius with 90″ conservative floor; "
            "not inferred from DESI seeing.")
        self.prior_source.setWordWrap(True)
        self.prior_source.setObjectName("psf0_prior_provenance")
        layout.addWidget(self.prior_source)
        self.gaia_source = QLabel(
            "No Gaia/Tractor source checked. Strict mode requires the original "
            "catalogue for this DR11 brick; Gaia identity does not prove "
            "an uncontaminated stellar cutout.")
        self.gaia_source.setObjectName("psf0_gaia_catalog_source")
        self.gaia_source.setWordWrap(True)
        layout.addWidget(self.gaia_source)

        actions = QHBoxLayout()
        self.run_btn = QPushButton("Build PSF preview")
        self.run_btn.setObjectName("psf0_run_button")
        self.run_btn.setToolTip(
            "Runs real FITS validation, star selection and ePSF calculation "
            "in the existing background worker; does not write a PSF.")
        self.run_btn.clicked.connect(self._request_run)
        actions.addWidget(self.run_btn)
        self.review_btn = QPushButton("Review training stars…")
        self.review_btn.setObjectName("psf0_review_stars_button")
        self.review_btn.setToolTip(
            "Open a separate mosaic of all train cutouts, aligned candidate "
            "PSF and square radial profiles. Exclude stars and recompute "
            "without changing independent holdouts.")
        self.review_btn.clicked.connect(self._open_star_review)
        actions.addWidget(self.review_btn)
        self.save_btn = QPushButton("Save diagnostic FITS + report")
        self.save_btn.setObjectName("psf0_save_button")
        self.save_btn.setEnabled(False)
        self.save_btn.setToolTip(
            "Only after a successful preview: verify original FITS hashes "
            "again and save an immutable non-adopted PSF derivative.")
        self.save_btn.clicked.connect(self._request_save)
        actions.addWidget(self.save_btn)
        self.copy_qa_btn = QPushButton("Copy QA report")
        self.copy_qa_btn.setObjectName("psf0_copy_qa")
        self.copy_qa_btn.setToolTip(
            "Copies the scientific result/rejection and star audit to the "
            "clipboard; no invalid PSF is saved or adopted.")
        self.copy_qa_btn.clicked.connect(self._copy_qa)
        actions.addWidget(self.copy_qa_btn)
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
        self.compare = QLabel(
            "Experiment comparison (same band): no completed trials.")
        self.compare.setObjectName("psf0_experiment_comparison")
        self.compare.setWordWrap(True)
        self.compare.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.compare.setToolTip(
            "Historical non-adopted diagnostic trials on the current band. "
            "Different cutout sizes can change the stellar sample; no "
            "scientific best model is selected automatically.")
        layout.addWidget(self.compare)

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
        self.stars = QTableWidget(0, 6)
        self.stars.setObjectName("psf0_star_audit")
        self.stars.setHorizontalHeaderLabels(
            ["Star ID", "Role", "X (pix)", "Y (pix)", "Cutout sky", "Gradient/RMS"])
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
        # Before a galaxy is bound, no action may appear runnable.
        self._update_buttons()

    @property
    def active_preview(self):
        return self._preview

    def bind_project(self, project_root):
        self._project_root = project_root
        self.refresh_sources()

    def refresh_sources(self):
        self._star_review_data = None
        self._active_exclusions = ()
        self._experiment_history = []
        self._history_band = None
        self._history_policy = None
        self._render_comparison()
        self._preview = None
        self._qa_payload = None
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
            self.gaia_source.setText(
                ("Original DR11 Tractor FITS found (not yet validated): "
                 + str(snapshot.get("gaia_tractor_path"))
                 if snapshot.get("gaia_tractor_present_unverified")
                 else "No original DR11 Tractor+Gaia catalogue locally; "
                      "strict Gaia selection needs a deliberate GUI download."))
            self.catalog_policy.blockSignals(True)
            self.catalog_policy.setCurrentIndex(
                self.catalog_policy.findData(
                    "gaia_strict" if snapshot.get("gaia_tractor_present_unverified")
                    else "dao_diagnostic"))
            self.catalog_policy.blockSignals(False)
            self.prior_source.setText(
                "FWHM: initial 4.80 px code default (not measured). "
                f"Target protection: max(saved project radius, 90″) = "
                f"{snapshot['target_protection_default_arcsec']:.1f}″. "
                "The galaxy protection radius cannot be inferred from PSF seeing.")
            self.status.setText(
                "Original FITS paths found, but not yet validated. "
                "Click Build PSF preview to verify FITS/WCS/units and run.")
            if not available:
                self.status.setText(
                    "No complete local native FITS quartet for PSF-0. "
                    "Use Manage sources; no automatic download.")
        self.band.blockSignals(False)
        self._history_band = self.band.currentData()
        self._history_policy = self.catalog_policy.currentData()
        self._update_buttons()

    def _parameters_changed(self, *_):
        # The FWHM/target radius changes the detected candidates,
        # spatial split and the PSF's scientific recipe. Historical
        # same-band trials from the prior parameter set are NOT
        # controlled comparisons; invalidate them immediately.
        self._experiment_history = []
        self._history_band = self.band.currentData()
        self._history_policy = self.catalog_policy.currentData()
        self._render_comparison()
        self._active_exclusions = ()
        if self._preview is not None or self._qa_payload is not None:
            self._band_changed()

    def _band_changed(self, *_):
        self._star_review_data = None
        self._active_exclusions = ()
        # Changing filter/policy changes the source sample. Clear old
        # comparisons when identities are not scientifically comparable.
        # Changing the experiment retains same-band trials, as before.
        # Changing filter invalidates former-band comparison. Changing the
        # experiment retains historical same-band trials in memory.
        band = self.band.currentData()
        policy = self.catalog_policy.currentData()
        if band != self._history_band or policy != self._history_policy:
            self._history_band = band
            self._history_policy = policy
            self._experiment_history = []
            self._render_comparison()
        self._preview = None
        self._qa_payload = None
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
        catalog_ready = (
            self.catalog_policy.currentData() != "gaia_strict" or
            bool(self._snapshot.get("gaia_tractor_present_unverified"))
            if self._snapshot is not None else False)
        self.run_btn.setEnabled(eligible and catalog_ready and not self._busy)
        self.save_btn.setEnabled(self._preview is not None and not self._busy
                                 and self._preview["result"].lineage_id is not None)
        self.band.setEnabled(not self._busy)
        self.fwhm.setEnabled(not self._busy)
        self.protection.setEnabled(not self._busy)
        self.experiment.setEnabled(not self._busy)
        self.catalog_policy.setEnabled(not self._busy)
        self.get_tractor_btn.setEnabled(
            not self._busy and bool(self._snapshot)
            and self._snapshot.get("dataset", ())[0:2] == ("dr11", "south"))
        self.auto_fwhm_btn.setEnabled(eligible and not self._busy)
        self.copy_qa_btn.setEnabled(self._qa_payload is not None and not self._busy)
        self.review_btn.setEnabled(self._star_review_data is not None and not self._busy)

    def set_busy(self, state: bool, message: str = ""):
        self._busy = bool(state)
        if message:
            self.status.setText(message)
        self._update_buttons()

    def _request_run(self, *, review_confirmed: bool = False):
        if self._project_root is None or self._busy or self.band.currentData() is None:
            return
        self._star_review_data = None
        self._preview = None
        self._qa_payload = None
        self.save_btn.setEnabled(False)
        self.stars.setRowCount(0)
        self.kernel_view.setPixmap(QPixmap())
        self.kernel_view.setText("PSF calculation running…")
        self.summary.setText("No accepted PSF; calculating provisional model.")
        self.qa.clear()
        self._update_buttons()
        self.runRequested.emit({
            "project_root": str(self._project_root),
            "band": self.band.currentData(),
            "fwhm_pix": self.fwhm.value(),
            "experiment": self.experiment.currentData(),
            "star_catalog_policy": self.catalog_policy.currentData(),
            "excluded_training_star_ids": tuple(self._active_exclusions),
            "manual_review_confirmed": review_confirmed,
            "target_protection_arcsec": self.protection.value(),
        })

    def _request_tractor_catalog(self):
        if self._busy or self._project_root is None or self._snapshot is None:
            return
        self.tractorRequested.emit({"project_root":str(self._project_root)})

    def show_tractor_acquired(self, payload):
        if (self._project_root is None or
                str(self._project_root) != str(payload.get("project_root"))):
            self.status.setText("Tractor catalogue belongs to another project")
            return
        self._snapshot = inspect_psf0_bands(self._project_root)
        self.catalog_policy.setCurrentIndex(
            self.catalog_policy.findData("gaia_strict"))
        self._band_changed()
        self.gaia_source.setText(
            f"Native DR11 Tractor+Gaia catalogue verified: "
            f"{payload['n_catalog_sources']} brick-primary sources; "
            f"SHA256 {payload['sha256_local'][:16]}… "
            "Gaia association is necessary but not sufficient for clean PSF.")
        self.status.setText(
            "Gaia-strict selection enabled. Build a fresh PSF preview "
            "to reject nonstellar and contaminated candidates.")
        self._update_buttons()

    def _request_auto_fwhm(self):
        if self._busy or self._project_root is None or self.band.currentData() is None:
            return
        self.autoFwhmRequested.emit({
            "project_root": str(self._project_root),
            "band": self.band.currentData(),
        })

    def show_fwhm_suggestion(self, payload):
        if (self._project_root is None or
                str(self._project_root) != str(payload.get("project_root")) or
                self.band.currentData() != payload.get("band")):
            self.status.setText(
                "PSFSIZE suggestion belongs to a different project/band; discarded.")
            return
        seeing = payload["seeing"]
        if seeing.get("available") is not True:
            message = seeing.get(
                "message", "Neither exact coadd map nor nearby Tractor "
                "PSFSIZE values can be verified")
            self.status.setText(message)
            self.prior_source.setText(
                "PSFSIZE lookup unsuccessful — " + message +
                " FWHM unchanged. Verify project native-image path and the "
                "matching original map/Tractor FITS; no silent default.")
            return
        pix = float(seeing["suggested_detection_fwhm_pix"])
        source_id = seeing.get("source", "UNKNOWN")
        via_tractor = source_id == (
            "LOCAL_DR11_TRACTOR_NEARBY_PSF_SIZE_NOT_PIXEL_MAP")
        source_name = (
            "nearby DR11 Tractor PSFSIZE band-column median"
            if via_tractor else "original local coadd PSFSIZE pixel map")
        self.fwhm.setValue(round(pix, 2))
        resolved = self.fwhm.value()
        caveat = (
            " Catalog positions approximate the local PSF prior; not "
            "a per-pixel seeing measurement at Leo P."
            if via_tractor else "")
        map_note = seeing.get("native_map_status")
        self.prior_source.setText(
            f"FWHM suggested from {source_name}: "
            f"{seeing['psfsize_arcsec']:.3f}″ / "
            f"{seeing['image_pixel_scale_arcsec']:.4f} arcsec/pix "
            f"= {pix:.3f} px (GUI {resolved:.2f} px). "
            f"Evidence: {seeing.get('sample', 'source sample unknown')}; "
            f"file {seeing.get('path', 'unknown')}. "
            + (f"Original map status: {map_note}. " if map_note else "")
            + caveat
            + " Not an empirical ePSF; target protection unchanged.")
        self.status.setText(
            f"Applied DESI {source_name} to the initial detection FWHM "
            f"guess: {resolved:.2f} px. Review and rebuild PSF.")

    def _open_star_review(self):
        if self._busy or self._star_review_data is None:
            return
        dialog = PSFStarReviewDialog(self._star_review_data, parent=self)
        if dialog.exec() != PSFStarReviewDialog.Accepted:
            return
        selected = dialog.excluded_ids()
        if selected == self._active_exclusions:
            return
        self._active_exclusions = selected
        self._request_run(review_confirmed=True)

    def _request_save(self):
        if self._busy or self._preview is None:
            return
        self.saveRequested.emit(self._preview)

    def show_failure(self, error: str):
        self._star_review_data = None
        self._preview = None
        self._qa_payload = None
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

    def _attach_star_review(self, payload: dict):
        """Expose saved checkboxes for review, but don't auto-alter PSF."""
        current = payload.get("star_review")
        if current is None:
            self._star_review_data = None
            return
        saved = payload.get("saved_manual_review")
        if saved and not current.get("excluded_training_star_ids"):
            # Data has already been SHA- and recipe-matched in the worker.
            # The researcher must still inspect and click Apply; baseline
            # remains exactly the default algorithm until then.
            self._star_review_data = {
                **current,
                "excluded_training_star_ids": tuple(
                    saved["excluded_training_star_ids"]),
            }
            note = (" Previous training-star exclusions were found for "
                    "these original FITS; inspect them in Review training "
                    "stars. They were NOT applied automatically.")
            if note.strip() not in self.prior_source.text():
                self.prior_source.setText(self.prior_source.text() + note)
        else:
            self._star_review_data = current

    @staticmethod
    def _catalog_screening_summary(catalog):
        # Full per-star score/iteration evidence is in Copy QA report.
        # The on-screen panel must show the short catalog provenance
        # separately, not dump 40–100 rows of nested trial statistics.
        omit = {
            "training_residual_screening_history",
            "automatic_training_rejections",
        }
        return str({key:value for key,value in catalog.items()
                    if key not in omit})

    @staticmethod
    def _training_sky_qa_summary(metrics):
        sky = metrics.get("sky_sample_qa", {})
        passes = metrics.get("training_residual_screening_history", [])
        stars = sky.get("per_star", [])
        example = next((x for x in stars
                        if "annulus_inner_in_guessed_fwhm" in x), {})
        ratios = (
            f"{example['annulus_inner_in_guessed_fwhm']:.2f}–"
            f"{example.get('annulus_outer_in_guessed_fwhm', 0):.2f}"
            if example else "unknown")
        notes = [
            "PSF background QA (not a joint PSF+sky fit):",
            f"Sky annulus starts/ends at {ratios} × guessed detection FWHM.",
            "Stellar wings can bias the fitted sky even after clipping; "
            "coadd noise covariance is not modeled.",
            f"Sky ring inside 4× guessed FWHM: "
            f"{sky.get('ring_starts_inside_four_assumed_fwhm','unknown')}.",
            "Training-only residual outlier screening (not χ²):",
        ]
        for trial in passes:
            notes.append(
                f"Pass {trial.get('iteration','?')}: median score "
                f"{trial.get('median_score','—')}, "
                f"cutoff {trial.get('cutoff','—')}, "
                f"auto-removed {len(trial.get('selected_for_exclusion',[]))}; "
                "holdouts excluded from the decision.")
            bias = trial.get("brightness_score_coupling_qa", {})
            rho = bias.get("score_amplitude_spearman_rho")
            if rho is not None:
                notes.append(
                    f"Star score versus model amplitude: Spearman "
                    f"rho={rho:+.3f} (descriptive; not a selection test).")
            if bias.get("auto_trim_blocked"):
                proposed = trial.get("proposed_before_brightness_guard", [])
                notes.append(
                    "AUTOMATIC TRAINING EXCLUSION BLOCKED: residual "
                    "score tracks star amplitude too strongly. "
                    f"Candidates {', '.join(proposed)} remain in training; "
                    "inspect individually and use a brightness-independent "
                    "selection model before retrying.")
                notes.append(
                    "The current model is the UNTRIMMED training result; "
                    "no iterative star-removal improvement is claimed.")
        if not passes:
            notes.append("Disabled in this experiment.")
        return "\n".join(notes)

    def show_quality_rejection(self, payload: dict):
        """Rejected ePSF: inspect its signed structure, never adopt a kernel."""
        if (self._project_root is None or
                str(self._project_root) != str(payload.get("project_root"))):
            self.show_failure("Rejected PSF belongs to another Galaxy Project")
            return
        self._preview = None
        self._attach_star_review(payload)
        self._qa_payload = {
            "science_status": "PSF0_QUALITY_REJECTED_NOT_ADOPTED",
            "source_project": str(self._project_root),
            "band": self.band.currentData(),
            "experiment": self.experiment.currentData(),
            "reason": payload.get("reason"),
            "diagnostics": {
                name: value for name, value in payload.get("diagnostics", {}).items()
                if name != "signed_kernel_values_for_qa_only"
            },
        }
        self.kernel_view.setPixmap(QPixmap())
        self.kernel_view.setText("No signed model available")
        self.status.setText(
            "PSF-0 QUALITY REJECTED · no PSF model adopted or saved")
        diagnostic = payload.get("diagnostics", {})
        catalog_qa = diagnostic.get("catalog_selection_qa", {})
        negative = diagnostic.get("negative_absolute_mass_fraction")
        if negative is not None:
            self.summary.setText(
                f"Negative absolute ePSF mass: {negative:.2%} "
                f"(experimental limit {diagnostic.get('blocking_threshold_fraction', .025):.2%}); "
                f"training stars {diagnostic.get('n_train', 0)}, "
                f"holdout stars {diagnostic.get('n_holdout', 0)}. "
                "Review the source selection and local backgrounds. "
                f"Policy: {catalog_qa.get('policy', 'unknown')}.")
        else:
            self.summary.setText("Provisional PSF rejected before adoption.")
        rows = [(item, "train") for item in diagnostic.get("training", [])]
        rows += [(item, "holdout") for item in diagnostic.get("holdout", [])]
        self.stars.setRowCount(len(rows))
        for index, (star, role) in enumerate(rows):
            for col, label in enumerate((
                    str(star.get("star_id", "unknown")),
                    role, f"{star['x_pix']:.2f}", f"{star['y_pix']:.2f}",
                    f"{star['cutout_sky_nMgy']:.4g}",
                    f"{star['gradient_sigma_span']:.2f}")):
                cell = QTableWidgetItem(label)
                if "ra_deg" in star and "dec_deg" in star:
                    cell.setToolTip(
                        f"ICRS RA={star['ra_deg']:.7f}°, Dec={star['dec_deg']:.7f}°")
                self.stars.setItem(index, col, cell)
        signed = diagnostic.get("signed_kernel_qa", {})
        holdouts = diagnostic.get("signed_holdout_summary", {})
        holdout_regions = holdouts.get("regions", {})
        annuli = signed.get("annuli", {})
        breakdown = []
        for label, key in (
                ("Core (r ≤ 1× guessed FWHM)", "core_r_le_1x_guessed_fwhm"),
                ("Intermediate (1–2× guessed FWHM)",
                 "intermediate_1x_to_2x_guessed_fwhm"),
                ("Outer (r > 2× guessed FWHM)",
                 "outer_r_gt_2x_guessed_fwhm")):
            if key in annuli:
                breakdown.append(
                    f"{label}: negative mass contribution "
                    f"{annuli[key]['negative_fraction_of_total_absolute_mass']:.2%} "
                    "of total |model|.")
        holdout_lines = []
        for label, key in (
                ("Core", "core_r_le_1x_guessed_fwhm"),
                ("Intermediate", "intermediate_1x_to_2x_guessed_fwhm"),
                ("Outer", "outer_r_gt_2x_guessed_fwhm")):
            zone = holdout_regions.get(key, {})
            absolute = zone.get("median_star_abs_residual_over_peak")
            signed_value = zone.get("median_star_signed_residual_over_peak")
            if absolute is not None and signed_value is not None:
                holdout_lines.append(
                    f"{label}: median |residual| / star peak "
                    f"{absolute:.4g}; signed {signed_value:+.4g}; "
                    f"stars {zone.get('n_stars', 0)}")
        raw = diagnostic.get("signed_kernel_values_for_qa_only")
        if raw is not None:
            self._show_signed_rejected_kernel(raw)
        solver = diagnostic.get("solver_diagnostics", {})
        stable = solver.get("final_converged_fraction_if_reported")
        stable_target = solver.get("convergence_fraction_target_if_reported")
        last_step = (solver.get("iteration_history_if_reported") or [{}])[-1]
        solver_text = (
            f"\n\nStable-centroid fraction: "
            f"{stable:.1%}" if stable is not None else
            "\n\nStable-centroid fraction: unavailable")
        solver_text += (
            f" (required ≥ {stable_target:.1%})"
            if stable_target is not None else "")
        solver_text += (
            f"\nLast relative center shift (pix): "
            f"{last_step['max_center_shift']:.4g}"
            if last_step.get("max_center_shift") is not None else "")
        solver_text += (
            f"\nPhotutils version: {solver.get('photutils_version', 'unknown')}"
            f"\nBuilder smoothing: "
            f"{solver.get('smoothing_kernel_final_if_reported', 'unknown')}"
            f"\nBuilder fit box: {solver.get('fit_shape_final_if_reported', 'unknown')}"
            f"\nBuilder wing smoothing: "
            f"{solver.get('wing_smoothing_requested', 'unknown')}"
            f"\nIterations run: "
            f"{solver.get('iterations_performed_if_reported', 'unknown')}"
            f"\nBuilder converged: {solver.get('converged_if_reported', 'unknown')}")
        self.qa.setPlainText(
            payload.get("reason", "PSF reconstruction quality failure")
            + "\n\nSigned ePSF preview: cyan = positive, red = negative; "
            "nonlinear contrast for diagnosis only."
            + "\n\nNegative contribution by radial zone "
            "(radii use an assumed detection FWHM, NOT a measured FWHM):\n"
            + ("\n".join(breakdown) if breakdown else "Unavailable")
            + "\n\nIndependent holdout residuals (core-fitted signed model; "
            "not covariance-calibrated):\n"
            + ("\n".join(holdout_lines)
               if holdout_lines else "Unavailable")
            + solver_text
            + "\n\nStar catalogue screening: "
            + self._catalog_screening_summary(catalog_qa)
            + "\n\n" + self._training_sky_qa_summary(diagnostic)
            + "\n\nRejected candidates: "
            + str(diagnostic.get("rejected_reasons", {}))
            + "\n\n"
            + "\n".join(diagnostic.get("notes", []))
            + "\n\nNo invalid PSF kernel or photometric measurement was saved.")
        self._record_experiment(
            n_train=diagnostic.get("n_train", 0),
            n_holdout=diagnostic.get("n_holdout", 0),
            negative=negative,
            converged=diagnostic.get(
                "solver_diagnostics", {}).get("converged_if_reported"),
            status="QUALITY REJECTED",
            outer_holdout_abs=holdout_regions.get(
                "outer_r_gt_2x_guessed_fwhm", {}
            ).get("median_star_abs_residual_over_peak"))
        self._update_buttons()

    def _show_signed_rejected_kernel(self, values):
        """Red = negative model, cyan = positive; descriptive, not accepted."""
        data = np.asarray(values, dtype=float)
        if (data.ndim != 2 or min(data.shape) < 9
                or not np.isfinite(data).all()):
            self.kernel_view.setText("Signed model grid unavailable")
            return
        peak = float(np.max(np.abs(data)))
        if peak <= 0:
            self.kernel_view.setText("Signed model grid is empty")
            return
        stretch = np.arcsinh(100. * np.abs(data) / peak) / np.arcsinh(100.)
        brightness = np.clip(np.rint(205. * stretch), 0, 205).astype(np.int16)
        rgb = np.zeros((*data.shape, 3), dtype=np.uint8)
        rgb[:] = (34, 39, 45)
        positive = data > 0
        negative = data < 0
        rgb[positive, 0] = 25
        rgb[positive, 1] = 50 + brightness[positive] * 3 // 4
        rgb[positive, 2] = 50 + brightness[positive]
        rgb[negative, 0] = 50 + brightness[negative]
        rgb[negative, 1] = 39
        rgb[negative, 2] = 45
        h, w = data.shape
        img = QImage(rgb.tobytes(), w, h, w * 3, QImage.Format_RGB888).copy()
        self.kernel_view.setPixmap(QPixmap.fromImage(img).scaled(
            230, 230, Qt.KeepAspectRatio, Qt.FastTransformation))
        self.kernel_view.setToolTip(
            "REJECTED signed ePSF · cyan positive; red negative. "
            "Asinh-stretched for visual diagnosis; no physically accepted PSF.")

    def _record_experiment(self, *, n_train, n_holdout, negative,
                           converged, status, outer_holdout_abs=None):
        entry = {
            "band": self.band.currentData(),
            "experiment": self.experiment.currentData(),
            "label": self.experiment.currentText(),
            "status": status,
            "negative_fraction": (
                float(negative) if negative is not None else None),
            "solver_converged": converged,
            "n_train": int(n_train),
            "n_holdout": int(n_holdout),
            "holdout_outer_median_abs_residual_over_peak": (
                float(outer_holdout_abs)
                if outer_holdout_abs is not None else None),
            "science_adopted": False,
        }
        self._experiment_history.append(entry)
        self._experiment_history = self._experiment_history[-10:]
        self._render_comparison()

    def _render_comparison(self):
        if not self._experiment_history:
            self.compare.setText(
                "Experiment comparison (same band): no completed trials.")
            return
        descriptions = []
        for entry in self._experiment_history:
            pct = (f"{entry['negative_fraction']:.2%}"
                   if entry["negative_fraction"] is not None else "N/A")
            converged = ("yes" if entry["solver_converged"] is True else
                         "no" if entry["solver_converged"] is False else "?")
            outer = entry.get("holdout_outer_median_abs_residual_over_peak")
            outer_text = f"{outer:.3g}" if outer is not None else "N/A"
            descriptions.append(
                f"{entry['experiment']}: neg {pct}, "
                f"conv {converged}, "
                f"outer holdout |resid|/peak {outer_text}, "
                f"stars {entry['n_train']}/{entry['n_holdout']}")
        self.compare.setText(
            "Non-adopted same-band trials · "
            + " | ".join(descriptions)
            + ". Different cutouts may select different stars.")

    def _copy_qa(self):
        if self._qa_payload is None or self._busy:
            return
        payload = dict(self._qa_payload)
        payload["comparison_trials_not_adopted"] = [
            dict(trial) for trial in self._experiment_history]
        QApplication.clipboard().setText(json.dumps(
            payload, ensure_ascii=False, indent=2,
            sort_keys=True, allow_nan=False))
        self.qa.appendPlainText("\nQA report copied to clipboard (no FITS written).")

    def show_result(self, preview: dict):
        if (self._project_root is None or
                str(self._project_root) != str(preview["project_root"])):
            self.show_failure("Result project is different from the open Galaxy Project")
            return
        result = preview["result"]
        self._preview = preview
        self._attach_star_review(preview)
        report = result.report()
        self._qa_payload = report
        self.status.setText(
            f"{result.dataset_id} · EXPERIMENTAL_PSF_NOT_ADOPTED")
        self.summary.setText(
            f"Star selection: {self.catalog_policy.currentText()} · "
            f"Experiment: {self.experiment.currentText()} · "
            f"Training stars: {len(result.training)} · Independent holdout: "
            f"{len(result.holdout)} · Detected candidates: {result.detection_count} · "
            f"Median heldout core |residual|/peak: "
            f"{result.metrics['holdout_median_abs_core_residual_over_peak']:.3g}. "
            "This is NOT a calibrated science PSF or FWHM.")
        rows = result.training + result.holdout
        self.stars.setRowCount(len(rows))
        for i, star in enumerate(rows):
            for j, text in enumerate((
                    star.star_id, star.role, f"{star.x_pix:.2f}",
                    f"{star.y_pix:.2f}", f"{star.cutout_sky_nMgy:.4g}",
                    f"{star.gradient_sigma_span:.2f}")):
                item = QTableWidgetItem(text)
                item.setToolTip(
                    f"ICRS RA={star.ra_deg:.7f}°, Dec={star.dec_deg:.7f}°")
                self.stars.setItem(i, j, item)
        diagnostics = result.metrics.get(
            "signed_epsf_diagnostics_before_clipping", {})
        self._record_experiment(
            n_train=len(result.training), n_holdout=len(result.holdout),
            negative=diagnostics.get("negative_absolute_mass_fraction"),
            converged=result.metrics.get(
                "solver_diagnostics", {}).get("converged_if_reported"),
            status="EXPERIMENTAL PREVIEW",
            outer_holdout_abs=result.metrics.get(
                "signed_holdout_summary", {}).get("regions", {}).get(
                "outer_r_gt_2x_guessed_fwhm", {}).get(
                    "median_star_abs_residual_over_peak"))
        self.qa.setPlainText(
            "Scientific limits:\n" + "\n".join(f"• {f}" for f in result.flags)
            + "\n\nCatalog screening: "
            + self._catalog_screening_summary(
                result.metrics.get("catalog_selection_qa", {}))
            + "\n\n" + self._training_sky_qa_summary(result.metrics)
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

    def show_review_saved(self, receipt_path: str):
        """Distinguish immutable manual star selection from PSF adoption."""
        if self._qa_payload is not None:
            self._qa_payload["manual_training_review_receipt"] = receipt_path
        self.qa.appendPlainText(
            "\nManual training-star review saved (NOT an adopted PSF): "
            + receipt_path)
        self.status.setText(
            self.status.text() + " · manual train exclusions saved")

    def show_saved(self, output_path: str):
        self.status.setText(
            "Diagnostic ePSF + QA saved (not adopted). "
            "Original FITS were rehashed and remain unchanged.")
        self.qa.appendPlainText("\nSaved diagnostic path: " + output_path)
        # Saving a duplicate is intentionally disallowed by the immutable
        # exporter, so discourage a repeated save click in the same session.
        self._preview = None
        self._update_buttons()
