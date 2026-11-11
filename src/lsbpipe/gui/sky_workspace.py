"""Stage 03: project-bound, read-only residual sky diagnostics.

Shares the already opened Scene & Masks dataset/draft/Tractor preparation.
No network, no additional bricks, no automatic sky subtraction or adoption.
"""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import json

from PySide6.QtCore import QTimer
from PySide6.QtWidgets import (
    QApplication, QGroupBox, QHBoxLayout, QLabel, QPlainTextEdit,
    QPushButton, QVBoxLayout, QWidget, QSizePolicy,
)

from lsbpipe.background import diagnose_current_project_residual_sky


class SkyWorkspace(QWidget):
    """Explicit native Stage 03 diagnostic, never an adopted background."""

    def __init__(self, scene_editor, parent=None):
        super().__init__(parent)
        self.setObjectName("photometry_sky_stage")
        # A third scientific tab must not impose its QWidget minimum
        # size hint on the whole Photometry Workspace at narrow widths.
        self.setMinimumWidth(0)
        self.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Expanding)
        self._scene = scene_editor
        self._executor = ThreadPoolExecutor(
            max_workers=1, thread_name_prefix="desidera-sky")
        self._future = None
        self._active_identity = None
        self._result_identity = None
        self._result = None
        self._timer = QTimer(self)
        self._timer.setInterval(80)
        self._timer.timeout.connect(self._poll_worker)

        root = QVBoxLayout(self)
        root.setContentsMargins(15, 13, 15, 13)
        root.setSpacing(10)
        heading = QLabel("03 · Residual background")
        heading.setObjectName("sky_stage_title")
        heading.setStyleSheet("font-size: 19px; font-weight: 600;")
        root.addWidget(heading)
        note = QLabel(
            "Diagnostic only. DESI coadds have already been sky-subtracted. "
            "This stage estimates the remaining local background but cannot "
            "restore flux removed upstream, calibrate PSF wings, or validate "
            "the galaxy's outer envelope.")
        note.setWordWrap(True)
        root.addWidget(note)

        setup = QGroupBox("Use current Scene & Masks evidence")
        layout = QVBoxLayout(setup)
        self.context_label = QLabel("Open a galaxy project in Scene & Masks.")
        self.context_label.setObjectName("sky_project_context")
        self.context_label.setWordWrap(True)
        layout.addWidget(self.context_label)
        help_label = QLabel(
            "In 02 Scene & Masks, load the existing original FITS and "
            "prepare neighbour masks. Desidera automatically excludes "
            "the protected target, isophote candidates, catalogued "
            "neighbours and bright Gaia stars, then selects spatially "
            "separated sky training and holdout tiles from the SAME "
            "native image. No manual sky-region drawing is required. "
            "Review failed geometry or PSF-wing warnings rather than "
            "forcing a sky fit.")
        help_label.setWordWrap(True)
        layout.addWidget(help_label)
        root.addWidget(setup)

        tools = QHBoxLayout()
        self.run_btn = QPushButton("Build automatic sky")
        self.run_btn.setObjectName("sky_diagnostic_run")
        self.run_btn.setToolTip(
            "Uses only the native Scene image/invvar, current draft and "
            "provisional Tractor masks in memory; background worker. "
            "No image subtraction, model adoption or new downloads.")
        self.run_btn.clicked.connect(self.run_diagnostic)
        tools.addWidget(self.run_btn)
        self.copy_btn = QPushButton("Copy background QA")
        self.copy_btn.setObjectName("sky_diagnostic_copy_qa")
        self.copy_btn.setEnabled(False)
        self.copy_btn.clicked.connect(self.copy_qa)
        tools.addWidget(self.copy_btn)
        tools.addStretch(1)
        root.addLayout(tools)

        self.status = QLabel("Not evaluated · sky and masks remain provisional.")
        self.status.setObjectName("sky_diagnostic_status")
        self.status.setWordWrap(True)
        root.addWidget(self.status)
        self.readout = QPlainTextEdit(self)
        self.readout.setObjectName("sky_diagnostic_readout")
        self.readout.setMinimumWidth(0)
        self.readout.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Expanding)
        self.readout.setReadOnly(True)
        self.readout.setPlainText(
            "No background result. Sky training/holdout regions will be "
            "selected automatically from the current masked Scene.")
        root.addWidget(self.readout, 1)
        self.refresh_context()

    def _snapshot(self):
        """Return immutable identity plus currently loaded read-only inputs."""
        scene = self._scene
        dataset = scene._loaded
        draft = scene._draft
        contamination = scene._contamination
        if dataset is None or draft is None or contamination is None:
            return None
        if (scene._future is not None or scene._project_root is None or
                scene._project_id != dataset["project_id"] or
                scene._band != dataset["band"] or
                draft.project_id != dataset["project_id"] or
                contamination.get("project_id") != dataset["project_id"] or
                contamination.get("band") != dataset["band"] or
                contamination.get("dataset_source") != dataset["source"] or
                contamination.get("draft_digest") != draft.digest or
                contamination.get("status") !=
                "CATALOGUE_CONTAMINATION_MASK_NOT_ADOPTED"):
            return None
        identity = (
            dataset["project_id"], dataset["source"], dataset["band"],
            draft.digest, contamination.get("catalog_sha256"),
            id(dataset), id(contamination))
        return identity, dataset, draft, contamination

    def refresh_context(self):
        scene = self._scene
        snap = self._snapshot()
        band = scene._band.upper() if scene._band else "?"
        if snap is None:
            self.context_label.setText(
                f"Band {band} · No current validated-identity Scene snapshot. "
                "Load original FITS and prepare masks in 02 Scene & Masks "
                "before starting a background diagnostic.")
        else:
            identity = snap[0]
            self.context_label.setText(
                f"Current project · band {band} · "
                f"{identity[1]} · draft {identity[3][:12]} · "
                "provisional Tractor masks available.")
        if self._result is not None and (
                snap is None or snap[0] != self._result_identity):
            self._result = None
            self._result_identity = None
            self.copy_btn.setEnabled(False)
            self.readout.setPlainText(
                "Previous sky result invalidated by a change in "
                "project, band, native pixels or scene mask draft.")
            self.status.setText("Stale input · evaluate again.")
        self.run_btn.setEnabled(snap is not None and self._future is None)
        return snap is not None

    def run_diagnostic(self):
        if self._future is not None:
            return
        snap = self._snapshot()
        if snap is None:
            self.refresh_context()
            self.status.setText(
                "Cannot run: load native FITS and prepare current "
                "Tractor neighbour masks in Scene & Masks first.")
            return
        identity, dataset, draft, contamination = snap
        self._result = None
        self._result_identity = None
        self.copy_btn.setEnabled(False)
        self._active_identity = identity
        self._future = self._executor.submit(
            diagnose_current_project_residual_sky,
            dataset, draft, contamination)
        self.run_btn.setEnabled(False)
        self.status.setText(
            "Automatically selecting safe Scene sky tiles, then "
            "fitting constant and planar residual sky candidates…")
        self._timer.start()

    def _poll_worker(self):
        future = self._future
        if future is None or not future.done():
            return
        self._timer.stop()
        identity = self._active_identity
        self._active_identity = None
        self._future = None
        snap = self._snapshot()
        if snap is None or snap[0] != identity:
            self._result = None
            self._result_identity = None
            self.status.setText(
                "Sky result discarded: project/Scene evidence changed "
                "while the worker was running.")
            self.refresh_context()
            return
        try:
            result = future.result()
        except (OSError, ValueError, KeyError, TypeError,
                ArithmeticError, RuntimeError) as exc:
            self.status.setText(f"Sky diagnostic failed closed: {exc}")
            self.readout.setPlainText(
                "No background candidate accepted or subtracted. "
                "Review prepared Scene evidence and mask coverage.")
            self.refresh_context()
            return
        if (result.get("project_id") != identity[0] or
                result.get("dataset_source") != identity[1] or
                result.get("band") != identity[2] or
                result.get("draft_sha256") != identity[3] or
                result.get("catalog_sha256") != identity[4]):
            self.status.setText("Sky result source identity mismatch; discarded.")
            self.refresh_context()
            return
        self._result = result
        self._result_identity = identity
        self.copy_btn.setEnabled(True)
        warning_count = len(
            (result.get("automatic_sky_spatial_qa") or {}).get(
                "warning_flags", ()))
        self.status.setText(
            f"{result['status']} · band {identity[2].upper()} · "
            f"{warning_count} spatial quality warnings · "
            "UNADOPTED · original FITS unchanged.")
        self.readout.setPlainText(self._format_readout(result))
        self.refresh_context()

    @staticmethod
    def _format_readout(result):
        lines = [
            "Residual-sky diagnostic · UNADOPTED",
            f"Source: {result['dataset_source']}",
            f"Scene draft: {result['draft_sha256']}",
            f"Automatically chosen training pixels: "
            f"{result['sky_train_effective_pixels']}",
            f"Automatically chosen holdout pixels: "
            f"{result['sky_holdout_effective_pixels']}",
            f"Scene tile policy: "
            f"{result.get('automatic_sky_selection_qa', {}).get('strategy')}",
            f"Spatial training / holdout tiles: "
            f"{result.get('automatic_sky_selection_qa', {}).get('training_tiles')} / "
            f"{result.get('automatic_sky_selection_qa', {}).get('holdout_tiles')}",
            "Original-image pixels ≠ independent noise realizations.",
            "",
        ]
        if not result["methods"]:
            lines += [
                "No numerical sky fit performed.",
                "Automatically eligible clean Scene tiles are "
                "insufficient or source inventory is unavailable. "
                "Review the masked field; no manual sky regions required.",
            ]
        for method, model in result["methods"].items():
            coeff = model["coefficients_nMgy"]
            lines += [
                f"{method}: fit={model['fit_status']}, "
                f"holdout={model['holdout_status']}",
                f"  sky b0, bx, by [nMgy/pixel, nMgy/pixel²]: {coeff}",
                f"  training residual RMS [nMgy/pixel]: "
                f"{model['training_rms_residual_nMgy']}",
                f"  heldout residual RMS [nMgy/pixel]: "
                f"{model['holdout_rms_residual_nMgy']}",
                "  holdout quadrant median residuals "
                "[nMgy/pixel]: "
                + ", ".join(
                    f"Q{tile['quadrant']}={tile['median_residual_nMgy']:.6g}"
                    for tile in model["tiles"]),
                "",
            ]
        spatial = result.get("automatic_sky_spatial_qa")
        if spatial is not None:
            lines += [
                "Spatial tile-quality diagnostic · same coadd, NOT science PASS",
            ]
            for role, geometry in spatial["geometry"].items():
                lines += [
                    f"  {role}: dominant tile row fraction = "
                    f"{geometry['dominant_single_row_fraction']}, "
                    f"partial edge tiles = {geometry['n_partial_edge_tiles']}/"
                    f"{geometry['n_tiles']}, missing quadrants = "
                    f"{geometry['missing_quadrants']}",
                ]
            for method, audit in spatial["methods"].items():
                lines += [
                    f"  {method}: training/holdout pixel RMS ratio = "
                    f"{audit['training_to_holdout_pixel_rms_ratio']}",
                    f"  {method}: holdout tile-median range [nMgy/pixel] = "
                    f"{audit['holdout']['tile_median_peak_to_peak_nMgy_per_pixel']}",
                ]
            if spatial["warning_flags"]:
                lines += [
                    "QUALITY FLAGS — spatial support not calibrated:",
                    *["  " + flag for flag in spatial["warning_flags"]],
                ]
            lines.append("")
        cmp = result.get("method_comparison")
        if cmp:
            lines += [
                "Plane − constant at target [nMgy/pixel]: "
                + f"{cmp['plane_minus_constant_at_target_nMgy_per_pixel']:.7g}",
                "This difference is NOT a galaxy-flux error estimate.",
            ]
        preferred = result.get("provisional_automatic_model")
        if preferred is not None:
            lines += [
                f"Automatically preferred diagnostic model: "
                f"{preferred['method']}",
                f"Predicted residual sky at target [nMgy/pixel]: "
                f"{preferred['predicted_residual_sky_at_target_nMgy_per_pixel']:.7g}",
                "Selection is provisional and conditional on current Scene masks.",
            ]
        lines += [
            "",
            "Status: SKY SCIENCE HOLD. No scientifically adopted "
            "model, background subtraction, calibrated covariance or "
            "outer-galaxy light recovery.",
            "Full provenance and warnings: Copy background QA.",
        ]
        return "\n".join(lines)

    def copy_qa(self):
        if self._result is None or self._snapshot() is None or (
                self._snapshot()[0] != self._result_identity):
            self.refresh_context()
            self.status.setText(
                "No current background QA to copy. Re-evaluate.")
            return
        payload = {
            **self._result,
            "export_generated_utc": datetime.now(timezone.utc).isoformat(),
        }
        QApplication.clipboard().setText(
            json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False))
        self.status.setText(
            "Copied complete current-band residual-sky QA. "
            "No image, mask or science model changed.")

    def shutdown(self):
        self._timer.stop()
        self._executor.shutdown(wait=False, cancel_futures=True)

    def closeEvent(self, event):
        self.shutdown()
        super().closeEvent(event)
