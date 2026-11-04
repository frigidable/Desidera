"""Project-bound Photometry Workspace: preparation and provenance only.

No photometric algorithms are exposed before their validation gates exist.
"""
from __future__ import annotations

from pathlib import Path
from datetime import datetime, timezone
import json
import platform
import sys

import numpy as np

from PySide6.QtCore import Qt, Signal, QSaveFile, QIODevice
from PySide6.QtWidgets import (
    QAbstractItemView, QApplication, QFileDialog, QGroupBox, QHBoxLayout,
    QHeaderView, QLabel, QPushButton, QScrollArea, QSplitter, QTableWidget,
    QTableWidgetItem,
    QVBoxLayout, QWidget,
)

from lsbpipe.photometry.readiness import inspect_photometry_inputs
from lsbpipe.gui.psf0_panel import PSF0Panel
from lsbpipe.gui.psf0_full_qa import detailed_band_report

_STAGES = (
    ("01", "Data preparation", "Original FITS, WCS, native masks and calibration"),
    ("02", "Scene initialization", "Tentative sources, PSF star candidates and target protection"),
    ("03", "Preliminary sky", "Detection background and local star-cutout sky"),
    ("04", "PSF reconstruction", "Empirical per-band PSF, spatial QA and extended wings"),
    ("05", "Scene & smart masks", "PSF-aware contaminants and reviewed manual masks"),
    ("06", "Sky reconstruction", "Galaxy residual sky, holdouts and diffuse-light degeneracy"),
    ("07", "Refinement & QA", "Conditionally revisit PSF, masks, sky and uncertainty"),
    ("08", "Galaxy photometry", "Isophotes, structural models, colors and science results"),
)


class PhotometryWorkspace(QWidget):
    """Dedicated screen; never presents input-file presence as science-ready."""

    backRequested = Signal()
    sourcesRequested = Signal()
    psf0RunRequested = Signal(object)
    psf0SaveRequested = Signal(object)
    psf0AutoFwhmRequested = Signal(object)
    psf0TractorRequested = Signal(object)
    psf0BatchRequested = Signal(object)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("photometry_workspace")
        self._project_root = None
        self._bound_project_id = None
        self._bound_source_identity = None
        root = QVBoxLayout(self)
        root.setContentsMargins(18, 14, 18, 14)
        root.setSpacing(12)

        header = QHBoxLayout()
        self.back_btn = QPushButton("← Galaxy Workspace")
        self.back_btn.setObjectName("photometry_back_to_workspace")
        self.back_btn.clicked.connect(self.backRequested.emit)
        header.addWidget(self.back_btn)
        header.addStretch(1)
        self.sources_btn = QPushButton("Manage sources")
        self.sources_btn.setObjectName("photometry_manage_sources")
        self.sources_btn.clicked.connect(self.sourcesRequested.emit)
        header.addWidget(self.sources_btn)
        self.refresh_btn = QPushButton("Refresh saved selection")
        self.refresh_btn.setObjectName("photometry_refresh")
        self.refresh_btn.setToolTip(
            "Read project metadata and local file existence, without "
            "opening FITS, connecting to the archive or running analysis.")
        self.refresh_btn.clicked.connect(lambda _checked=False: self.refresh(force=True))
        header.addWidget(self.refresh_btn)
        root.addLayout(header)

        self.title = QLabel("Photometry")
        self.title.setObjectName("photometry_heading")
        self.title.setStyleSheet("font-size: 21px; font-weight: 600;")
        root.addWidget(self.title)
        self.subtitle = QLabel(
            "PSF-0 diagnostics are experimental · no adopted PSF, sky or galaxy photometry")
        self.subtitle.setObjectName("photometry_not_adopted")
        self.subtitle.setWordWrap(True)
        root.addWidget(self.subtitle)

        splitter = QSplitter(Qt.Horizontal)
        splitter.setChildrenCollapsible(False)

        left_scroll = QScrollArea()
        left_scroll.setWidgetResizable(True)
        left_scroll.setFrameShape(QScrollArea.NoFrame)
        left_scroll.setMinimumWidth(235)
        left_scroll.setMaximumWidth(345)
        navigation = QWidget()
        left = QVBoxLayout(navigation)
        left.setContentsMargins(2, 2, 10, 2)
        left.setSpacing(9)
        left.addWidget(QLabel("Scientific workflow"))
        for number, name, description in _STAGES:
            block = QGroupBox(f"{number} · {name}")
            block.setObjectName(f"photometry_stage_{number}")
            stage = QVBoxLayout(block)
            caption = QLabel(description)
            caption.setWordWrap(True)
            caption.setStyleSheet("font-size: 11px;")
            stage.addWidget(caption)
            left.addWidget(block)
        left.addStretch(1)
        left_scroll.setWidget(navigation)
        splitter.addWidget(left_scroll)

        center_scroll = QScrollArea()
        center_scroll.setFrameShape(QScrollArea.NoFrame)
        center_scroll.setWidgetResizable(True)
        central = QWidget()
        body = QVBoxLayout(central)
        body.setContentsMargins(8, 2, 6, 2)
        body.setSpacing(12)

        source_box = QGroupBox("Observation selection · saved project state")
        source_layout = QVBoxLayout(source_box)
        self.source_identity = QLabel("No project bound")
        self.source_identity.setObjectName("photometry_source_identity")
        self.source_identity.setWordWrap(True)
        source_layout.addWidget(self.source_identity)
        self.source_path = QLabel("")
        self.source_path.setObjectName("photometry_source_path")
        self.source_path.setWordWrap(True)
        self.source_path.setTextInteractionFlags(Qt.TextSelectableByMouse)
        source_layout.addWidget(self.source_path)
        self.source_state = QLabel("No observational input checked")
        self.source_state.setObjectName("photometry_source_state")
        self.source_state.setWordWrap(True)
        source_layout.addWidget(self.source_state)
        self.table = QTableWidget(0, 3)
        self.table.setObjectName("photometry_quartet_table")
        self.table.setHorizontalHeaderLabels(["Original FITS", "Presence", "Path"])
        self.table.horizontalHeader().setSectionResizeMode(
            0, QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(
            1, QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(
            2, QHeaderView.Stretch)
        self.table.setSelectionMode(QAbstractItemView.NoSelection)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setTextElideMode(Qt.ElideMiddle)
        self.table.verticalHeader().setVisible(False)
        self.table.setMinimumHeight(175)
        source_layout.addWidget(self.table)

        # Four independent, vertically stacked scientific panels. A band
        # never shares its preview, star exclusions or QA with another band.
        # The saved project and its original science files are unchanged.
        batch = QGroupBox("PSF-0 · four-filter workbench")
        batch.setObjectName("psf0_four_band_toolbar")
        batch_layout = QVBoxLayout(batch)
        batch_note = QLabel(
            "Each g/r/i/z panel retains its own recipe and experimental "
            "result. Build all runs available bands sequentially, with each "
            "band's own original PSFSIZE FWHM unless manually overridden. "
            "No PSF is science-adopted.")
        batch_note.setWordWrap(True)
        batch_layout.addWidget(batch_note)
        batch_buttons = QHBoxLayout()
        self.build_all_psfs_btn = QPushButton("Build PSF previews · all available griz")
        self.build_all_psfs_btn.setObjectName("psf0_build_all_bands")
        self.build_all_psfs_btn.setToolTip(
            "Runs one real native FITS/Tractor PSF-0 diagnostic for each "
            "available band, sequentially, in the existing background "
            "worker. Errors are isolated by band. No automatic downloads, "
            "science adoption or FITS export.")
        self.build_all_psfs_btn.clicked.connect(self._request_all_psfs)
        batch_buttons.addWidget(self.build_all_psfs_btn)
        self.copy_full_qa_btn = QPushButton("Copy full griz debug report")
        self.copy_full_qa_btn.setObjectName("psf0_copy_full_griz_debug")
        self.copy_full_qa_btn.clicked.connect(self.copy_full_qa)
        self.copy_full_qa_btn.setToolTip(
            "Copy full settings, all band outcomes, source evidence, "
            "star IDs, Gaia identity, sky metrics and signed per-star core/"
            "wing residual diagnostics as JSON. No native image pixels.")
        batch_buttons.addWidget(self.copy_full_qa_btn)
        self.save_full_qa_btn = QPushButton("Save full QA JSON…")
        self.save_full_qa_btn.setObjectName("psf0_save_full_griz_debug")
        self.save_full_qa_btn.clicked.connect(self.save_full_qa)
        self.save_full_qa_btn.setToolTip(
            "Export the same full four-filter diagnostic JSON as an "
            "explicit user-selected new .json file, for ChatGPT attachment. "
            "Original FITS files are never modified.")
        batch_buttons.addWidget(self.save_full_qa_btn)
        batch_buttons.addStretch(1)
        batch_layout.addLayout(batch_buttons)
        self.batch_status = QLabel("No batch run started. Original FITS unchanged.")
        self.batch_status.setObjectName("psf0_griz_batch_status")
        self.batch_status.setWordWrap(True)
        batch_layout.addWidget(self.batch_status)
        body.addWidget(batch)
        self.psf0_panels = {}
        for native_band in ("g", "r", "i", "z"):
            panel = PSF0Panel(central, fixed_band=native_band)
            panel.setObjectName("psf0_panel_" + native_band)
            panel.runRequested.connect(self.psf0RunRequested.emit)
            panel.saveRequested.connect(self.psf0SaveRequested.emit)
            panel.autoFwhmRequested.connect(self.psf0AutoFwhmRequested.emit)
            panel.tractorRequested.connect(self.psf0TractorRequested.emit)
            self.psf0_panels[native_band] = panel
            body.addWidget(panel)
        # Compatibility alias for older callers; resolves to the first
        # available band on bind, not a second shared scientific panel.
        self.psf0_panel = self.psf0_panels["g"]
        body.addWidget(source_box)

        method_box = QGroupBox("Measurement design")
        method_layout = QVBoxLayout(method_box)
        for heading, detail in (
            ("Preparatory scene",
             "Identify preliminary sources and target-protected regions; "
             "pre-PSF user mask drafts are provisional and never science-adopted."),
            ("Star-cutout sky",
             "Estimate local residual backgrounds before building an empirical PSF. "
             "Keep this distinct from the final galaxy sky solution."),
            ("PSF before final masks",
             "Compare global, local and spatial models per filter with held-out "
             "stars; evaluate extended wings and survey-specific halo subtraction."),
            ("Final scene and sky",
             "Review PSF-aware contaminant masks, refine residual galaxy sky, "
             "and conditionally revisit upstream models when QA demands."),
            ("Photometry & final evidence",
             "Measure fixed/free isophotes and PSF-convolved models; propagate "
             "mask, sky and PSF uncertainty with immutable run provenance."),
        ):
            line = QLabel(f"<b>{heading}</b> — {detail}")
            line.setTextFormat(Qt.RichText)
            line.setWordWrap(True)
            line.setMargin(5)
            method_layout.addWidget(line)
        body.addWidget(method_box)

        self.blocker = QLabel(
            "PSF-0 preview is experimental and not adopted science. "
            "No final background, masks, photometry or SBF is implemented. "
            "File presence alone never establishes source-FITS validity.")
        self.blocker.setObjectName("photometry_blocker")
        self.blocker.setWordWrap(True)
        body.addWidget(self.blocker)
        body.addStretch(1)
        center_scroll.setWidget(central)
        splitter.addWidget(center_scroll)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([275, 840])
        root.addWidget(splitter, 1)

    def set_psf_busy(self, state):
        for panel in self.psf0_panels.values():
            panel.set_busy(state)
        self.build_all_psfs_btn.setEnabled(
            not state and any(panel.run_btn.isEnabled()
                              for panel in self.psf0_panels.values()))
        self.copy_full_qa_btn.setEnabled(not state)
        self.save_full_qa_btn.setEnabled(not state)

    def _request_all_psfs(self):
        if self._project_root is None:
            self.batch_status.setText("Open Leo P in Galaxy Projects first.")
            return
        requests, skipped = [], []
        for band, panel in self.psf0_panels.items():
            request = panel.make_run_request()
            if request is None:
                skipped.append({"band":band, "reason":"NATIVE_QUARTET_UNAVAILABLE"})
            elif not panel.run_btn.isEnabled():
                skipped.append({"band":band, "reason":"CATALOGUE_OR_WORKER_NOT_READY"})
            else:
                requests.append(request)
        if not requests:
            self.batch_status.setText(
                "No complete native FITS bands with required Gaia/Tractor "
                "catalogue. Manage sources; no science run started.")
            return
        self.batch_status.setText(
            f"Queued {len(requests)} available optical bands; "
            f"unavailable: {', '.join(x['band'] for x in skipped) or 'none'}. "
            "Diagnostic models only; no FITS exports.")
        for request in requests:
            self.psf0_panels[request["band"]].prepare_run_display()
        self.psf0BatchRequested.emit({
            "project_root":str(self._project_root),
            "requests":requests,
            "skipped_bands":skipped,
        })

    def full_qa_report(self):
        from lsbpipe.photometry.psf_gui_backend import (
            PSF0_DIAGNOSTIC_EXPERIMENTS,
        )
        return {
            "schema":"desidera.psf0-griz-comprehensive-debug/1",
            "created_utc":datetime.now(timezone.utc).isoformat(),
            "project_root":str(self._project_root) if self._project_root else None,
            "science_adopted":False,
            "python":sys.version,
            "platform":platform.platform(),
            "numpy_version":np.__version__,
            "supported_experiments":list(PSF0_DIAGNOSTIC_EXPERIMENTS),
            "batch_status":self.batch_status.text(),
            "saved_dataset":self.source_identity.text(),
            "original_fits_root":self.source_path.text(),
            "selected_source_status":self.source_state.text(),
            "bands":{band:detailed_band_report(panel)
                     for band, panel in self.psf0_panels.items()},
            "no_raw_fits_or_stamp_pixels_in_clipboard":True,
        }

    def _serialize_full_qa(self):
        """Deterministic serializable JSON; never embed raw FITS stamps."""
        def native_json(value):
            if isinstance(value, np.bool_):
                return bool(value)
            if isinstance(value, np.integer):
                return int(value)
            if isinstance(value, np.floating):
                return float(value)
            if isinstance(value, np.ndarray):
                return value.tolist()
            if isinstance(value, Path):
                return str(value)
            raise TypeError(f"Unsupported QA value: {type(value).__name__}")
        return json.dumps(
            self.full_qa_report(), ensure_ascii=False, sort_keys=True,
            indent=2, allow_nan=False, default=native_json)

    def copy_full_qa(self):
        try:
            content = self._serialize_full_qa()
        except (TypeError, ValueError, OverflowError) as exc:
            self.batch_status.setText(
                f"Full QA report serialization failed: {exc}. Nothing copied.")
            return
        QApplication.clipboard().setText(content)
        self.batch_status.setText(
            f"Copied complete non-adopted PSF-0 griz debug JSON: "
            f"{len(content):,} characters. Paste it into ChatGPT. "
            "No FITS or source pixels copied or modified.")

    def save_full_qa(self):
        """Explicit atomic JSON export, not a scientific PSF product."""
        if self._project_root is None:
            self.batch_status.setText("Open Leo P first; no QA file written.")
            return
        try:
            content = self._serialize_full_qa()
        except (TypeError, ValueError, OverflowError) as exc:
            self.batch_status.setText(
                f"Cannot serialize complete QA: {exc}. No file written.")
            return
        suffix = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
        proposed = self._project_root / (
            f"Desidera_LeoP_PSF0_griz_debug_{suffix}.json")
        selected, _ = QFileDialog.getSaveFileName(
            self, "Save non-adopted PSF-0 griz QA", str(proposed),
            "JSON diagnostic report (*.json)")
        if not selected:
            return
        path = Path(selected)
        if path.suffix.lower() != ".json" or path.is_symlink() or path.exists():
            self.batch_status.setText(
                "Choose a NEW .json filename. Existing files, symlinks "
                "and original FITS will not be overwritten.")
            return
        writer = QSaveFile(str(path))
        if not writer.open(QIODevice.WriteOnly):
            self.batch_status.setText(
                "Cannot open diagnostic JSON for writing: "
                + writer.errorString())
            return
        data = content.encode("utf-8")
        if writer.write(data) != len(data) or not writer.commit():
            writer.cancelWriting()
            self.batch_status.setText(
                "QA JSON export failed. Original FITS unchanged.")
            return
        self.batch_status.setText(
            f"Full non-adopted griz debug JSON saved: {path}. "
            "Attach this .json to ChatGPT; no original FITS modified.")

    def bind_project(self, project):
        """Read the immutable UUID from the existing project manifest."""
        from lsbpipe.surveys.object_projects import read_object_project
        saved = read_object_project(project.root)
        if saved.project_id != project.project_id:
            raise ValueError("Photometry workspace project identity mismatch")
        identity = (
            str(saved.data_root) if saved.data_root is not None else None,
            tuple(saved.selected_dataset[:3]) if saved.selected_dataset else None)
        same_project = (
            self._project_root is not None
            and Path(saved.root) == self._project_root
            and saved.project_id == self._bound_project_id
            and identity == self._bound_source_identity)
        self._project_root = Path(saved.root)
        self._bound_project_id = saved.project_id
        self._bound_source_identity = identity
        # A new UUID, release/region/brick or data root invalidates PSFs.
        # Changing only the saved default band does NOT change the four
        # native per-band input datasets, so existing PSF evidence survives.
        self.refresh(force=not same_project)

    def refresh(self, *, force=False):
        # Only a new project or explicit Refresh can discard a preview.
        # Re-entering the same project/photometry view is non-destructive.
        for panel in self.psf0_panels.values():
            if force and panel._project_root is not None:
                panel._project_root = self._project_root
                panel.refresh_sources()
            else:
                panel.bind_project(self._project_root)
        self.psf0_panel = next(
            (panel for panel in self.psf0_panels.values()
             if panel.run_btn.isEnabled()), self.psf0_panels["g"])
        self.build_all_psfs_btn.setEnabled(any(
            panel.run_btn.isEnabled() for panel in self.psf0_panels.values()))
        if self._project_root is None:
            self.title.setText("Photometry")
            self.source_state.setText("Open a galaxy project first.")
            self.table.setRowCount(0)
            return
        try:
            result = inspect_photometry_inputs(self._project_root)
        except (OSError, ValueError, KeyError, TypeError) as exc:
            self.source_state.setText(f"Cannot review saved observation data: {exc}")
            self.table.setRowCount(0)
            self.blocker.setText(
                "Input review failed. No FITS or science outputs were changed.")
            return
        self.title.setText(f"Photometry · {result['galaxy']}")
        selected = result["selected_dataset"]
        identity = " / ".join(selected) if selected else "No dataset selected"
        self.source_identity.setText(
            f"Selected native source: {identity}"
            + (" · ICRS target recorded" if result["target_defined"]
               else " · target coordinates not recorded"))
        self.source_path.setText("Original FITS root: " +
                                 (result["data_root"] or "not linked"))
        self.source_state.setText(
            f"{result['status']} · {result['present_count']}/4 files present. "
            f"{result['reason']}")
        self.table.setRowCount(0)
        for product in result["products"]:
            index = self.table.rowCount()
            self.table.insertRow(index)
            for column, value in enumerate((
                    product["kind"].upper(),
                    "File found" if product["present"] else "Missing / unsafe",
                    product["path"])):
                item = QTableWidgetItem(value)
                item.setToolTip(value)
                self.table.setItem(index, column, item)
        self.blocker.setText(
            "Science state: NOT_MEASURED_NOT_ADOPTED. PSF-0 remains diagnostic. "
            "The source table checks saved paths only. Building a PSF preview "
            "validates real FITS/WCS and performs an experimental stellar fit, "
            "but PSF wings/spatial representation, galaxy sky, correlated "
            "noise and photometric measurements remain unvalidated.")
