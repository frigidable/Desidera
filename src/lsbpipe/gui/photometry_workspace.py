"""Project-bound Photometry Workspace: experimental science stages.

Stage 04 exposes only finite-aperture observed-area diagnostics; no
scientific sky, PSF, total luminosity, or color is adopted.
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
    QComboBox, QHeaderView, QLabel, QPushButton, QScrollArea, QTabBar, QTabWidget, QTableWidget,
    QTableWidgetItem, QVBoxLayout, QWidget, QSizePolicy,
)

from lsbpipe.photometry.readiness import inspect_photometry_inputs
from lsbpipe.gui.psf0_panel import PSF0Panel
from lsbpipe.gui.psf0_full_qa import detailed_band_report
from lsbpipe.gui.scene_editor import SceneEditor
from lsbpipe.gui.sky_workspace import SkyWorkspace
from lsbpipe.gui.aperture_workspace import ApertureWorkspace
from lsbpipe.gui.project_qa_bundle import serialize_project_diagnostic_bundle



class PhotometryWorkspace(QWidget):
    """Dedicated screen; never presents input-file presence as science-ready."""

    backRequested = Signal()
    sourcesRequested = Signal()
    psf0RunRequested = Signal(object)
    psf0SaveRequested = Signal(object)
    psf0AutoFwhmRequested = Signal(object)
    psf0TractorRequested = Signal(object)
    psf0BatchRequested = Signal(object)
    navigationChanged = Signal(str, str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("photometry_workspace")
        self._project_root = None
        self._bound_project_id = None
        self._bound_source_identity = None
        self._suppress_navigation = False
        self._tab_source_identity = None
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
        self.copy_project_qa_btn = QPushButton("Copy QA")
        self.copy_project_qa_btn.setObjectName("photometry_copy_project_qa")
        self.copy_project_qa_btn.setToolTip(
            "Copy a single current-project JSON report: PSF-0 griz, Scene & "
            "Masks, residual background, aperture curves/systematics, stale "
            "evidence, blockers and environment. Missing analyses are marked "
            "NOT_RUN; this never opens FITS or starts workers.")
        self.copy_project_qa_btn.clicked.connect(self.copy_project_qa)
        root.addLayout(header)

        self.title = QLabel("Photometry")
        self.title.setObjectName("photometry_heading")
        self.title.setStyleSheet("font-size: 21px; font-weight: 600;")
        title_row = QHBoxLayout()
        title_row.addWidget(self.title)
        title_row.addStretch(1)
        title_row.addWidget(self.copy_project_qa_btn)
        root.addLayout(title_row)
        self.subtitle = QLabel(
            "PSF, sky and aperture fluxes are diagnostic only · "
            "no adopted galaxy photometry")
        self.subtitle.setObjectName("photometry_not_adopted")
        self.subtitle.setWordWrap(True)
        root.addWidget(self.subtitle)
        self.project_qa_status = QLabel("")
        self.project_qa_status.setObjectName("photometry_project_qa_status")
        self.project_qa_status.setWordWrap(True)
        root.addWidget(self.project_qa_status)
        self.project_qa_status.hide()

        # These are two real user workflows in ONE object-bound
        # Photometry Workspace, not eight non-executable roadmap cards.
        self.stage_tabs = QTabWidget(self)
        self.stage_tabs.setObjectName("photometry_science_stage_tabs")
        psf_stage = QWidget(self.stage_tabs)
        psf_stage.setObjectName("photometry_psf_stage")
        psf_stage_layout = QVBoxLayout(psf_stage)
        psf_stage_layout.setContentsMargins(0, 0, 0, 0)
        psf_stage_layout.setSpacing(7)
        self.stage_tabs.addTab(psf_stage, "1 · PSF-0")
        # One persistent selector, not a tab strip that disappears behind
        # the scroll area when the user examines PSF settings or QA.
        self.band_selector = QTabBar(self)
        self.band_selector.setObjectName("psf0_sticky_band_selector")
        self.band_selector.setDocumentMode(True)
        self.band_selector.setExpanding(False)
        self.band_selector.setToolTip(
            "Select g/r/i/z to inspect its independent PSF model, "
            "radial comparison and scientific rejection status.")
        psf_stage_layout.addWidget(self.band_selector)

        center_scroll = QScrollArea()
        center_scroll.setFrameShape(QScrollArea.NoFrame)
        center_scroll.setMinimumWidth(0)
        center_scroll.setWidgetResizable(True)
        center_scroll.setAlignment(Qt.AlignHCenter | Qt.AlignTop)
        central = QWidget()
        central.setObjectName("photometry_narrow_science_column")
        # Fill the actual space beside the persistent survey sidebar.
        # The former 980-px cap created wide, useless outer gutters on
        # standard monitors; inner controls retain their own width limits.
        central.setMaximumWidth(1400)
        central.setMinimumWidth(0)
        central.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Maximum)
        body = QVBoxLayout(central)
        body.setContentsMargins(8, 2, 8, 2)
        body.setSpacing(12)
        mode_row = QHBoxLayout()
        mode_caption = QLabel("Controls")
        mode_caption.setObjectName("psf0_workspace_view_caption")
        mode_row.addWidget(mode_caption)
        self.view_mode = QComboBox()
        self.view_mode.setObjectName("psf0_workspace_view_mode")
        self.view_mode.addItem("Standard", "guided")
        self.view_mode.addItem("Advanced diagnostics", "research")
        self.view_mode.setToolTip(
            "Changes which controls are visible, never the PSF recipe, "
            "source FITS, independent holdouts or scientific QA.")
        mode_row.addWidget(self.view_mode)
        mode_row.addStretch(1)
        body.addLayout(mode_row)
        source_box = QGroupBox("Observation selection · saved project state")
        self.source_box = source_box
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
        self.source_details_toggle = QPushButton(
            "Observation and original FITS details  ▸")
        self.source_details_toggle.setObjectName(
            "photometry_source_details_toggle")
        self.source_details_toggle.setCheckable(True)
        self.source_details_toggle.setMaximumWidth(390)
        self.source_details_toggle.toggled.connect(source_box.setVisible)
        self.source_details_toggle.toggled.connect(
            lambda checked: self.source_details_toggle.setText(
                "Observation and original FITS details  ▾" if checked else
                "Observation and original FITS details  ▸"))
        source_box.setVisible(False)
        body.addWidget(self.source_details_toggle, 0, Qt.AlignLeft)
        body.addWidget(source_box)

        # Four scientific objects persist, but only one filter is visible
        # at once. A band never shares its preview, exclusions or QA
        # with another band, and batch runs still visit all four.

        batch = QGroupBox("Multi-band PSF diagnostics")
        self.batch_box = batch
        batch.setObjectName("psf0_four_band_toolbar")
        batch_layout = QVBoxLayout(batch)
        batch_note = QLabel(
            "Independent per-band PSF diagnostics · not adopted")
        batch_note.setToolTip(
            "Runs one diagnostic for each available band. Original FITS "
            "are never modified and results are never adopted automatically.")
        batch_note.setWordWrap(True)
        batch_layout.addWidget(batch_note)
        batch_buttons = QHBoxLayout()
        self.build_all_psfs_btn = QPushButton("Build PSF for available bands")
        self.build_all_psfs_btn.setObjectName("psf0_build_all_bands")
        self.build_all_psfs_btn.setToolTip(
            "Runs one real native FITS/Tractor PSF-0 diagnostic for each "
            "available band, sequentially, in the existing background "
            "worker. Errors are isolated by band. No automatic downloads, "
            "science adoption or FITS export.")
        self.build_all_psfs_btn.clicked.connect(self._request_all_psfs)
        batch_buttons.addWidget(self.build_all_psfs_btn)
        self.copy_full_qa_btn = QPushButton("Copy complete QA")
        self.copy_full_qa_btn.setObjectName("psf0_copy_full_griz_debug")
        self.copy_full_qa_btn.clicked.connect(self.copy_full_qa)
        self.copy_full_qa_btn.setToolTip(
            "Copy full settings, all band outcomes, source evidence, "
            "star IDs, Gaia identity, sky metrics and signed per-star core/"
            "wing residual diagnostics as JSON. No native image pixels.")
        batch_buttons.addWidget(self.copy_full_qa_btn)
        self.save_full_qa_btn = QPushButton("Export QA JSON…")
        self.save_full_qa_btn.setObjectName("psf0_save_full_griz_debug")
        self.save_full_qa_btn.clicked.connect(self.save_full_qa)
        self.save_full_qa_btn.setToolTip(
            "Export the same full four-filter diagnostic JSON as an "
            "explicit user-selected new .json file, for ChatGPT attachment. "
            "Original FITS files are never modified.")
        batch_buttons.addWidget(self.save_full_qa_btn)
        batch_buttons.addStretch(1)
        batch_layout.addLayout(batch_buttons)
        self.batch_status = QLabel("")
        self.batch_status.setObjectName("psf0_griz_batch_status")
        self.batch_status.setWordWrap(True)
        batch_layout.addWidget(self.batch_status)
        self.psf0_panels = {}
        self.band_tabs = QTabWidget(central)
        self.band_tabs.setObjectName("psf0_single_band_tabs")
        self.band_tabs.setDocumentMode(True)
        self.band_tabs.setMinimumWidth(0)
        # The tab's initial three-card horizontal size hint must not
        # impose a hard minimum on the entire desktop window. Reflow is
        # decided by the parent before the tab gets its viewport.
        self.band_tabs.setSizePolicy(
            QSizePolicy.Ignored, QSizePolicy.Maximum)
        for native_band in ("g", "r", "i", "z"):
            panel = PSF0Panel(self.band_tabs, fixed_band=native_band)
            panel.setObjectName("psf0_panel_" + native_band)
            panel.runRequested.connect(self.psf0RunRequested.emit)
            panel.saveRequested.connect(self.psf0SaveRequested.emit)
            panel.autoFwhmRequested.connect(self.psf0AutoFwhmRequested.emit)
            panel.tractorRequested.connect(self.psf0TractorRequested.emit)
            self.psf0_panels[native_band] = panel
            index = self.band_tabs.addTab(panel, native_band.upper())
            panel.overviewChanged.connect(self._update_band_tab_statuses)
            self.band_tabs.setTabToolTip(
                index, "Independent " + native_band + "-band settings, "
                "scientific QA and non-adopted preview")
            self.band_selector.addTab(native_band.upper())
        # The fixed selector above the scroll replaces this inner tab bar;
        # both selectors use the same existing QTabWidget state.
        self.band_tabs.tabBar().hide()
        self.band_selector.currentChanged.connect(
            self.band_tabs.setCurrentIndex)
        self.band_tabs.currentChanged.connect(
            self.band_selector.setCurrentIndex)
        body.addWidget(self.band_tabs)
        body.addWidget(batch)
        # No empty multi-band action block when the selected galaxy has
        # no usable local PSF inputs. It reappears as soon as data exist.
        self.batch_box.setVisible(False)
        self.view_mode.currentIndexChanged.connect(
            self._set_view_mode)
        # Compatibility alias for older callers; resolves to the first
        # available band on bind, not a second shared scientific panel.
        self.psf0_panel = self.psf0_panels["g"]

        self.blocker = QLabel(
            "PSF-0 preview is experimental and not adopted science. "
            "No final background, masks, photometry or SBF is implemented. "
            "File presence alone never establishes source-FITS validity.")
        self.blocker.setObjectName("photometry_blocker")
        self.blocker.setWordWrap(True)
        self.blocker.setVisible(False)
        body.addWidget(self.blocker)
        center_scroll.setWidget(central)
        psf_stage_layout.addWidget(center_scroll, 1)
        self.scene_editor = SceneEditor(self.stage_tabs)
        self.stage_tabs.addTab(self.scene_editor, "2 · Scene and Masks")
        self.stage_tabs.setTabToolTip(
            1, "Original FITS/WCS editor; editable ICRS mask drafts. "
               "No final galaxy segmentation or adopted photometry.")
        self.sky_stage = SkyWorkspace(self.scene_editor, self.stage_tabs)
        self.stage_tabs.addTab(self.sky_stage, "3 · Background")
        self.stage_tabs.setTabToolTip(
            2, "Automatically selected same-project residual-sky controls; "
               "only diagnostic models, no scientific sky adoption.")
        self.aperture_stage = ApertureWorkspace(
            self.scene_editor, self.sky_stage, self.stage_tabs)
        self.stage_tabs.addTab(
            self.aperture_stage, "4 · Aperture Photometry")
        self.stage_tabs.setTabToolTip(
            3, "Exact observed finite-aperture fluxes, independent elliptical "
               "annuli and curve-of-growth diagnostics. No adopted total.")
        root.addWidget(self.stage_tabs, 1)
        self._last_stage_widget = psf_stage
        self.stage_tabs.currentChanged.connect(self._stage_changed)
        self.stage_tabs.currentChanged.connect(self._emit_navigation)
        self.band_tabs.currentChanged.connect(self._emit_navigation)
        self.scene_editor.band.currentIndexChanged.connect(self._emit_navigation)
        self._set_view_mode()
        self._update_band_tab_statuses()
        self._apply_science_reflow()

    def current_navigation(self):
        current = self.stage_tabs.currentWidget()
        stage = ("scene" if current is self.scene_editor
                 else "sky" if current is self.sky_stage
                 else "aperture" if current is self.aperture_stage else "psf0")
        if stage in ("scene", "sky", "aperture"):
            band = self.scene_editor.band.currentData()
        else:
            panel = self.band_tabs.currentWidget()
            band = panel.band.currentData() if panel is not None else None
        return stage, band if band in ("g", "r", "i", "z") else "g"

    def _emit_navigation(self, *_):
        if self._suppress_navigation or self._project_root is None:
            return
        stage, band = self.current_navigation()
        self.navigationChanged.emit(stage, band)

    def restore_navigation(self, *, stage, band):
        """Restore user location only; never run PSF/scene FITS workers."""
        if stage not in ("psf0", "scene", "sky", "aperture") or band not in ("g", "r", "i", "z"):
            raise ValueError("Unknown photometry resume location")
        self._suppress_navigation = True
        try:
            panel = self.psf0_panels.get(band)
            if panel is not None:
                index = self.band_tabs.indexOf(panel)
                if index >= 0:
                    self.band_tabs.setCurrentIndex(index)
            self.scene_editor.set_band(band)
            self.stage_tabs.setCurrentWidget(
                self.scene_editor if stage == "scene"
                else self.sky_stage if stage == "sky"
                else self.aperture_stage if stage == "aperture"
                else self.stage_tabs.widget(0))
            if stage == "sky":
                self.sky_stage.refresh_context()
            if stage == "aperture":
                self.aperture_stage.refresh_context()
        finally:
            self._suppress_navigation = False

    def _stage_changed(self, index):
        current = self.stage_tabs.widget(index)
        previous = self._last_stage_widget
        self._last_stage_widget = current
        psf = self.stage_tabs.widget(0)
        # Synchronize from PSF only when LEAVING the PSF tab. Returning
        # from 03 Background to 02 Scene must preserve the loaded native
        # band, annotations and prepared masks, regardless of PSF tab.
        if previous is psf and current in (
                self.scene_editor, self.sky_stage, self.aperture_stage):
            # Do not invalidate an already loaded scene, its masks, and
            # completed sky/photometry merely because the researcher
            # inspected a different PSF-0 band. An intentional change of
            # the native Scene band uses its own explicit band selector.
            if self.scene_editor._loaded is None:
                band = self.band_tabs.currentWidget().band.currentData()
                if band in ("g", "r", "i", "z"):
                    self.scene_editor.set_band(band)
        if current is self.sky_stage:
            self.sky_stage.refresh_context()
        if current is self.aperture_stage:
            self.aperture_stage.refresh_context()
        # Opening Scene is deliberately read-only: local FITS load needs
        # the explicit Load original FITS button. Mask revision saving
        # also requires a separate explicit researcher action.

    def _apply_science_reflow(self):
        if hasattr(self, "psf0_panels"):
            # Decide from the available workspace viewport BEFORE the
            # inner tab attempts to enforce a horizontal minimum width.
            narrow = self.width() < 1140
            for panel in self.psf0_panels.values():
                panel.set_narrow_visual_layout(narrow)

    def resizeEvent(self, event):
        self._apply_science_reflow()
        super().resizeEvent(event)

    def _set_view_mode(self, *_):
        research = self.view_mode.currentData() == "research"
        for panel in self.psf0_panels.values():
            panel.set_research_mode(research)
        # No speculative multi-stage roadmap or duplicated methodology
        # boxes in the operational PSF workbench. Scientific workflow
        # details live in repository docs and reproducible full QA.
        self.source_details_toggle.setVisible(research)
        self.source_box.setVisible(
            research and self.source_details_toggle.isChecked())
        self.copy_full_qa_btn.setVisible(research)
        # Scientific controls remain inspectable in explicit research mode,
        # including before any local FITS can be used.
        self.batch_box.setVisible(
            research or self.build_all_psfs_btn.isEnabled()
            or bool(self.batch_status.text().strip()))
        # Complete JSON export remains accessible in Guided mode.
        self._update_band_tab_statuses()

    def _update_band_tab_statuses(self):
        for band, panel in self.psf0_panels.items():
            index = self.band_tabs.indexOf(panel)
            if index < 0:
                continue
            state = panel.overview_state()
            label = f"{band.upper()} · {state}"
            tooltip = (
                f"{band}-band: {state}. Independent controls, "
                "science result and source provenance; not adopted.")
            self.band_tabs.setTabText(index, label)
            self.band_tabs.setTabToolTip(index, tooltip)
            self.band_selector.setTabText(index, label)
            self.band_selector.setTabToolTip(index, tooltip)

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

    def copy_project_qa(self):
        """One-click current evidence bundle; never alter stage state."""
        if self._project_root is None or self._bound_project_id is None:
            self.project_qa_status.setText(
                "Open a saved galaxy project before copying project QA.")
            self.project_qa_status.show()
            return
        try:
            content = serialize_project_diagnostic_bundle(self)
        except (TypeError, ValueError, OverflowError, RuntimeError,
                OSError) as exc:
            self.project_qa_status.setText(
                f"Project QA collection failed ({type(exc).__name__}): "
                f"{exc}. Clipboard not changed.")
            self.project_qa_status.show()
            return
        QApplication.clipboard().setText(content)
        self.project_qa_status.setText(
            f"Copied all-stage project QA: {len(content):,} JSON characters. "
            "Paste into ChatGPT or save the clipboard as UTF-8 .json. "
            "No worker, FITS or draft changed. NOT_RUN is not a pass.")
        self.project_qa_status.show()

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
        self._suppress_navigation = True
        self._project_root = Path(saved.root)
        self.scene_editor.bind_project(self._project_root)
        self._bound_project_id = saved.project_id
        self._bound_source_identity = identity
        # A new UUID, release/region/brick or data root invalidates PSFs.
        # Changing only the saved default band does NOT change the four
        # native per-band input datasets, so existing PSF evidence survives.
        try:
            self.refresh(force=not same_project)
        finally:
            self._suppress_navigation = False

    def refresh(self, *, force=False):
        # Only a new project or explicit Refresh can discard a preview.
        # Re-entering the same project/photometry view is non-destructive.
        self.scene_editor.bind_project(self._project_root)
        self.sky_stage.refresh_context()
        self.aperture_stage.refresh_context()
        for panel in self.psf0_panels.values():
            if force and panel._project_root is not None:
                panel._project_root = self._project_root
                panel.refresh_sources()
            else:
                panel.bind_project(self._project_root)
        # Select the first native band with a complete on-disk quartet,
        # even when Gaia/Tractor has not yet been downloaded. Otherwise a
        # project with only R-band FITS lands incorrectly on empty G rather
        # than explaining the missing star catalogue on the usable band.
        self.psf0_panel = next(
            (panel for panel in self.psf0_panels.values()
             if panel._snapshot is not None and any(
                 entry["band"] == panel.band.currentData() and
                 entry["paths_present_unverified"]
                 for entry in panel._snapshot["bands"])),
            self.psf0_panels["g"])
        source_key = (self._bound_project_id, self._bound_source_identity)
        if source_key != self._tab_source_identity:
            # A newly bound project opens on its first usable band.
            # Re-entering an unchanged project preserves the user's
            # currently inspected band and every in-memory PSF result.
            selected_index = self.band_tabs.indexOf(self.psf0_panel)
            self.band_tabs.setCurrentIndex(
                selected_index if selected_index >= 0 else 0)
            self._tab_source_identity = source_key
        self.build_all_psfs_btn.setEnabled(any(
            panel.run_btn.isEnabled() for panel in self.psf0_panels.values()))
        self.batch_box.setVisible(
            self.view_mode.currentData() == "research"
            or self.build_all_psfs_btn.isEnabled()
            or bool(self.batch_status.text().strip()))
        self._update_band_tab_statuses()
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
