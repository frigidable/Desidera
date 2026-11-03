"""Project-bound Photometry Workspace: preparation and provenance only.

No photometric algorithms are exposed before their validation gates exist.
"""
from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QAbstractItemView, QGroupBox, QHBoxLayout, QHeaderView, QLabel,
    QPushButton, QScrollArea, QSplitter, QTableWidget, QTableWidgetItem,
    QVBoxLayout, QWidget,
)

from lsbpipe.photometry.readiness import inspect_photometry_inputs
from lsbpipe.gui.psf0_panel import PSF0Panel

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

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("photometry_workspace")
        self._project_root = None
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
        self.refresh_btn.clicked.connect(self.refresh)
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

        # Keep the first genuinely runnable scientific stage above the
        # long original-FITS inventory: PSF controls remain discoverable
        # without an initial scroll on laptop-sized Windows viewports.
        self.psf0_panel = PSF0Panel(central)
        self.psf0_panel.runRequested.connect(self.psf0RunRequested.emit)
        self.psf0_panel.saveRequested.connect(self.psf0SaveRequested.emit)
        self.psf0_panel.autoFwhmRequested.connect(
            self.psf0AutoFwhmRequested.emit)
        body.addWidget(self.psf0_panel)
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

    def bind_project(self, project):
        """Read the immutable UUID from the existing project manifest."""
        from lsbpipe.surveys.object_projects import read_object_project
        saved = read_object_project(project.root)
        if saved.project_id != project.project_id:
            raise ValueError("Photometry workspace project identity mismatch")
        self._project_root = Path(saved.root)
        self.refresh()

    def refresh(self):
        # Always clear/rebind the scientific preview first: even a corrupt
        # project manifest must never leave a previous galaxy's PSF active.
        self.psf0_panel.bind_project(self._project_root)
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
