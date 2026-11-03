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

_STAGES = (
    ("01", "Observations", "Native quartet and observational provenance"),
    ("02", "Scene & masks", "Target protection, contaminants and valid pixels"),
    ("03", "Background", "Residual sky models and identifiability"),
    ("04", "PSF & noise", "Wings, resampling and covariance"),
    ("05", "Isophotes", "Fixed/free geometry, radial profiles and growth"),
    ("06", "Structural models", "Forward-modelled Sérsic and alternatives"),
    ("07", "Validation", "Uncertainties and injection-recovery"),
    ("08", "Results", "Immutable measurements, figures and QA"),
)


class PhotometryWorkspace(QWidget):
    """Dedicated screen; never presents input-file presence as science-ready."""

    backRequested = Signal()
    sourcesRequested = Signal()

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
            "Scientific preparation · no photometric measurements have been made")
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
            caption.setStyleSheet("color: palette(mid);")
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
        body.addWidget(source_box)

        method_box = QGroupBox("Measurement design")
        method_layout = QVBoxLayout(method_box)
        for heading, detail in (
            ("Coadd preparation",
             "Original pixels, mask semantics, WCS and units; no implicit "
             "background subtraction, interpolation or mosaic."),
            ("Non-parametric branch",
             "Masked elliptical aperture sums, free/fixed isophotes, "
             "growth curves and multi-band profiles."),
            ("Parametric branch",
             "PSF-convolved scene model with galaxy, neighbor and sky "
             "parameters; distinguish measured from extrapolated flux."),
            ("Hybrid extension",
             "Use coadd for initialization and then jointly fit distinct "
             "CCD exposures with explicit calibration and lineage."),
            ("Result gate",
             "Covariance-aware uncertainty, sky-model sensitivity, recovery "
             "tests and run provenance before adoption."),
        ):
            line = QLabel(f"<b>{heading}</b> — {detail}")
            line.setTextFormat(Qt.RichText)
            line.setWordWrap(True)
            line.setMargin(5)
            method_layout.addWidget(line)
        body.addWidget(method_box)

        self.blocker = QLabel(
            "No science run is available on this screen yet. "
            "Input presence alone never enables photometry.")
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
            "Science state: NOT_MEASURED_NOT_ADOPTED. "
            "The screen checks saved file paths only. FITS structure, WCS, "
            "pixel masks, sky, PSF, correlated noise and photometry have "
            "not been validated here. No science or catalogue value is adopted.")
