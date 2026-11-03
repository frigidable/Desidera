"""P0.8a: human-reviewed, immutable setup snapshots in Project Hub."""
from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView, QDialog, QDialogButtonBox, QHBoxLayout,
    QInputDialog, QLabel, QMessageBox, QPushButton, QTableWidget,
    QTableWidgetItem, QVBoxLayout,
)

from lsbpipe.surveys.object_projects import read_object_project
from lsbpipe.surveys.analysis_configurations import (
    list_analysis_configurations, freeze_analysis_configuration,
    compare_current_setup,
)


class AnalysisConfigurationsDialog(QDialog):
    def __init__(self, project_root: str | Path, parent=None):
        super().__init__(parent)
        self.project_root = Path(project_root)
        project = read_object_project(project_root)
        self.setWindowTitle(f"Analysis configurations · {project.display_name}")
        self.resize(870, 480)
        root = QVBoxLayout(self)
        intro = QLabel(
            "Immutable snapshots of target, dataset selection and recorded "
            "original-FITS identities. Freezing neither reruns QA nor verifies "
            "today's FITS bytes, and never adopts scientific measurements.")
        intro.setWordWrap(True)
        root.addWidget(intro)
        self.table = QTableWidget(0, 4)
        self.table.setObjectName("analysis_config_table")
        self.table.setHorizontalHeaderLabels([
            "Label", "Created (UTC)", "Input evidence", "SHA-256 (short)"])
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.itemSelectionChanged.connect(self._selected)
        root.addWidget(self.table, 1)
        self.status = QLabel("No setup selected.")
        self.status.setObjectName("analysis_config_status")
        self.status.setWordWrap(True)
        root.addWidget(self.status)
        row = QHBoxLayout()
        self.freeze_btn = QPushButton("Freeze current setup…")
        self.freeze_btn.setObjectName("analysis_config_freeze")
        self.freeze_btn.setToolTip(
            "Save a new immutable setup. No existing configuration, "
            "project, FITS or science output is overwritten.")
        self.freeze_btn.clicked.connect(self._freeze)
        row.addWidget(self.freeze_btn)
        row.addStretch(1)
        close = QDialogButtonBox(QDialogButtonBox.Close)
        close.rejected.connect(self.reject)
        row.addWidget(close)
        root.addLayout(row)
        self.refresh()

    def refresh(self):
        try:
            items = list_analysis_configurations(self.project_root)
        except (OSError, ValueError, KeyError, TypeError) as exc:
            self.table.setRowCount(0)
            self.status.setText(f"Invalid configuration evidence: {exc}")
            self.freeze_btn.setEnabled(False)
            return
        self.table.setRowCount(0)
        for config in items:
            row = self.table.rowCount()
            self.table.insertRow(row)
            values = (
                config["label"], config["created_at_utc"],
                config["input_evidence_state"], config["digest"][:12])
            for column, value in enumerate(values):
                item = QTableWidgetItem(value)
                item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                self.table.setItem(row, column, item)
            self.table.item(row, 0).setData(Qt.UserRole, config["digest"])
        if not items:
            self.status.setText(
                "No frozen configurations · the current project remains editable.")
        else:
            self.table.selectRow(0)
            self._selected()

    def _selected(self):
        row = self.table.currentRow()
        if row < 0 or self.table.item(row, 0) is None:
            return
        digest = self.table.item(row, 0).data(Qt.UserRole)
        try:
            report = compare_current_setup(self.project_root, digest)
            self.status.setText(
                f"Snapshot {digest[:12]} · {report['state']}. "
                "Original FITS not rehashed; this is not a completed science run.")
        except (OSError, ValueError, KeyError, TypeError) as exc:
            self.status.setText(f"Cannot verify snapshot: {exc}")

    def _freeze(self):
        label, confirmed = QInputDialog.getText(
            self, "Freeze analysis setup", "Configuration label:")
        if not confirmed:
            return
        try:
            digest, _ = freeze_analysis_configuration(self.project_root, label)
            self.refresh()
            self.status.setText(
                f"Saved immutable configuration {digest[:12]}. "
                "Run not executed; original FITS not rehashed.")
        except (OSError, ValueError, KeyError, TypeError) as exc:
            QMessageBox.warning(
                self, "Cannot freeze analysis configuration", str(exc))
