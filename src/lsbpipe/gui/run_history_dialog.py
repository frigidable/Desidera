"""Read-only presentation of immutable project run history."""
from __future__ import annotations

from pathlib import Path
from PySide6.QtCore import Qt, QUrl
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QLabel, QTableWidget, QTableWidgetItem,
    QAbstractItemView, QPlainTextEdit, QHBoxLayout, QPushButton,
)


class RunHistoryDialog(QDialog):
    """No science is executed or adopted by opening this window."""

    def __init__(self, records: list[dict], parent=None):
        super().__init__(parent)
        self.setWindowTitle("DESIDERA · Scientific run history")
        self.setObjectName("science_run_history_dialog")
        self.resize(940, 580)
        self.records = list(records)
        outer = QVBoxLayout(self)
        description = QLabel(
            "Completed diagnostics only · Old and un-rehashed records are NOT "
            "certified current. A saved run is not adopted photometry.")
        description.setWordWrap(True)
        outer.addWidget(description)
        self.table = QTableWidget(0, 5)
        self.table.setObjectName("science_run_history_table")
        self.table.setHorizontalHeaderLabels(
            ["Recorded run", "Action", "Dataset", "Review state", "Created (UTC)"])
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.verticalHeader().hide()
        self.table.setAlternatingRowColors(True)
        outer.addWidget(self.table, 2)
        self.details = QPlainTextEdit()
        self.details.setObjectName("science_run_history_detail")
        self.details.setReadOnly(True)
        outer.addWidget(self.details, 1)
        actions = QHBoxLayout()
        self.open_manifest = QPushButton("Open selected provenance JSON")
        self.open_manifest.setObjectName("science_run_open_receipt")
        self.open_manifest.setEnabled(False)
        self.open_manifest.clicked.connect(self._open_selected)
        close = QPushButton("Close")
        close.clicked.connect(self.accept)
        actions.addWidget(self.open_manifest)
        actions.addStretch()
        actions.addWidget(close)
        outer.addLayout(actions)
        for record in self.records:
            row = self.table.rowCount()
            self.table.insertRow(row)
            fields = [
                record.get("run_id") or "(unknown)",
                record.get("action") or "(unknown)",
                record.get("dataset") or "(unknown)",
                record.get("state") or "(unknown)",
                record.get("created_at_utc") or "(unknown)",
            ]
            for col, item in enumerate(fields):
                cell = QTableWidgetItem(str(item))
                cell.setFlags(cell.flags() & ~Qt.ItemIsEditable)
                self.table.setItem(row, col, cell)
        self.table.setColumnWidth(0, 250)
        self.table.setColumnWidth(1, 90)
        self.table.setColumnWidth(2, 215)
        self.table.setColumnWidth(3, 220)
        self.table.itemSelectionChanged.connect(self._selection_changed)
        if self.records:
            self.table.selectRow(0)
        else:
            self.details.setPlainText(
                "No completed QA receipts found in this project's folders.")

    def _selection_changed(self):
        row = self.table.currentRow()
        if row < 0 or row >= len(self.records):
            self.open_manifest.setEnabled(False)
            return
        entry = self.records[row]
        manifest = Path(entry["manifest_path"])
        self.open_manifest.setEnabled(manifest.is_file() and not manifest.is_symlink())
        self.details.setPlainText(
            f"Review status: {entry.get('state')}\n"
            f"Reason: {entry.get('reason', 'No additional reason')}\n"
            f"Manifest: {manifest}\n\n"
            "Recorded input hashes are historical until the original FITS "
            "have been rehashed for comparison. Never infer photometric "
            "readiness from a completed diagnostic.")

    def _open_selected(self):
        row = self.table.currentRow()
        if 0 <= row < len(self.records):
            path = Path(self.records[row]["manifest_path"])
            if path.is_file() and not path.is_symlink():
                QDesktopServices.openUrl(QUrl.fromLocalFile(str(path)))
