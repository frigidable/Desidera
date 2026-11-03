"""Explicit SIMBAD search/review UI: no network on project opening."""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from PySide6.QtCore import Qt, QTimer
from PySide6.QtWidgets import (
    QDialog, QDialogButtonBox, QDoubleSpinBox, QHBoxLayout, QLabel,
    QComboBox,
    QMessageBox, QPushButton, QTableWidget, QTableWidgetItem, QVBoxLayout,
    QAbstractItemView,
)

from lsbpipe.surveys.object_projects import read_object_project
from .catalogue_links import confirmed_name_link
from lsbpipe.surveys.galaxy_knowledge import (
    query_simbad_candidates, store_simbad_snapshot,
    latest_simbad_snapshot, latest_simbad_confirmation,
    read_simbad_snapshot, confirm_simbad_identity,
)

# Bounded shared worker pool: opening or closing a dialog never owns/destroys
# a running QThread; obsolete completions are never applied to another project.
_POOL = ThreadPoolExecutor(max_workers=2, thread_name_prefix="desidera-simbad")


class GalaxyKnowledgeDialog(QDialog):
    def __init__(self, project_root: str | Path, parent=None):
        super().__init__(parent)
        self.project_root = Path(project_root)
        self.project = read_object_project(self.project_root)
        self.setWindowTitle(f"Galaxy Knowledge · {self.project.display_name}")
        self.resize(850, 570)
        self._snapshot_digest = None
        self._future = None

        root = QVBoxLayout(self)
        self.intro = QLabel(
            f"SIMBAD identity candidates for {self.project.display_name}. "
            "Coordinates and names are unverified catalogue information. "
            "Choosing a candidate never changes the project's ICRS position.")
        self.intro.setWordWrap(True)
        root.addWidget(self.intro)

        controls = QHBoxLayout()
        controls.addWidget(QLabel("Search radius"))
        self.radius = QDoubleSpinBox()
        self.radius.setObjectName("knowledge_radius_arcsec")
        self.radius.setRange(1, 300)
        self.radius.setValue(60)
        self.radius.setDecimals(1)
        self.radius.setSuffix(" arcsec")
        self.radius.setEnabled(self.project.field is not None)
        controls.addWidget(self.radius)
        controls.addWidget(QLabel("SIMBAD server"))
        self.endpoint_policy = QComboBox()
        self.endpoint_policy.setObjectName("knowledge_endpoint_policy")
        self.endpoint_policy.addItem("Automatic · Strasbourg → Harvard", "auto")
        self.endpoint_policy.addItem("Automatic · Harvard → Strasbourg", "prefer_harvard")
        self.endpoint_policy.addItem("Strasbourg only", "strasbourg_only")
        self.endpoint_policy.addItem("Harvard only", "harvard_only")
        self.endpoint_policy.setToolTip(
            "Only the two official SIMBAD/astroquery mirrors. "
            "Failover occurs on network/service errors, not empty searches.")
        controls.addWidget(self.endpoint_policy)
        controls.addStretch(1)
        self.search_btn = QPushButton("Search SIMBAD…")
        self.search_btn.setObjectName("knowledge_search_simbad")
        self.search_btn.clicked.connect(self.search)
        controls.addWidget(self.search_btn)
        root.addLayout(controls)

        self.status = QLabel("Offline · no remote query has been made.")
        self.status.setObjectName("knowledge_search_status")
        self.status.setWordWrap(True)
        root.addWidget(self.status)
        self.endpoint_note = QLabel("Server provenance: no query yet.")
        self.endpoint_note.setObjectName("knowledge_endpoint_provenance")
        self.endpoint_note.setWordWrap(True)
        root.addWidget(self.endpoint_note)

        self.table = QTableWidget(0, 5)
        self.table.setObjectName("knowledge_candidate_table")
        self.table.setHorizontalHeaderLabels(
            ["SIMBAD name", "Separation (″)", "RA (ICRS)", "Dec (ICRS)", "Found via"])
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.itemSelectionChanged.connect(self._selection_changed)
        root.addWidget(self.table, 1)

        self.confirmed = QLabel("No confirmed SIMBAD catalogue association.")
        self.confirmed.setObjectName("knowledge_confirmed")
        self.confirmed.setWordWrap(True)
        self.confirmed.setTextFormat(Qt.RichText)
        self.confirmed.setTextInteractionFlags(
            Qt.TextBrowserInteraction | Qt.TextSelectableByMouse)
        self.confirmed.setOpenExternalLinks(True)
        root.addWidget(self.confirmed)

        lower = QHBoxLayout()
        self.confirm_btn = QPushButton("Confirm selected identity…")
        self.confirm_btn.setObjectName("knowledge_confirm_candidate")
        self.confirm_btn.clicked.connect(self._confirm)
        self.confirm_btn.setEnabled(False)
        lower.addWidget(self.confirm_btn)
        lower.addStretch(1)
        close = QDialogButtonBox(QDialogButtonBox.Close)
        close.rejected.connect(self.reject)
        lower.addWidget(close)
        root.addLayout(lower)

        self._poll = QTimer(self)
        self._poll.setInterval(100)
        self._poll.timeout.connect(self._check_result)
        self._load_offline()

    def _load_offline(self):
        self.confirmed.setText("No confirmed SIMBAD catalogue association.")
        try:
            existing = latest_simbad_snapshot(self.project_root)
            if existing:
                self._show_snapshot(*existing)
                self.status.setText(
                    "Stored SIMBAD snapshot · offline. Refresh explicitly to query again.")
            receipt = latest_simbad_confirmation(self.project_root)
            if receipt:
                c = receipt["candidate"]
                source = read_simbad_snapshot(
                    self.project_root, receipt["snapshot_sha256"])
                endpoint = source.get("endpoint_evidence") or {}
                host = (endpoint.get("endpoint_host")
                        if isinstance(endpoint, dict) else None)
                link = confirmed_name_link(
                    "SIMBAD", c["main_id"], endpoint_host=host)
                self.confirmed.setText(
                    "User-confirmed SIMBAD identity: " + link + ". "
                    "External reference only; ICRS target and science results unchanged.")
        except (OSError, ValueError, KeyError, TypeError) as exc:
            self.status.setText(f"Stored evidence cannot be trusted: {exc}")
            self.search_btn.setEnabled(False)
            self.confirm_btn.setEnabled(False)

    def _show_snapshot(self, digest: str, data: dict):
        self._snapshot_digest = digest
        provenance = data.get("endpoint_evidence") or {}
        host = provenance.get("endpoint_host")
        if host:
            failed = provenance.get("failed_attempts") or []
            attempted = provenance.get("attempted_hosts") or []
            self.endpoint_note.setText(
                f"SIMBAD server: {host} · "
                + ("fallback after " + ", ".join(
                    f["host"] for f in failed) if failed else "primary/selected endpoint")
                + (f" · tried {len(attempted)} endpoint(s)" if attempted else "")
                + " · external identity unverified")
        else:
            self.endpoint_note.setText(
                "SIMBAD server: not recorded in this older snapshot")
        candidates = data["candidates"]
        self.table.setRowCount(0)
        for c in candidates:
            row = self.table.rowCount()
            self.table.insertRow(row)
            values = (
                c["main_id"],
                "—" if c["separation_arcsec"] is None else
                f"{c['separation_arcsec']:.2f}",
                f"{c['ra_deg']:.7f}", f"{c['dec_deg']:+.7f}",
                ", ".join(c["discovered_via"]))
            for column, value in enumerate(values):
                self.table.setItem(row, column, QTableWidgetItem(value))
            self.table.item(row, 0).setData(Qt.UserRole, c["main_id"])
        self.table.clearSelection()
        self.confirm_btn.setEnabled(False)
        if data.get("truncated"):
            self.status.setText(
                f"{len(candidates)} candidates shown · list truncated. "
                "No automatic identification.")
        elif not candidates:
            self.status.setText(
                "SIMBAD returned no candidates. Project identity remains unresolved.")
        else:
            self.status.setText(
                f"{len(candidates)} candidate(s) found. Select and confirm "
                "only after checking the object's identity.")

    def _selection_changed(self):
        self.confirm_btn.setEnabled(
            self._future is None and self._snapshot_digest is not None
            and len(self.table.selectedItems()) > 0)

    def search(self):
        if self._future is not None:
            return
        project = self.project
        field = project.field
        self.search_btn.setEnabled(False)
        self.endpoint_policy.setEnabled(False)
        self.confirm_btn.setEnabled(False)
        self.status.setText("Contacting SIMBAD · no local scientific data will change…")
        self._future = _POOL.submit(
            query_simbad_candidates, name=project.display_name,
            ra_deg=field.ra_deg if field else None,
            dec_deg=field.dec_deg if field else None,
            radius_arcsec=self.radius.value(),
            mirror_policy=self.endpoint_policy.currentData())
        self._poll.start()

    def _check_result(self):
        if self._future is None or not self._future.done():
            return
        future, self._future = self._future, None
        self._poll.stop()
        self.search_btn.setEnabled(True)
        self.endpoint_policy.setEnabled(True)
        try:
            response = future.result()
            digest = store_simbad_snapshot(self.project_root, response)
            self._show_snapshot(digest, response)
        except Exception as exc:
            self.status.setText(
                f"SIMBAD request or evidence storage failed: {type(exc).__name__}: {exc}")
            # Do not let a failed new query silently validate a previous hit.
            self._snapshot_digest = None
            self.confirm_btn.setEnabled(False)

    def _confirm(self):
        items = self.table.selectedItems()
        if not items or self._snapshot_digest is None:
            return
        row = self.table.currentRow()
        if row < 0 or self.table.item(row, 0) is None:
            return
        name = self.table.item(row, 0).data(Qt.UserRole)
        answer = QMessageBox.question(
            self, "Confirm catalogue association",
            f"Associate this project with SIMBAD object '{name}'?\n\n"
            "This confirms a catalogue identity only. It will not change "
            "the project's coordinates, data, or science results.",
            QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if answer != QMessageBox.Yes:
            return
        try:
            confirm_simbad_identity(self.project_root, self._snapshot_digest, name)
            self._load_offline()
        except (OSError, ValueError, KeyError, TypeError) as exc:
            QMessageBox.warning(self, "Cannot confirm SIMBAD identity", str(exc))
