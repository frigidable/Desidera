"""NED ObjectLookup: explicit online query and offline immutable identity review."""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from PySide6.QtCore import Qt, QTimer
from PySide6.QtWidgets import (
    QDialog, QDialogButtonBox, QFormLayout, QHBoxLayout, QLabel,
    QLineEdit, QMessageBox, QPushButton, QVBoxLayout,
)

from lsbpipe.surveys.object_projects import read_object_project
from lsbpipe.surveys.ned_knowledge import (
    query_ned_identity, save_ned_snapshot, latest_ned_snapshot,
    confirm_ned_identity, latest_ned_confirmation,
)

_POOL = ThreadPoolExecutor(max_workers=2, thread_name_prefix="desidera-ned")


class NedKnowledgeDialog(QDialog):
    """Never queries NED on construction; all saved review works offline."""

    def __init__(self, project_root: str | Path, parent=None):
        super().__init__(parent)
        self.project_root = Path(project_root)
        self.project = read_object_project(project_root)
        self.setWindowTitle(f"NED Knowledge · {self.project.display_name}")
        self.resize(790, 510)
        self._future = None
        self._digest = None
        self._payload = None

        outer = QVBoxLayout(self)
        intro = QLabel(
            "NASA/IPAC NED ObjectLookup resolves a catalogue name, but a "
            "name match alone does not establish the astronomical identity. "
            "Review coordinates, object type, and source references before "
            "confirming. External data never overwrite science results.")
        intro.setWordWrap(True)
        outer.addWidget(intro)

        form = QFormLayout()
        self.name_input = QLineEdit(self.project.display_name)
        self.name_input.setObjectName("ned_query_name")
        form.addRow("Query name", self.name_input)
        outer.addLayout(form)

        toolbar = QHBoxLayout()
        self.search_btn = QPushButton("Search NED ObjectLookup…")
        self.search_btn.setObjectName("ned_manual_lookup")
        self.search_btn.clicked.connect(self.search)
        toolbar.addWidget(self.search_btn)
        toolbar.addStretch(1)
        outer.addLayout(toolbar)

        self.status = QLabel("Offline · no NED query started.")
        self.status.setObjectName("ned_lookup_status")
        self.status.setWordWrap(True)
        outer.addWidget(self.status)
        self.response = QLabel("No cached NED result.")
        self.response.setObjectName("ned_lookup_response")
        self.response.setWordWrap(True)
        self.response.setTextInteractionFlags(Qt.TextSelectableByMouse)
        outer.addWidget(self.response, 1)
        self.confirmed = QLabel("No user-confirmed NED association.")
        self.confirmed.setObjectName("ned_confirmed")
        self.confirmed.setWordWrap(True)
        outer.addWidget(self.confirmed)

        actions = QHBoxLayout()
        self.confirm_btn = QPushButton("Confirm this NED identity…")
        self.confirm_btn.setObjectName("ned_confirm_identity")
        self.confirm_btn.setEnabled(False)
        self.confirm_btn.clicked.connect(self._confirm)
        actions.addWidget(self.confirm_btn)
        actions.addStretch(1)
        close = QDialogButtonBox(QDialogButtonBox.Close)
        close.rejected.connect(self.reject)
        actions.addWidget(close)
        outer.addLayout(actions)

        self.poll = QTimer(self)
        self.poll.setInterval(100)
        self.poll.timeout.connect(self._check_result)
        self._load_offline()

    def _load_offline(self):
        try:
            query = latest_ned_snapshot(self.project_root)
            if query:
                self._show_result(*query)
                self.status.setText(
                    "Saved NED lookup · offline. Search again to refresh.")
            confirmed = latest_ned_confirmation(self.project_root)
            if confirmed:
                self.confirmed.setText(
                    "User-confirmed NED association: "
                    + confirmed["candidate"]["main_id"]
                    + " · project ICRS and science values unchanged.")
        except (ValueError, OSError, KeyError, TypeError) as exc:
            self.status.setText(f"Saved NED evidence cannot be trusted: {exc}")
            self.search_btn.setEnabled(False)
            self.confirm_btn.setEnabled(False)

    def _show_result(self, digest, result):
        self._digest = digest
        self._payload = result
        status = result["status"]
        item = result.get("candidate")
        self.confirm_btn.setEnabled(status == "KNOWN_OBJECT" and item is not None)
        if item:
            z = item["redshift"]["value"]
            zref = item["redshift"]["refcode"] or "no source refcode"
            details = [
                f"NED name: {item['main_id']}",
                f"Object type: {item.get('object_type') or 'unclassified'}",
                f"NED J2000: RA {item['ra_deg']:.8f}°, "
                f"Dec {item['dec_deg']:+.8f}°",
                "Coordinate source: "
                + (item.get("position_refcode") or "not supplied"),
                "Heliocentric redshift: "
                + (f"{z:.7g} (source: {zref})" if z is not None else "unavailable"),
                "Aliases (not adopted): " +
                ", ".join(result.get("aliases_unadopted", [])[:12]),
                "NED API version: " + result["api_version"],
                "Catalogue values are external, unadopted evidence.",
            ]
            self.response.setText("\n".join(details))
        elif status == "AMBIGUOUS_NAME":
            options = result.get("ambiguous_names", [])
            self.response.setText(
                "NED reports an ambiguous name. Enter a more specific "
                "identifier and search again. No confirmation is possible.\n"
                + "\n".join(options[:25]))
        else:
            self.response.setText(
                f"NED returned {status}. No known object association "
                "is available for confirmation.")
        self.status.setText(f"NED ObjectLookup · {status}")

    def search(self):
        if self._future is not None:
            return
        name = self.name_input.text().strip()
        if not name:
            self.status.setText("Enter a valid NED search name.")
            return
        self._future = _POOL.submit(query_ned_identity, name)
        self.search_btn.setEnabled(False)
        self.name_input.setEnabled(False)
        self.confirm_btn.setEnabled(False)
        self.status.setText(
            "Contacting NED manually · existing evidence remains unchanged.")
        self.poll.start()

    def _check_result(self):
        if self._future is None or not self._future.done():
            return
        task, self._future = self._future, None
        self.poll.stop()
        self.search_btn.setEnabled(True)
        self.name_input.setEnabled(True)
        try:
            result = task.result()
            digest = save_ned_snapshot(self.project_root, result)
            self._show_result(digest, result)
        except Exception as exc:
            self._digest = None
            self._payload = None
            self.confirm_btn.setEnabled(False)
            self.status.setText(
                f"NED query/evidence failure: {type(exc).__name__}: {exc}")

    def _confirm(self):
        if self._digest is None or self._payload is None:
            return
        item = self._payload.get("candidate")
        if not item or self._payload["status"] != "KNOWN_OBJECT":
            return
        decision = QMessageBox.question(
            self, "Confirm NED catalogue identity",
            f"Associate this project with NED object '{item['main_id']}'?\n\n"
            "This confirms only the catalogue identity. NED redshifts, "
            "other measurements, object coordinates, and aliases will not "
            "overwrite project data.",
            QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if decision != QMessageBox.Yes:
            return
        try:
            confirm_ned_identity(self.project_root, self._digest)
            self._load_offline()
        except (OSError, ValueError, TypeError, KeyError) as exc:
            QMessageBox.warning(self, "Cannot confirm NED identity", str(exc))
