"""GUI-only survey product table, filtering and download selection.

The controller never changes source FITS or executes HTTP calls; it handles
only the manifest-based selection presentation and preserves hidden rows.
SurveyExplorer retains legacy public methods as thin delegates.
"""
from __future__ import annotations

from PySide6.QtCore import QObject, Qt
from PySide6.QtWidgets import (
    QCheckBox, QHBoxLayout, QInputDialog, QMessageBox,
    QTableWidgetItem, QWidget)
from lsbpipe.surveys.selection import matching_product, summarize_selection


class SurveyDataController(QObject):
    def __init__(self, owner):
        super().__init__(owner)
        self.owner = owner

    @staticmethod
    def short_http_evidence(evidence):
        return {
            "HTTP_OK": "OK", "HTTP_MISSING": "404",
            "HTTP_ERROR": "ERR", "HTTP_UNVERIFIED": "?",
            "UNVERIFIED_URL": "?",
        }.get(evidence, "?")

    def _fill_brick_filters(self):
        w = self.owner
        previous = w.brick_filter.currentData()
        w.brick_filter.blockSignals(True)
        w.brick_filter.clear()
        w.brick_filter.addItem("All bricks", None)
        for release, region, brick in sorted({
                (hit.release, hit.region, hit.name) for hit in w.hits}):
            w.brick_filter.addItem(
                f"{release.upper()} · {region} · {brick}",
                (release, region, brick))
        index = next((
            i for i in range(w.brick_filter.count())
            if w.brick_filter.itemData(i) is not None and
            tuple(w.brick_filter.itemData(i)) == tuple(previous or ())
        ), 0)
        w.brick_filter.setCurrentIndex(index)
        w.brick_filter.blockSignals(False)

    def _fill_files(self, *, checked_paths=None):
        """Compact row identity, full provenance on hover, persistent selection."""
        w = self.owner
        entries = w.manifest["products"] if w.manifest else []
        w.files.setRowCount(len(entries))
        for row, entry in enumerate(entries):
            check = QCheckBox()
            check.setChecked(
                checked_paths is not None
                and entry["relative_path"] in checked_paths)
            check.toggled.connect(w._update_selection_summary)
            holder = QWidget()
            layout = QHBoxLayout(holder)
            layout.setContentsMargins(1, 0, 1, 0)
            layout.addWidget(check, alignment=Qt.AlignCenter)
            w.files.setCellWidget(row, 0, holder)
            suffix = f"{entry['kind']}-{entry['band']}" if entry["band"] else entry["kind"]
            name = (
                f"{entry['release'].upper()} {entry['region']} "
                f"{entry['brick']} · {suffix}")
            detail = entry["relative_path"] + chr(10) + entry["url"]
            values = (
                (1, name),
                (2, w._short_http_evidence(entry["evidence"])),
                (3, "—"))
            for col, value in values:
                item = QTableWidgetItem(str(value))
                item.setToolTip(detail + chr(10) + f"Status: {entry['evidence']}")
                w.files.setItem(row, col, item)
        w._filter_products()
        w._mark_local_files()

    def _mark_local_files(self):
        """Annotate local evidence separately from remote HTTP verification."""
        w = self.owner
        project = w._local_project
        if project is None or w.manifest is None:
            return
        for row, entry in enumerate(w.manifest["products"]):
            path = entry["relative_path"]
            local = project.data_root / path
            if path in project.on_disk and local.is_file():
                mib = local.stat().st_size / 1048576.
                w.files.item(row, 3).setText(f"{mib:.0f}")
                w.files.item(row, 3).setToolTip(
                    f"LOCAL PRESENT · {local} · {mib:.2f} MiB. "
                    "Existence is not a science/remote integrity check.")
                name = w.files.item(row, 1)
                name.setToolTip(name.toolTip() + "\nLOCAL PRESENT (unverified)")
            else:
                name = w.files.item(row, 1)
                name.setToolTip(name.toolTip() + "\nNOT FOUND in active local project")
        complete = sum(project.complete.values())
        w.project_status.setText(
            f"LOCAL · {project.root.name} · {len(project.on_disk)} original "
            f"products on disk · {complete}/{len(project.complete)} "
            "complete band quartets (not scientifically verified)")

    def _is_visible_product(self, entry):
        w = self.owner
        dataset = w.brick_filter.currentData()
        if dataset is not None and (
                entry["release"], entry["region"], entry["brick"]) != tuple(dataset):
            return False
        return matching_product(
            entry, brick="All bricks",
            band=w.band_filter.currentText(),
            kind=w.kind_filter.currentText())

    def _filter_products(self, *_):
        w = self.owner
        if w.manifest is None:
            return
        for row, entry in enumerate(w.manifest["products"]):
            w.files.setRowHidden(row, not w._is_visible_product(entry))
        w._update_selection_summary()

    def _set_visible_selection(self, checked: bool):
        w = self.owner
        if w.manifest is None:
            return
        for row in range(w.files.rowCount()):
            if not w.files.isRowHidden(row):
                w.files.cellWidget(
                    row, 0).findChild(QCheckBox).setChecked(checked)
        w._update_selection_summary()

    def select_science_quartet(self):
        """Prepare *exactly* four primary science products for a brick and band.

        Explicitly clears the whole manifest selection, not just visible
        rows, so choosing one science bundle cannot accidentally leave 80
        unrelated files selected.
        """
        w = self.owner
        if w.manifest is None:
            return
        row = w.bricks.currentRow()
        if not (0 <= row < len(w.hits)):
            QMessageBox.information(
                w, "Select a brick",
                "Select one release/region brick row above first.")
            return
        hit = w.hits[row]
        bands = tuple(hit.median_exposures)
        if not bands:
            QMessageBox.warning(w, "No band", "No optical bands listed.")
            return
        default = bands.index("r") if "r" in bands else 0
        band, ok = QInputDialog.getItem(
            w, "Science quartet", "Band:", list(bands), default, False)
        if not ok:
            return
        needed = {
            (hit.release, hit.region, hit.name, "image", band),
            (hit.release, hit.region, hit.name, "invvar", band),
            (hit.release, hit.region, hit.name, "nexp", band),
            (hit.release, hit.region, hit.name, "maskbits", None)}
        available = {
            (e["release"],e["region"],e["brick"],e["kind"],e["band"])
            for e in w.manifest["products"]}
        if not needed.issubset(available):
            QMessageBox.warning(
                w, "Incomplete manifest",
                "Not all four products are in the current download profile.")
            return
        for i, entry in enumerate(w.manifest["products"]):
            key = (entry["release"],entry["region"],entry["brick"],
                   entry["kind"],entry["band"])
            w.files.cellWidget(
                i,0).findChild(QCheckBox).setChecked(key in needed)
        w._update_selection_summary()
        w.message.setText(
            f"Science quartet selected: {hit.release.upper()} "
            f"{hit.region} {hit.name} {band}-band. "
            "Verify URLs before downloading.")

    def _clear_all_selection(self):
        w = self.owner
        """Explicitly deselect visible AND hidden products."""
        if w.manifest is None:
            return
        for row in range(w.files.rowCount()):
            w.files.cellWidget(
                row, 0).findChild(QCheckBox).setChecked(False)
        w._update_selection_summary()
        w.message.setText(
            "All products deselected, including hidden rows. "
            "Use Select visible or Science quartet to choose download files.")

    def _selected(self):
        w = self.owner
        if not w.manifest:
            return set()
        return {
            w.manifest["products"][row]["relative_path"]
            for row in range(w.files.rowCount())
            if w.files.cellWidget(row, 0).findChild(QCheckBox).isChecked()
        }

    def _update_selection_summary(self, *_):
        w = self.owner
        if not w.manifest:
            w.selected_label.setText("No files selected")
            w.download_btn.setEnabled(False)
            return
        selected = w._selected()
        report = summarize_selection(w.manifest["products"], selected)
        visible = sum(not w.files.isRowHidden(row)
                      for row in range(w.files.rowCount()))
        sizes = (f"{report.known_mib:.1f} MiB confirmed + unknown"
                 if report.unknown_size_files else
                 f"{report.known_mib:.1f} MiB")
        selection = (
            f"{report.selected_files} selected · " + sizes
            if report.selected_files else
            "0 selected · choose a quartet or individual FITS")
        w.selected_label.setText(
            selection + f" · {visible} candidate files shown (URLs unverified)")
        # Never preselect an entire survey product inventory for a
        # potentially enormous download: user explicitly chooses a quartet
        # or specific original FITS products.
        w.download_btn.setEnabled(bool(selected) and w.thread is None)
        w.selected_label.setToolTip(
            f"{visible} products currently visible; "
            f"{report.verified_files} URLs HTTP OK; "
            f"{report.unknown_size_files} selected files have unknown sizes; "
            f"{report.missing_files} known HTTP-missing.")

    def _apply_profile(self, _=None):
        w = self.owner
        if w.manifest is None:
            return
        previous_selection = w._selected()
        from lsbpipe.surveys.legacy import build_manifest, DEFAULT_PRODUCTS
        kinds = DEFAULT_PRODUCTS if w.profile.currentIndex() == 0 else (
            DEFAULT_PRODUCTS + ("psfsize", "depth", "galdepth", "model"))
        w.manifest = build_manifest(w._field(), w.hits, kinds=kinds)
        w._fill_files(checked_paths=previous_selection)
