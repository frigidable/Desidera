"""Explicit, project-bound review of an existing joint photometry archive."""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import zipfile

from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import (
    QAbstractItemView, QCheckBox, QComboBox, QDialog, QFileDialog, QHBoxLayout,
    QHeaderView, QLabel, QPlainTextEdit, QPushButton, QScrollArea, QSizePolicy,
    QTabWidget, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget,
)

from lsbpipe.photometry.joint_archive_review import (
    read_joint_archive, reopen_review_receipt, save_review_receipt,
)


class JointArchiveDialog(QDialog):
    """No scientific worker is launched; only user-selected files are read."""

    def __init__(self, *, project_name, project_root, context, parent=None):
        super().__init__(parent)
        self.setObjectName("joint_archive_dialog")
        self.setWindowTitle("Joint photometry · archive review")
        self.resize(900, 650)
        self.setMinimumSize(420, 360)
        self._root=Path(project_root)
        self._context=dict(context)
        self._review=None
        self._future=None
        self._executor=ThreadPoolExecutor(max_workers=1,thread_name_prefix="joint-archive-review")
        self._timer=QTimer(self)
        self._timer.setInterval(80)
        self._timer.timeout.connect(self._poll)
        self._pixmap=None
        layout=QVBoxLayout(self)
        heading=QLabel("Joint annular photometry")
        heading.setStyleSheet("font-size: 19px; font-weight: 600;")
        layout.addWidget(heading)
        self.identity=QLabel()
        self.identity.setTextFormat(Qt.PlainText)
        self.identity.setWordWrap(True)
        self.identity.setText(f"{project_name} · {' / '.join(context['source_triplet'])}")
        layout.addWidget(self.identity)
        note=QLabel("DIAGNOSTIC_NOT_ADOPTED · Review existing results only. "
                    "Original FITS are not reopened; archive checksums do not establish scientific acceptance.")
        note.setWordWrap(True)
        layout.addWidget(note)
        actions=QHBoxLayout()
        self.open_btn=QPushButton("Open diagnostic ZIP…")
        self.open_btn.setObjectName("joint_open_zip")
        self.open_btn.clicked.connect(self._choose_archive)
        actions.addWidget(self.open_btn)
        self.reopen_btn=QPushButton("Reopen review…")
        self.reopen_btn.setObjectName("joint_reopen_receipt")
        self.reopen_btn.clicked.connect(self._choose_receipt)
        actions.addWidget(self.reopen_btn)
        actions.addStretch(1)
        layout.addLayout(actions)
        self.status=QLabel("No archive selected. Open a diagnostic ZIP from the joint workflow.")
        self.status.setObjectName("joint_archive_status")
        self.status.setTextFormat(Qt.PlainText)
        self.status.setWordWrap(True)
        layout.addWidget(self.status)

        self.tabs=QTabWidget()
        self.tabs.setMinimumWidth(0)
        self.tabs.setSizePolicy(QSizePolicy.Ignored,QSizePolicy.Expanding)
        self.summary=QPlainTextEdit()
        self.summary.setReadOnly(True)
        self.summary.setAccessibleName("Archive provenance and structural fit summary")
        self.tabs.addTab(self.summary,"Summary")
        self.rings=QTableWidget(0,7)
        self.rings.setObjectName("joint_ring_table")
        self.rings.setHorizontalHeaderLabels(
            ["Band","a inner–outer [arcsec]","Observed [nMgy]","Model [nMgy]",
             "Conditional σ [nMgy]","Coverage","Role"])
        self.rings.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.rings.verticalHeader().setVisible(False)
        self.rings.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.rings.setMinimumWidth(0)
        self.tabs.addTab(self.rings,"Signed rings")
        page=QWidget()
        panel=QVBoxLayout(page)
        panel.setContentsMargins(4,4,4,4)
        self.figure_choice=QComboBox()
        self.figure_choice.setAccessibleName("Archived diagnostic figure")
        self.figure_choice.currentIndexChanged.connect(self._show_figure)
        panel.addWidget(self.figure_choice)
        self.fit_view=QCheckBox("Fit figure to view")
        self.fit_view.setChecked(True)
        self.fit_view.toggled.connect(self._scale_figure)
        panel.addWidget(self.fit_view)
        self.scroll=QScrollArea()
        self.scroll.setWidgetResizable(True)
        self.image=QLabel("No figures available.")
        self.image.setAlignment(Qt.AlignCenter)
        self.image.setWordWrap(True)
        self.image.setMinimumSize(0,0)
        self.image.setSizePolicy(QSizePolicy.Ignored,QSizePolicy.Ignored)
        self.scroll.setWidget(self.image)
        panel.addWidget(self.scroll,1)
        self.tabs.addTab(page,"Figures")
        layout.addWidget(self.tabs,1)
        self.confirm=QCheckBox("I confirm that this archive belongs to the current galaxy.")
        self.confirm.setToolTip("A matching brick and nearby anchor do not establish galaxy identity.")
        self.confirm.toggled.connect(self._update_actions)
        layout.addWidget(self.confirm)
        footer=QHBoxLayout()
        self.save_btn=QPushButton("Save review receipt…")
        self.save_btn.setObjectName("joint_save_receipt")
        self.save_btn.clicked.connect(self._choose_save)
        footer.addWidget(self.save_btn)
        footer.addStretch(1)
        close=QPushButton("Close")
        close.clicked.connect(self.reject)
        footer.addWidget(close)
        layout.addLayout(footer)
        self._clear()

    def _clear(self):
        self._review=None
        self.confirm.setChecked(False)
        self.summary.setPlainText("No verified archive loaded.")
        self.rings.setRowCount(0)
        self.figure_choice.clear()
        self._pixmap=None
        self.image.clear()
        self.image.setText("No figures available.")
        self.tabs.setCurrentIndex(0)
        self._update_actions()

    def _update_actions(self,*_):
        busy=self._future is not None
        self.open_btn.setEnabled(not busy)
        self.reopen_btn.setEnabled(not busy)
        self.confirm.setEnabled(not busy and self._review is not None)
        self.save_btn.setEnabled(not busy and self._review is not None and self.confirm.isChecked())
        self.figure_choice.setEnabled(self.figure_choice.count()>0 and not busy)

    def _choose_archive(self):
        path,_=QFileDialog.getOpenFileName(self,"Open joint diagnostic",str(self._root),"Diagnostic ZIP (*.zip)")
        if path:
            self.start_load(path)

    def _choose_receipt(self):
        path,_=QFileDialog.getOpenFileName(self,"Reopen joint review",str(self._root),"Review JSON (*.json)")
        if path:
            self.start_load(path,receipt=True)

    def start_load(self,path,*,receipt=False):
        """Called only after an explicit local file choice (or offline test)."""
        if self._future is not None:
            return
        self._clear()
        self.status.setText("Verifying archive and project association…")
        fn=reopen_review_receipt if receipt else read_joint_archive
        self._future=self._executor.submit(fn,path,**self._context)
        self._update_actions()
        self._timer.start()

    def _poll(self):
        if self._future is None or not self._future.done():
            return
        self._timer.stop()
        future=self._future
        self._future=None
        try:
            result=future.result()
            self._populate(result)
        except Exception as exc:
            # Worker errors are shown as failure, never a partially loaded review.
            self._clear()
            self.status.setText(f"Review failed: {type(exc).__name__}: {exc}. Choose another archive or receipt.")
        self._update_actions()

    def _populate(self,result):
        self._review=result
        fit=result["fit"]
        bands=result["manifest"]["band_order"]
        self.status.setText(
            f"Archive member hashes verified · {', '.join(bands)} · "
            f"anchor offset {result['anchor_separation_arcsec']:.3f} arcsec. "
            "Confirm galaxy association before saving a receipt.")
        lines=[
            "SCIENCE HOLD — diagnostic, not adopted",
            f"Archive: {result['archive_path']}",
            f"SHA-256: {result['archive_sha256']}",
            f"Archive anchor ICRS: {result['archive_anchor_icrs']}",
            f"Project target ICRS: {result['context']['target_icrs']}",
            f"Anchor separation: {result['anchor_separation_arcsec']:.4f} arcsec",
            "Original FITS reverified: no",
            "Covariance calibrated: no",
            "Full contour containment verified: no",
            "",
            f"Band order: {', '.join(bands)}",
            f"Conditional fit chi-square: {fit.get('chi2','not recorded')}",
            f"Degrees of freedom: {fit.get('dof','not recorded')}",
            f"Covariance hypothesis: {fit.get('covariance_source','not recorded')}",
        ]
        for entry in result["manifest"]["inputs"]:
            lines.extend(["",entry["dataset_id"],
                          "Native array SHA-256: "+entry["native_input_sha256"],
                          "Grid: "+entry["grid_id"]])
        held=fit.get("holdout")
        if isinstance(held,dict):
            lines.extend(["",f"Held-out conditional chi-square: {held.get('chi2','not recorded')}"])
        if result["asymptotic"] is not None:
            lines.extend(["","Asymptotic diagnostic present (model-dependent, not adopted)."])
        lines.extend(["",f"Available PNG figures: {len(result['figures'])}",
                      "Negative observed flux is retained in the signed ring table.",
                      "Archive review does not run the optimizer or update project science."])
        self.summary.setPlainText("\n".join(lines))
        self.rings.setRowCount(len(result["rows"]))
        for i,row in enumerate(result["rows"]):
            role="train" if row["train"] else "holdout" if row["test"] else "unused"
            values=(row["band"],f"{row['a_inner_arcsec']:.4g}–{row['a_outer_arcsec']:.4g}",
                    f"{row['observed_nMgy']:.7g}",f"{row['model_same_support_nMgy']:.7g}",
                    f"{row['conditional_sigma_nMgy']:.7g}",
                    "eligible" if row["eligible"] else "excluded",role)
            for j,v in enumerate(values):
                item=QTableWidgetItem(v)
                self.rings.setItem(i,j,item)
        self.figure_choice.clear()
        for name in result["figures"]:
            self.figure_choice.addItem(name.replace("_"," "),name)
        self._show_figure()
        self._update_actions()

    def _show_figure(self,*_):
        key=self.figure_choice.currentData()
        self._pixmap=None
        self.image.clear()
        if self._review is None or key not in self._review["figures"]:
            self.image.setText("No archived PNG figures available.")
            return
        pixmap=QPixmap()
        if not pixmap.loadFromData(self._review["figures"][key],"PNG"):
            self.image.setText("Archived PNG could not be decoded. Numerical review remains available.")
            return
        self._pixmap=pixmap
        self._scale_figure()

    def _scale_figure(self,*_):
        if self._pixmap is None:
            self.scroll.setWidgetResizable(True)
            self.image.setMinimumSize(0,0)
            return
        fit=self.fit_view.isChecked()
        self.scroll.setWidgetResizable(fit)
        if fit:
            self.image.setMinimumSize(0,0)
            size=self.scroll.viewport().size()
            if size.width()>0 and size.height()>0:
                self.image.setPixmap(self._pixmap.scaled(size,Qt.KeepAspectRatio,Qt.SmoothTransformation))
        else:
            self.image.setPixmap(self._pixmap)
            self.image.resize(self._pixmap.size())

    def resizeEvent(self,event):
        super().resizeEvent(event)
        if hasattr(self,"_pixmap"):
            self._scale_figure()

    def _choose_save(self):
        if self._review is None or not self.confirm.isChecked():
            return
        path,_=QFileDialog.getSaveFileName(self,"Save archive review receipt",
                                          str(self._root/"joint-review.json"),"Review JSON (*.json)")
        if path:
            try:
                save_review_receipt(path,self._review,confirmed=True)
                self.status.setText("Review receipt saved. Scientific results remain unadopted.")
            except (OSError,ValueError) as exc:
                self.status.setText(f"Receipt not saved: {exc}")

    def done(self,result):
        self._timer.stop()
        if self._future is not None:
            self._future.cancel()
            self._future=None
        self._executor.shutdown(wait=False,cancel_futures=True)
        super().done(result)
