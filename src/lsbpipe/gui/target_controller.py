"""Target-entry Qt widgets and isolated galaxy center controller.

This owns editable target input/validation only. It never downloads FITS,
chooses an observational dataset or changes the project's UUID.
"""
from __future__ import annotations

from PySide6.QtCore import QObject
from PySide6.QtWidgets import QApplication, QLineEdit, QMessageBox
from lsbpipe.surveys.legacy import SkyField
from lsbpipe.gui.target_input import parse_coordinate_pair

class TargetTextField(QLineEdit):
    """Ordinary pasteable text field with strict read-time parsing.

    Unlike QDoubleSpinBox it has no spin buttons, locale-dependent input
    filtering, or silently clamped values. Invalid drafts remain editable.
    """

    def __init__(self, initial: float, parser, *, decimals: int = 6, parent=None):
        super().__init__(parent)
        self._parser = parser
        self._decimals = decimals
        self.setProperty("targetInput", True)
        self.setClearButtonEnabled(True)
        self.setMinimumHeight(35)
        self.setValue(initial)
        self.editingFinished.connect(self._validate_quietly)
        self.textEdited.connect(self._clear_invalid)

    def value(self) -> float:
        return self._parser(self.text())

    def setValue(self, value: float):
        self.setText(f"{float(value):.{self._decimals}f}".rstrip("0").rstrip("."))
        self._clear_invalid()

    def _mark_invalid(self, reason: str):
        if not self.property("invalid"):
            self._normal_tooltip = self.toolTip()
        self.setProperty("invalid", True)
        self.setToolTip(reason)
        self._refresh_style()

    def _clear_invalid(self, *_):
        if self.property("invalid"):
            self.setProperty("invalid", False)
            self.setToolTip(getattr(self, "_normal_tooltip", ""))
            self._refresh_style()

    def _refresh_style(self):
        self.style().unpolish(self)
        self.style().polish(self)
        self.update()

    def _validate_quietly(self):
        try:
            self.value()
        except ValueError as exc:
            self._mark_invalid(str(exc))
        else:
            self._clear_invalid()

class TargetController(QObject):
    """Adapt target widgets without owning survey/science execution."""

    def __init__(self, owner):
        super().__init__(owner)
        self.owner = owner

    def field(self):
        w = self.owner
        return SkyField(w.ra.value(), w.dec.value(), w.field_width.value())

    def validated_field(self):
        """Reject invalid ICRS drafts without changing adopted coordinates."""
        w = self.owner
        for widget in (w.ra, w.dec, w.field_width):
            try:
                widget.value()
            except ValueError as exc:
                widget._mark_invalid(str(exc))
                widget.setFocus()
                widget.selectAll()
                w.target_hint.setText(str(exc))
                QMessageBox.warning(w, "Invalid target coordinates", str(exc))
                return None
        try:
            field = self.field()
        except ValueError as exc:
            w.target_hint.setText(str(exc))
            QMessageBox.warning(w, "Invalid target coordinates", str(exc))
            return None
        w.target_hint.setText("ICRS target ready · press Find data or Show target")
        return field

    def edited(self, *_):
        w = self.owner
        if w._restoring_local_project:
            return
        w._manual_target_edited = True
        if hasattr(w, "science_qa"):
            w._on_science_context_changed()
        w.target_hint.setText(
            "Target changed · run Find data again to refresh the bricks"
            if w.manifest is not None else
            "Target edited · press Find data or Show target")

    def paste_coordinate_pair(self):
        """Explicit clipboard action; never read the clipboard during editing."""
        w = self.owner
        text = QApplication.clipboard().text()
        try:
            ra, dec = parse_coordinate_pair(text)
        except ValueError as exc:
            w.target_hint.setText(str(exc))
            QMessageBox.warning(w, "Cannot paste RA / Dec", str(exc))
            return
        w.ra.setValue(ra)
        w.dec.setValue(dec)
        w.target_hint.setText(
            f"Pasted ICRS: RA {ra:.6f}°, Dec {dec:+.6f}°")

    def preferred_viewer_region(self, release: str):
        """North/south remote preview only for field wholly on one side."""
        field = self.field()
        half = field.width_arcmin / 120.0
        boundary = 32.375
        if field.dec_deg - half > boundary and release in ("dr9", "dr11"):
            return "north"
        if field.dec_deg + half < boundary:
            return "south"
        return None
