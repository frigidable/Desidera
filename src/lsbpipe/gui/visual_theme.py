"""Desidera visual language for Qt Widgets (presentation only).

No global QApplication palette changes, science-rendered pixels, or files.
Colors and vector symbols are scoped to the widget owning the stylesheet.
Keep semantic status text: color is supplemental, never science acceptance.
"""
from __future__ import annotations

from PySide6.QtCore import QPointF, QRectF, QSize, Qt
from PySide6.QtGui import QColor, QIcon, QPainter, QPainterPath, QPalette, QPen, QPixmap


def theme_colors(widget) -> dict[str, str]:
    """Two legible palettes selected from the current OS Window role."""
    dark = widget.palette().color(QPalette.Window).lightness() < 128
    if dark:
        return {
            "BG": "#0C1424",
            "SURFACE": "#172238",
            "RAISED": "#1D2C46",
            "TEXT": "#E9F1FD",
            "MUTED": "#AABAD1",
            "BORDER": "#344862",
            "ACCENT": "#65B8FA",
            "ACCENT_HOVER": "#8BCAFC",
            "PRIMARY_TEXT": "#071629",
            "FOCUS": "#74D8D1",
            "SOFT_ACCENT": "#1A3552",
            "WARN": "#EEC47A",
            "WARN_BG": "#33291D",
            "ERROR": "#F3A0A5",
            "ERROR_BG": "#402630",
        }
    return {
        "BG": "#F2F5FA",
        "SURFACE": "#FFFFFF",
        "RAISED": "#E9F0FA",
        "TEXT": "#1A2B45",
        "MUTED": "#4D617A",
        "BORDER": "#CED9E8",
        "ACCENT": "#215CAA",
        "ACCENT_HOVER": "#174B8D",
        "PRIMARY_TEXT": "#FFFFFF",
        "FOCUS": "#117F8C",
        "SOFT_ACCENT": "#E3F0FA",
        "WARN": "#865313",
        "WARN_BG": "#FFF3DE",
        "ERROR": "#A33848",
        "ERROR_BG": "#FFF0F1",
    }


def visual_stylesheet(widget) -> str:
    """Build contrast-aware QSS; never style FITS/QGraphicsView pixels."""
    css = """
QMainWindow#lsbpipe_project_hub,
QMainWindow#desidera_survey_explorer,
QWidget#photometry_workspace {
    background-color: $BG;
    color: $TEXT;
}
QWidget#survey_navigation {
    background-color: $SURFACE;
    border-right: 1px solid $BORDER;
}
QFrame#hub_project_details {
    background-color: $SURFACE;
    border: 1px solid $BORDER;
    border-radius: 12px;
}
QLabel { color: $TEXT; }
QLabel#hub_brand_title, QLabel#photometry_heading {
    color: $ACCENT;
    font-weight: 700;
}
QLabel#workflow_stage, QLabel#hub_detail_heading {
    color: $ACCENT;
    font-weight: 650;
}
QLabel#photometry_not_adopted, QLabel#photometry_blocker {
    background-color: $WARN_BG;
    border-left: 3px solid $WARN;
    border-radius: 5px;
    color: $WARN;
    padding: 8px 11px;
}
QLabel#sky_diagnostic_status, QLabel#aperture_stage_status,
QLabel#scene_status {
    color: $MUTED;
    background-color: $RAISED;
    border-left: 3px solid $ACCENT;
    border-radius: 5px;
    padding: 7px 10px;
}
QLabel#sky_diagnostic_status[visualState="error"],
QLabel#aperture_stage_status[visualState="error"],
QLabel#scene_status[visualState="error"] {
    color: $ERROR;
    background-color: $ERROR_BG;
    border-left: 3px solid $ERROR;
}
QGroupBox {
    background-color: $SURFACE;
    border: 1px solid $BORDER;
    border-radius: 10px;
    margin-top: 11px;
    padding: 10px 8px 7px 8px;
    font-weight: 550;
    color: $TEXT;
}
QGroupBox::title {
    subcontrol-origin: margin;
    subcontrol-position: top left;
    left: 12px;
    padding: 0px 6px;
    color: $ACCENT;
}
QPushButton {
    background-color: $SURFACE;
    color: $TEXT;
    border: 1px solid $BORDER;
    border-radius: 7px;
    padding: 6px 10px;
    min-height: 27px;
}
QPushButton:hover:!disabled {
    background-color: $SOFT_ACCENT;
    border-color: $ACCENT;
}
QPushButton:focus {
    border: 2px solid $FOCUS;
}
QPushButton:disabled {
    background-color: $BG;
    color: $MUTED;
    border-color: $BORDER;
}
QPushButton[visualRole="primary"], QPushButton[primaryAction="true"] {
    background-color: $ACCENT;
    color: $PRIMARY_TEXT;
    border: 1px solid $ACCENT;
    border-radius: 7px;
    font-weight: 650;
    min-height: 31px;
}
QPushButton[visualRole="primary"]:hover:!disabled,
QPushButton[primaryAction="true"]:hover:!disabled {
    background-color: $ACCENT_HOVER;
}
QPushButton[visualRole="primary"]:disabled,
QPushButton[primaryAction="true"]:disabled {
    background-color: $RAISED;
    color: $MUTED;
    border-color: $BORDER;
}
QPushButton[visualRole="caution"] {
    color: $WARN;
    border-color: $WARN;
    background-color: $WARN_BG;
}
QTabWidget::pane {
    border: 1px solid $BORDER;
    border-radius: 8px;
    background-color: $SURFACE;
}
QTabBar::tab {
    color: $MUTED;
    background-color: $BG;
    border: 1px solid $BORDER;
    border-bottom: 0px;
    border-top-left-radius: 6px;
    border-top-right-radius: 6px;
    padding: 7px 10px;
    margin-right: 3px;
}
QTabBar::tab:selected {
    background-color: $SURFACE;
    color: $ACCENT;
    border-top: 3px solid $ACCENT;
    padding-top: 5px;
    font-weight: 650;
}
QTabBar::tab:hover:!selected {
    background-color: $SOFT_ACCENT;
    color: $TEXT;
}
QLineEdit, QComboBox, QDoubleSpinBox, QPlainTextEdit,
QTableWidget, QListWidget {
    background-color: $SURFACE;
    color: $TEXT;
    border: 1px solid $BORDER;
    border-radius: 6px;
    selection-background-color: $ACCENT;
    selection-color: $PRIMARY_TEXT;
}
QLineEdit:focus, QComboBox:focus, QDoubleSpinBox:focus,
QPlainTextEdit:focus, QTableWidget:focus, QListWidget:focus {
    border: 2px solid $FOCUS;
}
QHeaderView::section {
    background-color: $RAISED;
    color: $TEXT;
    border: 0px;
    border-bottom: 1px solid $BORDER;
    padding: 6px 8px;
}
QProgressBar {
    background-color: $RAISED;
    color: $TEXT;
    border: 1px solid $BORDER;
    border-radius: 5px;
    text-align: center;
}
QProgressBar::chunk {
    background-color: $ACCENT;
    border-radius: 4px;
}
QToolTip {
    background-color: $SURFACE;
    color: $TEXT;
    border: 1px solid $BORDER;
    padding: 4px;
}
"""
    for key, color in theme_colors(widget).items():
        css = css.replace("$" + key, color)
    return css


def icon_for(name: str, *, color: str = "#215CAA", size: int = 18) -> QIcon:
    """A small platform-independent vector icon, no fonts or asset files."""
    known = {
        "galaxy", "folder", "plus", "back", "database", "chart",
        "download", "search", "refresh", "copy", "settings",
        "mask", "sky", "aperture", "open", "check", "save",
    }
    if name not in known:
        raise ValueError(f"Unknown Desidera icon: {name}")
    pixmap = QPixmap(size * 2, size * 2)
    pixmap.setDevicePixelRatio(2)
    pixmap.fill(Qt.transparent)
    painter = QPainter(pixmap)
    try:
        painter.setRenderHint(QPainter.Antialiasing, True)
        painter.scale(size / 24.0 * 2.0, size / 24.0 * 2.0)
        painter.setPen(QPen(QColor(color), 1.9, Qt.SolidLine,
                            Qt.RoundCap, Qt.RoundJoin))
        painter.setBrush(Qt.NoBrush)
        if name == "galaxy":
            painter.save()
            painter.translate(12, 12)
            painter.rotate(-32)
            painter.drawEllipse(QRectF(-10, -4, 20, 8))
            painter.restore()
            painter.drawEllipse(QRectF(9, 9, 6, 6))
            painter.drawLine(18, 3, 18, 7)
            painter.drawLine(16, 5, 20, 5)
        elif name == "folder":
            path = QPainterPath()
            path.moveTo(3, 7)
            path.lineTo(9, 7)
            path.lineTo(11, 9)
            path.lineTo(21, 9)
            path.lineTo(21, 19)
            path.lineTo(3, 19)
            path.closeSubpath()
            painter.drawPath(path)
        elif name == "plus":
            painter.drawEllipse(QRectF(3, 3, 18, 18))
            painter.drawLine(12, 7, 12, 17)
            painter.drawLine(7, 12, 17, 12)
        elif name == "back":
            painter.drawLine(19, 12, 5, 12)
            painter.drawLine(11, 6, 5, 12)
            painter.drawLine(5, 12, 11, 18)
        elif name == "database":
            painter.drawEllipse(QRectF(4, 4, 16, 6))
            painter.drawLine(4, 7, 4, 18)
            painter.drawLine(20, 7, 20, 18)
            painter.drawArc(QRectF(4, 13, 16, 8), 180 * 16, 180 * 16)
            painter.drawArc(QRectF(4, 8, 16, 8), 180 * 16, 180 * 16)
        elif name == "chart":
            painter.drawLine(4, 3, 4, 20)
            painter.drawLine(4, 20, 21, 20)
            path = QPainterPath()
            path.moveTo(6, 16)
            path.lineTo(11, 12)
            path.lineTo(15, 14)
            path.lineTo(20, 6)
            painter.drawPath(path)
        elif name == "download":
            painter.drawLine(12, 3, 12, 15)
            painter.drawLine(7, 10, 12, 15)
            painter.drawLine(12, 15, 17, 10)
            painter.drawLine(5, 20, 19, 20)
            painter.drawLine(5, 17, 5, 20)
            painter.drawLine(19, 17, 19, 20)
        elif name == "search":
            painter.drawEllipse(QRectF(3.5, 3.5, 13, 13))
            painter.drawLine(15, 15, 21, 21)
        elif name == "refresh":
            painter.drawArc(QRectF(4, 4, 16, 16), 24 * 16, 300 * 16)
            painter.drawLine(18, 4, 20, 9)
            painter.drawLine(20, 9, 15, 10)
        elif name == "save":
            painter.drawRoundedRect(QRectF(4, 3, 16, 18), 2, 2)
            painter.drawRect(QRectF(7, 3, 10, 7))
            painter.drawRoundedRect(QRectF(7, 13, 10, 8), 1, 1)
        elif name == "copy":
            painter.drawRoundedRect(QRectF(8, 7, 12, 14), 2, 2)
            painter.drawLine(5, 17, 4, 17)
            painter.drawLine(4, 17, 4, 3)
            painter.drawLine(4, 3, 16, 3)
        elif name == "settings":
            painter.drawLine(4, 6, 20, 6)
            painter.drawLine(4, 12, 20, 12)
            painter.drawLine(4, 18, 20, 18)
            painter.drawEllipse(QRectF(8, 3, 6, 6))
            painter.drawEllipse(QRectF(14, 9, 6, 6))
            painter.drawEllipse(QRectF(6, 15, 6, 6))
        elif name == "mask":
            painter.drawRoundedRect(QRectF(4, 4, 16, 16), 3, 3)
            painter.drawEllipse(QRectF(7, 7, 10, 10))
            painter.drawLine(5, 19, 19, 5)
        elif name == "sky":
            painter.drawLine(12, 3, 12, 7)
            painter.drawLine(12, 17, 12, 21)
            painter.drawLine(3, 12, 7, 12)
            painter.drawLine(17, 12, 21, 12)
            painter.drawEllipse(QRectF(8, 8, 8, 8))
            painter.drawLine(5, 5, 7, 7)
            painter.drawLine(17, 17, 19, 19)
        elif name == "aperture":
            painter.drawEllipse(QRectF(3, 3, 18, 18))
            painter.drawEllipse(QRectF(7, 7, 10, 10))
            painter.drawEllipse(QRectF(10, 10, 4, 4))
        elif name == "open":
            painter.drawRoundedRect(QRectF(4, 4, 16, 16), 3, 3)
            painter.drawLine(8, 12, 16, 12)
            painter.drawLine(13, 9, 16, 12)
            painter.drawLine(16, 12, 13, 15)
        elif name == "check":
            painter.drawEllipse(QRectF(3, 3, 18, 18))
            painter.drawLine(7, 12, 10, 15)
            painter.drawLine(10, 15, 17, 8)
    finally:
        painter.end()
    return QIcon(pixmap)


def adorn(button, icon_name: str, *, primary: bool = False,
          caution: bool = False) -> None:
    """Style a real action while keeping its label, tooltip and activation."""
    if primary and caution:
        raise ValueError("An action cannot be primary and caution together")
    tokens = theme_colors(button)
    color = (tokens["PRIMARY_TEXT"] if primary else
             tokens["WARN"] if caution else tokens["ACCENT"])
    button.setIcon(icon_for(icon_name, color=color))
    button.setIconSize(QSize(18, 18))
    if primary:
        button.setProperty("visualRole", "primary")
    if caution:
        button.setProperty("visualRole", "caution")


def set_status_tone(label, *, error: bool) -> None:
    """Repolish a semantic Qt state; visible text remains authoritative."""
    desired = "error" if error else "normal"
    if label.property("visualState") == desired:
        return
    label.setProperty("visualState", desired)
    style = label.style()
    style.unpolish(label)
    style.polish(label)
    label.update()
