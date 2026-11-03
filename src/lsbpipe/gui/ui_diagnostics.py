"""Opt-in Windows UI artifact source diagnostics, never science data.

Run with DESIDERA_UI_DIAG=1. Save *only the top-left crop of this application's*
workspace (two rendering paths) and Qt widget geometry to an isolated TEMP
folder. Do not collect input values, URLs, project data, other windows or FITS.
No diagnostics run by default.
"""
from __future__ import annotations

import json
from pathlib import Path
import tempfile
import traceback

from PySide6.QtCore import QPoint, QRect, Qt
from PySide6.QtWidgets import QApplication, QWidget


def _widget_record(widget, workspace):
    origin = widget.mapTo(workspace, QPoint(0, 0))
    return {
        "class": type(widget).__name__,
        "object_name": widget.objectName(),
        "geometry_in_workspace": [
            origin.x(), origin.y(), widget.width(), widget.height()],
        "visible": widget.isVisible(),
        "is_window": widget.isWindow(),
        "native_window_attribute": widget.testAttribute(Qt.WA_NativeWindow),
        "translucent_background": widget.testAttribute(
            Qt.WA_TranslucentBackground),
        "opacity": widget.windowOpacity(),
    }


def capture_workspace_corner(workspace: QWidget):
    """Snapshot just the app's top-left 250×120 pixels for paint diagnosis.

    Widget.grab() obtains Qt's QWidget composition. Screen.grabWindow(winId())
    captures a window compositor surface, when supported. Comparing them can
    implicate a native overlay versus a regular Qt child (not definitive).
    """
    base = Path(tempfile.mkdtemp(prefix="desidera-ui-", dir=None))
    info = {
        "diagnosis": "GUI-only, no scientific or data-quality acceptance",
        "window_class": type(workspace).__name__,
        "workspace_size": [workspace.width(), workspace.height()],
        "captures": {},
        "widget_at_points": [],
        "intersecting_qt_widgets": [],
        "other_visible_top_level_widgets": [],
    }
    try:
        app = QApplication.instance()
        rect = QRect(0, 0, min(250, workspace.width()),
                     min(120, workspace.height()))
        widget_image = workspace.grab(rect)
        if not widget_image.isNull():
            target = base / "qt_widget_corner.png"
            if widget_image.save(str(target)):
                info["captures"]["qt_widget"] = str(target)
        try:
            screen = workspace.screen()
            if screen is not None:
                native_image = screen.grabWindow(
                    int(workspace.winId()), 0, 0, rect.width(), rect.height())
                if not native_image.isNull():
                    target = base / "native_window_corner.png"
                    if native_image.save(str(target)):
                        info["captures"]["native_window"] = str(target)
        except Exception as exc:
            info["captures"]["native_error"] = (
                type(exc).__name__ + ": " + str(exc))
        for x, y in [(14, 15), (25, 35), (45, 35),
                     (65, 35), (85, 35), (170, 35)]:
            widget = app.widgetAt(workspace.mapToGlobal(QPoint(x, y)))
            info["widget_at_points"].append({
                "point": [x, y],
                "class": type(widget).__name__ if widget else None,
                "object_name": widget.objectName() if widget else None})
        for widget in [workspace] + workspace.findChildren(QWidget):
            if not widget.isVisible() or widget.isWindow():
                continue
            try:
                record = _widget_record(widget, workspace)
            except RuntimeError:
                continue
            x, y, width, height = record["geometry_in_workspace"]
            if (width > 0 and height > 0 and x < rect.width() and
                    y < rect.height() and x + width > 0 and y + height > 0):
                info["intersecting_qt_widgets"].append(record)
        for widget in app.topLevelWidgets():
            if widget is workspace or not widget.isVisible():
                continue
            info["other_visible_top_level_widgets"].append({
                "class": type(widget).__name__,
                "object_name": widget.objectName(),
                "size": [widget.width(), widget.height()],
                "relative_position": [
                    widget.frameGeometry().x() - workspace.frameGeometry().x(),
                    widget.frameGeometry().y() - workspace.frameGeometry().y()],
            })
    except Exception as exc:
        info["diagnostic_exception"] = (
            type(exc).__name__ + ": " + str(exc))
        info["diagnostic_traceback"] = traceback.format_exc()[-1800:]
    finally:
        destination = base / "ui_diagnostic.json"
        destination.write_text(
            json.dumps(info, indent=2, ensure_ascii=False), encoding="utf-8")
        print("DESIDERA UI DIAGNOSTICS:", str(base), flush=True)
    return info
