"""GUI-only top-left paint provenance diagnostics, never original FITS."""
import json

from PySide6.QtWidgets import QWidget, QPushButton, QVBoxLayout

from lsbpipe.gui.ui_diagnostics import capture_workspace_corner


def test_captured_qt_corner_and_intersecting_widget_tree(qtbot, tmp_path,
                                                         monkeypatch):
    from lsbpipe.gui import ui_diagnostics
    monkeypatch.setattr(
        ui_diagnostics.tempfile, "mkdtemp",
        lambda prefix, dir=None: str(tmp_path))
    root = QWidget()
    root.setObjectName("diagnostic_test_workspace")
    root.resize(360, 220)
    layout = QVBoxLayout(root)
    button = QPushButton("All galaxy projects")
    button.setObjectName("workspace_back_to_projects")
    layout.addWidget(button)
    layout.addStretch()
    qtbot.addWidget(root)
    root.show()
    info = capture_workspace_corner(root)
    file = tmp_path / "ui_diagnostic.json"
    assert file.exists()
    record = json.loads(file.read_text(encoding="utf-8"))
    assert info["workspace_size"] == [360, 220]
    assert record["captures"].get("qt_widget") is not None
    assert (tmp_path / "qt_widget_corner.png").exists()
    assert any(
        item["object_name"] == "workspace_back_to_projects"
        for item in record["intersecting_qt_widgets"])
    assert len(record["widget_at_points"]) == 6
    assert all("text" not in item for item in
               record["intersecting_qt_widgets"])
