"""P0.8a GUI: freeze setup explicitly, review history offline."""
import pytest

pytest.importorskip("pytestqt")
pytest.importorskip("PySide6.QtWidgets")

from PySide6.QtWidgets import QInputDialog
from lsbpipe.gui.analysis_configurations_dialog import AnalysisConfigurationsDialog
from lsbpipe.surveys.object_projects import create_object_project
from lsbpipe.surveys.legacy import SkyField


def test_analysis_config_dialog_freeze_and_reopen(qtbot, tmp_path, monkeypatch):
    project = create_object_project(
        "UGC 2162", tmp_path / "UGC_2162",
        field=SkyField(40.096208, 1.22925, 5))
    manifest = project.root / "lsbpipe-object.json"
    original = manifest.read_bytes()
    dlg = AnalysisConfigurationsDialog(project.root)
    qtbot.addWidget(dlg)
    dlg.show()
    assert dlg.table.rowCount() == 0
    monkeypatch.setattr(
        QInputDialog, "getText", lambda *args, **kwargs: ("Initial baseline", True))
    dlg.freeze_btn.click()
    assert dlg.table.rowCount() == 1
    assert dlg.table.item(0, 0).text() == "Initial baseline"
    assert "not rehashed" in dlg.status.text().lower()
    assert manifest.read_bytes() == original
    reopened = AnalysisConfigurationsDialog(project.root)
    qtbot.addWidget(reopened)
    assert reopened.table.rowCount() == 1
    assert reopened.table.item(0, 0).text() == "Initial baseline"
