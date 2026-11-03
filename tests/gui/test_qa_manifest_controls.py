"""GUI regression: immutable per-run provenance is directly discoverable."""
import pytest

pytest.importorskip("pytestqt")
pytest.importorskip("PySide6.QtWidgets")

from lsbpipe.gui.qa_workbench import ScienceQAControls


def test_manifest_receipt_visible_only_for_current_completed_science_run(qtbot):
    view = ScienceQAControls()
    qtbot.addWidget(view)
    assert not view.open_manifest.isEnabled()
    assert not view.open_report.isEnabled()
    view.show_result({
        "action": "square",
        "summary": "DR9 South · original FITS · 99.8% valid",
        "path": "C:/test/qa-runs/unique/square.json",
        "manifest_path": "C:/test/qa-runs/unique/manifest.json",
        "run_id": "unique",
        "is_figure": False,
    })
    assert view.open_manifest.isEnabled()
    assert view.open_report.isEnabled()
    assert "Run ID: unique" in view.result.toPlainText()
    assert "manifest.json" in view.result.toPlainText()
    view.directory.setText("C:/test/different-pilot")
    assert not view.open_manifest.isEnabled()
    assert not view.open_report.isEnabled()
    assert view._last_manifest_path is None


def test_failed_or_legacy_qa_cannot_open_stale_science_manifest(qtbot):
    view = ScienceQAControls()
    qtbot.addWidget(view)
    view.show_result({
        "action": "figure", "summary": "Earlier result",
        "path": "C:/temp/figure.png",
        "manifest_path": "C:/temp/manifest.json",
        "run_id": "old",
        "is_figure": True,
    })
    view.show_error("Invalid WCS")
    assert not view.open_manifest.isEnabled()
    assert not view.open_report.isEnabled()
    assert "Invalid WCS" in view.result.toPlainText()
    view.show_result({
        "action": "square", "summary": "Unversioned legacy diagnostic",
        "path": "C:/temp/old-square.json", "is_figure": False,
    })
    assert view.open_report.isEnabled()
    assert not view.open_manifest.isEnabled()
    assert "Run ID" not in view.result.toPlainText()



def test_scientific_run_history_button_emits_and_dialog_is_read_only(qtbot, tmp_path):
    from lsbpipe.gui.run_history_dialog import RunHistoryDialog
    view = ScienceQAControls()
    qtbot.addWidget(view)
    events = []
    view.historyRequested.connect(lambda: events.append("open"))
    assert view.history_btn.isEnabled()
    view.history_btn.click()
    assert events == ["open"]
    folder = tmp_path / "qa-runs" / "20261102T090000000000Z-abcdef1234"
    folder.mkdir(parents=True)
    manifest = folder / "manifest.json"
    manifest.write_text('{"schema":"legacy-only"}', encoding="utf-8")
    records = [{
        "run_id": folder.name, "action": "square",
        "dataset": "dr9/south/0401p012/r",
        "state": "RECORDED_NOT_REVALIDATED",
        "created_at_utc": "2026-11-02T09:00:00Z",
        "manifest_path": str(manifest),
        "reason": "No current local original-FITS hash comparison",
    }]
    dialog = RunHistoryDialog(records)
    qtbot.addWidget(dialog)
    dialog.show()
    assert dialog.table.rowCount() == 1
    assert "RECORDED_NOT_REVALIDATED" in (
        dialog.table.item(0, 3).text())
    assert "not automatically" not in dialog.details.toPlainText().lower() or (
        "not" in dialog.details.toPlainText().lower())
    assert dialog.open_manifest.isEnabled()
    assert dialog.details.isReadOnly()
