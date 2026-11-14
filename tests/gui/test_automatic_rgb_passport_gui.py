"""P10e: real Qt automatic selected-object colour fetch and compact passport."""
from pathlib import Path

import pytest

pytest.importorskip("pytestqt")
pytest.importorskip("PySide6.QtWebEngineWidgets")

from PySide6.QtCore import QPoint, QBuffer, QByteArray, QIODevice
from PySide6.QtGui import QColor, QImage

from scripts.offline_gui_fixture import offline_legacy_viewer
from lsbpipe.gui import project_hub as hub_module
from lsbpipe.gui.remote_card_previews import fetch_remote_rgb_preview
from lsbpipe.surveys.legacy import SkyField
from lsbpipe.surveys.object_projects import create_object_project, save_object_project


def _photo():
    img = QImage(384, 384, QImage.Format_RGB32)
    img.fill(QColor(56, 92, 150))
    data = QByteArray()
    buffer = QBuffer(data)
    assert buffer.open(QIODevice.WriteOnly)
    assert img.save(buffer, "JPEG", 92)
    buffer.close()
    return bytes(data)


def _project(tmp_path, name, ra):
    p = create_object_project(
        name, tmp_path / name.replace(" ", "_"),
        field=SkyField(ra, 18.088025, 4))
    p.selected_dataset = ("dr11", "south", "1554p180", "r")
    p.data_root = p.root
    save_object_project(p)
    return p


def _hub(qtbot):
    with offline_legacy_viewer():
        instance = hub_module.ProjectHub()
    qtbot.addWidget(instance)
    instance.resize(1300, 810)
    instance.show()
    return instance


def test_colour_automatically_loads_for_selected_leo_p_no_visible_fetch_button(
        qtbot, monkeypatch, tmp_path):
    monkeypatch.delenv("DESIDERA_DISABLE_AUTOMATIC_RGB")
    project = _project(tmp_path, "Leo P", 155.438012)
    original_manifest = (project.root / "lsbpipe-object.json").read_bytes()
    seen = []
    jpeg = _photo()

    def fake_fetch(root):
        seen.append(Path(root).name)
        return fetch_remote_rgb_preview(
            root, transport=lambda url: (jpeg, "image/jpeg", url))

    monkeypatch.setattr(hub_module, "fetch_remote_rgb_preview", fake_fetch)
    hub = _hub(qtbot)
    assert not hub.rgb_btn.isVisible()
    assert not hub.preview_btn.isVisible()
    assert hub.catalogues_menu_btn.isVisible()
    assert hub.tools_menu_btn.isVisible()
    assert seen == []  # debounce before first network access
    qtbot.waitUntil(
        lambda: (hub._rgb_future is None and
                 hub.galaxy_card.snapshot is not None and
                 hub.galaxy_card.snapshot.preview_source == "REMOTE_LEGACY_RGB"),
        timeout=12000)
    assert seen == ["Leo_P"]
    assert hub.galaxy_card.preview.pixmap() is not None
    assert not hub.galaxy_card.preview.pixmap().isNull()
    assert not hub.galaxy_card.preview_note.isVisible()
    assert (project.root / "lsbpipe-object.json").read_bytes() == original_manifest

    card = hub.galaxy_card
    assert card.preview.width() == card.preview.height()
    assert card.preview.mapToGlobal(QPoint(0, 0)).x() < (
        card.name.mapToGlobal(QPoint(0, 0)).x())
    assert card.dataset.mapToGlobal(QPoint(0, 0)).y() >= (
        card.preview.mapToGlobal(QPoint(0, 0)).y() + card.preview.height())
    assert not card.details_panel.isVisible()
    card.details_toggle.click()
    assert card.details_panel.isVisible()
    assert card.path.text().startswith("Project folder:")

    # A cached image must not open another request on project refresh/reopen.
    hub.refresh()
    qtbot.wait(950)
    assert seen == ["Leo_P"]
    reopened = _hub(qtbot)
    qtbot.wait(950)
    assert reopened.galaxy_card.snapshot.preview_source == "REMOTE_LEGACY_RGB"
    assert seen == ["Leo_P"]


def test_fast_project_switch_debounces_only_last_selected_galaxy(
        qtbot, monkeypatch, tmp_path):
    monkeypatch.delenv("DESIDERA_DISABLE_AUTOMATIC_RGB")
    first = _project(tmp_path, "Leo P", 155.438012)
    second = _project(tmp_path, "Other Galaxy", 155.538012)
    seen = []
    jpeg = _photo()

    def fake_fetch(root):
        seen.append(Path(root).name)
        return fetch_remote_rgb_preview(
            root, transport=lambda url: (jpeg, "image/jpeg", url))

    monkeypatch.setattr(hub_module, "fetch_remote_rgb_preview", fake_fetch)
    hub = _hub(qtbot)
    chosen = 0
    for i in range(hub.projects.count()):
        if hub.projects.item(i).data(0x0100)["project_id"] == first.project_id:
            chosen = i
    hub.projects.setCurrentRow(chosen)
    # Revert quickly before 650ms debounce can dispatch.
    for i in range(hub.projects.count()):
        if hub.projects.item(i).data(0x0100)["project_id"] == second.project_id:
            hub.projects.setCurrentRow(i)
    qtbot.waitUntil(
        lambda: (hub._rgb_future is None and
                 hub.galaxy_card.snapshot is not None and
                 hub.galaxy_card.snapshot.project_id == second.project_id and
                 hub.galaxy_card.snapshot.preview_source == "REMOTE_LEGACY_RGB"),
        timeout=12000)
    assert seen == ["Other_Galaxy"]
    assert not (first.root / ".desidera" /
                "card-preview-remote-rgb.jpg").exists()


def test_failed_auto_download_is_not_retried_on_selection_or_refresh(
        qtbot, monkeypatch, tmp_path):
    monkeypatch.delenv("DESIDERA_DISABLE_AUTOMATIC_RGB")
    _project(tmp_path, "Leo P", 155.438012)
    seen = []

    def broken(root):
        seen.append(str(root))
        raise TimeoutError("viewer network offline")

    monkeypatch.setattr(hub_module, "fetch_remote_rgb_preview", broken)
    hub = _hub(qtbot)
    qtbot.waitUntil(
        lambda: bool(hub._rgb_future is None and seen and
                     "TimeoutError" in hub.rgb_feedback.text()),
        timeout=12000)
    assert len(seen) == 1
    hub.refresh()
    qtbot.wait(900)
    assert len(seen) == 1
    assert "Refresh colour preview" in [a.text() for a in
                                         hub.tools_menu.actions()]
    hub.rgb_action.trigger()  # explicit retry still available in tools menu
    qtbot.waitUntil(lambda: hub._rgb_future is None and len(seen) == 2,
                    timeout=12000)
    assert len(seen) == 2

def test_remote_rgb_auto_fetch_can_be_disabled_and_restored_from_hub(
        qtbot, monkeypatch, tmp_path):
    """A persisted user preference prevents background network requests."""
    from PySide6.QtCore import QSettings
    monkeypatch.delenv("DESIDERA_DISABLE_AUTOMATIC_RGB", raising=False)
    prefs = str(tmp_path / "hub-previews.ini")
    monkeypatch.setattr(
        hub_module, "QSettings",
        lambda *_: QSettings(prefs, QSettings.IniFormat))
    _project(tmp_path, "Auto RGB toggle galaxy", 155.438012)
    seen = []
    monkeypatch.setattr(
        hub_module, "fetch_remote_rgb_preview",
        lambda root: seen.append(str(root)))
    hub = _hub(qtbot)
    assert hub.preview_settings_btn.isVisible()
    assert hub.auto_rgb_action.isChecked()  # Preserves existing default.
    hub.auto_rgb_action.setChecked(False)
    assert not hub._rgb_auto_timer.isActive()
    qtbot.wait(800)
    assert seen == []  # No automatic network action.
    assert QSettings(prefs, QSettings.IniFormat).value(
        "previews/automatic_remote_rgb", True, type=bool) is False
    reopened = _hub(qtbot)
    assert not reopened.auto_rgb_action.isChecked()
    assert not reopened._rgb_auto_timer.isActive()
    qtbot.wait(750)
    assert seen == []
    manual_calls = []
    monkeypatch.setattr(reopened, "_fetch_rgb_preview",
                        lambda *, manual=True: manual_calls.append(manual))
    assert reopened.rgb_action.isEnabled()
    reopened.rgb_action.trigger()
    assert manual_calls == [True]
    # Opting in again schedules, but does not immediately perform, a fetch.
    reopened.auto_rgb_action.setChecked(True)
    assert reopened._rgb_auto_timer.isActive()
    reopened.auto_rgb_action.setChecked(False)
    assert not reopened._rgb_auto_timer.isActive()


def test_rgb_environment_override_is_visible_and_prevents_automatic_fetch(
        qtbot, monkeypatch, tmp_path):
    from PySide6.QtCore import QSettings
    prefs = str(tmp_path / "override.ini")
    monkeypatch.setattr(
        hub_module, "QSettings",
        lambda *_: QSettings(prefs, QSettings.IniFormat))
    monkeypatch.setenv("DESIDERA_DISABLE_AUTOMATIC_RGB", "1")
    _project(tmp_path, "Environment override target", 155.438012)
    hub = _hub(qtbot)
    assert not hub.auto_rgb_action.isEnabled()
    assert "environment override" in hub.auto_rgb_action.toolTip()
    assert not hub._rgb_auto_timer.isActive()
