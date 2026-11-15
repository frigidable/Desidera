"""Real Qt Galaxy Hub: RGB on explicit click, square layout, offline fallback."""
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


def _jpeg() -> bytes:
    image = QImage(384, 384, QImage.Format_RGB32)
    image.fill(QColor(45, 92, 135))
    output = QByteArray()
    writer = QBuffer(output)
    assert writer.open(QIODevice.WriteOnly)
    assert image.save(writer, "JPEG", 90)
    writer.close()
    return bytes(output)


def _project(tmp_path):
    project = create_object_project(
        "Leo P", tmp_path / "Leo_P",
        field=SkyField(155.438012, 18.088025, 4.))
    project.selected_dataset = ("dr11", "south", "1554p180", "r")
    project.data_root = project.root
    save_object_project(project)
    return project


def test_square_passport_and_manual_colour_fetch_then_offline_reopen(
        qtbot, tmp_path, monkeypatch):
    project = _project(tmp_path)
    manifest = project.root / "lsbpipe-object.json"
    before = manifest.read_bytes()
    calls = []
    jpg = _jpeg()

    def fake_request(root):
        def transport(url):
            calls.append(url)
            return jpg, "image/jpeg", url
        return fetch_remote_rgb_preview(root, transport=transport)

    monkeypatch.setattr(hub_module, "fetch_remote_rgb_preview", fake_request)
    with offline_legacy_viewer():
        hub = hub_module.ProjectHub()
    qtbot.addWidget(hub)
    hub.resize(1300, 820)
    hub.show()
    assert not calls  # no network at application start or selection
    assert hub.rgb_btn.isEnabled()
    assert hub.galaxy_card.preview.width() == hub.galaxy_card.preview.height()
    assert (hub.galaxy_card.preview.mapToGlobal(QPoint(0, 0)).x()
            < hub.galaxy_card.name.mapToGlobal(QPoint(0, 0)).x())
    hub.rgb_action.trigger()
    qtbot.waitUntil(
        lambda: (hub._rgb_future is None and
                 hub.galaxy_card.snapshot is not None and
                 hub.galaxy_card.snapshot.preview_source == "REMOTE_LEGACY_RGB"),
        timeout=10000)
    assert len(calls) == 1
    assert "ls-dr11-south" in calls[0]
    assert "size=384" in calls[0]
    assert "Legacy Surveys Viewer RGB" in hub.galaxy_card.preview.toolTip()
    assert not hub.galaxy_card.preview_note.isVisible()
    assert not hub.galaxy_card.preview.pixmap().isNull()
    assert manifest.read_bytes() == before

    # Entirely offline reopen reads only cached receipt/JPEG.
    def forbidden(*args, **kwargs):
        pytest.fail("Unexpected network on Galaxy Projects reopen")
    monkeypatch.setattr(hub_module, "fetch_remote_rgb_preview", forbidden)
    with offline_legacy_viewer():
        second = hub_module.ProjectHub()
    qtbot.addWidget(second)
    second.resize(1300, 820)
    second.show()
    assert second._rgb_future is None
    assert second.galaxy_card.snapshot.preview_source == "REMOTE_LEGACY_RGB"
    assert not second.galaxy_card.preview.pixmap().isNull()
    assert manifest.read_bytes() == before


def test_failed_colour_download_retains_local_preview_and_no_overwrite(
        qtbot, tmp_path, monkeypatch):
    project = _project(tmp_path)
    fallback = QImage(184, 184, QImage.Format_RGB32)
    fallback.fill(QColor(128, 128, 128))
    cache = project.root / ".desidera" / "card-preview.png"
    cache.parent.mkdir(exist_ok=True)
    assert fallback.save(str(cache), "PNG")
    local_before = cache.read_bytes()
    monkeypatch.setattr(
        hub_module, "fetch_remote_rgb_preview",
        lambda root: (_ for _ in ()).throw(TimeoutError("viewer offline")))
    with offline_legacy_viewer():
        hub = hub_module.ProjectHub()
    qtbot.addWidget(hub)
    hub.show()
    assert hub.galaxy_card.snapshot.thumbnail == cache
    hub.rgb_action.trigger()
    qtbot.waitUntil(lambda: hub._rgb_future is None, timeout=10000)
    assert hub.galaxy_card.snapshot.thumbnail == cache
    assert hub.galaxy_card.snapshot.preview_source == "LOCAL_FITS_GRAYSCALE"
    assert "TimeoutError" in hub.rgb_feedback.text()
    assert cache.read_bytes() == local_before
    assert not (project.root / ".desidera" /
                "card-preview-remote-rgb.jpg").exists()


def test_responsive_preview_remains_square_when_inspector_narrows(qtbot):
    from lsbpipe.gui.galaxy_cards import GalaxyDetailCard
    card = GalaxyDetailCard()
    qtbot.addWidget(card)
    card.resize(420, 560)
    card.show()
    qtbot.waitUntil(lambda: card.preview.width() == 124)
    assert card.preview.height() == 124
    # Empty placeholders stay compact even at large window widths.
    card.resize(760, 650)
    qtbot.waitUntil(lambda: card.preview.width() == 124)
    assert card.preview.height() == 124
    # A real pixel-bearing preview expands to 184 logical pixels.
    from PySide6.QtGui import QPixmap
    from PySide6.QtGui import QColor
    card._pixmap = QPixmap(184, 184)
    card._pixmap.fill(QColor("#526f91"))
    card.resize(750, 650)
    card.resize(760, 650)
    qtbot.waitUntil(lambda: card.preview.width() == 184)
    assert card.preview.height() == 184
