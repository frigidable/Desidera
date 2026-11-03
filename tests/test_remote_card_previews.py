"""Official Legacy Viewer RGB preview: bounded offline cache and provider safety."""
from hashlib import sha256
import json
from pathlib import Path

import pytest

from lsbpipe.surveys.legacy import SkyField
from lsbpipe.surveys.object_projects import (
    create_object_project, read_object_project, save_object_project,
)
from lsbpipe.gui.galaxy_card_model import snapshot_from_record
from lsbpipe.gui.remote_card_previews import (
    SCHEMA, SIZE, STATUS, HOST, request_for_project,
    fetch_remote_rgb_preview, inspect_cached_remote_preview,
    _jpeg_dimensions,
)


def synthetic_jpeg(size=SIZE):
    """Minimal JPEG marker fixture for wire-format/cache contract tests."""
    sof = (b"\x08" + size.to_bytes(2, "big") + size.to_bytes(2, "big")
           + b"\x03\x01\x11\x00\x02\x11\x00\x03\x11\x00")
    block = b"\xff\xc0" + (len(sof) + 2).to_bytes(2, "big") + sof
    return (b"\xff\xd8\xff\xe0\x00\x04AB" + block
            + b"\xff\xda\x00\x08\x01\x01\x00\x00\x3f\x00"
            + b"\x00" * 24 + b"\xff\xd9")


def _setup(tmp_path, *, release="dr11", region="south"):
    p = create_object_project("Leo P", tmp_path / "Leo_P",
                              field=SkyField(155.438012, 18.088025, 4))
    p.selected_dataset = (release, region, "1554p180", "r")
    p.data_root = p.root
    save_object_project(p)
    return p


def _record(p):
    return {"project_id": p.project_id, "root": str(p.root),
            "display_name": p.display_name, "available": True,
            "legacy": False, "last_opened_at": ""}


def test_documented_release_specific_rgb_layers_and_actual_square(tmp_path):
    p = _setup(tmp_path)
    request = request_for_project(p)
    assert request["layer"] == "ls-dr11-south"
    assert request["size_px"] == 384
    assert request["pixel_scale_arcsec"] == pytest.approx(240 / 384)
    assert request["url"].startswith("https://www.legacysurvey.org/viewer/jpeg-cutout?")
    assert "layer=ls-dr11-south" in request["url"]
    assert "size=384" in request["url"]
    assert _jpeg_dimensions(synthetic_jpeg()) == (SIZE, SIZE)
    for release, region, expected in [
        ("dr9", "north", "ls-dr9-north"),
        ("dr9", "south", "ls-dr9-south"),
        ("dr10", "south", "ls-dr10-south"),
    ]:
        p.selected_dataset = (release, region, "1554p180", "r")
        assert request_for_project(p)["layer"] == expected
    p.selected_dataset = ("dr11", "north", "1554p180", "r")
    with pytest.raises(ValueError, match="release/region-specific"):
        request_for_project(p)


def test_offline_start_and_independent_rgb_receipt(tmp_path):
    p = _setup(tmp_path)
    fits = p.root / "original.fits.fz"
    fits.write_bytes(b"UNCHANGED")
    manifest = p.root / "lsbpipe-object.json"
    manifest_original = manifest.read_bytes()
    assert inspect_cached_remote_preview(p) == (None, "REMOTE_RGB_NOT_CACHED")
    assert snapshot_from_record(_record(p)).preview_source == "NONE"
    calls = []

    def pretend(url):
        calls.append(url)
        return synthetic_jpeg(), "image/jpeg", url

    result = fetch_remote_rgb_preview(p.root, transport=pretend)
    assert calls == [request_for_project(p)["url"]]
    assert result["status"] == STATUS
    assert result["schema"] == SCHEMA
    assert result["image_sha256"] == sha256(synthetic_jpeg()).hexdigest()
    jpeg = p.root / ".desidera" / "card-preview-remote-rgb.jpg"
    assert jpeg.read_bytes() == synthetic_jpeg()
    cached, status = inspect_cached_remote_preview(read_object_project(p.root))
    assert cached == jpeg and status == "REMOTE_RGB_CACHED"
    card = snapshot_from_record(_record(p))
    assert card.thumbnail == jpeg
    assert card.preview_source == "REMOTE_LEGACY_RGB"
    assert card.preview_status == status
    assert manifest.read_bytes() == manifest_original
    assert fits.read_bytes() == b"UNCHANGED"


def test_stale_config_corrupt_cache_and_invalid_media_fall_back(tmp_path):
    p = _setup(tmp_path)
    good = synthetic_jpeg()
    fetch_remote_rgb_preview(p.root, transport=lambda url: (good, "image/jpeg", url))
    jpg = p.root / ".desidera" / "card-preview-remote-rgb.jpg"
    jpg.write_bytes(b"not a JPEG")
    assert inspect_cached_remote_preview(p)[0] is None
    assert snapshot_from_record(_record(p)).preview_source != "REMOTE_LEGACY_RGB"
    jpg.write_bytes(good)
    p.field = SkyField(155.438012, 18.088025, 5)
    save_object_project(p)
    assert inspect_cached_remote_preview(p)[1] == "STALE_REMOTE_RGB_CACHE"

    with pytest.raises(ValueError, match="non-JPEG|not a JPEG"):
        fetch_remote_rgb_preview(
            p.root, transport=lambda url: (b"<html>server error</html>",
                                           "image/jpeg", url))
    with pytest.raises(ValueError, match="non-JPEG|not a JPEG"):
        fetch_remote_rgb_preview(
            p.root, transport=lambda url: (b"fake", "image/jpeg", url))
    with pytest.raises(ValueError, match="non-JPEG|not a JPEG|incomplete"):
        fetch_remote_rgb_preview(
            p.root, transport=lambda url: (
                b"\xff\xd8\xff" + b"broken" + b"\xff\xd9", "image/jpeg", url))


def test_fail_closed_on_external_provider_or_unapproved_mime(tmp_path):
    p = _setup(tmp_path)
    payload = synthetic_jpeg()
    with pytest.raises(ValueError, match="Unapproved"):
        fetch_remote_rgb_preview(
            p.root, transport=lambda url: (
                payload, "image/jpeg", "https://example.com/fake.jpg"))
    with pytest.raises(ValueError, match="not a JPEG"):
        fetch_remote_rgb_preview(
            p.root, transport=lambda url: (
                payload, "text/html", url))
    assert inspect_cached_remote_preview(p)[0] is None


def test_card_rgb_cache_never_mutates_local_grayscale(tmp_path):
    p = _setup(tmp_path)
    local = p.root / ".desidera" / "card-preview.png"
    local.parent.mkdir(exist_ok=True)
    local.write_bytes(b"OLD_GRAYSCALE")
    before = local.read_bytes()
    fetch_remote_rgb_preview(
        p.root, transport=lambda url: (synthetic_jpeg(), "image/jpeg", url))
    assert local.read_bytes() == before
    assert snapshot_from_record(_record(p)).preview_source == "REMOTE_LEGACY_RGB"


def test_no_dataset_or_position_rejects_request_without_side_effect(tmp_path):
    p = create_object_project("Empty", tmp_path / "empty")
    with pytest.raises(ValueError, match="Select a target"):
        request_for_project(p)
    assert not (p.root / ".desidera").exists()
