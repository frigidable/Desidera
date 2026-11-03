"""Manually fetched, display-only colour cutouts from the official Legacy Viewer.

This is a separate *visual cache*, not a science input or a substitute for
original Legacy survey FITS. Never performs I/O during Galaxy Hub startup.
"""
from __future__ import annotations

from datetime import datetime, timezone
from hashlib import sha256
import json
import math
from pathlib import Path
from urllib.parse import urlencode, urlparse
from urllib.request import Request, build_opener, HTTPRedirectHandler

from lsbpipe.surveys.object_projects import read_object_project
from .card_previews import _write_atomic

HOST = "www.legacysurvey.org"
ENDPOINT = f"https://{HOST}/viewer/jpeg-cutout"
SCHEMA = "desidera.remote-rgb-preview/1"
SIZE = 384  # official Viewer maximum is currently 512
MAX_BYTES = 4 * 1024 * 1024
STATUS = "DISPLAY_ONLY_REMOTE_RGB_NOT_SCIENTIFICALLY_ADOPTED"
LAYERS = {
    ("dr9", "north"): "ls-dr9-north",
    ("dr9", "south"): "ls-dr9-south",
    ("dr10", "south"): "ls-dr10-south",
    ("dr11", "south"): "ls-dr11-south",
}
# DR10 and DR11 in the northern footprint reuse DR9 North imagery. Avoid
# pretending this is an independent DR10/DR11 north-specific RGB release.
JPEG_SOFS = {0xC0, 0xC1, 0xC2, 0xC3, 0xC5, 0xC6, 0xC7,
             0xC9, 0xCA, 0xCB, 0xCD, 0xCE, 0xCF}


def request_for_project(project) -> dict:
    if project.field is None or project.selected_dataset is None:
        raise ValueError("Select a target position and a survey dataset first")
    release, region, brick, _ = project.selected_dataset
    layer = LAYERS.get((release, region))
    if layer is None:
        raise ValueError(
            "No release/region-specific official RGB layer is verified; "
            "use the local FITS thumbnail instead")
    ra, dec, width = (float(project.field.ra_deg),
                      float(project.field.dec_deg),
                      float(project.field.width_arcmin))
    if not (math.isfinite(ra) and math.isfinite(dec)
            and math.isfinite(width) and 0 <= ra < 360
            and -90 <= dec <= 90 and 0 < width <= 120):
        raise ValueError("Invalid project position or preview angular field")
    pixscale = round(width * 60 / SIZE, 7)  # arcseconds per pixel
    params = {
        "ra": f"{ra:.8f}", "dec": f"{dec:.8f}",
        "layer": layer, "pixscale": f"{pixscale:.7f}",
        "size": str(SIZE),
    }
    return {
        "project_id": project.project_id,
        "field": {"ra_deg": project.field.ra_deg,
                  "dec_deg": project.field.dec_deg,
                  "width_arcmin": project.field.width_arcmin},
        "dataset": [release, region, brick],
        "layer": layer,
        "pixel_scale_arcsec": pixscale,
        "size_px": SIZE,
        "url": ENDPOINT + "?" + urlencode(params),
    }


def _jpeg_dimensions(raw: bytes) -> tuple[int, int]:
    """Read SOF dimensions without trusting a server's MIME or file suffix."""
    if (len(raw) < 32 or len(raw) > MAX_BYTES or
            not raw.startswith(b"\xff\xd8\xff") or not raw.endswith(b"\xff\xd9")):
        raise ValueError("Legacy Viewer returned an incomplete/non-JPEG image")
    offset = 2
    while offset + 4 <= len(raw):
        if raw[offset] != 0xff:
            raise ValueError("Malformed JPEG segment marker")
        while offset < len(raw) and raw[offset] == 0xff:
            offset += 1
        if offset >= len(raw):
            break
        marker = raw[offset]
        offset += 1
        if marker == 0xda:  # start of compressed scan
            break
        if marker in (0xd8, 0xd9) or 0xd0 <= marker <= 0xd7:
            continue
        if offset + 2 > len(raw):
            break
        length = int.from_bytes(raw[offset:offset + 2], "big")
        if length < 2 or offset + length > len(raw):
            raise ValueError("Malformed JPEG segment length")
        if marker in JPEG_SOFS:
            if length < 8:
                raise ValueError("Truncated JPEG frame header")
            height = int.from_bytes(raw[offset + 3:offset + 5], "big")
            width = int.from_bytes(raw[offset + 5:offset + 7], "big")
            return width, height
        offset += length
    raise ValueError("JPEG image dimensions could not be verified")


class _OfficialRedirectsOnly(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        address = urlparse(newurl)
        if address.scheme != "https" or address.hostname != HOST:
            raise ValueError("RGB cutout redirect outside official HTTPS host")
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def _download(url: str) -> tuple[bytes, str, str]:
    opener = build_opener(_OfficialRedirectsOnly())
    req = Request(
        url, headers={"Accept": "image/jpeg",
                      "User-Agent": "DESIDERA Galaxy Projects visual preview"})
    with opener.open(req, timeout=20) as response:
        final_url = response.geturl()
        parsed = urlparse(final_url)
        if (parsed.scheme != "https" or parsed.hostname != HOST or
                parsed.path not in ("/viewer/jpeg-cutout",
                                    "/viewer/jpeg-cutout/")):
            raise ValueError("Unexpected RGB cutout source URL")
        mime = response.headers.get("Content-Type", "").split(";")[0].lower()
        if mime not in ("image/jpeg", "image/jpg"):
            raise ValueError(f"Legacy Viewer returned {mime or 'unknown'} instead of JPEG")
        data = response.read(MAX_BYTES + 1)
        declared = response.headers.get("Content-Length")
        if declared is not None and int(declared) != len(data):
            raise ValueError("Remote RGB cutout Content-Length mismatch")
        return data, mime, final_url


def _cache_paths(project):
    folder = project.root / ".desidera"
    return (folder / "card-preview-remote-rgb.jpg",
            folder / "card-preview-remote-rgb.json")


def inspect_cached_remote_preview(project) -> tuple[Path | None, str]:
    """Offline/readonly: validate cache against current accepted project field."""
    image, receipt = _cache_paths(project)
    if image.is_symlink() or receipt.is_symlink():
        return None, "INVALID_REMOTE_RGB_CACHE"
    if not receipt.is_file():
        return None, "REMOTE_RGB_NOT_CACHED"
    try:
        expected = request_for_project(project)
        if receipt.stat().st_size > 32768 or not image.is_file():
            return None, "INVALID_REMOTE_RGB_CACHE"
        record = json.loads(receipt.read_text(encoding="utf-8"))
        if (record.get("schema") != SCHEMA or
                record.get("status") != STATUS or
                record.get("request") != expected):
            return None, "STALE_REMOTE_RGB_CACHE"
        data = image.read_bytes() if image.stat().st_size <= MAX_BYTES else b""
        if (sha256(data).hexdigest() != record.get("image_sha256") or
                _jpeg_dimensions(data) != (SIZE, SIZE)):
            return None, "INVALID_REMOTE_RGB_CACHE"
        return image, "REMOTE_RGB_CACHED"
    except (OSError, ValueError, KeyError, TypeError, AttributeError):
        return None, "INVALID_REMOTE_RGB_CACHE"


def fetch_remote_rgb_preview(root: str | Path, *, transport=None) -> dict:
    """Explicit user action only. No automatic network calls on selection."""
    project = read_object_project(root)
    request = request_for_project(project)
    downloader = transport or _download
    raw, mime, final_url = downloader(request["url"])
    final = urlparse(final_url)
    if (final.scheme != "https" or final.hostname != HOST or
            final.path not in ("/viewer/jpeg-cutout", "/viewer/jpeg-cutout/")):
        raise ValueError("Unapproved remote RGB provider")
    if mime.lower().split(";")[0] not in ("image/jpeg", "image/jpg"):
        raise ValueError("Remote cutout is not a JPEG")
    if _jpeg_dimensions(raw) != (SIZE, SIZE):
        raise ValueError("Remote RGB size differs from requested square")
    if request_for_project(read_object_project(root)) != request:
        raise ValueError("Project target/dataset changed during RGB download")
    cache, receipt = _cache_paths(project)
    if cache.parent.is_symlink():
        raise ValueError("Refusing symlinked galaxy preview cache")
    record = {
        "schema": SCHEMA, "request": request,
        "provider": "Legacy Surveys Viewer",
        "response_url": final_url, "response_mime": mime,
        "image_sha256": sha256(raw).hexdigest(),
        "fetched_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": STATUS,
    }
    # Existing caches are disposable, but the original science FITS, project
    # manifest and QA receipts remain read-only.
    _write_atomic(cache, raw)
    _write_atomic(receipt, (
        json.dumps(record, ensure_ascii=False, sort_keys=True, indent=2)
        + "\n").encode("utf-8"))
    return record
