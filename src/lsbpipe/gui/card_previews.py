"""P10c: cached, display-only Galaxy Projects thumbnail from original FITS.

Original science IMAGE files are read *only* by an explicit/selected-card
background preview task. Reading the project list remains JSON-only. No
network, resampling of scientific datasets, sky fitting, science adoption,
FITS rewriting, or additional virtual environment is involved.

The preview is a disposable byte-image cache under the existing project's
.desidera directory, not a scientific measurement or QA receipt.
"""
from __future__ import annotations

from datetime import datetime, timezone
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import struct
import tempfile
import zlib

import numpy as np

from lsbpipe.surveys.legacy import BrickHit, product_for
from lsbpipe.surveys.object_projects import read_object_project

SCHEMA = "desidera.card-preview/1"
RENDER_VERSION = 1
MAX_SIDE = 512
MAX_CROP_PIXELS = 1800
MAX_PNG_BYTES = 4 * 1024 * 1024
PNG_SIG = b"\x89PNG\r\n\x1a\n"


def _image_paths(project):
    """Small deterministic candidate list. Never walk a FITS directory."""
    if project.field is None or project.selected_dataset is None:
        return []
    if project.data_root is None or not project.data_root.is_dir():
        return []
    release, region, brick, selected_band = project.selected_dataset
    supported = ("g", "r", "z") if region == "north" else ("g", "r", "i", "z")
    order = list(dict.fromkeys(
        ["r", selected_band, "g", "i", "z"]))
    paths = []
    for base in (project.data_root, project.data_root / "data"):
        for band in order:
            if band not in supported:
                continue
            hit = BrickHit(release, region, brick, 0, 0, 0, 0, {band: 1})
            product = product_for(hit, "image", band)
            target = base / product.relative_path
            if target.is_file() and not target.is_symlink():
                paths.append((target, band))
    return paths


def has_local_image_source(project) -> bool:
    return bool(_image_paths(project))


def _field_record(project) -> dict:
    if project.field is None or project.selected_dataset is None:
        raise ValueError("Galaxy target and selected survey dataset required")
    return {
        "ra_deg": project.field.ra_deg,
        "dec_deg": project.field.dec_deg,
        "width_arcmin": project.field.width_arcmin,
    }


def _chunk(tag: bytes, data: bytes) -> bytes:
    return (struct.pack(">I", len(data)) + tag + data
            + struct.pack(">I", zlib.crc32(tag + data) & 0xffffffff))


def _png_gray(gray: np.ndarray) -> bytes:
    """Write a bounded grayscale PNG without new image dependencies."""
    a = np.asarray(gray)
    if (a.ndim != 2 or a.dtype != np.uint8 or
            not 0 < a.shape[0] <= MAX_SIDE or
            not 0 < a.shape[1] <= MAX_SIDE):
        raise ValueError("Card preview must be a bounded uint8 grayscale image")
    h, w = a.shape
    raw = b"".join(
        b"\x00" + a[row, :].tobytes()
        for row in range(h))
    header = struct.pack(">IIBBBBB", w, h, 8, 0, 0, 0, 0)
    return (PNG_SIG + _chunk(b"IHDR", header)
            + _chunk(b"IDAT", zlib.compress(raw, 6))
            + _chunk(b"IEND", b""))


def _hash_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        while True:
            chunk = stream.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def _write_atomic(path: Path, data: bytes):
    """Atomic cache replacement, with no write access to original images."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.parent.is_symlink() or path.is_symlink():
        raise ValueError("Refusing to write through a symlinked preview cache")
    pending = None
    try:
        with tempfile.NamedTemporaryFile(
                mode="wb", prefix=".preview-", suffix=".tmp",
                dir=path.parent, delete=False) as fp:
            pending = Path(fp.name)
            fp.write(data)
            fp.flush()
            os.fsync(fp.fileno())
        os.replace(pending, path)
    finally:
        if pending is not None:
            pending.unlink(missing_ok=True)


def _receipt_path(project) -> Path:
    return project.root / ".desidera" / "card-preview.json"


def _png_path(project) -> Path:
    return project.root / ".desidera" / "card-preview.png"


def inspect_cached_preview(project) -> tuple[Path | None, str]:
    """Check cache receipt and source stat, without opening original FITS.

    Returns an existing compatible PNG or a refusal/status reason. Checking
    source size/mtime does NOT independently revalidate its SHA-256.
    """
    png, receipt = _png_path(project), _receipt_path(project)
    if not receipt.is_file():
        if png.is_file() and not png.is_symlink():
            return png, "LEGACY_PREVIEW_UNVERIFIED"
        return None, "NOT_GENERATED"
    if receipt.is_symlink() or png.is_symlink():
        return None, "INVALID_CACHE"
    try:
        if receipt.stat().st_size > 32768:
            return None, "INVALID_CACHE"
        item = json.loads(receipt.read_text(encoding="utf-8"))
        if (item["schema"] != SCHEMA or
                item["version"] != RENDER_VERSION or
                item["project_id"] != project.project_id or
                item["field"] != _field_record(project) or
                item["dataset"] != list(project.selected_dataset[:3]) or
                item["status"] != "DISPLAY_ONLY_NOT_SCIENTIFICALLY_ADOPTED"):
            return None, "STALE_CONFIGURATION"
        if project.data_root is None:
            return None, "STALE_SOURCE"
        band = item["source_band"]
        if band not in ("g", "r", "i", "z"):
            return None, "INVALID_CACHE"
        paths = _image_paths(project)
        source = next(
            (path for path, candidate_band in paths
             if candidate_band == band and
             str(path.resolve()) == item["source_path"]), None)
        if source is None:
            return None, "STALE_SOURCE"
        stat = source.stat()
        if (stat.st_size != item["source_size_bytes"] or
                stat.st_mtime_ns != item["source_mtime_ns"]):
            return None, "STALE_SOURCE"
        if not png.is_file() or png.stat().st_size > MAX_PNG_BYTES:
            return None, "INVALID_CACHE"
        image = png.read_bytes()
        if not image.startswith(PNG_SIG) or sha256(image).hexdigest() != item["png_sha256"]:
            return None, "INVALID_CACHE"
        return png, "READY_RECENT_STAT_NOT_FITS_REHASHED"
    except (OSError, ValueError, TypeError, KeyError, AttributeError):
        return None, "INVALID_CACHE"


def build_local_preview(root: str | Path) -> dict:
    """Create exactly one reusable, display-only local preview, on demand."""
    from astropy.io import fits
    from astropy.nddata import Cutout2D
    from astropy.wcs import WCS
    from astropy.wcs.utils import proj_plane_pixel_scales
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    from lsbpipe.surveys.coverage import _science_hdu_index
    from lsbpipe.gui.multiband_preview import display_preview
    from lsbpipe.gui.display_transfer import display_gray

    project = read_object_project(root)
    options = _image_paths(project)
    if not options:
        raise FileNotFoundError(
            "No original Legacy IMAGE FITS for the current dataset; "
            "download/link an image before generating its thumbnail")
    source, band = options[0]
    initial_stat = source.stat()
    with fits.open(source, memmap=True) as hdus:
        idx = _science_hdu_index(hdus, "image")
        hdu = hdus[idx]
        hdu.verify("exception")
        if "CHECKSUM" in hdu.header and hdu.verify_checksum() != 1:
            raise ValueError("IMAGE FITS CHECKSUM invalid")
        if "DATASUM" in hdu.header and hdu.verify_datasum() != 1:
            raise ValueError("IMAGE FITS DATASUM invalid")
        image = hdu.data
        if image is None or image.ndim != 2 or image.dtype.kind not in "iuf":
            raise ValueError("Invalid original IMAGE FITS pixels")
        wcs = WCS(hdu.header).celestial
        if not wcs.has_celestial or wcs.pixel_n_dim != 2:
            raise ValueError("Original IMAGE lacks valid celestial WCS")
        scale = float(np.min(proj_plane_pixel_scales(wcs))) * 3600
        if not math.isfinite(scale) or scale <= 0:
            raise ValueError("Invalid original FITS angular pixel scale")
        requested = max(8, math.ceil(project.field.width_arcmin * 60 / scale))
        side = min(requested, MAX_CROP_PIXELS)
        centre = SkyCoord(project.field.ra_deg * u.deg,
                          project.field.dec_deg * u.deg, frame="icrs")
        cutout = Cutout2D(image, centre, (side, side),
                          wcs=wcs, mode="trim", copy=True)
        data = cutout.data
        if data is None or not np.any(np.isfinite(data)):
            raise ValueError("No finite source pixels in preview field")
        overview, factor = display_preview(data, max_side=MAX_SIDE)
        gray = display_gray(overview, "image", contrast="faint")
        # FITS origin is lower left; PNG origin is upper left.
        png_bytes = _png_gray(np.ascontiguousarray(gray[::-1, :]))
        partial = (data.shape != (side, side) or requested > MAX_CROP_PIXELS)

    final_stat = source.stat()
    if (final_stat.st_size != initial_stat.st_size or
            final_stat.st_mtime_ns != initial_stat.st_mtime_ns):
        raise ValueError("Original FITS changed while building card preview")
    source_sha256 = _hash_file(source)
    verify_stat = source.stat()
    if (verify_stat.st_size != initial_stat.st_size or
            verify_stat.st_mtime_ns != initial_stat.st_mtime_ns):
        raise ValueError("Original FITS changed during preview hashing")
    receipt = {
        "schema": SCHEMA, "version": RENDER_VERSION,
        "project_id": project.project_id,
        "dataset": list(project.selected_dataset[:3]),
        "field": _field_record(project),
        "source_path": str(source.resolve()),
        "source_band": band, "source_size_bytes": initial_stat.st_size,
        "source_mtime_ns": initial_stat.st_mtime_ns,
        "source_sha256_at_generation": source_sha256,
        "png_sha256": sha256(png_bytes).hexdigest(),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "crop_pixels": list(data.shape),
        "display_block_factor": factor,
        "partial_display_or_size_cap": bool(partial),
        "status": "DISPLAY_ONLY_NOT_SCIENTIFICALLY_ADOPTED",
        "fits_integrity_certified": False,
    }
    cache = project.root / ".desidera"
    if cache.is_symlink():
        raise ValueError("Refusing symlinked preview cache directory")
    # Cache files only; never touch original FITS, project manifest or QA.
    _write_atomic(_png_path(project), png_bytes)
    _write_atomic(_receipt_path(project), (
        json.dumps(receipt, sort_keys=True, ensure_ascii=False,
                   indent=2) + "\n").encode("utf-8"))
    return receipt
