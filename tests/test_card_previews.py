"""P10c local preview tests: real compressed FITS, no science adoption."""
import hashlib
import json
import os
import struct

import numpy as np
import pytest

from lsbpipe.gui.card_previews import (
    PNG_SIG, MAX_SIDE, build_local_preview, inspect_cached_preview,
    has_local_image_source,
)
from lsbpipe.gui.galaxy_card_model import snapshot_from_record
from lsbpipe.surveys.legacy import SkyField, BrickHit, product_for
from lsbpipe.surveys.object_projects import (
    create_object_project, read_object_project, save_object_project,
)


def _real_project(tmp_path, *, bands=("r",), near_edge=False):
    fits = pytest.importorskip("astropy.io.fits")
    WCS = pytest.importorskip("astropy.wcs").WCS
    root = tmp_path / "Leo_P"
    project = create_object_project(
        "Leo P", root, field=SkyField(155.438012, 18.088025, 4))
    project.selected_dataset = ("dr11", "south", "1554p180", "r")
    project.data_root = project.root
    save_object_project(project)
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [35., 400.] if near_edge else [400., 400.]
    wcs.wcs.cdelt = [-1. / 3600, 1. / 3600]
    wcs.wcs.crval = [155.438012, 18.088025]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    yy, xx = np.mgrid[:800, :800]
    gx = 34 if near_edge else 399
    flux = (np.exp(-((xx - gx)**2 + (yy - 399)**2) / (2 * 30**2))
            - 0.007).astype("f4")
    flux[10, 10] = np.nan
    hit = BrickHit("dr11", "south", "1554p180",
                   0., 0., 0., 0., {band: 3 for band in bands})
    for band in bands:
        product = product_for(hit, "image", band)
        path = project.root / product.relative_path
        path.parent.mkdir(parents=True, exist_ok=True)
        fits.HDUList([
            fits.PrimaryHDU(),
            fits.CompImageHDU(data=flux, header=wcs.to_header()),
        ]).writeto(path)
    return project


def _record(project):
    return {
        "project_id": project.project_id,
        "root": str(project.root),
        "display_name": project.display_name, "available": True,
        "legacy": False, "last_opened_at": "",
    }


def test_real_fits_local_thumb_is_bounded_provenance_linked_and_read_only(tmp_path):
    project = _real_project(tmp_path, bands=("g", "r"))
    original = project.root / "lsbpipe-object.json"
    before = original.read_bytes()
    hit = BrickHit("dr11", "south", "1554p180",
                   0, 0, 0, 0, {"r": 3})
    image = project.root / product_for(hit, "image", "r").relative_path
    bytes_before = image.read_bytes()
    assert has_local_image_source(project)
    receipt = build_local_preview(project.root)
    assert receipt["status"] == "DISPLAY_ONLY_NOT_SCIENTIFICALLY_ADOPTED"
    assert receipt["fits_integrity_certified"] is False
    assert receipt["source_band"] == "r"
    assert receipt["partial_display_or_size_cap"] is False
    assert receipt["source_sha256_at_generation"] == hashlib.sha256(bytes_before).hexdigest()
    assert original.read_bytes() == before
    assert image.read_bytes() == bytes_before

    png = project.root / ".desidera" / "card-preview.png"
    assert png.read_bytes().startswith(PNG_SIG)
    width, height = struct.unpack(">II", png.read_bytes()[16:24])
    assert 0 < width <= MAX_SIDE
    assert 0 < height <= MAX_SIDE
    cached, state = inspect_cached_preview(read_object_project(project.root))
    assert cached == png
    assert state == "READY_RECENT_STAT_NOT_FITS_REHASHED"
    snapshot = snapshot_from_record(_record(project))
    assert snapshot.thumbnail == png
    assert snapshot.preview_status == state
    assert original.read_bytes() == before


def test_same_source_stat_can_be_reused_without_refetch_or_fits_decoding(tmp_path, monkeypatch):
    project = _real_project(tmp_path)
    build_local_preview(project.root)
    with monkeypatch.context() as patch:
        patch.setattr("astropy.io.fits.open", lambda *a, **kw:
                      pytest.fail("Cache inspection must not decode FITS"))
        candidate, status = inspect_cached_preview(
            read_object_project(project.root))
        assert candidate.is_file() and status.startswith("READY")
        assert snapshot_from_record(_record(project)).thumbnail == candidate


def test_source_mutation_invalidates_cache_without_auto_certification(tmp_path):
    project = _real_project(tmp_path)
    build_local_preview(project.root)
    hit = BrickHit("dr11", "south", "1554p180", 0, 0, 0, 0, {"r": 1})
    source = project.root / product_for(hit, "image", "r").relative_path
    old = source.stat()
    os.utime(source, ns=(old.st_atime_ns, old.st_mtime_ns + 10_000_000))
    cached, state = inspect_cached_preview(read_object_project(project.root))
    assert cached is None and state == "STALE_SOURCE"
    assert snapshot_from_record(_record(project)).thumbnail is None


def test_cache_corruption_is_not_displayed_as_original(tmp_path):
    project = _real_project(tmp_path)
    build_local_preview(project.root)
    cache = project.root / ".desidera" / "card-preview.png"
    cache.write_bytes(PNG_SIG + b"damaged")
    path, state = inspect_cached_preview(read_object_project(project.root))
    assert path is None and state == "INVALID_CACHE"


def test_config_change_invalidates_cache_without_rewriting_project(tmp_path):
    project = _real_project(tmp_path)
    build_local_preview(project.root)
    project.field = SkyField(155.438012, 18.088025, 3.)
    save_object_project(project)
    path, status = inspect_cached_preview(read_object_project(project.root))
    assert path is None and status == "STALE_CONFIGURATION"


def test_partial_original_frame_is_flagged_not_padded(tmp_path):
    project = _real_project(tmp_path, near_edge=True)
    receipt = build_local_preview(project.root)
    assert receipt["partial_display_or_size_cap"] is True
    assert receipt["crop_pixels"][1] < 240
    assert receipt["status"] == "DISPLAY_ONLY_NOT_SCIENTIFICALLY_ADOPTED"


def test_no_original_image_keeps_project_unchanged_and_creates_no_cache(tmp_path):
    p = create_object_project(
        "No FITS", tmp_path / "No_FITS",
        field=SkyField(40.096208, 1.22925, 4))
    p.selected_dataset = ("dr11", "south", "0401p012", "r")
    p.data_root = p.root
    save_object_project(p)
    manifest = (p.root / "lsbpipe-object.json").read_bytes()
    with pytest.raises(FileNotFoundError, match="No original"):
        build_local_preview(p.root)
    assert not (p.root / ".desidera").exists()
    assert (p.root / "lsbpipe-object.json").read_bytes() == manifest
