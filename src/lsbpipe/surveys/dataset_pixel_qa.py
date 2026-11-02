"""P0.5d: bounded, read-only Legacy per-band pixel coverage diagnostics.

Geometry is assessed separately from FITS/WCS validity and from LSB sky/PSF
suitability. Only the *same source_id / brick / filter* may form a quartet.
A shared MASKBITS applies only to that same source/brick. Nothing is coadded
or silently substituted across DR/regions. Coadd INVVAR is not an independent
noise or SBF-covariance estimate.

Official DR9 MASKBITS semantics: https://www.legacysurvey.org/dr9/bitmasks/
"""
from __future__ import annotations

from contextlib import ExitStack
import math

import numpy as np
from .dataset_registry import read_dataset_index, resolve_asset_path
from .dataset_integrity import _validate_open_hdus
from .task_control import checkpoint

PIXEL_QA_SCHEMA = "lsbpipe.pixel-coverage-diagnostic/1"
MAX_PIXEL_QA_QUARTETS = 16
MAX_PIXEL_QA_PIXELS = 20_000_000
STRIPE_ROWS = 128
REGIONS = ("target", "sky_annulus")
BANDS = ("g", "r", "i", "z")
MASK_FLAGS = {
    "NPRIMARY": 0,
    "BRIGHT": 1,
    "SATUR_G": 2, "SATUR_R": 3, "SATUR_Z": 4,
    "ALLMASK_G": 5, "ALLMASK_R": 6, "ALLMASK_Z": 7,
    "BAILOUT": 10, "MEDIUM": 11,
    "SATUR_I": 14, "ALLMASK_I": 15,
}


def _float_fraction(count, denominator):
    return count / denominator if denominator else None


def _metrics_template():
    return {
        region: {
            "sampled_pixels": 0,
            "finite_image": 0,
            "positive_invvar": 0,
            "positive_nexp": 0,
            "numerical_triplet": 0,
            "flagged_nprimary": 0,
            "flagged_bright": 0,
            "flagged_medium": 0,
            "flagged_bailout": 0,
            "flagged_band_satur": 0,
            "flagged_band_allmask": 0,
            "numerical_triplet_without_band_satur_allmask": 0,
        } for region in REGIONS
    }


def _safe_original_file(root, asset):
    path = resolve_asset_path(root, asset)
    if not path.is_file():
        raise FileNotFoundError("Required indexed original FITS missing: " + str(path))
    stat = path.stat()
    if stat.st_size == 0:
        raise ValueError("Required original FITS is empty")
    if (asset.get("size_bytes") is not None and
            stat.st_size != asset["size_bytes"]):
        raise ValueError("Indexed original FITS size changed; reindex source")
    if (asset.get("mtime_ns") is not None and
            stat.st_mtime_ns != asset["mtime_ns"]):
        raise ValueError("Indexed original FITS timestamp changed; reindex source")
    if asset.get("checksum_status") == "CHANGED_SINCE_CHECKSUM":
        raise ValueError("Changed source checksum; reimport a new asset version")
    return path


def _one_quartet(root, assets, project, band, callback):
    from astropy.io import fits
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    files = {kind: _safe_original_file(root, record)
             for kind, record in assets.items()}
    pre = {kind: (path.stat().st_size, path.stat().st_mtime_ns)
           for kind, path in files.items()}
    with ExitStack() as stack:
        hdus = {kind: stack.enter_context(
            fits.open(path, memmap=True, do_not_scale_image_data=False))
            for kind, path in files.items()}
        arrays, wcs, indices, _units = _validate_open_hdus(hdus, band)
        image, invvar, nexp, maskbits = (
            arrays[kind] for kind in ("image", "invvar", "nexp", "maskbits"))
        ny, nx = image.shape
        if image.size > MAX_PIXEL_QA_PIXELS:
            raise ValueError("Original FITS exceeds bounded pixel-QA resource limit")
        metrics = _metrics_template()
        center = SkyCoord(project.field.ra_deg * u.deg,
                          project.field.dec_deg * u.deg, frame="icrs")
        target = float(project.target_radius_arcsec)
        sky_inner = float(project.sky_inner_arcsec)
        sky_outer = float(project.sky_outer_arcsec)
        if not 0 < target <= sky_inner < sky_outer:
            raise ValueError("Invalid galaxy/sky aperture geometry")
        x = np.arange(nx, dtype=np.float64)
        flag_bit = MASK_FLAGS[f"SATUR_{band.upper()}"]
        all_bit = MASK_FLAGS[f"ALLMASK_{band.upper()}"]
        for top in range(0, ny, STRIPE_ROWS):
            checkpoint(callback)
            bottom = min(ny, top + STRIPE_ROWS)
            y, xx = np.broadcast_arrays(
                np.arange(top, bottom, dtype=np.float64)[:, None], x[None, :])
            ra, dec = wcs.all_pix2world(xx, y, 0)
            sky = SkyCoord(ra * u.deg, dec * u.deg, frame="icrs")
            radius = sky.separation(center).arcsec
            img = np.asarray(image[top:bottom])
            iv = np.asarray(invvar[top:bottom])
            ex = np.asarray(nexp[top:bottom])
            # Preserve high bits of signed int16 MASKBITS (do not treat
            # negative two's-complement values as absent flags).
            bits = np.asarray(maskbits[top:bottom]).astype(np.uint16).astype(
                np.uint32)
            finite = np.isfinite(img)
            ivgood = np.isfinite(iv) & (iv > 0)
            exgood = np.isfinite(ex) & (ex > 0)
            triplet = finite & ivgood & exgood
            flag = lambda bit: (bits & np.uint32(1 << bit)) != 0
            nprimary = flag(MASK_FLAGS["NPRIMARY"])
            bright = flag(MASK_FLAGS["BRIGHT"])
            medium = flag(MASK_FLAGS["MEDIUM"])
            bailout = flag(MASK_FLAGS["BAILOUT"])
            saturated = flag(flag_bit)
            allmasked = flag(all_bit)
            scopes = (
                ("target", np.isfinite(radius) & (radius <= target)),
                ("sky_annulus", np.isfinite(radius) &
                 (radius >= sky_inner) & (radius <= sky_outer)))
            for name, region in scopes:
                m = metrics[name]
                m["sampled_pixels"] += int(region.sum())
                for key, value in (
                    ("finite_image", finite), ("positive_invvar", ivgood),
                    ("positive_nexp", exgood), ("numerical_triplet", triplet),
                    ("flagged_nprimary", nprimary),
                    ("flagged_bright", bright),
                    ("flagged_medium", medium),
                    ("flagged_bailout", bailout),
                    ("flagged_band_satur", saturated),
                    ("flagged_band_allmask", allmasked),
                    ("numerical_triplet_without_band_satur_allmask",
                     triplet & ~saturated & ~allmasked)):
                    m[key] += int(np.count_nonzero(region & value))
        for kind, path in files.items():
            after = path.stat()
            if pre[kind] != (after.st_size, after.st_mtime_ns):
                raise RuntimeError("Original FITS changed during pixel QA")
    normalized = {}
    for name, m in metrics.items():
        denominator = m["sampled_pixels"]
        normalized[name] = {
            "sampled_pixels": denominator,
            "fractions": {
                key: _float_fraction(value, denominator)
                for key, value in m.items() if key != "sampled_pixels"
            },
        }
    return {"status": "PIXEL_DIAGNOSTIC_NOT_ADOPTED",
            "shape": [ny, nx], "regions": normalized,
            "input_asset_ids": {
                kind: record["asset_id"] for kind, record in assets.items()},
            "flags_are_individual_counts_not_default_rejection": True}


def inspect_pixel_coverage(project_root, *, max_quartets=MAX_PIXEL_QA_QUARTETS,
                           checkpoint_fn=None):
    """Returns independent, per-source/brick/band quantitative QA receipts.

    Missing products are BLOCKED, never borrowed from North/South. Reads are
    bounded; no silent result downsampling, no scientific source ranking.
    """
    from .object_projects import read_object_project
    if not isinstance(max_quartets, int) or not 1 <= max_quartets <= MAX_PIXEL_QA_QUARTETS:
        raise ValueError("Invalid bounded quartet limit")
    project = read_object_project(project_root)
    if project.field is None:
        raise ValueError("ICRS project coordinates needed for pixel QA")
    index = read_dataset_index(project.root)
    sources = {source["source_id"]: source for source in index["sources"]}
    by_group = {}
    for asset in index["assets"]:
        if sources[asset["source_id"]].get("provider") != "legacy_surveys":
            continue
        brick = asset.get("brick")
        if not brick:
            continue
        source = sources[asset["source_id"]]
        if asset.get("product_kind") == "maskbits":
            group = (asset["source_id"], brick, None)
        else:
            group = (asset["source_id"], brick, asset.get("band"))
        by_group.setdefault(group, []).append(asset)
    keys = sorted(
        (sid, brick, band) for sid, brick, band in by_group
        if band in BANDS)
    output = []
    for sid, brick, band in keys[:max_quartets]:
        checkpoint(checkpoint_fn)
        source = sources[sid]
        assets = {}
        for a in by_group[(sid, brick, band)] + by_group.get((sid, brick, None), []):
            kind = a.get("product_kind")
            if kind in assets:
                raise ValueError("Ambiguous duplicate original product: "+str(
                    (sid, brick, band, kind)))
            assets[kind] = a
        missing = sorted({"image", "invvar", "nexp", "maskbits"} - assets.keys())
        row = {"source_id": sid, "release": source.get("release"),
               "region": source.get("region"), "brick": brick, "band": band}
        if missing:
            row.update(status="INCOMPLETE_FILES", missing=missing)
        else:
            try:
                row.update(_one_quartet(
                    project.root, assets, project, band, checkpoint_fn))
            except Exception as exc:
                # Never invent valid coverage when a FITS/WCS check fails.
                from .task_control import OperationCancelled
                if isinstance(exc, OperationCancelled):
                    raise
                row.update(status="INVALID_INPUT", error=f"{type(exc).__name__}: {exc}")
        output.append(row)
    return {
        "schema": PIXEL_QA_SCHEMA,
        "status": "DIAGNOSTIC_NOT_ADOPTED",
        "rows": output,
        "omitted_quartets": max(0, len(keys)-max_quartets),
        "policy": (
            "Fractions are per native WCS-selected target and annulus pixels. "
            "Finite IMAGE, positive INVVAR/NEXP, and each MASKBITS condition "
            "are reported separately. A bit flag is NOT automatically rejected. "
            "No sky-background, PSF, pixel covariance, SBF, or photometric "
            "suitability certification; original FITS remain unchanged."),
    }
