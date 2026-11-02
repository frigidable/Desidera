"""Read-only, header-only FITS sky footprints for Dataset Manager.

This reports *geometric* coverage, not good pixels, photometric coverage,
background quality, PSF, matched exposures, mosaics, or scientific adoption.
No FITS data arrays are read or modified. WCS maps remain native.
"""
from __future__ import annotations

from pathlib import Path
import math

from .dataset_registry import read_dataset_index, resolve_asset_path

MAX_DRAWN_FOOTPRINTS = 240
EDGE_SAMPLES = 9


def _asset_footprint(root: Path, asset: dict, center):
    from astropy.io import fits
    from astropy.wcs import WCS
    from astropy.coordinates import SkyCoord, SkyOffsetFrame
    import astropy.units as u
    import numpy as np

    path = resolve_asset_path(root, asset)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError("FITS file missing or empty")
    with fits.open(path, memmap=False, lazy_load_hdus=True) as hdus:
        selected = None
        for hdu in hdus:
            header = hdu.header
            nx = int(header.get("NAXIS1") or header.get("ZNAXIS1") or 0)
            ny = int(header.get("NAXIS2") or header.get("ZNAXIS2") or 0)
            if nx < 2 or ny < 2:
                continue
            candidate = WCS(header).celestial
            if candidate.has_celestial and candidate.pixel_n_dim == 2:
                selected = (candidate, nx, ny)
                break
        if selected is None:
            raise ValueError("No 2D celestial FITS WCS header")
        wcs, nx, ny = selected
        # Trace sampled pixel-boundary edges, not an invented square brick.
        edge = np.linspace(-0.5, nx - 0.5, EDGE_SAMPLES)
        side = np.linspace(-0.5, ny - 0.5, EDGE_SAMPLES)
        xs = np.concatenate([edge, np.full(EDGE_SAMPLES, nx-.5),
                             edge[::-1], np.full(EDGE_SAMPLES, -.5)])
        ys = np.concatenate([np.full(EDGE_SAMPLES, -.5), side,
                             np.full(EDGE_SAMPLES, ny-.5), side[::-1]])
        ra, dec = wcs.all_pix2world(xs, ys, 0)
        if not (np.all(np.isfinite(ra)) and np.all(np.isfinite(dec))):
            raise ValueError("Nonfinite FITS sky footprint WCS")
        sky = SkyCoord(ra*u.deg, dec*u.deg, frame="icrs")
        offset = sky.transform_to(SkyOffsetFrame(origin=center))
        east_arcsec = offset.lon.to_value(u.arcsec)
        north_arcsec = offset.lat.to_value(u.arcsec)
        if not (np.all(np.isfinite(east_arcsec)) and
                np.all(np.isfinite(north_arcsec))):
            raise ValueError("Invalid tangent-plane footprint")
        # East left, north up; scene unit is an arcsecond on the sky.
        return {"points": [
            [-float(east), -float(north)]
            for east, north in zip(east_arcsec, north_arcsec)],
            "native_shape": [ny, nx],
            "geometry_status": "FITS_WCS_HEADER_ONLY",
        }


def survey_coverage_footprints(project_root, *, max_drawn=MAX_DRAWN_FOOTPRINTS,
                               checkpoint=None):
    """Inspect bounded number of original IMAGE headers with explicit omissions.

    No guessed brick geometry. Each source's original file status is retained,
    invalid WCS reported per asset; incomplete provenance is never certified.
    """
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from .object_projects import read_object_project
    root = Path(project_root).expanduser().resolve()
    project = read_object_project(root)
    if project.field is None:
        return {"status": "NO_OBJECT_COORDINATES", "footprints": [],
                "rejected": [], "omitted_count": 0}
    index = read_dataset_index(root)
    sources = {source["source_id"]: source for source in index["sources"]}
    images = [a for a in index["assets"]
              if a.get("product_kind") == "image"]
    if not isinstance(max_drawn, int) or max_drawn < 1:
        raise ValueError("Positive coverage display limit required")
    center = SkyCoord(project.field.ra_deg*u.deg,
                      project.field.dec_deg*u.deg, frame="icrs")
    footprints, rejected = [], []
    for asset in images[:max_drawn]:
        if checkpoint is not None:
            checkpoint()
        try:
            footprint = _asset_footprint(root, asset, center)
            source = sources[asset["source_id"]]
            footprints.append({
                "asset_id": asset["asset_id"], "source_id": asset["source_id"],
                "provider": source["provider"], "release": source.get("release"),
                "region": source.get("region"),
                "brick": asset.get("brick"), "band": asset.get("band"),
                **footprint})
        except (OSError, ValueError, KeyError, TypeError) as exc:
            rejected.append({"asset_id": asset["asset_id"],
                             "reason": f"{type(exc).__name__}: {exc}"})
    return {
        "status": "GEOMETRY_ONLY_NOT_ADOPTED",
        "field_width_arcmin": project.field.width_arcmin,
        "target_radius_arcsec": project.target_radius_arcsec,
        "sky_inner_arcsec": project.sky_inner_arcsec,
        "sky_outer_arcsec": project.sky_outer_arcsec,
        "footprints": footprints, "rejected": rejected,
        "omitted_count": max(0, len(images) - max_drawn),
        "warning": "Native FITS WCS header geometry only; no pixel-mask,"
                   " depth, PSF, sky or science mosaic validation",
    }
