"""Optional DESI Legacy PSFSIZE seeing prior, never an adopted ePSF.

Legacy coadd psfsize-<band>.fits.fz gives *weighted average PSF FWHM
in arcseconds at each native coadd pixel*. If an original map is already
present on disk, sample it near the project's ICRS coordinate using
verified map-vs-native-image WCS. Convert to native-image pixels from
actual WCS (do not assume 0.262 arcsec/pix).
Missing/invalid/unsupported maps are an explicit absence, NOT a reason
to silently download or invent an empirical measured seeing.

The galaxy-protection radius is conceptually separate and cannot be
inferred from the PSF seeing map.
"""
from __future__ import annotations

from hashlib import sha256
import math
from pathlib import Path

import numpy as np


def suggest_detection_fwhm_from_psfsize(frame, *, target_icrs_deg):
    from astropy.io import fits
    from astropy.wcs import WCS
    from astropy.wcs.utils import proj_plane_pixel_scales

    files = frame.source_evidence.get("files", {})
    source = files.get("image")
    if not source:
        raise ValueError("Native original FITS provenance is required")
    image_path = Path(source)
    marker = f"-image-{frame.band}.fits.fz"
    if not image_path.name.endswith(marker):
        raise ValueError("Unsupported image path for matching PSFSIZE map")
    path = image_path.with_name(
        image_path.name[:-len(marker)] + f"-psfsize-{frame.band}.fits.fz")
    if not path.is_file() or path.is_symlink():
        return {
            "available": False,
            "source": "OPTIONAL_DESI_LOCAL_PSF_SIZE_MAP_ABSENT",
            "message": ("The per-band PSFSIZE map is not stored in the "
                        "current Leo P native FITS project. Manual FWHM "
                        "guess remains an explicit heuristic."),
        }
    with fits.open(path, memmap=True) as hdul:
        images = [hdu for hdu in hdul if getattr(hdu, "data", None) is not None
                  and getattr(hdu.data, "ndim", None) == 2]
        if len(images) != 1:
            raise ValueError("PSFSIZE map needs exactly one native 2D image")
        hdu = images[0]
        if tuple(hdu.data.shape) != frame.shape:
            raise ValueError("PSFSIZE map shape differs from original IMAGE")
        map_wcs = WCS(hdu.header).celestial
        source_wcs = frame.wcs.celestial
        x, y = source_wcs.all_world2pix(
            float(target_icrs_deg[0]), float(target_icrs_deg[1]), 0)
        if not np.isfinite(x) or not np.isfinite(y):
            raise ValueError("Nonfinite target center in native IMAGE WCS")
        ix, iy = int(np.rint(x)), int(np.rint(y))
        ny, nx = frame.shape
        if ix < 2 or ix >= nx-2 or iy < 2 or iy >= ny-2:
            raise ValueError("Target outside safe PSFSIZE sample footprint")
        probe_x, probe_y = map_wcs.all_world2pix(
            float(target_icrs_deg[0]), float(target_icrs_deg[1]), 0)
        if (not np.isfinite(probe_x) or not np.isfinite(probe_y) or
                np.hypot(probe_x-x, probe_y-y) > .05):
            raise ValueError("PSFSIZE map WCS differs at Leo P target")
        # Check an additional off-center position to reject different scale
        # or orientation despite matching target CRVAL.
        world0, world1 = source_wcs.all_pix2world(
            np.array([x+35., x]), np.array([y, y+35.]), 0)
        test_x, test_y = map_wcs.all_world2pix(world0, world1, 0)
        offset = np.hypot(test_x-np.array([x+35., x]),
                          test_y-np.array([y, y+35.]))
        if not np.all(np.isfinite(offset)) or float(np.max(offset)) > .05:
            raise ValueError("PSFSIZE map WCS scale/orientation mismatch")
        patch = np.asarray(hdu.data[iy-2:iy+3, ix-2:ix+3], dtype=float)
        samples = patch[np.isfinite(patch) & (patch > 0)]
        if len(samples) < 12:
            raise ValueError("PSFSIZE map lacks positive finite data near Leo P")
        arcsec = float(np.median(samples))
    pixscales = np.asarray(proj_plane_pixel_scales(frame.wcs.celestial),
                           dtype=float)*3600.
    if (not np.isfinite(pixscales).all()
            or np.min(pixscales) <= 0 or
            np.max(pixscales)/np.min(pixscales) > 1.05):
        raise ValueError("Unexpected anisotropic/invalid native image WCS scale")
    scale = float(np.mean(pixscales))
    fwhm = arcsec/scale
    if not (0.2 <= arcsec <= 8. and 2.5 <= fwhm <= 20.):
        raise ValueError(
            f"PSFSIZE {arcsec:.3g} arcsec -> {fwhm:.3g} px is outside "
            "the initial DAOStarFinder recipe's validated range")
    digest = sha256()
    with path.open("rb") as reader:
        for block in iter(lambda: reader.read(4*1024*1024), b""):
            digest.update(block)
    return {
        "available": True,
        "source": "LOCAL_LEGACY_PSF_SIZE_MAP_NOT_EMPIRICAL_EPSF",
        "path": str(path),
        "source_sha256_local": digest.hexdigest(),
        "band": frame.band,
        "image_pixel_scale_arcsec": scale,
        "psfsize_arcsec": arcsec,
        "suggested_detection_fwhm_pix": fwhm,
        "sample": "5x5 median at project ICRS coordinate",
        "is_measured_empirical_psf": False,
        "requires_researcher_review": True,
    }
