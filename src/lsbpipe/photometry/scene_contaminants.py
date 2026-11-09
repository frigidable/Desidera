"""Provisional catalogue-assisted neighbour masks for Scene & Masks.

Uses same-brick, same-band DR11 Tractor positions and shape metadata as
*diagnostic* evidence, not a rendered Tractor model or validated PSF.
Never changes original FITS, ignores catalogue rows inside the protected
target footprint, and refuses to proceed when catalogue evidence is absent.
"""
from __future__ import annotations

import math
from types import SimpleNamespace

import numpy as np
from astropy.wcs.utils import proj_plane_pixel_scales

from lsbpipe.photometry.psf_tractor_gaia import read_native_tractor_gaia
from lsbpipe.photometry.sky_masks import rasterize_sky_shape
from lsbpipe.survey import named_bit


def conservative_catalogue_mask(shape, positions, rows, *,
                                arcsec_per_pixel, protected):
    """Bounded circular neighbour envelopes, not model subtraction.

    Deliberately skips every source whose *mask footprint* touches the
    protected galaxy. This may leave contamination near the target;
    near-target sources are explicitly counted for human review.
    """
    if (len(shape) != 2 or protected.shape != tuple(shape)
            or protected.dtype.kind != "b" or
            not math.isfinite(arcsec_per_pixel)
            or not .05 < arcsec_per_pixel < 10.):
        raise ValueError("Invalid native mask geometry")
    ny, nx = shape
    output = np.zeros((ny, nx), dtype=bool)
    accepted, near, outside, ignored = 0, 0, 0, 0
    if np.shape(positions) != (len(rows), 2):
        raise ValueError("Tractor positions and row count differ")
    for (x, y), row in zip(positions, rows):
        x, y = float(x), float(y)
        if not math.isfinite(x) or not math.isfinite(y):
            ignored += 1
            continue
        flux = row.get("flux")
        ivar = row.get("flux_ivar")
        if (flux is None or ivar is None or
                not math.isfinite(flux) or not math.isfinite(ivar)
                or flux <= 0 or ivar <= 0 or
                flux*math.sqrt(ivar) < 5.):
            ignored += 1
            continue
        source_type = row.get("type", "")
        if source_type not in ("PSF", "REX", "EXP", "DEV", "SER"):
            ignored += 1
            continue
        magnitude = 22.5 - 2.5*math.log10(flux)
        # A catalogue-only *screening envelope*. Size is not an estimate
        # of a validated stellar PSF wing nor a Tractor render footprint.
        arcsec = min(30., max(4., 4.+max(0., 20.-magnitude)*2.2))
        if source_type != "PSF":
            shape_r = row.get("shape_r_arcsec")
            if shape_r is not None and math.isfinite(float(shape_r)) and shape_r > 0:
                arcsec = min(35., max(arcsec, 3.*float(shape_r)))
            else:
                arcsec = max(arcsec, 7.)
        radius = arcsec/arcsec_per_pixel
        if (x+radius < 0 or x-radius >= nx or
                y+radius < 0 or y-radius >= ny):
            outside += 1
            continue
        x1, x2 = max(0, int(math.floor(x-radius))), min(nx, int(math.ceil(x+radius))+1)
        y1, y2 = max(0, int(math.floor(y-radius))), min(ny, int(math.ceil(y+radius))+1)
        yy, xx = np.ogrid[y1:y2, x1:x2]
        disk = (xx-x)**2+(yy-y)**2 <= radius**2
        if np.any(disk & protected[y1:y2, x1:x2]):
            near += 1
            continue
        output[y1:y2, x1:x2] |= disk
        accepted += 1
    output.setflags(write=False)
    return output, {
        "masked_catalog_sources": accepted,
        "skipped_target_overlap": near,
        "outside_cutout": outside,
        "unusable_catalog_sources": ignored,
    }


def prepare_scene_contamination(dataset, project, draft, anchor_icrs):
    """Prepare source-bound candidate mask in a worker, with no network.

    The current stage supports verified local DR11 South native Tractor
    catalogues. In other datasets or when the catalogue is absent, fail
    closed rather than call a raw threshold contour a cleaned isophote.
    """
    from lsbpipe.photometry.sky_masks import SkyMaskShape
    if (project.project_id != dataset["project_id"] or
            project.selected_dataset is None or
            project.data_root is None or draft.project_id != project.project_id):
        raise ValueError("Scene, draft and original project identities differ")
    release, region, brick, _ = project.selected_dataset
    band = dataset["band"]
    if (release, region) != ("dr11", "south"):
        raise ValueError(
            "Catalogued-neighbour masking currently supports DR11 South only; "
            "no unmasked isophote will be generated.")
    frame = SimpleNamespace(
        dataset_id=f"{release}/{region}/{brick}/{band}",
        band=band, wcs=dataset["wcs"])
    catalog = read_native_tractor_gaia(
        project.data_root, release=release, region=region,
        brick=brick, frame=frame)
    shape = tuple(dataset["pixel_shape"])
    hard = np.asarray(dataset["native_hard"], dtype=bool)
    if hard.shape != shape:
        raise ValueError("Native hard mask shape differs from image")
    scales = np.asarray(proj_plane_pixel_scales(dataset["wcs"].celestial))
    if scales.shape != (2,) or not np.isfinite(scales).all():
        raise ValueError("Invalid celestial WCS pixel scale")
    pxscale = float(math.sqrt(scales[0]*scales[1])*3600.)
    anchor = tuple(float(x) for x in anchor_icrs)
    protected = rasterize_sky_shape(
        SkyMaskShape("circle", center=anchor,
                     radius_arcsec=float(project.target_radius_arcsec)),
        dataset["wcs"], shape).astype(bool)
    manual_excluded = np.zeros(shape, dtype=bool)
    restored_auto = np.zeros(shape, dtype=bool)
    for edit in draft.edits:
        if not edit.applies_to(band):
            continue
        if edit.purpose not in ("protect_target", "exclude", "restore_auto"):
            continue
        selected = rasterize_sky_shape(edit.shape, dataset["wcs"], shape) > 0.
        if edit.purpose == "protect_target":
            protected |= selected
        elif edit.purpose == "exclude":
            manual_excluded |= selected
        else:
            restored_auto |= selected
    auto, counts = conservative_catalogue_mask(
        shape, catalog.pixel_positions, catalog.rows,
        arcsec_per_pixel=pxscale, protected=protected)
    auto = auto & ~restored_auto
    # Explicit manual exclusions are applied even if they overlap target
    # protection, whereas automated catalogue proposals never do.
    excluded = auto | manual_excluded
    usable = ~hard & ~excluded
    # QA warning bits are NOT all automatically invalid pixels: BRIGHT
    # marks a blob; GALAXY may refer to the target and must be protected.
    bits = dataset["data"]["maskbits"]
    caution = {key: int(np.count_nonzero(named_bit(bits, key)))
               for key in ("BRIGHT", "MEDIUM", "GALAXY", "BAILOUT")}
    if np.count_nonzero(usable) < 500:
        raise ValueError("Prepared source mask leaves too few usable pixels")
    flags = ["CATALOGUE_GEOMETRY_NOT_TRACTOR_MODEL_RENDER",
             "PSF_WINGS_NOT_VALIDATED",
             "PROVISIONAL_PROTECTED_TARGET_RADIUS"]
    if counts["skipped_target_overlap"]:
        flags.append("TARGET_NEIGHBOURS_REQUIRE_MANUAL_REVIEW")
    if caution["BAILOUT"]:
        flags.append("TRACTOR_BAILOUT_PIXELS_PRESENT")
    if caution["BRIGHT"]:
        flags.append("BRIGHT_STAR_MASKBITS_REQUIRE_WING_REVIEW")
    output = {
        "status": "CATALOGUE_CONTAMINATION_MASK_NOT_ADOPTED",
        "mask": excluded.copy(),
        "auto_mask": auto.copy(),
        "protected": protected.copy(),
        "hard": hard.copy(),
        "project_id": project.project_id,
        "draft_digest": draft.digest,
        "dataset_source": dataset["source"],
        "band": band,
        "anchor_icrs": anchor,
        "catalog_sha256": catalog.source_sha256,
        "catalog_path": catalog.source_path,
        "stats": {**counts,
                  "catalog_rows": len(catalog.rows),
                  "excluded_pixels": int(np.count_nonzero(excluded)),
                  "native_hard_pixels": int(np.count_nonzero(hard)),
                  "protected_pixels": int(np.count_nonzero(protected)),
                  "usable_pixels": int(np.count_nonzero(usable))},
        "maskbits_caution_counts": caution,
        "flags": tuple(flags),
    }
    for name in ("mask", "auto_mask", "protected", "hard"):
        output[name].setflags(write=False)
    return output
