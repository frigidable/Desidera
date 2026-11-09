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
                                arcsec_per_pixel, protected, wcs=None):
    """Source-type-aware *screening* footprints, never rendered Tractor models.

    Point sources get conservative round brightness-dependent envelopes.
    Where the original Tractor provides usable galaxy shape_r/e1/e2,
    extended models get ellipse-shaped envelopes on the original WCS
    tangent plane. Missing/degenerate morphology falls back to a circle,
    never an invented PA or q. All model footprints touching user or
    project target protection remain unresolved and are NOT removed.
    """
    if (len(shape) != 2 or protected.shape != tuple(shape)
            or protected.dtype.kind != "b" or
            not math.isfinite(arcsec_per_pixel)
            or not .05 < arcsec_per_pixel < 10.):
        raise ValueError("Invalid native mask geometry")
    ny, nx = shape
    output = np.zeros((ny, nx), dtype=bool)
    accepted, near, outside, ignored = 0, 0, 0, 0
    elliptical, missing_shape = 0, 0
    conflicts = []
    if np.shape(positions) != (len(rows), 2):
        raise ValueError("Tractor positions and row count differ")
    # Original WCS CD/PC matrix maps *native FITS pixels* to east/north
    # arcsec at the small selected field's tangent location.
    jacobian = None
    if wcs is not None:
        cd = np.asarray(wcs.celestial.pixel_scale_matrix, dtype=float)
        dec0 = float(wcs.celestial.wcs.crval[1])
        jacobian = np.array([
            [cd[0, 0]*math.cos(math.radians(dec0))*3600.,
             cd[0, 1]*math.cos(math.radians(dec0))*3600.],
            [cd[1, 0]*3600., cd[1, 1]*3600.]], dtype=float)
        if (jacobian.shape != (2, 2) or
                not np.isfinite(jacobian).all() or
                abs(float(np.linalg.det(jacobian))) < 1e-8):
            raise ValueError("Invalid local WCS for Tractor morphology")
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
        # A provisional magnitude-dependent neighbourhood, not predicted
        # wing flux or a calibrated mask-to-depth relation.
        baseline_arcsec = min(30., max(4., 4.+max(0.,20.-magnitude)*2.2))
        major = minor = baseline_arcsec
        pa = None
        if source_type != "PSF":
            radius = row.get("shape_r_arcsec")
            if radius is not None and math.isfinite(float(radius)) and float(radius) > 0.:
                radius = float(radius)
                major = min(35., max(baseline_arcsec, 3.*radius+3.))
                if source_type != "REX":
                    e1, e2 = row.get("shape_e1"), row.get("shape_e2")
                    if e1 is not None and e2 is not None and all(
                            math.isfinite(float(e)) for e in (e1,e2)):
                        ellipticity = math.hypot(float(e1),float(e2))
                        if 1e-3 < ellipticity < .95 and jacobian is not None:
                            # Official Legacy Tractor convention:
                            # epsilon=(a-b)/(a+b)*exp(2i phi); PA=(180-phi)%180.
                            q = (1.-ellipticity)/(1.+ellipticity)
                            minor = min(major, max(
                                3., 3.*radius*q+3., baseline_arcsec*.55))
                            pa = (-.5*math.atan2(float(e2),float(e1))) % math.pi
                        elif ellipticity >= .95:
                            missing_shape += 1
                    else:
                        missing_shape += 1
            else:
                major = max(major, 7.)
                if source_type != "REX":
                    missing_shape += 1
        # The geometric-mean scale is insufficient for a rotated or
        # anisotropic native WCS: bound *both* axes by the smallest
        # singular value of the local east/north pixel Jacobian.
        min_scale = (float(np.min(np.linalg.svd(jacobian, compute_uv=False)))
                     if pa is not None else arcsec_per_pixel)
        r_pix = major/min_scale
        if (x+r_pix < 0 or x-r_pix >= nx or
                y+r_pix < 0 or y-r_pix >= ny):
            outside += 1
            continue
        x1, x2 = max(0, int(math.floor(x-r_pix))), min(nx, int(math.ceil(x+r_pix))+1)
        y1, y2 = max(0, int(math.floor(y-r_pix))), min(ny, int(math.ceil(y+r_pix))+1)
        yy, xx = np.ogrid[y1:y2, x1:x2]
        dx, dy = xx-x, yy-y
        if pa is not None:
            east = jacobian[0,0]*dx+jacobian[0,1]*dy
            north = jacobian[1,0]*dx+jacobian[1,1]*dy
            along = east*math.sin(pa)+north*math.cos(pa)
            across = east*math.cos(pa)-north*math.sin(pa)
            disk = (along/major)**2+(across/minor)**2 <= 1.
        else:
            disk = dx*dx+dy*dy <= r_pix**2
        if np.any(disk & protected[y1:y2, x1:x2]):
            near += 1
            if len(conflicts) < 128:
                conflicts.append({
                    "x": x, "y": y, "objid": row.get("objid"),
                    "type": source_type, "radius_pix": r_pix,
                })
            continue
        output[y1:y2, x1:x2] |= disk
        accepted += 1
        if pa is not None:
            elliptical += 1
    output.setflags(write=False)
    return output, {
        "masked_catalog_sources": accepted,
        "elliptical_catalog_sources": elliptical,
        "missing_shape_fallbacks": missing_shape,
        "skipped_target_overlap": near,
        "near_target_sources": tuple(conflicts),
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
        arcsec_per_pixel=pxscale, protected=protected, wcs=dataset["wcs"])
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
