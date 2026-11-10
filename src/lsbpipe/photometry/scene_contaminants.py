"""Provisional catalogue-assisted neighbour masks for Scene & Masks.

Uses same-brick, same-band DR11 Tractor positions and shape metadata as
*diagnostic* evidence, not a rendered Tractor model or validated PSF.
Never changes original FITS. Preliminary target protection is a
source-screening VETO only, not physical membership. Explicit manual
exclusions can work inside target protection, but automatic foreground
PSF/galaxy decomposition is NOT implemented. Refuse absent evidence.
"""
from __future__ import annotations

import math
from types import SimpleNamespace

import numpy as np
from astropy.wcs.utils import proj_plane_pixel_scales
from scipy.spatial import cKDTree

from lsbpipe.photometry.psf_tractor_gaia import read_native_tractor_gaia
from lsbpipe.photometry.scene_source_inventory import (
    build_scene_source_inventory, gaia_astrometric_foreground_review,
    conservative_red_extended_candidate,
)
from lsbpipe.photometry.sky_masks import rasterize_sky_shape
from lsbpipe.survey import named_bit


def _has_gaia_catalogue_reference(row):
    """An identified Gaia DR2/DR3 reference, not merely a placeholder mag.

    Native Tractor FITS stores gaia_g_mag=0.0 for sources without Gaia
    references. The non-None test mislabeled *every* such source as a
    Gaia counterpart in SC1 QA. IDs denote catalogue associations,
    not verified astrometry or foreground membership.
    """
    ref_cat = str(row.get("ref_cat") or "").strip().upper()
    if ref_cat not in ("G2", "G3"):
        return False
    ref_id = row.get("ref_id")
    try:
        return ref_id is not None and int(ref_id) > 0
    except (TypeError, ValueError, OverflowError):
        return False


def conservative_catalogue_mask(shape, positions, rows, *,
                                arcsec_per_pixel, protected, wcs=None,
                                sga_core=None,
                                gaia_inner_preview=False,
                                red_extended_preview=False,
                                anchor_pixel=None, target_radius_arcsec=None):
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
    if sga_core is not None and (
            np.shape(sga_core) != tuple(shape) or
            np.asarray(sga_core).dtype.kind != "b"):
        raise ValueError("SGA core region must match native Boolean grid")
    ny, nx = shape
    output = np.zeros((ny, nx), dtype=bool)
    accepted, near, outside, ignored = 0, 0, 0, 0
    elliptical, missing_shape = 0, 0
    possible_shreds = 0
    inner_preview = np.zeros((ny, nx), dtype=bool)
    inner_preview_centres = []
    red_preview = np.zeros((ny,nx), dtype=bool)
    red_preview_centres = []
    conflicts = []
    masked_centres = []
    inventory_truncated = 0
    if np.shape(positions) != (len(rows), 2):
        raise ValueError("Tractor positions and row count differ")
    # Nearest-neighbour isolation is an independent veto. Do not use
    # red/extended heuristics to discard blended target light in dense
    # galaxy cores. Only compute this layer if explicitly requested.
    near_pix = None
    if red_extended_preview:
        if (anchor_pixel is None or len(anchor_pixel) != 2 or
                not all(math.isfinite(float(v)) for v in anchor_pixel) or
                target_radius_arcsec is None or
                not math.isfinite(float(target_radius_arcsec)) or
                float(target_radius_arcsec) <= 0):
            raise ValueError("Red-extended QA needs current target geometry")
        p = np.asarray(positions, dtype=float)
        good = np.where(np.isfinite(p).all(axis=1))[0]
        near_pix = np.full(len(rows), np.inf)
        if len(good) > 1:
            near_pix[good] = cKDTree(p[good]).query(p[good], k=2)[0][:, 1]
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
    for index, ((x, y), row) in enumerate(zip(positions, rows)):
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
            # Red extended interlopers need three independent signals:
            # genuine extended fit, multi-band red SED and isolation.
            # Cutoffs are conservative HEURISTICS pending cross-galaxy
            # purity/galaxy-light loss validation; never science adoption.
            is_red_preview = False
            red_evidence = None
            if red_extended_preview:
                offset = (math.hypot(x-float(anchor_pixel[0]),
                                     y-float(anchor_pixel[1])) *
                          arcsec_per_pixel)
                is_red_preview, red_evidence = (
                    conservative_red_extended_candidate(
                        row, offset_arcsec=offset,
                        neighbour_arcsec=near_pix[index]*arcsec_per_pixel,
                        target_radius_arcsec=target_radius_arcsec))
            if is_red_preview:
                output[y1:y2, x1:x2] |= disk
                red_preview[y1:y2, x1:x2] |= disk
                accepted += 1
                candidate = {
                    "x":x, "y":y, "objid":row.get("objid"),
                    "type":source_type, "radius_pix":r_pix,
                    "gaia_counterpart":False,
                    "diagnostic_only":True, "rule":"red_extended_v1",
                    "evidence":red_evidence,
                }
                red_preview_centres.append(candidate)
                if len(masked_centres) < 4096:
                    masked_centres.append(candidate)
                else:
                    inventory_truncated += 1
                if pa is not None:
                    elliptical += 1
                continue
            # Explicit opt-in DIAGNOSTIC only. Strong reference astrometry
            # nominates possible foreground stars, not verified PSF models.
            # Even here we require PSF-type, valid Gaia G3 ID, and no
            # Gaia duplicate/strongly negative parallax quality alert.
            # Preserve the *entire* geometric envelope as a candidate
            # exclusion, not as a subtraction or science mask.
            is_gaia_preview = (
                gaia_inner_preview and source_type == "PSF"
                and not row.get("gaia_duplicated_source", False)
                and gaia_astrometric_foreground_review(row)[0])
            if is_gaia_preview:
                output[y1:y2, x1:x2] |= disk
                inner_preview[y1:y2, x1:x2] |= disk
                accepted += 1
                candidate = {
                    "x": x, "y": y, "objid": row.get("objid"),
                    "type": source_type, "radius_pix": r_pix,
                    "gaia_counterpart": True,
                    "diagnostic_only": True,
                }
                inner_preview_centres.append(candidate)
                if len(masked_centres) < 4096:
                    masked_centres.append(candidate)
                else:
                    inventory_truncated += 1
                continue
            near += 1
            in_sga = (sga_core is not None and
                      0 <= int(round(x)) < nx and
                      0 <= int(round(y)) < ny and
                      bool(sga_core[int(round(y)), int(round(x))]))
            has_gaia = _has_gaia_catalogue_reference(row)
            # SGA priors and DR11 Tractor fits share input imagery; this
            # is a CANDIDATE explanation, never probabilistic evidence.
            # Gaia counterparts require individual stellar assessment:
            # they are not proof of foreground contamination.
            classification = (
                "POSSIBLE_GALAXY_SHRED" if in_sga and not has_gaia
                else "GAIA_COUNTERPART_REVIEW" if has_gaia
                else "UNRESOLVED_TARGET_BLEND")
            if classification == "POSSIBLE_GALAXY_SHRED":
                possible_shreds += 1
            if len(conflicts) < 256:
                conflicts.append({
                    "x": x, "y": y, "objid": row.get("objid"),
                    "type": source_type, "radius_pix": r_pix,
                    "classification": classification,
                    "sga_core_overlap": bool(in_sga),
                    "ref_cat": row.get("ref_cat"),
                    "ref_id": row.get("ref_id"),
                })
            continue
        output[y1:y2, x1:x2] |= disk
        accepted += 1
        if len(masked_centres) < 4096:
            masked_centres.append({
                "x":x,"y":y,"objid":row.get("objid"),
                "type":source_type,
                "gaia_counterpart":_has_gaia_catalogue_reference(row),
            })
        else:
            inventory_truncated += 1
        if pa is not None:
            elliptical += 1
    output.setflags(write=False)
    inner_preview.setflags(write=False)
    red_preview.setflags(write=False)
    return output, {
        "masked_catalog_sources": accepted,
        "gaia_inner_preview_mask": inner_preview,
        "gaia_inner_preview_sources": len(inner_preview_centres),
        "gaia_inner_preview_centres": tuple(inner_preview_centres),
        "red_extended_preview_mask": red_preview,
        "red_extended_preview_sources": len(red_preview_centres),
        "red_extended_preview_centres": tuple(red_preview_centres),
        "elliptical_catalog_sources": elliptical,
        "missing_shape_fallbacks": missing_shape,
        "skipped_target_overlap": near,
        "possible_sga_shreds": possible_shreds,
        # Catalogue source CENTRES, not model footprints. Use only as
        # unresolved potential galaxy-component inventory for QA after
        # an independently measured initial radial extent exists.
        "masked_catalog_centres": tuple(masked_centres),
        "masked_centres_truncated": inventory_truncated,
        "near_target_sources": tuple(conflicts),
        "outside_cutout": outside,
        "unusable_catalog_sources": ignored,
    }



def prepare_scene_contamination(dataset, project, draft, anchor_icrs,
                                *, sga_prior=None, gaia_inner_preview=False,
                                red_extended_preview=False):
    """Prepare source-bound candidate mask in a worker, with no network.

    The current stage supports verified local DR11 South native Tractor
    catalogues. In other datasets or when the catalogue is absent, fail
    closed rather than call a raw threshold contour a cleaned isophote.
    """
    from lsbpipe.photometry.sky_masks import (
        SkyMaskShape, is_legacy_isophote_protection,
    )
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
    sga_core = None
    if sga_prior is not None:
        from lsbpipe.photometry.sga_priors import SGAGeometryPrior
        if (not isinstance(sga_prior, SGAGeometryPrior) or
                not isinstance(sga_prior.source_sha256, str) or
                len(sga_prior.source_sha256) != 64 or
                sga_prior.core_shape() is None):
            raise ValueError("Valid selected SGA geometry prior required")
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        separation = SkyCoord(*anchor, unit="deg").separation(
            SkyCoord(sga_prior.ra_deg, sga_prior.dec_deg, unit="deg")
        ).arcsec
        if separation > 180.:
            raise ValueError("SGA association exceeds 180 arcsec; review match")
        sga_core = rasterize_sky_shape(
            sga_prior.core_shape(), dataset["wcs"], shape).astype(bool)
        if not np.any(sga_core):
            raise ValueError("Chosen SGA prior outside native FITS cutout")
    protected = rasterize_sky_shape(
        SkyMaskShape("circle", center=anchor,
                     radius_arcsec=float(project.target_radius_arcsec)),
        dataset["wcs"], shape).astype(bool)
    if sga_core is not None:
        # Extends protection of plausible TARGET light, never substitutes
        # SGA D26 for a final galaxy radius/outer isophote.
        protected |= sga_core
    manual_excluded = np.zeros(shape, dtype=bool)
    restored_auto = np.zeros(shape, dtype=bool)
    legacy_isophote_count = 0
    for edit in draft.edits:
        if not edit.applies_to(band):
            continue
        if is_legacy_isophote_protection(edit):
            # Historical SC1 candidate previews were serialized as
            # protect_target. Fail-safe: never let them change *any*
            # subsequent source-screening decision.
            legacy_isophote_count += 1
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
        arcsec_per_pixel=pxscale, protected=protected,
        wcs=dataset["wcs"], sga_core=sga_core,
        gaia_inner_preview=bool(gaia_inner_preview),
        red_extended_preview=bool(red_extended_preview),
        anchor_pixel=(dataset["wcs"].all_world2pix(*anchor, 0)
                      if red_extended_preview else None),
        target_radius_arcsec=(float(project.target_radius_arcsec)
                              if red_extended_preview else None))
    inner_gaia_mask = counts.pop("gaia_inner_preview_mask")
    red_extended_mask = counts.pop("red_extended_preview_mask")
    source_inventory = build_scene_source_inventory(
        catalog, pixel_shape=shape, protected=protected,
        screen_stats=counts, release=release, region=region,
        brick=brick, band=band, sga_core=sga_core,
        anchor_icrs=anchor)
    foreground_review_ids = tuple(
        e["source_id"] for e in source_inventory["entries"]
        if "GAIA_ASTROMETRIC_SIGNAL_FOREGROUND_REVIEW"
           in e["evidence_flags"]
        and (e["centre_in_target_protection"] or
             e["screening_disposition"] in (
                 "PROTECTED_OVERLAP", "GAIA_FOREGROUND_MASK_PREVIEW"))
    )
    auto = auto & ~restored_auto
    # The displayed orange footprint must show EFFECTIVE proposed pixels,
    # not pixels explicitly restored by the researcher.
    inner_gaia_mask = inner_gaia_mask & ~restored_auto
    red_extended_mask = red_extended_mask & ~restored_auto
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
    if foreground_review_ids and counts["gaia_inner_preview_sources"] < (
            len(foreground_review_ids)):
        flags.append("PROTECTED_GAIA_ASTROMETRY_STARS_NOT_REMOVED")
    if red_extended_preview:
        flags.extend((
            "EXPLICIT_RED_EXTENDED_MASK_PREVIEW_ONLY",
            "COLOR_MORPHOLOGY_CLASSIFIER_PURITY_NOT_VALIDATED",
            "NO_AUTOMATIC_MEMBERSHIP_PROOF",
        ))
    if gaia_inner_preview:
        flags.extend((
            "EXPLICIT_GAIA_INTERIOR_MASK_PREVIEW_ONLY",
            "GAIA_ASTROMETRY_NOT_FIDELITY_VALIDATED",
            "GAIA_WING_AND_GALAXY_FLUX_LOSS_NOT_QUANTIFIED",
        ))
    if legacy_isophote_count:
        flags.append("LEGACY_ISOPHOTE_PROTECT_EDITS_INACTIVE")
    if sga_core is not None:
        flags.extend(("SGA_CORE_ONLY_NOT_OUTER_LSB_BOUNDARY",
                      "SGA_TRACTOR_EVIDENCE_NOT_STATISTICALLY_INDEPENDENT"))
    if caution["BAILOUT"]:
        flags.append("TRACTOR_BAILOUT_PIXELS_PRESENT")
    if caution["BRIGHT"]:
        flags.append("BRIGHT_STAR_MASKBITS_REQUIRE_WING_REVIEW")
    output = {
        "status": "CATALOGUE_CONTAMINATION_MASK_NOT_ADOPTED",
        "mask": excluded.copy(),
        "auto_mask": auto.copy(),
        "gaia_inner_preview_mask": inner_gaia_mask.copy(),
        "gaia_inner_preview_enabled": bool(gaia_inner_preview),
        "red_extended_preview_mask": red_extended_mask.copy(),
        "red_extended_preview_enabled": bool(red_extended_preview),
        "protected": protected.copy(),
        "hard": hard.copy(),
        "project_id": project.project_id,
        "draft_digest": draft.digest,
        "dataset_source": dataset["source"],
        "band": band,
        "anchor_icrs": anchor,
        "sga_prior_id": sga_prior.sga_id if sga_prior is not None else None,
        "sga_prior_sha256": (
            sga_prior.source_sha256 if sga_prior is not None else None),
        "sga_core_pixels": (
            int(np.count_nonzero(sga_core)) if sga_core is not None else 0),
        # This provisional circular exclusion/protection scale affects
        # which Tractor sources can be masked. Preserve it for independent
        # radial-boundary sensitivity QA; not a science galaxy radius.
        "target_protection_radius_arcsec": float(
            project.target_radius_arcsec),
        "catalog_sha256": catalog.source_sha256,
        "catalog_path": catalog.source_path,
        "foreground_astrometry_protected_review_ids": foreground_review_ids,
        "stats": {**counts,
                  "legacy_isophote_protect_edits_ignored":
                      legacy_isophote_count,
                  "protected_gaia_astrometry_review":
                      len(foreground_review_ids),
                  "gaia_inner_preview_pixels":
                      int(np.count_nonzero(inner_gaia_mask)),
                  "red_extended_preview_pixels":
                      int(np.count_nonzero(red_extended_mask)),
                  "catalog_rows": len(catalog.rows),
                  "excluded_pixels": int(np.count_nonzero(excluded)),
                  "native_hard_pixels": int(np.count_nonzero(hard)),
                  "protected_pixels": int(np.count_nonzero(protected)),
                  "usable_pixels": int(np.count_nonzero(usable))},
        "maskbits_caution_counts": caution,
        # Read-only raw Tractor-component evidence. It MUST NEVER be
        # interpreted as galaxy membership or modify screening masks.
        "source_inventory": source_inventory,
        "flags": tuple(flags),
    }
    for name in ("mask", "auto_mask", "gaia_inner_preview_mask",
                 "red_extended_preview_mask", "protected", "hard"):
        output[name].setflags(write=False)
    return output
