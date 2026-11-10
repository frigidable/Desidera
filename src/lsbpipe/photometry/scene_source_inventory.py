"""Read-only, source-level evidence registry for Desidera Scene & Masks.

A Tractor component is NOT a physical object classification. This module
records all native-catalog sources projected into the original FITS cutout
(including Gaia references with NO fitted Legacy band flux), and associates
each with the actual stage-SC1 screening decision. It neither changes
mask pixels nor estimates physical membership or star formation.

The only classifications here are *screening dispositions* and observable
catalogue flags, with explicit uncertainty.
"""
from __future__ import annotations

import math

import numpy as np


SOURCE_EVIDENCE_SCHEMA = "desidera.scene-source-evidence/1"
# DR11 official MASKBITS/FITBITS; masks denote processing context, not origin.
MASKBIT_NAMES = {
    1: "BRIGHT", 2: "SATUR_G", 3: "SATUR_R", 4: "SATUR_Z",
    5: "ALLMASK_G", 6: "ALLMASK_R", 7: "ALLMASK_Z",
    10: "BAILOUT", 11: "MEDIUM", 12: "GALAXY",
    13: "CLUSTER", 14: "SATUR_I", 15: "ALLMASK_I",
    16: "SUB_BLOB", 17: "RESOLVED", 18: "MCLOUDS",
}
FITBIT_NAMES = {
    0: "FORCED_POINTSOURCE", 1: "FIT_BACKGROUND",
    2: "HIT_RADIUS_LIMIT", 3: "HIT_SERSIC_LIMIT",
    4: "FROZEN", 5: "BRIGHT", 6: "MEDIUM", 7: "GAIA",
}


def _flags(value, names):
    try:
        v = int(value or 0)
    except (TypeError, ValueError):
        return ()
    return tuple(name for bit, name in names.items() if v & (1 << bit))


def _sig(value, invvar):
    try:
        if (value is None or invvar is None or
                not math.isfinite(float(value)) or
                not math.isfinite(float(invvar)) or float(invvar) <= 0):
            return None
        return abs(float(value))*math.sqrt(float(invvar))
    except (ValueError, TypeError, OverflowError):
        return None


def _separation_arcsec(ra, dec, centre):
    if centre is None or ra is None or dec is None:
        return None
    ra0, dec0 = map(float, centre)
    ra, dec = float(ra), float(dec)
    if not all(math.isfinite(v) for v in (ra0,dec0,ra,dec)):
        return None
    # Haversine remains stable across the RA=0 meridian and at small
    # dwarf-galaxy angular separations.
    da = math.radians((ra-ra0+180.)%360.-180.)
    dd = math.radians(dec-dec0)
    a = (math.sin(dd/2.)**2 +
         math.cos(math.radians(dec))*math.cos(math.radians(dec0))
         * math.sin(da/2.)**2)
    return math.degrees(2.*math.asin(math.sqrt(min(1.,max(0.,a)))))*3600.


def _point_key(x, y, objid):
    return (round(float(x), 4), round(float(y), 4),
            int(objid) if objid is not None else None)


def _described_source(row, *, x, y, disposition, release, region, brick,
                      band, centre_protected, centre_sga, angular_arcsec):
    ref_cat = str(row.get("ref_cat") or "").upper()
    ref_id = row.get("ref_id")
    is_gaia = ref_cat == "G3" and ref_id is not None and int(ref_id) > 0
    mflags = _flags(row.get("maskbits"), MASKBIT_NAMES)
    fflags = _flags(row.get("fitbits"), FITBIT_NAMES)
    flux_ivar = row.get("flux_ivar")
    fitted_flux = (flux_ivar is not None and
                   math.isfinite(float(flux_ivar)) and float(flux_ivar) > 0)
    par_sig = _sig(row.get("parallax_mas"), row.get("parallax_ivar"))
    pmra_sig = _sig(row.get("pmra_masyr"), row.get("pmra_ivar"))
    pmdec_sig = _sig(row.get("pmdec_masyr"), row.get("pmdec_ivar"))
    evidence = []
    if is_gaia:
        evidence.append("GAIA_DR3_REFERENCE_ID")
    elif ref_cat == "G3":
        evidence.append("GAIA_REF_WITHOUT_VALID_ID")
    elif (row.get("gaia_g_mag") is not None and
          0. < float(row["gaia_g_mag"]) < 35.):
        evidence.append("GAIA_MAG_ONLY_REFERENCE_NOT_VERIFIED")
    if ref_cat == "L4":
        evidence.append("SGA_REFERENCE_NOT_TARGET_MEMBERSHIP")
    if is_gaia and (
            (par_sig is not None and par_sig >= 5) or
            (pmra_sig is not None and pmra_sig >= 5) or
            (pmdec_sig is not None and pmdec_sig >= 5)):
        evidence.append("GAIA_ASTROMETRIC_SIGNAL_FOREGROUND_REVIEW")
    if row.get("gaia_duplicated_source"):
        evidence.append("GAIA_DUPLICATED_SOURCE")
    if not fitted_flux:
        evidence.append("NO_FITTED_TRACTOR_BAND_FLUX")
        if is_gaia:
            evidence.append("RETAINED_GAIA_SOURCE_POSSIBLE")
    if "FORCED_POINTSOURCE" in fflags or (
            "PSF" == row.get("type") and
            any(b in mflags for b in ("BRIGHT","GALAXY","CLUSTER"))):
        evidence.append("PSF_TYPE_NOT_STELLAR_MEMBERSHIP_PROOF")
    if "RESOLVED" in mflags:
        evidence.append("RESOLVED_DWARF_CONTEXT_NOT_MEMBERSHIP_PROOF")
    if "GALAXY" in mflags:
        evidence.append("SGA_GALAXY_CONTEXT_NOT_MEMBERSHIP_PROOF")
    if "BAILOUT" in mflags:
        evidence.append("TRACTOR_BAILOUT_REVIEW")
    if "BRIGHT" in mflags or "MEDIUM" in mflags:
        evidence.append("BRIGHT_STAR_WINGS_REVIEW")
    if any(f in mflags for f in ("SATUR_G","SATUR_R","SATUR_I",
                                 "SATUR_Z","ALLMASK_G","ALLMASK_R",
                                 "ALLMASK_I","ALLMASK_Z")):
        evidence.append("SATURATION_OR_ALLMASK_REVIEW")
    if centre_protected or disposition == "PROTECTED_OVERLAP":
        evidence.append("INSIDE_OR_TOUCHES_PROVISIONAL_TARGET_PROTECTION")
    if centre_sga:
        evidence.append("SGA_CENTRAL_REGION_NOT_MEMBERSHIP_PROOF")
    if disposition == "PROPOSED_EXCLUDE":
        evidence.append("GEOMETRIC_MASK_ONLY_NO_TRACTOR_MODEL")
    objid = row.get("objid")
    sid = f"{release}/{region}/{brick}/OBJID={objid}"
    if objid is None:
        sid += f"/pixel={float(x):.3f},{float(y):.3f}"
    model_fluxes = dict(row.get("model_fluxes_nMgy") or {})
    model_ivars = dict(row.get("model_flux_ivars") or {})
    return {
        "source_id":sid, "objid":objid,
        "x":float(x), "y":float(y),
        "ra_deg":row.get("ra_deg"), "dec_deg":row.get("dec_deg"),
        "offset_arcsec":angular_arcsec,
        "tractor_type":str(row.get("type") or ""),
        "screening_disposition":disposition,
        "physical_membership":"UNDETERMINED",
        "photometry_action":"NO_NEW_ACTION_READ_ONLY",
        "centre_in_target_protection":bool(centre_protected),
        "centre_in_sga_prior":bool(centre_sga),
        "reference_catalog":ref_cat, "reference_id":ref_id,
        "gaia_dr3_source_id":int(ref_id) if is_gaia else None,
        "gaia_duplicated_source":bool(row.get("gaia_duplicated_source")),
        "gaia_g_mag":row.get("gaia_g_mag"),
        "parallax_mas":row.get("parallax_mas"),
        "parallax_significance":par_sig,
        "pmra_masyr":row.get("pmra_masyr"),
        "pmra_significance":pmra_sig,
        "pmdec_masyr":row.get("pmdec_masyr"),
        "pmdec_significance":pmdec_sig,
        "ref_epoch":row.get("ref_epoch"),
        "band":band, "band_model_flux_nMgy":row.get("flux"),
        "band_model_flux_ivar":flux_ivar,
        "model_fluxes_nMgy":model_fluxes,
        "model_flux_ivars":model_ivars,
        "maskbits":row.get("maskbits"),
        "maskbits_names":mflags,
        "fitbits":row.get("fitbits"),
        "fitbits_names":fflags,
        "shape_r_arcsec":row.get("shape_r_arcsec"),
        "shape_e1":row.get("shape_e1"),
        "shape_e2":row.get("shape_e2"),
        "evidence_flags":tuple(evidence),
    }


def build_scene_source_inventory(catalog, *, pixel_shape, protected,
                                 screen_stats, release, region, brick, band,
                                 sga_core=None, anchor_icrs=None,
                                 max_records=12000):
    """Join native catalog identities with actual SC1 disposition.

    All sources with centres in the native cutout plus protected-border
    overlap sources are eligible, even if they lack a fitted flux.
    The selected-protected and excluded sets are DIFFERENT policies,
    not physical hypotheses. Any truncation is explicitly reported.
    """
    ny, nx = (int(v) for v in pixel_shape)
    if (protected.shape != (ny,nx) or protected.dtype.kind != "b" or
            sga_core is not None and (
                np.shape(sga_core) != (ny,nx) or
                np.asarray(sga_core).dtype.kind != "b")):
        raise ValueError("Source inventory native grid mismatch")
    points = np.asarray(catalog.pixel_positions)
    if points.shape != (len(catalog.rows), 2):
        raise ValueError("Source inventory rows and positions differ")
    if not isinstance(max_records, int) or max_records < 1:
        raise ValueError("Source inventory bounded maximum required")
    protected_index = {
        _point_key(c["x"],c["y"],c.get("objid"))
        for c in screen_stats.get("near_target_sources", ())
    }
    masked_index = {
        _point_key(c["x"],c["y"],c.get("objid"))
        for c in screen_stats.get("masked_catalog_centres", ())
    }
    candidates = []
    for (x,y),row in zip(points,catalog.rows):
        if not np.isfinite(x) or not np.isfinite(y):
            continue
        x,y=float(x),float(y)
        inside = 0 <= x < nx and 0 <= y < ny
        key = _point_key(x,y,row.get("objid"))
        if not inside and key not in protected_index:
            continue
        i,j=int(round(x)),int(round(y))
        centre_protected = (0 <= i < nx and 0 <= j < ny
                            and bool(protected[j,i]))
        centre_sga = (sga_core is not None and
                      0 <= i < nx and 0 <= j < ny
                      and bool(sga_core[j,i]))
        if key in protected_index:
            disposition = "PROTECTED_OVERLAP"
        elif key in masked_index:
            disposition = "PROPOSED_EXCLUDE"
        else:
            disposition = "NOT_SCREENED_OR_UNFITTED"
        candidates.append((0 if disposition == "PROTECTED_OVERLAP"
                           else 1 if disposition == "PROPOSED_EXCLUDE"
                           else 2, x,y,row,disposition,
                           centre_protected,centre_sga))
    # Stable identity and focus ordering, independent of FITS row order.
    candidates.sort(key=lambda item:(
        item[0], str(item[3].get("objid")), item[1], item[2]))
    n_total = len(candidates)
    full_counts = {
        name:sum(item[4] == name for item in candidates)
        for name in ("PROTECTED_OVERLAP","PROPOSED_EXCLUDE",
                     "NOT_SCREENED_OR_UNFITTED")
    }
    if n_total > max_records:
        candidates = candidates[:max_records]
    entries = tuple(_described_source(
        row,x=x,y=y,disposition=disposition,
        release=release,region=region,brick=brick,band=band,
        centre_protected=cp,centre_sga=cs,
        angular_arcsec=_separation_arcsec(
            row.get("ra_deg"),row.get("dec_deg"),anchor_icrs))
        for _,x,y,row,disposition,cp,cs in candidates)
    source_screening_truncated = (
        screen_stats.get("skipped_target_overlap",0) >
        len(screen_stats.get("near_target_sources",())) or
        screen_stats.get("masked_catalog_sources",0) >
        len(screen_stats.get("masked_catalog_centres",()))
    )
    return {
        "schema":SOURCE_EVIDENCE_SCHEMA,
        "catalog_sha256":catalog.source_sha256,
        "catalog_source_path":catalog.source_path,
        "dataset_id":catalog.dataset_id if hasattr(catalog,"dataset_id")
                     else f"{release}/{region}/{brick}/{band}",
        "total_in_native_frame":n_total,
        "stored_count":len(entries),
        "truncated":max(0,n_total-len(entries)),
        "screening_counts_all_native":full_counts,
        "screening_evidence_incomplete":bool(source_screening_truncated),
        "screening_stats_legacy":{
            "skipped_target_overlap":screen_stats.get(
                "skipped_target_overlap",0),
            "masked_catalog_sources":screen_stats.get(
                "masked_catalog_sources",0),
            "unusable_catalog_sources":screen_stats.get(
                "unusable_catalog_sources",0),
        },
        "entries":entries,
        "status":"SOURCE_EVIDENCE_REVIEW_ONLY",
        "cautions":(
            "TRACTOR_COMPONENT_NOT_PHYSICAL_OBJECT",
            "GAIA_DR3_REFERENCE_NOT_FOREGROUND_MEMBERSHIP",
            "NO_HII_OR_YOUNG_CLUSTER_CLASSIFICATION",
            "NO_NEW_MASK_OR_PHOTOMETRY_ACTION",
            "SC1_MASKED_CENTRE_INVENTORY_MAY_BE_CAPPED",
        ) + (("INVENTORY_TRUNCATED",) if n_total>len(entries) else ())
          + (("UPSTREAM_SCREENING_OUTCOMES_TRUNCATED",)
             if source_screening_truncated else ()),
    }
