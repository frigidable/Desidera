"""P16B: masked, *original-pixel* outer stellar-profile feasibility audit.

This observes azimuthal medians of selected stars on larger native DR11
coadd cutouts, with two DISJOINT conditional outer-sky hypotheses and
conservative catalog-neighbor cores masked. It is emphatically NOT a
wide-wing PSF reconstruction, sky-free halo detection, statistical
significance test or calibration for LSB science. No model is fitted or
adopted; no pixel/PSF/star-sample array is changed.
"""
from __future__ import annotations

import math
import numpy as np

SCHEMA = "desidera.psf0-original-native-wide-wing-profile/1"
PROFILE_EDGES_IN_MEASURED_FWHM = (1., 2., 3., 4., 6., 8., 12., 16., 20., 22.)
SKY_RING_A = (22., 26.)
SKY_RING_B = (26., 30.)
NEIGHBOR_MASK_RADIUS_IN_MEASURED_FWHM = 3.


def _star_outer_profile(frame, star, positions, fwhm_pix, protected_mask):
    """Read-only original-pixel annular medians with explicit limitations."""
    center_x, center_y = float(star.x_pix), float(star.y_pix)
    record = {
        "source_id":str(star.star_id),
        "role":str(star.role),
        "status":"NOT_EVALUATED",
        "source_selected_before_outer_wing_measurement":True,
        "original_image_and_inverse_variance_unchanged":True,
        "no_model_fit_or_psf_adoption":True,
    }
    half = int(math.ceil(SKY_RING_B[1]*fwhm_pix)) + 1
    ix, iy = int(np.rint(center_x)), int(np.rint(center_y))
    ny, nx = frame.shape
    x0, y0 = ix-half, iy-half
    if x0 < 0 or y0 < 0 or ix+half >= nx or iy+half >= ny:
        return dict(
            record,status="FULL_OUTER_NATIVE_CUTOUT_OUTSIDE_IMAGE",
            half_width_native_pix=half)
    rows=slice(y0,iy+half+1)
    cols=slice(x0,ix+half+1)
    pixels=np.asarray(frame.image_nMgy[rows,cols],dtype=float)
    inv=np.asarray(frame.invvar_nMgy2[rows,cols],dtype=float)
    usable=np.asarray(frame.valid[rows,cols],dtype=bool).copy()
    for key in (f"SATUR_{frame.band.upper()}",
                f"ALLMASK_{frame.band.upper()}"):
        if key in frame.mask_reasons:
            usable &= ~np.asarray(frame.mask_reasons[key][rows,cols],dtype=bool)
    usable &= np.isfinite(pixels)&np.isfinite(inv)&(inv>0)
    if protected_mask is not None:
        usable &= ~protected_mask[rows,cols]
    gy,gx=np.indices(pixels.shape,dtype=float)
    gx+=x0
    gy+=y0
    r=np.hypot(gx-center_x,gy-center_y)
    angle=np.arctan2(gy-center_y,gx-center_x)%(2*np.pi)
    sectors=np.floor(angle*(8/(2*np.pi))).astype(int)
    neighbor_radius=NEIGHBOR_MASK_RADIUS_IN_MEASURED_FWHM*fwhm_pix
    points=np.asarray(positions,dtype=float)
    dist=np.hypot(points[:,0]-center_x,points[:,1]-center_y)
    central=np.flatnonzero(dist<=2.5)
    if len(central)!=1:
        return dict(record,status="CATALOG_CENTRAL_MATCH_NOT_UNIQUE",
                    n_possible_catalog_centers=int(len(central)))
    neighbors=np.flatnonzero(
        (dist>2.5)&(dist<=half*math.sqrt(2)+neighbor_radius))
    n_intersecting=0
    for j in neighbors:
        px,py=points[j]
        if not (x0-neighbor_radius<=px<=ix+half+neighbor_radius
                and y0-neighbor_radius<=py<=iy+half+neighbor_radius):
            continue
        n_intersecting+=1
        usable &= ((gx-px)**2+(gy-py)**2)>neighbor_radius**2

    def summary(a,b):
        band=(r>=a*fwhm_pix)&(r<b*fwhm_pix)
        values=pixels[band&usable]
        present=int(values.size)
        available=int(np.count_nonzero(band))
        sectors_present=[
            int(k) for k in range(8) if np.any(band&usable&(sectors==k))]
        return {
            "radius_inner_in_measured_fwhm":float(a),
            "radius_outer_in_measured_fwhm":float(b),
            "n_original_pixels_geometrically_in_ring":available,
            "n_unmasked_native_pixels":present,
            "unmasked_coverage_fraction":float(present/available)
                if available else None,
            "n_azimuthal_octants_with_usable_pixels":len(sectors_present),
            "median_original_nMgy_per_native_pixel":(
                float(np.median(values)) if present else None),
            "median_is_not_flux_weighted_or_covariance_calibrated":True,
        }

    sky_a=summary(*SKY_RING_A)
    sky_b=summary(*SKY_RING_B)
    med_a=sky_a["median_original_nMgy_per_native_pixel"]
    med_b=sky_b["median_original_nMgy_per_native_pixel"]
    peak=float(star.peak_nMgy)
    if (med_a is None or med_b is None
            or not math.isfinite(peak) or peak<=0):
        return dict(
            record,status="INSUFFICIENT_OR_INVALID_CONDITIONAL_OUTER_SKY",
            sky_hypothesis_inner=sky_a,sky_hypothesis_outer=sky_b,
            n_catalog_neighbors_whose_cores_were_masked=n_intersecting)
    profiles=[]
    for a,b in zip(PROFILE_EDGES_IN_MEASURED_FWHM[:-1],
                   PROFILE_EDGES_IN_MEASURED_FWHM[1:]):
        entry=summary(a,b)
        raw=entry["median_original_nMgy_per_native_pixel"]
        entry.update(
            sky_a_subtracted_median_over_original_star_peak=(
                float((raw-med_a)/peak) if raw is not None else None),
            sky_b_subtracted_median_over_original_star_peak=(
                float((raw-med_b)/peak) if raw is not None else None),
            ring_is_a_direct_robust_native_pixel_statistic_not_wing_fit=True,
        )
        profiles.append(entry)
    return dict(
        record,
        status="MASKED_OUTER_PROFILE_DESCRIPTIVE_SKY_DEGENERATE",
        half_width_native_pix=half,
        original_native_cutout_shape_pix=list(pixels.shape),
        measured_core_fwhm_native_pix=float(fwhm_pix),
        n_catalog_neighbors_whose_cores_were_masked=n_intersecting,
        catalog_neighbor_exclusion_radius_native_pix=float(neighbor_radius),
        catalog_mask_only_covers_neighbor_cores_not_their_wings=True,
        original_target_protection_mask_applied=protected_mask is not None,
        conditional_sky_inner_ring=sky_a,
        conditional_sky_outer_ring=sky_b,
        sky_hypotheses_from_disjoint_original_native_annuli=True,
        sky_hypothesis_difference_in_original_rms_units_not_calibrated=True,
        sky_inner_minus_outer_nMgy_per_native_pixel=float(med_a-med_b),
        sky_difference_over_star_peak=float((med_a-med_b)/peak),
        original_star_peak_nMgy_per_native_pixel=float(peak),
        radial_profile=profiles,
        individual_annular_medians_are_correlated_and_sky_degenerate=True,
        stellar_wing_detection_not_claimed=True,
        bright_stars_independent_of_existing_psf_training_not_identified=True,
        no_psf_kernel_change=True,
        science_adopted=False,
    )


def audit_native_broad_radius_star_profiles(frame, stars, catalog_matcher,
                                            core_geometry, *,
                                            protected_mask=None):
    """Probe broad-radius feasibility on unchanged selected-star coadd pixels.

    Does not choose or reject stars from residuals. In particular, fixed
    holdout stars remain held out in the PSF fit even though read-only
    wing data are shown to diagnose the current model.
    """
    out={
        "schema":SCHEMA,
        "status":"NOT_RUN",
        "science_adopted":False,
        "all_selected_stars_fixed_before_broad_radius_probe":True,
        "not_an_independent_bright_star_wing_calibration":True,
        "not_an_outer_sky_or_wide_angle_psf_solution":True,
        "no_original_pixels_or_psf_kernel_modified":True,
        "radial_edges_in_measured_core_fwhm":list(
            PROFILE_EDGES_IN_MEASURED_FWHM),
        "two_disjoint_outer_sky_hypothesis_rings_in_measured_fwhm":[
            list(SKY_RING_A),list(SKY_RING_B)],
        "neighbor_mask_radius_in_measured_fwhm":
            NEIGHBOR_MASK_RADIUS_IN_MEASURED_FWHM,
        "no_automatic_star_selection_or_psf_adoption":True,
        "per_star":[],
    }
    if catalog_matcher is None:
        return dict(out,status="NO_VERIFIED_NATIVE_CATALOG_FOR_NEIGHBOR_MASK")
    positions=np.asarray(catalog_matcher.pixel_positions,dtype=float)
    if (positions.ndim!=2 or positions.shape[1]!=2
            or not len(positions) or not np.isfinite(positions).all()):
        return dict(out,status="INVALID_NATIVE_CATALOG_POSITIONS")
    if protected_mask is not None:
        protected_mask=np.asarray(protected_mask,dtype=bool)
        if protected_mask.shape!=frame.shape:
            return dict(out,status="INVALID_TARGET_PROTECTION_MASK_SHAPE")
    if not isinstance(core_geometry,dict):
        return dict(out,status="NO_CORE_GEOMETRY")
    width=core_geometry.get(
        "measured_core_fwhm_axis_geometric_mean_native_pix")
    if (core_geometry.get("status")!="FINITE_NATIVE_CORE_GEOMETRY_MEASURED"
            or not isinstance(width,(float,int))
            or not math.isfinite(width) or width<=0):
        return dict(out,status="NO_MEASURED_NATIVE_CORE_FWHM")
    rows=[
        _star_outer_profile(frame, star, positions, float(width),
                            protected_mask)
        for star in stars
    ]
    n_ok=sum(x["status"]=="MASKED_OUTER_PROFILE_DESCRIPTIVE_SKY_DEGENERATE"
             for x in rows)
    return dict(
        out,
        status="ORIGINAL_NATIVE_OUTER_PROFILE_FEASIBILITY_ONLY",
        n_selected_stars=len(rows),
        n_scenes_with_two_outer_sky_hypotheses=n_ok,
        n_scenes_not_measurable=len(rows)-n_ok,
        original_target_protection_mask_applied=protected_mask is not None,
        n_fixed_holdouts_in_original_selected_inventory=sum(
            x["role"]=="holdout" for x in rows),
        n_fixed_holdouts_with_outer_profile=sum(
            x["role"]=="holdout"
            and x["status"]=="MASKED_OUTER_PROFILE_DESCRIPTIVE_SKY_DEGENERATE"
            for x in rows),
        per_star=rows,
    )


# P16C: these fixed catalog-only rules must never be optimized on the P16B
# stellar profile residuals. Prospective means ranked before P16C pixel reads;
# its stars are disjoint from the *existing* ePSF training and holdouts, but
# they are not independent of the DR11 coadd, sky or instrument.
PROSPECTIVE_SCHEMA = "desidera.psf0-disjoint-bright-tractor-gaia-native-profile/1"
PROSPECTIVE_LIMIT = 16
PROSPECTIVE_MIN_DIAGONAL_TRACTOR_SNR = 20.


def audit_disjoint_bright_gaia_wing_probes(
    frame, stars, catalog_matcher, core_geometry, *, protected_mask=None,
    max_catalog_probes=None
):
    """Catalog-preregistered bright-star wing *feasibility* sample.

    Selection uses only the original native Tractor row, coordinates, the
    pre-existing target exclusion geometry and the *frozen* core FWHM. Never
    select/refill on observed flux, radial profile, measured background, mask
    coverage or perceived scientific merit. A selected star with bad original
    pixels stays selected but is returned as failed, never replaced.
    """
    if max_catalog_probes is None:
        max_catalog_probes = PROSPECTIVE_LIMIT
    if type(max_catalog_probes) is not int or not 1 <= max_catalog_probes <= 48:
        raise ValueError("Catalog-only fixed probe quota must be 1..48")
    from collections import Counter
    from types import SimpleNamespace
    from scipy.spatial import cKDTree
    from lsbpipe.photometry.psf_tractor_gaia import (
        MASK_REJECT_BITS, FIT_REJECT_BITS, SATURATION_BITS,
    )

    report = {
        "schema": PROSPECTIVE_SCHEMA,
        "status": "NOT_RUN",
        "science_adopted": False,
        "catalog_selection_uses_image_or_wing_residuals": False,
        "catalog_rank_uses_native_tractor_band_flux_descending": True,
        "predeclared_max_catalog_star_probes": max_catalog_probes,
        "predeclared_min_diagonal_tractor_snr_not_covariance_calibrated":
            PROSPECTIVE_MIN_DIAGONAL_TRACTOR_SNR,
        "existing_kernel_training_and_fixed_holdout_stars_excluded": True,
        "prospective_star_identities_disjoint_not_independent_coadd": True,
        "all_selected_sources_audited_without_replacement_or_resplitting": True,
        "no_correction_for_uncataloged_sources_or_neighbor_wings": True,
        "no_original_pixels_kernel_training_or_science_gates_changed": True,
        "predeclared_catalog_rejection_counts": {},
        "per_star": [],
    }
    if catalog_matcher is None or not hasattr(catalog_matcher, "rows"):
        return dict(report, status="NO_FULL_NATIVE_GAIA_TRACTOR_CATALOG_METADATA")
    positions = np.asarray(catalog_matcher.pixel_positions, dtype=float)
    rows = catalog_matcher.rows
    if (positions.shape != (len(rows), 2) or
            not np.isfinite(positions).all() or not len(rows)):
        return dict(report, status="INVALID_NATIVE_CATALOG_FOR_NEW_COHORT")
    if not isinstance(core_geometry, dict):
        return dict(report, status="NO_MEASURED_NATIVE_CORE_FWHM")
    width = core_geometry.get(
        "measured_core_fwhm_axis_geometric_mean_native_pix")
    if (core_geometry.get("status") != "FINITE_NATIVE_CORE_GEOMETRY_MEASURED"
            or not isinstance(width, (int, float))
            or not math.isfinite(width) or width <= 0):
        return dict(report, status="NO_MEASURED_NATIVE_CORE_FWHM")
    if protected_mask is not None:
        protected_mask = np.asarray(protected_mask, dtype=bool)
        if protected_mask.shape != frame.shape:
            return dict(report, status="INVALID_TARGET_PROTECTION_MASK_SHAPE")

    if frame.band not in SATURATION_BITS:
        return dict(report, status="UNSUPPORTED_NATIVE_CATALOG_BAND")
    excluded_ids = {
        int(star.gaia_source_id) for star in stars
        if getattr(star, "gaia_source_id", None) is not None
    }
    excluded_positions = np.asarray(
        [(float(s.x_pix), float(s.y_pix)) for s in stars], dtype=float
    ).reshape(-1, 2)
    tree = cKDTree(positions)
    heldout_tree = (cKDTree(excluded_positions) if len(excluded_positions)
                    else None)
    cutoff_half = int(math.ceil(SKY_RING_B[1] * width)) + 1
    minimum_other_catalog_separation = (
        NEIGHBOR_MASK_RADIUS_IN_MEASURED_FWHM * width
    )
    # REF_ID is an identity, not a guarantee that the catalog itself has
    # no duplicated REF_ID across separate Tractor rows. Refuse every
    # ambiguous Gaia identity before catalog-only flux ranking.
    gaia_identity_counts = Counter(
        (row.get("ref_id") for row in rows
         if isinstance(row, dict) and row.get("ref_cat") == "G3"
         and isinstance(row.get("ref_id"), int) and row["ref_id"] > 0))
    candidates = []
    rejected = Counter()
    for index, (row, xy) in enumerate(zip(rows, positions)):
        if not isinstance(row, dict):
            rejected["MALFORMED_ORIGINAL_TRACTOR_ROW"] += 1
            continue
        flux, ivar = row.get("flux"), row.get("flux_ivar")
        if (row.get("ref_cat") != "G3" or
                not isinstance(row.get("ref_id"), int) or
                row["ref_id"] <= 0 or row.get("type") != "PSF"):
            rejected["NOT_VERIFIED_GAIA_PSF"] += 1
            continue
        if gaia_identity_counts[row["ref_id"]] != 1:
            rejected["DUPLICATE_GAIA_ID_IN_NATIVE_TRACTOR_CATALOG"] += 1
            continue
        if (not isinstance(flux, (float, int)) or
                not isinstance(ivar, (float, int)) or
                not math.isfinite(flux) or not math.isfinite(ivar) or
                flux <= 0 or ivar <= 0 or
                flux * math.sqrt(ivar) <
                PROSPECTIVE_MIN_DIAGONAL_TRACTOR_SNR):
            rejected["INSUFFICIENT_POSITIVE_TRACTOR_BAND_FLUX"] += 1
            continue
        try:
            fitbits = int(row["fitbits"])
            maskbits = int(row["maskbits"])
        except (KeyError, TypeError, ValueError, OverflowError):
            rejected["MISSING_CATALOG_QUALITY_FLAGS"] += 1
            continue
        if (fitbits & FIT_REJECT_BITS or
                maskbits & MASK_REJECT_BITS or
                any(maskbits & (1 << bit)
                    for bit in SATURATION_BITS[frame.band]) or
                row.get("gaia_duplicated_source") is not False):
            rejected["UNSAFE_CATALOG_QUALITY_FLAGS"] += 1
            continue
        if row["ref_id"] in excluded_ids or (
                heldout_tree is not None and
                heldout_tree.query(xy)[0] <= 2.5):
            rejected["EXISTING_TRAIN_OR_HOLDOUT_STAR"] += 1
            continue
        x, y = map(float, xy)
        ix, iy = int(np.rint(x)), int(np.rint(y))
        if (ix-cutoff_half < 0 or iy-cutoff_half < 0 or
                ix+cutoff_half >= frame.shape[1] or
                iy+cutoff_half >= frame.shape[0]):
            rejected["NO_COMPLETE_30FWHM_NATIVE_STAMP"] += 1
            continue
        if protected_mask is not None and protected_mask[iy, ix]:
            rejected["TARGET_PROTECTED_STELLAR_CENTER"] += 1
            continue
        overlaps = tree.query_ball_point(xy, r=2.5)
        if len(overlaps) != 1:
            rejected["AMBIGUOUS_CATALOG_STELLAR_CENTER"] += 1
            continue
        if len(rows)>1:
            sep = float(tree.query(xy, k=2)[0][1])
            if sep <= minimum_other_catalog_separation:
                rejected["NEIGHBOR_CATALOG_CORE_OVERLAP"] += 1
                continue
        try:
            objid = int(row["objid"])
        except (KeyError, TypeError, ValueError, OverflowError):
            rejected["MISSING_TRACTOR_OBJECT_ID"] += 1
            continue
        candidates.append((float(flux), int(row["ref_id"]),
                           objid, int(index)))

    # The quota is frozen and filled BEFORE measuring the original pixels.
    # No backfilling based on the scientific result of a selected star.
    candidates.sort(key=lambda entry: (-entry[0], entry[1], entry[2]))
    fixed = candidates[:max_catalog_probes]
    original_native = np.asarray(frame.image_nMgy, dtype=float)
    ivar = np.asarray(frame.invvar_nMgy2, dtype=float)
    valid = np.asarray(frame.valid, dtype=bool)
    source_mask = np.array(valid & np.isfinite(original_native)
                           & np.isfinite(ivar) & (ivar > 0), copy=True)
    for key in (f"SATUR_{frame.band.upper()}", f"ALLMASK_{frame.band.upper()}"):
        if key in frame.mask_reasons:
            source_mask &= ~np.asarray(frame.mask_reasons[key], dtype=bool)
    if protected_mask is not None:
        source_mask &= ~protected_mask
    profiles = []
    for rank, (flux, gaia_id, objid, idx) in enumerate(fixed, start=1):
        x, y = map(float, positions[idx])
        ix, iy = int(np.rint(x)), int(np.rint(y))
        core_pixels = original_native[iy-1:iy+2, ix-1:ix+2]
        core_valid = source_mask[iy-1:iy+2, ix-1:ix+2]
        identity = {
            "source_id": f"wing-cat-{rank:04d}",
            "cohort": "PROSPECTIVE_GAIA_TRACTOR_BRIGHT_NON_EPSF",
            "role": "catalog_wing_probe",
            "gaia_source_id": gaia_id,
            "tractor_objid": objid,
            "native_tractor_flux_nMgy_catalog_not_direct_observed_peak": flux,
            "catalog_flux_rank_in_fixed_new_cohort": rank,
            "native_xy_pix": [x, y],
            "catalog_selected_before_original_pixel_profile_read": True,
            "not_frozen_psf_training_or_holdout": True,
        }
        if not np.any(core_valid):
            profiles.append(dict(identity,
                status="NO_UNMASKED_ORIGINAL_NATIVE_CORE_PEAK",
                science_adopted=False))
            continue
        measured_peak = float(np.max(core_pixels[core_valid]))
        if not math.isfinite(measured_peak) or measured_peak <= 0:
            profiles.append(dict(identity,
                status="NONPOSITIVE_ORIGINAL_NATIVE_CORE_PEAK",
                science_adopted=False))
            continue
        star = SimpleNamespace(
            star_id=identity["source_id"],
            x_pix=x, y_pix=y, role="catalog_wing_probe",
            peak_nMgy=measured_peak)
        result = _star_outer_profile(
            frame, star, positions, float(width), protected_mask)
        result.update(identity)
        result["native_peak_from_postselection_3x3_original_pixels"] = True
        # P16B has no independent bright-star sample; P16C does.
        result["bright_stars_independent_of_existing_psf_training_not_identified"] = False
        result["science_adopted"] = False
        profiles.append(result)
    return dict(
        report,
        status="DISJOINT_BRIGHT_GAIA_NATIVE_PROFILE_FEASIBILITY_ONLY",
        native_catalog_sha256=catalog_matcher.source_sha256,
        n_native_catalog_entries_examined=len(rows),
        n_catalog_only_geometry_eligible_before_brightness_cap=len(candidates),
        n_catalog_probes_fixed_before_original_pixel_reads=len(fixed),
        n_profiles_with_both_conditional_outer_skies=sum(
            row["status"] ==
            "MASKED_OUTER_PROFILE_DESCRIPTIVE_SKY_DEGENERATE"
            for row in profiles),
        n_probe_profiles_not_measurable=sum(
            row["status"] !=
            "MASKED_OUTER_PROFILE_DESCRIPTIVE_SKY_DEGENERATE"
            for row in profiles),
        predeclared_catalog_rejection_counts=dict(sorted(rejected.items())),
        excluded_existing_gaia_star_ids_count=len(excluded_ids),
        full_native_footprint_width_pix=2*cutoff_half+1,
        minimum_other_catalog_separation_native_pix=float(
            minimum_other_catalog_separation),
        catalog_qualified_source_inventory_count=len(candidates),
        no_profile_based_star_replacement=True,
        per_star=profiles,
    )


# P16D: read-only stress test on the EXACT preselected P16C identities.
# Do not use outcomes to re-rank, replace, retrain or adopt any star/PSF.
P16D_SCHEMA = "desidera.psf0-fixed-bright-gaia-neighbor-sky-stress/1"
P16D_MASK_RADII_FWHM = (3., 6.)
P16D_PROFILE_RINGS_FWHM = ((4., 6.), (6., 8.), (8., 12.),
                           (12., 16.), (16., 20.), (20., 22.))
P16D_MIN_SKY_PIXELS = 100
P16D_MIN_OCTANT_PIXELS = 25


def _p16d_conditional_sky_plane(pixels, dx, dy, fit_mask):
    """Clipped planar background *hypothesis*, never a calibrated sky."""
    if np.count_nonzero(fit_mask) < P16D_MIN_SKY_PIXELS:
        return None
    yy = pixels[fit_mask]
    design = np.column_stack((np.ones(yy.size), dx[fit_mask],
                              dy[fit_mask]))
    keep = np.ones(yy.size, dtype=bool)
    if np.linalg.matrix_rank(design) != 3:
        return None
    for _ in range(3):
        if np.count_nonzero(keep) < P16D_MIN_SKY_PIXELS:
            return None
        coeff = np.linalg.lstsq(design[keep], yy[keep], rcond=None)[0]
        residual = yy - design @ coeff
        center = np.median(residual[keep])
        scatter = 1.4826 * np.median(np.abs(residual[keep] - center))
        if not np.isfinite(scatter) or scatter <= 0:
            break
        next_keep = np.abs(residual - center) <= 3. * scatter
        if (np.count_nonzero(next_keep) < P16D_MIN_SKY_PIXELS
                or np.array_equal(next_keep, keep)):
            break
        keep = next_keep
    coeff = np.linalg.lstsq(design[keep], yy[keep], rcond=None)[0]
    if not np.all(np.isfinite(coeff)):
        return None
    return tuple(float(v) for v in coeff)


def audit_fixed_bright_gaia_neighbor_sky_stress(
    frame, catalog_matcher, p16c_qa, core_geometry, *, protected_mask=None
):
    """P16D: fixed-cohort mask-radius, octant and sky-plane sensitivity.

    Evaluates catalog neighbors at 3 and 6 measured FWHM; compares two
    disjoint sky-ring constants with a fitted outer (22--30 FWHM) plane.
    These are descriptive, highly correlated coadd-pixel statistics.
    They cannot certify stellar wings, a sky model or a photometric PSF.
    """
    report = {
        "schema": P16D_SCHEMA,
        "status": "NOT_RUN",
        "science_adopted": False,
        "fixed_p16c_cohort_never_reranked_or_replaced": True,
        "pixels_invvar_and_psf_kernel_unchanged": True,
        "mask_radii_in_measured_fwhm": list(P16D_MASK_RADII_FWHM),
        "radial_rings_in_measured_fwhm": [
            list(v) for v in P16D_PROFILE_RINGS_FWHM],
        "sky_rings_in_measured_fwhm": [
            list(SKY_RING_A), list(SKY_RING_B)],
        "plane_fit_uses_pooled_22_to_30_fwhm_pixels": True,
        "plane_is_conditional_and_can_absorb_stellar_wings": True,
        "uncataloged_sources_and_neighbor_wings_unmodeled": True,
        "no_covariance_calibration_or_statistical_detection_test": True,
        "no_core_wing_stitching_or_science_adoption": True,
        "per_star": [],
    }
    if not isinstance(p16c_qa, dict) or (
            p16c_qa.get("schema") != PROSPECTIVE_SCHEMA or
            p16c_qa.get("status") !=
            "DISJOINT_BRIGHT_GAIA_NATIVE_PROFILE_FEASIBILITY_ONLY"):
        return dict(report, status="NO_FIXED_P16C_COHORT")
    if catalog_matcher is None or not hasattr(catalog_matcher, "pixel_positions"):
        return dict(report, status="NO_VERIFIED_NATIVE_CATALOG")
    if (p16c_qa.get("native_catalog_sha256") !=
            getattr(catalog_matcher, "source_sha256", None)):
        return dict(report, status="P16C_CATALOG_LINEAGE_MISMATCH")
    positions = np.asarray(catalog_matcher.pixel_positions, dtype=float)
    if (positions.ndim != 2 or positions.shape[1] != 2 or
            len(positions) == 0 or not np.isfinite(positions).all()):
        return dict(report, status="INVALID_NATIVE_CATALOG_POSITIONS")
    width = (core_geometry.get(
        "measured_core_fwhm_axis_geometric_mean_native_pix")
        if isinstance(core_geometry, dict) and core_geometry.get("status") ==
        "FINITE_NATIVE_CORE_GEOMETRY_MEASURED" else None)
    if (not isinstance(width, (float, int)) or
            not math.isfinite(width) or width <= 0):
        return dict(report, status="NO_MEASURED_NATIVE_CORE_FWHM")
    if protected_mask is not None:
        protected_mask = np.asarray(protected_mask, dtype=bool)
        if protected_mask.shape != frame.shape:
            return dict(report, status="INVALID_TARGET_PROTECTION_MASK")
    original = np.asarray(frame.image_nMgy, dtype=float)
    inverse = np.asarray(frame.invvar_nMgy2, dtype=float)
    half = int(math.ceil(SKY_RING_B[1] * width)) + 1
    ny, nx = frame.shape
    output = []
    for old in p16c_qa.get("per_star", []):
        item = {
            "source_id": old.get("source_id"),
            "gaia_source_id": old.get("gaia_source_id"),
            "catalog_flux_rank_in_fixed_new_cohort": old.get(
                "catalog_flux_rank_in_fixed_new_cohort"),
            "p16c_status": old.get("status"),
            "status": "NOT_EVALUATED",
            "mask_variants": [],
            "science_adopted": False,
        }
        if old.get("status") != (
                "MASKED_OUTER_PROFILE_DESCRIPTIVE_SKY_DEGENERATE"):
            item["status"] = "P16C_PROFILE_NOT_AVAILABLE"
            output.append(item)
            continue
        x, y = map(float, old["native_xy_pix"])
        ix, iy = int(np.rint(x)), int(np.rint(y))
        x0, y0 = ix - half, iy - half
        if (x0 < 0 or y0 < 0 or ix + half >= nx or iy + half >= ny):
            item["status"] = "FULL_ORIGINAL_STAMP_OUTSIDE_IMAGE"
            output.append(item)
            continue
        rows, cols = slice(y0, iy + half + 1), slice(x0, ix + half + 1)
        pixels = original[rows, cols]
        ivar = inverse[rows, cols]
        usable_base = (np.asarray(frame.valid[rows, cols], dtype=bool).copy()
                       & np.isfinite(pixels) & np.isfinite(ivar) & (ivar > 0))
        for key in (f"SATUR_{frame.band.upper()}",
                    f"ALLMASK_{frame.band.upper()}"):
            if key in frame.mask_reasons:
                usable_base &= ~np.asarray(
                    frame.mask_reasons[key][rows, cols], dtype=bool)
        if protected_mask is not None:
            usable_base &= ~protected_mask[rows, cols]
        gy, gx = np.indices(pixels.shape, dtype=float)
        gx += x0
        gy += y0
        dx, dy = (gx - x) / width, (gy - y) / width
        native_radius = np.hypot(gx - x, gy - y)
        octants = np.minimum(
            np.floor((np.arctan2(dy, dx) % (2 * np.pi)) *
                     (8. / (2 * np.pi))).astype(int), 7)
        peak = float(old["original_star_peak_nMgy_per_native_pixel"])
        if not math.isfinite(peak) or peak <= 0:
            item["status"] = "INVALID_P16C_NATIVE_PEAK"
            output.append(item)
            continue
        dist = np.hypot(positions[:, 0] - x, positions[:, 1] - y)
        if np.count_nonzero(dist <= 2.5) != 1:
            item["status"] = "CATALOG_CENTRAL_MATCH_CHANGED"
            output.append(item)
            continue
        for mask_fwhm in P16D_MASK_RADII_FWHM:
            mask_radius = mask_fwhm * width
            usable = usable_base.copy()
            eligible = ((dist > 2.5) &
                        (dist <= half * math.sqrt(2) + mask_radius) &
                        (positions[:, 0] >= x0 - mask_radius) &
                        (positions[:, 0] <= ix + half + mask_radius) &
                        (positions[:, 1] >= y0 - mask_radius) &
                        (positions[:, 1] <= iy + half + mask_radius))
            neighbor_indices = np.flatnonzero(eligible)
            for j in neighbor_indices:
                px, py = positions[j]
                usable &= ((gx - px) ** 2 + (gy - py) ** 2 >
                           mask_radius ** 2)
            def median_in(selection):
                if not np.any(selection):
                    return None
                return float(np.median(pixels[selection]))
            sky_a_mask = (usable & (native_radius >= 22 * width)
                          & (native_radius < 26 * width))
            sky_b_mask = (usable & (native_radius >= 26 * width)
                          & (native_radius < 30 * width))
            sky_a = median_in(sky_a_mask)
            sky_b = median_in(sky_b_mask)
            sky_pixels = sky_a_mask | sky_b_mask
            plane = _p16d_conditional_sky_plane(pixels, dx, dy, sky_pixels)
            sky_model = (plane[0] + plane[1] * dx + plane[2] * dy
                         if plane is not None else None)
            record = {
                "mask_radius_in_measured_fwhm": mask_fwhm,
                "n_catalog_neighbors_whose_cores_masked": int(
                    len(neighbor_indices)),
                "sky_ring_a_n_unmasked": int(np.count_nonzero(sky_a_mask)),
                "sky_ring_b_n_unmasked": int(np.count_nonzero(sky_b_mask)),
                "sky_ring_a_median_nMgy_per_pixel": sky_a,
                "sky_ring_b_median_nMgy_per_pixel": sky_b,
                "sky_a_minus_b_over_native_peak": (
                    (sky_a - sky_b) / peak if sky_a is not None and
                    sky_b is not None else None),
                "conditional_outer_sky_plane_coefficients_nMgy_per_pix": (
                    list(plane) if plane is not None else None),
                "conditional_plane_center_minus_sky_b_over_peak": (
                    (plane[0] - sky_b) / peak if plane is not None
                    and sky_b is not None else None),
                "radial_profile": [],
            }
            for inner, outer in P16D_PROFILE_RINGS_FWHM:
                area = ((native_radius >= inner * width)
                        & (native_radius < outer * width))
                ring = usable & area
                total = int(np.count_nonzero(area))
                n = int(np.count_nonzero(ring))
                octant_values = []
                for sector in range(8):
                    subset = ring & (octants == sector)
                    octant_values.append(
                        float(np.median(pixels[subset] - sky_b) / peak)
                        if sky_b is not None and
                        np.count_nonzero(subset) >= P16D_MIN_OCTANT_PIXELS
                        else None)
                populated = [v for v in octant_values if v is not None]
                raw = median_in(ring)
                plane_residual = (
                    float(np.median((pixels - sky_model)[ring]) / peak)
                    if sky_model is not None and n else None)
                record["radial_profile"].append({
                    "radius_inner_in_measured_fwhm": inner,
                    "radius_outer_in_measured_fwhm": outer,
                    "n_unmasked_native_pixels": n,
                    "unmasked_fraction": n / total if total else None,
                    "sky_a_subtracted_median_over_peak": (
                        (raw - sky_a) / peak if raw is not None
                        and sky_a is not None else None),
                    "sky_b_subtracted_median_over_peak": (
                        (raw - sky_b) / peak if raw is not None
                        and sky_b is not None else None),
                    "conditional_plane_subtracted_median_over_peak":
                        plane_residual,
                    "sky_b_subtracted_octant_medians_over_peak":
                        octant_values,
                    "n_octants_with_at_least_25_pixels": len(populated),
                    "max_minus_min_octant_median_over_peak": (
                        max(populated) - min(populated)
                        if len(populated) >= 4 else None),
                })
            item["mask_variants"].append(record)
        baseline, stress = item["mask_variants"]
        base_sky_a = old["conditional_sky_inner_ring"][
            "median_original_nMgy_per_native_pixel"]
        base_sky_b = old["conditional_sky_outer_ring"][
            "median_original_nMgy_per_native_pixel"]
        # Reproduce P16C at 3 FWHM before interpreting mask stress.
        agrees = (
            baseline["sky_ring_a_median_nMgy_per_pixel"] is not None
            and baseline["sky_ring_b_median_nMgy_per_pixel"] is not None
            and np.isclose(baseline["sky_ring_a_median_nMgy_per_pixel"],
                           base_sky_a, rtol=0, atol=1.e-10)
            and np.isclose(baseline["sky_ring_b_median_nMgy_per_pixel"],
                           base_sky_b, rtol=0, atol=1.e-10))
        p16c_rings = {
            (row["radius_inner_in_measured_fwhm"],
             row["radius_outer_in_measured_fwhm"]): row
            for row in old["radial_profile"]
        }
        for measured in baseline["radial_profile"]:
            key = (measured["radius_inner_in_measured_fwhm"],
                   measured["radius_outer_in_measured_fwhm"])
            orig = p16c_rings.get(key)
            raw = (orig.get(
                "sky_b_subtracted_median_over_original_star_peak")
                if orig is not None else None)
            new = measured["sky_b_subtracted_median_over_peak"]
            agrees = bool(
                agrees and raw is not None and new is not None
                and np.isclose(raw, new, rtol=0, atol=1.e-10))
        item["reproduces_original_p16c_three_fwhm_sky_and_profile"] = bool(agrees)
        for original_ring, stress_ring in zip(
                baseline["radial_profile"], stress["radial_profile"]):
            a = original_ring["sky_b_subtracted_median_over_peak"]
            b = stress_ring["sky_b_subtracted_median_over_peak"]
            stress_ring["six_minus_three_fwhm_mask_residual_over_peak"] = (
                b - a if a is not None and b is not None else None)
        item["status"] = (
            "BASELINE_P16C_REPRODUCTION_MISMATCH" if not agrees
            else "DESCRIPTIVE_SKY_NEIGHBOR_STRESS_ONLY"
            if all(v["sky_ring_a_n_unmasked"] >= P16D_MIN_SKY_PIXELS
                   and v["sky_ring_b_n_unmasked"] >= P16D_MIN_SKY_PIXELS
                   and v["sky_ring_a_median_nMgy_per_pixel"] is not None
                   and v["sky_ring_b_median_nMgy_per_pixel"] is not None
                   and v["conditional_outer_sky_plane_coefficients_nMgy_per_pix"]
                   is not None for v in item["mask_variants"])
            else "INSUFFICIENT_SKY_OR_PLANE_COVERAGE")
        output.append(item)
    return dict(
        report,
        status="FIXED_COHORT_SKY_NEIGHBOR_STRESS_QA_ONLY",
        n_fixed_p16c_star_identities=len(output),
        n_baseline_p16c_reproduced=sum(
            row.get("reproduces_original_p16c_three_fwhm_sky_and_profile")
            is True for row in output),
        n_stress_profiles_with_both_sky_and_plane=sum(
            row["status"] == "DESCRIPTIVE_SKY_NEIGHBOR_STRESS_ONLY"
            for row in output),
        per_star=output,
    )



# P17: fixed BEFORE inspecting any additional Leo P wing profiles.
# This is a single catalog-constrained sensitivity experiment, not an
# optimizer for isolation radius, sample size or science-ready PSF.
P17_SCHEMA = "desidera.psf0-fixed-48-gaia-catalog-neighbor-masks/1"
P17_COHORT_LIMIT = 48
P17_FROZEN_P16C_PREFIX = 16
P17_CONTROL_EVERY_NEW_RANK = 4
P17_MASK_CHOICES = ("fixed_3_fwhm", "fixed_4p5_fwhm",
                    "flux_tiered_catalog", "fixed_6_fwhm")
P17_ADAPTIVE_FLUX_THRESHOLDS = (0.01, 0.1)


def _p17_neighbor_mask(pixels_shape, x0, y0, positions, indices, radii):
    """Mask *local subrectangles*: avoid stamp-sized masks per neighbor."""
    result = np.ones(pixels_shape, dtype=bool)
    height, width = pixels_shape
    for j in indices:
        x, y = map(float, positions[j])
        radius = float(radii[j])
        left = max(0, int(math.floor(x - radius)) - x0)
        right = min(width, int(math.ceil(x + radius)) - x0 + 1)
        low = max(0, int(math.floor(y - radius)) - y0)
        high = min(height, int(math.ceil(y + radius)) - y0 + 1)
        if left >= right or low >= high:
            continue
        gy, gx = np.ogrid[y0+low:y0+high, x0+left:x0+right]
        result[low:high, left:right] &= ((gx - x)**2 + (gy - y)**2 >
                                            radius**2)
    return result


def _p17_adaptive_neighbor_radii(rows, center_flux, fwhm):
    """Predeclared Tractor-flux tiers, not a physically validated wing model.

    Missing/nonpositive neighboring flux, or a non-PSF neighbor type,
    conservatively receives 6 FWHM. Valid PSF neighbor flux ratio:
    below 0.01 -> 3, [0.01, 0.1) -> 4.5, >=0.1 -> 6 FWHM.
    """
    result = []
    for row in rows:
        flux = row.get("flux") if isinstance(row, dict) else None
        source_type = row.get("type") if isinstance(row, dict) else None
        if (source_type != "PSF" or not isinstance(flux, (int, float))
                or not math.isfinite(flux) or flux <= 0):
            tier = 6.
        elif flux / center_flux >= P17_ADAPTIVE_FLUX_THRESHOLDS[1]:
            tier = 6.
        elif flux / center_flux >= P17_ADAPTIVE_FLUX_THRESHOLDS[0]:
            tier = 4.5
        else:
            tier = 3.
        result.append(float(tier * fwhm))
    return np.asarray(result, dtype=float)


def audit_extended_bright_gaia_mask_experiment(
    frame, stars, catalog_matcher, core_geometry, p16c_qa, *,
    protected_mask=None
):
    """P17: 16 legacy + 24 new exploratory + 8 new identity-held-out probes.

    The control cohort is fixed by the Tractor-flux rank BEFORE native-pixel
    profiling, but still shares the same coadd and correlated sky.
    All mask strategies are preregistered and computed for every fixed star.
    No masking strategy is selected based on these residuals.
    """
    report = {
        "schema": P17_SCHEMA,
        "status": "NOT_RUN",
        "science_adopted": False,
        "target_catalog_rank_limit": P17_COHORT_LIMIT,
        "legacy_first_16_p16c_identities_preserved": True,
        "prospective_control_every_fourth_rank_from_17_to_48": True,
        "selection_uses_only_catalog_and_frozen_core_geometry": True,
        "selection_and_masking_do_not_depend_on_observed_wing_profile": True,
        "native_catalog_flux_used_for_neighbor_mask_hypotheses_not_psf_fit":
            True,
        "adaptive_flux_ratio_thresholds": list(P17_ADAPTIVE_FLUX_THRESHOLDS),
        "adaptive_non_psf_or_missing_flux_receives_6_fwhm": True,
        "predeclared_mask_variants": list(P17_MASK_CHOICES),
        "mask_tier_3_fwhm_for_psf_neighbor_below_flux_ratio": 0.01,
        "mask_tier_4p5_fwhm_for_psf_neighbor_flux_ratio": [0.01, 0.1],
        "mask_tier_6_fwhm_for_psf_neighbor_flux_ratio_at_least": 0.1,
        "sky_plane_is_not_an_independent_calibration": True,
        "octants_and_coadd_pixels_not_statistically_independent": True,
        "new_control_stars_not_independent_of_same_coadd": True,
        "no_wing_detection_or_model_fit_or_psf_adoption": True,
        "per_star": [],
    }
    if (not isinstance(p16c_qa, dict) or
            p16c_qa.get("schema") != PROSPECTIVE_SCHEMA or
            p16c_qa.get("status") !=
            "DISJOINT_BRIGHT_GAIA_NATIVE_PROFILE_FEASIBILITY_ONLY" or
            len(p16c_qa.get("per_star", [])) != P17_FROZEN_P16C_PREFIX):
        return dict(report, status="NO_COMPLETE_FROZEN_16_STAR_P16C_COHORT")
    if (catalog_matcher is None or
            not hasattr(catalog_matcher, "rows") or
            not hasattr(catalog_matcher, "pixel_positions") or
            getattr(catalog_matcher, "source_sha256", None) !=
            p16c_qa.get("native_catalog_sha256")):
        return dict(report, status="CATALOG_LINEAGE_OR_METADATA_MISMATCH")
    positions = np.asarray(catalog_matcher.pixel_positions, dtype=float)
    if (positions.shape != (len(catalog_matcher.rows), 2) or
            not len(positions) or not np.isfinite(positions).all()):
        return dict(report, status="INVALID_NATIVE_CATALOG_POSITIONS")
    if (not isinstance(core_geometry, dict) or
            core_geometry.get("status") !=
            "FINITE_NATIVE_CORE_GEOMETRY_MEASURED"):
        return dict(report, status="NO_FROZEN_NATIVE_CORE_GEOMETRY")
    fwhm = core_geometry.get(
        "measured_core_fwhm_axis_geometric_mean_native_pix")
    if (not isinstance(fwhm, (int, float)) or
            not math.isfinite(fwhm) or fwhm <= 0):
        return dict(report, status="NO_FROZEN_NATIVE_CORE_FWHM")
    if protected_mask is not None:
        protected_mask = np.asarray(protected_mask, dtype=bool)
        if protected_mask.shape != frame.shape:
            return dict(report, status="INVALID_TARGET_PROTECTION_MASK_SHAPE")

    # Invoke EXACT P16C catalog eligibility and brightness ordering rules
    # with a larger, predeclared quota, preserving the original top sixteen.
    expanded = audit_disjoint_bright_gaia_wing_probes(
        frame, stars, catalog_matcher, core_geometry,
        protected_mask=protected_mask,
        max_catalog_probes=P17_COHORT_LIMIT)
    if (expanded.get("status") !=
            "DISJOINT_BRIGHT_GAIA_NATIVE_PROFILE_FEASIBILITY_ONLY" or
            expanded.get("native_catalog_sha256") !=
            p16c_qa["native_catalog_sha256"]):
        return dict(report, status="EXPANDED_CATALOG_SELECTION_FAILED")
    fixed = expanded["per_star"]
    prefix_keys = ("gaia_source_id", "tractor_objid", "native_xy_pix",
                   "native_tractor_flux_nMgy_catalog_not_direct_observed_peak",
                   "catalog_flux_rank_in_fixed_new_cohort")
    for a, b in zip(fixed[:P17_FROZEN_P16C_PREFIX],
                    p16c_qa["per_star"]):
        if any(a.get(key) != b.get(key) for key in prefix_keys):
            return dict(report, status="P16C_FIRST_16_IDENTITIES_CHANGED")
    if len(fixed) < P17_FROZEN_P16C_PREFIX:
        return dict(report, status="P16C_FIRST_16_IDENTITIES_MISSING")
    original = np.asarray(frame.image_nMgy, dtype=float)
    inverse = np.asarray(frame.invvar_nMgy2, dtype=float)
    half = int(math.ceil(SKY_RING_B[1] * fwhm)) + 1
    result_stars = []
    n_rows = len(catalog_matcher.rows)

    for selected in fixed:
        rank = int(selected["catalog_flux_rank_in_fixed_new_cohort"])
        category = ("P16C_LEGACY_EXPLORATORY" if rank <= 16
                    else "NEW_CATALOG_HOLDOUT" if rank % 4 == 0
                    else "NEW_CATALOG_EXPLORATORY")
        item = {
            "source_id": selected["source_id"],
            "gaia_source_id": selected["gaia_source_id"],
            "tractor_objid": selected["tractor_objid"],
            "catalog_rank": rank,
            "cohort_role": category,
            "native_xy_pix": selected["native_xy_pix"],
            "native_tractor_flux_nMgy": selected[
                "native_tractor_flux_nMgy_catalog_not_direct_observed_peak"],
            "predeclared_before_new_profile_reads": True,
            "status": "NOT_MEASURABLE",
            "mask_variants": [],
            "science_adopted": False,
        }
        if selected["status"] != (
                "MASKED_OUTER_PROFILE_DESCRIPTIVE_SKY_DEGENERATE"):
            item["status"] = "ORIGINAL_NATIVE_PEAK_OR_SKY_NOT_MEASURABLE"
            result_stars.append(item)
            continue
        peak = selected["original_star_peak_nMgy_per_native_pixel"]
        center_flux = item["native_tractor_flux_nMgy"]
        x, y = map(float, selected["native_xy_pix"])
        ix, iy = int(np.rint(x)), int(np.rint(y))
        x0, y0 = ix - half, iy - half
        if (x0 < 0 or y0 < 0 or ix + half >= frame.shape[1] or
                iy + half >= frame.shape[0]):
            item["status"] = "FULL_NATIVE_STAMP_OUTSIDE_IMAGE"
            result_stars.append(item)
            continue
        sl_y, sl_x = slice(y0, iy+half+1), slice(x0, ix+half+1)
        pixels = original[sl_y, sl_x]
        ivar = inverse[sl_y, sl_x]
        base = (np.asarray(frame.valid[sl_y, sl_x], dtype=bool).copy()
                & np.isfinite(pixels) & np.isfinite(ivar) & (ivar > 0))
        for key in (f"SATUR_{frame.band.upper()}",
                    f"ALLMASK_{frame.band.upper()}"):
            if key in frame.mask_reasons:
                base &= ~np.asarray(
                    frame.mask_reasons[key][sl_y, sl_x], dtype=bool)
        if protected_mask is not None:
            base &= ~protected_mask[sl_y, sl_x]
        gy, gx = np.indices(pixels.shape, dtype=float)
        gx += x0
        gy += y0
        dx, dy = (gx-x)/fwhm, (gy-y)/fwhm
        radius = np.hypot(gx-x, gy-y)
        sector = np.minimum(np.floor(
            (np.arctan2(dy, dx) % (2*np.pi))*8/(2*np.pi)
        ).astype(int), 7)
        dcat = np.hypot(positions[:, 0]-x, positions[:, 1]-y)
        if np.count_nonzero(dcat <= 2.5) != 1:
            item["status"] = "CATALOG_CENTRAL_MATCH_NOT_UNIQUE"
            result_stars.append(item)
            continue
        adaptive = _p17_adaptive_neighbor_radii(
            catalog_matcher.rows, center_flux, fwhm)
        item["n_catalog_rows_in_original_native_inventory"] = n_rows
        item["mask_variant_summaries_are_not_independent"] = True
        for label in P17_MASK_CHOICES:
            radii = (
                adaptive if label == "flux_tiered_catalog"
                else np.full(n_rows, {
                    "fixed_3_fwhm": 3.,
                    "fixed_4p5_fwhm": 4.5,
                    "fixed_6_fwhm": 6.,
                }[label] * fwhm))
            # Only catalog centroids close enough to intersect the stamp
            # are traversed; no star is excluded from the fixed cohort.
            maximum = 6.*fwhm
            nearby = np.flatnonzero(
                (dcat > 2.5) &
                (positions[:, 0] >= x0-maximum) &
                (positions[:, 0] <= ix+half+maximum) &
                (positions[:, 1] >= y0-maximum) &
                (positions[:, 1] <= iy+half+maximum))
            mask = _p17_neighbor_mask(
                pixels.shape, x0, y0, positions, nearby, radii)
            usable = base & mask
            sky_masks = [
                usable & (radius >= lo*fwhm) & (radius < hi*fwhm)
                for lo, hi in (SKY_RING_A, SKY_RING_B)
            ]
            sky_n = [int(np.count_nonzero(v)) for v in sky_masks]
            skies = [
                float(np.median(pixels[v])) if np.any(v) else None
                for v in sky_masks
            ]
            plane = _p16d_conditional_sky_plane(
                pixels, dx, dy, sky_masks[0] | sky_masks[1])
            plane_pixels = (plane[0] + plane[1]*dx + plane[2]*dy
                            if plane is not None else None)
            record = {
                "strategy": label,
                "n_catalog_neighbors_considered_in_cutout": int(
                    len(nearby)),
                "n_neighbor_masks_intersecting_cutout": int(
                    sum(
                        (max(x0, math.floor(positions[j, 0]-radii[j])) <=
                         min(ix+half, math.ceil(positions[j, 0]+radii[j])))
                        and
                        (max(y0, math.floor(positions[j, 1]-radii[j])) <=
                         min(iy+half, math.ceil(positions[j, 1]+radii[j])))
                        for j in nearby)),
                "sky_ring_n_unmasked": sky_n,
                "sky_ring_median_nMgy_per_native_pixel": skies,
                "sky_plane_coefficients_nMgy_per_native_pixel": (
                    list(plane) if plane is not None else None),
                "radial_profile": [],
            }
            for low, high in P16D_PROFILE_RINGS_FWHM:
                full = (radius >= low*fwhm) & (radius < high*fwhm)
                ring = full & usable
                n = int(np.count_nonzero(ring))
                total = int(np.count_nonzero(full))
                raw = float(np.median(pixels[ring])) if n else None
                octant_values = []
                for k in range(8):
                    area = ring & (sector == k)
                    octant_values.append(
                        float(np.median(pixels[area]-skies[1])/peak)
                        if skies[1] is not None and
                        np.count_nonzero(area) >= P16D_MIN_OCTANT_PIXELS
                        else None)
                populated = [v for v in octant_values if v is not None]
                record["radial_profile"].append({
                    "radius_inner_in_measured_fwhm": low,
                    "radius_outer_in_measured_fwhm": high,
                    "n_unmasked_native_pixels": n,
                    "unmasked_fraction": n / total if total else None,
                    "n_octants_with_at_least_25_pixels": len(populated),
                    "sky_a_subtracted_median_over_peak": (
                        (raw-skies[0])/peak if raw is not None and
                        skies[0] is not None else None),
                    "sky_b_subtracted_median_over_peak": (
                        (raw-skies[1])/peak if raw is not None and
                        skies[1] is not None else None),
                    "conditional_plane_subtracted_median_over_peak": (
                        float(np.median((pixels-plane_pixels)[ring])/peak)
                        if plane_pixels is not None and n else None),
                    "sky_b_subtracted_octant_medians_over_peak":
                        octant_values,
                    "max_minus_min_octant_median_over_peak": (
                        max(populated)-min(populated)
                        if len(populated) >= 4 else None),
                })
            item["mask_variants"].append(record)
        # Enforce EXACT legacy baseline continuity; do not accept a
        # "better" mask when the old 3-FWHM result was not reproduced.
        if rank <= P17_FROZEN_P16C_PREFIX:
            legacy = p16c_qa["per_star"][rank-1]
            first = item["mask_variants"][0]
            pairs = zip(first["sky_ring_median_nMgy_per_native_pixel"],
                        (legacy["conditional_sky_inner_ring"][
                            "median_original_nMgy_per_native_pixel"],
                         legacy["conditional_sky_outer_ring"][
                            "median_original_nMgy_per_native_pixel"]))
            agree = all(a is not None and np.isclose(
                a, b, rtol=0, atol=1.e-10) for a, b in pairs)
            old_rings = {
                (r["radius_inner_in_measured_fwhm"],
                 r["radius_outer_in_measured_fwhm"]): r
                for r in legacy["radial_profile"]}
            for r in first["radial_profile"]:
                old = old_rings.get((
                    r["radius_inner_in_measured_fwhm"],
                    r["radius_outer_in_measured_fwhm"]))
                measured = r["sky_b_subtracted_median_over_peak"]
                expected = (old.get(
                    "sky_b_subtracted_median_over_original_star_peak")
                    if old is not None else None)
                agree = bool(agree and measured is not None and
                             expected is not None and np.isclose(
                                 measured, expected, atol=1.e-10, rtol=0))
            item["reproduces_frozen_p16c_three_fwhm"] = bool(agree)
            if not agree:
                item["status"] = "P16C_BASELINE_MISMATCH"
                result_stars.append(item)
                continue
        item["status"] = "MULTIMASK_DESCRIPTIVE_ONLY"
        result_stars.append(item)

    return dict(
        report,
        status="FIXED_EXPANDED_COHORT_MASK_FEASIBILITY_ONLY",
        native_catalog_sha256=p16c_qa["native_catalog_sha256"],
        n_catalog_only_eligible_before_limit=expanded[
            "n_catalog_only_geometry_eligible_before_brightness_cap"],
        n_catalog_selected=len(fixed),
        n_legacy_p16c_stars=sum(
            v["cohort_role"] == "P16C_LEGACY_EXPLORATORY"
            for v in result_stars),
        n_new_exploratory=sum(
            v["cohort_role"] == "NEW_CATALOG_EXPLORATORY"
            for v in result_stars),
        n_new_identity_holdout=sum(
            v["cohort_role"] == "NEW_CATALOG_HOLDOUT"
            for v in result_stars),
        n_legacy_baseline_reproduced=sum(
            v.get("reproduces_frozen_p16c_three_fwhm") is True
            for v in result_stars),
        n_multimask_stars_measured=sum(
            v["status"] == "MULTIMASK_DESCRIPTIVE_ONLY"
            for v in result_stars),
        n_legacy_p16c_expected=P17_FROZEN_P16C_PREFIX,
        planned_new_cohort_has_32_stars=(
            len(fixed) == P17_COHORT_LIMIT),
        n_new_catalog_holdout_stars_is_not_cross_field_validation=True,
        no_outcome_driven_refilling=True,
        per_star=result_stars,
    )



# P18 is a new, predeclared single *diagnostic* morphology hypothesis.
# It does not retroactively select the most favorable P17 mask or modify
# the already frozen PSF training / 16-star P16C or 48-star P17 cohorts.
P18_SCHEMA = "desidera.psf0-catalog-shape-neighbor-sector-sky/2"
P18_POINT_CORE_FWHM = (3., 3.5, 4.5)   # ratio <0.1; 0.1..1; >=1
P18_EXTENDED_HALF_LIGHT_FACTOR = 4.
P18_MISSING_EXTENDED_SHAPE_FWHM = 6.
P18_SKY_SECTOR_MIN_PIXELS = 100


def _p18_neighbor_footprint(row, center_flux, fwhm, arcsec_per_pix):
    """Build a catalog-only bounded *hypothetical* ellipse in sky arcsec.

    The uniform PSF-sized core is combined in quadrature with four
    Tractor half-light radii of an extended model. This is NOT a
    direct PSF-convolved source model or a guarantee of zero flux
    beyond the exclusion ellipse.
    """
    flux = row.get("flux")
    typ = row.get("type")
    if (not isinstance(flux, (int, float)) or
            not math.isfinite(flux) or flux <= 0):
        return (P18_MISSING_EXTENDED_SHAPE_FWHM * fwhm *
                arcsec_per_pix,) * 2 + (0., "MISSING_FLUX_FALLBACK")
    ratio = flux / center_flux
    base = (P18_POINT_CORE_FWHM[2] if ratio >= 1.
            else P18_POINT_CORE_FWHM[1] if ratio >= .1
            else P18_POINT_CORE_FWHM[0])
    core = base * fwhm * arcsec_per_pix
    if typ == "PSF":
        return core, core, 0., "CATALOG_PSF"
    if typ not in ("REX", "EXP", "DEV", "SER"):
        radius = P18_MISSING_EXTENDED_SHAPE_FWHM * fwhm * arcsec_per_pix
        return radius, radius, 0., "UNKNOWN_MORPHOLOGY_FALLBACK"
    radius = row.get("shape_r_arcsec")
    e1 = row.get("shape_e1")
    e2 = row.get("shape_e2")
    if (not isinstance(radius, (int, float)) or
            not math.isfinite(radius) or radius <= 0):
        fallback = P18_MISSING_EXTENDED_SHAPE_FWHM*fwhm*arcsec_per_pix
        return fallback, fallback, 0., "MISSING_SHAPE_FALLBACK"
    if typ == "REX":
        e1, e2 = 0., 0.
    if (not isinstance(e1, (int, float)) or
            not isinstance(e2, (int, float)) or
            not math.isfinite(e1) or not math.isfinite(e2)):
        fallback = P18_MISSING_EXTENDED_SHAPE_FWHM*fwhm*arcsec_per_pix
        return fallback, fallback, 0., "INVALID_ELLIPTICITY_FALLBACK"
    eps = math.hypot(e1, e2)
    if eps >= 1.:
        fallback = P18_MISSING_EXTENDED_SHAPE_FWHM*fwhm*arcsec_per_pix
        return fallback, fallback, 0., "INVALID_ELLIPTICITY_FALLBACK"
    q = (1. - eps) / (1. + eps)
    theta = .5*math.atan2(e2, e1)
    wing = P18_EXTENDED_HALF_LIGHT_FACTOR*radius
    return (math.hypot(core, wing), math.hypot(core, q*wing),
            theta, "EXTENDED_SHAPE_ELLIPSE")


def _p18_sky_to_pixel_matrix(frame, x, y):
    """Local native-pixel to east/north tangent Jacobian (arcsec/pixel).

    For TAN images the WCS CD/PC matrix is already defined in the
    intermediate projection plane, NOT in raw RA/DEC sky increments.
    Multiplying that matrix by cos(dec) again biases the east scale.
    Instead finite-difference the full native WCS near the source,
    including SIP distortion when present, then convert the local RA
    derivative to east with cos(dec) exactly once. At high declination,
    near WCS discontinuities or with ill-conditioned Jacobians, fail.
    """
    if not hasattr(frame, "wcs"):
        return None
    try:
        coords = np.asarray(frame.wcs.all_pix2world(
            np.array([x, x+1., x-1., x, x]),
            np.array([y, y, y, y+1., y-1.]), 0),dtype=float)
        if coords.shape != (2,5) or not np.isfinite(coords).all():
            return None
        dec0 = float(coords[1,0])
        cosdec = math.cos(math.radians(dec0))
        if abs(cosdec) < 1.e-5:
            return None
        d_ra = (coords[0,[1,3]]-coords[0,[2,4]]+180.)%360.-180.
        d_dec = coords[1,[1,3]]-coords[1,[2,4]]
        mat = np.stack([d_ra*cosdec,d_dec])*1800.
    except (AttributeError, TypeError, ValueError, IndexError):
        return None
    if (mat.shape != (2, 2) or not np.isfinite(mat).all() or
            np.linalg.det(mat) == 0):
        return None
    sig = np.linalg.svd(mat, compute_uv=False)
    if sig[-1] <= 0 or sig[0]/sig[-1] > 2:
        return None
    return mat


def _p18_shape_mask(shape, x0, y0, positions, rows, center_index,
                    center_flux, fwhm, tangent_matrix):
    """Compute local shape+brightness masks without star-pixel decisions."""
    usable = np.ones(shape, dtype=bool)
    height, width = shape
    singular = np.linalg.svd(tangent_matrix, compute_uv=False)
    min_scale = float(singular[-1])
    scale = float(np.sqrt(abs(np.linalg.det(tangent_matrix))))
    counts = {}
    n_intersecting = 0
    for j, (row, xy) in enumerate(zip(rows, positions)):
        if j == center_index:
            continue
        cx, cy = map(float, xy)
        a, b, theta, reason = _p18_neighbor_footprint(
            row, center_flux, fwhm, scale)
        margin = a/min_scale
        left = max(0, int(math.floor(cx-margin))-x0)
        right = min(width, int(math.ceil(cx+margin))-x0+1)
        bottom = max(0, int(math.floor(cy-margin))-y0)
        top = min(height, int(math.ceil(cy+margin))-y0+1)
        if left >= right or bottom >= top:
            continue
        n_intersecting += 1
        counts[reason] = counts.get(reason, 0) + 1
        yy, xx = np.ogrid[y0+bottom:y0+top, x0+left:x0+right]
        east = tangent_matrix[0,0]*(xx-cx)+tangent_matrix[0,1]*(yy-cy)
        north = tangent_matrix[1,0]*(xx-cx)+tangent_matrix[1,1]*(yy-cy)
        cc, ss = math.cos(theta), math.sin(theta)
        major = cc*east + ss*north
        minor = -ss*east + cc*north
        usable[bottom:top, left:right] &= (
            (major/a)**2 + (minor/b)**2 > 1.)
    return usable, counts, n_intersecting


def _p18_sky_sector_deletion(pixels, usable, radius, octants, fwhm,
                             peak):
    """Eight leave-one-octant-out sky-B medians: correlated sensitivity."""
    ring_b = usable & (radius >= SKY_RING_B[0]*fwhm) & (
        radius < SKY_RING_B[1]*fwhm)
    full = (float(np.median(pixels[ring_b]))
            if np.count_nonzero(ring_b) >= P16D_MIN_SKY_PIXELS
            else None)
    deletion = []
    for s in range(8):
        remaining = ring_b & (octants != s)
        value = (float(np.median(pixels[remaining]))
                 if np.count_nonzero(remaining) >= P16D_MIN_SKY_PIXELS
                 else None)
        deletion.append(value)
    if full is None or any(v is None for v in deletion):
        worst = None
    else:
        worst = float(max(abs(v-full) for v in deletion)/peak)
    return full, deletion, worst


def audit_catalog_shape_neighbor_and_sky_sensitivity(
    frame, catalog_matcher, core_geometry, p17_qa, *,
    protected_mask=None
):
    """P18: no re-ranking, no PSF fit, no sky error bars.

    Uses exact P17 fixed 48 star IDs and the native Tractor morphology
    already read from that same DR11 brick. Records an additional
    conditional mask profile plus disjoint sky and leave-one-octant-out
    diagnostics; selection and mask geometry never use wing residuals.
    """
    report = {
        "schema": P18_SCHEMA,
        "status": "NOT_RUN",
        "science_adopted": False,
        "cohort_fixed_from_p17": True,
        "geometry_uses_original_catalog_only": True,
        "center_fwhm_uses_frozen_measured_core": True,
        "point_core_fwhm_by_flux_ratio": list(P18_POINT_CORE_FWHM),
        "point_flux_ratio_edges": [0.1, 1.],
        "extended_shape_half_light_radius_multiplier":
            P18_EXTENDED_HALF_LIGHT_FACTOR,
        "missing_or_uninterpretable_shape_fallback_fwhm":
            P18_MISSING_EXTENDED_SHAPE_FWHM,
        "ellipticity_convention": "(1-q)/(1+q), angle in tangent east/north",
        "local_wcs_jacobian": "FULL_ALL_PIX2WORLD_CENTRAL_DIFFERENCE_EAST_NORTH_V2",
        "sky_sector_deletions_correlated_not_error_bars": True,
        "sky_model_not_validated": True,
        "coadd_noise_covariance_not_calibrated": True,
        "extended_model_not_psf_convolved_or_flux_free_outside_mask": True,
        "no_posthoc_mask_selection_no_psf_or_science_adoption": True,
        "per_star": [],
    }
    if (not isinstance(p17_qa, dict) or
            p17_qa.get("schema") != P17_SCHEMA or
            p17_qa.get("status") !=
            "FIXED_EXPANDED_COHORT_MASK_FEASIBILITY_ONLY" or
            len(p17_qa.get("per_star", [])) != P17_COHORT_LIMIT):
        return dict(report, status="NO_VALID_FIXED_P17_COHORT")
    if (catalog_matcher is None or
            getattr(catalog_matcher, "source_sha256", None) !=
            p17_qa.get("native_catalog_sha256") or
            not hasattr(catalog_matcher, "rows") or
            not hasattr(catalog_matcher, "pixel_positions")):
        return dict(report, status="CATALOG_LINEAGE_MISMATCH")
    positions = np.asarray(catalog_matcher.pixel_positions, dtype=float)
    rows = catalog_matcher.rows
    if (positions.shape != (len(rows), 2) or
            not np.isfinite(positions).all()):
        return dict(report, status="INVALID_NATIVE_CATALOG_POSITIONS")
    fwhm = (core_geometry.get(
        "measured_core_fwhm_axis_geometric_mean_native_pix")
        if isinstance(core_geometry, dict) and core_geometry.get("status") ==
        "FINITE_NATIVE_CORE_GEOMETRY_MEASURED" else None)
    if (not isinstance(fwhm, (int, float)) or
            not math.isfinite(fwhm) or fwhm <= 0):
        return dict(report, status="NO_FROZEN_MEASURED_FWHM")
    if protected_mask is not None:
        protected_mask = np.asarray(protected_mask, dtype=bool)
        if protected_mask.shape != frame.shape:
            return dict(report, status="INVALID_TARGET_PROTECTION_MASK")
    original = np.asarray(frame.image_nMgy, dtype=float)
    inverse = np.asarray(frame.invvar_nMgy2, dtype=float)
    half = int(math.ceil(SKY_RING_B[1]*fwhm))+1
    output = []
    n_native_sky = 0
    for star in p17_qa["per_star"]:
        item = {
            "source_id": star.get("source_id"),
            "gaia_source_id": star.get("gaia_source_id"),
            "catalog_rank": star.get("catalog_rank"),
            "cohort_role": star.get("cohort_role"),
            "status": "NOT_MEASURABLE",
            "science_adopted": False,
        }
        # A P17 error MUST remain an error; no sample refilling.
        if star.get("status") != "MULTIMASK_DESCRIPTIVE_ONLY":
            item["status"] = "P17_ORIGINAL_PROFILE_NOT_MEASURABLE"
            output.append(item)
            continue
        x, y = map(float, star["native_xy_pix"])
        ix, iy = int(np.rint(x)), int(np.rint(y))
        x0, y0 = ix-half, iy-half
        if (x0 < 0 or y0 < 0 or ix+half >= frame.shape[1] or
                iy+half >= frame.shape[0]):
            item["status"] = "FULL_NATIVE_STAMP_OUTSIDE_IMAGE"
            output.append(item)
            continue
        tangent = _p18_sky_to_pixel_matrix(frame, x, y)
        if tangent is None:
            item["status"] = "MISSING_OR_UNSAFE_WCS_TANGENT_MATRIX"
            output.append(item)
            continue
        matching = [
            j for j, row in enumerate(rows)
            if row.get("ref_cat") == "G3" and
            row.get("ref_id") == star["gaia_source_id"] and
            row.get("objid") == star["tractor_objid"]]
        if len(matching) != 1 or np.linalg.norm(
                positions[matching[0]]-np.array([x,y])) > 1.e-6:
            item["status"] = "P17_CATALOG_STAR_IDENTITY_MISMATCH"
            output.append(item)
            continue
        sl_y, sl_x = slice(y0, iy+half+1), slice(x0, ix+half+1)
        pixels = original[sl_y, sl_x]
        inverse_cut = inverse[sl_y, sl_x]
        usable = (np.asarray(frame.valid[sl_y,sl_x], dtype=bool).copy()
                  & np.isfinite(pixels) & np.isfinite(inverse_cut)
                  & (inverse_cut > 0))
        for key in (f"SATUR_{frame.band.upper()}",
                    f"ALLMASK_{frame.band.upper()}"):
            if key in frame.mask_reasons:
                usable &= ~np.asarray(frame.mask_reasons[key][
                    sl_y,sl_x], dtype=bool)
        if protected_mask is not None:
            usable &= ~protected_mask[sl_y,sl_x]
        catalog_mask, reasons, n_neighbors = _p18_shape_mask(
            pixels.shape, x0, y0, positions, rows, matching[0],
            float(star["native_tractor_flux_nMgy"]), float(fwhm), tangent)
        usable &= catalog_mask
        yy, xx = np.indices(pixels.shape, dtype=float)
        xx += x0
        yy += y0
        dx, dy = (xx-x)/fwhm, (yy-y)/fwhm
        radius = np.hypot(xx-x, yy-y)
        octants = np.minimum(np.floor(
            (np.arctan2(dy,dx)%(2*np.pi))*8/(2*np.pi)
        ).astype(int), 7)
        sky_masks = [
            usable & (radius >= a*fwhm) & (radius < b*fwhm)
            for a,b in (SKY_RING_A, SKY_RING_B)
        ]
        sky_counts = [int(np.count_nonzero(v)) for v in sky_masks]
        skies = [
            float(np.median(pixels[v]))
            if np.count_nonzero(v) >= P16D_MIN_SKY_PIXELS else None
            for v in sky_masks
        ]
        plane = _p16d_conditional_sky_plane(
            pixels, dx, dy, sky_masks[0] | sky_masks[1])
        plane_field = (plane[0]+plane[1]*dx+plane[2]*dy
                       if plane is not None else None)
        # P17's center is defined by exactly the same unchanged core
        # and 3x3 native peak, not a newly optimized amplitude.
        peak = float(star["native_tractor_flux_nMgy"])
        # The P17 residual denominators use native peak, not Tractor
        # integrated flux. Derive it from the *unchanged* original 3x3.
        center = original[iy-1:iy+2,ix-1:ix+2]
        center_valid = (
            np.asarray(frame.valid[iy-1:iy+2,ix-1:ix+2], dtype=bool)
            & np.isfinite(center)
            & np.isfinite(inverse[iy-1:iy+2,ix-1:ix+2])
            & (inverse[iy-1:iy+2,ix-1:ix+2] > 0))
        for key in (f"SATUR_{frame.band.upper()}",
                    f"ALLMASK_{frame.band.upper()}"):
            if key in frame.mask_reasons:
                center_valid &= ~np.asarray(frame.mask_reasons[key][
                    iy-1:iy+2,ix-1:ix+2], dtype=bool)
        if protected_mask is not None:
            center_valid &= ~protected_mask[iy-1:iy+2,ix-1:ix+2]
        if not np.any(center_valid):
            item["status"] = "INVALID_ORIGINAL_NATIVE_PEAK"
            output.append(item)
            continue
        peak = float(np.max(center[center_valid]))
        if not math.isfinite(peak) or peak <= 0:
            item["status"] = "INVALID_ORIGINAL_NATIVE_PEAK"
            output.append(item)
            continue
        sky_b, deletion, maximum = _p18_sky_sector_deletion(
            pixels, usable, radius, octants, fwhm, peak)
        item.update({
            "n_neighbor_footprints_intersecting_stamp": n_neighbors,
            "catalog_footprint_class_counts": reasons,
            "native_peak_nMgy_per_pixel": peak,
            "sky_ring_n_unmasked": sky_counts,
            "sky_ring_medians_nMgy_per_pixel": skies,
            "sky_plane_coefficients_nMgy_per_pixel": (
                list(plane) if plane is not None else None),
            "sky_b_leave_one_octant_out_medians_nMgy": deletion,
            "sky_b_max_abs_leave_one_octant_change_over_peak": maximum,
            "sky_a_minus_sky_b_over_peak": (
                (skies[0]-skies[1])/peak if skies[0] is not None and
                skies[1] is not None else None),
            "radial_profile": [],
        })
        for a,b in P16D_PROFILE_RINGS_FWHM:
            full = (radius >= a*fwhm)&(radius < b*fwhm)
            usable_ring = full & usable
            n = int(np.count_nonzero(usable_ring))
            total = int(np.count_nonzero(full))
            raw = float(np.median(pixels[usable_ring])) if n else None
            octant_values = []
            for k in range(8):
                sel = usable_ring & (octants == k)
                octant_values.append(
                    float(np.median((pixels-skies[1])[sel])/peak)
                    if skies[1] is not None and
                    np.count_nonzero(sel) >= P16D_MIN_OCTANT_PIXELS
                    else None)
            item["radial_profile"].append({
                "radius_inner_in_measured_fwhm": a,
                "radius_outer_in_measured_fwhm": b,
                "n_unmasked_native_pixels": n,
                "unmasked_fraction": n/total if total else None,
                "n_octants_with_at_least_25_pixels": sum(
                    v is not None for v in octant_values),
                "sky_b_subtracted_median_over_peak": (
                    (raw-skies[1])/peak
                    if raw is not None and skies[1] is not None else None),
                "sky_a_subtracted_median_over_peak": (
                    (raw-skies[0])/peak
                    if raw is not None and skies[0] is not None else None),
                "conditional_plane_subtracted_median_over_peak": (
                    float(np.median((pixels-plane_field)[usable_ring])/peak)
                    if plane_field is not None and n else None),
                "sky_b_subtracted_octant_medians_over_peak": octant_values,
            })
        if (sky_counts[0] < P16D_MIN_SKY_PIXELS or
                sky_counts[1] < P16D_MIN_SKY_PIXELS or
                plane is None):
            item["status"] = "INSUFFICIENT_OR_UNSAFE_CONDITIONAL_SKY"
        else:
            item["status"] = "MORPHOLOGY_MASK_AND_SKY_DESCRIPTIVE_ONLY"
            n_native_sky += 1
        output.append(item)
    return dict(
        report,
        status="FIXED_P17_COHORT_MORPHOLOGY_SKY_DIAGNOSTIC_ONLY",
        native_catalog_sha256=p17_qa["native_catalog_sha256"],
        n_fixed_star_identities=len(output),
        n_new_identity_holdout=sum(
            v["cohort_role"] == "NEW_CATALOG_HOLDOUT" for v in output),
        n_profiles_with_both_sky_and_plane=n_native_sky,
        n_profiles_not_qualifying_sky_or_geometry=len(output)-n_native_sky,
        no_star_replacement_or_psf_adoption=True,
        per_star=output,
    )
