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
    frame, stars, catalog_matcher, core_geometry, *, protected_mask=None
):
    """Catalog-preregistered bright-star wing *feasibility* sample.

    Selection uses only the original native Tractor row, coordinates, the
    pre-existing target exclusion geometry and the *frozen* core FWHM. Never
    select/refill on observed flux, radial profile, measured background, mask
    coverage or perceived scientific merit. A selected star with bad original
    pixels stays selected but is returned as failed, never replaced.
    """
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
        "predeclared_max_catalog_star_probes": PROSPECTIVE_LIMIT,
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
    fixed = candidates[:PROSPECTIVE_LIMIT]
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
        radial = np.hypot(dx, dy)
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
            sky_a_mask = usable & (radial >= 22) & (radial < 26)
            sky_b_mask = usable & (radial >= 26) & (radial < 30)
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
                area = (radial >= inner) & (radial < outer)
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
            if all(v["sky_ring_b_median_nMgy_per_pixel"] is not None
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
