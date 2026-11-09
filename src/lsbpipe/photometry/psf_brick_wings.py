"""P20: same-native-brick stellar wings, blank-sky controls and spatial checks.

Read-only fit of a common *conditional* radial profile plus one local
sky nuisance offset per bright star. It is not a calibrated PSF tail,
a physical scattering kernel, or a science-adoption mechanism.
All pixel controls and both fixed identity splits stay within one brick.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.optimize import least_squares

from lsbpipe.photometry.psf_wing_profile import (
    P16D_MIN_SKY_PIXELS, P16D_MIN_OCTANT_PIXELS,
    P16D_PROFILE_RINGS_FWHM, P18_SCHEMA,
    P17_SCHEMA, SKY_RING_A, SKY_RING_B,
    _p18_shape_mask, _p18_sky_to_pixel_matrix,
)

SCHEMA = "desidera.psf0-same-brick-joint-wing-sky-spatial/1"
BRIGHT_TRAIN_RANKS = tuple(x for x in range(1, 17) if x % 4)
BRIGHT_RANK_HOLDOUT = (4, 8, 12, 16)
BLANK_SPACING_FACTOR = 2.
BLANK_MAX_PER_QUADRANT = 4
BLANK_MIN_CENTRAL_CATALOG_SEPARATION_FWHM = 3.
BLANK_MIN_BRIGHT_STAR_SEPARATION_FWHM = 40.
MIN_RADIAL_PIXELS = 100
MIN_RADIAL_OCTANTS = 4


def _med(values, *, absolute=False):
    clean = [abs(float(v)) if absolute else float(v)
             for v in values if v is not None and np.isfinite(v)]
    return float(np.median(clean)) if clean else None


def _quadrant(xy, shape):
    return int(xy[0] >= shape[1]/2.) + 2*int(xy[1] >= shape[0]/2.)


def _within_brick_blank_controls(frame, matcher, stars, fwhm,
                                 protected_mask):
    """Choose *only by catalog geometry* on a fixed, nonoverlap grid.

    Invalid sky/mask data is recorded without replacing that center.
    The grid, fixed quota and eligibility are independent of native
    image residuals. The same coadd couples controls and stars.
    """
    positions = np.asarray(matcher.pixel_positions, dtype=float)
    rows = matcher.rows
    bright = np.asarray([
        s["native_xy_pix"] for s in stars
        if s["catalog_rank"] <= 16], dtype=float)
    height, width = frame.shape
    half = int(math.ceil(SKY_RING_B[1]*fwhm))+1
    step = int(math.ceil(BLANK_SPACING_FACTOR*half))+3
    x_grid = range(half+1, width-half-1, step)
    y_grid = range(half+1, height-half-1, step)
    quotas = [0]*4
    controls = []
    candidates = 0
    eligible = 0
    # Fixed reference from the preselected bright *catalog* (not pixel flux).
    reference_flux = float(np.median([
        s["native_tractor_flux_nMgy"] for s in stars
        if s["catalog_rank"] <= 16]))
    image = np.asarray(frame.image_nMgy)
    ivar = np.asarray(frame.invvar_nMgy2)
    for y in y_grid:
        for x in x_grid:
            candidates += 1
            q = _quadrant((x,y), frame.shape)
            if quotas[q] >= BLANK_MAX_PER_QUADRANT:
                continue
            if np.any(np.sum((positions-(x,y))**2, axis=1) <
                      (BLANK_MIN_CENTRAL_CATALOG_SEPARATION_FWHM*fwhm)**2):
                continue
            if np.any(np.sum((bright-(x,y))**2, axis=1) <
                      (BLANK_MIN_BRIGHT_STAR_SEPARATION_FWHM*fwhm)**2):
                continue
            if protected_mask is not None and protected_mask[y,x]:
                continue
            # Do not use pixels, SKY medians, or faint residual amplitudes
            # to choose or backfill a site.
            eligible += 1
            quotas[q] += 1
            entry = {
                "center_native_xy_pix": [int(x),int(y)],
                "brick_quadrant": q,
                "status": "NOT_MEASURABLE",
                "radial_residual_nMgy_per_pixel": [None]*6,
                "n_radial_unmasked_native_pixels": [0]*6,
                "radial_unmasked_fraction": [None]*6,
            }
            x0,y0 = x-half,y-half
            sx,sy = slice(x0,x+half+1), slice(y0,y+half+1)
            tangent = _p18_sky_to_pixel_matrix(frame, x, y)
            if tangent is None:
                entry["status"] = "WCS_INVALID"
                controls.append(entry)
                continue
            pixels = image[sy,sx]
            inverse = ivar[sy,sx]
            usable = (
                np.asarray(frame.valid[sy,sx],dtype=bool).copy()
                & np.isfinite(pixels) & np.isfinite(inverse)
                & (inverse>0))
            for key in (f"SATUR_{frame.band.upper()}",
                        f"ALLMASK_{frame.band.upper()}"):
                if key in frame.mask_reasons:
                    usable &= ~np.asarray(frame.mask_reasons[key][sy,sx],
                                          dtype=bool)
            if protected_mask is not None:
                usable &= ~protected_mask[sy,sx]
            neighbor_mask, reasons, _n = _p18_shape_mask(
                pixels.shape, x0, y0, positions, rows, -1,
                reference_flux, fwhm, tangent)
            usable &= neighbor_mask
            yy,xx = np.indices(pixels.shape,dtype=float)
            radius = np.hypot(xx+x0-x,yy+y0-y)
            sky_a = usable & (radius>=SKY_RING_A[0]*fwhm) & (
                radius<SKY_RING_A[1]*fwhm)
            sky_b = usable & (radius>=SKY_RING_B[0]*fwhm) & (
                radius<SKY_RING_B[1]*fwhm)
            na,nb = int(sky_a.sum()),int(sky_b.sum())
            entry["sky_n_unmasked_native_pixels"] = [na,nb]
            entry["catalog_mask_reasons"] = reasons
            if na<P16D_MIN_SKY_PIXELS or nb<P16D_MIN_SKY_PIXELS:
                entry["status"] = "NO_BOTH_NATIVE_SKY_RINGS"
                controls.append(entry)
                continue
            sa,sb = float(np.median(pixels[sky_a])),float(np.median(pixels[sky_b]))
            entry["sky_a_minus_sky_b_nMgy_per_pixel"] = sa-sb
            # One blank control has no central star by definition.
            # Its false "wing" is the same radial-minus-outer-sky
            # measurement in physical image units as at bright stars.
            residuals,counts,fractions = [],[],[]
            for low,high in P16D_PROFILE_RINGS_FWHM:
                full = (radius>=low*fwhm)&(radius<high*fwhm)
                selected = full & usable
                n = int(selected.sum())
                counts.append(n)
                fractions.append(float(n/full.sum()) if full.any() else None)
                residuals.append(
                    float(np.median(pixels[selected])-sb)
                    if n>=MIN_RADIAL_PIXELS else None)
            entry["radial_residual_nMgy_per_pixel"] = residuals
            entry["n_radial_unmasked_native_pixels"] = counts
            entry["radial_unmasked_fraction"] = fractions
            entry["status"] = "BLANK_CONDITIONAL_PROFILE_ONLY"
            controls.append(entry)
    return {
        "grid_step_native_pix": step,
        "n_grid_centers_visited": candidates,
        "n_catalog_geometry_eligible_before_quota": eligible,
        "n_catalog_fixed_blank_controls": len(controls),
        "n_measured_blank_controls": sum(
            x["status"]=="BLANK_CONDITIONAL_PROFILE_ONLY" for x in controls),
        "n_selected_by_quadrant": quotas,
        "reference_flux_nMgy_is_bright_catalog_median": reference_flux,
        "centers_chosen_before_any_pixel_or_residual_tests": True,
        "blank_controls_retain_unknown_uncatalogued_objects": True,
        "controls_same_coadd_not_statistically_independent": True,
        "per_control": controls,
        "per_ring": [
            {
                "radii_in_measured_fwhm": list(ring),
                "n_available_controls": sum(
                    row["radial_residual_nMgy_per_pixel"][idx] is not None
                    for row in controls),
                "median_signed_false_wing_nMgy_per_pixel":
                    _med([row["radial_residual_nMgy_per_pixel"][idx]
                          for row in controls]),
                "median_abs_false_wing_nMgy_per_pixel":
                    _med([row["radial_residual_nMgy_per_pixel"][idx]
                          for row in controls],absolute=True),
                "max_abs_false_wing_nMgy_per_pixel":
                    max((abs(row["radial_residual_nMgy_per_pixel"][idx])
                         for row in controls
                         if row["radial_residual_nMgy_per_pixel"][idx] is not None),
                        default=None),
            }
            for idx,ring in enumerate(P16D_PROFILE_RINGS_FWHM)
        ],
    }


def _joint_nuisance_fit(stars, blank_ring_abs, hypothesis):
    """Robust finite-bin common profile + constrained per-star sky offsets.

    Nuisance scales encode sensitivity across *pre-existing* background
    hypotheses and blank controls. They are not Gaussian 1-sigma errors.
    A unit penalty anchors otherwise unidentifiable additive degeneracy;
    the fitted absolute offset is accordingly prior-dependent.
    """
    observations,scales,star_indices,ring_indices = [],[],[],[]
    star_scales=[]
    for i,star in enumerate(stars):
        peak=star.get("native_peak_nMgy_per_pixel")
        sky_a=star.get("sky_a_minus_sky_b_over_peak")
        sector=star.get("sky_b_max_abs_leave_one_octant_change_over_peak")
        floor=max(
            abs(float(sky_a)) if sky_a is not None else 0.,
            float(sector) if sector is not None else 0.,
            1.e-8)
        for k,entry in enumerate(star.get("radial_profile", [])):
            amount=entry.get(hypothesis)
            if amount is None or not np.isfinite(amount):
                continue
            if (entry.get("n_unmasked_native_pixels",0)<MIN_RADIAL_PIXELS
                    or entry.get("n_octants_with_at_least_25_pixels",0)
                    <MIN_RADIAL_OCTANTS):
                continue
            null=blank_ring_abs[k]
            scale=max(floor,float(null)/peak if null is not None else 0.)
            observations.append(float(amount))
            scales.append(scale)
            star_indices.append(i)
            ring_indices.append(k)
        star_scales.append(floor)
    nstar=len(stars)
    if len(set(star_indices))<3 or len(set(ring_indices))!=6:
        return {"status":"INSUFFICIENT_FIXED_BRIGHT_STAR_RADIAL_COVERAGE",
                "n_measured_points":len(observations),
                "n_active_stars":len(set(star_indices))}
    y=np.asarray(observations)
    w=np.asarray(scales)
    si=np.asarray(star_indices)
    ri=np.asarray(ring_indices)
    prior=np.asarray(star_scales)
    seed=np.array([np.median(y[ri==j]) for j in range(6)]+[0.]*nstar)
    def residuals(parameters):
        return np.concatenate([
            (y-parameters[ri]-parameters[6+si])/w,
            parameters[6:]/prior,
        ])
    fit=least_squares(residuals,seed,loss="soft_l1",f_scale=1.,
                      max_nfev=300,ftol=1.e-9,xtol=1.e-9,gtol=1.e-9)
    if not fit.success or not np.isfinite(fit.x).all():
        return {"status":"JOINT_NUMERICAL_SOLVER_DID_NOT_CONVERGE",
                "n_measured_points":len(y)}
    return {
        "status":"ROBUST_COMMON_PROFILE_PLUS_SKY_NUISANCE_DESCRIPTIVE_ONLY",
        "n_input_stars":nstar,
        "n_active_stars":len(set(star_indices)),
        "n_measured_points":len(y),
        "radial_common_wing_candidate_over_native_peak":
            [float(v) for v in fit.x[:6]],
        "per_star_additive_background_nuisance_over_native_peak":
            [float(v) for v in fit.x[6:]],
        "nuisance_penalty_anchors_nonidentifiable_absolute_level":True,
        "weight_scales_are_systematic_sensitivities_not_sigmas":True,
        "coadd_noise_covariance_not_calibrated":True,
    }


def _predictions_against_stars(model, stars):
    if model.get("status") != (
        "ROBUST_COMMON_PROFILE_PLUS_SKY_NUISANCE_DESCRIPTIVE_ONLY"):
        return {"status":"NO_FIXED_TRAINING_MODEL"}
    amplitudes=model["radial_common_wing_candidate_over_native_peak"]
    rows=[]
    for star in stars:
        differences=[]
        conditional_scales=[]
        for k,r in enumerate(star["radial_profile"]):
            y=r.get("sky_b_subtracted_median_over_peak")
            if y is None or r.get("n_unmasked_native_pixels",0)<MIN_RADIAL_PIXELS:
                continue
            if r.get("n_octants_with_at_least_25_pixels",0)<MIN_RADIAL_OCTANTS:
                continue
            differences.append(float(y)-amplitudes[k])
            conditional_scales.append(max(
                abs(star.get("sky_a_minus_sky_b_over_peak") or 0.),
                star.get("sky_b_max_abs_leave_one_octant_change_over_peak") or 0.,
            ))
        rows.append({
            "catalog_rank":star["catalog_rank"],
            "gaia_source_id":star["gaia_source_id"],
            "n_comparable_rings":len(differences),
            "median_abs_prediction_residual_over_native_peak":
                _med(differences,absolute=True),
            "median_conditional_sky_sensitivity_over_native_peak":
                _med(conditional_scales),
            "no_prediction_based_sample_replacement":True,
        })
    return {
        "status":"PREDICTED_FIXED_WITHHELD_IDENTITIES_DESCRIPTIVE_ONLY",
        "n_withheld_stars":len(rows),
        "n_with_at_least_one_comparable_ring":sum(
            x["n_comparable_rings"]>0 for x in rows),
        "median_abs_prediction_residual_over_native_peak":_med([
            r["median_abs_prediction_residual_over_native_peak"] for r in rows]),
        "per_star":rows,
        "withheld_stars_share_same_coadd_and_historical_exploration":True,
    }


def audit_brick_joint_wings_sky_spatial(frame,catalog_matcher,
                                        core_geometry,p17_qa,p18_qa,*,
                                        protected_mask=None):
    """P20 within-brick controls; never mutates image, PSF or science gate."""
    out={
        "schema":SCHEMA,
        "status":"NOT_RUN",
        "science_adopted":False,
        "psf_wings_validated":False,
        "wide_angle_kernel_constructed":False,
        "dataset_id":getattr(frame,"dataset_id",None),
        "native_brick_only":True,
        "no_galaxy_or_star_injection_performed":True,
        "predeclared_legacy_training_ranks":list(BRIGHT_TRAIN_RANKS),
        "predeclared_legacy_identity_holdout_ranks":list(BRIGHT_RANK_HOLDOUT),
        "bright_identity_holdout_retrospective_not_blind":True,
        "no_mask_radius_or_sky_hypothesis_selected_by_residuals":True,
        "no_psf_tail_without_photometric_validation":True,
    }
    if (not isinstance(p17_qa,dict) or p17_qa.get("schema")!=P17_SCHEMA
            or p17_qa.get("status")!=
            "FIXED_EXPANDED_COHORT_MASK_FEASIBILITY_ONLY"
            or not isinstance(p18_qa,dict)
            or p18_qa.get("schema")!=P18_SCHEMA
            or p18_qa.get("status")!=
            "FIXED_P17_COHORT_MORPHOLOGY_SKY_DIAGNOSTIC_ONLY"):
        return dict(out,status="P17_P18_FIXED_COHORT_UNAVAILABLE")
    if (p17_qa.get("native_catalog_sha256")!=
            p18_qa.get("native_catalog_sha256") or
            catalog_matcher is None or
            getattr(catalog_matcher,"source_sha256",None)!=
            p17_qa.get("native_catalog_sha256")):
        return dict(out,status="NATIVE_CATALOG_LINEAGE_MISMATCH")
    old=p17_qa.get("per_star",[])
    new=p18_qa.get("per_star",[])
    if len(old)!=48 or len(new)!=48 or any(
        (a.get("catalog_rank"),a.get("gaia_source_id"),a.get("source_id"),
         a.get("cohort_role")) !=
        (b.get("catalog_rank"),b.get("gaia_source_id"),b.get("source_id"),
         b.get("cohort_role")) for a,b in zip(old,new)
    ):
        return dict(out,status="FROZEN_48_STAR_IDENTITY_MISMATCH")
    if (not isinstance(core_geometry,dict) or core_geometry.get("status") !=
            "FINITE_NATIVE_CORE_GEOMETRY_MEASURED"):
        return dict(out,status="FROZEN_CORE_UNAVAILABLE")
    fwhm=core_geometry.get(
        "measured_core_fwhm_axis_geometric_mean_native_pix")
    if (not isinstance(fwhm,(int,float)) or
            not np.isfinite(fwhm) or fwhm<=0):
        return dict(out,status="FROZEN_FWHM_UNAVAILABLE")
    if protected_mask is not None:
        protected_mask=np.asarray(protected_mask,dtype=bool)
        if protected_mask.shape!=frame.shape:
            return dict(out,status="TARGET_MASK_SHAPE_MISMATCH")
    positions=np.asarray(catalog_matcher.pixel_positions,dtype=float)
    if positions.shape!=(len(catalog_matcher.rows),2) or not np.isfinite(
            positions).all():
        return dict(out,status="NATIVE_CATALOG_POSITIONS_INVALID")
    if not hasattr(frame,"wcs"):
        return dict(out,status="NATIVE_WCS_MISSING")
    if set(s["catalog_rank"] for s in old)!=set(range(1,49)):
        return dict(out,status="NONCONTIGUOUS_OR_DUPLICATE_RANKS")
    by_rank={s["catalog_rank"]:s for s in new}
    if any(by_rank[r]["cohort_role"]!="P16C_LEGACY_EXPLORATORY"
           for r in range(1,17)):
        return dict(out,status="BRIGHT_COHORT_ROLE_MISMATCH")
    blanks=_within_brick_blank_controls(
        frame,catalog_matcher,old,fwhm,protected_mask)
    blank_abs=[r["median_abs_false_wing_nMgy_per_pixel"]
               for r in blanks["per_ring"]]
    training=[by_rank[r] for r in BRIGHT_TRAIN_RANKS]
    heldout=[by_rank[r] for r in BRIGHT_RANK_HOLDOUT]
    hypotheses=(
        "sky_b_subtracted_median_over_peak",
        "sky_a_subtracted_median_over_peak",
        "conditional_plane_subtracted_median_over_peak",
    )
    variants={key:_joint_nuisance_fit(training,blank_abs,key)
              for key in hypotheses}
    primary=variants[hypotheses[0]]
    heldout_results=_predictions_against_stars(primary,heldout)
    # Spatial LOO on *the same fixed bright cohort*, with no posterior
    # exclusion of a poorly fit star.
    spatial=[]
    for q in range(4):
        train=[by_rank[r] for r in BRIGHT_TRAIN_RANKS
               if _quadrant(old[r-1]["native_xy_pix"],frame.shape)!=q]
        test=[by_rank[r] for r in range(1,17)
              if _quadrant(old[r-1]["native_xy_pix"],frame.shape)==q]
        fit=_joint_nuisance_fit(
            train,blank_abs,"sky_b_subtracted_median_over_peak")
        spatial.append({
            "quadrant":q,
            "n_bright_training_outside_quadrant":len(train),
            "n_bright_stars_inside_quadrant":len(test),
            "n_bright_identity_holdouts_inside_quadrant":
                sum(s["catalog_rank"] in BRIGHT_RANK_HOLDOUT for s in test),
            "fitted_conditional_profile":fit.get(
                "radial_common_wing_candidate_over_native_peak"),
            "prediction":_predictions_against_stars(fit,test),
            "spatial_holdout_retrospective_within_one_coadd":True,
        })
    by_ring=[]
    for idx,ring in enumerate(P16D_PROFILE_RINGS_FWHM):
        estimates=[
            variants[v].get(
                "radial_common_wing_candidate_over_native_peak")
            for v in hypotheses
        ]
        values=[v[idx] for v in estimates if v is not None]
        by_ring.append({
            "radius_in_measured_fwhm":list(ring),
            "joint_sky_b_candidate_over_peak":(
                values[0] if estimates[0] is not None else None),
            "joint_hypotheses_range_over_peak":(
                [min(values),max(values)] if values else None),
            "all_three_sky_hypotheses_measured":len(values)==3,
            "blank_median_abs_false_wing_nMgy_per_pixel":blank_abs[idx],
            "blank_max_abs_false_wing_nMgy_per_pixel":blanks[
                "per_ring"][idx]["max_abs_false_wing_nMgy_per_pixel"],
            "sky_and_blank_controls_are_not_calibrated_sigma":True,
        })
    return dict(
        out,
        status="SAME_BRICK_WING_SKY_AND_SPATIAL_FEASIBILITY_ONLY",
        native_catalog_sha256=p17_qa["native_catalog_sha256"],
        n_frozen_stars=48,
        n_bright_joint_fit_train=12,
        n_bright_fixed_identity_holdout=4,
        n_new_faint_identity_holdout_not_promoted=8,
        blank_sky_controls=blanks,
        joint_profile_hypotheses=variants,
        heldout_bright_star_prediction=heldout_results,
        leave_one_brick_quadrant_out=spatial,
        radial_nuisance_and_blank_comparison=by_ring,
        blank_sky_independent_of_frozen_star_profiles_not_exposures=True,
        absolute_outer_wing_level_unidentified_without_sky_prior=True,
        scientific_wing_model_not_adopted=True,
    )
