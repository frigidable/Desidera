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
