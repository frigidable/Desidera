"""P15 exploratory native-pixel PSF-convolved Tractor galaxy scene model.

Fixed original coadd pixels and already-fitted experimental PSF. REX,
EXP and DEV shapes use actual same-brick Tractor catalog metadata,
not fabricated radii. The analytic oversampled Sersic + FFT model is
NOT exactly the Tractor Gaussian-mixture, pixelization or PSF renderer.
Never use these diagnostics for PSF adoption or calibrated photometry.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.signal import fftconvolve

from lsbpipe.photometry.psf_joint_scene import (
    _projected_nnls, _shifted_native_kernel,
)

SCHEMA = "desidera.psf0-native-tractor-extended-scene/1"
SERSIC_B = {1.: 1.6783469900166612, 4.: 7.669249442500804}
PROFILE_N = {"REX": 1., "EXP": 1., "DEV": 4.}
SUBPIXEL_GRID = 3
MAX_NUISANCE_SOURCES = 28  # compute cap: fail whole fit, never truncate


def finite(value):
    try:
        number = float(value)
    except (ValueError, TypeError, OverflowError):
        return None
    return number if math.isfinite(number) else None


def audit_tractor_shape(row):
    kind = str(row.get("tractor_type", row.get("type", ""))).strip().upper()
    radius = finite(row.get("shape_r_arcsec"))
    e1, e2 = finite(row.get("shape_e1")), finite(row.get("shape_e2"))
    out = {
        "type": kind, "shape_r_arcsec": radius,
        "shape_e1": e1, "shape_e2": e2,
        "catalog_sersic_n": finite(row.get("sersic_n")),
        "flux_not_used_as_amplitude_prior": True,
        "science_adopted": False,
    }
    if kind not in PROFILE_N:
        out["status"] = "UNSUPPORTED_TRACTOR_PROFILE_TYPE"
        return out
    if radius is None or radius <= 0:
        out["status"] = "MISSING_OR_INVALID_HALF_LIGHT_RADIUS_ARCSEC"
        return out
    if kind == "REX":
        # REX is round exponential by catalog model definition.
        q, pa = 1., None
    else:
        if e1 is None or e2 is None:
            out["status"] = "MISSING_TRACTOR_ELLIPTICITY"
            return out
        e = math.hypot(e1, e2)
        if e >= 1.:
            out["status"] = "INVALID_TRACTOR_ELLIPTICITY"
            return out
        q = (1.-e)/(1.+e)
        # DR11 convention PA = (180 - phi) % 180,
        # phi = 0.5 atan2(e2,e1). PA from North to East.
        pa = (180.-.5*math.degrees(math.atan2(e2,e1)))%180.
    out.update({
        "status": "VALID_APPROXIMATE_EXTENDED_TEMPLATE",
        "sersic_n": PROFILE_N[kind], "axis_ratio": q,
        "position_angle_east_of_north_deg": pa,
        "rex_round_by_definition": kind == "REX",
        "not_exact_tractor_renderer": True,
    })
    return out


def local_tangent_jacobian(wcs, xy):
    """Arcsec east/north per native image x/y pixel, including sky rotation."""
    x, y = map(float, xy)
    ra, dec = wcs.all_pix2world(
        np.array([x,x+1.,x]), np.array([y,y,y+1.]), 0)
    ra, dec = np.asarray(ra,float), np.asarray(dec,float)
    if not np.isfinite(ra).all() or not np.isfinite(dec).all():
        raise ValueError("NONFINITE_WCS_PROJECTION")
    dr = ((ra[1:]-ra[0]+180.)%360.-180.)
    east = dr*np.cos(np.deg2rad(dec[0]))*3600.
    north = (dec[1:]-dec[0])*3600.
    J = np.array([[east[0],east[1]],[north[0],north[1]]])
    if (not np.isfinite(J).all()
            or abs(np.linalg.det(J))<1e-9
            or np.linalg.cond(J)>1e4):
        raise ValueError("DEGENERATE_LOCAL_WCS")
    return J


def convolved_native_sersic(*, kernel, wcs, source_xy,
                             origin_xy, width, morphology):
    """Finite-domain, subpixel integrated intrinsic Sersic convolved with PSF.

    Total-flux normalization is intentionally NOT claimed: intrinsic
    DEV wings beyond the finite modeled domain are unmeasured.
    """
    if morphology.get("status") != "VALID_APPROXIMATE_EXTENDED_TEMPLATE":
        raise ValueError("NO_VALID_NATIVE_TRACTOR_MORPHOLOGY")
    k = np.asarray(kernel,float)
    if (k.ndim!=2 or k.shape[0]!=k.shape[1] or
            k.shape[0]%2!=1 or type(width) is not int or width<=0
            or not np.isfinite(k).all() or np.min(k)<0
            or not np.isclose(k.sum(),1.,rtol=1e-8,atol=1e-8)):
        raise ValueError("INVALID_FROZEN_NATIVE_PSF")
    n = morphology["sersic_n"]
    q = morphology["axis_ratio"]
    r_eff = morphology["shape_r_arcsec"]
    pa = morphology["position_angle_east_of_north_deg"]
    if n not in SERSIC_B or not 0<q<=1:
        raise ValueError("INVALID_SERSIC_GEOMETRY")
    theta = math.radians(0. if pa is None else pa)
    J = local_tangent_jacobian(wcs,source_xy)
    pad = k.shape[0]//2
    yy,xx = np.indices((width+2*pad,width+2*pad),dtype=float)
    xx += float(origin_xy[0])-pad
    yy += float(origin_xy[1])-pad
    intrinsic = np.zeros_like(xx)
    for iy in range(SUBPIXEL_GRID):
        for ix in range(SUBPIXEL_GRID):
            dx = xx+(ix+.5)/SUBPIXEL_GRID-.5-float(source_xy[0])
            dy = yy+(iy+.5)/SUBPIXEL_GRID-.5-float(source_xy[1])
            east=J[0,0]*dx+J[0,1]*dy
            north=J[1,0]*dx+J[1,1]*dy
            major=east*math.sin(theta)+north*math.cos(theta)
            minor=east*math.cos(theta)-north*math.sin(theta)
            elliptical_R=np.hypot(major,minor/q)
            intrinsic += np.exp(-SERSIC_B[n]*
                                (elliptical_R/r_eff)**(1./n))
    intrinsic /= SUBPIXEL_GRID**2
    if not np.isfinite(intrinsic).all() or intrinsic.max()<=0.:
        raise ValueError("NONFINITE_INTRINSIC_PROFILE")
    img=fftconvolve(intrinsic,k,mode="same")
    out=np.maximum(img[pad:pad+width,pad:pad+width],0.)
    if not np.isfinite(out).all():
        raise ValueError("NONFINITE_CONVOLVED_PROFILE")
    return out



def _predeclared_ring_prediction(z, B, P, weights, radius_pix, *,
                                guessed_fwhm_pix, n_base_columns):
    """Withhold 2–4 detection FWHM annulus from BOTH model fits.

    Fixed, predeclared annular split on ORIGINAL same-pixel coadds.
    Templates/morphology/PSF are frozen before the split. This is
    predictive *pixel* QA only, NOT an independent star validation:
    resampled coadd pixel covariances cross the train/test boundary.
    """
    ring = ((radius_pix > 2.*guessed_fwhm_pix)
            & (radius_pix <= 4.*guessed_fwhm_pix))
    train = ~ring
    out = {
        "schema": "desidera.psf0-p15-native-annular-predictive-check/1",
        "heldout_radius_inner_guessed_fwhm": 2.,
        "heldout_radius_outer_guessed_fwhm": 4.,
        "n_fitting_native_pixels": int(np.count_nonzero(train)),
        "n_reserved_native_pixels": int(np.count_nonzero(ring)),
        "both_fits_exclude_all_reserved_ring_pixels": True,
        "fixed_templates_morphology_and_psf_prior_to_pixel_split": True,
        "predictive_test_uses_same_original_native_pixel_mask": True,
        "pixel_resampling_covariance_across_ring_boundary_unmodeled": True,
        "psf_itself_not_independently_trained_for_every_star": True,
        "no_model_selection_or_science_adoption": True,
    }
    if (not np.any(ring) or
            np.count_nonzero(train) <= P.shape[1]+3):
        return dict(out, status="INSUFFICIENT_RING_OR_TRAINING_PIXELS")
    try:
        fitted, sky, _, condition = _projected_nnls(
            z[train], B[train], P[train], weights[train])
        baseline, sky_base, _, condition_base = _projected_nnls(
            z[train], B[train], P[train, :n_base_columns],
            weights[train])
    except (ValueError, RuntimeError) as exc:
        return dict(
            out, status="RING_EXCLUDED_FITS_NOT_IDENTIFIABLE",
            reason=str(exc))
    prediction = (
        z[ring] - B[ring]@sky - P[ring]@fitted)
    baseline_prediction = (
        z[ring] - B[ring]@sky_base
        - P[ring,:n_base_columns]@baseline)
    score = float(np.dot(
        weights[ring]*prediction, weights[ring]*prediction))
    score_base = float(np.dot(
        weights[ring]*baseline_prediction,
        weights[ring]*baseline_prediction))
    return dict(
        out, status="PREDECLARED_RING_OUT_OF_SAMPLE_PIXEL_PREDICTION",
        weighted_reserved_ring_rss_with_extended=score,
        weighted_reserved_ring_rss_without_extended=score_base,
        signed_predictive_ring_rss_improvement_fraction=(
            float((score_base-score)/score_base)
            if score_base>0 else None),
        training_extended_design_condition_number=float(condition),
        training_point_only_design_condition_number=float(condition_base),
        negative_ring_gain_retained_instead_of_clamped=True,
        not_chi_squared_without_coadd_covariance=True,
    )

def pilot_extended_scene(image, invvar, valid, *, wcs, kernel, target_xy,
                         point_sources, extended_sources,
                         guessed_fwhm_pix, scene_width=None):
    """Fit same-pixel star+PSF-points+sky versus adding typed galaxies."""
    k=np.asarray(kernel,dtype=float)
    width=2*k.shape[0]-1 if scene_width is None else scene_width
    if type(width) is not int or width<k.shape[0] or width%2!=1:
        raise ValueError("INVALID_SCENE_WIDTH")
    image, invvar, valid = (np.asarray(image,float),
                              np.asarray(invvar,float),
                              np.asarray(valid,bool))
    if image.ndim!=2 or image.shape!=invvar.shape or image.shape!=valid.shape:
        raise ValueError("INCONSISTENT_ORIGINAL_IMAGE_ARRAYS")
    cx,cy=map(float,target_xy)
    half=width//2
    x0,y0=int(np.rint(cx))-half,int(np.rint(cy))-half
    base={"schema":SCHEMA,"science_adopted":False,
          "new_psf_reconstruction_performed":False,
          "native_pixels_or_catalog_modified":False}
    if x0<0 or y0<0 or x0+width>image.shape[1] or y0+width>image.shape[0]:
        return dict(base,status="EXPANDED_SCENE_OUTSIDE_NATIVE_FITS")
    sy,sx=slice(y0,y0+width),slice(x0,x0+width)
    z, inv=image[sy,sx],invvar[sy,sx]
    mask=valid[sy,sx]&np.isfinite(z)&np.isfinite(inv)&(inv>0)
    nmask=int(mask.sum())
    yy,xx=np.indices((width,width),dtype=float)
    xx+=x0;yy+=y0
    if not mask[half,half]:
        return dict(base,status="INVALID_NATIVE_CENTRAL_PIXEL")
    templates=[_shifted_native_kernel(k,xx,yy,cx,cy)]
    points=[]
    for source in point_sources:
        src=(float(source["x_pix"]),float(source["y_pix"]))
        t=_shifted_native_kernel(k,xx,yy,*src)
        if np.any(t[mask]>0.):
            templates.append(t)
            points.append(int(source["tractor_objid"]))
    n_point=len(points)
    extended_audit=[]
    admitted=[]
    for item in extended_sources:
        morphology=audit_tractor_shape(item)
        record={
            "tractor_objid":int(item["tractor_objid"]),
            "tractor_type":str(item["tractor_type"]),
            "separation_pix":float(math.hypot(
                item["x_pix"]-cx,item["y_pix"]-cy)),
            "catalog_flux_nMgy":finite(item.get("flux_nMgy")),
            "morphology":morphology}
        if morphology["status"]!="VALID_APPROXIMATE_EXTENDED_TEMPLATE":
            record["status"]="SKIPPED_NO_VERIFIED_MORPHOLOGY"
            extended_audit.append(record)
            continue
        if record["catalog_flux_nMgy"] is None or record["catalog_flux_nMgy"]<=0.:
            record["status"]="SKIPPED_MISSING_OR_NONPOSITIVE_FLUX"
            extended_audit.append(record)
            continue
        try:
            t=convolved_native_sersic(
                kernel=k,wcs=wcs,
                source_xy=(item["x_pix"],item["y_pix"]),
                origin_xy=(x0,y0),width=width,morphology=morphology)
        except (ArithmeticError,ValueError) as exc:
            record["status"]="SKIPPED_INVALID_PROFILE_RENDERING"
            record["reason"]=str(exc)
            extended_audit.append(record)
            continue
        if not np.any(t[mask]>0.):
            record["status"]="SKIPPED_NO_FINITE_TEMPLATE_OVERLAP"
            extended_audit.append(record)
            continue
        record["status"]="APPROXIMATE_TYPED_PROFILE_INCLUDED"
        extended_audit.append(record)
        templates.append(t)
        admitted.append(int(item["tractor_objid"]))
    base.update({
        "n_catalog_extended_candidates":len(extended_sources),
        "n_extended_profiles_rendered":len(admitted),
        "n_extended_catalog_unmodeled":len(extended_sources)-len(admitted),
        "n_point_templates_rendered":n_point,
        "n_valid_original_native_pixels":nmask,
        "extended_source_morphology_audit":extended_audit,
        "all_catalog_extended_profiles_accounted_for":(
            len(admitted)==len(extended_sources)),
        "not_equivalent_to_tractor_gaussian_mixture_renderer":True,
        "finite_intrinsic_and_empirical_psf_support_unvalidated":True,
        "no_coadd_covariance_or_parameter_uncertainty_model":True,
        "profile_amplitudes_not_catalog_total_fluxes":True,
        "subpixel_integration_grid":SUBPIXEL_GRID,
        "source_count_budget":MAX_NUISANCE_SOURCES,
    })
    if n_point+len(admitted)>MAX_NUISANCE_SOURCES:
        return dict(base,status="TOO_MANY_RENDERED_SOURCES_NO_TRUNCATION")
    if not admitted:
        return dict(base,status="NO_VALID_EXTENDED_SOURCE_PROFILE_TO_TEST")
    if nmask<=len(templates)+3:
        return dict(base,status="INSUFFICIENT_NATIVE_PIXEL_DOF")
    if any(np.max(t[mask])<=0 for t in templates):
        return dict(base,status="DEGENERATE_NATIVE_PROFILE_SUPPORT")
    P=np.column_stack([t[mask]/np.max(t[mask]) for t in templates])
    B=np.column_stack((np.ones(nmask),
        ((xx-cx)/max(1.,half))[mask],((yy-cy)/max(1.,half))[mask]))
    w=np.sqrt(inv[mask])
    baseline_columns=1+n_point
    try:
        a,sky,rss,cond=_projected_nnls(z[mask],B,P,w)
        a0,sky0,rss0,cond0=_projected_nnls(
            z[mask],B,P[:,:baseline_columns],w)
    except (RuntimeError, ValueError) as exc:
        return dict(base,status="NUMERICALLY_UNIDENTIFIABLE_EXTENDED_SCENE",
                    reason=str(exc))
    residual=(z[mask]-B@sky-P@a)*w
    residual0=(z[mask]-B@sky0-P[:,:baseline_columns]@a0)*w
    distance=np.hypot(xx-cx,yy-cy)[mask]
    predictive_ring = _predeclared_ring_prediction(
        z[mask], B, P, w, distance,
        guessed_fwhm_pix=guessed_fwhm_pix,
        n_base_columns=baseline_columns)
    regions={}
    for name,radius in [("one_guessed_fwhm",guessed_fwhm_pix),
                        ("four_guessed_fwhm",4*guessed_fwhm_pix)]:
        subset=distance<=radius
        u=float(np.dot(residual0[subset],residual0[subset]))
        v=float(np.dot(residual[subset],residual[subset]))
        regions[name]={
            "n_pixels":int(subset.sum()),
            "baseline_point_star_sky_rss":u,
            "typed_extended_star_point_sky_rss":v,
            "signed_rss_improvement_fraction": (
                float((u-v)/u) if u>0 else None),
            "same_original_native_pixels_globally_fit_models":True,
            "not_a_significance_or_psf_validation":True,
        }
    return dict(base,
        status="PARTIAL_APPROXIMATE_EXTENDED_SCENE_DIAGNOSTIC",
        fitted_extended_source_objids=admitted,
        weighted_rss_without_extended=float(rss0),
        weighted_rss_with_extended=float(rss),
        full_scene_fractional_rss_change=(
            float((rss0-rss)/rss0) if rss0>0 else None),
        central_aperture_comparison=regions,
        predeclared_2_to_4_fwhm_ring_prediction=predictive_ring,
        sky_with_minus_without_extended_at_star_nMgy_pix=
            float(sky[0]-sky0[0]),
        main_star_fit_with_extended_peak_nMgy_pix=float(a[0]),
        main_star_fit_without_extended_peak_nMgy_pix=float(a0[0]),
        full_design_condition_number=float(cond),
        baseline_design_condition_number=float(cond0),
        nested_global_rss_gain_not_evidence_of_better_psf=True,
        omitted_catalog_sources_prevent_complete_deblending=(
            len(admitted)<len(extended_sources)),
    )
