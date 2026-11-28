"""Conditional asymptotic flux decomposition on exact annular support.

Every estimate remains diagnostic. The hybrid total assumes a Sérsic
correction for all unobserved light and PSF redistribution; it is not an
empirical asymptote or a bound on arbitrary diffuse envelopes.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import cholesky, solve_triangular

from .native_annular_diagnostic import annular_record_sha256, numeric_sha256
from .radial_models import (EllipseNode, NestedEllipse, sersic_intensity,
                            sersic_model_total, conditional_plateau)
from .uncertainty_budget import _covariance, STATUS

COMPONENTS = (
    "observed_valid_flux", "model_missing_in_frame_aperture",
    "hybrid_in_frame_aperture", "intrinsic_model_enclosed",
    "intrinsic_model_outer_tail", "model_psf_and_sampling_transfer",
    "pure_model_total", "hybrid_model_corrected_total",
)


def _derivative(fn, z, step=1e-5):
    return np.column_stack([(fn(z+v*step)-fn(z-v*step))/(2*step)
                            for v in np.eye(len(z))])


def linearized_parameter_response(jacobian, covariance, train):
    """d fitted parameters / d full measured flux, at fixed covariance.

    Uses the Gauss-Newton approximation of the interior training fit.
    Retains correlations with measurements outside training through C.
    """
    j=np.asarray(jacobian,dtype=float)
    t=np.asarray(train)
    if (j.ndim!=2 or not j.shape[1] or not np.isfinite(j).all()
            or t.shape!=(len(j),) or t.dtype.kind!="b"
            or t.sum()<=j.shape[1]):
        raise ValueError("Invalid training Jacobian")
    c=_covariance(covariance,len(j))
    ix=np.flatnonzero(t)
    chol=cholesky(c[np.ix_(ix,ix)],lower=True)
    w=solve_triangular(chol,j[ix],lower=True)
    u,s,vt=np.linalg.svd(w,full_matrices=False)
    if np.linalg.matrix_rank(w)!=j.shape[1]:
        raise ValueError("Training Jacobian is rank deficient")
    pinv=(vt.T/s)@u.T
    response=np.zeros((j.shape[1],len(j)))
    response[:,ix]=solve_triangular(chol.T,pinv.T,lower=False).T
    return response


def asymptotic_flux_diagnostic(records, covariance, fit, *, aperture_count,
                               plateau_inner_count=None,
                               fractional_tolerance=None, normal_multiplier=None):
    """Separate measured, model-corrected and extrapolated signed flux.

    Aperture_count is a predeclared number of inner rings. The full
    in-frame aperture operator is required, even for masked pixels.
    Conditional propagation includes the covariance between measured
    flux and the parameters fitted to it. Geometry, masks, source
    corrections, PSF and C stay fixed. No extra independent parameter
    error is added to the measured flux uncertainty.
    """
    records=tuple(records)
    if (not records or fit.get("status")!=STATUS or not fit.get("success")
            or fit.get("boundary") or fit.get("conditional_parameter_covariance") is None):
        raise ValueError("Successful interior diagnostic fit required")
    if (tuple(annular_record_sha256(r) for r in records)!=tuple(fit["record_sha256"])
            or numeric_sha256(covariance)!=fit["covariance_sha256"]):
        raise ValueError("Fit, covariance and exact records must match")
    bands=tuple(fit["band_order"])
    n=len(records[0]["ring_edges_arcsec"])-1
    if (len(bands)!=len(records) or not isinstance(aperture_count,int)
            or isinstance(aperture_count,bool) or not 1<=aperture_count<=n):
        raise ValueError("Explicit aperture ring count required")
    for b,r in zip(bands,records):
        if (r["dataset_id"].split("/")[-1]!=b
                or len(r["ring_edges_arcsec"])!=n+1
                or not np.array_equal(r["ring_edges_arcsec"],records[0]["ring_edges_arcsec"])
                or not np.all(r["coverage"]["eligible"][:aperture_count])
                or not r.get("full_aperture_finite_psf_support",False)
                or "full_aperture_operator" not in r):
            raise ValueError("Eligible common aperture and complete finite PSF footprint required")
        op=np.asarray(r["full_aperture_operator"])
        if (op.shape!=r["operator"].shape or op.shape[0]!=n
                or not np.isfinite(op).all()):
            raise ValueError("Full aperture operator shape mismatch")
    c=_covariance(covariance,len(bands)*n)
    d=np.concatenate([r["observed_flux_nMgy"] for r in records])
    z=np.asarray(fit["parameters"],dtype=float)
    if z.shape!=(3*len(bands),) or not np.isfinite(z).all() or not np.isfinite(d).all():
        raise ValueError("Finite signed fluxes and log Sérsic parameters required")
    geo=records[0]["geometry"]
    g=NestedEllipse(EllipseNode(**geo["inner"]),EllipseNode(**geo["outer"]),
                    geo["transition_start"],geo["transition_end"])
    a=float(records[0]["ring_edges_arcsec"][aperture_count])
    def predict(p):
        return np.concatenate([r["operator"]@sersic_intensity(
            r["radial_nodes_arcsec"],*np.exp(p[3*i:3*i+3]))
            for i,r in enumerate(records)])
    m=predict(z)
    if not np.allclose(m,fit["prediction"],rtol=1e-8,atol=1e-9):
        raise ValueError("Stored model prediction does not match fitted parameters")
    tr=np.asarray(fit["train_rings"])
    if tr.shape!=(n,) or tr.dtype.kind!="b":
        raise ValueError("Invalid training ring selection")
    j=_derivative(predict,z)
    response=linearized_parameter_response(j,c,np.tile(tr,len(bands)))
    def model_terms(p):
        rows=[]
        for i,r in enumerate(records):
            amp,re,idx=np.exp(p[3*i:3*i+3])
            profile=sersic_intensity(r["radial_nodes_arcsec"],amp,re,idx)
            valid=float(np.sum(r["operator"][:aperture_count]@profile))
            full=float(np.sum(r["full_aperture_operator"][:aperture_count]@profile))
            integ=sersic_model_total(g,amp/r["pixel_area_arcsec2"],re,idx,
                                     aperture_radius=a)
            inside,tail,total=(integ[k] for k in
                ("intrinsic_enclosed_flux","intrinsic_outer_tail","model_total"))
            missing=full-valid
            if missing < -1e-8*max(1.,full):
                raise ValueError("Full model aperture has less light than valid support")
            rows.extend([0.,missing,missing,inside,tail,inside-full,total,total-valid])
        return np.asarray(rows)
    terms=model_terms(z)
    deriv=_derivative(model_terms,z)
    direct=np.zeros((len(terms),len(d)))
    for i in range(len(bands)):
        for k in (0,2,7):
            direct[i*len(COMPONENTS)+k,i*n:i*n+aperture_count]=1
    values=direct@d+terms
    influence=direct+deriv@response
    cov=influence@c@influence.T
    cov=_covariance((cov+cov.T)/2,len(values))
    out={}
    plateau_requested=any(x is not None for x in
                          (plateau_inner_count,fractional_tolerance,normal_multiplier))
    if plateau_requested and any(x is None for x in
                          (plateau_inner_count,fractional_tolerance,normal_multiplier)):
        raise ValueError("Declare all plateau settings together")
    for i,b in enumerate(bands):
        ix=i*len(COMPONENTS)
        row={}
        for k,name in enumerate(COMPONENTS):
            flux=float(values[ix+k])
            sigma=float(np.sqrt(max(0.,cov[ix+k,ix+k])))
            row[name]={"flux_nMgy":flux,"conditional_sigma_nMgy":sigma}
            if k in (0,2,6,7) and flux>0:
                row[name].update({"AB_mag_before_extinction":float(22.5-2.5*np.log10(flux)),
                                  "conditional_sigma_mag_linearized":float(2.5/np.log(10)*sigma/flux)})
        row["intrinsic_tail_fraction_of_model"]=float(values[ix+4]/values[ix+6])
        if plateau_requested:
            sl=slice(i*n,(i+1)*n)
            row["observed_plateau"]=conditional_plateau(d[sl],c[sl,sl],
                records[i]["coverage"]["eligible"],inner_count=plateau_inner_count,
                outer_count=aperture_count,fractional_tolerance=fractional_tolerance,
                normal_multiplier=normal_multiplier)
        out[b]=row
    return {"status":STATUS,"aperture_radius_arcsec":a,"aperture_count":aperture_count,
            "band_order":bands,"component_order":COMPONENTS,"bands":out,
            "conditional_component_covariance_nMgy2":cov,
            "covariance_order":"BAND_MAJOR_COMPONENT_MINOR",
            "scope":"IN_FRAME_APERTURE_AND_UNTRUNCATED_MODEL_ONLY",
            "uncertainty_kind":"FIXED_GEOMETRY_GAUSS_NEWTON_CONDITIONAL",
            "record_sha256":fit["record_sha256"],"covariance_sha256":fit["covariance_sha256"],
            "full_geometry_containment_verified":False,
            "covariance_calibrated":False,"asymptotic_flux_adopted":False}
