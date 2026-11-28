"""Experimental signed annular inference and conditional GLS diagnostics.

All array operands refer to one validated native grid, or to an explicitly
ordered vector of native-band annuli. No missing pixels are reconstructed.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import cholesky, solve_triangular, cho_solve
from scipy.optimize import least_squares
from scipy.signal import fftconvolve

from .uncertainty_budget import _covariance, STATUS


def annular_support(rings, sectors, valid, *, n_rings, n_sectors,
                    ring_min, cumulative_min, sector_min):
    """In-frame coverage only; -1 denotes outside the declared domain.

    Caller must separately demonstrate the full geometry is inside the frame.
    A missing sector fails the gate. Eligibility is contiguous from ring zero.
    """
    r, s, v = map(np.asarray, (rings, sectors, valid))
    if (r.ndim != 2 or r.shape != s.shape or r.shape != v.shape or
            r.dtype.kind not in "iu" or s.dtype.kind not in "iu" or v.dtype.kind != "b"):
        raise ValueError("Expected same-grid integer labels and Boolean validity")
    if not isinstance(n_rings, int) or not isinstance(n_sectors, int) or min(n_rings, n_sectors) < 1:
        raise ValueError("Positive explicit ring/sector counts required")
    thresholds = np.asarray([ring_min, cumulative_min, sector_min], dtype=float)
    if not np.isfinite(thresholds).all() or np.any(thresholds <= 0) or np.any(thresholds > 1):
        raise ValueError("Coverage thresholds must lie in (0,1]")
    if np.any(r < -1) or np.any(r >= n_rings):
        raise ValueError("Invalid ring label")
    inside = r >= 0
    if np.any((s[inside] < 0) | (s[inside] >= n_sectors)):
        raise ValueError("Invalid sector label")
    idx = r[inside].astype(np.int64)*n_sectors+s[inside].astype(np.int64)
    total = np.bincount(idx, minlength=n_rings*n_sectors).reshape(n_rings,n_sectors)
    kept = np.bincount(idx[v[inside]], minlength=n_rings*n_sectors).reshape(n_rings,n_sectors)
    divide = lambda a,b: np.divide(a,b,out=np.zeros_like(a,dtype=float),where=b>0)
    sector = divide(kept,total)
    ring = divide(kept.sum(axis=1),total.sum(axis=1))
    cumulative = divide(kept.sum(axis=1).cumsum(),total.sum(axis=1).cumsum())
    gate = (ring >= ring_min) & (cumulative >= cumulative_min) & (sector.min(axis=1) >= sector_min)
    return {"status": STATUS, "ring_fraction": ring, "cumulative_fraction": cumulative,
            "sector_fraction": sector, "eligible": np.logical_and.accumulate(gate),
            "coverage_scope": "IN_FRAME_ONLY"}


def adjoint_annular_operator(radius, radial_nodes, ring_masks, psf, *,
                             grid_id, mask_grid_id):
    """Project a linearly interpolated radial intensity through a finite PSF.

    Pixel-center rings, exact supplied Boolean support; zero exterior boundary
    for the finite image. All radii in the same units. Flux is nMgy if input
    profile is nMgy/pixel. No interpolation over missing data.
    """
    if not isinstance(grid_id,str) or not grid_id or grid_id != mask_grid_id:
        raise ValueError("Grid identity mismatch")
    r, nodes, k = map(lambda x: np.asarray(x,dtype=float), (radius,radial_nodes,psf))
    masks = np.asarray(ring_masks)
    if (r.ndim != 2 or not r.size or nodes.ndim != 1 or len(nodes)<2 or
            masks.ndim!=3 or masks.shape[1:]!=r.shape or masks.dtype.kind!="b"):
        raise ValueError("Incompatible native ring masks or radial grid")
    if not all(np.isfinite(x).all() for x in (r,nodes,k)) or np.any(np.diff(nodes)<=0):
        raise ValueError("Nonfinite or nonmonotonic grid")
    if nodes[0]<0 or r.min()<nodes[0] or r.max()>nodes[-1]:
        raise ValueError("Radial grid must cover every native pixel")
    if (k.ndim!=2 or not k.size or any(n%2!=1 for n in k.shape) or
            np.any(k<0) or not np.isclose(k.sum(),1,rtol=1e-6,atol=1e-9)):
        raise ValueError("Odd finite normalized nonnegative PSF required")
    if not len(masks):
        raise ValueError("Nonempty disjoint ring masks required")
    used=np.zeros(r.shape,dtype=bool)
    result=[]
    py,px=np.array(k.shape)//2
    for mask in masks:
        if not mask.any() or np.any(mask&used):
            raise ValueError("Nonempty disjoint ring masks required")
        used |= mask
        iy,ix=np.nonzero(mask)
        sl=np.s_[max(0,iy.min()-py):min(r.shape[0],iy.max()+py+1),
                 max(0,ix.min()-px):min(r.shape[1],ix.max()+px+1)]
        rr=r[sl]
        left=np.clip(np.searchsorted(nodes,rr,side="right")-1,0,len(nodes)-2)
        f=(rr-nodes[left])/(nodes[left+1]-nodes[left])
        # The finite adjoint support is the observed mask dilated by the
        # kernel half-size. Cropping here preserves the full zero-boundary sum.
        weight=fftconvolve(mask[sl].astype(float),k[::-1,::-1],mode="same")
        a=np.bincount(left.ravel(),weights=(weight*(1-f)).ravel(),minlength=len(nodes))
        a+=np.bincount((left+1).ravel(),weights=(weight*f).ravel(),minlength=len(nodes))
        result.append(a)
    return np.asarray(result)


def propagate_linear_covariance(operator, input_covariance, independent_covariance):
    """L C L^T + D, conditional on fixed operators and no cross term with D."""
    l = np.asarray(operator,dtype=float)
    if l.ndim!=2 or min(l.shape)<1 or not np.isfinite(l).all():
        raise ValueError("Invalid linear operator")
    c = _covariance(input_covariance,l.shape[1])
    d = _covariance(independent_covariance,l.shape[0])
    return _covariance(l@c@l.T+d,l.shape[0])


def _indices(mask, n):
    m=np.asarray(mask)
    if m.shape!=(n,) or m.dtype.kind!="b" or not m.any():
        raise ValueError("Explicit nonempty Boolean selection required")
    return np.flatnonzero(m)


def fit_conditional_gls(observed, covariance, predictor, initial, bounds, *, train):
    """Bounded GLS; no automatic rescaling by reduced chi-square.

    Predictor returns full band-ordered signed ring fluxes. Conditional
    parameter covariance is returned only for full-rank interior solutions.
    """
    d=np.asarray(observed,dtype=float)
    p=np.asarray(initial,dtype=float)
    if d.ndim!=1 or p.ndim!=1 or not np.isfinite(d).all() or not np.isfinite(p).all():
        raise ValueError("Invalid flux or initial parameter vector")
    ix=_indices(train,len(d))
    if len(ix)<=len(p):
        raise ValueError("Insufficient degrees of freedom")
    c=_covariance(covariance,len(d))
    lower=cholesky(c[np.ix_(ix,ix)],lower=True)
    def evaluate(z):
        a=np.asarray(predictor(z),dtype=float)
        if a.shape!=d.shape or not np.isfinite(a).all():
            raise ValueError("Predictor returned invalid flux vector")
        return a
    fit=least_squares(lambda z:solve_triangular(lower,(evaluate(z)-d)[ix],lower=True),
                      p,bounds=bounds,x_scale="jac")
    rank=int(np.linalg.matrix_rank(fit.jac))
    boundary=bool(np.any(fit.active_mask))
    cov=None
    if fit.success and rank==len(p) and not boundary:
        _, singular, vt=np.linalg.svd(fit.jac,full_matrices=False)
        cov=(vt.T/singular**2)@vt
    return {"status":STATUS,"success":bool(fit.success),"boundary":boundary,
            "rank":rank,"parameters":fit.x,"prediction":evaluate(fit.x),
            "chi2":float(fit.fun@fit.fun),"dof":len(ix)-len(p),
            "conditional_parameter_covariance":cov,"covariance_calibrated":False}


def conditional_holdout(observed, prediction, covariance, jacobian,
                         parameter_covariance, *, train, test):
    """Correlated predictive residual with first-order parameter uncertainty.

    J is the derivative of full predicted flux w.r.t. the fitted parameters.
    C_theta must refer to an interior training-only fit. No independence
    assumption between training and test noise; no calibrated p-value.
    """
    d,m,j=map(lambda x:np.asarray(x,dtype=float),(observed,prediction,jacobian))
    if (d.ndim!=1 or m.shape!=d.shape or j.ndim!=2 or j.shape[0]!=len(d) or
            not all(np.isfinite(x).all() for x in (d,m,j))):
        raise ValueError("Invalid observed/model/Jacobian dimensions")
    f,t=_indices(train,len(d)),_indices(test,len(d))
    if np.intersect1d(f,t).size:
        raise ValueError("Training and test overlap")
    c=_covariance(covariance,len(d))
    cp=_covariance(parameter_covariance,j.shape[1])
    ff=cholesky(c[np.ix_(f,f)],lower=True)
    tf=c[np.ix_(t,f)]
    a=cho_solve((ff,True),tf.T).T
    residual=(d-m)[t]-a@(d-m)[f]
    jj=j[t]-a@j[f]
    cc=c[np.ix_(t,t)]-a@tf.T+jj@cp@jj.T
    cc=_covariance((cc+cc.T)/2,len(t))
    w=solve_triangular(cholesky(cc,lower=True),residual,lower=True)
    return {"status":STATUS,"residual":residual,"predictive_covariance":cc,
            "chi2":float(w@w),"n_test":len(t),"covariance_calibrated":False}


def control_model_comparison(observed, nominal_model_flux, control_residuals,
                             covariance, models, *, null_model, train,
                             reference_support_id, control_support_ids,
                             control_ids, overlap, control_provenance):
    """Ring-domain control injections with predeclared comparable support.

    Each model maps to (predictor, initial, bounds). Fit the same full set
    for observations and nominal model plus each *measured* control vector.
    This routine does not choose positions or certify exchangeability,
    identical detector noise, or a calibrated empirical false-positive rate.
    Overlap describes the ENTIRE control footprint including sky training.
    """
    obs=np.asarray(observed,dtype=float)
    base=np.asarray(nominal_model_flux,dtype=float)
    controls=np.asarray(control_residuals,dtype=float)
    if (obs.ndim!=1 or base.shape!=obs.shape or controls.ndim!=2
            or controls.shape[1:]!=obs.shape or not len(controls)
            or not all(np.isfinite(x).all() for x in (obs,base,controls))):
        raise ValueError("Finite observed and matched control vectors required")
    count=len(controls)
    ids=tuple(control_ids)
    supports=tuple(control_support_ids)
    if (not isinstance(reference_support_id,str) or not reference_support_id
            or len(ids)!=count or len(set(ids))!=count or any(not isinstance(i,str) or not i for i in ids)
            or len(supports)!=count or any(i!=reference_support_id for i in supports)
            or not isinstance(control_provenance,str) or not control_provenance.strip()):
        raise ValueError("Unique controls with exact declared support and provenance required")
    ov=np.asarray(overlap)
    if ov.shape!=(count,count) or ov.dtype.kind!="b" or not np.array_equal(ov,ov.T):
        raise ValueError("Explicit symmetric footprint overlap matrix required")
    if np.any(ov & ~np.eye(count,dtype=bool)):
        raise ValueError("Controls overlap; do not count them as separate controls")
    if null_model not in models or len(models)<2:
        raise ValueError("Explicit null and alternative model set required")
    def score(values):
        results={}
        for name,(predict,initial,bounds) in models.items():
            fit=fit_conditional_gls(values,covariance,predict,initial,bounds,train=train)
            if not fit["success"]:
                raise ValueError("Model comparison contains a failed fit")
            results[name]={"chi2":fit["chi2"],"boundary":fit["boundary"],"rank":fit["rank"]}
        null=results[null_model]["chi2"]
        alternative=min((name for name in results if name!=null_model),
                        key=lambda name:results[name]["chi2"])
        return {"delta_chi2":max(0.,null-results[alternative]["chi2"]),
                "best_alternative":alternative,"fits":results}
    target=score(obs)
    rows=[dict(control_id=cid,**score(base+v)) for cid,v in zip(ids,controls)]
    exceed=sum(row["delta_chi2"]>=target["delta_chi2"] for row in rows)
    return {"status":STATUS,"target":target,"controls":rows,
            "n_controls":count,"n_equal_or_exceed":exceed,
            "control_provenance":control_provenance,"reference_support_id":reference_support_id,
            "calibrated_false_positive_rate":None,"exchangeability_verified":False,
            "kind":"RING_DOMAIN_CONTROL_STRESS_TEST"}
