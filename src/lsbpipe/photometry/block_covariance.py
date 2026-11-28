"""Conditional spatial block covariance, including fitted sky/halo influence.

The block model is a sensitivity hypothesis, not a calibrated detector-noise
model. All bands use the same explicitly ordered tile identities and coordinate
system; no missing tile measurements are synthesized here.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial.distance import cdist

from .uncertainty_budget import _covariance, STATUS


def correction_influence(ring_counts, sky_design_sums, halo_template_sums,
                         sky_influence, halo_influence, *, training_tiles):
    """Linear response of signed ring flux to the mean of every block.

    Influence arrays are derivatives of fitted sky coefficients and halo
    amplitude with respect to all block means, including sky/halo coupling.
    Observed ring counts and training tiles must have disjoint footprints.
    Inputs use the same ordered tile IDs; sky design sums and halo sums
    are evaluated on exactly the valid pixels used in ring_counts.
    """
    c,x,h,b,a=map(lambda v:np.asarray(v,dtype=float),
                (ring_counts,sky_design_sums,halo_template_sums,
                 sky_influence,halo_influence))
    if c.ndim!=2 or min(c.shape)<1:
        raise ValueError("Expected nonempty ring by tile counts")
    nr,nt=c.shape
    if (x.ndim!=2 or x.shape[0]!=nr or b.shape!=(x.shape[1],nt)
            or h.shape!=(nr,) or a.shape!=(nt,)):
        raise ValueError("Inconsistent correction influence dimensions")
    if (not all(np.isfinite(v).all() for v in (c,x,h,b,a))
            or np.any(c<0) or np.any(c!=np.floor(c)) or np.any(h<0)):
        raise ValueError("Invalid count, design or influence arrays")
    t=np.asarray(training_tiles)
    if t.shape!=(nt,) or t.dtype.kind!="b":
        raise ValueError("Explicit Boolean training tile selection required")
    if np.any(c[:,t]!=0) or np.any(b[:,~t]!=0) or np.any(a[~t]!=0):
        raise ValueError("Training and observed footprints must be disjoint")
    return c-x@b-np.outer(h,a)


def within_tile_diagonal_remainder(ring_tile_variance):
    """Conditional IVAR component orthogonal to each tile sum.

    For independent pixel errors with variances summed over disjoint rings,
    subtract the covariance explained by the total of those same pixels
    in each tile. This is NOT the covariance of an arbitrary block estimator.
    """
    v=np.asarray(ring_tile_variance,dtype=float)
    if v.ndim!=2 or min(v.shape)<1 or not np.isfinite(v).all() or np.any(v<0):
        raise ValueError("Expected finite nonnegative ring by tile IVAR variances")
    total=v.sum(axis=0)
    use=total>0
    w=v[:,use]/np.sqrt(total[use])
    # Pairwise Laplacian form avoids subtracting almost equal diagonal
    # terms, including the exact single-ring/zero-remainder case.
    d=-(w@w.T)
    np.fill_diagonal(d,0.)
    np.fill_diagonal(d,-d.sum(axis=1))
    return (d+d.T)/2


def conditional_block_covariance(tile_xy_arcsec, operators, band_scales,
                                 band_correlation, diagonal_remainders, *,
                                 band_order, length_arcsec, mix, max_tiles=4096):
    """C_ij = scale_i scale_j rho_ij L_i K L_j^T plus within-tile terms.

    K = mix exp(-distance/length) + (1-mix) I. The band correlation
    and robust block-mean scales are explicit caller inputs; no inference
    of spatial stationarity or calibration is claimed.
    """
    xy=np.asarray(tile_xy_arcsec,dtype=float)
    bands=tuple(band_order)
    scales=np.asarray(band_scales,dtype=float)
    length,m=float(length_arcsec),float(mix)
    if (xy.ndim!=2 or xy.shape[1]!=2 or not len(xy) or
            not isinstance(max_tiles,int) or len(xy)>max_tiles
            or not np.isfinite(xy).all() or len(np.unique(xy,axis=0))!=len(xy)):
        raise ValueError("Unique finite tile coordinates and explicit memory budget required")
    if (not bands or len(set(bands))!=len(bands) or
            any(not isinstance(b,str) or not b for b in bands)
            or scales.shape!=(len(bands),) or not np.isfinite(scales).all()
            or np.any(scales<=0) or not np.isfinite(length) or length<=0
            or not np.isfinite(m) or not 0<=m<=1):
        raise ValueError("Invalid band order, block scales or kernel hypothesis")
    rho=_covariance(band_correlation,len(bands))
    if not np.allclose(np.diag(rho),1.,rtol=0,atol=1e-12):
        raise ValueError("Band correlation diagonal must equal one")
    ops=tuple(np.asarray(v,dtype=float) for v in operators)
    rem=tuple(diagonal_remainders)
    if len(ops)!=len(bands) or len(rem)!=len(bands):
        raise ValueError("Each ordered band requires an operator and remainder")
    if any(v.ndim!=2 or not v.shape[0] or v.shape[1]!=len(xy)
           or not np.isfinite(v).all() for v in ops):
        raise ValueError("Incompatible tile operator")
    rem=tuple(_covariance(v,o.shape[0]) for v,o in zip(rem,ops))
    k=cdist(xy,xy)
    k/=(-length)
    np.exp(k,out=k)
    k*=m
    k.flat[::len(xy)+1]+=1-m
    rows=[]
    projected=[o@k for o in ops]
    for i,oi in enumerate(ops):
        row=[]
        for j,oj in enumerate(ops):
            block=scales[i]*scales[j]*rho[i,j]*(projected[i]@oj.T)
            if i==j:
                block+=rem[i]
            row.append(block)
        rows.append(row)
    c=np.block(rows)
    c=_covariance((c+c.T)/2,len(c))
    return {"status":STATUS,"covariance":c,"band_order":bands,
            "ring_counts_by_band":tuple(o.shape[0] for o in ops),
            "length_arcsec":length,"mix":m,"n_tiles":len(xy),
            "covariance_calibrated":False}


def paired_control_statistics(residuals, *, band_order, tile_ids, minimum_tiles, clip_sigma):
    """Robust marginal scales and PSD sample correlation on matched far tiles.

    Rows must represent the same disjoint control footprints in every band.
    Clip and minimum count are explicit diagnostic choices, not calibration.
    No assumption of independent samples or effective degrees of freedom.
    """
    r=np.asarray(residuals,dtype=float)
    bands=tuple(band_order)
    ids=tuple(tile_ids)
    clip=float(clip_sigma)
    if (r.ndim!=2 or r.shape[1]!=len(bands) or not bands or len(set(bands))!=len(bands)
            or len(ids)!=len(r) or len(set(ids))!=len(ids) or not np.isfinite(r).all()
            or not isinstance(minimum_tiles,int) or minimum_tiles<max(3,len(bands)+1)
            or len(r)<minimum_tiles or not np.isfinite(clip) or clip<=0):
        raise ValueError("Insufficient matched finite control tiles or invalid configuration")
    center=np.median(r,axis=0)
    scale=1.4826*np.median(np.abs(r-center),axis=0)
    if np.any(scale<=0) or not np.isfinite(scale).all():
        raise ValueError("Degenerate robust control scale")
    z=np.clip((r-center)/scale,-clip,clip)
    z-=z.mean(axis=0)
    norm=np.linalg.norm(z,axis=0)
    if np.any(norm<=0) or not np.isfinite(norm).all():
        raise ValueError("Degenerate matched control correlation")
    unit=z/norm
    rho=_covariance(unit.T@unit,len(bands))
    return {"status":STATUS,"band_order":bands,"n_matched_tiles":len(r),
            "scales":scale,"correlation":rho,"clip_sigma":clip,
            "covariance_calibrated":False}


def split_control_block_noise(scales, correlation, formal_mean_variances):
    """Separate an explicit independent-IVAR component from control scatter.

    A conditional additive model: total control covariance = coherent
    block covariance + independent formal mean variance. Subtract the
    median formal variance per band, and reject a nonpositive/indefinite
    result rather than clipping eigenvalues or silently using a zero term.
    The empirical estimates remain sensitivity hypotheses.
    """
    s=np.asarray(scales,dtype=float)
    v=np.asarray(formal_mean_variances,dtype=float)
    if (s.ndim!=1 or not len(s) or not np.isfinite(s).all() or np.any(s<=0)
            or v.ndim!=2 or v.shape[1]!=len(s) or not len(v)
            or not np.isfinite(v).all() or np.any(v<=0)):
        raise ValueError("Positive finite control scales and mean variances required")
    rho=_covariance(correlation,len(s))
    if not np.allclose(np.diag(rho),1.,rtol=0,atol=1e-12):
        raise ValueError("Expected a correlation matrix")
    white=np.median(v,axis=0)
    coherent=_covariance(rho*np.outer(s,s)-np.diag(white),len(s))
    if np.any(np.diag(coherent)<=0):
        raise ValueError("Coherent block variance not resolved above formal IVAR")
    sig=np.sqrt(np.diag(coherent))
    return {"status":STATUS,"coherent_scales":sig,
            "coherent_correlation":coherent/np.outer(sig,sig),
            "formal_mean_variance_median":white,"covariance_calibrated":False}


def independent_ivar_covariance(ring_variances, correction_operator, training_mean_variances):
    """Exact linear propagation conditional on independent native pixels.

    Ring sums and fitted training tile means have disjoint support. The
    correction operator contains sky and halo influence, including their
    shared training data. No cross-band white covariance is assumed.
    """
    v=np.asarray(ring_variances,dtype=float)
    q=np.asarray(correction_operator,dtype=float)
    t=np.asarray(training_mean_variances,dtype=float)
    if (v.ndim!=1 or not len(v) or q.ndim!=2 or q.shape[0]!=len(v)
            or t.shape!=(q.shape[1],) or not all(np.isfinite(a).all() for a in (v,q,t))
            or np.any(v<=0) or np.any(t<0)):
        raise ValueError("Invalid formal ring/training variance arrays")
    d=np.diag(v)+(q*t)@q.T
    return _covariance((d+d.T)/2,len(v))
