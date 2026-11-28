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
    d=np.diag(v.sum(axis=1))-w@w.T
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
