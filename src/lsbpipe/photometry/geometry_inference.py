"""Opt-in PSF-forward geometry candidates on frozen native block support.

Scores use explicitly supplied block-intensity scales, not calibrated
pixel likelihoods. No geometry, correction, or photometry is adopted.
"""
from __future__ import annotations

from dataclasses import asdict

import numpy as np
from scipy.ndimage import distance_transform_cdt
from scipy.optimize import least_squares
from scipy.signal import fftconvolve

from .radial_models import EllipseNode, NestedEllipse, sersic_intensity
from .uncertainty_budget import STATUS


def _node_vector(node):
    return np.array([node.east_arcsec, node.north_arcsec, node.q, node.pa_deg])


def _radius(geometry, east, north):
    """Exact endpoint ellipses; bracket only pixels within the transition.

    NestedEllipse has already certified non-crossing. The 40 bisection
    iterations bound the transition radius error by its width / 2**40.
    """
    def rho(node, e, n):
        x, y = e-node.east_arcsec, n-node.north_arcsec
        p = np.deg2rad(node.pa_deg)
        return np.hypot(x*np.sin(p)+y*np.cos(p),
                        (x*np.cos(p)-y*np.sin(p))/node.q)
    inner = rho(geometry.inner, east, north)
    if geometry.inner == geometry.outer:
        return inner
    outer = rho(geometry.outer, east, north)
    lo, hi = geometry.transition_start, geometry.transition_end
    use = (inner > lo) & (outer < hi)
    a = np.where(inner <= lo, inner, outer)
    if use.any():
        e, n = east[use], north[use]
        low, high = np.full(len(e),lo), np.full(len(e),hi)
        for _ in range(40):
            mid = (low+high)/2
            x,y,q,p = geometry.at(mid)
            dx,dy=e-x,n-y
            outside = np.hypot(dx*np.sin(p)+dy*np.cos(p),
                               (dx*np.cos(p)-dy*np.sin(p))/q)>mid
            low=np.where(outside,mid,low)
            high=np.where(outside,high,mid)
        a[use]=(low+high)/2
    return a


def fit_block_geometry(east_arcsec, north_arcsec, images, kernels, *, band_order,
                       valid, region, train_tiles, holdout_tiles, block_size,
                       min_block_fraction, buffer_tiles, block_scales,
                       seed, transition_start, transition_end, midpoint_bounds,
                       delta_limits, delta_scales, profile_initial,
                       profile_bounds, max_nfev=100, max_work_pixels=500000):
    """Compare fixed, constant and variable common ellipse hypotheses.

    Inputs are already prepared on one native pixel grid. The region,
    validity, block split and noise scales are frozen before optimization.
    Whole disjoint tiles are assigned to training or holdout. Model images
    include masked pixels before PSF convolution, but only observed usable
    pixels contribute to the block means. Data are never reconstructed.

    Geometry coordinates: east, north (arcsec), q, unwrapped PA (degrees).
    Midpoint bounds and maximum endpoint differences must guarantee nesting
    everywhere in the optimization box. Per-band profile parameters are
    positive (central nMgy/pixel, radial scale arcsec, Sersic index).
    """
    e,n=map(lambda x:np.asarray(x,dtype=float),(east_arcsec,north_arcsec))
    if (e.ndim!=2 or not e.size or e.shape!=n.shape
            or not np.isfinite(e).all() or not np.isfinite(n).all()):
        raise ValueError("Finite native tangent-plane coordinates required")
    if (not isinstance(max_work_pixels,int) or max_work_pixels<1
            or e.size>max_work_pixels or not isinstance(max_nfev,int) or max_nfev<1):
        raise ValueError("Explicit geometry work budget exceeded or invalid")
    bands=tuple(band_order)
    if not bands or len(set(bands))!=len(bands):
        raise ValueError("Unique explicit band order required")
    data=tuple(np.asarray(x,dtype=float) for x in images)
    psfs=tuple(np.asarray(x,dtype=float) for x in kernels)
    if (len(data)!=len(bands) or len(psfs)!=len(bands)
            or any(x.shape!=e.shape for x in data)):
        raise ValueError("One same-grid image and PSF per band required")
    for k in psfs:
        if (k.ndim!=2 or not k.size or any(v%2!=1 for v in k.shape)
                or not np.isfinite(k).all() or np.any(k<0)
                or not np.isclose(k.sum(),1.,rtol=1e-6,atol=1e-9)):
            raise ValueError("Odd normalized finite nonnegative PSF required")
    v,reg=np.asarray(valid),np.asarray(region)
    if any(x.shape!=e.shape or x.dtype.kind!="b" for x in (v,reg)):
        raise ValueError("Explicit native Boolean valid and fixed region masks required")
    if (not isinstance(block_size,int) or isinstance(block_size,bool) or block_size<2
            or not isinstance(buffer_tiles,int) or buffer_tiles<0
            or not 0<float(min_block_fraction)<=1):
        raise ValueError("Invalid block configuration")
    b=block_size
    ny,nx=e.shape
    ty,tx=ny//b,nx//b
    tr,te=np.asarray(train_tiles),np.asarray(holdout_tiles)
    if any(x.shape!=(ty,tx) or x.dtype.kind!="b" for x in (tr,te)):
        raise ValueError("Explicit Boolean masks on complete native tiles required")
    if not tr.any() or not te.any() or np.any(tr&te):
        raise ValueError("Nonempty disjoint training and holdout tiles required")
    if buffer_tiles and np.any(distance_transform_cdt(~te,metric="chessboard")[tr]<=buffer_tiles):
        raise ValueError("Training tiles intersect the holdout buffer")
    # Ignore incomplete edge tiles, but preserve the actual pixel support
    # within every retained complete tile.
    domain=v[:ty*b,:tx*b]&reg[:ty*b,:tx*b]
    y,x=np.indices(domain.shape)
    ids=y//b*tx+x//b
    counts=np.bincount(ids[domain],minlength=ty*tx)
    eligible=counts>=float(min_block_fraction)*b*b
    train=tr.ravel()&eligible
    test=te.ravel()&eligible
    used=np.flatnonzero(train|test)
    if not train.any() or not test.any():
        raise ValueError("No eligible training or holdout tiles")
    if int(train.sum())*len(bands)<=8+3*len(bands):
        raise ValueError("Insufficient training blocks for variable geometry")
    lut=np.full(ty*tx,-1,dtype=int)
    lut[used]=np.arange(len(used))
    pick=domain&(lut[ids]>=0)
    py,px=np.nonzero(pick)
    mapped=lut[ids[pick]]
    cnt=counts[used]
    # Exact finite convolution footprint, not an arbitrary data cutout.
    hy=max(k.shape[0]//2 for k in psfs)
    hx=max(k.shape[1]//2 for k in psfs)
    y0,y1=int(py.min())-hy,int(py.max())+hy+1
    x0,x1=int(px.min())-hx,int(px.max())+hx+1
    if min(y0,x0)<0 or y1>ny or x1>nx:
        raise ValueError("Selected blocks lack full finite PSF support at frame edge")
    sl=np.s_[y0:y1,x0:x1]
    ee,nn=e[sl],n[sl]
    iy,ix=py-y0,px-x0
    def mean(image):
        return np.bincount(mapped,weights=image[iy,ix],minlength=len(used))/cnt
    observed=[]
    for image in data:
        if not np.isfinite(image[py,px]).all():
            raise ValueError("Nonfinite observed intensity on frozen support")
        observed.append(mean(image[sl]))
    observed=np.asarray(observed)
    scales=np.asarray(block_scales,dtype=float)
    if scales.shape!=(len(bands),) or not np.isfinite(scales).all() or np.any(scales<=0):
        raise ValueError("Explicit positive block intensity scales required")
    if not isinstance(seed,EllipseNode):
        raise ValueError("Explicit constant geometry seed required")
    lo,hi=np.asarray(midpoint_bounds,dtype=float)
    limit=np.asarray(delta_limits,dtype=float)
    ds=np.asarray(delta_scales,dtype=float)
    base=_node_vector(seed)
    if (any(a.shape!=(4,) for a in (lo,hi,limit,ds))
            or not all(np.isfinite(a).all() for a in (lo,hi,limit,ds))
            or np.any(lo>=hi) or np.any(base<lo) or np.any(base>hi)
            or np.any(limit<=0) or np.any(ds<=0)):
        raise ValueError("Invalid midpoint or endpoint-difference bounds/scales")
    a0,a1=float(transition_start),float(transition_end)
    NestedEllipse(seed,seed,a0,a1)
    qmin=lo[2]-limit[2]/2
    qmax=hi[2]+limit[2]/2
    nesting=1.5/(a1-a0)*(np.hypot(*limit[:2])+a1*(limit[2]+np.deg2rad(limit[3])))
    if qmin<=nesting or qmax>1:
        raise ValueError("Entire variable-geometry box must certify nested ellipses")
    ini=np.asarray(profile_initial,dtype=float)
    pl,ph=np.asarray(profile_bounds,dtype=float)
    if (any(a.shape!=(len(bands),3) for a in (ini,pl,ph))
            or not all(np.isfinite(a).all() for a in (ini,pl,ph))
            or np.any(pl<=0) or np.any(pl>=ph) or np.any(ini<pl) or np.any(ini>ph)):
        raise ValueError("Positive bounded profile triples per band required")
    ini,pl,ph=(np.log(x.ravel()) for x in (ini,pl,ph))
    train,test=train[used],test[used]
    results={}
    for mode,ng in (("fixed",0),("constant",4),("variable",8)):
        mid0=base
        prof0=ini
        if mode!="fixed" and results["fixed"]["success"]:
            prof0=results["fixed"]["log_profile_parameters"].ravel()
        if mode=="variable" and results["constant"]["success"]:
            mid0=np.asarray(results["constant"]["midpoint"])
            prof0=results["constant"]["log_profile_parameters"].ravel()
        start=np.r_[mid0[:ng] if ng else [], np.zeros(4) if ng==8 else [],prof0]
        # For variable mode mid0 has four entries, then four deltas.
        lower=np.r_[lo if ng else [],-limit if ng==8 else [],pl]
        upper=np.r_[hi if ng else [],limit if ng==8 else [],ph]
        cache=[None,None]
        def geometry(z):
            mid=base if ng==0 else z[:4]
            delta=z[4:8] if ng==8 else np.zeros(4)
            return NestedEllipse(EllipseNode(*(mid-delta/2)),
                                 EllipseNode(*(mid+delta/2)),a0,a1)
        def predict(z):
            geo=geometry(z)
            key=z[:ng].copy()
            if cache[0] is None or not np.array_equal(cache[0],key):
                cache[:]=[key,_radius(geo,ee,nn)]
            radius=cache[1]
            return np.array([mean(fftconvolve(sersic_intensity(radius,*np.exp(
                z[ng+3*i:ng+3*i+3])),k,mode="same")) for i,k in enumerate(psfs)])
        def residual(z):
            pred=predict(z)
            r=((pred[:,train]-observed[:,train])/scales[:,None]).ravel()
            return np.r_[r,z[4:8]/ds] if ng==8 else r
        f=least_squares(residual,start,bounds=(lower,upper),x_scale="jac",
                        max_nfev=max_nfev,ftol=1e-8,xtol=1e-8,gtol=1e-8)
        pred=predict(f.x)
        diff=pred-observed
        weighted=(diff/scales[:,None])**2
        penalty=float(np.sum((f.x[4:8]/ds)**2)) if ng==8 else 0.
        results[mode]={"status":STATUS,"success":bool(f.success),
            "message":str(f.message),"boundary":bool(np.any(f.active_mask)),
            "jacobian_rank_with_penalty":int(np.linalg.matrix_rank(f.jac)),
            "geometry":asdict(geometry(f.x)),
            "midpoint":(base if ng==0 else f.x[:4]).copy(),
            "log_profile_parameters":f.x[ng:].reshape(len(bands),3).copy(),
            "training_weighted_score":float(weighted[:,train].sum()),
            "holdout_weighted_score":float(weighted[:,test].sum()),
            "training_mae_by_band":np.abs(diff[:,train]).mean(axis=1),
            "holdout_mae_by_band":np.abs(diff[:,test]).mean(axis=1),
            "regularization_penalty":penalty,"prediction":pred,"nfev":int(f.nfev)}
    return {"status":STATUS,"band_order":bands,"candidates":results,
            "observed_block_means":observed,"block_ids":used,
            "block_counts":cnt,"train":train,"holdout":test,"block_scales":scales.copy(),
            "block_size":b,"tile_shape":(ty,tx),"buffer_tiles":buffer_tiles,
            "configuration":{"seed":asdict(seed),"transition_start":a0,"transition_end":a1,
                "midpoint_bounds":[lo.tolist(),hi.tolist()],
                "delta_limits":limit.tolist(),"delta_scales":ds.tolist(),
                "profile_initial":np.exp(ini).reshape(len(bands),3).tolist(),
                "profile_bounds":[np.exp(pl).reshape(len(bands),3).tolist(),
                                  np.exp(ph).reshape(len(bands),3).tolist()],
                "min_block_fraction":float(min_block_fraction),"max_nfev":max_nfev,
                "max_work_pixels":max_work_pixels},
            "score_kind":"CONDITIONAL_DIAGONAL_BLOCK_SCORE_NOT_CALIBRATED_CHI2",
            "selected_geometry":None,"covariance_calibrated":False}


def fit_native_geometry_candidates(frames, scenes, *, anchor_icrs,
                                    corrections_by_band, psfs, psf_identities,
                                    region, region_grid_id, **configuration):
    """Native adapter; all masks and correction/PSF identities must match.

    psf_identities are ordered (native input SHA256, grid fingerprint) pairs.
    Region and tile selections are predeclared; this function does not infer
    them from image brightness or change the scene. Candidate geometries
    must be explicitly reviewed before reuse in annular photometry.
    """
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    from .contracts import NativeBandFrame, PSFDescription
    from .native_annular_diagnostic import NativeCorrection, common_native_validity
    from .sky_masks import native_boolean_mask_sha256
    from hashlib import sha256

    frames,scenes,psfs,psf_identities=map(tuple,(frames,scenes,psfs,psf_identities))
    corrections_by_band=tuple(tuple(c) for c in corrections_by_band)
    if not frames or any(len(x)!=len(frames) for x in
            (scenes,psfs,psf_identities,corrections_by_band)):
        raise ValueError("One explicit scene, correction set and PSF per native band required")
    if any(not isinstance(f,NativeBandFrame) for f in frames):
        raise ValueError("Validated native frames required")
    first=frames[0]
    grid=first.native_grid_fingerprint
    if region_grid_id!=grid:
        raise ValueError("Region grid identity mismatch")
    budget=configuration.get("max_work_pixels",500000)
    if not isinstance(budget,int) or budget<1 or np.prod(first.shape)>budget:
        raise ValueError("Native geometry exceeds explicit work-pixel budget")
    valid=common_native_validity(frames,scenes)
    bands=tuple(f.band for f in frames)
    ra,dec=map(float,anchor_icrs)
    if not np.isfinite([ra,dec]).all() or not 0<=ra<360 or not -90<=dec<=90:
        raise ValueError("Invalid ICRS anchor")
    center=SkyCoord(ra*u.deg,dec*u.deg,frame="icrs")
    yy,xx=np.indices(first.shape)
    sky=first.wcs.pixel_to_world(xx,yy).icrs
    separation=center.separation(sky)
    if not np.isfinite(separation.deg).all() or np.max(separation.deg)>.5:
        raise ValueError("Small-field geometry restricted to 0.5 degree")
    angle=center.position_angle(sky).rad
    east=separation.arcsec*np.sin(angle)
    north=separation.arcsec*np.cos(angle)
    images=[]; kernels=[]; inputs=[]
    for f,s,cs,psf,identity in zip(frames,scenes,corrections_by_band,psfs,psf_identities):
        if (not isinstance(psf,PSFDescription) or psf.state!="EMPIRICAL_KERNEL"
                or tuple(identity)!=(f.native_input_sha256,grid)):
            raise ValueError("PSF must refer to this native observation")
        if not cs or any(not isinstance(c,NativeCorrection) for c in cs):
            raise ValueError("Explicit native corrections required")
        kinds=[c.kind for c in cs]
        if "residual_sky" not in kinds or len(set(kinds))!=len(kinds):
            raise ValueError("One residual sky and no repeated correction kind required")
        image=np.array(f.image_nMgy,copy=True)
        for c in cs:
            if (c.grid_id!=grid or c.input_sha256!=f.native_input_sha256
                    or c.values.shape!=f.shape or not np.isfinite(c.values[valid]).all()):
                raise ValueError("Correction identity or usable-pixel mismatch")
            image[valid]-=c.values[valid]
        images.append(image)
        kernels.append(psf.kernel)
        inputs.append({"dataset_id":f.dataset_id,"native_input_sha256":f.native_input_sha256,
                       "scene_revision":s.revision,
                       "corrections":[{"kind":c.kind,"sha256":c.sha256,"source":c.source} for c in cs],
                       "psf_sha256":sha256(psf.kernel.tobytes()).hexdigest()})
    result=fit_block_geometry(east,north,images,kernels,band_order=bands,
                               valid=valid,region=region,**configuration)
    result.update({"grid_id":grid,"anchor_icrs":[ra,dec],"inputs":inputs,
                   "common_valid_sha256":native_boolean_mask_sha256(valid),
                   "region_sha256":native_boolean_mask_sha256(region),
                   "train_tiles_sha256":native_boolean_mask_sha256(configuration["train_tiles"]),
                   "holdout_tiles_sha256":native_boolean_mask_sha256(configuration["holdout_tiles"]),
                   "geometry_adopted":False})
    return result


def reviewed_geometry_candidate(result, mode, *, confirmed):
    """Return an explicit candidate for further diagnostic measurement only.

    Does not adopt photometry, change a scene, or automatically select a
    winner. Preserve the full result separately as fit provenance.
    """
    if confirmed is not True or result.get("status")!=STATUS:
        raise ValueError("Explicit researcher confirmation of diagnostic geometry required")
    if mode not in ("fixed","constant","variable"):
        raise ValueError("Unknown geometry candidate")
    candidate=result["candidates"][mode]
    if not candidate["success"] or candidate["boundary"]:
        raise ValueError("Failed or boundary-constrained geometry requires separate investigation")
    g=candidate["geometry"]
    return NestedEllipse(EllipseNode(**g["inner"]),EllipseNode(**g["outer"]),
                         g["transition_start"],g["transition_end"])
