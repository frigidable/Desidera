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


def native_scene_block_covariance(frames, records, halo_fits, *,
                                  length_arcsec, mix, minimum_control_tiles,
                                  clip_sigma, max_tiles=4096):
    """Build an additive conditional covariance from exact native ring support.

    Each halo fit must be the independent sky/halo fit whose corrections
    were used in the corresponding annular record. Far held-out controls
    must have identical tile locations AND valid-pixel footprints across
    bands. The common spatial component is assumed uniform within a tile;
    formal IVAR is treated as independent pixel noise. Neither assumption
    constitutes a calibrated coadd covariance. Compact source corrections
    and geometry remain fixed.
    """
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    from .contracts import NativeBandFrame
    from .native_annular_diagnostic import numeric_sha256, annular_record_sha256
    from .stellar_contamination import compact_halo_template

    frames,records,halo_fits=map(tuple,(frames,records,halo_fits))
    if not frames or len(frames)!=len(records) or len(frames)!=len(halo_fits):
        raise ValueError("One native record and halo fit per band required")
    if any(not isinstance(f,NativeBandFrame) for f in frames):
        raise ValueError("Validated native frames required")
    first=records[0]
    bands=tuple(f.band for f in frames)
    if len(set(bands))!=len(bands):
        raise ValueError("Repeated band")
    size=halo_fits[0]["block_size"]
    if not isinstance(size,int) or size<2:
        raise ValueError("Invalid block size")
    if not isinstance(max_tiles,int) or max_tiles<1:
        raise ValueError("Invalid tile budget")
    nr=len(first["ring_edges_arcsec"])-1
    tile_sets=[]; band_inputs=[]; control_maps=[]
    for frame,rec,fit in zip(frames,records,halo_fits):
        if (rec["status"]!=STATUS or fit["status"]!=STATUS
                or rec["dataset_id"]!=frame.dataset_id or fit["dataset_id"]!=frame.dataset_id
                or rec["grid_id"]!=frame.native_grid_fingerprint or fit["grid_id"]!=rec["grid_id"]
                or rec["native_input_sha256"]!=frame.native_input_sha256
                or fit["native_input_sha256"]!=rec["native_input_sha256"]
                or frame.dataset_id.rsplit("/",1)[0]!=frames[0].dataset_id.rsplit("/",1)[0]
                or rec["grid_id"]!=first["grid_id"] or rec["valid_sha256"]!=first["valid_sha256"]
                or rec["ring_index_sha256"]!=first["ring_index_sha256"]
                or rec["geometry"]!=first["geometry"] or rec["anchor_icrs"]!=first["anchor_icrs"]
                or not np.array_equal(rec["ring_edges_arcsec"],first["ring_edges_arcsec"])
                or fit["block_size"]!=size):
            raise ValueError("Native grid, band, geometry, support or block identity mismatch")
        corr={c["kind"]:c["sha256"] for c in rec["corrections"]}
        if (corr.get("residual_sky")!=fit["sky"].sha256
                or corr.get("stellar_halo")!=fit["halo"].sha256
                or corr.get("compact_sources")!=fit["compact_sources_sha256"]):
            raise ValueError("Annular record used different contamination corrections")
        ring=np.asarray(rec.get("ring_index"))
        if (ring.shape!=frame.shape or ring.dtype.kind not in "iu"
                or numeric_sha256(ring)!=rec["ring_index_sha256"]
                or np.any(ring < -1) or np.any(ring>=nr)):
            raise ValueError("Retain the exact native pixel ring support first")
        yy,xx=np.nonzero(ring>=0)
        rr=ring[yy,xx]
        if not len(rr) or not np.array_equal(np.bincount(rr,minlength=nr),rec["valid_pixel_count"]):
            raise ValueError("Ring support count mismatch")
        iv=frame.invvar_nMgy2[yy,xx]
        if not np.isfinite(iv).all() or np.any(iv<=0):
            raise ValueError("Invalid native IVAR on ring support")
        variances=np.bincount(rr,weights=1/iv,minlength=nr)
        if not np.allclose(variances,rec["formal_diagonal_variance_nMgy2"],rtol=1e-10,atol=1e-15):
            raise ValueError("Ring support variance mismatch")
        # The same sky-offset coordinates and analytic basis as the fit.
        star=SkyCoord(*fit["star_icrs"],unit=u.deg,frame="icrs")
        sky=frame.wcs.pixel_to_world(xx,yy).icrs
        rad=star.separation(sky).arcsec
        pa=star.position_angle(sky).rad
        h=compact_halo_template(rad,cutoff_arcsec=fit["cutoff_arcsec"],
            power=fit["power"],pivot_arcsec=fit["pivot_arcsec"],core_arcsec=fit["core_arcsec"])
        if not np.isfinite(h).all():
            raise ValueError("Observed ring intersects excluded stellar core")
        sums=[np.bincount(rr,minlength=nr)]
        if fit["sky_degree"]==1:
            for val in (rad*np.sin(pa),rad*np.cos(pa)):
                sums.append(np.bincount(rr,weights=val/fit["far_star_guard_arcsec"],minlength=nr))
        design=np.asarray(sums,dtype=float).T
        halo=np.bincount(rr,weights=h,minlength=nr)
        # Check the reconstructed bases against the exact correction rasters.
        if not np.allclose(design@fit["sky_coefficients"],
                np.bincount(rr,weights=fit["sky"].values[yy,xx],minlength=nr),rtol=1e-9,atol=1e-10):
            raise ValueError("Sky design does not reproduce the applied correction")
        if not np.allclose(halo*fit["amplitude"],
                np.bincount(rr,weights=fit["halo"].values[yy,xx],minlength=nr),rtol=1e-9,atol=1e-10):
            raise ValueError("Halo template does not reproduce the applied correction")
        nt=len(fit["block_xy"])
        bxy=np.asarray(fit["block_xy"])
        if (bxy.shape!=(nt,2) or bxy.dtype.kind not in "iu" or np.any(bxy<0)
                or np.any(bxy%size) or np.any(bxy[:,0]+size>frame.shape[1])
                or np.any(bxy[:,1]+size>frame.shape[0])):
            raise ValueError("Invalid fitted block coordinates")
        keys=[tuple(map(int,row)) for row in bxy]
        if len(set(keys))!=nt:
            raise ValueError("Repeated fitted block")
        train=np.asarray(fit["far"])|np.asarray(fit["near"])
        held=np.asarray(fit["far_control"])
        if train.shape!=(nt,) or held.shape!=(nt,) or train.dtype.kind!="b" or held.dtype.kind!="b" or np.any(train&held):
            raise ValueError("Invalid disjoint block roles")
        b=np.asarray(fit["sky_influence"],dtype=float)
        a=np.asarray(fit["halo_influence"],dtype=float)
        if b.shape!=(design.shape[1],nt) or a.shape!=(nt,) or np.any(b[:,~train]!=0) or np.any(a[~train]!=0):
            raise ValueError("Invalid fitted influence support")
        q= -design@b-np.outer(halo,a)
        v=np.asarray(fit["block_formal_mean_variance"],dtype=float)
        if v.shape!=(nt,) or not np.isfinite(v).all() or np.any(v<=0):
            raise ValueError("Missing formal block variance")
        # Integer tile indices include incomplete edge tiles for observations.
        obs_xy=np.column_stack((xx//size*size,yy//size*size))
        obs_keys,inv=np.unique(obs_xy,axis=0,return_inverse=True)
        obs_keys=[tuple(map(int,row)) for row in obs_keys]
        if set(obs_keys)&{keys[i] for i in np.flatnonzero(train|held)}:
            raise ValueError("Observed ring intersects training or far-control tile")
        counts=np.bincount(rr*len(obs_keys)+inv,minlength=nr*len(obs_keys)).reshape(nr,-1)
        used=set(obs_keys)|{keys[i] for i in np.flatnonzero(train)}
        tile_sets.append(used)
        band_inputs.append((obs_keys,counts,keys,q,train,independent_ivar_covariance(variances,q,v)))
        residual=np.asarray(fit["block_means"])-np.asarray(fit["block_sky_design"])@fit["sky_coefficients"]
        fingerprints=tuple(fit["block_footprint_sha256"])
        if residual.shape!=(nt,) or len(fingerprints)!=nt:
            raise ValueError("Incomplete control block provenance")
        control_maps.append({keys[i]:(fingerprints[i],residual[i],v[i])
                             for i in np.flatnonzero(held)})
    matched=set(control_maps[0])
    for m in control_maps[1:]:
        matched &= set(m)
    matched=sorted(k for k in matched if len({m[k][0] for m in control_maps})==1)
    if not matched:
        raise ValueError("No far-control tiles with identical native support")
    control=np.array([[m[k][1] for m in control_maps] for k in matched])
    formal=np.array([[m[k][2] for m in control_maps] for k in matched])
    stats=paired_control_statistics(control,band_order=bands,tile_ids=matched,
                                    minimum_tiles=minimum_control_tiles,clip_sigma=clip_sigma)
    split=split_control_block_noise(stats["scales"],stats["correlation"],formal)
    tiles=sorted(set.union(*tile_sets))
    if len(tiles)>max_tiles:
        raise ValueError("Native covariance exceeds tile budget")
    lookup={k:i for i,k in enumerate(tiles)}
    ops=[]; white=[]
    for obs,counts,keys,q,train,d in band_inputs:
        l=np.zeros((nr,len(tiles)))
        l[:,[lookup[k] for k in obs]]=counts
        ix=np.flatnonzero(train)
        l[:,[lookup[keys[i]] for i in ix]]+=q[:,ix]
        ops.append(l);white.append(d)
    xy=np.asarray(tiles,dtype=float)+(size-1)/2
    # Incomplete observed edge tiles use their nominal tile centers.
    sky=frames[0].wcs.pixel_to_world(xy[:,0],xy[:,1]).icrs
    anchor=SkyCoord(*first["anchor_icrs"],unit=u.deg,frame="icrs")
    rad=anchor.separation(sky).arcsec
    pa=anchor.position_angle(sky).rad
    result=conditional_block_covariance(np.column_stack((rad*np.sin(pa),rad*np.cos(pa))),
        ops,split["coherent_scales"],split["coherent_correlation"],white,
        band_order=bands,length_arcsec=length_arcsec,mix=mix,max_tiles=max_tiles)
    result.update({"record_sha256":tuple(annular_record_sha256(r) for r in records),
                   "matched_far_control_tiles":matched,"n_matched_far_controls":len(matched),
                   "tile_xy_native":tiles,"block_size":size,"control_statistics":stats,
                   "noise_split":split,"independent_ivar_covariances":white,
                   "covariance_model":"ADDITIVE_COHERENT_TILE_AND_INDEPENDENT_IVAR",
                   "control_sky_prediction_uncertainty_not_removed":True})
    return result
