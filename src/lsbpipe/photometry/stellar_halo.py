"""SH1-A: native-pixel sector measurements for a saturated foreground star.

No PSF extrapolation, science adoption, sky subtraction or pixel modification.
Optional SH1-B compares explicit, non-adopted halo + background hypotheses.
Target-facing and protected pixels cannot train the future halo model.
The angular sectors are defined relative to the fixed star-target vector.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.ndimage import binary_dilation

EDGES_ARCSEC = (12., 16., 20., 25., 30., 40., 50., 60., 70., 85., 100., 115.)
EXPANDED_EDGES_ARCSEC = EDGES_ARCSEC + (130., 145., 160., 180.)
EXPANDED_ANCHORS_ARCSEC = (70., 85., 100., 120., 150.)
MIN_PIXELS = 100
MIN_COVERAGE = .25


def _stats(image, select):
    v = image[select]
    if not v.size:
        return {"n": 0, "median_nMgy_pix": None, "mad_nMgy_pix": None,
                "q16_nMgy_pix": None, "q84_nMgy_pix": None}
    m = float(np.median(v))
    q = np.quantile(v, [.16, .84])
    return {"n": int(v.size), "median_nMgy_pix": m,
            "mad_nMgy_pix": float(np.median(np.abs(v-m))),
            "q16_nMgy_pix": float(q[0]), "q84_nMgy_pix": float(q[1])}


def measure_saturated_halo_sectors(image, invvar, native_hard, *,
        source_xy, target_xy, pixel_scale_arcsec, source_mask, protected,
        unresolved_sky_veto, neighbor_centres=(), compare_models=False,
        expanded=False):
    a = np.asarray(image, dtype=float)
    iv = np.asarray(invvar, dtype=float)
    masks = [np.asarray(m, dtype=bool) for m in
             (native_hard, source_mask, protected, unresolved_sky_veto)]
    if a.ndim != 2 or iv.shape != a.shape or any(m.shape != a.shape for m in masks):
        raise ValueError("SH1 requires aligned native image, INVVAR and scene masks")
    sx, sy = map(float, source_xy)
    tx, ty = map(float, target_xy)
    scale = float(pixel_scale_arcsec)
    if (not all(math.isfinite(v) for v in (sx,sy,tx,ty,scale))
            or not .1 <= scale <= 1.
            or not (0 <= sx < a.shape[1] and 0 <= sy < a.shape[0])
            or math.hypot(tx-sx,ty-sy)*scale < 1.):
        raise ValueError("SH1 invalid source/target/scale")
    yy, xx = np.ogrid[:a.shape[0], :a.shape[1]]
    dx, dy = xx-sx, yy-sy
    distance = np.hypot(dx,dy)*scale
    direction = math.atan2(ty-sy, tx-sx)
    angle = (np.arctan2(dy,dx)-direction+np.pi/8)%(2*np.pi)
    sector = np.floor(angle/(np.pi/4)).astype(np.int8)
    hard, source, protection, veto = masks
    # Explicit engineering guard around invalid pixels. Do not allow
    # the protection mask to restore saturation/bleed pixels.
    guarded_hard = binary_dilation(hard, iterations=3)
    native_valid = (~guarded_hard & np.isfinite(a) & np.isfinite(iv) & (iv > 0))
    # Mask catalogue cores even when automatic photometry screening
    # left them protected or unfitted. This is only a halo measurement mask.
    cores = np.zeros(a.shape, dtype=bool)
    n_neighbours = 0
    for cx,cy in neighbor_centres:
        cx,cy = float(cx),float(cy)
        if not math.isfinite(cx) or not math.isfinite(cy):
            continue
        radius = 4./scale
        x0=max(0, math.floor(cx-radius)); x1=min(a.shape[1],math.ceil(cx+radius)+1)
        y0=max(0, math.floor(cy-radius)); y1=min(a.shape[0],math.ceil(cy+radius)+1)
        if x1 <= x0 or y1 <= y0:
            continue
        cyy,cxx=np.ogrid[y0:y1,x0:x1]
        cores[y0:y1,x0:x1] |= (cxx-cx)**2+(cyy-cy)**2 <= radius**2
        n_neighbours += 1
    descriptive = native_valid & ~source & ~cores
    train_valid = descriptive & ~protection & ~veto
    edges = EXPANDED_EDGES_ARCSEC if expanded else EDGES_ARCSEC
    rows=[]
    for lo,hi in zip(edges[:-1], edges[1:]):
        annulus = (distance >= lo) & (distance < hi)
        nominal = math.pi*(hi**2-lo**2)/(8*scale**2)
        for k in range(8):
            geom = annulus & (sector == k)
            # only opposite three sectors train; side sectors 2,6
            # remain independent checks; 0,1,7 face the galaxy.
            role = "AWAY_TRAIN_CANDIDATE" if k in (3,4,5) else (
                "SIDE_CONTROL" if k in (2,6) else "TARGET_FACING_DESCRIPTIVE")
            desc = _stats(a, geom & descriptive)
            safe = _stats(a, geom & train_valid)
            fraction = min(1.,safe["n"]/nominal)
            rows.append({"r_inner_arcsec":lo, "r_outer_arcsec":hi,
                         "sector_index":k,"relative_angle_center_deg":45*k,
                         "role":role,"nominal_full_sector_pixels_approx":nominal,
                         "in_frame_pixels":int(geom.sum()),
                         "native_valid_pixels":int((geom & native_valid).sum()),
                         "descriptive":desc,"outside_protection_and_veto":safe,
                         "safe_full_sector_fraction_approx":fraction,
                         "adequate_for_descriptive_comparison":
                             safe["n"] >= MIN_PIXELS and fraction >= MIN_COVERAGE})
    comparison = (_compare_block_models(a, train_valid, distance, sector, rows,
                  sx=sx, sy=sy, scale=scale, target_distance=math.hypot(tx-sx,ty-sy)*scale,
                  edges=edges, expanded=expanded)
                  if compare_models else None)
    if expanded and comparison is not None:
        comparison["diagnostic_tile_map"] = _diagnostic_tile_map(
            a, descriptive, protection, distance, sx=sx,sy=sy,scale=scale,
            candidates=comparison["candidates"])
    return {"schema":"desidera.sh1-saturated-halo-sector-measurement/1",
            "status":"SH1_SECTORS_MEASURED_NOT_ADOPTED",
            "source_xy_native_pix":[sx,sy],"target_xy_native_pix":[tx,ty],
            "pixel_scale_arcsec":scale,"star_target_distance_arcsec":
                math.hypot(tx-sx,ty-sy)*scale,
            "radial_edges_arcsec":list(edges),"expanded_radial_experiment":bool(expanded),
            "native_invalid_guard_pix":3,"catalogue_core_radius_arcsec":4.,
            "catalogue_neighbour_count":n_neighbours,
            "min_pixels":MIN_PIXELS,"min_full_sector_fraction_approx":MIN_COVERAGE,
            "rows":rows,
            "catalog_flux_used":False,"compact_psf_used":False,
            "sky_subtracted":False,"halo_fitted":bool(comparison and comparison["candidates"]),
            "sh1b_model_comparison":comparison,"science_adopted":False,
            "original_pixels_modified":False,"noise_covariance_calibrated":False,
            "warning":"Signed coadd medians include unknown sky; coverage and MAD are descriptive, not calibrated uncertainties. Target-facing values may include galaxy light. No adopted halo amplitude or galaxy correction inferred."}


def _spatial_blocks(image, valid, distance, sector, rows, *, sx, sy, scale, step, edges=EDGES_ARCSEC):
    """Disjoint native tiles; only accepted sector/ring combinations enter.
    Within a tile, aggregate each ring/sector subset separately. Side subsets
    sharing a tile with ANY training subset are not independent holdouts.
    """
    allowed = {(r["r_inner_arcsec"], r["sector_index"])
               for r in rows if r["adequate_for_descriptive_comparison"]
               and r["sector_index"] in (2,3,4,5,6)}
    blocks = []
    for y0 in range(0, image.shape[0], step):
        for x0 in range(0, image.shape[1], step):
            sl = np.s_[y0:y0+step,x0:x0+step]
            y,x = np.nonzero(valid[sl])
            if len(x) < 32:
                continue
            x=x+x0; y=y+y0
            radii=distance[y,x]; sectors=sector[y,x]
            rings=np.searchsorted(edges,radii,side="right")-1
            group=8*rings+sectors
            for code in np.unique(group):
                ring,k=divmod(int(code),8)
                if not 0<=ring<len(edges)-1:
                    continue
                lo=edges[ring]
                if (lo,k) not in allowed:
                    continue
                q=group==code
                if int(q.sum()) < max(32, int(.15*step*step)):
                    continue
                v=image[y[q],x[q]]
                blocks.append({"tile":[x0,y0],"sector":k,"r_inner_arcsec":lo,
                               "n":int(q.sum()),"r_arcsec":float(np.median(radii[q])),
                               "x_arcsec":float(np.median(x[q]-sx)*scale),
                               "y_arcsec":float(np.median(y[q]-sy)*scale),
                               "intensity_nMgy_pix":float(np.median(v))})
    train_tiles={tuple(b["tile"]) for b in blocks if b["sector"] in (3,4,5)}
    for b in blocks:
        b["role"]=("TRAIN" if b["sector"] in (3,4,5) else
                   "SHARED_TILE_EXCLUDED" if tuple(b["tile"]) in train_tiles else "SIDE_HOLDOUT")
    return blocks


def _radial_basis(r, knots):
    # Positive coefficients are drops between successive radial knots.
    # The imposed zero at the final knot is a hypothesis, not a measurement.
    r=np.asarray(r,dtype=float)
    return np.clip((knots[1:]-r[...,None])/np.diff(knots),0.,1.)


def _block_design(blocks, knots, plane):
    r=np.array([b["r_arcsec"] for b in blocks])
    cols=[_radial_basis(r,knots),np.ones((len(blocks),1))]
    if plane:
        cols.append(np.array([[b["x_arcsec"]/100.,b["y_arcsec"]/100.]
                              for b in blocks]).reshape(-1,2))
    return np.column_stack(cols)


def _compare_block_models(image, valid, distance, sector, rows, *, sx,sy,scale,target_distance,
                          edges=EDGES_ARCSEC, expanded=False):
    """Hypothesis/sensitivity comparison, never an adopted correction.
    The additive halo/sky degeneracy is exposed via several outer anchors.
    Equal spatial-subset weights, signed block medians, robust residual loss.
    """
    from scipy.optimize import least_squares
    candidates=[]; rejected=[]; block_sets=[]
    for step in (24,40):
        blocks=_spatial_blocks(image,valid,distance,sector,rows,sx=sx,sy=sy,scale=scale,step=step,edges=edges)
        train=[b for b in blocks if b["role"]=="TRAIN"]
        side=[b for b in blocks if b["role"]=="SIDE_HOLDOUT"]
        block_sets.append({"tile_size_pix":step,"blocks":blocks,
                           "n_training":len(train),"n_side_holdout":len(side)})
        for anchor in (EXPANDED_ANCHORS_ARCSEC if expanded else (70.,85.,100.)):
            for plane in (False,True):
                meta={"tile_size_pix":step,"assumed_zero_halo_radius_arcsec":anchor,
                      "sky_model":"plane" if plane else "constant"}
                outer=[b for b in train if b["r_arcsec"]>=anchor]
                inner=[b for b in train if b["r_arcsec"]<anchor]
                if (len(train)<25 or len(side)<6 or len(outer)<6
                    or len({b["sector"] for b in outer})<2 or len(inner)<12):
                    rejected.append(dict(meta,reason="INSUFFICIENT_BLOCK_OR_OUTER_ANGULAR_SUPPORT"))
                    continue
                rmin=min(b["r_arcsec"] for b in train)
                rmax=max(b["r_arcsec"] for b in train)
                interior = (16.,20.,25.,30.,40.,50.,60.,70.,85.,100.,115.,130.,145.) if expanded else (16.,20.,25.,30.,40.,50.,60.,70.,85.)
                knots=np.array([rmin]+[k for k in interior if rmin+2.<k<anchor]+[anchor])
                X=_block_design(train,knots,plane)
                y=np.array([b["intensity_nMgy_pix"] for b in train])*1000.
                if np.linalg.matrix_rank(X)<X.shape[1] or len(train)<3*X.shape[1]:
                    rejected.append(dict(meta,reason="RANK_OR_SAMPLE_COUNT"))
                    continue
                cond=float(np.linalg.cond(X))
                if not np.isfinite(cond) or cond>1e5:
                    rejected.append(dict(meta,reason="ILL_CONDITIONED_DESIGN"))
                    continue
                nd=len(knots)-1
                lower=np.r_[np.zeros(nd),np.full(X.shape[1]-nd,-np.inf)]
                init=np.zeros(X.shape[1]); init[nd]=np.median(y)
                fit=least_squares(lambda p:X@p-y,init,bounds=(lower,np.full(X.shape[1],np.inf)),
                                  loss="soft_l1",f_scale=.5,max_nfev=300)
                if not fit.success or not np.isfinite(fit.x).all():
                    rejected.append(dict(meta,reason="ROBUST_SOLVER_FAILED"))
                    continue
                p=fit.x/1000.
                residual=X@p-y/1000.
                controls=[b for b in side if rmin<=b["r_arcsec"]<=rmax]
                if len(controls)<6 or (expanded and any(
                        sum(b['sector']==k for b in controls)<3 for k in (2,6))):
                    rejected.append(dict(meta,reason="INSUFFICIENT_IN_RANGE_HOLDOUT"))
                    continue
                sr=_block_design(controls,knots,plane)@p-np.array([b["intensity_nMgy_pix"] for b in controls])
                radii=sorted({float(r) for r in (*edges,target_distance) if rmin<=r<=rmax})
                profile=[{"radius_arcsec":r,"halo_nMgy_pix":float(_radial_basis(r,knots)@p[:nd])} for r in radii]
                target=(float(_radial_basis(target_distance,knots)@p[:nd])
                        if rmin<=target_distance<=rmax else None)
                candidates.append(dict(meta,radial_knots_arcsec=knots.tolist(),
                    halo_at_knots_nMgy_pix=(_radial_basis(knots,knots)@p[:nd]).tolist(),
                    sky_constant_nMgy_pix=float(p[nd]),
                    sky_x_nMgy_pix_per_arcsec=float(p[nd+1]/100) if plane else 0.,
                    sky_y_nMgy_pix_per_arcsec=float(p[nd+2]/100) if plane else 0.,
                    training_r_range_arcsec=[rmin,rmax],design_condition=cond,
                    n_training=len(train),n_outer=len(outer),n_side_holdout=len(controls),
                    train_median_abs_residual_nMgy_pix=float(np.median(np.abs(residual))),
                    side_median_abs_residual_nMgy_pix=float(np.median(np.abs(sr))),
                    side_signed_median_residual_nMgy_pix=float(np.median(sr)),
                    side_by_sector=[{"sector":k,"n":len(v),
                        "signed_median_residual_nMgy_pix":float(np.median(v)),
                        "median_abs_residual_nMgy_pix":float(np.median(np.abs(v)))}
                        for k in (2,6) if (v:=sr[np.array([b["sector"]==k for b in controls])]).size],
                    side_radial_residuals=[{"sector":k,"r_inner_arcsec":lo,
                        "n":len(v),"observed_minus_model_median_nMgy_pix":float(-np.median(v)),
                        "median_abs_residual_nMgy_pix":float(np.median(np.abs(v)))}
                        for k in (2,6) for lo in edges[:-1]
                        if (v:=sr[np.array([b["sector"]==k and b["r_inner_arcsec"]==lo
                                            for b in controls],dtype=bool)]).size],
                    halo_at_target_distance_nMgy_pix=target,profile=profile))
    values=[c["halo_at_target_distance_nMgy_pix"] for c in candidates
            if c["halo_at_target_distance_nMgy_pix"] is not None]
    return {"schema":"desidera.sh1b-block-model-comparison/1",
            "status":"SH1B_HYPOTHESES_COMPARED_NOT_ADOPTED" if candidates else "SH1B_INSUFFICIENT_SUPPORT",
            "candidates":candidates,"rejected":rejected,"block_sets":block_sets,
            "target_distance_hypothesis_range_nMgy_pix":[min(values),max(values)] if values else None,
            "science_adopted":False,"source_light_subtracted":False,
            "outer_zero_is_assumption":True,"uncertainty_calibrated":False,
            "warning":"No unique halo/sky separation. Outer zero anchors are unverified hypotheses, not measurements. Block medians use median coordinates and equal weights; coadd covariance is not calibrated. Side sectors are disjoint-tile checks, not independent exposures. No model selected, extrapolated correction, or pixels modified."}


def _diagnostic_tile_map(image, valid, protected, distance, *, sx,sy,scale,candidates):
    """Native-pixel evaluation aggregated to 48px tiles, not a science FITS.
    Every entry uses common radial training support of ALL candidates.
    Spread is the range of conditional hypotheses, not a calibrated error.
    """
    step=48
    out={"schema":"desidera.sh1b2-diagnostic-tile-map/1","tile_size_pix":step,
         "native_shape":list(image.shape),"common_training_radius_arcsec":None,
         "cells":[],"science_adopted":False,"source_light_subtracted":False,
         "warning":"Tile summaries, not a corrected image. Target-facing residuals include galaxy light. Hypothesis spread is not a confidence interval; no extrapolation outside common radial support."}
    if not candidates:
        return out
    lo=max(c["training_r_range_arcsec"][0] for c in candidates)
    hi=min(c["training_r_range_arcsec"][1] for c in candidates)
    if hi<=lo:
        return out
    out["common_training_radius_arcsec"]=[lo,hi]
    for y0 in range(0,image.shape[0],step):
        for x0 in range(0,image.shape[1],step):
            sl=np.s_[y0:y0+step,x0:x0+step]
            rr=distance[sl]
            good=valid[sl] & (rr>=lo) & (rr<=hi)
            n=int(good.sum())
            if n<max(32,.15*good.size):
                continue
            y,x=np.nonzero(good)
            r=rr[good]; x=(x+x0-sx)*scale; y=(y+y0-sy)*scale
            values=image[sl][good]
            halo=[]; sky=[]; residual=[]
            for c in candidates:
                h=np.interp(r,c["radial_knots_arcsec"],c["halo_at_knots_nMgy_pix"],right=0.)
                b=c["sky_constant_nMgy_pix"]+x*c["sky_x_nMgy_pix_per_arcsec"]+y*c["sky_y_nMgy_pix_per_arcsec"]
                halo.append(float(np.median(h)))
                sky.append(float(np.median(b)))
                residual.append(float(np.median(values-h-b)))
            out["cells"].append({"x0":x0,"y0":y0,"n":n,
                "valid_fraction":n/good.size,
                "protected_fraction_of_used":float(protected[sl][good].mean()),
                "observed_median_nMgy_pix":float(np.median(values)),
                "halo_median_range_nMgy_pix":[min(halo),max(halo)],
                "sky_median_range_nMgy_pix":[min(sky),max(sky)],
                "residual_median_range_nMgy_pix":[min(residual),max(residual)],
                "residual_median_across_hypotheses_nMgy_pix":float(np.median(residual))})
    return out
