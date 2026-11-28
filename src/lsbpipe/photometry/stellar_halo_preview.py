"""SH1-C conditional native image previews and signed aperture sums."""
import numpy as np

APERTURES=(10.,15.,20.,25.,30.,35.,40.,45.)


def conditional_halo_images(image, invvar, native_hard, source_mask, *,
                            source_xy,target_xy,pixel_scale_arcsec,comparison,radial_edges):
    image=np.asarray(image,dtype=float)
    iv=np.asarray(invvar,dtype=float)
    hard=np.asarray(native_hard,dtype=bool)
    excluded=np.asarray(source_mask,dtype=bool)
    if image.ndim!=2 or any(v.shape!=image.shape for v in (iv,hard,excluded)):
        raise ValueError("SH1-C requires aligned native arrays")
    candidates=comparison["candidates"]
    if not candidates:
        return {"status":"SH1C_NO_MODEL_SUPPORT_NOT_ADOPTED","science_adopted":False}
    sx,sy=map(float,source_xy); tx,ty=map(float,target_xy)
    scale=float(pixel_scale_arcsec)
    sep=np.hypot(tx-sx,ty-sy)
    if not np.isfinite([sx,sy,tx,ty,scale]).all() or scale<=0 or sep<=0:
        raise ValueError("SH1-C invalid geometry")
    if not (0<=tx<image.shape[1] and 0<=ty<image.shape[0]):
        raise ValueError("SH1-C target outside image")
    lo=max(c["training_r_range_arcsec"][0] for c in candidates)
    hi=min(c["training_r_range_arcsec"][1] for c in candidates)
    if not np.isfinite([lo,hi]).all() or hi<=lo:
        return {"status":"SH1C_NO_COMMON_SUPPORT_NOT_ADOPTED","science_adopted":False}
    yy,xx=np.ogrid[:image.shape[0],:image.shape[1]]
    x=(xx-sx)*scale; y=(yy-sy)*scale
    r=np.hypot(x,y); rt=np.hypot(xx-tx,yy-ty)*scale
    native=~hard & np.isfinite(image) & np.isfinite(iv) & (iv>0)
    valid=native & ~excluded
    support=(r>=lo)&(r<=hi)
    use=valid&support
    if not use.any():
        return {"status":"SH1C_NO_COMMON_PIXELS_NOT_ADOPTED","science_adopted":False}
    idx=min(range(len(candidates)),key=lambda i:(candidates[i]["side_median_abs_residual_nMgy_pix"],i))
    masks=[rt<=a for a in APERTURES]
    edge=min(tx,ty,image.shape[1]-1-tx,image.shape[0]-1-ty)*scale
    rows=[]
    for radius,ap in zip(APERTURES,masks):
        v=ap&valid; u=ap&use
        rows.append(dict(radius_arcsec=radius,in_frame_pixels=int(ap.sum()),
            aperture_fully_in_frame=bool(radius<=edge),valid_pixels=int(v.sum()),
            common_support_pixels=int(u.sum()),unsupported_valid_pixels=int((v&~support).sum()),
            raw_valid_flux_nMgy=float(image[v].sum()) if v.any() else None,
            raw_matched_flux_nMgy=float(image[u].sum()) if u.any() else None,
            model_fluxes=[]))
    hmin=np.full(image.shape,np.inf); hmax=np.full(image.shape,-np.inf)
    bmin=hmin.copy(); bmax=hmax.copy()
    edges=np.asarray(radial_edges,dtype=float)
    if edges.ndim!=1 or len(edges)<2 or not np.isfinite(edges).all() or np.any(np.diff(edges)<=0):
        raise ValueError("SH1-C invalid radial edges")
    cosphi=np.divide(((tx-sx)*x+(ty-sy)*y)/sep,r,out=np.zeros_like(r),where=r>0)
    for i,c in enumerate(candidates):
        k=np.asarray(c["radial_knots_arcsec"],dtype=float)
        h=np.asarray(c["halo_at_knots_nMgy_pix"],dtype=float)
        coeff=np.array([c["sky_constant_nMgy_pix"],c["sky_x_nMgy_pix_per_arcsec"],c["sky_y_nMgy_pix_per_arcsec"]])
        if (k.ndim!=1 or len(k)<2 or h.shape!=k.shape or not np.isfinite(k).all()
            or not np.isfinite(h).all() or np.any(np.diff(k)<=0) or np.any(h<0)
            or np.any(np.diff(h)>1e-10) or abs(h[-1])>1e-10 or not np.isfinite(coeff).all()):
            raise ValueError("SH1-C malformed profile")
        halo=np.interp(r,k,h,left=np.nan,right=0.)
        sky=coeff[0]+coeff[1]*x+coeff[2]*y
        if not (np.isfinite(halo[use]).all() and np.isfinite(sky[use]).all()):
            raise ValueError("SH1-C declared radial support mismatch")
        np.minimum(hmin,halo,out=hmin,where=use); np.maximum(hmax,halo,out=hmax,where=use)
        np.minimum(bmin,sky,out=bmin,where=use); np.maximum(bmax,sky,out=bmax,where=use)
        lookup={(v["sector"],v["r_inner_arcsec"]):v for v in c.get("side_radial_residuals",[])}
        amp=np.full(image.shape,np.nan)
        for a,b in zip(edges[:-1],edges[1:]):
            l=lookup.get((2,a)); z=lookup.get((6,a))
            if l is None or z is None or min(l["n"],z["n"])<3:
                continue
            delta=.5*abs(l["observed_minus_model_median_nMgy_pix"]-z["observed_minus_model_median_nMgy_pix"])
            if np.isfinite(delta):
                amp[(r>=a)&(r<b)]=delta
        hp=np.maximum(0.,halo+amp*cosphi); hm=np.maximum(0.,halo-amp*cosphi)
        for ap,row in zip(masks,rows):
            u=ap&use; stress=u&np.isfinite(amp)
            n=int(u.sum()); ns=int(stress.sum())
            row["model_fluxes"].append(dict(candidate_index=i,
                halo_flux_nMgy=float(halo[u].sum()) if n else None,
                sky_flux_nMgy=float(sky[u].sum()) if n else None,
                halo_corrected_flux_nMgy=float((image[u]-halo[u]).sum()) if n else None,
                halo_sky_residual_flux_nMgy=float((image[u]-halo[u]-sky[u]).sum()) if n else None,
                stress_supported_pixels=ns,stress_common_pixel_fraction=ns/n if n else None,
                stress_radial_residual_flux_nMgy=float((image[stress]-halo[stress]-sky[stress]).sum()) if ns else None,
                stress_plus_residual_flux_nMgy=float((image[stress]-hp[stress]-sky[stress]).sum()) if ns else None,
                stress_minus_residual_flux_nMgy=float((image[stress]-hm[stress]-sky[stress]).sum()) if ns else None))
        if i==idx:
            layers=dict(halo_model=halo,sky_model=sky,halo_corrected=image-halo,residual=image-halo-sky)
    for row in rows:
        for f in ("halo_flux_nMgy","sky_flux_nMgy","halo_corrected_flux_nMgy","halo_sky_residual_flux_nMgy"):
            v=[m[f] for m in row["model_fluxes"] if m[f] is not None]
            row[f+"_range"]=[min(v),max(v)] if v else None
        row["full_valid_pixel_model_support"]=bool(row["aperture_fully_in_frame"] and
            row["valid_pixels"]>0 and row["unsupported_valid_pixels"]==0)
    layers={k:np.where(use,v,np.nan) for k,v in layers.items()}
    layers["original"]=np.where(native,image,np.nan)
    spread=np.full(image.shape,np.nan); np.subtract(hmax,hmin,out=spread,where=use)
    layers["halo_model_spread"]=spread
    spread=np.full(image.shape,np.nan); np.subtract(bmax,bmin,out=spread,where=use)
    layers["sky_model_spread"]=spread
    for v in layers.values():
        v.setflags(write=False)
    return dict(schema="desidera.sh1c-conditional-native-preview/1",
        status="SH1C_CONDITIONAL_PREVIEW_NOT_ADOPTED",science_adopted=False,
        original_pixels_modified=False,diagnostic_subtraction_only=True,
        uncertainty_calibrated=False,reference_candidate_index=idx,
        reference_rule="MIN_SIDE_MEDIAN_ABS_RESIDUAL_DISPLAY_ONLY",
        candidate_count=len(candidates),common_radial_support_arcsec=[lo,hi],
        source_xy_native_pix=[sx,sy],target_xy_native_pix=[tx,ty],
        pixel_scale_arcsec=scale,valid_pixel_count=int(valid.sum()),
        common_supported_pixel_count=int(use.sum()),apertures=rows,layers=layers,
        stress_rule="PLUS_MINUS_HALF_SIDE_RESIDUAL_CONTRAST_TIMES_COS_TARGET_ANGLE; NONNEGATIVE_HALO; PAIRED_BINS_MIN_3_BLOCKS_PER_SIDE",
        warning="Conditional signed sums on matched observed pixels only. No masked-pixel fill, total flux, calibrated sky or error. Stress transfer is a hypothesis, not a bound.")
