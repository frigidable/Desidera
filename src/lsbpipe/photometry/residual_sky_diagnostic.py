"""SKY-1: conditional remaining coadd background on expanded native window.

All training/validation excludes existing protection, source masks and
additional target guard. Spatial holdouts are independent of THIS sky
fit, not necessarily of the preceding halo fit. No science adoption.
"""
import numpy as np
from scipy.spatial import Delaunay, QhullError
from ..background import BackgroundSettings, FitStatus, _fit_robust
from .stellar_halo_spikes import safe_pixels, model_arrays, key
from .galaxy_geometry_diagnostic import local_grid, block_design, fit_geometry, fitted, radial_profiles
from .sky_masks import native_boolean_mask_sha256

GUARDS=(45.,60.,75.)
TILES=(32,48)
METHODS={"constant":1,"plane":3,"quadratic":6}
AP_RADII=(5.,10.,15.)
NORM=100.


def features(x,y):
    x,y=np.broadcast_arrays(x/NORM,y/NORM)
    return np.stack((np.ones_like(x),x,y,x*x,x*y,y*y),axis=-1)


def quadrant(x,y):
    return (x>=0).astype(int)+2*(y>=0).astype(int)


def tile_rows(res,safe,tx,ty,scale,tile):
    rows=[]
    ny,nx=res.shape
    for y in range(0,ny-tile+1,tile):
        for x in range(0,nx-tile+1,tile):
            # Exclude blocks crossing validation quadrant boundaries.
            if (x<tx<=x+tile-1) or (y<ty<=y+tile-1): continue
            sl=np.s_[y:y+tile,x:x+tile]; m=safe[sl]
            if m.sum()<.5*tile*tile: continue
            yy,xx=np.nonzero(m)
            dx=(xx+x-tx)*scale; dy=(yy+y-ty)*scale
            # Block means of basis on exactly the observed pixels.
            v=res[sl][m]
            rows.append(dict(tile_x_pix=x,tile_y_pix=y,x_pix=float(xx.mean()+x),y_pix=float(yy.mean()+y),
                quadrant=int(quadrant(np.asarray(dx.mean()),np.asarray(dy.mean()))),
                n=int(m.sum()),mean_nMgy_pix=float(v.mean()),
                median_nMgy_pix=float(np.median(v)),design=features(dx,dy).mean(axis=0).tolist()))
    return rows


def fit_blocks(rows,method):
    k=METHODS[method]
    if len(rows)<max(12,3*k): return None
    X=np.asarray([r["design"][:k] for r in rows])
    z=np.asarray([r["mean_nMgy_pix"] for r in rows])
    scale=max(float(np.max(abs(z))),1e-12)
    # Equal tile weights; robust scatter, not per-pixel formal IVAR.
    f,cond,it,status=_fit_robust(X,z/scale,np.ones(len(z)),
        "constant" if k==1 else "tilted_plane",BackgroundSettings(max_design_condition=1e5))
    if status!=FitStatus.FITTED or f is None: return None
    return dict(coefficients_nMgy_pix=(f*scale).tolist(),design_condition=float(cond),iterations=it)


def predict(f,x,y):
    c=np.asarray(f["coefficients_nMgy_pix"])
    return np.sum(features(x,y)[...,:len(c)]*c,axis=-1)


def empty_samples(res,safe,tx,ty,scale):
    """Nonoverlapping circles within each radius, FULL safe pixel coverage.

    Different radii and halo scenarios may overlap; no independent noise
    claim. Samples evaluated by quadrant-held-out model later.
    """
    out=[]
    ny,nx=res.shape
    for radius in AP_RADII:
        r=radius/scale; pad=int(np.ceil(r)); step=2*pad+1
        yy,xx=np.mgrid[-pad:pad+1,-pad:pad+1]
        disk=xx*xx+yy*yy<=r*r
        for y in range(pad,ny-pad,step):
            for x in range(pad,nx-pad,step):
                if abs(x-tx)<=pad or abs(y-ty)<=pad: continue
                sl=np.s_[y-pad:y+pad+1,x-pad:x+pad+1]
                if not safe[sl][disk].all(): continue
                dx=(x+xx[disk]-tx)*scale; dy=(y+yy[disk]-ty)*scale
                out.append(dict(radius_arcsec=radius,x_pix=x,y_pix=y,
                    quadrant=int(quadrant(np.asarray(x-tx),np.asarray(y-ty))),
                    n=int(disk.sum()),raw_sum_nMgy=float(res[sl][disk].sum()),
                    design_sum=features(dx,dy).sum(axis=0).tolist()))
    return out


def summarize_apertures(samples,folds):
    rows=[]
    for radius in AP_RADII:
        vals=[]
        for s in samples:
            if s["radius_arcsec"]!=radius or s["quadrant"] not in folds: continue
            c=np.asarray(folds[s["quadrant"]]["coefficients_nMgy_pix"])
            val=s["raw_sum_nMgy"]-np.dot(s["design_sum"][:len(c)],c)
            vals.append(dict(x_pix=s["x_pix"],y_pix=s["y_pix"],quadrant=s["quadrant"],
                n=s["n"],residual_sum_nMgy=float(val)))
        v=np.asarray([a["residual_sum_nMgy"] for a in vals])
        rows.append(dict(radius_arcsec=radius,count=len(vals),
            quadrants=sorted(set(a["quadrant"] for a in vals)),
            median_nMgy=float(np.median(v)) if v.size else None,
            rms_nMgy=float(np.sqrt(np.mean(v*v))) if v.size else None,
            mad_sigma_nMgy=float(1.4826*np.median(abs(v-np.median(v)))) if v.size>=5 else None,
            samples=vals))
    return rows


def candidate(rows,method,samples,tx,ty,scale,required_radius_arcsec=48.):
    q=np.asarray([r["quadrant"] for r in rows],dtype=int)
    counts=[int((q==i).sum()) for i in range(4)]
    out=dict(method=method,quadrant_block_counts=counts,required_target_support_radius_arcsec=float(required_radius_arcsec),
             status="SKY1_INSUFFICIENT_SPATIAL_SUPPORT")
    if min(counts)<3: return out
    f=fit_blocks(rows,method)
    if f is None: return dict(out,status="SKY1_FIT_FAILED")
    xy=np.asarray([[(r["x_pix"]-tx)*scale,(r["y_pix"]-ty)*scale] for r in rows])
    try:
        hull=Delaunay(xy)
        # Reference target region is inside training convex hull, not
        # extrapolated beyond its exterior. Interior interpolation is
        # still an unvalidated assumption across the protected region.
        theta=np.linspace(0,2*np.pi,32,endpoint=False)
        # Circumscribed polygon: its interior contains the entire required disk.
        radius=required_radius_arcsec/np.cos(np.pi/32)
        test=np.vstack(([0.,0.],radius*np.column_stack((np.cos(theta),np.sin(theta)))))
        inside=bool(np.all(hull.find_simplex(test)>=0))
    except QhullError:
        inside=False
    folds={}; cv=[]
    for i in range(4):
        train=[r for r in rows if r["quadrant"]!=i]
        test=[r for r in rows if r["quadrant"]==i]
        ff=fit_blocks(train,method)
        if ff is None: continue
        c=np.asarray(ff["coefficients_nMgy_pix"])
        d=np.asarray([r["mean_nMgy_pix"] for r in test])-np.asarray([r["design"][:len(c)] for r in test])@c
        folds[i]=ff
        cv.append(dict(quadrant=i,training_blocks=len(train),heldout_blocks=len(test),
            median_nMgy_pix=float(np.median(d)),mae_nMgy_pix=float(np.mean(abs(d))),
            rms_nMgy_pix=float(np.sqrt(np.mean(d*d))),fit=ff))
    out.update(fit=f,target_region_inside_training_hull=inside,heldout_quadrants=cv,
        empty_apertures=summarize_apertures(samples,folds))
    if len(cv)==4:
        zero=float(np.mean([np.mean([abs(r["mean_nMgy_pix"]) for r in rows if r["quadrant"]==i]) for i in range(4)]))
        out["zero_equal_quadrant_mae_nMgy_pix"]=zero
        out["cv_equal_quadrant_mae_nMgy_pix"]=float(np.mean([r["mae_nMgy_pix"] for r in cv]))
        out["cv_max_abs_quadrant_median_nMgy_pix"]=max(abs(r["median_nMgy_pix"]) for r in cv)
    out["status"]="SKY1_CONDITIONAL_NOT_ADOPTED" if inside and len(cv)==4 else "SKY1_VALIDATION_OR_HULL_BLOCKED"
    return out


def sky_diagnostics(image,iv,hard,*,source_xy,target_xy,scale,source_mask,protected,veto,
                    neighbours,comparison,preview,wcs):
    a=np.asarray(image,dtype=float); iv=np.asarray(iv,dtype=float)
    hard,source_mask,protected,veto=[np.asarray(m,dtype=bool) for m in (hard,source_mask,protected,veto)]
    if a.ndim!=2 or any(m.shape!=a.shape for m in (iv,hard,source_mask,protected,veto)):
        raise ValueError("SKY-1 expanded grid mismatch")
    if not np.isfinite([*source_xy,*target_xy,scale]).all() or scale<=0:
        raise ValueError("SKY-1 invalid geometry")
    p=preview; tx,ty=target_xy; sx,sy=source_xy
    ox,oy=p["target_xy_native_pix"]; osx,osy=p["source_xy_native_pix"]
    if not np.allclose([sx-tx,sy-ty],[osx-ox,osy-oy],rtol=0,atol=.01):
        raise ValueError("SKY-1 expanded/native relative geometry mismatch")
    safe=safe_pixels(a,iv,hard,source_mask,protected,veto,neighbours,scale)
    yy,xx=np.ogrid[:a.shape[0],:a.shape[1]]
    rt=np.hypot(xx-tx,yy-ty)*scale
    safe &= abs(xx-sx)*scale>=4.
    ref=comparison["candidates"][p["reference_candidate_index"]]
    variants={"SH1C":(ref,p["layers"]["residual"])}
    e=p.get("sh1e") or {}
    for t in e.get("half_width_trials",[]):
        if t["half_width_arcsec"]!=4.: continue
        match=next((c for c in t.get("refit_candidates",[]) if key(c)==key(ref)),None)
        im=e.get("layers",{}).get("spike_4_refit_residual")
        if match is not None and im is not None: variants["REFIT_4"]=(match,im)
    # Shared expanded training support and same original photometric pixels
    # for all halo alternatives. No mirrored-spike model in sky training.
    res={}
    for name,(c,im) in variants.items():
        h,b=model_arrays(a.shape,source_xy,scale,c)
        res[name]=a-h-b
        safe &= np.isfinite(res[name])
    g=local_grid(p,wcs)
    native={name:im[g["sl"]] for name,(c,im) in variants.items()}
    valid=np.logical_and.reduce([np.isfinite(im) for im in native.values()])
    valid &= np.isfinite(p["layers"]["original"][g["sl"]])
    refgeo=(p.get("geo1") or {}).get("reference") or {}
    # A 40 arcsec ellipse around the fitted centre fits inside this disk.
    # Also cover the fixed 35 arcsec geometry fit window about the anchor.
    shift=float(np.hypot(refgeo["east_offset_arcsec"],refgeo["north_offset_arcsec"])) if fitted(refgeo) else 0.
    required=max(48.,40.+shift+1.)
    d=block_design(g,valid,35.,3.)
    out=dict(schema="desidera.sky1-expanded-residual/1",status="SKY1_NOT_ADOPTED",science_adopted=False,
        original_pixels_modified=False,scene_mask_modified=False,expanded_shape=list(a.shape),
        expanded_target_xy=list(target_xy),expanded_source_xy=list(source_xy),
        shared_safe_mask_sha256=native_boolean_mask_sha256(safe),
        shared_original_crop_mask_sha256=native_boolean_mask_sha256(valid),
        shared_original_crop_origin_xy=[int(g["x"][0,0]),int(g["y"][0,0])],
        trials=[],layers={},reference_geometry=refgeo,required_target_support_radius_arcsec=required,
        warning="Residual correction on top of SH1 sky, not upstream sky. Holdout quadrants exclude this fit's training only; prior halo fit may overlap. No missing-pixel filling or calibrated uncertainty; quadratic is a stress scenario. Masked target interior is interpolation, not observed sky.")
    # Expanded display layers, not replacement science images.
    out["expanded_layers"]={"original":np.where(np.isfinite(a)&~hard&(iv>0),a,np.nan)}
    for name,v in res.items():
        out["expanded_layers"][name]=np.where(np.isfinite(a)&~hard&(iv>0),v,np.nan)
    out["expanded_safe_mask"]=safe.copy()
    for v in list(out["expanded_layers"].values())+[out["expanded_safe_mask"]]:
        v.setflags(write=False)
    out["pixel_scale_arcsec"]=float(scale)
    out["native_reference_residual"]=p["layers"]["residual"]
    # Compute profiles on the fixed GEO-1 ellipse before optional geometry refit.
    for guard in GUARDS:
        m=safe&(rt>=guard)
        for name,v in res.items():
            samples=empty_samples(v,m,tx,ty,scale)
            for tile in TILES:
                rows=tile_rows(v,m,tx,ty,scale,tile)
                trial=dict(halo_scenario=name,halo_hypothesis=list(key(variants[name][0])),
                    target_guard_arcsec=guard,tile_size_pix=tile,
                    safe_pixels=int(m.sum()),safe_mask_sha256=native_boolean_mask_sha256(m),
                    blocks=rows,candidates=[])
                for method in METHODS:
                    c=candidate(rows,method,samples,tx,ty,scale,required_radius_arcsec=required)
                    trial["candidates"].append(c)
                    if c["status"]!="SKY1_CONDITIONAL_NOT_ADOPTED": continue
                    # Evaluate sky only in the local geometry crop. Full
                    # display is restricted to the training convex hull.
                    delta=predict(c["fit"],(g["x"]-ox)*scale,(g["y"]-oy)*scale)
                    corrected=native[name]-delta
                    if fitted(refgeo):
                        c["fixed_geometry_profiles"]=radial_profiles(g,
                            {"baseline":native[name],"corrected":corrected},valid,refgeo)
                    # Predeclared reference configuration, not winner selection.
                    if guard==60. and tile==32:
                        c["geometry_fixed_zero_background"]=fit_geometry(d,corrected,g,
                            refgeo.get("robust_loss_scale_nMgy_pix"),refgeo if fitted(refgeo) else None,
                            fixed_background_nMgy_pix=0.)
                        if name=="SH1C" and method in ("constant","plane","quadratic"):
                            fy,fx=np.ogrid[:p["layers"]["original"].shape[0],:p["layers"]["original"].shape[1]]
                            X,Y=np.broadcast_arrays((fx-ox)*scale,(fy-oy)*scale)
                            points=np.asarray([[(r["x_pix"]-tx)*scale,(r["y_pix"]-ty)*scale] for r in rows])
                            hull=Delaunay(points)
                            supported=hull.find_simplex(np.column_stack((X.ravel(),Y.ravel()))).reshape(X.shape)>=0
                            baseline=variants[name][1]
                            sky=np.where(supported&np.isfinite(baseline),predict(c["fit"],X,Y),np.nan)
                            residual=baseline-sky
                            for v2 in (sky,residual): v2.setflags(write=False)
                            out["layers"][f"sky1_{method}_model"]=sky
                            out["layers"][f"sky1_{method}_corrected"]=residual
                if guard==60. and tile==32:
                    # Symmetric edge stress, not a change to the scene.
                    pad=25./scale
                    ny,nx=a.shape
                    for side in ("left","right","bottom","top","all"):
                        def keep(x,y,w,h):
                            return ((side not in ("left","all") or x>=pad)
                                and (side not in ("right","all") or x+w<=nx-pad)
                                and (side not in ("bottom","all") or y>=pad)
                                and (side not in ("top","all") or y+h<=ny-pad))
                        rr=[r for r in rows if keep(r["tile_x_pix"],r["tile_y_pix"],tile,tile)]
                        ss=[]
                        for sample in samples:
                            radius=sample["radius_arcsec"]/scale
                            if keep(sample["x_pix"]-radius,sample["y_pix"]-radius,2*radius+1,2*radius+1):
                                ss.append(sample)
                        trial.setdefault("edge_stress",[]).append(dict(side=side,width_arcsec=25.,
                            retained_blocks=len(rr),removed_blocks=len(rows)-len(rr),
                            candidates=[candidate(rr,method,ss,tx,ty,scale,required_radius_arcsec=required)
                                        for method in METHODS]))
                out["trials"].append(trial)
    return out


def paired_sky_diagnostics(g,r):
    """Matched observed sky pixels, blocks and empty-aperture centres.

    Same-session snapshots only; no new file reads, no PSF matching.
    """
    if g.get("band")!="g" or r.get("band")!="r":
        raise ValueError("SKY-2 pair requires g and r")
    for k in ("project_id","draft_sha256","catalog_sha256","foreground_source_id",
              "source_triplet","original_grid_fingerprint","sh1_policy"):
        if g.get(k) is None or g.get(k)!=r.get(k):
            raise ValueError("SKY-2 pair identity mismatch: "+k)
    gg,rr=g["expanded_sky"],r["expanded_sky"]
    for k in ("expanded_grid_fingerprint","expanded_shape"):
        if gg.get(k) is None or gg.get(k)!=rr.get(k):
            raise ValueError("SKY-2 expanded grid mismatch: "+k)
    for k in ("expanded_target_xy","expanded_source_xy","pixel_scale_arcsec"):
        if not np.allclose(gg[k],rr[k],rtol=0,atol=1e-6):
            raise ValueError("SKY-2 pair geometry mismatch: "+k)
    tx,ty=gg["expanded_target_xy"]; scale=gg["pixel_scale_arcsec"]
    shape=tuple(gg["expanded_shape"])
    images={b:d["expanded_layers"]["SH1C"] for b,d in (("g",gg),("r",rr))}
    masks=[d["expanded_safe_mask"] for d in (gg,rr)]
    if any(m.shape!=shape for m in masks+list(images.values())):
        raise ValueError("SKY-2 pair array shape mismatch")
    y,x=np.ogrid[:shape[0],:shape[1]]
    safe=masks[0]&masks[1]&np.isfinite(images["g"])&np.isfinite(images["r"])
    safe &= np.hypot(x-tx,y-ty)*scale>=60.
    required=max(gg["required_target_support_radius_arcsec"],rr["required_target_support_radius_arcsec"])
    results={}
    for band,image in images.items():
        rows=tile_rows(image,safe,tx,ty,scale,32)
        samples=empty_samples(image,safe,tx,ty,scale)
        results[band]=dict(blocks=rows,candidates=[candidate(rows,m,samples,tx,ty,scale,required)
                                                   for m in METHODS],
                           empty_aperture_centres=[[s["radius_arcsec"],s["x_pix"],s["y_pix"],s["n"]] for s in samples])
    def loc(b): return [(v["tile_x_pix"],v["tile_y_pix"],v["n"]) for v in results[b]["blocks"]]
    if loc("g")!=loc("r") or results["g"]["empty_aperture_centres"]!=results["r"]["empty_aperture_centres"]:
        raise ValueError("SKY-2 paired sampling mismatch")
    return dict(schema="desidera.sky2-matched-controls/1",status="SKY2_PAIRED_NOT_ADOPTED",
        science_adopted=False,psf_matched=False,common_mask_sha256=native_boolean_mask_sha256(safe),
        common_pixels=int(safe.sum()),expanded_grid_fingerprint=gg["expanded_grid_fingerprint"],
        guard_arcsec=60.,tile_size_pix=32,halo_scenario="SH1C",bands=results,
        band_fits_sha256={b:v["original_fits_sha256_local"] for b,v in (("g",g),("r",r))},
        warning="Same-session snapshots, prior-band file hashes not refreshed. Spatial validation conditional on halo fit; no adopted colour, sky or covariance.")
