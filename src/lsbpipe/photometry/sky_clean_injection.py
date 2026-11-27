"""SKY-CLEAN / INJ-1: conditional sky-only source screening and galaxy recovery.

Never modifies the Scene photometric mask, original FITS or saved draft.
Uses fixed SH1C residuals; does NOT validate the upstream halo/PSF fit.
All models in a comparison use identical pixels. Cross-validation after
image-based screening is conditional on that screening.
"""
from types import SimpleNamespace
import numpy as np
from scipy.ndimage import gaussian_filter, label, distance_transform_edt
from astropy.wcs import WCS
from .psf_tractor_gaia import read_native_tractor_gaia, SATURATION_BITS
from .sky_masks import native_boolean_mask_sha256
from .residual_sky_diagnostic import tile_rows, empty_samples, candidate, fit_blocks, predict, METHODS
from .galaxy_geometry_diagnostic import local_grid, block_design, fit_geometry, fitted

SAT_RADIUS=20.
PEAK_SIGMA=8.
PEAK_GROW=4.
GUARD=60.
TILE=32
# Observed exponential profiles, not intrinsic PSF-convolved models.
SHAPES=((15.,24.5,.56,44.5),(25.,25.5,.70,75.))
INJ_METHODS=("zero","constant","plane","quadratic")


def mad(v):
    v=np.asarray(v)
    return float(1.4826*np.median(abs(v-np.median(v)))) if v.size else 0.


def saturations(shape,points,rows,band,scale, radius=SAT_RADIUS):
    """Use saturation bit, not flux or IVAR, to veto sky training."""
    if band not in SATURATION_BITS: raise ValueError("Unsupported band")
    mask=np.zeros(shape,bool)
    sources=[]
    ny,nx=shape
    r=radius/scale
    bit=1<<SATURATION_BITS[band][0]
    for (x,y),row in zip(points,rows):
        if not row["maskbits"]&bit: continue
        if not np.isfinite([x,y]).all(): continue
        x0=max(0,int(np.floor(x-r))); x1=min(nx,int(np.ceil(x+r))+1)
        y0=max(0,int(np.floor(y-r))); y1=min(ny,int(np.ceil(y+r))+1)
        if x1<=x0 or y1<=y0: continue
        yy,xx=np.ogrid[y0:y1,x0:x1]
        disk=(xx-x)**2+(yy-y)**2<=r*r
        if not disk.any(): continue
        mask[y0:y1,x0:x1]|=disk
        sources.append(dict(objid=int(row["objid"]),x_pix=float(x),y_pix=float(y),
            brick_primary=bool(row.get("brick_primary",True)),
            flux_nMgy=row.get("flux"),flux_ivar=row.get("flux_ivar"),
            radius_arcsec=radius))
    return mask,sources


def compact_peaks(res,safe,scale):
    """Positive compact excess relative to 10-pixel smoothed background.

    Normalized convolutions avoid zero-filled gaps entering the estimate.
    Negative components provide a diagnostic control, not a flux correction.
    Detection operates only on pre-existing eligible off-target pixels.
    """
    m=safe&np.isfinite(res)
    if m.sum()<100: return np.zeros_like(m),dict(status="INSUFFICIENT_DETECTION_SUPPORT")
    z=np.where(m,res,0.)
    def smooth(s):
        den=gaussian_filter(m.astype(float),s,mode="constant",cval=0.)
        num=gaussian_filter(z,s,mode="constant",cval=0.)
        return np.divide(num,den,out=np.zeros_like(num),where=den>0),den
    small,ds=smooth(1.5); large,dl=smooth(10.)
    use=m&(ds>=.75)&(dl>=.5)
    v=small-large
    if use.sum()<100:
        return np.zeros_like(m),dict(status="INSUFFICIENT_DETECTION_SUPPORT")
    med=float(np.median(v[use])); noise=mad(v[use])
    if noise<=0 or not np.isfinite(noise):
        return np.zeros_like(m),dict(status="NO_MEASURABLE_DETECTION_SCATTER")
    def components(which):
        lab,n=label(which,structure=np.ones((3,3),int))
        count=np.bincount(lab.ravel(),minlength=n+1)
        keep=count>=3; keep[0]=False
        return keep[lab],int(keep.sum())
    pos,npk=components(use&(v-med>PEAK_SIGMA*noise))
    neg,nn=components(use&(v-med< -PEAK_SIGMA*noise))
    grown=(distance_transform_edt(~pos)<=PEAK_GROW/scale) if pos.any() else pos
    return grown,dict(status="CONDITIONAL_COMPACT_EXCESS",threshold_sigma=PEAK_SIGMA,
        highpass_median_nMgy_pix=med,highpass_mad_nMgy_pix=noise,
        positive_components=npk,negative_components=nn,
        positive_core_pixels=int(pos.sum()),negative_core_pixels=int(neg.sum()),
        grown_pixels=int(grown.sum()),grow_arcsec=PEAK_GROW)


def quality_controls(res,mask,tx,ty,scale,required):
    rows=tile_rows(res,mask,tx,ty,scale,TILE)
    samples=empty_samples(res,mask,tx,ty,scale)
    # Do not remove outliers in the main fit. Show an explicitly
    # data-selected, two-sided block-heterogeneity stress experiment.
    d=np.asarray([r["mean_nMgy_pix"]-r["median_nMgy_pix"] for r in rows])
    loc=float(np.median(d)) if d.size else 0.
    sig=mad(d)
    stress=[r for r,v in zip(rows,d) if abs(v-loc)<=6*sig] if sig>0 else rows
    main=[candidate(rows,m,samples,tx,ty,scale,required) for m in METHODS]
    for c in main:
        for a in c.get("empty_apertures",[]):
            v=np.asarray([s["residual_sum_nMgy"] for s in a["samples"]])
            a["max_abs_nMgy"]=float(np.max(abs(v))) if v.size else None
            a["tail_beyond_5mad_count"]=int((abs(v-np.median(v))>5*mad(v)).sum()) if v.size>=10 and mad(v)>0 else None
    return dict(block_count=len(rows),empty_aperture_count=len(samples),
        tile_ids=[[r["tile_x_pix"],r["tile_y_pix"],r["n"]] for r in rows],
        empty_aperture_centres=[[s["radius_arcsec"],s["x_pix"],s["y_pix"],s["n"]] for s in samples],
        candidates=main,mean_median_stress=dict(removed_blocks=len(rows)-len(stress),
            centre_nMgy_pix=loc,mad_nMgy_pix=sig,
            candidates=[candidate(stress,m,[],tx,ty,scale,required) for m in METHODS],
            warning="Two-sided, data-selected stress; CV is conditional on retained blocks; no noise calibration."))


def select_sites(mask,tx,ty,scale,max_sites=2):
    """Deterministic coverage-only selection, separated from real target."""
    ny,nx=mask.shape
    pad=int(np.ceil(50./scale)); r=int(np.ceil(40./scale))
    step=max(1,int(round(20./scale)))
    yy,xx=np.ogrid[-r:r+1,-r:r+1]
    disk=(xx*xx+yy*yy)*scale**2<=40.**2
    sites=[]
    for y in range(pad,ny-pad,step):
        for x in range(pad,nx-pad,step):
            if np.hypot(x-tx,y-ty)*scale<100.: continue
            cover=float(mask[y-r:y+r+1,x-r:x+r+1][disk].mean())
            if cover>=.75: sites.append((cover,x,y))
    chosen=[]
    for cover,x,y in sorted(sites,key=lambda z:(-z[0],z[2],z[1])):
        if all(np.hypot(x-s["x_pix"],y-s["y_pix"])*scale>=90. for s in chosen):
            chosen.append(dict(x_pix=x,y_pix=y,coverage_r40=cover))
        if len(chosen)>=max_sites: break
    return chosen


def galaxy_image(shape,wcs,site,params,scale):
    """Evaluate an observed exponential on native pixel centres.

    The local tangent-plane Jacobian uses the same ICRS convention as GEO-1.
    Extend the exponential across the entire available field, including
    its wings in the sky-training domain. No re-normalization at edges.
    """
    re,mu,q,pa=params
    x,y=site["x_pix"],site["y_pix"]
    from ..sky_geometry import ellipse_from_wcs
    sky=wcs.celestial.pixel_to_world(x,y).icrs
    geo=ellipse_from_wcs(wcs,sky.ra.deg,sky.dec.deg,sky_q=1.,sky_pa_deg=0.)
    k=np.linalg.inv(geo.local_pixel_per_arcsec_en)
    yy,xx=np.ogrid[:shape[0],:shape[1]]
    e=k[0,0]*(xx-x)+k[0,1]*(yy-y)
    n=k[1,0]*(xx-x)+k[1,1]*(yy-y)
    t=np.deg2rad(pa)
    rho=np.hypot(e*np.sin(t)+n*np.cos(t),(e*np.cos(t)-n*np.sin(t))/q)
    amp=geo.pixel_area_arcsec2*10**(.4*(22.5-mu))
    return amp*np.exp(-1.67834699001666*rho/re),rho


def inject(res,mask,wcs,site,params,tx,ty,scale,required):
    x,y=site["x_pix"],site["y_pix"]
    model,rho=galaxy_image(res.shape,wcs,site,params,scale)
    yy,xx=np.ogrid[:res.shape[0],:res.shape[1]]
    base=mask&(np.hypot(xx-tx,yy-ty)*scale>=GUARD)
    # Refit on exactly the same pixels with and without injection.
    # Exclude 45 arcsec around the injected centre; the exponential wings
    # beyond this radius remain in the injected training image.
    train=base&(np.hypot(xx-x,yy-y)*scale>=45.)
    injected=res+model
    before=tile_rows(res,train,tx,ty,scale,TILE)
    after=tile_rows(injected,train,tx,ty,scale,TILE)
    p=dict(target_xy_native_pix=[x,y],layers={"original":res})
    grid=local_grid(p,wcs)
    valid=mask[grid["sl"]]
    d=block_design(grid,valid,35.,3.)
    area=grid["area"]
    re,mu,q,pa=params
    result=dict(site=site,truth=dict(re_arcsec=re,mu0_mag_arcsec2=mu,q=q,pa_deg=pa,
        profile="observed exponential, sampled at native pixel centres",
        analytic_total_nMgy=float(2*np.pi*q*(re/1.67834699001666)**2*10**(.4*(22.5-mu)))),
        sky_training_blocks=len(before),sky_training_mask_sha256=native_boolean_mask_sha256(train),
        models={})
    # Check the actual injection aperture is within the training hull.
    # The helper tests a disk around its input centre; fit rows have
    # coordinates/basis relative to the REAL target, so pass the real
    # centre and check the displaced region separately.
    from scipy.spatial import Delaunay,QhullError
    points=np.asarray([[r["x_pix"],r["y_pix"]] for r in before])
    angle=np.linspace(0,2*np.pi,32,endpoint=False)
    test=np.column_stack((x+40.5/scale/np.cos(np.pi/32)*np.cos(angle),
                          y+40.5/scale/np.cos(np.pi/32)*np.sin(angle)))
    try: inside=bool(np.all(Delaunay(points).find_simplex(test)>=0))
    except (QhullError,ValueError): inside=False
    result["injection_aperture_inside_training_hull"]=inside
    dx=(grid["x"]-tx)*scale; dy=(grid["y"]-ty)*scale
    for method in INJ_METHODS:
        f0=f1=None
        if method!="zero":
            if not inside:
                result["models"][method]=dict(status="INJ1_SKY_HULL_BLOCKED");continue
            f0=fit_blocks(before,method); f1=fit_blocks(after,method)
            if f0 is None or f1 is None:
                result["models"][method]=dict(status="INJ1_SKY_FIT_BLOCKED");continue
        bg0=predict(f0,dx,dy) if f0 else 0.
        bg1=predict(f1,dx,dy) if f1 else 0.
        control=res[grid["sl"]]-bg0
        recovered=injected[grid["sl"]]-bg1
        geo=fit_geometry(d,recovered,grid,fixed_background_nMgy_pix=0.)
        record=dict(status="INJ1_CONDITIONAL_NOT_ADOPTED",geometry=geo,sky_before=f0,sky_after=f1,apertures=[])
        if fitted(geo):
            record["geometry_errors"]=dict(centre_arcsec=float(np.hypot(geo["east_offset_arcsec"],geo["north_offset_arcsec"])),
                delta_q=float(geo["q"]-q),delta_pa_deg=float((geo["pa_deg"]-pa+90)%180-90))
        for rad in (20.,35.,40.):
            domain=rho[grid["sl"]]<=rad
            use=domain&valid
            true=float(model[grid["sl"]][use].sum())
            if true<=0: continue
            obs=float(recovered[use].sum()); ctrl=float(control[use].sum())
            record["apertures"].append(dict(a_arcsec=rad,pixels=int(use.sum()),
                coverage=float(use.sum()/domain.sum()),true_observed_flux_nMgy=true,
                absolute_recovered_flux_nMgy=obs,control_flux_nMgy=ctrl,
                differential_recovered_flux_nMgy=obs-ctrl,
                absolute_fractional_error=obs/true-1.,
                differential_fractional_error=(obs-ctrl)/true-1.))
        result["models"][method]=record
    return result


def expanded_wcs(probe,ds,sky):
    meta=probe["sh1b2_native_window"]
    if meta["grid_fingerprint"]!=sky["expanded_grid_fingerprint"] or meta["native_shape"]!=sky["expanded_shape"]:
        raise ValueError("INJ1 expanded grid identity mismatch")
    w=WCS(meta["wcs_fits_header"])
    x0,y0=meta["actual_bounds_xyxy_exclusive"][:2]
    ox,oy=meta["scene_origin_full_brick_xy"]
    ny,nx=ds["pixel_shape"]
    points=np.array([[0.,0.],[nx-1.,0.],[0.,ny-1.],[nx-1.,ny-1.],[(nx-1)/2,(ny-1)/2]])
    shifted=w.all_world2pix(ds["wcs"].all_pix2world(points,0),0)
    if not np.isfinite(shifted).all() or np.max(abs(shifted-points-[ox-x0,oy-y0]))>1e-3:
        raise ValueError("INJ1 expanded WCS translation mismatch")
    return w


def run_validation(project,products,datasets,paired):
    """Independent downstream diagnostic, not a science-ready photometry stage."""
    out=dict(schema="desidera.sky-clean-inj1/1",status="INJ1_NOT_ADOPTED",
        science_adopted=False,bands={},masks={},sites=[],
        warning="Fixed SH1 residuals and source masks; no upstream halo refit, PSF convolution, injected Poisson noise, mask re-detection after injection, or calibrated uncertainties. Absolute recovery includes actual sky residuals; differential recovery isolates changes in the sky fit on identical pixels.")
    if len(products)>1 and (set(products)!={"g","r"} or paired.get("sky",{}).get("status")!="SKY2_PAIRED_NOT_ADOPTED"):
        return dict(out,status="INJ1_PAIR_IDENTITY_BLOCKED")
    images={}; masks={}; cleaned={}; grids={}; geom=None
    for band,b in products.items():
        pr=b["probe"]; s=b["sky"]["expanded_sky"]; ds=datasets[band]
        if "sh1b2_native_window" not in pr or "expanded_safe_mask" not in s:
            return dict(out,status="INJ1_MISSING_EXPANDED_EVIDENCE")
        w=expanded_wcs(pr,ds,s)
        res=s["expanded_layers"]["SH1C"]; safe=np.asarray(s["expanded_safe_mask"],bool)&np.isfinite(res)
        tx,ty=s["expanded_target_xy"]; scale=s["pixel_scale_arcsec"]
        ident=(tuple(res.shape),s["expanded_grid_fingerprint"],tx,ty,scale)
        if geom is not None and ident!=geom: raise ValueError("INJ1 matched expanded grid mismatch")
        geom=ident; grids[band]=w; images[band]=res
        release,region,brick,_=ds["source_quartet_key"]
        catalog=read_native_tractor_gaia(project.data_root,release=release,region=region,brick=brick,
            frame=SimpleNamespace(dataset_id="/".join(ds["source_quartet_key"]),band=band,wcs=w),include_nonprimary=True)
        if catalog.source_sha256!=b["scene"]["catalog_sha256"]: raise ValueError("INJ1 catalogue hash mismatch")
        sat,sources=saturations(res.shape,catalog.pixel_positions,catalog.rows,band,scale)
        yy,xx=np.ogrid[:res.shape[0],:res.shape[1]]
        outside=np.hypot(xx-tx,yy-ty)*scale>=GUARD
        peaks,audit=compact_peaks(res,safe&~sat&outside,scale)
        masks[band]=safe
        cleaned[band]=safe&~sat&~peaks
        out["masks"][band]=dict(saturation=(safe&sat),compact=(safe&~sat&peaks))
        out["bands"][band]=dict(saturated_sources=sources,compact_detection=audit,
            saturation_removed_pixels=int((safe&sat&outside).sum()),
            compact_removed_pixels=int((safe&~sat&peaks&outside).sum()))
    _,_,tx,ty,scale=geom
    common0=np.logical_and.reduce(list(masks.values()))
    common=np.logical_and.reduce(list(cleaned.values()))
    yy,xx=np.ogrid[:common.shape[0],:common.shape[1]]
    guard=np.hypot(xx-tx,yy-ty)*scale>=GUARD
    out["common_before_pixels"]=int((common0&guard).sum())
    out["common_after_pixels"]=int((common&guard).sum())
    out["common_mask_sha256"]=native_boolean_mask_sha256(common&guard)
    out["masks"]["common_removed"]=(common0&~common&guard)
    required=max(b["sky"]["expanded_sky"]["required_target_support_radius_arcsec"] for b in products.values())
    sites=select_sites(common,tx,ty,scale)
    out["sites"]=sites
    for band,res in images.items():
        z=out["bands"][band]
        z["before"]=quality_controls(res,common0&guard,tx,ty,scale,required)
        z["after"]=quality_controls(res,common&guard,tx,ty,scale,required)
        z["injections"]=[inject(res,common,grids[band],site,params,tx,ty,scale,required)
                         for site in sites for params in SHAPES]
    if len(images)==2:
        for key in ("before","after"):
            g,r=(out["bands"][b][key] for b in ("g","r"))
            if g["tile_ids"]!=r["tile_ids"] or g["empty_aperture_centres"]!=r["empty_aperture_centres"]:
                raise ValueError("INJ1 common sampling mismatch")
    if not sites: out["status"]="INJ1_NO_SUITABLE_BLANK_FIELD"
    return out
