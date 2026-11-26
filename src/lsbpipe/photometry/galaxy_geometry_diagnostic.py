"""GEO-1: provisional geometry from signed native residuals.

Robust block-mean fits of a generalized exponential plus free constant.
Neither the shape nuisance profile nor its amplitude is science photometry.
Same observed pixels for scenario comparisons; no masked-pixel fill.
"""
import hashlib
import numpy as np
from scipy.optimize import least_squares
from ..sky_geometry import ellipse_from_wcs, ellipse_from_jacobian

RADII=(20.,30.,35.,40.)
PROFILE_EDGES=(0.,5.,10.,15.,20.,25.,30.,35.,40.)
BLOCK=6
MAX_SHIFT=8.


def local_grid(p,wcs):
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    tx,ty=map(float,p["target_xy_native_pix"])
    center=wcs.celestial.pixel_to_world(tx,ty).icrs
    geo=ellipse_from_wcs(wcs,center.ra.deg,center.dec.deg,sky_q=1.,sky_pa_deg=0.)
    k=geo.local_pixel_per_arcsec_en
    shape=p["layers"]["original"].shape
    # Crop fixed region, including every permitted shifted 40 arcsec ellipse.
    margin=int(np.ceil(50.*np.linalg.norm(k,ord=2)))+2
    x0=max(0,int(np.floor(tx))-margin); x1=min(shape[1],int(np.ceil(tx))+margin+1)
    y0=max(0,int(np.floor(ty))-margin); y1=min(shape[0],int(np.ceil(ty))+margin+1)
    yy,xx=np.mgrid[y0:y1,x0:x1]
    en=np.linalg.inv(k)@np.stack([(xx-tx).ravel(),(yy-ty).ravel()])
    return dict(x=xx,y=yy,e=en[0].reshape(xx.shape),n=en[1].reshape(xx.shape),
                k=k,center=center,anchor=(tx,ty),sl=np.s_[y0:y1,x0:x1],
                area=geo.pixel_area_arcsec2,wcs=wcs,shape=shape)


def _hash(mask):
    return hashlib.sha256(np.ascontiguousarray(mask,dtype=np.uint8).tobytes()).hexdigest()


def scenario_images(p):
    s={"SH1C":p["layers"]["residual"]}
    e=p.get("sh1e") or {}
    for width in (2,4):
        key=f"spike_{width}_refit_residual"
        if key in e.get("layers",{}): s[f"REFIT_{width}"]=e["layers"][key]
    decay=e.get("decay_extent") or {}
    if decay.get("status")=="SH1E3_CONDITIONAL_GRID_NOT_ADOPTED" and "decay_spike_corrected" in decay.get("layers",{}):
        s["DECAY_CONDITIONAL"]=decay["layers"]["decay_spike_corrected"]
    return s


def crop_images(p,g,images):
    if any(a.shape!=g["shape"] for a in images.values()):
        raise ValueError("GEO-1 scenario grid mismatch")
    return {k:np.asarray(a[g["sl"]],dtype=float) for k,a in images.items()}


def block_design(g,valid,radius,inner=3.,strip=0.,source_x=None,scale=None):
    rad=np.hypot(g["e"],g["n"])
    domain=(rad>=inner)&(rad<=radius)
    if strip:
        if source_x is None or scale is None: raise ValueError("Missing spike strip geometry")
        domain &= abs(g["x"]-source_x)*scale>=strip
    use=domain&valid
    # Each 6x6 native cell is represented by mean signed intensity and
    # centroid of the SAME observed pixels. No imputation or clipping.
    bx=g["x"]//BLOCK; by=g["y"]//BLOCK
    ids=by*(g["shape"][1]//BLOCK+1)+bx
    _,indices=np.unique(ids,return_inverse=True)
    ids=indices.reshape(ids.shape)
    size=int(ids.max())+1
    nominal=np.bincount(ids[domain],minlength=size)
    observed=np.bincount(ids[use],minlength=size)
    kept=(observed>=8)&(observed>=.5*nominal)
    use &= kept[ids]
    counts=np.bincount(ids[use],minlength=size)
    sel=counts>0
    def mean(a):
        return np.bincount(ids[use],weights=a[use],minlength=size)[sel]/counts[sel]
    angle=np.mod(np.arctan2(g["e"],g["n"]),2*np.pi)
    coverage=[]
    for i in range(8):
        sec=domain&(angle>=i*np.pi/4)&(angle<(i+1)*np.pi/4)
        coverage.append(float((sec&use).sum()/sec.sum()) if sec.any() else 0.)
    return dict(use=use,mean=mean,e=mean(g["e"]),n=mean(g["n"]),
        block_count=int(sel.sum()),used_pixels=int(use.sum()),domain_pixels=int(domain.sum()),
        coverage=coverage,mask_sha256=_hash(use),radius=radius,inner=inner,strip=strip)


def _model(par,e,n):
    ce,cn,q,theta,log_a,log_h,index,bg=par
    de=e-ce; dn=n-cn
    x=de*np.sin(theta)+dn*np.cos(theta)
    y=de*np.cos(theta)-dn*np.sin(theta)
    r=np.hypot(x,y/q)
    return np.exp(log_a)*np.exp(-np.log(2)*(r/np.exp(log_h))**(1/index))+bg


def fit_geometry(design,image,g,noise=None,seed=None):
    d=design
    info={k:d[k] for k in ("radius","inner","strip","block_count","used_pixels","domain_pixels","mask_sha256","coverage")}
    if d["block_count"]<40 or d["used_pixels"]<.5*d["domain_pixels"] or min(d["coverage"])<.25:
        return dict(info,status="GEO1_INSUFFICIENT_COVERAGE")
    z=d["mean"](image)
    if not np.isfinite(z).all(): raise ValueError("Nonfinite GEO-1 common block")
    spread=max(float(np.quantile(z,.9)-np.quantile(z,.1)),1e-12)
    norm=spread
    z=z/norm
    # Shared per-comparison noise scale is an empirical robust loss scale,
    # not a calibrated uncertainty. Default uses robust outer block scatter.
    if noise is None:
        outer=np.hypot(d["e"],d["n"])>=.7*d["radius"]
        v=z[outer] if outer.any() else z
        noise=max(.02,float(1.4826*np.median(abs(v-np.median(v)))))
    else:
        noise=max(.02,float(noise)/norm)
    e,n=d["e"],d["n"]
    bg=float(np.median(z[np.hypot(e,n)>=.7*d["radius"]]))
    amp=max(.1,float(np.quantile(z,.95)-bg))
    lower=np.array([-MAX_SHIFT,-MAX_SHIFT,.25,-2*np.pi,-8.,np.log(1.5),.35,float(z.min()-5)])
    upper=np.array([MAX_SHIFT,MAX_SHIFT,1.,3*np.pi,5.,np.log(70.),2.5,float(z.max()+5)])
    starts=[np.array([0.,0.,.65,np.deg2rad(pa),np.log(amp),np.log(12.),1.,bg]) for pa in (0,60,120)]
    if seed is not None:
        starts=[np.array([seed["east_offset_arcsec"],seed["north_offset_arcsec"],seed["q"],
            np.deg2rad(seed["pa_deg"]),np.log(max(seed["nuisance_amplitude_nMgy_pix"]/norm,1e-8)),
            np.log(seed["nuisance_half_intensity_radius_arcsec"]),seed["nuisance_index"],
            seed["nuisance_constant_nMgy_pix"]/norm])]+starts[:1]
    candidates=[]
    for start in starts:
        f=least_squares(lambda par:_model(par,e,n)-z,np.clip(start,lower+1e-8,upper-1e-8),
            bounds=(lower,upper),loss="soft_l1",f_scale=noise,x_scale="jac",max_nfev=120)
        if f.success and np.isfinite(f.x).all() and np.isfinite(f.cost):
            candidates.append(f)
    if not candidates: return dict(info,status="GEO1_OPTIMIZER_UNCONVERGED")
    f=min(candidates,key=lambda t:t.cost)
    ce,cn,q,theta,la,lh,index,bg=f.x
    pa=float(np.degrees(theta)%180)
    dx,dy=g["k"]@np.array([ce,cn])
    x,y=g["anchor"][0]+dx,g["anchor"][1]+dy
    sky=g["wcs"].celestial.pixel_to_world(x,y).icrs
    null_cost=float(np.sum(noise**2*(np.sqrt(1+((z-np.median(z))/noise)**2)-1)))
    improvement=1-f.cost/null_cost if null_cost>0 else 0.
    flags=[]
    if q<=.26: flags.append("Q_LOWER_BOUND")
    if q>=.95: flags.append("PA_WEAK_NEAR_CIRCULAR")
    if max(abs(ce),abs(cn))>=MAX_SHIFT-.1: flags.append("CENTRE_BOUND")
    if index<=.36 or index>=2.49 or np.exp(lh)>=69. or np.exp(lh)<=1.52: flags.append("PROFILE_BOUND")
    if improvement<=.05: flags.append("WEAK_PROFILE_CONTRAST")
    if len(candidates)<len(starts): flags.append("SOME_STARTS_UNCONVERGED")
    # Close-cost multistart solutions can reveal shape degeneracy.
    for alt in candidates:
        if alt.cost<=f.cost*1.01+1e-10:
            da=abs((np.degrees(alt.x[3])-pa+90)%180-90)
            if abs(alt.x[2]-q)>.1 or np.hypot(alt.x[0]-ce,alt.x[1]-cn)>2 or (q<.9 and da>15):
                flags.append("MULTISTART_GEOMETRY_DISAGREEMENT");break
    return dict(info,status="GEO1_CONDITIONAL_FIT_NOT_ADOPTED",
        east_offset_arcsec=float(ce),north_offset_arcsec=float(cn),centre_xy_native_pix=[float(x),float(y)],
        centre_icrs_deg=[float(sky.ra.deg),float(sky.dec.deg)],q=float(q),pa_deg=pa,
        nuisance_amplitude_nMgy_pix=float(np.exp(la)*norm),
        nuisance_half_intensity_radius_arcsec=float(np.exp(lh)),nuisance_index=float(index),
        nuisance_constant_nMgy_pix=float(bg*norm),robust_loss_scale_nMgy_pix=float(noise*norm),
        relative_robust_loss_improvement=float(improvement),converged_starts=len(candidates),
        flags=flags)


def fitted(f):
    return f.get("status")=="GEO1_CONDITIONAL_FIT_NOT_ADOPTED"


def radial_profiles(g,images,valid,f):
    de=g["e"]-f["east_offset_arcsec"]; dn=g["n"]-f["north_offset_arcsec"]
    a=np.deg2rad(f["pa_deg"])
    x=de*np.sin(a)+dn*np.cos(a); y=de*np.cos(a)-dn*np.sin(a)
    r=np.hypot(x,y/f["q"])
    rows=[]
    # Half-planes along fitted major/minor axes, plus total. Native signed
    # sums, with no missing-area or PSF correction.
    for lo,hi in zip(PROFILE_EDGES[:-1],PROFILE_EDGES[1:]):
        for sector,sel in (("ALL",True),("MAJOR_PLUS",x>=0),("MAJOR_MINUS",x<0),
                           ("MINOR_PLUS",y>=0),("MINOR_MINUS",y<0),("N",dn>=0),("S",dn<0)):
            domain=(r>=lo)&(r<hi)&sel
            use=domain&valid; n=int(use.sum())
            ap=(r<hi)&sel&valid
            row=dict(inner_arcsec=lo,outer_arcsec=hi,sector=sector,
                used_pixels=n,in_frame_pixels=int(domain.sum()),cumulative_used_pixels=int(ap.sum()),scenarios={})
            for name,im in images.items():
                flux=float(im[use].sum()) if n else None
                row["scenarios"][name]=dict(ring_flux_nMgy=flux,
                    mean_nMgy_arcsec2=flux/(n*g["area"]) if n else None,
                    cumulative_flux_nMgy=float(im[ap].sum()) if ap.any() else None)
            rows.append(row)
    return rows


def overlay_points(g,f,radius):
    geo=ellipse_from_jacobian(*f["centre_xy_native_pix"],f["q"],f["pa_deg"],g["k"])
    z=np.linspace(0,2*np.pi,73)
    a=radius*geo.pixel_a_per_sky_arcsec; b=a*geo.pixel.q
    c,s=np.cos(geo.pixel.theta_rad),np.sin(geo.pixel.theta_rad)
    return np.column_stack([geo.pixel.x_pix+a*np.cos(z)*c-b*np.sin(z)*s,
                           geo.pixel.y_pix+a*np.cos(z)*s+b*np.sin(z)*c]).tolist()


def geometry_diagnostics(p,comparison,wcs):
    g=local_grid(p,wcs)
    images=crop_images(p,g,scenario_images(p))
    valid=np.logical_and.reduce([np.isfinite(v) for v in images.values()])
    valid &= np.isfinite(p["layers"]["original"][g["sl"]])
    scale=p["pixel_scale_arcsec"]; sx,sy=p["source_xy_native_pix"]
    trials=[]; reference=None; reference_design=None
    configs=[(r,3.,0.) for r in RADII]+[(35.,8.,0.),(35.,3.,2.),(35.,3.,4.)]
    for rad,inner,strip in configs:
        d=block_design(g,valid,rad,inner,strip,sx,scale)
        b=fit_geometry(d,images["SH1C"],g)
        if (rad,inner,strip)==(35.,3.,0.):
            reference=b; reference_design=d
        fits={"SH1C":b}
        for name,im in images.items():
            if name!="SH1C":
                fits[name]=fit_geometry(d,im,g,b.get("robust_loss_scale_nMgy_pix"),b if fitted(b) else None)
        trials.append(dict(radius_arcsec=rad,inner_excluded_arcsec=inner,strip_half_width_arcsec=strip,fits=fits))
    ensemble=[]; profiles=[]; circles=[]; overlay=[]
    if reference is not None and fitted(reference):
        d=reference_design; raw=p["layers"]["original"][g["sl"]]
        dx=(g["x"]-sx)*scale; dy=(g["y"]-sy)*scale; rr=np.hypot(dx,dy)
        for i,c in enumerate(comparison["candidates"]):
            h=np.interp(rr,c["radial_knots_arcsec"],c["halo_at_knots_nMgy_pix"],left=np.nan,right=0.)
            b=c["sky_constant_nMgy_pix"]+dx*c["sky_x_nMgy_pix_per_arcsec"]+dy*c["sky_y_nMgy_pix_per_arcsec"]
            im=raw-h-b
            if not np.isfinite(im[d["use"]]).all():
                ensemble.append(dict(index=i,status="GEO1_MODEL_SUPPORT_MISMATCH"));continue
            f=fit_geometry(d,im,g,reference["robust_loss_scale_nMgy_pix"],reference)
            ensemble.append(dict(index=i,**f))
        profiles=radial_profiles(g,images,valid,reference)
        circles=radial_profiles(g,images,valid,dict(reference,q=1.,pa_deg=0.))
        overlay=overlay_points(g,reference,35.)
    return dict(schema="desidera.geo1-native-geometry/1",status="GEO1_NOT_ADOPTED",
        science_adopted=False,psf_matched=False,geometry_method="ROBUST_BLOCK_MEAN_GENERALIZED_EXPONENTIAL_PLUS_CONSTANT",
        angle_convention="LOCAL_ICRS_NORTH_THROUGH_EAST_AT_ANCHOR",block_size_pix=BLOCK,
        common_crop_shape=list(valid.shape),common_crop_origin_xy=[int(g["x"][0,0]),int(g["y"][0,0])],
        common_mask_sha256=_hash(valid),common_pixels=int(valid.sum()),scenario_names=list(images),
        trials=trials,reference_rule="SH1C_35ARCSEC_3ARCSEC_CORE_EXCLUSION_NO_STRIP",
        reference=reference,ensemble_geometry=ensemble,elliptical_profiles=profiles,
        circular_profiles_same_centre=circles,overlay_35arcsec_native_xy=overlay,
        warning="Geometry nuisance fit only; no PSF convolution, no calibrated sky or uncertainty. Cumulative fit windows are not local isophotes. Profiles retain signed observed pixels only; no area fill or total flux. No scene/FITS changes.")


def paired_geometry(g_result,r_result,wcs):
    # Caller must already have a successful paired_sector_diagnostics result;
    # repeat identity guards so this helper cannot bypass them.
    for k in ("project_id","source_id","tractor_catalog_sha256","draft_digest",
              "original_scene_grid_fingerprint","sh1d_policy"):
        if g_result.get(k) is None or g_result.get(k)!=r_result.get(k):
            raise ValueError("GEO-1 pair mismatch: "+k)
    if g_result["band"]!="g" or r_result["band"]!="r": raise ValueError("GEO-1 requires g/r")
    gp,rp=g_result["sh1c_preview"],r_result["sh1c_preview"]
    for k in ("target_xy_native_pix","source_xy_native_pix","pixel_scale_arcsec"):
        if not np.allclose(gp[k],rp[k],rtol=0,atol=1e-6): raise ValueError("GEO-1 pair geometry mismatch: "+k)
    g=local_grid(gp,wcs)
    full={band+"_"+name:im for band,p in (("g",gp),("r",rp)) for name,im in scenario_images(p).items()}
    images=crop_images(gp,g,full)
    valid=np.logical_and.reduce([np.isfinite(v) for v in images.values()])
    for p in (gp,rp):
        if p["layers"]["original"].shape!=g["shape"]: raise ValueError("GEO-1 pair shape mismatch")
        valid &= np.isfinite(p["layers"]["original"][g["sl"]])
    d=block_design(g,valid,35.,3.)
    fits={band:fit_geometry(d,images[band+"_SH1C"],g) for band in ("g","r")}
    profiles={band:radial_profiles(g,images,valid,f) for band,f in fits.items() if fitted(f)}
    return dict(schema="desidera.geo1-paired-native/1",status="GEO1_PAIRED_NOT_ADOPTED",
        science_adopted=False,psf_matched=False,common_mask_sha256=_hash(valid),common_pixels=int(valid.sum()),
        common_crop_shape=list(valid.shape),common_crop_origin_xy=[int(g["x"][0,0]),int(g["y"][0,0])],
        grid_fingerprint=g_result["original_scene_grid_fingerprint"],draft_digest=g_result["draft_digest"],
        fits=fits,profiles_by_geometry_band=profiles,
        warning="Same-session snapshots, no fresh rehash of earlier band. All scenario profiles share the same observed mask for each geometry. No PSF matching or adopted colour.")
