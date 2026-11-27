"""INJ-2: target-scale finite-core PSF/Sersic recovery and conditional totals.

No adopted photometry, no original FITS or Scene edits. Signed observed
aperture flux stays separate from analytic model integrals.
"""
import numpy as np
from scipy.optimize import least_squares
from scipy.special import gammaincinv, gammainc, gammaln
from scipy.signal import fftconvolve
from scipy.spatial import Delaunay, QhullError
from .sky_clean_injection import expanded_wcs, GUARD, TILE
from .galaxy_geometry_diagnostic import local_grid, block_design, fit_geometry, fitted
from .residual_sky_diagnostic import tile_rows, fit_blocks, predict
from .stellar_halo import measure_saturated_halo_sectors
from .stellar_halo_spikes import key, model_arrays
from .sky_masks import native_boolean_mask_sha256
from ..sky_geometry import ellipse_from_wcs

METHODS=("zero","constant","plane","quadratic")
# Target-informed stress grid, not fitted truths. Ratios refer to
# measured target F40, not its (unknown) asymptotic magnitude.
SHAPES=((23.,.55,.56,44.5,.8),(30.,.85,.56,44.5,1.),
        (38.,1.15,.70,75.,1.2))
RADII=(5.,10.,20.,30.,35.,40.)


def psf_evidence(project,products,datasets):
    """Verified saved native PSF0 only; no invented seeing kernel."""
    if not all("original_fits_sha256_local" in datasets[b] and
               "source_quartet_key" in datasets[b] for b in products):
        return {},{"status":"PSF_SOURCE_EVIDENCE_UNAVAILABLE","bands":{}}
    from .psf_gui_backend import read_saved_psf0_bands
    restored=read_saved_psf0_bands(project.root)
    kernels={}; meta={}
    for b in products:
        p=restored["restored"].get(b)
        if not p:
            meta[b]={"status":"SAVED_PSF0_UNAVAILABLE"};continue
        r=p["result"]
        a=np.asarray(r.kernel,dtype=float)
        if (p["source_hashes"]!=datasets[b]["original_fits_sha256_local"] or
                tuple(p["dataset"])!=tuple(datasets[b]["source_quartet_key"]) or
                r.recipe.psf_oversampling!=1 or
                a.ndim!=2 or any(v%2!=1 for v in a.shape) or
                not np.isfinite(a).all() or a.sum()<=0. or min(a.shape)<3):
            meta[b]={"status":"SAVED_PSF0_IDENTITY_OR_KERNEL_BLOCKED"};continue
        kernels[b]=a/a.sum()
        report=r.report()
        meta[b]=dict(status="FINITE_CORE_PSF0_NOT_ADOPTED",lineage_id=r.lineage_id,
            kernel_sha256=report["kernel_sha256"],kernel_shape=list(a.shape),
            input_sum=float(a.sum()),negative_mass_fraction=float(-np.minimum(kernels[b],0).sum()),
            wings_validated=False,spatial_representativity_validated=False,
            source_hashes=p["source_hashes"])
    return kernels,dict(status="FINITE_CORE_PSF0_DIAGNOSTIC",bands=meta,
        warnings=restored["warnings"])


def total_flux(i0,re,n,q):
    """Sersic central intensity I0 in nMgy/arcsec2; a_e in arcsec."""
    bn=gammaincinv(2*n,.5)
    return float(2*np.pi*q*i0*re**2*n*np.exp(gammaln(2*n)-2*n*np.log(bn)))


def ellipse_radius(e,north,q,pa):
    t=np.deg2rad(pa)
    return np.hypot(e*np.sin(t)+north*np.cos(t),
                    (e*np.cos(t)-north*np.sin(t))/q)


def radius_pixels(shape,centre,k,q,pa,origin=(0,0),pad=(0,0)):
    """Local ICRS tangent-plane approximation, including convolution pad."""
    py,px=pad
    yy,xx=np.ogrid[-py:shape[0]+py,-px:shape[1]+px]
    x=xx+origin[0]-centre[0]; y=yy+origin[1]-centre[1]
    inv=np.linalg.inv(k)
    return ellipse_radius(inv[0,0]*x+inv[0,1]*y,inv[1,0]*x+inv[1,1]*y,q,pa)


def convolved_profile(rho,re,n,kernel):
    # rho is padded by half the finite kernel support, so valid convolution
    # does not impose a zero boundary inside the returned footprint.
    bn=gammaincinv(2*n,.5)
    image=np.exp(-bn*(rho/re)**(1/n))
    if kernel is not None:
        image=fftconvolve(image,kernel,mode="valid")
    return image


def strict_convolve(image,valid,kernel):
    """No interpolation across masked samples; reject the full kernel box."""
    if image.shape!=valid.shape or not np.isfinite(image[valid]).all():
        raise ValueError("Invalid convolution support")
    support=fftconvolve(valid.astype(float),np.ones(kernel.shape),mode="same")
    good=support>=kernel.size-.25
    out=fftconvolve(np.where(valid,image,0.),kernel,mode="same")
    return np.where(good,out,np.nan),good


def fit_profile(image,valid,grid,geom,kernel,outer=40.,inner=3.):
    """Fit signed annular means and forward-model the same observed pixels.

    Equal ring weights, robust loss, no diagonal/covariance likelihood claim.
    Fitted intensity/total are conditional model estimates, not observed sums.
    """
    q,pa=float(geom["q"]),float(geom["pa_deg"])
    ce,cn=geom["east_offset_arcsec"],geom["north_offset_arcsec"]
    centre=np.asarray(grid["anchor"])+grid["k"]@np.array([ce,cn])
    rho=ellipse_radius(grid["e"]-ce,grid["n"]-cn,q,pa)
    pad=(0,0) if kernel is None else tuple(s//2 for s in kernel.shape)
    rpad=radius_pixels(rho.shape,centre,grid["k"],q,pa,
        origin=(grid["x"][0,0],grid["y"][0,0]),pad=pad)
    edges=[inner]+[a for a in range(5,int(outer)+1,5) if a>inner]
    rows=[]; selected=[]; obs=[]
    for lo,hi in zip(edges[:-1],edges[1:]):
        domain=(rho>=lo)&(rho<hi); use=domain&valid
        count=int(use.sum())
        ideal=np.pi*q*(hi**2-lo**2)/grid["area"]
        cover=count/ideal
        row=dict(lo_arcsec=lo,hi_arcsec=hi,observed_pixels=count,
            geometric_pixels_approx=float(ideal),coverage=float(cover),
            signed_flux_nMgy=float(image[use].sum()) if count else None,
            signed_mean_nMgy_arcsec2=float(image[use].mean()/grid["area"]) if count else None)
        rows.append(row)
        if cover>=.70 and count>=20:
            selected.append(use); obs.append(row["signed_mean_nMgy_arcsec2"])
    info=dict(status="SERSIC_INSUFFICIENT_RADIAL_SUPPORT",outer_arcsec=outer,inner_arcsec=inner,
        q=q,pa_deg=pa,rows=rows,used_rings=len(obs),science_adopted=False)
    # Missing outer bins can otherwise produce a deceptively constrained tail.
    if len(obs)<5 or rows[-1]["coverage"]<.70: return info
    obs=np.asarray(obs)
    norm=float(np.max(abs(obs)))
    if norm<=1e-12 or not np.isfinite(norm):
        return dict(info,status="SERSIC_NO_MEASURABLE_CONTRAST")
    def model(p):
        return np.exp(p[0])*convolved_profile(rpad,np.exp(p[1]),p[2],kernel)
    def residual(p):
        m=model(p)
        return np.asarray([m[u].mean() for u in selected])-obs/norm
    low=[-10.,np.log(3.),.35]; high=[8.,np.log(100.),2.5]
    fits=[]
    for re,n in ((20.,.55),(35.,1.0),(50.,1.6)):
        f=least_squares(residual,[0.,np.log(re),n],bounds=(low,high),
            loss="soft_l1",f_scale=.05,x_scale="jac",max_nfev=100)
        if f.success and np.isfinite(f.x).all() and np.isfinite(f.cost): fits.append(f)
    if not fits: return dict(info,status="SERSIC_OPTIMIZER_BLOCKED")
    best=min(fits,key=lambda f:f.cost)
    i0=float(np.exp(best.x[0])*norm); re=float(np.exp(best.x[1])); n=float(best.x[2])
    full=total_flux(i0,re,n,q); bn=gammaincinv(2*n,.5)
    predicted=model(best.x)*norm*grid["area"]
    flags=[]
    if re<=3.05 or re>=99. or n<=.355 or n>=2.49: flags.append("PROFILE_BOUND")
    if not np.isfinite(full) or full<=0: return dict(info,status="SERSIC_INVALID_TOTAL")
    tail=1-float(gammainc(2*n,bn*(outer/re)**(1/n)))
    if tail>.30: flags.append("MODEL_TAIL_GT_30_PERCENT")
    if tail>.50: flags.append("MODEL_TAIL_GT_50_PERCENT")
    if len(fits)<3: flags.append("SOME_STARTS_UNCONVERGED")
    for f in fits:
        if f.cost<=best.cost*1.05+1e-8:
            alt=total_flux(float(np.exp(f.x[0])*norm),float(np.exp(f.x[1])),float(f.x[2]),q)
            if abs(alt/full-1)>.10:
                flags.append("MULTISTART_TOTAL_DISAGREEMENT");break
    # Native observed sums and separate model predictions on exactly
    # the same observed support; no replacement of missing pixels.
    aperture=[]
    for rad in RADII:
        if rad>outer: continue
        use=(rho<rad)&valid
        aperture.append(dict(a_arcsec=rad,observed_pixels=int(use.sum()),
            coverage_approx=float(use.sum()*grid["area"]/(np.pi*q*rad**2)),
            observed_signed_flux_nMgy=float(image[use].sum()),
            model_on_observed_pixels_nMgy=float(predicted[use].sum()),
            intrinsic_model_enclosed_nMgy=float(full*gammainc(2*n,bn*(rad/re)**(1/n)))))
    for row in rows:
        use=(rho>=row["lo_arcsec"])&(rho<row["hi_arcsec"])&valid
        row["model_mean_on_observed_pixels_nMgy_arcsec2"]=float(predicted[use].mean()/grid["area"]) if use.any() else None
    return dict(info,status="SERSIC_CONDITIONAL_NOT_ADOPTED",re_arcsec=re,n=n,
        mu0_mag_arcsec2=float(22.5-2.5*np.log10(i0)),i0_nMgy_arcsec2=i0,
        analytic_total_nMgy=full,analytic_total_mag=float(22.5-2.5*np.log10(full)),
        intrinsic_tail_fraction=tail,robust_cost=float(best.cost),flags=flags,
        converged_starts=len(fits),apertures=aperture,
        warning="Finite-core PSF model only; no calibrated covariance or empirical asymptote. Analytic total assumes Sersic continuation to infinity.")


def coverage_sites(mask,wcs,tx,ty,scale,max_sites=3):
    """Predeclared geometry/coverage only; never select on flux residuals."""
    ny,nx=mask.shape; step=max(1,int(round(15./scale)))
    pad=int(np.ceil(50./scale)); r=int(np.ceil(42./scale))
    candidates=[]
    for y in range(pad,ny-pad,step):
        for x in range(pad,nx-pad,step):
            if np.hypot(x-tx,y-ty)*scale<100.: continue
            sky=wcs.celestial.pixel_to_world(x,y).icrs
            g=ellipse_from_wcs(wcs,sky.ra.deg,sky.dec.deg,sky_q=1.,sky_pa_deg=0.)
            cut=mask[y-r:y+r+1,x-r:x+r+1]
            scores=[]
            for _,_,q,pa,_ in SHAPES:
                rho=radius_pixels(cut.shape,(x,y),g.local_pixel_per_arcsec_en,q,pa,
                    origin=(x-r,y-r))
                for lo,hi in ((0,5),(0,20),(0,40),(5,20),(20,35),(35,40)):
                    area=(rho>=lo)&(rho<hi)
                    cover=float(cut[area].mean()) if area.any() else 0.
                    scores.append(cover)
            if min(scores)>=.75:
                candidates.append(dict(x_pix=x,y_pix=y,minimum_elliptical_coverage=min(scores),
                    mean_elliptical_coverage=float(np.mean(scores))))
    sites=[]
    for s in sorted(candidates,key=lambda a:(-a["minimum_elliptical_coverage"],
            -a["mean_elliptical_coverage"],a["y_pix"],a["x_pix"])):
        if all(np.hypot(s["x_pix"]-v["x_pix"],s["y_pix"]-v["y_pix"])*scale>=90. for v in sites):
            sites.append(s)
        if len(sites)>=max_sites: break
    return sites


def halo_refit(image,inputs,extra_protected,extra_veto,reference):
    result=measure_saturated_halo_sectors(image,inputs["iv"],inputs["hard"],
        source_xy=inputs["source_xy"],target_xy=inputs["target_xy"],
        pixel_scale_arcsec=inputs["scale"],source_mask=inputs["source_mask"],
        protected=inputs["protected"]|extra_protected,
        unresolved_sky_veto=inputs["veto"]|extra_veto,
        neighbor_centres=inputs["neighbours"],compare_models=True,
        expanded=True,diagnostic_tile_map=False)
    c=next((c for c in result["sh1b_model_comparison"]["candidates"] if key(c)==key(reference)),None)
    if c is None: return None,None
    h,b=model_arrays(image.shape,inputs["source_xy"],inputs["scale"],c)
    return image-h-b,dict(hypothesis=list(key(c)),halo_at_knots_nMgy_pix=c["halo_at_knots_nMgy_pix"],
        radial_knots_arcsec=c["radial_knots_arcsec"],training_r_range_arcsec=c["training_r_range_arcsec"],
        sky_constant_nMgy_pix=c["sky_constant_nMgy_pix"],
        sky_x_nMgy_pix_per_arcsec=c["sky_x_nMgy_pix_per_arcsec"],
        sky_y_nMgy_pix_per_arcsec=c["sky_y_nMgy_pix_per_arcsec"])


def injection_case(inputs,wcs,mask,site,params,kernel,reference,f40,extra_veto,control=None):
    x,y=site["x_pix"],site["y_pix"]; re,n,q,pa,ratio=params
    scale=inputs["scale"]; tx,ty=inputs["target_xy"]
    raw=np.asarray(inputs["image"],float)
    sky=wcs.celestial.pixel_to_world(x,y).icrs
    geo=ellipse_from_wcs(wcs,sky.ra.deg,sky.dec.deg,sky_q=1.,sky_pa_deg=0.)
    k=geo.local_pixel_per_arcsec_en
    rho=radius_pixels(raw.shape,(x,y),k,q,pa)
    padded=radius_pixels(raw.shape,(x,y),k,q,pa,pad=tuple(v//2 for v in kernel.shape))
    image=convolved_profile(padded,re,n,kernel)
    norm=float(image[rho<40].sum())
    if norm<=0: return dict(status="INJ2_MODEL_FLUX_BLOCKED"),control
    amp=ratio*f40/norm
    model=image*amp
    i0=amp/geo.pixel_area_arcsec2
    truth=dict(re_arcsec=re,n=n,q=q,pa_deg=pa,mu0_mag_arcsec2=float(22.5-2.5*np.log10(i0)),
        analytic_total_nMgy=total_flux(i0,re,n,q),target_f40_ratio=ratio,
        full_sampled_convolved_f40_nMgy=float(model[rho<40].sum()),
        negative_model_mass_fraction=float(-np.minimum(model,0).sum()/max(model.sum(),1e-30)))
    result=dict(status="INJ2_NOT_ADOPTED",site=site,truth=truth,models={})
    if truth["negative_model_mass_fraction"]>1e-5:
        return dict(result,status="INJ2_SIGNED_PSF_MODEL_BLOCKED"),control
    yy,xx=np.ogrid[:raw.shape[0],:raw.shape[1]]
    guard=np.hypot(xx-x,yy-y)*scale<60.
    if control is None:
        control=halo_refit(raw,inputs,guard,extra_veto,reference)
    before,h0=control
    after,h1=halo_refit(raw+model,inputs,guard,extra_veto,reference)
    if before is None or after is None:
        return dict(result,status="INJ2_MATCHED_HALO_REFIT_BLOCKED"),control
    result["halo_control"]=h0; result["halo_injected"]=h1
    common=mask&np.isfinite(before)&np.isfinite(after)
    train=common&~guard&(np.hypot(xx-tx,yy-ty)*scale>=GUARD)
    b0=tile_rows(before,train,tx,ty,scale,TILE)
    b1=tile_rows(after,train,tx,ty,scale,TILE)
    result["training_blocks"]=len(b0)
    result["training_mask_sha256"]=native_boolean_mask_sha256(train)
    if [(b["tile_x_pix"],b["tile_y_pix"],b["n"]) for b in b0]!=[(b["tile_x_pix"],b["tile_y_pix"],b["n"]) for b in b1]:
        raise ValueError("INJ2 training pixel mismatch")
    points=np.asarray([[r["x_pix"],r["y_pix"]] for r in b0])
    angle=np.linspace(0,2*np.pi,32,endpoint=False)
    check=np.column_stack((x+45./scale/np.cos(np.pi/32)*np.cos(angle),
                           y+45./scale/np.cos(np.pi/32)*np.sin(angle)))
    try: inside=bool(np.all(Delaunay(points).find_simplex(check)>=0))
    except (QhullError,ValueError): inside=False
    result["injection_hull_supported"]=inside
    grid=local_grid(dict(target_xy_native_pix=[x,y],layers={"original":raw}),wcs)
    valid=common[grid["sl"]]
    d=block_design(grid,valid,35.,3.)
    true_rho=rho[grid["sl"]]
    # Re-check support after the upstream halo refit on the same truth
    # ellipse, including the central region and outer ring.
    covers=[]
    for lo,hi in ((0,5),(0,20),(0,40),(35,40)):
        region=(true_rho>=lo)&(true_rho<hi)
        covers.append(float(valid[region].mean()) if region.any() else 0.)
    result["post_refit_elliptical_coverage"]=covers
    if min(covers)<.75:
        return dict(result,status="INJ2_POST_REFIT_COVERAGE_BLOCKED"),control
    dx=(grid["x"]-tx)*scale; dy=(grid["y"]-ty)*scale
    for method in METHODS:
        f0=f1=None
        if method!="zero":
            if not inside or min(sum(r["quadrant"]==i for r in b0) for i in range(4))<3:
                result["models"][method]=dict(status="INJ2_SKY_SUPPORT_BLOCKED");continue
            f0=fit_blocks(b0,method);f1=fit_blocks(b1,method)
            if f0 is None or f1 is None:
                result["models"][method]=dict(status="INJ2_SKY_FIT_BLOCKED");continue
        ctrl=before[grid["sl"]]-(predict(f0,dx,dy) if f0 else 0.)
        recovered=after[grid["sl"]]-(predict(f1,dx,dy) if f1 else 0.)
        rec=dict(status="INJ2_CONDITIONAL_NOT_ADOPTED",sky_before=f0,sky_after=f1,apertures=[])
        for radius in (20.,35.,40.):
            use=valid&(true_rho<radius)
            f=float(model[grid["sl"]][use].sum())
            if f<=0:continue
            a=float(recovered[use].sum()); c=float(ctrl[use].sum())
            rec["apertures"].append(dict(a_arcsec=radius,truth_observed_nMgy=f,
                absolute_recovered_nMgy=a,control_nMgy=c,
                absolute_fractional_error=a/f-1,differential_fractional_error=(a-c)/f-1))
        geom=fit_geometry(d,recovered,grid,fixed_background_nMgy_pix=0.)
        rec["geometry"]=geom
        if fitted(geom):
            rec["geometry_errors"]=dict(centre_arcsec=float(np.hypot(geom["east_offset_arcsec"],geom["north_offset_arcsec"])),
                delta_q=float(geom["q"]-q),delta_pa_deg=float((geom["pa_deg"]-pa+90)%180-90))
            fit=fit_profile(recovered,valid,grid,geom,kernel)
            rec["sersic"]=fit
            if fit["status"]=="SERSIC_CONDITIONAL_NOT_ADOPTED":
                rec["sersic_errors"]=dict(delta_n=fit["n"]-n,
                    re_fractional_error=fit["re_arcsec"]/re-1,
                    delta_mu0=fit["mu0_mag_arcsec2"]-truth["mu0_mag_arcsec2"],
                    total_fractional_error=fit["analytic_total_nMgy"]/truth["analytic_total_nMgy"]-1)
        result["models"][method]=rec
    return result,control


def target_analysis(products,datasets,validation,kernels):
    """Native forward-convolved totals plus strict matched-PSF apertures."""
    refband="r" if "r" in products else next(iter(products))
    p=products[refband]["probe"]["sh1c_preview"]
    ref=validation["target"]["reference_geometry"]
    grid=local_grid(p,datasets[refband]["wcs"])
    anchor=p["target_xy_native_pix"]
    images={}; fits={}; layers={}
    valid=np.ones(grid["shape"],bool)
    for band,b in products.items():
        other=b["probe"]["sh1c_preview"]
        if (datasets[band]["grid_fingerprint"]!=datasets[refband]["grid_fingerprint"] or
            other["layers"]["residual"].shape!=grid["shape"] or
            not np.allclose(other["target_xy_native_pix"],anchor,rtol=0,atol=1e-6)):
            raise ValueError("INJ2 target grid mismatch")
        raw=other["layers"]["residual"]
        valid &= np.isfinite(raw)&np.isfinite(other["layers"]["original"])
        yy,xx=np.ogrid[:raw.shape[0],:raw.shape[1]]
        scale=other["pixel_scale_arcsec"]
        images[band]={"zero":raw}
        for c in validation["bands"][band]["after"]["candidates"]:
            if c["status"]=="SKY1_CONDITIONAL_NOT_ADOPTED":
                images[band][c["method"]]=raw-predict(c["fit"],(xx-anchor[0])*scale,(yy-anchor[1])*scale)
    # All full-image corrections outside the target hull are internal
    # scratch arrays only. Matching is restricted to the central domain
    # whose entire convolution footprint is supported by the sky fit.
    local=valid[grid["sl"]]
    for band,variants in images.items():
        fits[band]={}
        for name,im in variants.items():
            if band not in kernels:
                fits[band][name]=dict(status="SERSIC_PSF_UNAVAILABLE");continue
            fits[band][name]=[fit_profile(im[grid["sl"]],local,grid,ref,kernels[band],outer=a,inner=i)
                for a,i in ((30.,3.),(35.,3.),(40.,3.),(40.,5.))]
    matched=dict(status="PSF_MATCH_UNAVAILABLE",science_adopted=False)
    if set(images)=={"g","r"} and set(kernels)=={"g","r"}:
        names=sorted(set(images["g"])&set(images["r"]))
        yy,xx=np.ogrid[:valid.shape[0],:valid.shape[1]]
        scale=p["pixel_scale_arcsec"]
        # Every accepted SKY-CLEAN candidate guarantees this disk is in
        # its training hull. Erode by both native kernel boxes.
        radius=min(products[b]["sky"]["expanded_sky"]["required_target_support_radius_arcsec"] for b in images)
        safe=valid&(np.hypot(xx-anchor[0],yy-anchor[1])*scale<=radius)
        rho=ellipse_radius(grid["e"]-ref["east_offset_arcsec"],grid["n"]-ref["north_offset_arcsec"],ref["q"],ref["pa_deg"])
        rows=[]
        for name in names:
            gm,gv=strict_convolve(images["g"][name],safe,kernels["r"])
            rm,rv=strict_convolve(images["r"][name],safe,kernels["g"])
            common=(gv&rv)[grid["sl"]]
            a=gm[grid["sl"]]; b=rm[grid["sl"]]
            for radius0 in RADII:
                use=common&(rho<radius0); original=local&(rho<radius0)
                fg=float(a[use].sum()); fr=float(b[use].sum())
                coverage=float(use.sum()*grid["area"]/(np.pi*ref["q"]*radius0**2))
                rows.append(dict(method=name,a_arcsec=radius0,common_pixels=int(use.sum()),
                    coverage_approx=coverage,native_common_pixels=int(original.sum()),
                    g_signed_nMgy=fg,r_signed_nMgy=fr,
                    native_g_same_pixels_nMgy=float(images["g"][name][grid["sl"]][use].sum()),
                    native_r_same_pixels_nMgy=float(images["r"][name][grid["sl"]][use].sum()),
                    g_minus_r_mag=float(-2.5*np.log10(fg/fr)) if fg>0 and fr>0 and coverage>=.75 else None))
            if name=="plane":
                layers={"g_matched_plane":np.where(common,a,np.nan),
                        "r_matched_plane":np.where(common,b,np.nan)}
        effective=fftconvolve(kernels["g"],kernels["r"],mode="full")
        matched=dict(status="FINITE_CORE_CROSS_CONVOLUTION_NOT_ADOPTED",science_adopted=False,
            effective_kernel_shape=list(effective.shape),effective_kernel_sum=float(effective.sum()),
            effective_kernel_negative_mass=float(-np.minimum(effective,0).sum()),apertures=rows,
            warning="Cross convolution equalizes stored finite core kernels only. Entire kernel box must be observed; no mask filling. Colour is conditional; covariance and far wings unvalidated.")
    return dict(status="TARGET_SERSIC_NOT_ADOPTED",science_adopted=False,reference_band=refband,
        native_profile_fits=fits,matched_psf=matched,layers=layers)


def run_inj2(products,datasets,validation,kernels,evidence):
    out=dict(schema="desidera.inj2-sersic/1",status="INJ2_NOT_ADOPTED",
        science_adopted=False,psf=evidence,sites=[],bands={},
        warning="Fixed source masks, injected smooth expectation only (no additional Poisson noise); radial SH1 halo and sky refit with the same hypothesis, no directional spike re-estimation. Conditional finite-core PSF and Sersic extrapolation, not calibrated errors.")
    if validation.get("target",{}).get("status")!="TARGET_COMMON_PROFILES_NOT_ADOPTED":
        return dict(out,status="INJ2_TARGET_EVIDENCE_BLOCKED")
    out["target"]=target_analysis(products,datasets,validation,kernels)
    if set(kernels)!=set(products):
        return dict(out,status="INJ2_PSF_EVIDENCE_BLOCKED")
    clean=[]; old=[]; grids={}
    for band,b in products.items():
        s=b["sky"]["expanded_sky"]
        grids[band]=expanded_wcs(b["probe"],datasets[band],s)
        base=s["expanded_safe_mask"]
        m=validation["masks"][band]
        old.append(base)
        clean.append(base&~m["saturation"]&~m["compact"])
    common=np.logical_and.reduce(clean); common0=np.logical_and.reduce(old)
    first=next(iter(products))
    s=products[first]["sky"]["expanded_sky"]; tx,ty=s["expanded_target_xy"]; scale=s["pixel_scale_arcsec"]
    yy,xx=np.ogrid[:common.shape[0],:common.shape[1]]
    if native_boolean_mask_sha256(common&(np.hypot(xx-tx,yy-ty)*scale>=GUARD))!=validation["common_mask_sha256"]:
        raise ValueError("INJ2 clean mask identity mismatch")
    sites=coverage_sites(common,grids[first],tx,ty,scale)
    out["sites"]=sites
    target=next(r for r in validation["target"]["profiles"] if r["sector"]=="ALL" and r["outer_arcsec"]==40.)
    for band,b in products.items():
        p=b["probe"]["sh1c_preview"]; inputs=p["sky1_inputs"]
        ref=b["probe"]["sh1_sector_measurement"]["sh1b_model_comparison"]["candidates"][p["reference_candidate_index"]]
        f40=target["scenarios"].get(band+"_plane",target["scenarios"][band+"_zero"])["cumulative_flux_nMgy"]
        if f40 is None or f40<=0:
            out["bands"][band]=dict(status="INJ2_TARGET_FLUX_BLOCKED",cases=[]);continue
        cases=[]
        for site in sites:
            ctrl=None
            for params in SHAPES:
                case,ctrl=injection_case(inputs,grids[band],common,site,params,kernels[band],ref,f40,common0&~common,ctrl)
                cases.append(case)
        out["bands"][band]=dict(target_normalization_f40_nMgy=f40,cases=cases)
    if not sites: out["status"]="INJ2_NO_SUITABLE_ELLIPTICAL_SUPPORT"
    return out
