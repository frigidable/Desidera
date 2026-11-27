"""OUTER-PHOT / INJ-3, conditional expanded-grid profile and sensitivity.
No mask filling, empirical-total adoption, or calibrated error claim.
"""
import numpy as np
from scipy.spatial import Delaunay, QhullError
from .sky_clean_injection import expanded_wcs
from .sersic_recovery import (ellipse_radius, radius_pixels, halo_refit,
    fit_profile, injection_case)
from .stellar_halo_spikes import model_arrays
from .galaxy_geometry_diagnostic import local_grid
from .residual_sky_diagnostic import tile_rows, candidate, predict
from .sky_masks import native_boolean_mask_sha256
from ..sky_geometry import ellipse_from_wcs

GUARDS=(60.,75.,90.)
OUTERS=(40.,50.,60.,70.)
SHAPES=((22.,.4,.56,44.5,.85),(26.,.5,.56,44.5,1.),
        (29.,.7,.60,49.,1.15))
METHODS=("constant","plane","quadratic")
TILE=32


def hull_support(rows,centres,radii,scale):
    """Largest fully enclosed circular domain, conservative 32-gon."""
    pts=np.asarray([[r["x_pix"],r["y_pix"]] for r in rows])
    try:
        hull=Delaunay(pts)
    except (ValueError,QhullError):
        return []
    t=np.linspace(0,2*np.pi,32,endpoint=False)
    unit=np.column_stack((np.cos(t),np.sin(t)))/np.cos(np.pi/32)
    return [float(rad) for rad in radii if all(np.all(hull.find_simplex(
        np.asarray(c)+unit*rad/scale)>=0) for c in centres)]


def outer_rows(grid,valid,images,geom,source_x,scale,max_radius=70.):
    """Compare signed sums on common support; nominal areas include off-frame pixels."""
    de=grid["e"]-geom["east_offset_arcsec"]; dn=grid["n"]-geom["north_offset_arcsec"]
    t=np.deg2rad(geom["pa_deg"]); q=geom["q"]
    x=de*np.sin(t)+dn*np.cos(t); y=de*np.cos(t)-dn*np.sin(t)
    rho=np.hypot(x,y/q)
    angle=np.mod(np.arctan2(y/q,x),2*np.pi)
    rows=[]
    for hi in range(5,int(max_radius)+1,5):
        lo=hi-5; ring=(rho>=lo)&(rho<hi); ap=rho<hi
        ideal=np.pi*q*(hi**2-lo**2)/grid["area"]
        ideal_ap=np.pi*q*hi**2/grid["area"]
        use=valid&ring; ua=valid&ap
        sectors=[]
        for j in range(8):
            u=use&(angle>=j*np.pi/4)&(angle<(j+1)*np.pi/4)
            sectors.append(dict(index=j,pixels=int(u.sum()),coverage_approx=float(u.sum()/(ideal/8)),
                means={k:float(im[u].mean()/grid["area"]) if u.any() else None for k,im in images.items()}))
        row=dict(inner_arcsec=lo,outer_arcsec=hi,ring_pixels=int(use.sum()),
            ring_coverage_approx=float(use.sum()/ideal),aperture_coverage_approx=float(ua.sum()/ideal_ap),
            in_frame_ring_coverage_approx=float(ring.sum()/ideal),
            sectors=sectors,scenarios={})
        row["support_screen"]=bool(row["ring_coverage_approx"]>=.70 and
            row["aperture_coverage_approx"]>=.70 and min(s["coverage_approx"] for s in sectors)>=.5)
        for k,im in images.items():
            v=dict(ring_flux_nMgy=float(im[use].sum()) if use.any() else None,
                mean_nMgy_arcsec2=float(im[use].mean()/grid["area"]) if use.any() else None,
                cumulative_flux_nMgy=float(im[ua].sum()) if ua.any() else None,
                strip_trials=[])
            for width in (2.,4.):
                kept=use&(abs(grid["x"]-source_x)*scale>=width)
                ak=ua&(abs(grid["x"]-source_x)*scale>=width)
                v["strip_trials"].append(dict(half_width_arcsec=width,ring_pixels=int(kept.sum()),
                    mean_nMgy_arcsec2=float(im[kept].mean()/grid["area"]) if kept.any() else None,
                    cumulative_pixels=int(ak.sum()),cumulative_flux_nMgy=float(im[ak].sum()) if ak.any() else None))
            row["scenarios"][k]=v
        rows.append(row)
    return rows


def plateau(rows,name):
    """Require contiguous radial support; last increments are descriptive, not a noise gate."""
    good=[]
    for r in rows:
        if not r["support_screen"]: break
        good.append(r)
    if len(good)<3:
        return dict(status="OUTER_INSUFFICIENT_CONTIGUOUS_COVERAGE",
                    supported_outer_arcsec=good[-1]["outer_arcsec"] if good else 0.)
    a=good[-1]["scenarios"][name]["cumulative_flux_nMgy"]
    increments=[r["scenarios"][name]["ring_flux_nMgy"] for r in good[-2:]]
    if a is None or a<=0 or any(v is None for v in increments):
        return dict(status="OUTER_NONPOSITIVE_OR_MISSING_FLUX",supported_outer_arcsec=good[-1]["outer_arcsec"])
    v=[float(x/a) for x in increments]
    return dict(status="OUTER_SMALL_INCREMENTS_UNCALIBRATED" if max(map(abs,v))<=.02 else "OUTER_GROWTH_CONTINUES",
        supported_outer_arcsec=good[-1]["outer_arcsec"],last_two_signed_increments_over_flux=v,
        science_adopted=False)


def preflight_sites(mask,res,wcs,tx,ty,scale,outer,guard,max_sites):
    """Coverage and training-hull screening before any injection, no flux ranking."""
    ny,nx=mask.shape; pad=int(np.ceil((outer+10.)/scale))
    r=int(np.ceil((outer+2.)/scale)); step=max(1,round(15./scale))
    candidates=[]
    for y in range(pad,ny-pad,step):
        for x in range(pad,nx-pad,step):
            if np.hypot(x-tx,y-ty)*scale<100.: continue
            world=wcs.celestial.pixel_to_world(x,y).icrs
            g=ellipse_from_wcs(wcs,world.ra.deg,world.dec.deg,sky_q=1.,sky_pa_deg=0.)
            cut=mask[y-r:y+r+1,x-r:x+r+1]
            coverage=[]
            for _,_,q,pa,_ in SHAPES:
                rho=radius_pixels(cut.shape,(x,y),g.local_pixel_per_arcsec_en,q,pa,origin=(x-r,y-r))
                for lo,hi in ((0,5),(0,20),(0,40),(0,outer),(outer-5.,outer)):
                    u=(rho>=lo)&(rho<hi)
                    coverage.append(float(cut[u].mean()) if u.any() else 0.)
            if min(coverage)>=.75:
                candidates.append(dict(x_pix=x,y_pix=y,minimum_elliptical_coverage=min(coverage)))
    # Precompute conservative tiles fully outside the artificial guard.
    yy,xx=np.ogrid[:ny,:nx]
    base=mask&(np.hypot(xx-tx,yy-ty)*scale>=guard)
    rows=tile_rows(res,base,tx,ty,scale,TILE)
    selected=[]; blocked=0; passed=0
    for s in sorted(candidates,key=lambda v:(-v["minimum_elliptical_coverage"],v["y_pix"],v["x_pix"])):
        x,y=s["x_pix"],s["y_pix"]
        # Entire tile disjoint from injection guard, stricter than actual
        # masked tile support; no use of block means in site selection.
        rr=[b for b in rows if np.hypot(
            b["tile_x_pix"]+(TILE-1)/2-x,b["tile_y_pix"]+(TILE-1)/2-y)*scale
            >=guard+np.sqrt(2)*TILE*scale/2]
        if min(sum(b["quadrant"]==j for b in rr) for j in range(4))<3 or not hull_support(rr,[(x,y)],[outer+5.],scale):
            blocked+=1;continue
        passed+=1
        if len(selected)<max_sites and all(np.hypot(x-v["x_pix"],y-v["y_pix"])*scale>=90. for v in selected):
            selected.append(dict(s,preflight_training_blocks=len(rr)))
    return dict(outer_arcsec=outer,guard_arcsec=guard,coverage_candidates=len(candidates),
        hull_blocked=blocked,hull_passed=passed,sites=selected)


def _range(values):
    return [float(min(values)),float(max(values))] if values else None


def sensitivity(fits,inj):
    """Factor-specific ranges, no quadrature or statistical error bars."""
    good=lambda f:f.get("status")=="SERSIC_CONDITIONAL_NOT_ADOPTED"
    metrics=("analytic_total_mag","re_arcsec","n","mu0_mag_arcsec2")
    def spans(seq):
        seq=[f for f in seq if good(f)]
        return dict(count=len(seq),ranges={k:_range([f[k] for f in seq]) for k in metrics},
                    flagged_count=sum(bool(f.get("flags")) for f in seq))
    result={}
    for band,seq in fits.items():
        method="constant" if band=="g" else "plane"
        base=[f for f in seq if f["halo_scenario"]=="reference" and f["guard_arcsec"]==75. and f["method"]==method and f["outer_arcsec"]==40. and f["inner_arcsec"]==3. and good(f)]
        ref=base[0] if base else None
        b=dict(reference=ref,reference_rule="a priori reporting convention, not adopted",
            sky=spans([f for f in seq if f["halo_scenario"]=="reference" and f["guard_arcsec"]==75. and f["outer_arcsec"]==40. and f["inner_arcsec"]==3.]),
            radius=spans([f for f in seq if f["halo_scenario"]=="reference" and f["guard_arcsec"]==75. and f["method"]==method and f["inner_arcsec"]==3.]),
            inner=spans([f for f in seq if f["halo_scenario"]=="reference" and f["guard_arcsec"]==75. and f["method"]==method and f["outer_arcsec"]==40.]),
            guard_and_halo=spans([f for f in seq if f["method"]==method and f["outer_arcsec"]==40. and f["inner_arcsec"]==3.]),
            all_models=spans(seq),injection_errors={})
        for group in inj:
            cases=group.get("bands",{}).get(band,[])
            for m in ("zero",*METHODS):
                records=[c.get("models",{}).get(m,{}) for c in cases]
                errors=[v["sersic_errors"] for v in records if "sersic_errors" in v]
                b["injection_errors"][f"{group['outer_arcsec']:g}/{m}"]=dict(count=len(errors),
                    total_fractional_error_range=_range([e["total_fractional_error"] for e in errors]),
                    re_fractional_error_range=_range([e["re_fractional_error"] for e in errors]),
                    delta_n_range=_range([e["delta_n"] for e in errors]))
        result[band]=b
    return result


def run_outer(products,datasets,validation,kernels,psf):
    out=dict(schema="desidera.outer-phot-inj3/1",status="OUTER_INJ3_NOT_ADOPTED",
        science_adopted=False,psf=psf,profiles={},fits={},injections=[],layers={},
        warning="Observed signed sums, no missing flux reconstruction. Spatial models conditional on frozen masks and radial halo hypothesis; spike-strip exclusions change support. Model ranges and injection errors are sensitivity diagnostics, not calibrated confidence intervals.")
    if validation.get("target",{}).get("status")!="TARGET_COMMON_PROFILES_NOT_ADOPTED":
        return dict(out,status="OUTER_TARGET_EVIDENCE_BLOCKED")
    ref=validation["target"]["reference_geometry"]
    first="r" if "r" in products else next(iter(products))
    s0=products[first]["sky"]["expanded_sky"]; tx,ty=s0["expanded_target_xy"]; scale=s0["pixel_scale_arcsec"]
    shape=tuple(s0["expanded_shape"])
    yy,xx=np.ogrid[:shape[0],:shape[1]]
    rt=np.hypot(xx-tx,yy-ty)*scale
    grids={}; clean=[]; old=[]; valid=[]; refs={}; ins={}; residuals={}
    for band,b in products.items():
        s=b["sky"]["expanded_sky"]; p=b["probe"]["sh1c_preview"]
        if (s["expanded_grid_fingerprint"]!=s0["expanded_grid_fingerprint"] or tuple(s["expanded_shape"])!=shape
            or not np.allclose(s["expanded_target_xy"],[tx,ty],rtol=0,atol=1e-6)
            or not np.allclose(s["expanded_source_xy"],s0["expanded_source_xy"],rtol=0,atol=1e-6)
            or not np.isclose(s["pixel_scale_arcsec"],scale,rtol=0,atol=1e-9)):
            raise ValueError("OUTER expanded identity mismatch")
        grids[band]=expanded_wcs(b["probe"],datasets[band],s)
        i=p["sky1_inputs"]; ins[band]=i
        c=b["probe"]["sh1_sector_measurement"]["sh1b_model_comparison"]["candidates"][p["reference_candidate_index"]]
        refs[band]=c
        h,sky=model_arrays(shape,i["source_xy"],scale,c)
        image=np.asarray(i["image"],float)
        residuals[band]={"reference":image-h-sky}
        # Match the SH1-C ensemble's common radial support, not just the reference.
        lo,hi=p["common_radial_support_arcsec"]; sx,sy=i["source_xy"]
        rs=np.hypot(xx-sx,yy-sy)*scale
        valid.append(~i["hard"]&~i["source_mask"]&np.isfinite(image)&np.isfinite(i["iv"])&(i["iv"]>0)&(rs>=lo)&(rs<=hi)&np.isfinite(h))
        m=validation["masks"][band]; safe=s["expanded_safe_mask"]
        old.append(safe); clean.append(safe&~m["saturation"]&~m["compact"])
    common=np.logical_and.reduce(clean)
    common0=np.logical_and.reduce(old)
    if native_boolean_mask_sha256(common&(rt>=60.))!=validation["common_mask_sha256"]:
        raise ValueError("OUTER clean mask identity mismatch")
    extra=common0&~common
    # Refit only the original hypothesis; failures remain visible.
    out["halo_refits"]={}
    for band,i in ins.items():
        out["halo_refits"][band]={}
        for guard in GUARDS[1:]:
            im,qa=halo_refit(i["image"],i,rt<guard,extra,refs[band])
            name=f"protect_{int(guard)}"
            out["halo_refits"][band][name]=qa if qa is not None else {"status":"OUTER_HALO_REFIT_BLOCKED"}
            if im is not None: residuals[band][name]=im
    names=sorted(set.intersection(*(set(r) for r in residuals.values())))
    # Hold common radial/photometric support fixed across every retained
    # halo scenario and band. Sky-only extra veto never enters this mask.
    photo=np.logical_and.reduce(valid)
    for band in residuals:
        for name in names: photo &= np.isfinite(residuals[band][name])
    grid=local_grid(dict(target_xy_native_pix=[tx,ty],layers={"original":ins[first]["image"]}),
        grids[first],margin_arcsec=80.)
    mask=photo[grid["sl"]]
    out["common_photo_mask_sha256"]=native_boolean_mask_sha256(mask)
    out["reference_geometry"]=ref
    out["fits"]={b:[] for b in products}
    source_x=s0["expanded_source_xy"][0]
    # Halo and guard vary jointly for refits; reference permits separate
    # sky-guard sensitivity. A single fixed common support across models
    # is further restricted to a fully enclosed disk for sky correction.
    for name in names:
        for guard in GUARDS:
            if name!="reference" and name!=f"protect_{int(guard)}": continue
            train=common&(rt>=guard)
            for band in products: train &= np.isfinite(residuals[band][name])
            rows_by_band={b:tile_rows(residuals[b][name],train,tx,ty,scale,TILE) for b in products}
            if len({tuple((v["tile_x_pix"],v["tile_y_pix"],v["n"]) for v in rr) for rr in rows_by_band.values()})!=1:
                raise ValueError("OUTER sky block mismatch")
            shift=float(np.hypot(ref["east_offset_arcsec"],ref["north_offset_arcsec"]))
            support=hull_support(rows_by_band[first],[(tx,ty)],[r+shift+1. for r in OUTERS],scale)
            max_a=max([r for r in OUTERS if any(abs(v-(r+shift+1.))<1e-6 for v in support)],default=0.)
            trial=f"{name}_{int(guard)}"
            item=dict(halo_scenario=name,guard_arcsec=guard,sky_supported_outer_arcsec=max_a,
                sky_training_mask_sha256=native_boolean_mask_sha256(train),models={},profiles=[],plateau={})
            out["profiles"][trial]=item
            if max_a<40.: continue
            images={}; dx=(grid["x"]-tx)*scale; dy=(grid["y"]-ty)*scale
            for b,rr in rows_by_band.items():
                raw=residuals[b][name][grid["sl"]]
                images[b+"_zero"]=raw
                item["models"][b]=[]
                for method in METHODS:
                    c=candidate(rr,method,[],tx,ty,scale,max_a+shift+1.)
                    item["models"][b].append(c)
                    if c["status"]=="SKY1_CONDITIONAL_NOT_ADOPTED":
                        images[b+"_"+method]=raw-predict(c["fit"],dx,dy)
            # Even if a sky model fails, preserve the zero correction.
            rows=outer_rows(grid,mask,images,ref,source_x,scale,max_a)
            item["profiles"]=rows
            item["plateau"]={k:plateau(rows,k) for k in images}
            for k,image in images.items():
                band,method=k.split("_",1)
                for outer in OUTERS:
                    if outer>max_a: continue
                    supported=item["plateau"][k]["supported_outer_arcsec"]
                    for inner in ((3.,5.) if outer==40. else (3.,)):
                        tags=dict(halo_scenario=name,guard_arcsec=guard,method=method,outer_arcsec=outer,inner_arcsec=inner)
                        if band not in kernels or outer>supported:
                            fit=dict(status="OUTER_FIT_PSF_OR_COVERAGE_BLOCKED")
                        else:
                            fit=fit_profile(image,mask,grid,ref,kernels[band],outer=outer,inner=inner)
                        out["fits"][band].append(dict(fit,**tags))
                if trial=="reference_75" and method in ("constant","plane"):
                    out["layers"][k]=np.where(mask,image,np.nan)
    # Raw baseline outer rings up to 70 even if sky cannot be interpolated.
    out["raw_reference_profiles"]=outer_rows(grid,mask,{b:residuals[b]["reference"][grid["sl"]] for b in products},
        ref,source_x,scale)
    if set(kernels)==set(products):
        target=next(r for r in validation["target"]["profiles"] if r["sector"]=="ALL" and r["outer_arcsec"]==40.)
        for outer,guard,max_sites in ((40.,60.,2),(60.,75.,1)):
            group=preflight_sites(common,residuals[first]["reference"],grids[first],tx,ty,scale,outer,guard,max_sites)
            group["bands"]={}
            for band,i in ins.items():
                f40=target["scenarios"].get(band+"_plane",target["scenarios"][band+"_zero"])["cumulative_flux_nMgy"]
                cases=[]
                if f40 is not None and f40>0:
                    for site in group["sites"]:
                        ctrl=None
                        for params in SHAPES:
                            case,ctrl=injection_case(i,grids[band],common,site,params,kernels[band],refs[band],f40,extra,ctrl,
                                outer=outer,guard_arcsec=guard)
                            cases.append(case)
                group["bands"][band]=cases
            out["injections"].append(group)
    else: out["injection_status"]="INJ3_PSF_EVIDENCE_BLOCKED"
    out["sensitivity"]=sensitivity(out["fits"],out["injections"])
    return out
