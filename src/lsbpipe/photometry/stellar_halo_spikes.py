"""SH1-E: native-axis strip ablation and away-only mirrored spike stress test.
No target-facing fit, science adoption or scene edits.
"""
import numpy as np
from scipy.ndimage import binary_dilation
from .stellar_halo import measure_saturated_halo_sectors

EDGES=(12.,16.,20.,25.,30.,40.,50.,60.,70.,85.,100.,115.,130.,145.,160.,180.)


def coords(shape,star,target,scale,angle=0.):
    yy,xx=np.ogrid[:shape[0],:shape[1]]
    x=(xx-star[0])*scale; y=(yy-star[1])*scale
    a=np.deg2rad(angle)
    t=y*np.cos(a)+x*np.sin(a); q=x*np.cos(a)-y*np.sin(a)
    if (target[1]-star[1])*np.cos(a)+(target[0]-star[0])*np.sin(a)<0:
        t=-t; q=-q
    return t,q


def safe_pixels(a,iv,hard,source,protected,veto,neighbours,scale):
    safe=np.isfinite(a)&np.isfinite(iv)&(iv>0)&~binary_dilation(hard,iterations=3)&~source&~protected&~veto
    for x,y in neighbours:
        if not np.isfinite([x,y]).all(): continue
        r=4./scale
        x0=max(0,int(np.floor(x-r))); x1=min(a.shape[1],int(np.ceil(x+r))+1)
        y0=max(0,int(np.floor(y-r))); y1=min(a.shape[0],int(np.ceil(y+r))+1)
        if x1<=x0 or y1<=y0: continue
        yy,xx=np.ogrid[y0:y1,x0:x1]
        safe[y0:y1,x0:x1] &= (xx-x)**2+(yy-y)**2>r*r
    return safe


def model_arrays(shape,star,scale,c):
    yy,xx=np.ogrid[:shape[0],:shape[1]]
    x=(xx-star[0])*scale; y=(yy-star[1])*scale
    r=np.hypot(x,y)
    h=np.interp(r,c["radial_knots_arcsec"],c["halo_at_knots_nMgy_pix"],left=np.nan,right=0.)
    b=c["sky_constant_nMgy_pix"]+c["sky_x_nMgy_pix_per_arcsec"]*x+c["sky_y_nMgy_pix_per_arcsec"]*y
    lo,hi=c["training_r_range_arcsec"]
    return np.where((r>=lo)&(r<=hi),h,np.nan),b


def transverse_rows(res,safe,t,q,width):
    """Signed residual excess relative to local flanks, with split consistency."""
    rows=[]
    for lo,hi in zip(EDGES[:-1],EDGES[1:]):
        base=safe&np.isfinite(res)&(-t>=lo)&(-t<hi)
        left=base&(q>=-width-4)&(q<-width-1)
        right=base&(q>width+1)&(q<=width+4)
        def median(mask,n):
            v=res[mask]; return float(np.median(v)) if v.size>=n else None
        lm,rm=median(left,24),median(right,24)
        for j in range(int(width)):
            for sign in (-1,1):
                cell=base&(sign*q>=j)&(sign*q<j+1)
                if sign==-1: cell &= q<0
                med=median(cell,24)
                excess=med-(lm+rm)/2 if None not in (med,lm,rm) else None
                splits=[]
                for half in (-t<(lo+hi)/2,-t>=(lo+hi)/2):
                    m,l,r=median(cell&half,12),median(left&half,12),median(right&half,12)
                    splits.append(m-(l+r)/2 if None not in (m,l,r) else None)
                rows.append(dict(longitudinal_inner_arcsec=lo,longitudinal_outer_arcsec=hi,
                    transverse_inner_arcsec=j,transverse_outer_arcsec=j+1,transverse_sign=sign,
                    n=int(cell.sum()),left_flank_n=int(left.sum()),right_flank_n=int(right.sum()),
                    median_residual_nMgy_pix=med,left_flank_median_nMgy_pix=lm,right_flank_median_nMgy_pix=rm,
                    signed_excess_nMgy_pix=excess,split_excess_nMgy_pix=splits,
                    template_supported=bool(excess is not None and all(v is not None for v in splits))))
    return rows


def mirrored_template(t,q,width,rows):
    """Rotate away-side cells by 180 deg. Zero outside the assumed strip."""
    strip=np.abs(q)<width
    m=np.where(strip,np.nan,0.)
    for row in rows:
        if not row["template_supported"]: continue
        lo,hi=row["longitudinal_inner_arcsec"],row["longitudinal_outer_arcsec"]
        a,b=row["transverse_inner_arcsec"],row["transverse_outer_arcsec"]
        sign=row["transverse_sign"]
        for side in (-1,1):
            tq=side*t; cq=-side*q
            v=sign*cq
            cell=strip&(tq>=lo)&(tq<hi)&(v>=a)&(v<b)
            if sign==-1: cell &= cq<0
            m[cell]=max(0.,row["signed_excess_nMgy_pix"])
    return m


def scenario_sums(original,base,refit,spike,strip,target,scale):
    yy,xx=np.ogrid[:original.shape[0],:original.shape[1]]
    rt=np.hypot(xx-target[0],yy-target[1])*scale
    common=np.isfinite(original)&np.isfinite(base)&np.isfinite(refit)
    rows=[]
    regions=[("aperture",0.,a) for a in (10.,15.,20.,25.,30.,35.,40.,45.)]
    regions += [("ring",float(a),float(a+5)) for a in range(0,45,5)]
    for kind,lo,hi in regions:
        ann=(rt>=lo)&(rt<=hi if kind=="aperture" or hi==45. else rt<hi)
        for sector,direction in (("ALL",True),("NATIVE_Y_PLUS",yy>=target[1]),("NATIVE_Y_MINUS",yy<target[1])):
            u=ann&direction&common; off=u&~strip; modeled=u&np.isfinite(spike)
            def f(a,m): return float(a[m].sum()) if m.any() else None
            rows.append(dict(kind=kind,inner_arcsec=lo,outer_arcsec=hi,sector=sector,
                common_pixels=int(u.sum()),off_strip_pixels=int(off.sum()),template_supported_pixels=int(modeled.sum()),
                raw_common_flux_nMgy=f(original,u),baseline_common_flux_nMgy=f(base,u),refit_common_flux_nMgy=f(refit,u),
                baseline_off_strip_flux_nMgy=f(base,off),refit_off_strip_flux_nMgy=f(refit,off),
                baseline_template_matched_flux_nMgy=f(base,modeled),refit_template_matched_flux_nMgy=f(refit,modeled),
                spike_template_flux_nMgy=f(spike,modeled),spike_corrected_flux_nMgy=f(refit-spike,modeled)))
    return rows


def key(c):
    return (c["tile_size_pix"],c["assumed_zero_halo_radius_arcsec"],c["sky_model"])


def ensemble_refit_deltas(p,baseline,refits,width):
    """Common pixels for ALL matched hypotheses, including both model supports."""
    pairs=[(c,refits[key(c)]) for c in baseline if key(c) in refits]
    if not pairs: return []
    image=p["layers"]["original"]; base=p["layers"]["residual"]
    sx,sy=p["source_xy_native_pix"]; tx,ty=p["target_xy_native_pix"]
    scale=p["pixel_scale_arcsec"]
    yy,xx=np.ogrid[:image.shape[0],:image.shape[1]]
    rt=np.hypot(xx-tx,yy-ty)*scale; rs=np.hypot(xx-sx,yy-sy)*scale
    lo=max(c["training_r_range_arcsec"][0] for pair in pairs for c in pair)
    hi=min(c["training_r_range_arcsec"][1] for pair in pairs for c in pair)
    valid=np.isfinite(image)&np.isfinite(base)&(rs>=lo)&(rs<=hi)&(rt<=45.)
    y,x=np.nonzero(valid); r=rs[y,x]; radius=rt[y,x]
    dx=(x-sx)*scale; dy=(y-sy)*scale
    off=np.abs(dx)>=width; values=image[y,x]
    rows=[]; selections=[]
    for ap in (20.,30.,35.,40.,45.):
        for region,mask in (("ALL",radius<=ap),("OFF_STRIP",(radius<=ap)&off)):
            selections.append(mask)
            rows.append(dict(aperture_arcsec=ap,region=region,used_pixels=int(mask.sum()),
                common_radial_support_arcsec=[lo,hi],matched_hypotheses=[]))
    def residual(c):
        h=np.interp(r,c["radial_knots_arcsec"],c["halo_at_knots_nMgy_pix"],left=np.nan,right=0.)
        b=c["sky_constant_nMgy_pix"]+dx*c["sky_x_nMgy_pix_per_arcsec"]+dy*c["sky_y_nMgy_pix_per_arcsec"]
        return values-h-b
    for c,d in pairs:
        b=residual(c); f=residual(d)
        if not (np.isfinite(b).all() and np.isfinite(f).all()):
            raise ValueError("SH1-E ensemble support mismatch")
        for row,m in zip(rows,selections):
            row["matched_hypotheses"].append(dict(hypothesis=list(key(c)),
                baseline_flux_nMgy=float(b[m].sum()) if m.any() else None,
                refit_flux_nMgy=float(f[m].sum()) if m.any() else None,
                refit_minus_baseline_nMgy=float((f[m]-b[m]).sum()) if m.any() else None))
    for row in rows:
        vals=[v["refit_minus_baseline_nMgy"] for v in row["matched_hypotheses"]
              if v["refit_minus_baseline_nMgy"] is not None]
        row["refit_delta_range_nMgy"]=[min(vals),max(vals)] if vals else None
    return rows


def directional_halo_experiment(image,iv,hard,*,source_xy,target_xy,scale,
        source_mask,protected,veto,neighbours,baseline_comparison,original_preview):
    a=np.asarray(image,dtype=float); iv=np.asarray(iv,dtype=float)
    hard,source_mask,protected,veto=[np.asarray(m,dtype=bool) for m in (hard,source_mask,protected,veto)]
    if a.ndim!=2 or any(m.shape!=a.shape for m in (iv,hard,source_mask,protected,veto)):
        raise ValueError("SH1-E native grid mismatch")
    if not np.isfinite([*source_xy,*target_xy,scale]).all() or scale<=0:
        raise ValueError("SH1-E invalid geometry")
    if original_preview.get("status")!="SH1C_CONDITIONAL_PREVIEW_NOT_ADOPTED":
        return {"status":"SH1E_BASELINE_UNAVAILABLE","science_adopted":False}
    safe=safe_pixels(a,iv,hard,source_mask,protected,veto,neighbours,scale)
    t,q=coords(a.shape,source_xy,target_xy,scale)
    candidates=baseline_comparison["candidates"]
    ref=candidates[original_preview["reference_candidate_index"]]
    bh,bb=model_arrays(a.shape,source_xy,scale,ref)
    scans=[]
    for angle in (-8.,-4.,0.,4.,8.):
        st,sq=coords(a.shape,source_xy,target_xy,scale,angle)
        scans.append(dict(angle_from_native_y_deg=angle,rows=transverse_rows(a-bh-bb,safe,st,sq,2.)))
    p=original_preview; original=p["layers"]["original"]; base=p["layers"]["residual"]
    ot,oq=coords(original.shape,p["source_xy_native_pix"],p["target_xy_native_pix"],scale)
    trials=[]; layers={}; adaptive=None
    for width in (2.,4.):
        strip=np.abs(q)<width
        fit=measure_saturated_halo_sectors(a,iv,hard,source_xy=source_xy,target_xy=target_xy,
            pixel_scale_arcsec=scale,source_mask=source_mask|strip,protected=protected,
            unresolved_sky_veto=veto,neighbor_centres=neighbours,compare_models=True,
            expanded=True,diagnostic_tile_map=False)
        comp=fit["sh1b_model_comparison"]
        by_key={key(c):c for c in comp["candidates"]}
        matches=[dict(hypothesis=list(key(c)),baseline_target_halo_nMgy_pix=c["halo_at_target_distance_nMgy_pix"],
            refit_target_halo_nMgy_pix=by_key[key(c)]["halo_at_target_distance_nMgy_pix"],
            baseline_side_residual_nMgy_pix=c["side_median_abs_residual_nMgy_pix"],
            refit_side_residual_nMgy_pix=by_key[key(c)]["side_median_abs_residual_nMgy_pix"])
            for c in candidates if key(c) in by_key]
        trial=dict(half_width_arcsec=width,excluded_safe_pixels_all_window=int((safe&strip).sum()),
            candidate_count=len(comp["candidates"]),matched_hypotheses=matches,
            rejected=comp["rejected"],refit_candidates=comp["candidates"],
            ensemble_refit_aperture_deltas=ensemble_refit_deltas(p,candidates,by_key,width),
            control_residuals_same_holdout_pixels=False)
        c=by_key.get(key(ref))
        if c is None:
            trial["status"]="SH1E_MATCHED_REFERENCE_REFIT_UNAVAILABLE"; trials.append(trial); continue
        h,b=model_arrays(a.shape,source_xy,scale,c)
        rows=transverse_rows(a-h-b,safe,t,q,width)
        oh,ob=model_arrays(original.shape,p["source_xy_native_pix"],scale,c)
        common=np.isfinite(base)&np.isfinite(oh)&np.isfinite(ob)
        res=np.where(common,original-oh-ob,np.nan)
        if width==4.:
            from .stellar_halo_spike_extent import adaptive_spike_experiment
            adaptive=adaptive_spike_experiment(a-h-b,safe,source_xy=source_xy,
                target_xy=target_xy,scale=scale,preview=p,refit_residual=res)
            adaptive["reference_hypothesis"]=list(key(c))
        spike=np.where(common,mirrored_template(ot,oq,width,rows),np.nan)
        masked=np.where(np.abs(oq)<width,np.nan,res)
        prefix=f"spike_{int(width)}"
        yy,xx=np.ogrid[:original.shape[0],:original.shape[1]]
        inside=np.hypot(xx-p["target_xy_native_pix"][0],yy-p["target_xy_native_pix"][1])*scale<=45.
        supported_target=int((inside&common&(np.abs(oq)<width)&np.isfinite(spike)).sum())
        trial["template_supported_in_strip_pixels_within_45arcsec"]=supported_target
        output_layers=[("refit_residual",res),("masked_residual",masked)]
        if supported_target:
            output_layers += [("model",spike),("halo_plus_spike",np.where(common,oh+spike,np.nan)),("corrected",res-spike)]
        for name,val in output_layers:
            val.setflags(write=False); layers[prefix+"_"+name]=val
        trial.update(status=("SH1E_CONDITIONAL_REFIT_NOT_ADOPTED" if supported_target else
                            "SH1E_REFIT_ONLY_NO_TARGET_STRIP_SUPPORT"),reference_hypothesis=list(key(ref)),
            transverse_rows=rows,template_supported_cell_count=sum(r["template_supported"] for r in rows),
            aperture_and_ring_sums=scenario_sums(original,base,res,spike,np.abs(oq)<width,p["target_xy_native_pix"],scale))
        trials.append(trial)
    return dict(schema="desidera.sh1e-directional-experiment/1",status="SH1E_NOT_ADOPTED",
        science_adopted=False,original_pixels_modified=False,scene_mask_modified=False,
        axis_definition="NATIVE_Y_FIXED; ANGLE_SCAN_DESCRIPTIVE_ONLY",axis_scan=scans,
        half_width_trials=trials,layers=layers,adaptive_extent=adaptive,
        warning="Away-only residual contrast; mirrored nonnegative template is conditional, not measured target-facing spike. Split checks descriptive, not independent covariance validation. Off-strip zero is assumed. No unsupported in-strip extrapolation or total photometry. Sector signs refer to native y.")
