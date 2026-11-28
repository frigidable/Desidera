"""SH1-E3: fixed-axis, mask-aware decay grid, diagnostic only.

Fit away-side signed residuals with local transverse linear backgrounds
estimated from flanks. Equal longitudinal-block weights, no independent
pixel noise assumption. Grid ranges are sensitivity ranges, not confidence
intervals. Reflected model is a hypothesis, never adopted.
"""
import numpy as np

LENGTHS=(3.,5.,8.,12.,20.,35.,60.,100.)
SIGMAS=(.262,.524,.786,1.048,1.572)
GRID=[(l,s) for l in LENGTHS for s in SIGMAS]
EDGES=np.arange(12.,103.,5.)
OFFSETS=(-30.,-20.,-10.,10.,20.,30.)
ANCHOR=17.


def basis(l,q,length,sigma):
    return np.exp(-(l-ANCHOR)/length-.5*(q/sigma)**2)


def prepare_blocks(residual,safe,t,q,scale,offset=0.):
    """Flank regression and common signed core vectors for every hypothesis."""
    selected=safe&np.isfinite(residual)&(-t>=EDGES[0])&(-t<EDGES[-1])&(np.abs(q-offset)<=5.)
    l=-t[selected]; x=q[selected]-offset; v=residual[selected]
    blocks=[]; coverage=[]
    for i,(lo,hi) in enumerate(zip(EDGES[:-1],EDGES[1:])):
        z=(l>=lo)&(l<hi)
        ll,xx,yy=l[z],x[z],v[z]
        c=np.abs(xx)<=2.5; left=xx<=-3.; right=xx>=3.; f=left|right
        n=[int(m.sum()) for m in (c,left,right)]
        need=[max(10,int(.25*(hi-lo)*w/scale**2)) for w in (5.,2.,2.)]
        row=dict(index=i,inner_arcsec=float(lo),outer_arcsec=float(hi),counts=n,
                 minimum_counts=need,status="UNMEASURED")
        coverage.append(row)
        if any(a<b for a,b in zip(n,need)): continue
        # Fit the nuisance background only on flanks; remove its projection
        # from both measured core and each candidate model, including tails.
        X=np.column_stack([np.ones(xx.size),xx])
        inv=np.linalg.pinv(X[f])
        def contrast(y):
            return y[c]-X[c]@(inv@y[f])
        templates=np.column_stack([basis(ll,xx,*g) for g in GRID])
        y=contrast(yy); B=contrast(templates)
        weight=1./np.sqrt(n[0])
        blocks.append(dict(index=i,y=y*weight,B=B*weight,l=ll[c],q=xx[c]))
        row.update(status="MEASURED",signed_mean_contrast_nMgy_pix=float(np.mean(y)))
    return blocks,coverage


def fit_grid(blocks):
    if not blocks: return None
    y=np.concatenate([b["y"] for b in blocks])
    B=np.concatenate([b["B"] for b in blocks],axis=0)
    denom=np.sum(B*B,axis=0)
    numerator=y@B
    amp=np.maximum(0.,np.divide(numerator,denom,out=np.zeros_like(denom),where=denom>0))
    null=float(y@y)
    loss=np.maximum(0.,null-2*amp*numerator+amp**2*denom)
    best=int(np.argmin(loss))
    return dict(best=best,amplitudes=amp,loss=loss,null=null)


def evaluate(blocks):
    fit=fit_grid(blocks)
    if fit is None: return None
    best=fit["best"]; null=fit["null"]; loss=fit["loss"]
    improvement=max(0.,null-float(loss[best]))
    folds=[]
    for parity in (0,1):
        train=[b for b in blocks if b["index"]%2!=parity]
        test=[b for b in blocks if b["index"]%2==parity]
        f=fit_grid(train)
        if f is None or len(train)<2 or len(test)<2: continue
        j=f["best"]; amp=f["amplitudes"][j]
        n=sum(float(b["y"]@b["y"]) for b in test)
        s=sum(float(np.sum((b["y"]-amp*b["B"][:,j])**2)) for b in test)
        folds.append(dict(heldout_parity=parity,train_blocks=len(train),test_blocks=len(test),
            length_arcsec=GRID[j][0],sigma_arcsec=GRID[j][1],
            relative_improvement=float((n-s)/n) if n>0 else None))
    # Retain shapes within 10% of the maximal improvement over null.
    # This is deliberately NOT a statistical confidence set.
    near=np.flatnonzero(loss<=loss[best]+.1*improvement+1e-15)
    rows=[dict(length_arcsec=l,sigma_arcsec=s,amplitude_at_17arcsec_nMgy_pix=float(fit["amplitudes"][j]),
        loss_per_block=float(loss[j]/len(blocks)),near_best=bool(j in near)) for j,(l,s) in enumerate(GRID)]
    return dict(best_index=best,best=rows[best],grid=rows,block_count=len(blocks),
        null_loss_per_block=null/len(blocks),relative_improvement=improvement/null if null>0 else 0.,
        near_best_length_range_arcsec=[min(GRID[j][0] for j in near),max(GRID[j][0] for j in near)],
        near_best_sigma_range_arcsec=[min(GRID[j][1] for j in near),max(GRID[j][1] for j in near)],
        folds=folds)


def positive_cv(result):
    return bool(result and len(result["folds"])==2 and
                all(f["relative_improvement"] is not None and f["relative_improvement"]>0 for f in result["folds"]))


def injection_trials(blocks,result):
    """Conditional recovery on observed residual minus fitted grid component.
    Inject into projected core vectors (linear flank operation commutes);
    this does not repeat halo fitting or calibrate independent uncertainties.
    """
    if not result or len(blocks)<4: return []
    j=result["best_index"]
    amp=result["best"]["amplitude_at_17arcsec_nMgy_pix"]
    if amp<=0: return []
    out=[]
    for length,sigma in ((5.,.524),(20.,.786),(60.,1.048)):
        k=GRID.index((length,sigma))
        injected=[dict(b,y=b["y"]-amp*b["B"][:,j]+amp*b["B"][:,k]) for b in blocks]
        r=evaluate(injected)
        out.append(dict(injected_length_arcsec=length,injected_sigma_arcsec=sigma,
            injected_amplitude_at_17arcsec_nMgy_pix=amp,recovered=r["best"],
            recovered_near_best_length_range_arcsec=r["near_best_length_range_arcsec"],
            recovered_near_best_sigma_range_arcsec=r["near_best_sigma_range_arcsec"],
            recovered_positive_cv=positive_cv(r)))
    return out


def decay_experiment(residual,safe,*,source_xy,target_xy,scale,preview,refit_residual):
    from .stellar_halo_spikes import coords,scenario_sums
    if residual.ndim!=2 or safe.shape!=residual.shape or not np.isfinite(scale) or scale<=0:
        raise ValueError("SH1-E3 invalid aligned input")
    t,q=coords(residual.shape,source_xy,target_xy,scale)
    blocks,coverage=prepare_blocks(residual,safe,t,q,scale)
    fit=evaluate(blocks)
    controls=[]
    for off in OFFSETS:
        cb,_=prepare_blocks(residual,safe,t,q,scale,off)
        ids=set(b["index"] for b in cb)&set(b["index"] for b in blocks)
        f=evaluate([b for b in blocks if b["index"] in ids])
        g=evaluate([b for b in cb if b["index"] in ids])
        comparable=bool(f and g and len(ids)>=4 and len(f["folds"])==len(g["folds"])==2)
        def score(r):
            v=[x["relative_improvement"] for x in r["folds"]]
            return min(v) if v and None not in v else None
        fs,gs=(score(f),score(g)) if comparable else (None,None)
        controls.append(dict(offset_arcsec=off,matched_block_indices=sorted(ids),
            target_min_cv_improvement=fs,control_min_cv_improvement=gs,
            target_beats_control=bool(fs is not None and gs is not None and fs>max(0.,gs))))
    usable=[c for c in controls if c["target_min_cv_improvement"] is not None and c["control_min_cv_improvement"] is not None]
    passed=positive_cv(fit) and len(usable)>=3 and all(c["target_beats_control"] for c in usable)
    status="SH1E3_CONDITIONAL_GRID_NOT_ADOPTED" if passed else "SH1E3_INSUFFICIENT_PREDICTIVE_OR_CONTROL_EVIDENCE"
    out=dict(schema="desidera.sh1e3-fixed-axis-decay/1",status=status,science_adopted=False,
        fixed_native_axis_deg=0.,reference_arcsec=ANCHOR,coverage=coverage,fit=fit,
        control_comparisons=controls,usable_control_count=len(usable),positive_cv=positive_cv(fit),
        injections=injection_trials(blocks,fit),layers={},
        warning="Exploratory Gaussian transverse/exponential longitudinal grid; length is e-folding scale, not physical endpoint. Control blocks share exposure; no calibrated significance or confidence intervals. Injection is conditional on fixed halo/sky and observed mask. Mirroring is unvalidated.")
    if not passed: return out
    p=preview; original=p["layers"]["original"]; base=p["layers"]["residual"]
    ot,oq=coords(original.shape,p["source_xy_native_pix"],p["target_xy_native_pix"],scale)
    f=fit["best"]; l=np.abs(ot); sigma=f["sigma_arcsec"]
    # The finite fitting domain is explicit. Model predictions across
    # missing away blocks are conditional parametric extrapolations.
    # No prediction in unmeasured central region or beyond search range.
    m=f["amplitude_at_17arcsec_nMgy_pix"]*basis(l,oq,f["length_arcsec"],sigma)
    strip=np.abs(oq)<=5.
    common=np.isfinite(original)&np.isfinite(base)&np.isfinite(refit_residual)
    m=np.where(strip,m,0.)
    common &= ~(strip&((l<EDGES[0])|(l>=EDGES[-1])))
    m=np.where(common,m,np.nan)
    r=np.where(common,refit_residual,np.nan); corrected=r-m
    for v in (m,corrected): v.setflags(write=False)
    out["layers"]={"decay_spike_model":m,"decay_spike_corrected":corrected}
    out["aperture_and_ring_sums"]=scenario_sums(original,base,r,m,strip,p["target_xy_native_pix"],scale)
    return out
