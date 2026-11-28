"""SH1-E2: contiguous directional excess, empirical control-strip screening.

Detection is exploratory: control strips are not independent exposures,
and no calibrated significance, physical cutoff or target-side symmetry
is inferred. All outputs are conditional diagnostics.
"""
import numpy as np

ANGLES=(-8.,-4.,0.,4.,8.)
EDGES=(12.,17.,22.,27.,32.,37.,42.,47.,52.,57.,62.,67.,72.,77.,82.,87.,92.,97.,102.,107.,112.,117.,122.,127.,132.,137.,142.,147.,152.,157.,162.,167.,172.,177.,180.)
OFFSETS=(-30.,-20.,-10.,10.,20.,30.)
# Approx. two native pixels at Legacy 0.262 arcsec/px.
STEP=.524
MAX_WIDTH=4*STEP


def _contrast(v,t,q,lo,hi,offset,width):
    """Core minus mean of two flank medians, requiring all three regions."""
    dt=hi-lo
    # Count relative to nominal geometric area, never just surviving pixels.
    # t,q are projected native pixel centres; caller supplies scale.
    core=np.abs(q-offset)<width
    left=(q-offset>=-5.)&(q-offset<-3.)
    right=(q-offset>3.)&(q-offset<=5.)
    masks=(core,left,right)
    def med(m):
        return float(np.median(v[m])) if m.any() else None
    return [med(m) for m in masks],[int(m.sum()) for m in masks]


def _measurement(v,t,q,lo,hi,offset,width,scale):
    parts=[]; counts=[]
    for a,b in ((lo,(lo+hi)/2),((lo+hi)/2,hi)):
        m=(t>=a)&(t<b)
        vals,n=_contrast(v[m],t[m],q[m],a,b,offset,width)
        need=[max(8,int(.45*(b-a)*w/scale**2)) for w in (2*width,2.,2.)]
        parts.append(vals[0]-(vals[1]+vals[2])/2 if all(x>=y for x,y in zip(n,need)) else None)
        counts.append(n)
    return parts,counts


def scan_axis(residual,safe,t,q,scale):
    # Work on a narrow vector rather than repeatedly masking the full brick.
    m=safe&np.isfinite(residual)&(-t>=EDGES[0])&(-t<EDGES[-1])&(np.abs(q)<=36.)
    v=residual[m]; l=-t[m]; cross=q[m]
    rows=[]
    for lo,hi in zip(EDGES[:-1],EDGES[1:]):
        z=(l>=lo)&(l<hi); vv,ll,qq=v[z],l[z],cross[z]
        controls=[]
        for offset in OFFSETS:
            p,n=_measurement(vv,ll,qq,lo,hi,offset,2*STEP,scale)
            controls.append(dict(offset_arcsec=offset,split_excess_nMgy_pix=p,counts=n))
        paired=[c["split_excess_nMgy_pix"] for c in controls if None not in c["split_excess_nMgy_pix"]]
        # Max absolute control contrast per half: deliberately conservative,
        # not a p-value; no sqrt(Npix) assumption.
        threshold=[max(1e-12,max(abs(c[j]) for c in paired)) for j in (0,1)] if len(paired)>=3 else None
        widths=[]
        for k in (1,2,3,4):
            width=k*STEP
            p,n=_measurement(vv,ll,qq,lo,hi,0.,width,scale)
            # Same width controls are needed for each width.
            cs=[_measurement(vv,ll,qq,lo,hi,o,width,scale)[0] for o in OFFSETS]
            cs=[c for c in cs if None not in c]
            th=[max(1e-12,max(abs(c[j]) for c in cs)) for j in (0,1)] if len(cs)>=3 else None
            detect=bool(th is not None and None not in p and all(x>y for x,y in zip(p,th)))
            widths.append(dict(half_width_arcsec=width,split_excess_nMgy_pix=p,counts=n,
                control_count=len(cs),control_split_excesses_nMgy_pix=cs,
                control_envelope_nMgy_pix=th,detected=detect))
        # Use a fixed inner width for longitudinal detection, avoiding
        # selection of a different width solely to prolong the ray.
        fixed=widths[1]
        status=("UNMEASURED" if fixed["control_envelope_nMgy_pix"] is None or None in fixed["split_excess_nMgy_pix"]
                else "DETECTED" if fixed["detected"] else "NOT_DETECTED")
        # Signed fine transverse profile relative to local flank levels,
        # evaluated separately in each longitudinal half.
        cells=[]
        for j in range(-4,4):
            parts=[]; counts=[]
            for a,b in ((lo,(lo+hi)/2),((lo+hi)/2,hi)):
                h=(ll>=a)&(ll<b)
                regions=(h&(qq>=j*STEP)&(qq<(j+1)*STEP),
                         h&(qq>=-5)&(qq<-3),h&(qq>3)&(qq<=5))
                n=[int(r.sum()) for r in regions]
                need=[max(6,int(.45*(b-a)*w/scale**2)) for w in (STEP,2.,2.)]
                vals=[float(np.median(vv[r])) for r in regions] if all(x>=y for x,y in zip(n,need)) else None
                parts.append(vals[0]-(vals[1]+vals[2])/2 if vals is not None else None)
                counts.append(n)
            cells.append(dict(transverse_inner_arcsec=j*STEP,transverse_outer_arcsec=(j+1)*STEP,
                split_excess_nMgy_pix=parts,counts=counts))
        # Consecutive symmetric fine cells around the axis must both be
        # positive in both halves. Width is a support estimate, not FWHM.
        accepted=[]
        for k in range(4):
            pair=[cells[3-k],cells[4+k]]
            if any(None in c["split_excess_nMgy_pix"] or min(c["split_excess_nMgy_pix"])<=0 for c in pair):
                break
            accepted.extend(pair)
        # No width may be adopted beyond a detected cumulative aperture.
        detected_widths=[w["half_width_arcsec"] for w in widths if w["detected"]]
        limit=max(detected_widths,default=0.)
        accepted=[c for c in accepted if max(abs(c["transverse_inner_arcsec"]),abs(c["transverse_outer_arcsec"]))<=limit+1e-10]
        rows.append(dict(inner_arcsec=lo,outer_arcsec=hi,status=status,
            fixed_width_control_strips=controls,control_envelope_nMgy_pix=threshold,
            width_trials=widths,transverse_cells=cells,accepted_cells=accepted,
            detected_half_width_arcsec=max([max(abs(c["transverse_inner_arcsec"]),abs(c["transverse_outer_arcsec"])) for c in accepted],default=0.)))
    return rows


def contiguous_support(rows):
    # Seed near the saturated core, require two consecutive detections.
    # An unmeasured or non-detected bin stops the run; isolated far
    # detections are retained in QA but never extend the model.
    start=None
    for i in range(min(3,len(rows)-1)):
        if all(r["status"]=="DETECTED" and r["accepted_cells"] for r in rows[i:i+2]):
            start=i; break
    if start is None:
        return dict(status="NO_CONTIGUOUS_INNER_DETECTION",run=[],detected_outer_edge_arcsec=None,
                    endpoint_status="UNCONSTRAINED",next_interval_arcsec=None)
    end=start
    while end<len(rows) and rows[end]["status"]=="DETECTED" and rows[end]["accepted_cells"]:
        end+=1
    nextrow=rows[end] if end<len(rows) else None
    reason=("SEARCH_LIMIT" if nextrow is None else
            "COVERAGE_GAP" if nextrow["status"]=="UNMEASURED" else
            "TRANSVERSE_SUPPORT_GAP" if nextrow["status"]=="DETECTED" else "NEXT_BIN_NOT_DETECTED")
    return dict(status="CONTIGUOUS_DETECTION_NOT_ADOPTED",run=list(range(start,end)),
        detected_inner_edge_arcsec=rows[start]["inner_arcsec"],
        detected_outer_edge_arcsec=rows[end-1]["outer_arcsec"],endpoint_status=reason,
        next_interval_arcsec=[nextrow["inner_arcsec"],nextrow["outer_arcsec"]] if nextrow else None)


def compact_template(t,q,rows,support):
    """Conditional truncation scenario, not evidence of zero beyond run.

    Only detected fine cells contribute. Within consecutive measured
    intervals interpolate each cell's amplitude between adjacent centres
    when that same cell is supported in both. No bridging unsupported bins.
    """
    result=np.zeros(t.shape,dtype=float)
    for i in support["run"]:
        r=rows[i]; lo,hi=r["inner_arcsec"],r["outer_arcsec"]; mid=(lo+hi)/2
        for c in r["accepted_cells"]:
            a,b=c["transverse_inner_arcsec"],c["transverse_outer_arcsec"]
            value=float(np.mean(c["split_excess_nMgy_pix"]))
            for side in (-1,1):
                l=side*t; w=-side*q
                m=(l>=lo)&(l<hi)&(w>=a)&(w<b)
                if not m.any(): continue
                vv=np.full(int(m.sum()),value); ll=l[m]
                for neighbor,half in ((i-1,ll<mid),(i+1,ll>=mid)):
                    if neighbor not in support["run"]: continue
                    other=rows[neighbor]
                    match=next((x for x in other["accepted_cells"] if x["transverse_inner_arcsec"]==a),None)
                    if match is None: continue
                    omid=(other["inner_arcsec"]+other["outer_arcsec"])/2
                    ov=float(np.mean(match["split_excess_nMgy_pix"]))
                    vv[half]=value+(ov-value)*(ll[half]-mid)/(omid-mid)
                result[m]=vv
    return result


def adaptive_spike_experiment(residual,safe,*,source_xy,target_xy,scale,preview,refit_residual):
    from .stellar_halo_spikes import coords,scenario_sums
    scans=[]
    for angle in ANGLES:
        t,q=coords(residual.shape,source_xy,target_xy,scale,angle)
        rows=scan_axis(residual,safe,t,q,scale)
        support=contiguous_support(rows)
        # Rank only on the fixed near-core search interval, not ray length.
        # This is data-driven exploratory axis selection, not independent validation.
        scores=[]
        for r in rows[:3]:
            w=r["width_trials"][1]
            if None not in w["split_excess_nMgy_pix"] and w["control_envelope_nMgy_pix"] is not None:
                scores.append(min(x-y for x,y in zip(w["split_excess_nMgy_pix"],w["control_envelope_nMgy_pix"])))
        score=float(np.median(scores)) if len(scores)>=2 else None
        scans.append(dict(angle_deg=angle,axis_score_nMgy_pix=score,rows=rows,support=support))
    eligible=[s for s in scans if s["axis_score_nMgy_pix"] is not None]
    best=max(eligible,key=lambda s:(s["axis_score_nMgy_pix"],-abs(s["angle_deg"]))) if eligible else None
    out=dict(schema="desidera.sh1e2-contiguous-spike/1",science_adopted=False,
        status="NO_AXIS_SUPPORT" if best is None else best["support"]["status"],
        axis_scan=scans,selected_axis_deg=best["angle_deg"] if best else None,
        support=best["support"] if best else None,
        transverse_step_arcsec=STEP,search_edges_arcsec=list(EDGES),
        control_offsets_arcsec=list(OFFSETS),layers={},
        warning="Exploratory axis/width/extent; no calibrated significance or physical endpoint. Controls are spatial comparisons, not independent exposures. Mirrored compact model assumes zero outside the detected run, not an upper limit. No scene or FITS edits.")
    if best is None or not best["support"]["run"]:
        return out
    p=preview; original=p["layers"]["original"]; base=p["layers"]["residual"]
    t,q=coords(original.shape,p["source_xy_native_pix"],p["target_xy_native_pix"],scale,best["angle_deg"])
    m=compact_template(t,q,best["rows"],best["support"])
    common=np.isfinite(refit_residual)&np.isfinite(base)&np.isfinite(original)
    # Preserve the unsupported central strip: no extrapolation into saturated core.
    common &= ~((np.abs(t)<EDGES[0])&(np.abs(q)<MAX_WIDTH))
    m=np.where(common,m,np.nan)
    r=np.where(common,refit_residual,np.nan)
    corrected=r-m
    for a in (m,corrected): a.setflags(write=False)
    out["layers"]={"adaptive_spike_model":m,"adaptive_spike_corrected":corrected}
    out["aperture_and_ring_sums"]=scenario_sums(original,base,r,m,np.isfinite(m)&(m>0),
                                                p["target_xy_native_pix"],scale)
    return out
