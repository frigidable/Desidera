"""SH1-D: signed sector profiles and source-mask sensitivity, no adoption."""
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

EDGES=(0.,5.,10.,15.,20.,25.,30.,35.,40.,45.)
SECTORS=("N","S","E","W","NE","NW","SE","SW","ALL")


def _samples(p,wcs,common=None):
    from astropy.wcs.utils import wcs_to_celestial_frame
    original=p["layers"]["original"]
    residual=p["layers"]["residual"]
    tx,ty=p["target_xy_native_pix"]; sx,sy=p["source_xy_native_pix"]
    scale=p["pixel_scale_arcsec"]
    yy,xx=np.ogrid[:original.shape[0],:original.shape[1]]
    rt=np.hypot(xx-tx,yy-ty)*scale
    y,x=np.nonzero(rt<=EDGES[-1])
    frame=wcs_to_celestial_frame(wcs.celestial)
    ra,dec=wcs.all_pix2world(x,y,0)
    ra0,dec0=wcs.all_pix2world(tx,ty,0)
    center=SkyCoord(ra0*u.deg,dec0*u.deg,frame=frame).icrs
    points=SkyCoord(ra*u.deg,dec*u.deg,frame=frame).icrs
    east,north=center.spherical_offsets_to(points)
    n=north.arcsec>=0; e=east.arcsec>=0
    groups=(n,~n,e,~e,n&e,n&~e,~n&e,~n&~e,np.ones(x.size,dtype=bool))
    valid=np.isfinite(residual[y,x])&np.isfinite(original[y,x])
    if common is not None:
        if common.shape!=original.shape:
            raise ValueError("SH1-D common mask shape mismatch")
        valid &= common[y,x]
    return dict(x=x,y=y,rt=rt[y,x],groups=groups,valid=valid,
        raw=original[y,x],ref=residual[y,x],scale=scale,dx=(x-sx)*scale,dy=(y-sy)*scale)


def _summarize(s,candidates):
    rows=[]; selections=[]
    for lo,hi in zip(EDGES[:-1],EDGES[1:]):
        ann=(s["rt"]>=lo)&(s["rt"]<hi if hi!=EDGES[-1] else s["rt"]<=hi)
        for sector,group in zip(SECTORS,s["groups"]):
            geom=ann&group; sel=geom&s["valid"]; n=int(sel.sum())
            selections.append(sel)
            rows.append(dict(r_inner_arcsec=lo,r_outer_arcsec=hi,sector=sector,
                in_frame_pixels=int(geom.sum()),used_pixels=n,
                raw_flux_nMgy=float(s["raw"][sel].sum()) if n else None,
                reference_residual_flux_nMgy=float(s["ref"][sel].sum()) if n else None,
                candidate_residual_fluxes_nMgy=[]))
    dx,dy=s["dx"],s["dy"]; r=np.hypot(dx,dy)
    for c in candidates:
        h=np.interp(r,c["radial_knots_arcsec"],c["halo_at_knots_nMgy_pix"],left=np.nan,right=0.)
        b=c["sky_constant_nMgy_pix"]+c["sky_x_nMgy_pix_per_arcsec"]*dx+c["sky_y_nMgy_pix_per_arcsec"]*dy
        v=s["raw"]-h-b
        if not np.isfinite(v[s["valid"]]).all():
            raise ValueError("SH1-D model support mismatch")
        for row,sel in zip(rows,selections):
            row["candidate_residual_fluxes_nMgy"].append(float(v[sel].sum()) if row["used_pixels"] else None)
    for row in rows:
        n=row["used_pixels"]; vals=row["candidate_residual_fluxes_nMgy"]
        row["mean_intensity_range_nMgy_arcsec2"]=(
            [min(vals)/(n*s["scale"]**2),max(vals)/(n*s["scale"]**2)] if n else None)
        row["used_fraction_in_frame"]=n/max(1,row["in_frame_pixels"])
    return rows


def sector_diagnostics(p,comparison,wcs,entries):
    s=_samples(p,wcs); tx,ty=p["target_xy_native_pix"]; scale=s["scale"]
    sources=[]
    for e in entries:
        if not e.get("centre_in_native_frame") or e.get("x") is None or e.get("y") is None:
            continue
        x,y=float(e["x"]),float(e["y"])
        if not np.isfinite([x,y]).all() or np.hypot(x-tx,y-ty)*scale>49.:
            continue
        ex,ey=wcs.all_world2pix(float(e["ra_deg"]),float(e["dec_deg"]),0)
        if np.hypot(ex-x,ey-y)>.75:
            raise ValueError("SH1-D catalogue WCS mismatch")
        sources.append(dict(source_id=e["source_id"],objid=e["objid"],x=x,y=y,
            distance_arcsec=float(np.hypot(x-tx,y-ty)*scale),
            tractor_type=e.get("tractor_type"),gaia_dr3_source_id=e.get("gaia_dr3_source_id"),
            foreground_evidence_assessment=e.get("foreground_evidence_assessment"),
            evidence_flags=e.get("evidence_flags",[])))
    sources.sort(key=lambda e:e["distance_arcsec"])
    trials=[]
    for e in sources:
        d=np.hypot(s["x"]-e["x"],s["y"]-e["y"])*scale
        for mr in (2.,4.):
            for a in (20.,30.,35.,40.,45.):
                base=s["valid"]&(s["rt"]<=a)
                removed=base&(d<=mr); kept=base&~removed
                n=int(base.sum()); nr=int(removed.sum())
                trials.append(dict(objid=e["objid"],mask_radius_arcsec=mr,aperture_radius_arcsec=a,
                    base_pixels=n,removed_pixels=nr,remaining_pixels=int(kept.sum()),
                    baseline_reference_flux_nMgy=float(s["ref"][base].sum()) if n else None,
                    excluded_reference_flux_nMgy=float(s["ref"][kept].sum()) if kept.any() else None,
                    removed_signed_flux_nMgy=float(s["ref"][removed].sum()) if n else None))
    return dict(schema="desidera.sh1d-sector-diagnostic/1",status="SH1D_NOT_ADOPTED",
        science_adopted=False,geometry="NATIVE_CIRCULAR_RINGS_ICRS_PIXEL_CENTRE_SECTORS",
        radial_edges_arcsec=list(EDGES),rows=_summarize(s,comparison["candidates"]),
        sources=sources,individual_source_exclusion_trials=trials,source_exclusions_reference_only=True,
        warning="Signed matched-pixel fluxes; no fill, PSF matching or area correction. Individual mask trials remove galaxy light as well as possible contaminants. No membership inference.")


def paired_sector_diagnostics(g,r,wcs):
    """Same-session snapshots; exact native grid only."""
    for key in ("project_id","source_id","tractor_catalog_sha256","draft_digest",
                "original_scene_grid_fingerprint","sh1d_policy"):
        if g.get(key) is None or g.get(key)!=r.get(key):
            raise ValueError("SH1-D pair mismatch: "+key)
    if g["band"]!="g" or r["band"]!="r":
        raise ValueError("SH1-D requires g and r")
    gp,rp=g["sh1c_preview"],r["sh1c_preview"]
    for key in ("target_xy_native_pix","source_xy_native_pix","pixel_scale_arcsec"):
        if not np.allclose(gp[key],rp[key],rtol=0,atol=1e-6):
            raise ValueError("SH1-D pair geometry mismatch: "+key)
    ga,ra=gp["layers"]["residual"],rp["layers"]["residual"]
    if ga.shape!=ra.shape:
        raise ValueError("SH1-D pair native shape mismatch")
    common=np.isfinite(ga)&np.isfinite(ra)&np.isfinite(gp["layers"]["original"])&np.isfinite(rp["layers"]["original"])
    bands={}
    for band,res,p in (("g",g,gp),("r",r,rp)):
        rows=_summarize(_samples(p,wcs,common),res["sh1_sector_measurement"]["sh1b_model_comparison"]["candidates"])
        bands[band]=dict(rows=rows,original_fits_sha256_local=res["original_fits_sha256_local"],
                         effective_source_mask_sha256=p["effective_source_mask_sha256"])
    if any(a["used_pixels"]!=b["used_pixels"] for a,b in zip(bands["g"]["rows"],bands["r"]["rows"])):
        raise ValueError("SH1-D common support mismatch")
    import hashlib
    return dict(schema="desidera.sh1d-paired-sector-diagnostic/1",
        status="SH1D_PAIRED_NATIVE_PIXELS_NOT_ADOPTED",science_adopted=False,
        grid_fingerprint=g["original_scene_grid_fingerprint"],draft_digest=g["draft_digest"],
        common_mask_sha256=hashlib.sha256(np.ascontiguousarray(common,dtype=np.uint8).tobytes()).hexdigest(),
        common_mask_shape=list(common.shape),common_pixel_count=int(common.sum()),
        bands=bands,psf_matched=False,
        warning="Same-session snapshots, original hashes at measurement; no fresh rehash of previous band. Exact native WCS/grid only. No PSF homogenization, adopted colour, or common total flux.")
