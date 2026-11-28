"""Read-only CONTAM-QA inventory from an existing g/r diagnostics archive.
No source subtraction, mask edits, science adoption, or FITS access.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import zipfile


def ellipse_coordinate(ra, dec, geometry):
    a,d,a0,d0=map(math.radians, (ra,dec,*geometry['centre_icrs_deg']))
    da=a-a0
    den=math.sin(d)*math.sin(d0)+math.cos(d)*math.cos(d0)*math.cos(da)
    if den<=0: raise ValueError('Source outside tangent hemisphere')
    scale=180/math.pi*3600/den
    e=scale*math.cos(d)*math.sin(da)
    n=scale*(math.sin(d)*math.cos(d0)-math.cos(d)*math.sin(d0)*math.cos(da))
    t=math.radians(geometry['pa_deg'])
    return math.hypot(e*math.sin(t)+n*math.cos(t),(e*math.cos(t)-n*math.sin(t))/geometry['q'])


def audit(path, radius=120.):
    with zipfile.ZipFile(path) as z:
        if len(z.namelist())!=len(set(z.namelist())): raise ValueError('Duplicate archive member')
        load=lambda n:json.loads(z.read(n))
        manifest=load('manifest.json')
        outer=load('outer-qa.json')
        scenes={b:load(b+'-qa.json')['scene'] for b in ('g','r')}
    geom=outer['reference_geometry']
    if not (0<geom['q']<=1 and math.isfinite(radius) and radius>0): raise ValueError('Invalid geometry or radius')
    inv={b:s['source_inventory'] for b,s in scenes.items()}
    if any(v['truncated'] for v in inv.values()): raise ValueError('Inventory truncated; cannot claim complete audit')
    if scenes['g']['catalog_sha256']!=scenes['r']['catalog_sha256']: raise ValueError('Catalogue mismatch')
    sources={}
    for b,v in inv.items():
        seen=set()
        for e in v['entries']:
            sid=e['source_id']
            if sid in seen: raise ValueError('Duplicate source ID')
            seen.add(sid)
            sources.setdefault(sid,{})[b]=e
    rows=[]
    for sid,entries in sources.items():
        e=entries.get('g',entries.get('r'))
        if any((x['ra_deg'],x['dec_deg'])!=(e['ra_deg'],e['dec_deg']) for x in entries.values()): raise ValueError('Source coordinate mismatch')
        rho=ellipse_coordinate(e['ra_deg'],e['dec_deg'],geom)
        protected=any(x['screening_disposition']=='PROTECTED_OVERLAP' for x in entries.values())
        if rho>radius and not protected: continue
        flux=e['model_fluxes_nMgy']; g,r=flux.get('g'),flux.get('r')
        gaia=any(x['foreground_evidence_assessment']=='GAIA_ASTROMETRY_FOREGROUND_CANDIDATE_REVIEW' for x in entries.values())
        flags=sorted(set(f for x in entries.values() for f in x.get('maskbits_names',[])))
        halo=e['objid']==manifest['foreground_objid']
        if halo: review='EXISTING_SH1_BRIGHT_STAR'
        elif gaia: review='GAIA_FOREGROUND_CANDIDATE'
        elif e['tractor_type']=='PSF': review='POINT_SOURCE_OR_BLEND_UNRESOLVED'
        else: review='EXTENDED_TARGET_OR_BACKGROUND_UNRESOLVED'
        if rho<=40: zone='CORE_0_40'
        elif rho<=65: zone='OUTER_40_65'
        elif rho<=radius: zone='EXTENSION_65_LIMIT'
        else: zone='PROTECTION_OVERLAP_ONLY'
        rows.append(dict(source_id=sid,objid=e['objid'],ra_deg=e['ra_deg'],dec_deg=e['dec_deg'],elliptical_radius_arcsec=rho,zone=zone,tractor_type=e['tractor_type'],review_class=review,protected_overlap=protected,inventory_bands=''.join(sorted(entries)),g_model_nMgy=g,r_model_nMgy=r,model_g_minus_r=-2.5*math.log10(g/r) if g is not None and r is not None and g>0 and r>0 else None,gaia_dr3_source_id=e.get('gaia_dr3_source_id'),shape_r_arcsec=e.get('shape_r_arcsec'),maskbits='|'.join(flags),g_disposition=entries.get('g',{}).get('screening_disposition'),r_disposition=entries.get('r',{}).get('screening_disposition'),g_action=entries.get('g',{}).get('photometry_action'),r_action=entries.get('r',{}).get('photometry_action'),decision='REVIEW_NOT_ADOPTED'))
    rows.sort(key=lambda x:(x['elliptical_radius_arcsec'],x['source_id']))
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('archive',type=Path);p.add_argument('--radius',type=float,default=120.)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();rows=audit(a.archive,a.radius)
    if not rows: raise SystemExit('No relevant sources in inventory')
    with a.output.open('x',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print(f'{len(rows)} catalogue entries; REVIEW_NOT_ADOPTED; {a.output}')


if __name__=='__main__': main()
