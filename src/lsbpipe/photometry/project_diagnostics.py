"""Explicit headless (no windows) orchestration of existing experimental stages.

Uses the same Scene/Sky worker functions as Qt. Requires GUI extra because
those worker entry points currently live in GUI modules. No QApplication,
network, new branches, adopted products or original FITS changes.
"""
from __future__ import annotations
from datetime import datetime,timezone
import importlib.metadata as md
import json
from pathlib import Path
import platform
import struct
import zlib
import numpy as np

from .diagnostic_cache import atomic_zip,cached,digest,file_hash

RECIPE="desidera.project-diagnostics/1"
SKY_ONLY={"background.py","photometry/residual_sky_diagnostic.py","gui/sky_workspace.py"}
VALIDATION_ONLY={"photometry/sky_clean_injection.py","photometry/sersic_recovery.py"}


def code_identity():
    root=Path(__file__).resolve().parent.parent
    files={}
    for p in sorted(root.rglob("*.py")):
        files[p.relative_to(root).as_posix()]=file_hash(p)
    versions={"python":platform.python_version()}
    for package in ("numpy","scipy","astropy","PySide6","photutils"):
        try: versions[package]=md.version(package)
        except md.PackageNotFoundError: versions[package]=None
    # Conservative SH1 dependency set: all Python modules except the three
    # explicitly downstream sky modules. Other code edits can invalidate SH1.
    upstream={k:v for k,v in files.items() if k not in SKY_ONLY|VALIDATION_ONLY}
    sky={k:v for k,v in files.items() if k not in VALIDATION_ONLY}
    return {"sh1":digest([upstream,versions]),"sky":digest([sky,versions]),
            "validation":digest([files,versions]),"versions":versions}


def project_identity(project,draft,width,anchor,policy):
    return dict(project_id=project.project_id,source=list(project.selected_dataset[:3]),
        data_root=str(project.data_root.resolve()),field=[float(project.field.ra_deg),float(project.field.dec_deg),float(project.field.width_arcmin)],
        target_radius_arcsec=float(project.target_radius_arcsec),draft_digest=draft.digest,
        width_arcmin=float(width),anchor_icrs=list(anchor),policy=policy)


def stage_keys(identity,ds,cont,entry,code,upstream_artifact=None):
    from .sky_masks import native_boolean_mask_sha256
    base=dict(recipe=RECIPE,project=identity,band=ds["band"],fits=ds["original_fits_sha256_local"],
        hdu_indices=ds["hdu_indices"],grid=ds["grid_fingerprint"],catalog=cont["catalog_sha256"],
        source=entry["source_id"],objid=int(entry["objid"]),
        masks={k:native_boolean_mask_sha256(cont[k]) for k in
            ("mask","protected","hard","unresolved_target_overlap_sky_veto_mask")})
    sh1=digest([base,code["sh1"]])
    sky=digest([sh1,upstream_artifact,code["sky"]]) if upstream_artifact else None
    return sh1,sky


def verify_sources(project,ds,cont):
    from ..gui.science_workflows import science_data_root
    from ..surveys.legacy import BrickHit
    from ..surveys.dataset_integrity import _expected_paths,_sha256_file
    release,region,brick,band=ds["source_quartet_key"]
    root=science_data_root(project.data_root,hit=BrickHit(release,region,brick,0.,0.,0.,0.,{band:1}),band=band)
    hashes={k:_sha256_file(v,data_root=root)["sha256_local"] for k,v in _expected_paths(root,(release,region,brick,band)).items()}
    if hashes!=ds["original_fits_sha256_local"]: raise RuntimeError("Original FITS changed during diagnostics")
    catalog=Path(cont["catalog_path"])
    if _sha256_file(catalog,data_root=project.data_root)["sha256_local"]!=cont["catalog_sha256"]:
        raise RuntimeError("Tractor catalogue changed during diagnostics")


def png_bytes(a,band,limits=None):
    from ..gui.fits_raster import rgb_for_layer
    a=np.asarray(a,dtype=float)
    v=a[np.isfinite(a)]
    if limits is None: limits=tuple(np.quantile(v,[.01,.99])) if v.size else (0.,1.)
    if limits[1]<=limits[0]: limits=(limits[0],limits[0]+1e-10)
    rgb=rgb_for_layer(a,"image",band,contrast="faint",scaling="asinh",limits=limits)
    rgb[~np.isfinite(a)]=(135,49,139)
    def chunk(kind,data):
        return struct.pack(">I",len(data))+kind+data+struct.pack(">I",zlib.crc32(kind+data)&0xffffffff)
    raw=b"".join(b"\0"+row.tobytes() for row in rgb[::-1])
    return b"\x89PNG\r\n\x1a\n"+chunk(b"IHDR",struct.pack(">IIBBBBB",a.shape[1],a.shape[0],8,2,0,0,0))+chunk(b"IDAT",zlib.compress(raw))+chunk(b"IEND",b""),list(map(float,limits))


def compact_summary(bands,pairs):
    result={}
    for band,v in bands.items():
        sky=v["sky"]["expanded_sky"]
        rows=[]
        for t in sky["trials"]:
            def scores(candidates):
                return [dict(method=c["method"],status=c["status"],
                    zero_mae=c.get("zero_equal_quadrant_mae_nMgy_pix"),
                    cv_mae=c.get("cv_equal_quadrant_mae_nMgy_pix")) for c in candidates]
            rows.append(dict(halo=t["halo_scenario"],guard_arcsec=t["target_guard_arcsec"],
                tile_pix=t["tile_size_pix"],models=scores(t["candidates"]),
                edge_stress=[dict(side=e["side"],removed_blocks=e["removed_blocks"],models=scores(e["candidates"]))
                    for e in t.get("edge_stress",[])]))
        p=v["probe"]["sh1c_preview"]
        result[band]=dict(geometry=(p.get("geo1") or {}).get("reference"),
            spike_decay_status=((p.get("sh1e") or {}).get("decay_extent") or {}).get("status"),
            sky=rows)
        edge=sky.get("scene_edge") or {}
        if edge:
            result[band]["scene_edge"]={k:edge[k] for k in (
                "primary_safe_pixels","inclusive_safe_pixels","removed_safe_pixels","removed_right25_pixels",
                "removed_signed_residual_flux_nMgy","removed_right25_signed_residual_flux_nMgy","warning")
                if k in edge}
            meta=edge.get("scene") or {}
            result[band]["scene_edge"]["scene"]={k:v for k,v in meta.items() if k!="nonprimary_entries"}
            result[band]["scene_edge"]["controls"]={k:dict(block_count=v["block_count"],
                models=scores(v["candidates"]),brightest_blocks=v["brightest_blocks"])
                for k,v in edge.get("controls",{}).items()}
    return dict(bands=result,paired={k:{"status":v.get("status"),
        **{x:v[x] for x in ("reason","common_pixels","common_pixel_count","common_mask_sha256") if x in v}}
        for k,v in pairs.items()})


def report_archive(path,manifest,bands,pairs,validation=None,inj2=None):
    from ..gui.scene_editor import _scene_qa_json_value
    def write(z):
        images=[]
        for band,b in bands.items():
            p=b["probe"]["sh1c_preview"]; sky=b["sky"]["expanded_sky"]
            payload=dict(scene=b["scene"],probe=b["probe"],sky=b["sky"])
            z.writestr(f"{band}-qa.json",json.dumps(_scene_qa_json_value(payload),ensure_ascii=False,allow_nan=False,separators=(",",":")))
            layers={"native-original":p["layers"]["original"],"native-sh1c":p["layers"]["residual"]}
            layers.update({f"expanded-{k}":v for k,v in sky["expanded_layers"].items()})
            if "expanded_layers_edge_added" in sky:
                layers["expanded-edge-newly-excluded"]=sky["expanded_layers_edge_added"].astype(float)
            layers.update({k:v for k,v in sky["layers"].items()})
            t=next((t for t in sky["trials"] if t["halo_scenario"]=="SH1C" and t["target_guard_arcsec"]==60 and t["tile_size_pix"]==32),None)
            if t:
                for stat in ("mean","median"):
                    arr=np.full(sky["expanded_shape"],np.nan)
                    for block in t["blocks"]:
                        x,y=block["tile_x_pix"],block["tile_y_pix"]; tile=t["tile_size_pix"]
                        arr[y:y+tile,x:x+tile]=block[stat+"_nMgy_pix"]
                    layers["expanded-block-"+stat]=arr
            shared=p["layers"]["original"]; v=shared[np.isfinite(shared)]
            limits=tuple(np.quantile(v,[.01,.99])) if v.size else (0.,1.)
            for name,a in layers.items():
                same=name.startswith("native-") or name.endswith("_corrected")
                display=(0.,1.) if name=="expanded-edge-newly-excluded" else None
                if "block-" in name:
                    v=a[np.isfinite(a)]
                    lim=max(float(np.quantile(abs(v),.98)) if v.size else 0.,1e-12)
                    display=(-lim,lim)
                data,bounds=png_bytes(a,band,limits if same else display)
                filename=f"{band}-{name}.png"; z.writestr(filename,data)
                images.append(dict(file=filename,limits_nMgy_pix=bounds,orientation="FITS_Y_UP",purple="unavailable"))
        if validation is not None:
            z.writestr("validation-qa.json",json.dumps(_scene_qa_json_value(validation),ensure_ascii=False,allow_nan=False,separators=(",",":")))
            for band,masks in validation.get("masks",{}).items():
                if band=="common_removed": masks={"common_removed":masks}
                for name,mask in masks.items():
                    filename=f"validation-{band}-{name}.png"
                    data,bounds=png_bytes(mask.astype(float),"g",(0.,1.)); z.writestr(filename,data)
                    images.append(dict(file=filename,limits_nMgy_pix=bounds,orientation="FITS_Y_UP",meaning="1=additional sky-only veto; not a photometric mask"))
        if inj2 is not None:
            z.writestr("inj2-qa.json",json.dumps(_scene_qa_json_value(inj2),ensure_ascii=False,allow_nan=False,separators=(",",":")))
            for name,arr in inj2.get("target",{}).get("layers",{}).items():
                filename="inj2-"+name+".png"
                data,bounds=png_bytes(arr,name[0]);z.writestr(filename,data)
                images.append(dict(file=filename,limits_nMgy_pix=bounds,orientation="FITS_Y_UP",purple="unavailable"))
        manifest["images"]=images
        z.writestr("manifest.json",json.dumps(manifest,ensure_ascii=False,allow_nan=False,indent=2))
        z.writestr("paired-qa.json",json.dumps(_scene_qa_json_value(pairs),ensure_ascii=False,allow_nan=False,separators=(",",":")))
        summary=compact_summary(bands,pairs)
        if validation is not None:
            summary["validation"]={k:v for k,v in validation.items() if k not in ("bands","masks","target")}
            if validation.get("target"):
                target=validation["target"]
                summary["validation"]["target"]={k:v for k,v in target.items() if k!="profiles"}
                summary["validation"]["target"]["profiles_all"]=[r for r in target.get("profiles",[]) if r["sector"]=="ALL"]
            summary["validation"]["bands"]={}
            for band,b in validation.get("bands",{}).items():
                v={k:x for k,x in b.items() if k not in ("before","after")}
                for key in ("before","after"):
                    if key not in b: continue
                    v[key]={k:b[key][k] for k in ("block_count","empty_aperture_count")}
                    v[key]["models"]=[{k:c.get(k) for k in ("method","status","zero_equal_quadrant_mae_nMgy_pix","cv_equal_quadrant_mae_nMgy_pix")} for c in b[key]["candidates"]]
                summary["validation"]["bands"][band]=v
        if inj2 is not None:
            summary["inj2"]={k:v for k,v in inj2.items() if k!="target"}
            if "target" in inj2:
                summary["inj2"]["target"]={k:v for k,v in inj2["target"].items() if k!="layers"}
        z.writestr("summary.json",json.dumps(_scene_qa_json_value(summary),ensure_ascii=False,allow_nan=False,indent=2))
    atomic_zip(path,write)


def run_project(root,*,bands=("g","r"),objid,force_stage=None,width=None,anchor=None,
                gaia_inner=False,red_extended=False,progress=print):
    from ..surveys.object_projects import read_object_project
    from ..surveys.project_locking import exclusive_metadata_lock
    from .mask_drafts import read_active_mask_draft
    from ..gui.scene_editor import load_scene_original_cutout,prepare_scene_neighbour_masks,run_scene_psf0_foreground_probe
    from ..gui.sky_workspace import run_expanded_sky
    from .stellar_halo_sectors import paired_sector_diagnostics
    from .galaxy_geometry_diagnostic import paired_geometry
    from .residual_sky_diagnostic import paired_sky_diagnostics
    from .sky_clean_injection import run_validation
    from .sersic_recovery import psf_evidence,run_inj2
    project=read_object_project(root)
    if project.field is None or project.selected_dataset is None or project.data_root is None:
        raise ValueError("Project needs local native source, field and data root")
    if not bands or len(set(bands))!=len(bands) or any(b not in ("g","r","i","z") for b in bands):
        raise ValueError("Choose distinct optical bands")
    if force_stage not in (None,"sh1","sky","all"): raise ValueError("Unknown force stage")
    draft,_=read_active_mask_draft(root)
    width=float(width if width is not None else min(6.,project.field.width_arcmin))
    anchor=tuple(anchor if anchor is not None else (project.field.ra_deg,project.field.dec_deg))
    policy=dict(gaia_inner_preview=bool(gaia_inner),red_extended_preview=bool(red_extended),sga_prior=None)
    identity=project_identity(project,draft,width,anchor,policy)
    code=code_identity()
    folder=project.root/"analyses"
    if folder.is_symlink(): raise ValueError("Unsafe analyses folder")
    folder.mkdir(exist_ok=True)
    manifest=dict(schema=RECIPE,generated_utc=datetime.now(timezone.utc).isoformat(),
        status="PROJECT_DIAGNOSTICS_NOT_ADOPTED",science_adopted=False,identity=identity,
        code=code,requested_bands=list(bands),foreground_objid=objid,stages={},
        limitations=["Saved mask draft only; no unsaved GUI state or SGA prior.",
            "SH1 includes GEO-1 and directional diagnostics; no fresh PSF-0 fit, calibrated sky or final photometry. INJ-1 is fixed-halo, observed-profile recovery only.",
            "Caches preserve numerical arrays; report contains metadata and display PNGs, not new science FITS."])
    products={}; datasets={}; contours={}
    with exclusive_metadata_lock(folder/"project-diagnostics.lock"):
        for band in bands:
            progress(f"{band}: verify original FITS and prepare saved Scene masks")
            ds=load_scene_original_cutout(str(project.root),band,width)
            cont=prepare_scene_neighbour_masks(ds,str(project.root),draft,anchor,
                gaia_inner_preview=gaia_inner,red_extended_preview=red_extended)
            inventory=cont["source_inventory"]
            matches=[e for e in inventory["entries"] if e["objid"]==objid and
                e.get("tractor_type")=="PSF" and e.get("gaia_dr3_source_id") is not None and
                e.get("foreground_evidence_assessment")=="GAIA_ASTROMETRY_FOREGROUND_CANDIDATE_REVIEW" and
                f"SATUR_{band.upper()}" in (e.get("maskbits_names") or ())]
            if len(matches)!=1: raise ValueError(f"{band}: explicit saturated foreground OBJID {objid} absent or ambiguous")
            entry=matches[0]; key,_=stage_keys(identity,ds,cont,entry,code)
            def probe():
                v=run_scene_psf0_foreground_probe(ds,str(project.root),entry,inventory,anchor,
                    scene_masks=cont,scene_draft=draft)
                if (v.get("status")!="SH1_SECTORS_MEASURED_NOT_ADOPTED" or
                    not (v.get("sh1c_preview") or {}).get("sky1_inputs")):
                    raise RuntimeError(f"{band}: SH1 did not produce expanded support: {v.get('status')}")
                return v
            progress(f"{band}: SH1 + GEO-1 (cache check)")
            pr,mode,artifact=cached(folder/f"project-diagnostics-{band}-sh1.zip",key,probe,
                force=force_stage in ("sh1","all"))
            progress(f"{band}: SH1 {mode}")
            _,skykey=stage_keys(identity,ds,cont,entry,code,artifact)
            progress(f"{band}: SKY-2 (cache check)")
            sk,sm,sa=cached(folder/f"project-diagnostics-{band}-sky.zip",skykey,
                lambda:run_expanded_sky(ds,draft,cont,pr),force=force_stage in ("sky","sh1","all"))
            progress(f"{band}: SKY-2 {sm}")
            products[band]=dict(probe=pr,sky=sk,scene=cont)
            datasets[band]=ds; contours[band]=cont
            manifest["stages"][band]=dict(sh1=dict(cache=mode,key=key,archive_sha256=artifact,status=pr["status"]),
                sky=dict(cache=sm,key=skykey,archive_sha256=sa,status=sk["status"]),
                fits_sha256=ds["original_fits_sha256_local"],catalog_sha256=cont["catalog_sha256"],
                reference_geometry=(pr["sh1c_preview"].get("geo1") or {}).get("reference"))
        pairs={}
        if {"g","r"}<=set(products):
            g,r=products["g"],products["r"]
            for name,fn,args in (("sectors",paired_sector_diagnostics,(g["probe"],r["probe"],datasets["g"]["wcs"])),
                    ("geometry",paired_geometry,(g["probe"],r["probe"],datasets["g"]["wcs"])),
                    ("sky",paired_sky_diagnostics,(g["sky"],r["sky"]))):
                try: pairs[name]=fn(*args)
                except ValueError as exc: pairs[name]=dict(status="PAIR_BLOCKED",reason=str(exc),science_adopted=False)
        progress("SKY-CLEAN + INJ-1 (cache check)")
        vk=digest([RECIPE,code.get("validation",code["sky"]),list(bands),
            {b:manifest["stages"][b]["sky"]["archive_sha256"] for b in bands}])
        validation,vm,va=cached(folder/"project-diagnostics-validation.zip",vk,
            lambda:run_validation(project,products,datasets,pairs),force=force_stage is not None)
        manifest["validation"]=dict(cache=vm,key=vk,archive_sha256=va,status=validation["status"])
        progress(f"SKY-CLEAN + INJ-1 {vm}: {validation['status']}")
        progress("INJ-2 + conditional Sersic totals (cache check)")
        kernels,psf=psf_evidence(project,products,datasets)
        ik=digest([RECIPE,code.get("validation",code["sky"]),va,psf])
        inj2,im,ia=cached(folder/"project-diagnostics-inj2.zip",ik,
            lambda:run_inj2(products,datasets,validation,kernels,psf),force=force_stage is not None)
        manifest["inj2"]=dict(cache=im,key=ik,archive_sha256=ia,status=inj2["status"])
        progress(f"INJ-2 {im}: {inj2['status']}")
        # No completed report if inputs changed mid-run. The cache keys
        # still describe the old input snapshot and cannot match new data.
        current=read_object_project(root); current_draft,_=read_active_mask_draft(root)
        if project_identity(current,current_draft,width,anchor,policy)!=identity or code_identity()!=code:
            raise RuntimeError("Project, saved draft or code changed during diagnostics; rerun")
        for band in bands: verify_sources(current,datasets[band],contours[band])
        manifest["paired_status"]={k:v["status"] for k,v in pairs.items()}
        # Replace only the diagnostic report; earlier scientific artifacts untouched.
        path=folder/"project-diagnostics-qa.zip"
        report_archive(path,manifest,products,pairs,validation,inj2)
    progress(f"QA: {path}")
    return path
