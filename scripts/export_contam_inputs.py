"""Export verified numerical CONTAM-QA inputs from an existing diagnostic run.

No recomputation, network, original FITS writes or adopted products.
Run from the same installed Desidera environment used for diagnose-project.
"""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace
import zipfile

import numpy as np
from lsbpipe.photometry.diagnostic_cache import digest, file_hash, read_cache, write_cache
from lsbpipe.photometry.sersic_recovery import psf_evidence


def verified_cache(path, identity):
    expected=identity["archive_sha256"]
    if file_hash(path)!=expected:
        raise ValueError(f"Cache does not match QA snapshot: {path.name}")
    value=read_cache(path,identity["key"])
    if value is None or file_hash(path)!=expected:
        raise ValueError(f"Cache key, arrays or archive changed: {path.name}")
    return value


def extract_arrays(probe, sky):
    p=probe["sh1c_preview"]
    inputs=p["sky1_inputs"]
    s=sky["expanded_sky"]
    if inputs["expanded_grid_fingerprint"]!=s["expanded_grid_fingerprint"]:
        raise ValueError("Expanded grid mismatch")
    arrays={name:np.asarray(inputs[name]) for name in
            ("image","iv","hard","source_mask","protected","veto")}
    arrays["reference_halo_residual"]=np.asarray(s["expanded_layers"]["SH1C"])
    arrays["sky_training_safe"]=np.asarray(s["expanded_safe_mask"])
    shape=tuple(s["expanded_shape"])
    if len(shape)!=2 or any(x.shape!=shape for x in arrays.values()):
        raise ValueError("Expanded array shape mismatch")
    for name in ("hard","source_mask","protected","veto","sky_training_safe"):
        if arrays[name].dtype.kind!="b":
            raise ValueError(f"Expected boolean mask: {name}")
    return dict(arrays=arrays,expanded_grid_fingerprint=s["expanded_grid_fingerprint"],
                expanded_shape=shape,source_xy=inputs["source_xy"],
                target_xy=inputs["target_xy"],scale_arcsec=inputs["scale"],
                expanded_wcs_header=probe["sh1b2_native_window"]["wcs_fits_header"],
                expanded_window=probe["sh1b2_native_window"],
                native_reference_geometry=p["geo1"]["reference"],
                source_entries=s["source_entries"])


def export(root,output):
    root=Path(root)
    folder=root/"analyses"
    qa=folder/"project-diagnostics-qa.zip"
    output=Path(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    if output.parent.resolve()!=folder.resolve() or folder.is_symlink():
        raise ValueError("Output must be inside the existing analyses folder")
    qa_hash=file_hash(qa)
    with zipfile.ZipFile(qa) as z:
        if len(z.namelist())!=len(set(z.namelist())): raise ValueError("Duplicate QA member")
        load=lambda name:json.loads(z.read(name))
        m=load("manifest.json")
        outer=load("outer-qa.json")
        scenes={b:load(b+"-qa.json")["scene"] for b in ("g","r")}
    if not {"g","r"}<=set(m["stages"]): raise ValueError("g/r snapshot required")
    data={}; ds={}
    for b in ("g","r"):
        stage=m["stages"][b]
        scene=scenes[b]
        probe=verified_cache(folder/f"project-diagnostics-{b}-sh1.zip",stage["sh1"])
        sky=verified_cache(folder/f"project-diagnostics-{b}-sky.zip",stage["sky"])
        if (probe["original_fits_sha256_local"]!=stage["fits_sha256"] or
            sky["original_fits_sha256_local"]!=stage["fits_sha256"] or
            probe["tractor_catalog_sha256"]!=stage["catalog_sha256"] or
            sky["catalog_sha256"]!=stage["catalog_sha256"] or
            scene["catalog_sha256"]!=stage["catalog_sha256"]):
            raise ValueError(f"{b}: source provenance mismatch")
        data[b]=extract_arrays(probe,sky)
        data[b]["source_inventory"]=scene["source_inventory"]
        ds[b]={"original_fits_sha256_local":stage["fits_sha256"],
               "source_quartet_key":scene["source_quartet_key"]}
    if data["g"]["expanded_grid_fingerprint"]!=data["r"]["expanded_grid_fingerprint"]:
        raise ValueError("g/r expanded grids differ; no silent reprojection")
    kernels,meta=psf_evidence(SimpleNamespace(root=root),data,ds)
    for b in ("g","r"):
        want=outer["psf"]["bands"][b]
        got=meta["bands"].get(b,{})
        if b not in kernels or any(got.get(k)!=want.get(k) for k in
                ("kernel_sha256","lineage_id","source_hashes","kernel_shape")):
            raise ValueError(f"{b}: saved PSF no longer matches QA")
        data[b]["psf_kernel"]=kernels[b]
    if file_hash(qa)!=qa_hash: raise ValueError("QA changed during export")
    for b in ("g","r"):
        for stage in ("sh1","sky"):
            if file_hash(folder/f"project-diagnostics-{b}-{stage}.zip")!=m["stages"][b][stage]["archive_sha256"]:
                raise ValueError("Cache changed during export")
    payload=dict(schema="desidera.contam-inputs/1",science_adopted=False,
                 status="CONTAM_INPUTS_NOT_ADOPTED",qa_sha256=qa_hash,
                 stages=m["stages"],identity=m["identity"],
                 reference_geometry=outer["reference_geometry"],psf=meta,bands=data,
                 warning="Snapshot arrays only, not verified current FITS, no calibrated PSF wings, no new photometry.")
    write_cache(output,digest(["desidera.contam-inputs/1",qa_hash]),payload)
    print(output)
    return output


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--project",type=Path,required=True)
    a=p.parse_args()
    export(a.project,a.project/"analyses"/"contam-inputs.zip")


if __name__=="__main__": main()
