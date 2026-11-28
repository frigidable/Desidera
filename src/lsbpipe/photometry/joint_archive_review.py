"""Bounded read-only inspection of DESIDERA joint annular diagnostic ZIPs."""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
from pathlib import Path
import re
import struct
import zipfile

STATUS = "DIAGNOSTIC_NOT_ADOPTED"
SCHEMA = "desidera.annular-diagnostic/1"
RECEIPT_SCHEMA = "desidera.joint-archive-review/1"
HEX = re.compile(r"[0-9a-f]{64}")
FIGURES = (
    "01_stellar_halo", "02_prepared_images", "03_curve_of_growth",
    "04_surface_brightness", "05_model_residual", "06_radial_model",
    "curve_of_growth", "radial_intensity",
)
REQUIRED = {"manifest.json", "rings.csv", "arrays.npz", "fit.json"}
ALLOWED = REQUIRED | {"asymptotic.json", "publication_arrays.npz", "publication_qa.json"} | {
    f"{name}.{ext}" for name in FIGURES for ext in ("png", "svg")
}
MAX_ZIP = 512 * 1024**2
MAX_EXPANDED = 1024**3
MAX_TEXT = 16 * 1024**2
MAX_PNG = 16 * 1024**2
MAX_RECEIPT = 64 * 1024
MAX_ROWS = 20000


def _json(data):
    def unique(pairs):
        out = {}
        for k,v in pairs:
            if k in out:
                raise ValueError("Duplicate JSON key")
            out[k]=v
        return out
    def nonfinite(value):
        raise ValueError("Nonfinite JSON value")
    return json.loads(data, object_pairs_hook=unique, parse_constant=nonfinite)


def _sha_ok(value):
    return isinstance(value,str) and HEX.fullmatch(value) is not None


def _digest(stream):
    h=hashlib.sha256()
    for chunk in iter(lambda:stream.read(1024*1024),b""):
        h.update(chunk)
    return h.hexdigest()


def _context(project_id, source_triplet, target_icrs, target_radius_arcsec):
    source=tuple(source_triplet)
    target=tuple(map(float,target_icrs))
    radius=float(target_radius_arcsec)
    if (not isinstance(project_id,str) or not project_id or len(source)!=3
            or any(not isinstance(x,str) or not x for x in source)
            or len(target)!=2 or not all(map(math.isfinite,target))
            or not 0<=target[0]<360 or not -90<=target[1]<=90
            or not math.isfinite(radius) or radius<=0):
        raise ValueError("A project, native source and valid target are required")
    return {"project_id":project_id,"source_triplet":list(source),
            "target_icrs":list(target),"target_radius_arcsec":radius}


def _separation(a,b):
    r1,d1,r2,d2=map(math.radians,(*a,*b))
    h=math.sin((d2-d1)/2)**2+math.cos(d1)*math.cos(d2)*math.sin((r2-r1)/2)**2
    return math.degrees(2*math.asin(math.sqrt(min(1.,max(0.,h)))))*3600


def _png(data):
    if (len(data)<33 or data[:8]!=b"\x89PNG\r\n\x1a\n"
            or data[8:16]!=b"\x00\x00\x00\x0dIHDR"):
        raise ValueError("Invalid PNG header")
    w,h=struct.unpack(">II",data[16:24])
    if not 0<w<=10000 or not 0<h<=10000 or w*h>25_000_000:
        raise ValueError("PNG dimensions exceed review budget")
    return data


def read_joint_archive(path, *, project_id, source_triplet, target_icrs,
                       target_radius_arcsec, expected_sha256=None):
    """Read locally without extraction, pickle, network, or original FITS access.

    A matching source and nearby anchor are association checks, not proof of
    original science identity. Member hashes verify internal archive integrity,
    not authenticity or current FITS contents.
    """
    context=_context(project_id,source_triplet,target_icrs,target_radius_arcsec)
    path=Path(path).resolve(strict=True)
    if not path.is_file() or path.stat().st_size>MAX_ZIP:
        raise ValueError("Archive is not a regular file or exceeds review budget")
    keep={}
    with path.open("rb") as raw:
        digest=_digest(raw)
        if expected_sha256 is not None and digest!=expected_sha256:
            raise ValueError("Archive changed since the saved review receipt")
        raw.seek(0)
        with zipfile.ZipFile(raw) as z:
            infos=z.infolist()
            names=[i.filename for i in infos]
            if (len(names)!=len(set(names)) or not REQUIRED.issubset(names)
                    or set(names)-ALLOWED or any(i.is_dir() or i.flag_bits&1 for i in infos)
                    or sum(i.file_size for i in infos)>MAX_EXPANDED):
                raise ValueError("Unexpected archive members or expanded-size budget exceeded")
            for info in infos:
                limit=MAX_PNG if info.filename.endswith(".png") else MAX_TEXT
                if not info.filename.endswith(".npz") and info.file_size>limit:
                    raise ValueError("Archive member exceeds review budget")
            manifest=_json(z.read("manifest.json"))
            if (not isinstance(manifest,dict) or manifest.get("schema")!=SCHEMA
                    or manifest.get("status")!=STATUS
                    or manifest.get("covariance_calibrated") is not False
                    or manifest.get("full_geometry_containment_verified") is not False):
                raise ValueError("Unsupported or non-diagnostic manifest")
            hashes=manifest.get("files_sha256")
            if not isinstance(hashes,dict) or set(hashes)!=set(names)-{"manifest.json"}:
                raise ValueError("Manifest does not cover every archive member")
            for name,expected in hashes.items():
                if not _sha_ok(expected):
                    raise ValueError("Invalid member checksum")
                with z.open(name) as f:
                    if name in ("rings.csv","fit.json","asymptotic.json","publication_qa.json") or name.endswith(".png"):
                        data=f.read()
                        actual=hashlib.sha256(data).hexdigest()
                        keep[name]=_png(data) if name.endswith(".png") else data
                    else:
                        actual=_digest(f)
                if actual!=expected:
                    raise ValueError(f"Member checksum mismatch: {name}")
        raw.seek(0)
        if _digest(raw)!=digest:
            raise ValueError("Archive changed during review")
    bands=manifest.get("band_order")
    inputs=manifest.get("inputs")
    if (not isinstance(bands,list) or not bands or len(set(bands))!=len(bands)
            or any(b not in ("g","r","i","z") for b in bands)
            or not isinstance(inputs,list) or len(inputs)!=len(bands)):
        raise ValueError("Invalid diagnostic band order")
    anchors=[]
    for b,entry in zip(bands,inputs):
        if (not isinstance(entry,dict)
                or entry.get("dataset_id")!="/".join(context["source_triplet"]+[b])
                or not _sha_ok(entry.get("native_input_sha256"))
                or not isinstance(entry.get("grid_id"),str)
                or not entry["grid_id"].startswith("wcs-sha256:")
                or not _sha_ok(entry["grid_id"][len("wcs-sha256:"):])
                or "synthetic_injection_parent_sha256" in entry):
            raise ValueError("Archive does not match the selected native source")
        anchor=entry.get("anchor_icrs")
        if (not isinstance(anchor,list) or len(anchor)!=2
                or any(not isinstance(x,(float,int)) or isinstance(x,bool) or not math.isfinite(x) for x in anchor)
                or not 0<=anchor[0]<360 or not -90<=anchor[1]<=90):
            raise ValueError("Archive lacks a valid target anchor")
        anchors.append(anchor)
    if any(a!=anchors[0] for a in anchors):
        raise ValueError("Archive bands do not share an anchor")
    separation=_separation(context["target_icrs"],anchors[0])
    if separation>context["target_radius_arcsec"]:
        raise ValueError("Archive anchor is outside the current project target radius")
    fit=_json(keep["fit.json"])
    if (not isinstance(fit,dict) or fit.get("status")!=STATUS or fit.get("success") is not True
            or fit.get("boundary") is not False or fit.get("band_order")!=bands
            or fit.get("record_sha256")!=manifest.get("record_sha256")
            or fit.get("covariance_sha256")!=manifest.get("covariance_sha256")
            or not _sha_ok(fit.get("covariance_sha256"))
            or not isinstance(fit.get("record_sha256"),list)
            or len(fit["record_sha256"])!=len(bands)
            or not all(_sha_ok(v) for v in fit["record_sha256"])):
        raise ValueError("Fit provenance or status does not match manifest")
    if fit.get("covariance_calibrated") is not False or fit.get("science_adopted") is True:
        raise ValueError("Unexpected adopted or calibrated fit status")
    reader=csv.DictReader(io.StringIO(keep["rings.csv"].decode("utf-8")))
    cols=("band","a_inner_arcsec","a_outer_arcsec","observed_nMgy",
          "model_same_support_nMgy","conditional_sigma_nMgy","eligible","train","test")
    if reader.fieldnames is None or not set(cols).issubset(reader.fieldnames):
        raise ValueError("Incomplete ring table")
    rows=[]
    for row in reader:
        if len(rows)>=MAX_ROWS or row["band"] not in bands:
            raise ValueError("Unexpected ring table size or band")
        result={"band":row["band"]}
        for k in cols[1:6]:
            v=float(row[k])
            if not math.isfinite(v):
                raise ValueError("Nonfinite ring value")
            result[k]=v
        if (not 0<=result["a_inner_arcsec"]<result["a_outer_arcsec"]
                or result["conditional_sigma_nMgy"]<0):
            raise ValueError("Invalid ring bounds or uncertainty")
        for k in cols[6:]:
            if row[k] not in ("True","False"):
                raise ValueError("Invalid ring flags")
            result[k]=row[k]=="True"
        if (result["train"] and result["test"]) or ((result["train"] or result["test"]) and not result["eligible"]):
            raise ValueError("Inconsistent ring support")
        rows.append(result)
    if not rows:
        raise ValueError("Empty ring table")
    edges=[]
    for b in bands:
        subset=[r for r in rows if r["band"]==b]
        if not subset or subset[0]["a_inner_arcsec"]!=0 or any(
                l["a_outer_arcsec"]!=r["a_inner_arcsec"] for l,r in zip(subset,subset[1:])):
            raise ValueError("Noncontiguous ring table")
        edges.append([r["a_outer_arcsec"] for r in subset])
    if any(v!=edges[0] for v in edges):
        raise ValueError("Band ring radii differ")
    if [r["band"] for r in rows]!=[b for b in bands for _ in edges[0]]:
        raise ValueError("Ring table is not in declared band order")
    for flag in ("train","test"):
        expected=fit.get(flag+"_rings")
        if not isinstance(expected,list) or len(expected)!=len(edges[0]) or any(type(x) is not bool for x in expected):
            raise ValueError("Missing fit ring selection")
        if [r[flag] for r in rows]!=expected*len(bands):
            raise ValueError("Ring selection differs from fit")
    figures={name:keep[name+".png"] for name in FIGURES if name+".png" in keep}
    if manifest.get("six_figure_diagnostic") is True and (
            any(f"{name}.{ext}" not in names for name in FIGURES[:6] for ext in ("png","svg"))
            or not {"publication_arrays.npz","publication_qa.json"}.issubset(names)):
        raise ValueError("Incomplete declared six-figure package")
    asym=_json(keep["asymptotic.json"]) if "asymptotic.json" in keep else None
    qa=_json(keep["publication_qa.json"]) if "publication_qa.json" in keep else None
    if asym is not None and (not isinstance(asym,dict) or asym.get("status")!=STATUS
            or asym.get("record_sha256")!=fit["record_sha256"]
            or asym.get("covariance_sha256")!=fit["covariance_sha256"]
            or asym.get("asymptotic_flux_adopted") is not False
            or asym!=fit.get("asymptotic_diagnostic")):
        raise ValueError("Asymptotic report does not match the diagnostic fit")
    if qa is not None and (not isinstance(qa,dict) or qa.get("status")!=STATUS
            or qa.get("band_order")!=bands or qa.get("record_sha256")!=fit["record_sha256"]
            or qa.get("covariance_sha256")!=fit["covariance_sha256"]
            or qa.get("masked_pixels_reconstructed") is not False
            or qa.get("covariance_calibrated") is not False):
        raise ValueError("Publication QA does not match the diagnostic fit")
    return {"status":STATUS,"archive_path":str(path),"archive_sha256":digest,
            "context":context,"archive_anchor_icrs":anchors[0],
            "anchor_separation_arcsec":separation,"manifest":manifest,
            "fit":fit,"rows":rows,"figures":figures,
            "asymptotic":asym,"publication_qa":qa,
            "original_fits_reverified":False}


def save_review_receipt(path, review, *, confirmed):
    """Save an explicit association receipt, never overwrite or adopt science."""
    if confirmed is not True or review.get("status")!=STATUS:
        raise ValueError("Explicit archive-to-project association required")
    payload={"schema":RECEIPT_SCHEMA,"status":STATUS,"context":review["context"],
             "archive_path":review["archive_path"],"archive_sha256":review["archive_sha256"],
             "association_confirmed":True,"original_fits_reverified":False}
    data=(json.dumps(payload,indent=2,allow_nan=False)+"\n").encode("utf-8")
    with Path(path).open("xb") as f:
        f.write(data)


def reopen_review_receipt(path, **context):
    p=Path(path)
    if p.stat().st_size>MAX_RECEIPT:
        raise ValueError("Review receipt exceeds size budget")
    with p.open("rb") as f:
        receipt=_json(f.read(MAX_RECEIPT+1))
    if (not isinstance(receipt,dict) or receipt.get("schema")!=RECEIPT_SCHEMA
            or receipt.get("status")!=STATUS or receipt.get("association_confirmed") is not True
            or receipt.get("context")!=_context(**context)
            or not _sha_ok(receipt.get("archive_sha256"))
            or not isinstance(receipt.get("archive_path"),str)
            or not Path(receipt["archive_path"]).is_absolute()):
        raise ValueError("Receipt does not match the current project and source")
    return read_joint_archive(receipt["archive_path"],expected_sha256=receipt["archive_sha256"],**context)
