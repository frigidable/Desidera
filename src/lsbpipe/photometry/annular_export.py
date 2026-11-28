"""Portable, non-adopted annular diagnostics: signed CSV, JSON, NPZ and plots."""
from __future__ import annotations

import csv
import hashlib
import io
import json
from pathlib import Path
import zipfile

import numpy as np

from .native_annular_diagnostic import annular_record_sha256, numeric_sha256
from .uncertainty_budget import _covariance, STATUS


def _jsonable(obj):
    if isinstance(obj,np.ndarray):
        return _jsonable(obj.tolist())
    if isinstance(obj,np.generic):
        return _jsonable(obj.item())
    if isinstance(obj,dict):
        return {k:_jsonable(v) for k,v in obj.items()}
    if isinstance(obj,(list,tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj,float) and not np.isfinite(obj):
        return None
    return obj


def export_annular_diagnostic(path, records, covariance, fit, *, render=False):
    """Write a NEW archive for an exact joint fit and its measured support.

    No calibrated confidence intervals, reconstructed masked flux or physical
    outer component claims. Plots show signed valid-pixel flux; local pixel
    area is only an approximate projected area.
    """
    records=tuple(records)
    if not records or fit.get("status")!=STATUS:
        raise ValueError("Expected a non-adopted diagnostic fit")
    if tuple(annular_record_sha256(r) for r in records)!=tuple(fit["record_sha256"]):
        raise ValueError("Fit and annular records do not match")
    if numeric_sha256(covariance)!=fit["covariance_sha256"]:
        raise ValueError("Fit and covariance do not match")
    bands=tuple(fit["band_order"])
    if len(bands)!=len(records):
        raise ValueError("Band order mismatch")
    n=len(records[0]["ring_edges_arcsec"])-1
    c=_covariance(covariance,len(bands)*n)
    prediction=np.asarray(fit["prediction"],dtype=float)
    if prediction.shape!=(len(c),) or not np.isfinite(prediction).all():
        raise ValueError("Invalid fit prediction")
    rows=[]
    arrays={"covariance_nMgy2":c,"prediction_nMgy":prediction,
            "train_rings":fit["train_rings"],"test_rings":fit["test_rings"]}
    for i,(b,r) in enumerate(zip(bands,records)):
        sl=slice(i*n,(i+1)*n)
        observed=r["observed_flux_nMgy"]
        counts=r["valid_pixel_count"]
        area=counts*r["pixel_area_arcsec2"]
        errors=np.sqrt(np.maximum(0,np.diag(c)[sl]))
        for j in range(n):
            rows.append({"band":b,"a_inner_arcsec":r["ring_edges_arcsec"][j],
                         "a_outer_arcsec":r["ring_edges_arcsec"][j+1],
                         "observed_nMgy":observed[j],"model_same_support_nMgy":prediction[sl][j],
                         "conditional_sigma_nMgy":errors[j],
                         "valid_pixels":counts[j],"local_valid_area_arcsec2":area[j],
                         "eligible":bool(r["coverage"]["eligible"][j]),
                         "ring_fraction_in_frame":r["coverage"]["ring_fraction"][j],
                         "worst_sector_fraction_in_frame":r["coverage"]["sector_fraction"][j].min(),
                         "train":bool(fit["train_rings"][j]),"test":bool(fit["test_rings"][j])})
        for name in ("operator","radial_nodes_arcsec","ring_edges_arcsec",
                     "observed_flux_nMgy","valid_pixel_count"):
            arrays[f"band{i}_{name}"]=r[name]
        if "full_aperture_operator" in r:
            arrays[f"band{i}_full_aperture_operator"]=r["full_aperture_operator"]
    csvbuf=io.StringIO()
    writer=csv.DictWriter(csvbuf,fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
    npbuf=io.BytesIO()
    np.savez_compressed(npbuf,**arrays)
    contents={"rings.csv":csvbuf.getvalue().encode("utf-8"),"arrays.npz":npbuf.getvalue()}
    fit_meta={k:v for k,v in fit.items() if k!="prediction"}
    contents["fit.json"]=json.dumps(_jsonable(fit_meta),indent=2,allow_nan=False).encode("utf-8")
    if "asymptotic_diagnostic" in fit:
        asym=fit["asymptotic_diagnostic"]
        if (asym.get("status")!=STATUS or tuple(asym["record_sha256"])!=tuple(fit["record_sha256"])
                or asym["covariance_sha256"]!=fit["covariance_sha256"]
                or asym.get("asymptotic_flux_adopted") is not False):
            raise ValueError("Asymptotic diagnostic provenance mismatch")
        contents["asymptotic.json"]=json.dumps(_jsonable(asym),indent=2,allow_nan=False).encode("utf-8")
    if render:
        # Optional matplotlib dependency, no GUI/backend or external inputs.
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        for kind in ("curve_of_growth","radial_intensity"):
            fig=Figure(figsize=(8,5),layout="constrained")
            FigureCanvasAgg(fig)
            ax=fig.add_subplot()
            for i,(b,r) in enumerate(zip(bands,records)):
                sl=slice(i*n,(i+1)*n)
                observed=r["observed_flux_nMgy"]
                model=prediction[sl]
                a=r["ring_edges_arcsec"][1:]
                block=c[sl,sl]
                if kind=="curve_of_growth":
                    l=np.tril(np.ones((n,n)))
                    observed=l@observed
                    model=l@model
                    sigma=np.sqrt(np.maximum(0,np.diag(l@block@l.T)))
                else:
                    area=r["valid_pixel_count"]*r["pixel_area_arcsec2"]
                    if np.any(area<=0):
                        raise ValueError("Missing measured area")
                    observed=observed/area
                    model=model/area
                    sigma=np.sqrt(np.maximum(0,np.diag(block)))/area
                good=r["coverage"]["eligible"]
                line=ax.errorbar(a[good],observed[good],yerr=sigma[good],
                                 fmt="o",label=f"{b}: observed eligible")
                color=line[0].get_color()
                ax.plot(a,model,color=color,linestyle="--",label=f"{b}: PSF model on same support")
                if np.any(~good):
                    ax.plot(a[~good],observed[~good],marker="x",linestyle="none",color=color,
                            label=f"{b}: outside coverage gate")
            ax.set_xlabel("Semimajor radius (arcsec)")
            ax.set_ylabel("Signed observed cumulative flux (nMgy)" if kind=="curve_of_growth"
                          else "Signed mean intensity (nMgy / arcsec²)")
            ax.set_title("Conditional diagnostic — not adopted")
            ax.legend(fontsize=8)
            for ext in ("png","svg"):
                buf=io.BytesIO()
                fig.savefig(buf,format=ext,dpi=160)
                contents[f"{kind}.{ext}"]=buf.getvalue()
            fig.clear()
    manifest={"status":STATUS,"schema":"desidera.annular-diagnostic/1",
              "band_order":bands,"record_sha256":fit["record_sha256"],
              "covariance_sha256":fit["covariance_sha256"],
              "covariance_source":fit["covariance_source"],
              "inputs":[{k:r[k] for k in ("dataset_id","grid_id","native_input_sha256",
                    "scene_revision","valid_sha256","psf_sha256","corrections",
                    "anchor_icrs","geometry","coverage_thresholds")} for r in records],
              "covariance_calibrated":False,"full_geometry_containment_verified":False,
              "model_total_kind":"UNTRUNCATED_MODEL_EXTRAPOLATION",
              "files_sha256":{name:hashlib.sha256(data).hexdigest() for name,data in contents.items()}}
    contents["manifest.json"]=json.dumps(_jsonable(manifest),indent=2,allow_nan=False).encode("utf-8")
    path=Path(path)
    with zipfile.ZipFile(path,"x",compression=zipfile.ZIP_DEFLATED) as z:
        for name,data in contents.items():
            z.writestr(name,data)
    with zipfile.ZipFile(path) as z:
        if z.testzip() is not None:
            raise IOError("Diagnostic archive failed CRC verification")
    return manifest
