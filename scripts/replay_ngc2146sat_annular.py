"""Replay the archived NGC2146sat ring-domain reference through DESIDERA.

Requires local historical archives. No downloads, no FITS edits, no adoption.
Run with Python 3.10+ and an installed checkout. Output is a new JSON file.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import zipfile

import numpy as np

from lsbpipe.photometry.annular_inference import fit_conditional_gls, conditional_holdout
from lsbpipe.photometry.radial_models import sersic_intensity

TEMPLATE_MEMBER = "NGC2146sat_outer_control_template.npz"
TEMPLATE_HASH = "8ce18207e5c9adc9ffc4e31070c07e946708767d88c556feae4eb502d4caa355"
REFERENCE_HASH = "b677b20fa4d53beea80bf315fbeceeea1db0e00c36c2c7cc7cb125f684e38cb7"
HOLDOUT_HASH = "2c695197c2f679fb569383cc4211b6dcf31f81538e44a23c34eb2331798b775f"
KEY = "core26_compact_refit_50_L15_A50"
HOLDOUT_KEY = "core26_compact_refit_45_L15_A50"


def digest(data):
    return hashlib.sha256(data).hexdigest()


def checked_archive(path, expected):
    data=Path(path).read_bytes()
    if digest(data)!=expected:
        raise ValueError(f"Reference SHA-256 mismatch: {path}")
    return zipfile.ZipFile(io.BytesIO(data))


def csv_rows(z,name):
    return list(csv.DictReader(io.StringIO(z.read(name).decode("utf-8"))))


def replay(template_archive, reference_archive, holdout_archive):
    with zipfile.ZipFile(template_archive) as z:
        data=z.read(TEMPLATE_MEMBER)
    if digest(data)!=TEMPLATE_HASH:
        raise ValueError("Archived operator template SHA-256 mismatch")
    with np.load(io.BytesIO(data),allow_pickle=False) as z:
        t={name:z[name] for name in ("rg","ops_g","ops_r","d","C","init","lo","hi")}
    with checked_archive(reference_archive,REFERENCE_HASH) as z:
        rows=[r for r in csv_rows(z,"rings.csv") if r["key"]==KEY]
        fits=[r for r in csv_rows(z,"fits.csv") if r["key"]==KEY]
        with np.load(io.BytesIO(z.read("covariances.npz")),allow_pickle=False) as c:
            reference_cov=c[KEY]
    with checked_archive(holdout_archive,HOLDOUT_HASH) as z:
        held=[r for r in csv_rows(z,"holdout_45_50.csv") if r["key"]==HOLDOUT_KEY]
    ordered=[r for b in "gr" for r in sorted(
        (r for r in rows if r["band"]==b),key=lambda row:float(row["a"]))]
    if len(ordered)!=20 or len(fits)!=2 or len(held)!=1:
        raise ValueError("Incomplete archived reference")
    for i,b in enumerate("gr"):
        if [float(r["a"]) for r in ordered[i*10:(i+1)*10]]!=list(range(5,51,5)):
            raise ValueError("Reference ring order mismatch")
    d,c=t["d"],t["C"]
    np.testing.assert_allclose(d,[float(r["observed"]) for r in ordered],rtol=0,atol=1e-10)
    np.testing.assert_allclose(c,reference_cov,rtol=0,atol=1e-10)
    # Original archived vector uses log I0, log scale, linear index.
    # The new fitting convention uses log for all three parameters.
    initial,lo,hi=(t[k].copy() for k in ("init","lo","hi"))
    for a in (initial,lo,hi):
        a[2::3]=np.log(a[2::3])
    def predict(z):
        return np.concatenate([
            t["ops_"+b]@sersic_intensity(t["rg"],*np.exp(z[3*i:3*i+3]))
            for i,b in enumerate("gr")])
    train=np.tile(np.array([False]+[True]*9),2)
    full=fit_conditional_gls(d,c,predict,initial,(lo,hi),train=train)
    if not full["success"] or full["boundary"] or full["conditional_parameter_covariance"] is None:
        raise RuntimeError("Reference fit failed or is not an interior full-rank solution")
    np.testing.assert_allclose(full["chi2"],float(fits[0]["chi2"]),rtol=0,atol=1e-5)
    ref_model=np.array([float(r["model"]) for r in ordered])
    model_error=float(np.max(np.abs(full["prediction"]-ref_model)))
    if model_error>5e-4:
        raise AssertionError("Archived ring prediction mismatch")
    train=np.tile(np.array([False]+[True]*8+[False]),2)
    test=np.tile(np.array([False]*9+[True]),2)
    fitted=fit_conditional_gls(d,c,predict,initial,(lo,hi),train=train)
    if not fitted["success"] or fitted["boundary"] or fitted["conditional_parameter_covariance"] is None:
        raise RuntimeError("Held-out reference fit failed")
    p=fitted["parameters"]
    j=np.column_stack([(predict(p+v*1e-5)-predict(p-v*1e-5))/(2e-5)
                       for v in np.eye(len(p))])
    result=conditional_holdout(d,fitted["prediction"],c,j,
                              fitted["conditional_parameter_covariance"],
                              train=train,test=test)
    np.testing.assert_allclose(result["chi2"],float(held[0]["chi2"]),rtol=0,atol=5e-4)
    expected=np.array([float(held[0]["observed_"+b])-float(held[0]["prediction_"+b])
                       for b in "gr"])
    np.testing.assert_allclose(result["residual"],expected,rtol=0,atol=5e-4)
    return {"status":"DIAGNOSTIC_NOT_ADOPTED",
            "scope":"ARCHIVED_RING_DOMAIN_ONLY_NOT_NATIVE_FITS_REPLAY",
            "reference_key":KEY,"band_order":["g","r"],
            "source_sha256":{"template_member":TEMPLATE_HASH,"reference":REFERENCE_HASH,
                             "holdout":HOLDOUT_HASH},
            "full_fit":{"chi2":full["chi2"],"dof":full["dof"],
                        "physical_parameter_vector":np.exp(full["parameters"]).tolist(),
                        "parameter_order_per_band":["I0_nMgy_per_pixel","scale_arcsec","index"],
                        "max_model_difference_nMgy":model_error},
            "holdout":{"chi2":result["chi2"],"residual_nMgy":result["residual"].tolist(),
                       "covariance_nMgy2":result["predictive_covariance"].tolist()},
            "covariance_calibrated":False,"native_fits_replayed":False}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--template-archive",required=True,type=Path)
    p.add_argument("--reference-archive",required=True,type=Path)
    p.add_argument("--holdout-archive",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    result=replay(a.template_archive,a.reference_archive,a.holdout_archive)
    with a.output.open("x",encoding="utf-8") as f:
        json.dump(result,f,indent=2,allow_nan=False)
        f.write("\n")
    print(json.dumps(result,indent=2,allow_nan=False))


if __name__=="__main__":
    main()
