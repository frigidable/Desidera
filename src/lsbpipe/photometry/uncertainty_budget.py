"""Conditional calibration and tested-scenario summaries; never science adoption."""
from __future__ import annotations

import numpy as np

STATUS = "DIAGNOSTIC_NOT_ADOPTED"


def _covariance(value, size):
    c = np.asarray(value, dtype=float)
    if c.shape != (size, size) or not np.isfinite(c).all():
        raise ValueError("Expected finite square covariance in parameter order")
    scale = max(float(np.max(np.abs(c))), np.finfo(float).tiny)
    if not np.allclose(c, c.T, rtol=0, atol=1e-12 * scale):
        raise ValueError("Covariance must be symmetric")
    c = (c + c.T) / 2
    if np.linalg.eigvalsh(c).min() < -1e-12 * scale:
        raise ValueError("Covariance must be positive semidefinite")
    return c


def conditional_magnitude_budget(stat_cov, extinction_coefficients, sigma_ebv,
                                 zeropoint_cov):
    """AB magnitude covariance, assuming independent statistics, EBV and ZP.

    ZP covariance is caller-supplied; EBV is common to all bands. This does
    not calibrate either input or include model-family/systematic shifts.
    """
    k = np.asarray(extinction_coefficients, dtype=float)
    s = float(sigma_ebv)
    if k.ndim != 1 or not k.size or not np.isfinite(k).all() or not np.isfinite(s) or s < 0:
        raise ValueError("Invalid extinction inputs")
    stat = _covariance(stat_cov, k.size)
    zp = _covariance(zeropoint_cov, k.size)
    cal = zp + s**2 * np.outer(k, k)
    combined = _covariance(stat + cal, k.size)
    return {"status": STATUS, "stat_cov": stat, "assumed_cal_cov": cal,
            "conditional_plus_assumed_cal_cov": combined}


def linear_sigma(covariance, weights):
    """Standard deviation of a linear combination in the given band order."""
    w = np.asarray(weights, dtype=float)
    if w.ndim != 1 or not w.size or not np.isfinite(w).all():
        raise ValueError("Expected finite weights")
    c = _covariance(covariance, w.size)
    return float(np.sqrt(max(0., w @ c @ w)))


def scenario_envelope(values, reference, eligible):
    """Finite tested extrema, NOT a confidence interval or a probability law."""
    x = np.asarray(values, dtype=float)
    ok = np.asarray(eligible)
    ref = float(reference)
    if x.ndim != 1 or ok.shape != x.shape or ok.dtype.kind != 'b' or not np.isfinite(ref):
        raise ValueError("Expected a reference and explicit Boolean eligibility")
    x = x[ok]
    if not x.size or not np.isfinite(x).all():
        raise ValueError("No finite eligible scenario set")
    low, high = float(x.min()), float(x.max())
    return {"status": STATUS, "kind": "TESTED_EXTREMA_NOT_CONFIDENCE_INTERVAL",
            "n": int(x.size), "reference": ref, "minimum": low, "maximum": high,
            "delta_min": low - ref, "delta_max": high - ref}
