"""Preliminary *star-level* ePSF residual screening, distinct from builder pixel clipping.

Unlike EPSFBuilder's pixelwise 3-sigma rejection, this evaluates entire
native star cutouts after amplitude fitting and flags gross profile
outliers for an explicitly requested *experimental* rebuild.

The dimensionless error scores are NOT chi-square: DR11 coadds have
correlated noise and an uncertain residual sky. Selection must not use
holdout stars or adjust scientific adoption gates to obtain a pass.
"""
from __future__ import annotations

import math
import numpy as np
from scipy.ndimage import shift


def training_star_residual(entry, signed_kernel, *, guessed_fwhm_pix):
    """Evaluate one training cutout with an amplitude fitted to its core."""
    ident = entry[0]
    c = entry[-1]
    data = np.asarray(c["image"], dtype=float)
    weights = np.asarray(c["weights"], dtype=float)
    kernel = np.asarray(signed_kernel, dtype=float)
    if data.shape != kernel.shape or weights.shape != data.shape:
        raise ValueError("Training star and ePSF grids disagree")
    x0, y0 = c["center_in_cutout"]
    cy, cx = (np.array(kernel.shape)-1)/2.
    model = shift(kernel, (y0-cy, x0-cx), order=3,
                  mode="constant", cval=0., prefilter=True)
    yy, xx = np.indices(data.shape)
    radial = np.hypot(xx-x0, yy-y0)
    valid = np.isfinite(data) & np.isfinite(weights) & (weights > 0)
    inner = valid & (radial <= guessed_fwhm_pix)
    if np.count_nonzero(inner) < 12:
        raise ValueError("Too few pixels for core amplitude")
    invvar = weights[inner]**2
    denominator = float(np.sum(invvar*model[inner]**2))
    numerator = float(np.sum(invvar*model[inner]*data[inner]))
    if denominator <= 0 or numerator <= 0:
        raise ValueError("Invalid PSF core amplitude")
    amplitude = numerator/denominator
    predicted = amplitude*model
    # This is a deliberately conservative score scale, not an
    # inverse covariance: per-pixel IVAR cannot describe coadd noise.
    noise = np.maximum(
        np.where(valid, 1. / np.maximum(weights,1e-12), np.inf),
        max(float(c["rms_nMgy"]), 1e-12))
    scale = np.hypot(noise, 0.03*np.abs(predicted))
    residual = np.abs(data-predicted) / scale
    core = valid & (radial <= 1.5*guessed_fwhm_pix)
    outer = valid & (radial > 1.5*guessed_fwhm_pix) & (
        radial <= min(3.5*guessed_fwhm_pix, min(data.shape)/2.-1))
    if np.count_nonzero(core) < 15 or np.count_nonzero(outer) < 25:
        raise ValueError("Insufficient masked-free star pixels for outlier QA")
    core_p75 = float(np.quantile(residual[core], .75))
    outer_p90 = float(np.quantile(residual[outer], .90))
    tail_p99 = float(np.quantile(residual[valid], .99))
    score = max(core_p75, outer_p90, tail_p99/2.)
    if not np.isfinite(score):
        raise ValueError("Nonfinite star residual screening statistic")
    return {
        "star_id":ident,
        "score":score,
        "core_p75":core_p75,
        "outer_p90":outer_p90,
        "all_p99_div_2":tail_p99/2.,
        "amplitude_nMgy_diagnostic":float(amplitude),
        "metric":"PERCENTILE_WHITENED_RESIDUAL_HEURISTIC_NOT_CHI_SQUARE",
    }


def choose_training_outliers(rows, *, min_train, sigma_threshold=4.,
                             max_remove_fraction=.12):
    """Robust, capped trimming of one training-only scoring distribution.

    Median and MAD are estimated anew for each training iteration.
    A 4*MAD outlier must also exceed 1.8 times the median to be
    excluded. A small training sample never gets auto-trimmed.
    """
    if len(rows) < max(12, min_train+3):
        return (), {"reason":"SMALL_SAMPLE_NOT_AUTOTRIMMED",
                    "n_scored":len(rows)}
    scores = np.asarray([r["score"] for r in rows], dtype=float)
    if not np.isfinite(scores).all() or np.any(scores < 0):
        raise ValueError("Invalid training residual scores")
    median = float(np.median(scores))
    mad = float(1.4826*np.median(np.abs(scores-median)))
    cutoff = median + max(sigma_threshold*mad, 1.)
    cutoff = max(cutoff, 1.8*median)
    candidates = sorted(
        (r for r in rows if r["score"] > cutoff),
        key=lambda r:(-r["score"],r["star_id"]))
    max_remove = min(
        len(rows)-min_train,
        max(1, int(math.floor(max_remove_fraction*len(rows)))))
    chosen = tuple(sorted(r["star_id"] for r in candidates[:max_remove]))
    return chosen, {
        "n_scored":len(rows),"median_score":median,
        "robust_mad_sigma":mad,"cutoff":cutoff,
        "n_above_cutoff":len(candidates),
        "n_removed_this_pass":len(chosen),
        "max_removed_this_pass":max_remove,
        "metric":"TRAIN_ONLY_PROVISIONAL_NONCALIBRATED",
    }
