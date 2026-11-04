"""Non-adopted, held-out-star residual QA for *signed* ePSF candidates.

A model that passes a global pixel-sign test is not necessarily a
better physical PSF. Here a signed, unmodified ePSF is amplitude-fitted
ONLY in each independent held-out star's core and evaluated in three
radial zones, including the outer cutout used for background modelling.

These are diagnostic residuals divided by observed stellar peak, NOT
noise-whitened likelihoods, PSF-wing validation or a calibrated error
budget. Original pixels and signed negative model components are intact.
"""
from __future__ import annotations

import math
import numpy as np
from scipy.ndimage import shift as shift_image


def heldout_signed_residual(
    *,
    source_id: str,
    image: np.ndarray,
    sqrt_ivar: np.ndarray,
    signed_kernel: np.ndarray,
    local_center_xy: tuple[float, float],
    guessed_fwhm_pix: float,
    measured_peak: float,
) -> dict:
    """Fit core-only amplitude and inspect signed model residual by annulus."""
    science = np.asarray(image, dtype=float)
    weights = np.asarray(sqrt_ivar, dtype=float)
    kernel = np.asarray(signed_kernel, dtype=float)
    if (science.ndim != 2 or science.shape != kernel.shape
            or science.shape != weights.shape or min(science.shape) < 15
            or not np.all(np.isfinite(science))
            or not np.all(np.isfinite(weights))
            or np.any(weights < 0)
            or not np.all(np.isfinite(kernel))):
        raise ValueError("Signed holdout inputs must share a finite native-grid shape")
    if (len(local_center_xy) != 2
            or not all(math.isfinite(float(v)) for v in local_center_xy)
            or not math.isfinite(guessed_fwhm_pix)
            or guessed_fwhm_pix <= 0
            or not math.isfinite(measured_peak)
            or measured_peak <= 0):
        raise ValueError("Signed holdout requires valid center, FWHM guess and peak")

    cx0, cy0 = (np.array(kernel.shape[::-1])-1)/2
    cx, cy = map(float, local_center_xy)
    model = shift_image(
        kernel, (cy-cy0, cx-cx0), order=3, mode="constant",
        cval=0., prefilter=True)
    rr_y, rr_x = np.indices(science.shape, dtype=float)
    radius = np.hypot(rr_x-cx, rr_y-cy)
    good = weights > 0
    core = good & (radius <= float(guessed_fwhm_pix))
    if np.count_nonzero(core) < 5:
        raise ValueError("Too few independent holdout core pixels")
    w2 = weights[core] ** 2
    denominator = float(np.sum(w2 * model[core]**2))
    if not math.isfinite(denominator) or denominator <= 0:
        raise ValueError("Cannot amplitude-fit signed model in holdout core")
    amplitude = float(np.sum(w2 * model[core] * science[core]) / denominator)
    if not math.isfinite(amplitude) or amplitude <= 0:
        raise ValueError("Signed holdout core amplitude is not positive")
    residual = science - amplitude*model
    annuli = {
        "core_r_le_1x_guessed_fwhm": core,
        "intermediate_1x_to_2x_guessed_fwhm": (
            good & (radius > guessed_fwhm_pix) &
            (radius <= 2*guessed_fwhm_pix)),
        "outer_r_gt_2x_guessed_fwhm": (
            good & (radius > 2*guessed_fwhm_pix)),
    }
    metrics = {}
    for name, region in annuli.items():
        values = residual[region]/measured_peak
        metrics[name] = {
            "n_valid_pixels": int(values.size),
            "median_signed_residual_over_peak": (
                float(np.median(values)) if values.size else None),
            "median_absolute_residual_over_peak": (
                float(np.median(np.abs(values))) if values.size else None),
        }
    return {
        "source_id": str(source_id),
        "fitted_amplitude_model_units": amplitude,
        "radial_residuals": metrics,
        "fwhm_guess_not_measured": float(guessed_fwhm_pix),
        "signed_kernel_unchanged": True,
        "correlated_coadd_noise_not_whitened": True,
    }


def summarize_signed_holdouts(reports: list[dict]) -> dict:
    """Report independent-star summaries without inventing precision."""
    names = (
        "core_r_le_1x_guessed_fwhm",
        "intermediate_1x_to_2x_guessed_fwhm",
        "outer_r_gt_2x_guessed_fwhm",
    )
    # Repeated reports from one physical star are not independent trials.
    # This is an integrity guard, not a certification of distinct source IDs
    # as astrophysically independent stars or survey exposures.
    seen_ids: set[str] = set()
    for row in reports:
        if not isinstance(row, dict):
            raise TypeError("Each signed holdout report must be a mapping")
        source_id = row.get("source_id")
        if not isinstance(source_id, str) or not source_id.strip():
            raise ValueError("Each holdout report requires a nonempty source_id")
        source_id = source_id.strip()
        if source_id in seen_ids:
            raise ValueError(
                f"Duplicate holdout source_id cannot count as an independent star: {source_id}")
        seen_ids.add(source_id)
    result = {"n_usable_holdout_stars": len(seen_ids), "regions": {}}
    for name in names:
        signed = []
        absolute = []
        for row in reports:
            zone = row["radial_residuals"][name]
            if zone["n_valid_pixels"] < 4:
                continue
            if zone["median_signed_residual_over_peak"] is None:
                continue
            signed.append(zone["median_signed_residual_over_peak"])
            absolute.append(zone["median_absolute_residual_over_peak"])
        result["regions"][name] = {
            "n_stars": len(signed),
            "median_star_signed_residual_over_peak": (
                float(np.median(signed)) if signed else None),
            "median_star_abs_residual_over_peak": (
                float(np.median(absolute)) if absolute else None),
        }
    return result
