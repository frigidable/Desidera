"""Serializable PSF-0 debugging evidence: complete metadata, no FITS/pixel grids."""
from __future__ import annotations

import numpy as np


def _per_star(star, fwhm):
    data = np.asarray(star["cutout"], dtype=float)
    model = np.asarray(star["model"], dtype=float)
    valid = np.asarray(star["valid"], dtype=bool)
    cx, cy = star["center_in_cutout"]
    yy, xx = np.indices(data.shape)
    radius = np.hypot(xx - float(cx), yy - float(cy))
    residual = data - model
    peak = max(float(star["peak_nMgy"]), 1e-10)

    def region(mask):
        pixels = residual[mask & valid]
        if len(pixels) == 0:
            return None
        return {
            "n_valid": int(len(pixels)),
            "median_abs_over_star_peak": float(np.median(np.abs(pixels))/peak),
            "median_signed_over_star_peak": float(np.median(pixels)/peak),
        }
    return {
        "star_id": star["star_id"], "role": star["role"],
        "gaia": star.get("gaia", {}),
        "shape_qa": star.get("shape_qa", {}),
        "native_pixel_xy": [star["x_pix"], star["y_pix"]],
        "icrs_deg": [star["ra_deg"], star["dec_deg"]],
        "subpixel_center_xy": list(star["center_in_cutout"]),
        "peak_nMgy_per_pix": star["peak_nMgy"],
        "sky_nMgy_per_pix": star["sky_nMgy"],
        "gradient_sigma_span": star["gradient_sigma_span"],
        "sky_qa": star.get("sky_qa", {}),
        "n_invalid_stamp_pixels": int(np.count_nonzero(~valid)),
        "core": region(radius <= fwhm),
        "intermediate": region((radius > fwhm) & (radius <= 2*fwhm)),
        "outer": region(radius > 2*fwhm),
        "statistics_are_not_covariance_calibrated": True,
    }


def detailed_band_report(panel):
    """No raw image arrays, no network, no scientific adoption."""
    review = panel._star_review_data or {}
    fwhm = float(review.get("guessed_fwhm_pix", panel.fwhm.value()))
    return {
        "schema": "desidera.psf0-band-debug-report/1",
        "science_status": "DIAGNOSTIC_NOT_SCIENCE_ADOPTED",
        "band": panel.band.currentData(),
        "project_root": str(panel._project_root) if panel._project_root else None,
        "experiment": panel.experiment.currentData(),
        "star_catalog_policy": panel.catalog_policy.currentData(),
        "detection_fwhm_guess_pix_not_measured": panel.fwhm.value(),
        "fwhm_strategy": ("manual" if panel._fwhm_manually_edited
                          else "local_psfsize"),
        "last_verified_band_seeing": panel._latest_seeing_payload,
        "exact_result_recipe": (
            panel._last_reported_payload.get("effective_recipe")
            if panel._last_reported_payload else None),
        "exact_result_recipe_sha256": (
            panel._last_reported_payload.get("recipe_sha256")
            if panel._last_reported_payload else None),
        "original_source_sha256": (
            panel._last_reported_payload.get("source_hashes")
            if panel._last_reported_payload else None),
        "source_project_uuid": (
            panel._last_reported_payload.get("project_id")
            if panel._last_reported_payload else None),
        "source_dataset": (
            panel._last_reported_payload.get("dataset")
            if panel._last_reported_payload else None),
        "target_protection_arcsec": panel.protection.value(),
        "seeing_provenance": panel.prior_source.text(),
        "gaia_catalogue_provenance": panel.gaia_source.text(),
        "status": panel.status.text(),
        "summary": panel.summary.text(),
        "full_visible_qa": panel.qa.toPlainText(),
        "input_quartet_paths_presence_unverified": panel._snapshot,
        "science_report_or_rejection": panel._qa_payload,
        "historical_trials_not_adopted": list(panel._experiment_history),
        "manual_excluded_training_ids": list(panel._active_exclusions),
        "stars_in_training_review": [
            _per_star(star, fwhm) for star in review.get("stars", [])],
        "limitations": [
            "Core/wing residuals not covariance-calibrated chi-square",
            "Coadd noise correlations unmodeled",
            "Sky annulus can absorb wings; no joint star/sky fit",
            "Global PSF lacks spatial and extended-wing validation",
            "No raw stamps or input FITS bytes are included",
        ],
    }
