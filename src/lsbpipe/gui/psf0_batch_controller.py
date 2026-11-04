"""Sequential griz PSF-0 previews on the existing worker; no file writes."""
from dataclasses import replace
import traceback

from lsbpipe.photometry.psf_gui_backend import (
    psf0_experiment_recipe, run_gui_psf0, suggest_gui_psf0_fwhm,
)


def launch(window, batch):
    workspace = window.photometry_workspace
    root = batch.get("project_root")
    requests = list(batch.get("requests", ()))
    bands = [r.get("band") for r in requests]
    if window.thread is not None:
        workspace.batch_status.setText("Another task active; batch not started.")
        return
    if (window._object_project is None or
            root != str(window._object_project.root) or
            not 1 <= len(requests) <= 4 or len(set(bands)) != len(bands) or
            any(b not in ("g", "r", "i", "z") for b in bands) or
            any(r.get("project_root") != root for r in requests)):
        workspace.batch_status.setText("Invalid project/bands; nothing started.")
        return

    def calculate(progress):
        outcomes = []
        for index, req in enumerate(requests):
            band = req["band"]
            progress((index, len(requests), f"PSF-0 {band}: checking PSFSIZE/FITS"))
            try:
                if req.get("manual_review_confirmed"):
                    raise ValueError("Manual review requires separate confirmation")
                seeing = None
                strategy = req.get("fwhm_strategy", "manual")
                if strategy == "local_psfsize":
                    seeing = suggest_gui_psf0_fwhm(root, band=band)
                    info = seeing["seeing"]
                    if not info.get("available"):
                        raise ValueError(
                            "No valid PSFSIZE for "+band+": "+str(info.get("message")))
                    fwhm = float(info["suggested_detection_fwhm_pix"])
                elif strategy == "manual":
                    fwhm = float(req["fwhm_pix"])
                else:
                    raise ValueError("Unsupported FWHM source")
                cfg = replace(
                    psf0_experiment_recipe(fwhm, req["experiment"]),
                    star_catalog_policy=req["star_catalog_policy"],
                    excluded_training_star_ids=tuple(
                        req.get("excluded_training_star_ids", ())))
                payload = run_gui_psf0(
                    root, band=band, config=cfg,
                    protection_arcsec=float(req["target_protection_arcsec"]))
                outcomes.append({"band":band, "payload":payload, "seeing":seeing})
            except Exception as exc:
                outcomes.append({
                    "band":band, "error":type(exc).__name__+": "+str(exc),
                    "traceback":traceback.format_exc(limit=8)})
            progress((index+1,len(requests),
                      f"PSF-0 {band} finished {index+1}/{len(requests)}"))
        return outcomes

