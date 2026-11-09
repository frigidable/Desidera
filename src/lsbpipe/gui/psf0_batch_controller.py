"""Sequential griz PSF-0 diagnostics with immutable project checkpoints."""
from dataclasses import replace
import traceback

from lsbpipe.photometry.psf_gui_backend import (
    psf0_experiment_recipe, run_gui_psf0, suggest_gui_psf0_fwhm,
    checkpoint_gui_psf0,
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
                if (req.get("excluded_training_star_ids") and
                        abs(fwhm - float(req["fwhm_pix"])) > .015):
                    raise ValueError(
                        "Stored star-XXXX exclusions were reviewed with "
                        f"FWHM {float(req['fwhm_pix']):.2f} px but current "
                        f"{band}-band PSFSIZE gives {fwhm:.2f} px. The "
                        "source detection inventory can change. Apply "
                        "the new seeing prior, rebuild the PSF, then "
                        "review train stars again. No stale exclusions "
                        "were applied.")
                cfg = replace(
                    psf0_experiment_recipe(fwhm, req["experiment"]),
                    star_catalog_policy=req["star_catalog_policy"],
                    excluded_training_star_ids=tuple(
                        req.get("excluded_training_star_ids", ())))
                extended_qa = req.get("extended_research_qa", False)
                if type(extended_qa) is not bool:
                    raise ValueError("Full research QA option must be boolean")
                payload = run_gui_psf0(
                    root, band=band, config=cfg,
                    protection_arcsec=float(req["target_protection_arcsec"]),
                    extended_research_qa=extended_qa)
                if payload.get("status") == "EXPERIMENTAL_PSF_NOT_ADOPTED":
                    try:
                        payload["auto_saved_path"] = checkpoint_gui_psf0(payload)
                    except (OSError, ValueError, KeyError, TypeError) as exc:
                        payload["checkpoint_error"] = (
                            type(exc).__name__ + ": " + str(exc))
                outcomes.append({"band":band, "payload":payload, "seeing":seeing})
            except Exception as exc:
                outcomes.append({
                    "band":band, "error":type(exc).__name__+": "+str(exc),
                    "traceback":traceback.format_exc(limit=8)})
            progress((index+1,len(requests),
                      f"PSF-0 {band} finished {index+1}/{len(requests)}"))
        return outcomes

    def finished(outcomes):
        notes = []
        for entry in outcomes:
            band = entry["band"]
            panel = workspace.psf0_panels[band]
            if "error" in entry:
                panel.show_failure(entry["error"]+"\n\n"+entry["traceback"])
                notes.append(band+": ERROR")
                continue
            if entry["seeing"] is not None:
                panel.show_fwhm_suggestion(entry["seeing"])
            payload = entry["payload"]
            if payload.get("status") == "PSF0_QUALITY_REJECTED_NOT_ADOPTED":
                panel.show_quality_rejection(payload)
                notes.append(band+": QUALITY REJECTED")
            else:
                panel.show_result(payload)
                if payload.get("checkpoint_error"):
                    panel.status.setText(
                        "PSF-0 checkpoint FAILED: " + payload["checkpoint_error"]
                        + ". The preview remains in memory; save before exit.")
                    notes.append(band+": CHECKPOINT FAILED")
                else:
                    notes.append(band+": SAVED EXPERIMENTAL QA")
        for skipped in batch.get("skipped_bands", ()):
            notes.append(str(skipped["band"])+": SKIPPED ("+
                         str(skipped["reason"])+")")
        workspace.batch_status.setText(
            "griz batch completed · "+" | ".join(notes)
            + ". All PSFs remain NOT ADOPTED. Original FITS unchanged. "
              "Successful checkpoints are immutable derivatives.")
        window.message.setText(
            "Sequential griz PSF-0 diagnostics finished · no PSF adopted.")

    def failed(message):
        workspace.batch_status.setText(
            "PSF batch worker failed: "+message+". No PSF adopted.")

    window._start(calculate, finished, failed, operation_name="psf0-griz-batch")
    if window.thread is not None:
        workspace.set_psf_busy(True)
        workspace.batch_status.setText(
            f"Sequential PSF-0 running for {len(requests)} native bands, "
            "original FITS read-only.")
