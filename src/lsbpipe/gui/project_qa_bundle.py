"""One-click, project-bound, read-only diagnostic bundle for Photometry.

Collects the currently available PSF-0, Scene, residual-sky and finite
aperture evidence without running analyses, opening FITS, changing a draft,
or silently presenting a stale stage result as current. Large repeated
isophote radial traces are encoded losslessly for core numerical columns;
the full original stage reports remain available via their own Copy QA.
"""
from __future__ import annotations

from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
import json
import platform
import sys

import numpy as np

_SCHEMA = "desidera.gui.project-diagnostic-bundle/1"


def _package_version(name):
    try:
        return version(name)
    except PackageNotFoundError:
        return None


def _safe_stage(name, action):
    try:
        return action()
    except Exception as exc:
        return {
            "state": "COLLECTION_FAILED",
            "stage": name,
            "error_type": type(exc).__name__,
            "error_message": str(exc)[:1600],
            "no_science_result_inferred": True,
        }


def _compact_radial_bins(bins):
    """Keep every radial coordinate/status, remove repeated verbose keys."""
    columns = (
        "a_mid_arcsec", "outer_arcsec", "pixels", "valid_pixels",
        "coverage", "sectors_with_data", "sectors_positive",
        "median_smoothed_rms", "supported",
    )
    samples = [[b.get(k) for k in columns] for b in bins]
    # Preserve octant contrast at important transitions and edge radii,
    # without repeating 8 floats at EVERY radial sample in every policy.
    special = set()
    if bins:
        special.update((0, len(bins) - 1))
    for i, bin_data in enumerate(bins):
        if i > 0 and bin_data.get("supported") != bins[i-1].get("supported"):
            special.update((i-1, i))
        if bin_data.get("coverage", 1.) < .7:
            special.add(i)
    contrasts = [{
        "radial_index": i,
        "a_mid_arcsec": bins[i].get("a_mid_arcsec"),
        "sector_contrast_by_octant": bins[i].get("sector_contrast_by_octant"),
    } for i in sorted(special)]
    return {
        "columns": list(columns),
        "rows": samples,
        "octant_contrasts_at_boundaries": contrasts,
        "all_radial_rows_included": True,
        "octant_contrasts_are_selected_not_complete": True,
    }


def _compact_ablation(result):
    if result is None:
        return None
    out = {k: v for k, v in result.items() if k != "scenarios"}
    out["scenarios"] = []
    for scenario in result.get("scenarios", []):
        item = dict(scenario)
        if "sga_radial_bins" in item:
            item["sga_radial_trace"] = _compact_radial_bins(
                item.pop("sga_radial_bins"))
        out["scenarios"].append(item)
    return out


def _compact_scene(workspace):
    stage = workspace.scene_editor
    result = stage._complete_scene_qa_snapshot()
    if result["project"]["project_id"] != workspace._bound_project_id:
        return {
            "state": "STALE_PROJECT",
            "expected_project_id": workspace._bound_project_id,
            "reported_project_id": result["project"]["project_id"],
            "visible_status": result.get("display_messages"),
        }
    result["mask_policy_ablation"] = _compact_ablation(
        result.get("mask_policy_ablation"))
    # Source-level candidates are essential for diagnosis: keep the full
    # source evidence metadata, not images, cuts, or pixel masks.
    result["state"] = (
        "WORKER_IN_PROGRESS" if result["worker"]["in_progress"] else
        "CURRENT_LOADED" if result.get("native_fits") else
        "NOT_LOADED")
    result["source_band_is_current"] = (
        result["controls"]["band"] == result["native_fits"]["band"]
        if result.get("native_fits") else None)
    result["ablation_radial_format"] = (
        "EVERY_BIN_COMPACT_COLUMNS_WITH_SELECTED_OCTANT_CONTRAST; "
        "full eight-octant rows available in Copy Scene QA")
    return result


def _stage_result(workspace, name):
    stage = workspace.sky_stage if name == "background" else workspace.aperture_stage
    # Do NOT call refresh_context: that method mutates controls and may
    # discard stale evidence. A snapshot is a read-only identity check.
    current = stage._snapshot()
    outcome = {
        "state": "NOT_RUN",
        "visible_status": stage.status.text(),
        "worker_in_progress": stage._future is not None,
        "current_source_ready": current is not None,
        "current_result_provenance_matches": None,
    }
    result = stage._result
    if stage._future is not None:
        outcome["state"] = "WORKER_IN_PROGRESS"
    elif result is not None and (
            current is None or current[0] != stage._result_identity):
        outcome["state"] = "STALE_RESULT_WITHHELD"
        outcome["current_result_provenance_matches"] = False
    elif result is not None:
        outcome["state"] = "CURRENT_DIAGNOSTIC"
        outcome["current_result_provenance_matches"] = True
        if (result.get("project_id") != workspace._bound_project_id or
                result.get("band") != current[0][2] or
                result.get("dataset_source") != current[0][1] or
                result.get("draft_sha256", result.get(
                    "scene_draft_sha256")) != current[0][3] or
                result.get("catalog_sha256") != current[0][4]):
            outcome["state"] = "STALE_RESULT_WITHHELD"
            outcome["current_result_provenance_matches"] = False
        else:
            outcome["result"] = (_compact_background(result) if name == "background"
                                 else _compact_apertures(result))
    if outcome["state"] == "NOT_RUN" and current is None:
        outcome["state"] = "BLOCKED_BY_SCENE_INPUTS"
    return outcome


def _compact_background(result):
    # The sky worker exports tile summaries and fitted nuisance parameters,
    # never sky-subtracted frames. Retain its science-veto evidence intact.
    # A projection is used only for verbose per-tile arrays, when present.
    return _guard_raw_arrays(result)


def _compact_apertures(result):
    out = {k: v for k, v in result.items() if k != "curves"}
    out["curves"] = []
    for curve in result.get("curves", []):
        c = {k: v for k, v in curve.items()
             if k not in ("cumulative", "independently_integrated_annuli")}
        c["cumulative"] = [
            {k: v for k, v in aperture.items()
             if k in (
                 "geometry_icrs", "geometry_native", "area_native_pixel2",
                 "area_by_reason_native_pixel2", "frame_invalid_by_reason_native_pixel2",
                 "signed_valid_flux_nMgy", "formal_diagonal_sigma_nMgy",
                 "ab_magnitude_valid_area_only",
                 "measured_mean_intensity_nMgy_per_arcsec2",
                 "science_state", "warnings", "recipe_sha256")}
            for aperture in curve["cumulative"]]
        c["independently_integrated_annuli"] = (
            curve["independently_integrated_annuli"])
        c["cumulative_sector_qa"] = [
            aperture.get("sector_qa") for aperture in curve["cumulative"]]
        out["curves"].append(c)
    out["omitted_from_each_cumulative_record"] = (
        "Repeated project/source identifiers (reported in parent); "
        "full original per-aperture records in Copy aperture QA")
    return out


def _guard_raw_arrays(value, depth=0):
    """Recursively serialize diagnostics, never export full pixel grids."""
    if depth > 24:
        return {"state": "NESTED_DEPTH_LIMIT", "full_report_available": True}
    if isinstance(value, np.ndarray):
        if value.ndim >= 2 or value.size > 256:
            return {"omitted_raw_array": True, "shape": list(value.shape),
                    "dtype": str(value.dtype)}
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(k): _guard_raw_arrays(v, depth+1)
                for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_guard_raw_arrays(v, depth+1) for v in value]
    if isinstance(value, Path):
        return str(value)
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return {"unserializable_type": type(value).__name__,
            "full_report_available": True}


def _compact_psf(workspace):
    from lsbpipe.gui.psf0_full_qa import detailed_band_report
    bands = {}
    for band, panel in workspace.psf0_panels.items():
        def one():
            report = detailed_band_report(panel)
            status = report["result_provenance_state"]
            is_bound = (
                report["source_project_uuid"] in (None, workspace._bound_project_id)
                and (report["project_root"] is None or
                     workspace._project_root is not None and
                     Path(report["project_root"]) == workspace._project_root))
            # Controls and star-review metadata may remain useful, but a
            # result from another UUID must never be labeled CURRENT.
            if not is_bound:
                return {"state": "STALE_RESULT_WITHHELD",
                        "reason": "PROJECT_OR_UUID_MISMATCH",
                        "result_project_id": report["source_project_uuid"],
                        "selected_project_id": workspace._bound_project_id}
            if status == "NO_SCIENCE_RESULT":
                state = "NOT_RUN"
            elif status == "RESULT_RECIPE_POLICY_MISMATCH":
                state = "STALE_RESULT_WITHHELD"
            else:
                state = "CURRENT_DIAGNOSTIC"
            if state == "STALE_RESULT_WITHHELD":
                report.pop("science_report_or_rejection", None)
                report.pop("stars_in_training_review", None)
            return {"state": state, "report": _guard_raw_arrays(report)}
        bands[band] = _safe_stage("psf0." + band, one)
    return {
        "state": "CURRENT_PROJECT_BAND_INDEPENDENT",
        "bands": bands,
        "batch_status": workspace.batch_status.text(),
    }


def _environment():
    return {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "desidera_version": _package_version("desidera"),
        "numpy": np.__version__,
        "pyside6": _package_version("PySide6"),
        "astropy": _package_version("astropy"),
        "photutils": _package_version("photutils"),
        "git_revision": _git_revision(),
    }


def _git_revision():
    """Read only tiny Git pointers when launched from a source checkout."""
    try:
        root = Path(__file__).resolve().parents[3]
        dotgit = root / ".git"
        if not dotgit.is_dir():
            return None
        head = (dotgit / "HEAD").read_text(encoding="utf-8").strip()
        if not head.startswith("ref: "):
            return head if len(head) == 40 else None
        ref = dotgit / head[5:].strip()
        sha = ref.read_text(encoding="utf-8").strip() if ref.is_file() else None
        return sha if sha is not None and len(sha) == 40 else None
    except (OSError, ValueError):
        return None


def build_project_diagnostic_bundle(workspace):
    """Capture current in-memory evidence without triggering any worker."""
    project_id = workspace._bound_project_id
    scene = workspace.scene_editor
    loaded = scene._loaded
    draft = scene._draft
    contaminant = scene._contamination
    navigation = workspace.current_navigation()
    report = {
        "schema": _SCHEMA,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scientific_adoption": False,
        "science_hold": True,
        "read_only_snapshot": True,
        "original_fits_modified": False,
        "workers_started_by_copy": False,
        "snapshot_scope": (
            "CURRENT_PROJECT_LOADED_SCENE_BAND_PLUS_INDEPENDENT_PSF_GRIZ"),
        "limitations": [
            "The copy button does not run missing or stale calculations",
            "Only the currently loaded Scene band has sky and aperture diagnostics",
            "Previously invalidated bands cannot be reconstructed from memory",
            "All stages remain provisional and unadopted",
            "Scene radial octants are selectively retained; full detail via Copy Scene QA",
            "No native FITS, CCD, PSF stamp, or masked pixel grids are copied",
        ],
        "environment": _environment(),
        "project": {
            "project_id": project_id,
            "project_directory_name": (
                workspace._project_root.name if workspace._project_root else None),
            "bound_dataset_identity": workspace._bound_source_identity,
            "current_stage": navigation[0],
            "navigation_band": navigation[1],
            "scene_band": scene._band,
            "loaded_dataset_band": loaded.get("band") if loaded else None,
            "loaded_dataset_source": loaded.get("source") if loaded else None,
            "loaded_scene_shape": (
                list(loaded.get("pixel_shape")) if loaded else None),
            "target_anchor_icrs": scene._project_target,
            "scene_draft_digest": draft.digest if draft else None,
            "catalogue_sha256": (
                contaminant.get("catalog_sha256") if contaminant else None),
        },
        "stages": {
            "psf0": _safe_stage("psf0", lambda: _compact_psf(workspace)),
            "scene_masks": _safe_stage("scene_masks",
                                       lambda: _compact_scene(workspace)),
            "background": _safe_stage("background",
                                      lambda: _stage_result(workspace, "background")),
            "aperture": _safe_stage("aperture",
                                    lambda: _stage_result(workspace, "aperture")),
        },
    }
    # Only fully JSON-serializable, finite numerical values are permissible.
    # Failing serialization must NEVER produce a misleading partial copy.
    report = _guard_raw_arrays(report)
    json.dumps(report, ensure_ascii=False, allow_nan=False)
    return report


def serialize_project_diagnostic_bundle(workspace):
    return json.dumps(build_project_diagnostic_bundle(workspace),
                      ensure_ascii=False, allow_nan=False,
                      separators=(",", ":"), sort_keys=True)
