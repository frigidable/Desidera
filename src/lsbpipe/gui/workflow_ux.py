"""Presentation-only workflow guidance for the three essential galaxy steps.

Does not grant science acceptance, mutate a project, infer a catalog
selection, or start work. It reflects explicit scientific readiness.
"""
from __future__ import annotations
from dataclasses import dataclass
import json


@dataclass(frozen=True)
class WorkflowCue:
    step: int
    title: str
    instruction: str

    def __post_init__(self):
        if self.step not in (1, 2, 3):
            raise ValueError("Invalid workflow stage")


def workflow_cue(*, has_target: bool, has_catalog: bool,
                 selected_brick: bool, has_local_fits: bool,
                 has_interactive_image: bool) -> WorkflowCue:
    """Derive an honest next action from current UI, without inference."""
    if not has_target:
        return WorkflowCue(1, "01 / 03 · Target",
                           "Enter ICRS coordinates, then find survey data.")
    if not has_catalog and not has_local_fits:
        return WorkflowCue(1, "01 / 03 · Target",
                           "Find data to see which survey bricks cover this field.")
    if not selected_brick:
        return WorkflowCue(2, "02 / 03 · Data",
                           "Select an available brick to continue.")
    if not has_local_fits:
        return WorkflowCue(2, "02 / 03 · Data",
                           "Choose the four science FITS products, then download.")
    if not has_interactive_image:
        return WorkflowCue(3, "03 / 03 · Inspect",
                           "Original FITS found. Open the selected band to inspect.")
    return WorkflowCue(3, "03 / 03 · Inspect",
                       "Explore FITS pixels and masks; diagnostics are optional.")


# Copyable support reports intentionally use an allowlist. Worker exception
# text, traceback, source file paths, project labels/UUIDs, and credentials
# can appear in raw Qt status or private _last_failure diagnostics. NONE of
# those unreviewed fields belong on a shared clipboard by default.
_GUI_ERROR_CAUSES = {
    "FileNotFoundError": "Required local source file is missing.",
    "PermissionError": "Local filesystem access was denied.",
    "IsADirectoryError": "Expected a file, but the path resolves to a directory.",
    "OSError": "Local I/O or FITS access failed.",
    "ValueError": "The source evidence or numerical inputs failed validation.",
    "KeyError": "Required structured input fields are missing.",
    "TypeError": "An input or result has an unexpected type.",
    "MemoryError": "Available memory was insufficient.",
    "ImportError": "A required Python dependency could not be imported.",
    "ArithmeticError": "A numerical operation failed.",
    "RuntimeError": "The background task or its execution environment failed.",
}
_GUI_ERROR_STAGES = {
    "scene": (
        "02 / Scene and Masks", "load_original_fits",
        "Check the selected band and local IMAGE/INVVAR/NEXP/MASKBITS "
        "quartet; then choose Retry FITS load. Draft regions are retained."),
    "background": (
        "03 / Residual background", "evaluate_residual_sky",
        "Check Stage 02 native FITS, prepared masks and sky sampling; "
        "then choose Evaluate residual sky again."),
    "aperture": (
        "04 / Aperture Photometry", "build_aperture_curves",
        "Check Stage 02 source, geometry, masks and photutils dependency; "
        "then choose Build aperture curves again."),
}
_GUI_FAILURE_STATUSES = frozenset((
    "FITS_LOAD_FAILED_CLOSED", "WORKER_FAILED_CLOSED",
    "WORKER_START_FAILED_CLOSED",
))


def format_gui_failure_report(*, stage: str, error_type: str,
                              failure_status: str, band=None) -> str:
    """Share a deterministic, path-free diagnostic suitable for support.

    Only the four explicitly safe parameters contribute to the report.
    Never accept an exception, traceback or arbitrary message as input.
    This report is NOT a log dump or a scientific-quality assessment.
    """
    if stage not in _GUI_ERROR_STAGES:
        raise ValueError("Unrecognized GUI stage for support report")
    if failure_status not in _GUI_FAILURE_STATUSES:
        raise ValueError("Unrecognized failure status for support report")
    title, operation, recovery = _GUI_ERROR_STAGES[stage]
    safe_type = error_type if error_type in _GUI_ERROR_CAUSES else "Exception"
    reason = _GUI_ERROR_CAUSES.get(
        safe_type, "Unclassified operation error; review local diagnostics.")
    if failure_status == "WORKER_START_FAILED_CLOSED":
        reason = "The background worker could not be started."
    payload = {
        "report_format": "DESIDERA_GUI_FAILURE_REPORT_V1",
        "stage": title,
        "operation": operation,
        "status": failure_status,
        "error_type": safe_type,
        "possible_cause": reason,
        "recommended_action": recovery,
        "band": band if band in ("g", "r", "i", "z") else None,
        "science_status": "UNADOPTED_DIAGNOSTIC_ONLY",
        "raw_error_message_included": False,
        "traceback_included": False,
        "local_paths_or_project_identifiers_included": False,
    }
    return json.dumps(payload, ensure_ascii=False, indent=2)
