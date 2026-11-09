"""Read-only readiness snapshot for the future object-bound photometry workspace.

Presence is not FITS validation, WCS integrity, scientifically usable pixels,
a background estimate, or permission to publish a photometric result.
"""
from __future__ import annotations

from pathlib import Path

from lsbpipe.surveys.dataset_integrity import _expected_paths
from lsbpipe.surveys.legacy import BRICK_RE, regions_for_release
from lsbpipe.surveys.object_projects import read_object_project

KINDS = ("image", "invvar", "nexp", "maskbits")
SCIENCE_STATE = "NOT_MEASURED_NOT_ADOPTED"


def inspect_photometry_inputs(project_root: str | Path) -> dict:
    """Inspect the *saved* selected native quartet without opening or writing FITS.

    This intentionally does not use asset-index presence as scientific
    validation. Non-native and cross-brick observations are future adapters.
    """
    project = read_object_project(project_root)
    result = {
        "project_id": project.project_id,
        "galaxy": project.display_name,
        "target_defined": project.field is not None,
        "selected_dataset": project.selected_dataset,
        "data_root": str(project.data_root) if project.data_root else None,
        "products": (),
        "present_count": 0,
        "science_state": SCIENCE_STATE,
        "status": "NO_SELECTED_DATASET",
        "reason": "Select an observational dataset in Galaxy Workspace.",
    }
    key = project.selected_dataset
    if key is None:
        return result
    release, region, brick, band = key
    if (release not in ("dr9", "dr10", "dr11")
            or region not in regions_for_release(release)
            or not BRICK_RE.fullmatch(brick)
            or band not in ("g", "r", "i", "z")
            or (region == "north" and band == "i")):
        result.update(status="UNSUPPORTED_SELECTION",
                      reason="Selected dataset is not a supported native optical quartet.")
        return result
    root = project.data_root
    if root is None or not root.is_dir():
        result.update(status="DATA_ROOT_UNAVAILABLE",
                      reason="The saved original-FITS data directory is unavailable.")
        return result
    paths = _expected_paths(root, key)
    products = tuple({
        "kind": kind,
        "path": str(paths[kind]),
        "present": paths[kind].is_file() and not paths[kind].is_symlink(),
    } for kind in KINDS)
    count = sum(item["present"] for item in products)
    result.update(products=products, present_count=count)
    if count != len(KINDS):
        result.update(
            status="INCOMPLETE_QUARTET",
            reason="Some original FITS are missing or symlinks. "
                   "No photometric preparation was performed.")
    else:
        result.update(
            status="FILES_PRESENT_UNVERIFIED",
            reason="All four original paths exist. FITS data, WCS, units, "
                   "pixel quality, PSF, sky and noise remain unvalidated.")
    return result
