"""Lightweight, read-only observational readiness for Galaxy Projects.

A directory, cached JPEG, dataset registry record or saved QA reference is NOT
proof of usable calibrated FITS. This projection checks only a bounded quartet
of *specific* paths for the currently selected original Legacy brick/band.
It never scans directories, opens FITS, hashes assets or adopts science.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .legacy import BrickHit, product_for

KINDS = ("image", "invvar", "nexp", "maskbits")


@dataclass(frozen=True)
class PassportReadiness:
    code: str
    observation_label: str
    explanation: str
    next_action: str
    next_action_label: str
    file_count: int = 0
    expected_count: int = 0
    missing_kinds: tuple[str, ...] = ()


def _status(code, label, explanation, action, button, *, file_count=0,
            expected_count=0, missing_kinds=()):
    return PassportReadiness(
        code, label, explanation, action, button, file_count,
        expected_count, tuple(missing_kinds))


def inspect_passport_readiness(project) -> PassportReadiness:
    """Inspect bounded selected-brick file presence, not FITS/QA integrity."""
    if project.field is None:
        return _status(
            "TARGET_NOT_SET", "Target position not defined",
            "Define accepted project ICRS coordinates before searching surveys.",
            "DEFINE_TARGET", "Open project · define target")
    if project.selected_dataset is None:
        return _status(
            "DATASET_NOT_SELECTED", "No survey dataset selected",
            "The project exists, but no survey release, brick and band have "
            "been selected. Select observations in the workspace.",
            "FIND_DATA", "Open project · find data")
    if project.data_root is None:
        return _status(
            "DATA_ROOT_NOT_LINKED", "Survey selected · no local FITS location",
            "A survey brick is selected, but original FITS are not linked to "
            "this galaxy project. A colour Viewer preview is not FITS evidence.",
            "GET_FITS", "Open project · get original FITS")
    root = Path(project.data_root)
    if not root.is_dir():
        return _status(
            "DATA_ROOT_MISSING", "Saved FITS location unavailable",
            "The stored data-directory path is missing or inaccessible. "
            "Relink the original files; no file integrity claim is made.",
            "RELINK_DATA", "Open project · relink FITS")
    release, region, brick, band = project.selected_dataset
    try:
        hit = BrickHit(release, region, brick, 0, 0, 0, 0, {band: None})
        paths = [
            (kind, product_for(
                hit, kind, None if kind == "maskbits" else band).relative_path)
            for kind in KINDS
        ]
    except (ValueError, KeyError, TypeError) as exc:
        return _status(
            "DATASET_INVALID", "Selected dataset needs attention",
            "The selected survey key cannot be mapped to original FITS "
            "products; inspect the project configuration.",
            "REVIEW_DATASET", "Open project · review dataset")
    available = []
    for kind, relative in paths:
        # Imported archives may reside under the project root or its data/
        # child. Never glob unknown paths or follow symlinks.
        candidates = (root / relative, root / "data" / relative)
        try:
            exists = any(
                not item.is_symlink() and item.is_file()
                and item.stat().st_size > 0
                for item in candidates)
        except OSError:
            exists = False
        if exists:
            available.append(kind)
    missing = tuple(kind for kind in KINDS if kind not in available)
    count = len(available)
    if count == 0:
        return _status(
            "NO_MATCHING_ORIGINAL_FITS",
            "No original FITS found for the selected brick/band",
            "The stored data path exists, but this survey selection has no "
            "matching original IMAGE/INVVAR/NEXP/MASKBITS files.",
            "GET_FITS", "Open project · get original FITS",
            file_count=0, expected_count=4, missing_kinds=missing)
    if count < len(KINDS):
        return _status(
            "INCOMPLETE_ORIGINAL_FITS",
            f"Original FITS present: {count}/4 products",
            "Missing " + ", ".join(missing).upper() + ". File existence "
            "does not establish data integrity or photometric suitability.",
            "COMPLETE_FITS", "Open project · complete FITS",
            file_count=count, expected_count=4, missing_kinds=missing)
    return _status(
        "FOUR_ORIGINAL_FITS_PRESENT_UNVERIFIED",
        "Original FITS present: 4/4 products",
        "All four files exist for the selected brick and filter. FITS/WCS "
        "validity, correlated noise, PSF and sky have not been established "
        "by this passport.",
        "REVIEW_FITS", "Open project · inspect FITS",
        file_count=4, expected_count=4)
