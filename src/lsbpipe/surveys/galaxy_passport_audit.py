"""Conservative offline view of a saved *selected-dataset* FITS/WCS audit.

This is a read-only navigator to past diagnostics, NOT an independent current
SHA-256 verification. A quick stat match can only say the file metadata still
matches a receipt made when the original files were hashed and decoded.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
import math
from pathlib import Path
import re

from .dataset_integrity import SCHEMA as AUDIT_SCHEMA
from .legacy import BrickHit, product_for

_MAX_AUDIT_BYTES = 8 * 1024 * 1024
_SHA256 = re.compile(r"^[a-f0-9]{64}$")
_KINDS = ("image", "invvar", "nexp", "maskbits")


@dataclass(frozen=True)
class SavedAuditEvidence:
    state: str
    label: str
    explanation: str
    receipt_relative_path: str = ""
    created_utc: str = ""
    instrumentally_usable_fraction: float | None = None


def _result(state: str, label: str, explanation: str, *, path="",
            date="", usable=None) -> SavedAuditEvidence:
    return SavedAuditEvidence(state, label, explanation, path, date, usable)


def _receipt_path(project) -> tuple[Path, str]:
    pointer = project.last_dataset_audit
    if not isinstance(pointer, dict):
        raise ValueError("Invalid project audit pointer")
    rel = pointer.get("relative_path")
    if not isinstance(rel, str) or not rel or "\\" in rel:
        raise ValueError("Invalid audit receipt path")
    parts = Path(rel).parts
    # Fixed, existing on-disk layout. No traversal or alternative files.
    if (len(parts) != 3 or parts[0] != "dataset-audits"
            or not parts[1] or parts[1].startswith(".")
            or parts[2] != "audit.json" or ".." in parts):
        raise ValueError("Unexpected audit receipt path")
    path = project.root / Path(rel)
    if path.is_symlink() or path.parent.is_symlink():
        raise ValueError("Symlinked audit receipt not accepted")
    if not path.is_file() or path.stat().st_size > _MAX_AUDIT_BYTES:
        raise ValueError("Missing or oversized audit receipt")
    return path, rel


def _audit_data_root(project, report: dict) -> Path:
    if project.data_root is None:
        raise ValueError("No current survey data root")
    source = report.get("data_root")
    pointer_source = project.last_dataset_audit.get("data_root")
    if not isinstance(source, str) or not isinstance(pointer_source, str):
        raise ValueError("Invalid audit source directory")
    if source != pointer_source:
        raise ValueError("Audit pointer and receipt source differ")
    original = Path(source)
    if not original.is_absolute() or original.is_symlink() or not original.is_dir():
        raise ValueError("Original audited data root unavailable")
    current = project.data_root
    if not current.is_dir() or current.is_symlink():
        raise ValueError("Current data root unavailable")
    candidates = {current.resolve()}
    data = current / "data"
    if data.is_dir() and not data.is_symlink():
        candidates.add(data.resolve())
    if original.resolve() not in candidates:
        raise ValueError("Audit source is not the selected project's data root")
    return original


def _same_position(project, recorded: dict) -> bool:
    field = project.field
    if field is None or not isinstance(recorded, dict):
        return False
    for key, current in (("ra_deg", field.ra_deg),
                         ("dec_deg", field.dec_deg),
                         ("width_arcmin", field.width_arcmin)):
        value = recorded.get(key)
        if isinstance(value, bool) or not isinstance(value, (float, int)):
            return False
        if not math.isfinite(value) or abs(value - current) > 1e-9:
            return False
    return True


def _checked_product_stats(data_root: Path, dataset: dict) -> bool:
    """Allow *metadata agreement* only; never claim hashes rechecked."""
    band = dataset["band"]
    hit = BrickHit(dataset["release"], dataset["region"],
                   dataset["brick"], 0, 0, 0, 0, {band: None})
    products = dataset.get("products")
    if not isinstance(products, dict):
        return False
    for kind in _KINDS:
        recorded = products.get(kind)
        if not isinstance(recorded, dict):
            return False
        relative = product_for(
            hit, kind, None if kind == "maskbits" else band).relative_path
        if recorded.get("relative_path") != relative.as_posix():
            return False
        size, mtime, digest = (
            recorded.get("size_bytes"), recorded.get("mtime_ns"),
            recorded.get("sha256_local"))
        if (type(size) is not int or size <= 0 or type(mtime) is not int
                or mtime <= 0 or not isinstance(digest, str)
                or _SHA256.fullmatch(digest) is None):
            return False
        candidate = data_root / relative
        try:
            if (candidate.is_symlink() or not candidate.is_file()
                    or (candidate.stat().st_size, candidate.stat().st_mtime_ns)
                    != (size, mtime)):
                return False
        except OSError:
            return False
    return True


def inspect_saved_dataset_audit(project) -> SavedAuditEvidence:
    """Return historical audit status without reading original FITS pixels."""
    if not project.last_dataset_audit:
        return _result(
            "NO_SAVED_AUDIT", "No saved FITS/WCS audit",
            "Run Dataset QA in the galaxy workspace to inspect original "
            "FITS structure, WCS and numerical masks.")
    try:
        path, rel = _receipt_path(project)
        # A local JSON receipt is not cryptographically authenticated. It
        # remains *historical diagnostic evidence* even with valid shape.
        report = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(report, dict):
            raise ValueError("Invalid audit JSON document")
        if (report.get("schema") != AUDIT_SCHEMA
                or report.get("project_id") != project.project_id
                or report.get("science_status") != "DIAGNOSTIC_NOT_ADOPTED"):
            raise ValueError("Audit has wrong schema/project/science status")
        pointer = project.last_dataset_audit
        if (report.get("created_utc") != pointer.get("created_utc")
                or not isinstance(report.get("datasets"), list)):
            raise ValueError("Audit timestamp or dataset list differs")
        data_root = _audit_data_root(project, report)
        if not _same_position(project, report.get("field")):
            return _result(
                "AUDIT_SELECTION_STALE", "Saved audit · target changed",
                "The current project target/field differs from the audit's "
                "sampled field. Re-run QA; old coverage cannot be adopted.",
                path=rel)
        key = project.selected_dataset
        if key is None:
            return _result(
                "AUDIT_SELECTION_STALE", "Saved audit · no dataset selected",
                "Select a dataset and run an audit again.", path=rel)
        matching = [d for d in report["datasets"]
                    if isinstance(d, dict) and
                    tuple(d.get(k) for k in (
                        "release", "region", "brick", "band")) == tuple(key)]
        if len(matching) != 1:
            return _result(
                "AUDIT_SELECTION_STALE", "Selected band not in saved audit",
                "The saved audit does not match the currently selected "
                "release, region, brick and band.", path=rel)
        selected = matching[0]
        status = selected.get("status")
        created = report.get("created_utc", "")
        if status != "FITS_WCS_VALIDATED_DIAGNOSTIC":
            # Invalid/incomplete reports contain no independent complete
            # input-stat evidence: show only their historical outcome.
            return _result(
                "AUDIT_REPORTED_ISSUE",
                "Previous audit · " + str(status or "unrecognized outcome"),
                "The last saved diagnostic did not validate this quartet. "
                "This is historical, not a verified current-data failure.",
                path=rel, date=created)
        if not _checked_product_stats(data_root, selected):
            return _result(
                "AUDIT_INPUTS_CHANGED", "Saved FITS/WCS audit · needs rerun",
                "At least one original file is missing, moved, or has "
                "different file metadata; the old validation cannot "
                "describe current inputs.", path=rel, date=created)
        fractions = selected.get("fractions", {})
        usable = fractions.get("instrumentally_usable") if isinstance(
            fractions, dict) else None
        if (isinstance(usable, bool) or not isinstance(usable, (float, int))
                or not math.isfinite(usable) or not 0 <= usable <= 1):
            usable = None
        return _result(
            "AUDIT_RECORDED_STAT_MATCH",
            "FITS/WCS checked in saved audit · metadata unchanged",
            "The selected quartet passed the saved FITS/WCS diagnostic. "
            "Size and modification-time metadata still match, but content "
            "hashes were NOT rechecked at startup. Sky, PSF, coadd noise "
            "and scientific suitability are NOT certified.",
            path=rel, date=created, usable=usable)
    except (OSError, ValueError, TypeError, KeyError, AttributeError):
        return _result(
            "AUDIT_EVIDENCE_UNAVAILABLE", "Saved audit · cannot be verified",
            "The audit receipt is unavailable, malformed or belongs to a "
            "different data source. Do not infer FITS validity; run QA again.")
