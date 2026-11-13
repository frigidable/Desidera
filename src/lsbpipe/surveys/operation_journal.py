"""Crash-tolerant, non-scientific journal for workspace jobs (P0.4).

A journal entry says what the UI attempted. It is NEVER evidence of a
completed scientific calculation; only QA manifests and input hashes can
establish that. Unknown jobs left by process death are shown on next open
without silently replaying, deleting or adopting their outputs.
"""
from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import json
import re
import secrets
from typing import Literal

from .local_project import _atomic_json
from .project_locking import exclusive_metadata_lock

SCHEMA = "lsbpipe.operation-journal/1"
LIVE = frozenset(("RUNNING", "STOP_REQUESTED"))
TERMINAL = frozenset(("COMPLETED", "CANCELLED", "FAILED", "IGNORED"))
NAME = re.compile(r"^[A-Za-z][A-Za-z0-9_.:-]{0,79}$")
ID = re.compile(r"^[0-9]{8}T[0-9]{12}Z-[0-9a-f]{10}$")


def _utc():
    return datetime.now(timezone.utc).isoformat()


def _folder(project_root: str | Path) -> Path:
    root = Path(project_root).expanduser().resolve()
    if not root.is_dir():
        raise FileNotFoundError(root)
    path = root / "operations"
    if path.is_symlink():
        raise ValueError("Unsafe symlink at project operation journal")
    return path


def _journal_transaction(directory: Path):
    """One stable per-project lock for all job state transitions.

    The non-JSON sidecar is persistent and never interpreted as a job.
    """
    if directory.is_symlink():
        raise ValueError("Unsafe project operation journal directory")
    return exclusive_metadata_lock(directory / ".operation-journal.lock")


def _read(path: Path) -> dict:
    if path.is_symlink():
        raise ValueError("Unsafe operation journal symlink")
    data = json.loads(path.read_text(encoding="utf-8"))
    if (not isinstance(data, dict) or data.get("schema") != SCHEMA or
            data.get("operation_id") != path.stem or
            not ID.fullmatch(path.stem) or data.get("state") not in
            LIVE | TERMINAL or not isinstance(data.get("project_id"), str)):
        raise ValueError("Corrupt or foreign operation journal entry")
    return data


def start_operation(project_root: str | Path, project_id: str,
                    action: str) -> Path:
    """Create a unique independent journal; existing uncertain jobs survive."""
    from .object_projects import read_object_project
    if not isinstance(action, str) or not NAME.fullmatch(action):
        raise ValueError("Invalid operation kind")
    root = Path(project_root).expanduser().resolve()
    passport = read_object_project(root)
    if passport.project_id != project_id:
        raise ValueError("Operation journal UUID does not match project")
    directory = _folder(root)
    directory.mkdir(parents=True, exist_ok=True)
    with _journal_transaction(directory):
        ident = (datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
                 + "-" + secrets.token_hex(5))
        path = directory / f"{ident}.json"
        if path.exists():
            raise FileExistsError("Operation ID collision")
        _atomic_json(path, {
            "schema": SCHEMA, "operation_id": ident,
            "project_id": project_id, "action": action,
            "state": "RUNNING", "started_utc": _utc(),
            "updated_utc": _utc(),
            "scientific_status": "NOT_SCIENTIFIC_EVIDENCE",
            "message": "UI operation attempted; may be interrupted by process exit",
        })
        return path


def request_stop(path: str | Path) -> None:
    file = Path(path)
    with _journal_transaction(file.parent):
        info = _read(file)
        if info["state"] in TERMINAL:
            return
        info["state"] = "STOP_REQUESTED"
        info["updated_utc"] = _utc()
        _atomic_json(file, info)


def finish_operation(path: str | Path, outcome: str) -> None:
    """Record observed GUI terminal delivery, not scientific completion."""
    if outcome not in TERMINAL:
        raise ValueError("Invalid GUI operation outcome")
    file = Path(path)
    with _journal_transaction(file.parent):
        info = _read(file)
        if info["state"] in TERMINAL:
            raise ValueError("A terminal operation journal is immutable")
        info["state"] = outcome
        info["updated_utc"] = _utc()
        info["message"] = (
            "GUI callback handled; inspect immutable QA manifests separately"
            if outcome == "COMPLETED" else
            "No new scientific result was automatically adopted")
        _atomic_json(file, info)


def interrupted_operations(project_root: str | Path) -> list[dict]:
    """Read-only, bounded scan. Does not change an uncertain state."""
    directory = _folder(project_root)
    if not directory.exists():
        return []
    matches = []
    for path in directory.glob("*.json"):
        if path.is_symlink():
            continue
        try:
            info = _read(path)
        except (OSError, ValueError, TypeError):
            continue
        if info["state"] in LIVE:
            matches.append(info)
    return sorted(matches, key=lambda x: x["operation_id"], reverse=True)[:100]
