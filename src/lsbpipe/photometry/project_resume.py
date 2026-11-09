"""Small atomic, project-scoped UI resume pointer; never a science result.

Separate from the immutable scientific snapshots and object manifest to
avoid rewriting older project choices when only the active tab changes.
No new directories and no network activity.
"""
from __future__ import annotations

import json
from pathlib import Path

from lsbpipe.surveys.local_project import _atomic_json
from lsbpipe.surveys.object_projects import read_object_project

SCHEMA = "desidera.photometry-workspace-state/1"
FILENAME = "photometry-workspace-state.json"
STAGES = frozenset(("psf0", "scene"))
LOCATIONS = frozenset(("survey", "photometry"))
BANDS = frozenset(("g", "r", "i", "z"))


def _default(project):
    return {
        "location": "survey",
        "stage": "psf0",
        "band": (project.selected_dataset[3]
                 if project.selected_dataset is not None else "g"),
        "source_matches": True,
    }


def read_project_resume(root: str | Path) -> dict:
    """Fail closed on corrupt state, foreign project or changed brick."""
    project = read_object_project(root)
    path = project.root / FILENAME
    if not path.exists() and not path.is_symlink():
        return _default(project)
    if path.is_symlink() or not path.is_file():
        raise ValueError("Unsafe project resume pointer")
    try:
        record = json.loads(path.read_text(encoding="utf-8"))
    except (ValueError, OSError) as exc:
        raise ValueError("Invalid project resume pointer") from exc
    if (not isinstance(record, dict) or record.get("schema") != SCHEMA
            or record.get("project_id") != project.project_id
            or record.get("location") not in LOCATIONS
            or record.get("stage") not in STAGES
            or record.get("band") not in BANDS
            or not isinstance(record.get("source"), list)
            or any(not isinstance(v, str) for v in record["source"])):
        raise ValueError("Invalid or foreign project resume state")
    recorded = record["source"]
    selected = (list(project.selected_dataset[:3])
                if project.selected_dataset is not None else [])
    if recorded != selected:
        return {**_default(project), "source_matches": False}
    if project.selected_dataset is not None and (
            project.selected_dataset[1] == "north"
            and record["band"] == "i"):
        return {**_default(project), "source_matches": False}
    return {
        "location": record["location"],
        "stage": record["stage"],
        "band": record["band"],
        "source_matches": True,
    }


def save_project_resume(root: str | Path, *, stage: str, band: str,
                        location: str = "photometry") -> None:
    """Write one small atomic UI pointer, without mutating science/project."""
    project = read_object_project(root)
    if stage not in STAGES or band not in BANDS or location not in LOCATIONS:
        raise ValueError("Unsupported stage or optical band")
    if (project.selected_dataset is not None
            and project.selected_dataset[1] == "north" and band == "i"):
        raise ValueError("DR11 North has no native i-band")
    path = project.root / FILENAME
    if path.is_symlink() or (path.exists() and not path.is_file()):
        raise ValueError("Unsafe project resume path")
    source = (list(project.selected_dataset[:3])
              if project.selected_dataset is not None else [])
    payload = {
        "schema": SCHEMA,
        "project_id": project.project_id,
        "source": source,
        "location": location,
        "stage": stage,
        "band": band,
    }
    if path.is_file():
        try:
            old = json.loads(path.read_text(encoding="utf-8"))
        except (ValueError, OSError):
            old = None
        if old == payload:
            return
    _atomic_json(path, payload)
