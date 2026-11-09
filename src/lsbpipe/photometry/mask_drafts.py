"""Explicit, crash-safe, revision-preserving storage of manual mask drafts.

Only user-initiated saving changes project storage. Immutable historical
snapshots remain; the mutable HEAD pointer changes atomically. No FITS
or adopted scientific measurements are written.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import tempfile

from lsbpipe.surveys.object_projects import read_object_project
from .sky_masks import SkyMaskDraft

HEAD_SCHEMA = "desidera.sky-mask-head/1"


def _folder(project_root: str | Path) -> tuple[Path, str]:
    project = read_object_project(project_root)
    folder = project.root / "analyses" / "mask-drafts"
    if folder.is_symlink():
        raise ValueError("Mask draft folder may not be a symlink")
    return folder, project.project_id


def _write_atomic(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp = tempfile.mkstemp(prefix=".mask-partial-", dir=str(path.parent))
    try:
        with os.fdopen(fd, "wb") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temp, path)
    finally:
        if os.path.exists(temp):
            os.unlink(temp)


def _json_bytes(data: dict) -> bytes:
    return json.dumps(data, ensure_ascii=False, sort_keys=True,
                      allow_nan=False, separators=(",", ":")).encode("utf-8")


def read_active_mask_draft(project_root: str | Path) -> tuple[SkyMaskDraft, str | None]:
    folder, project_id = _folder(project_root)
    head = folder / "HEAD.json"
    if not head.exists():
        return SkyMaskDraft(project_id), None
    if head.is_symlink() or not head.is_file():
        raise ValueError("Unsafe mask draft HEAD")
    record = json.loads(head.read_text("utf-8"))
    digest = record.get("digest")
    if (record.get("schema") != HEAD_SCHEMA or
            record.get("project_id") != project_id or
            not isinstance(digest, str) or len(digest) != 64 or
            any(c not in "0123456789abcdef" for c in digest)):
        raise ValueError("Invalid mask draft HEAD")
    snapshot = folder / ("draft-" + digest + ".json")
    if snapshot.is_symlink() or not snapshot.is_file():
        raise ValueError("Mask draft snapshot missing or unsafe")
    draft = SkyMaskDraft.from_dict(json.loads(snapshot.read_text("utf-8")))
    if draft.project_id != project_id or draft.digest != digest:
        raise ValueError("Modified/cross-project mask draft snapshot")
    return draft, digest


def save_mask_draft(project_root: str | Path, draft: SkyMaskDraft,
                    *, expected_digest: str | None) -> str:
    """Save a new immutable draft and update HEAD only on expected revision.

    Caller must supply last read digest. This prevents silent lost updates
    from simultaneous windows/restarts and never overwrites historical edits.
    """
    folder, project_id = _folder(project_root)
    if draft.project_id != project_id:
        raise ValueError("Mask draft belongs to a different galaxy")
    _, current = read_active_mask_draft(project_root)
    if current != expected_digest:
        raise ValueError("Mask draft modified elsewhere: refresh before saving")
    folder.mkdir(parents=True, exist_ok=True)
    digest = draft.digest
    snapshot = folder / ("draft-" + digest + ".json")
    payload = _json_bytes(draft.to_dict())
    if snapshot.exists():
        if snapshot.is_symlink() or snapshot.read_bytes() != payload:
            raise ValueError("Existing mask snapshot differs from its recorded identity")
    else:
        _write_atomic(snapshot, payload)
    # Avoid churn when unchanged, but still validate the current HEAD.
    if digest != current:
        _write_atomic(folder / "HEAD.json", _json_bytes({
            "schema": HEAD_SCHEMA, "project_id": project_id, "digest": digest,
            "science_state": "MASK_DRAFT_NOT_ADOPTED",
        }))
    return digest
