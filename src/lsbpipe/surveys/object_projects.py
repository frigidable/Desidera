"""Portable, object-centric research projects and a user-local project library.

The original FITS, legacy local-project metadata and immutable QA runs are
never overwritten by this module. All operations are offline and atomic.
"""
from __future__ import annotations

from dataclasses import dataclass, field as dc_field
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import uuid

from .legacy import SkyField
from .project_locking import exclusive_metadata_lock
from .local_project import (
    _atomic_json, open_local_project, session_pointer_path,
    last_project_directory,
)

SCHEMA = "lsbpipe.object-project/1"
REGISTRY_SCHEMA = "lsbpipe.project-registry/1"
PROJECT_FILENAME = "lsbpipe-object.json"
REGISTRY_FILENAME = "projects.json"


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def library_index_path() -> Path:
    """Dedicated registry; does not replace the legacy last-project pointer."""
    return session_pointer_path().parent / REGISTRY_FILENAME


def _registry_transaction():
    """One stable cross-process lock for each user-local project library.

    Do not delete the lock sidecar. Its inode is the synchronization
    identity; cleanup while another process waits could split the lock.
    """
    index = library_index_path()
    return exclusive_metadata_lock(
        index.with_name(index.name + ".lock"))


def _read_json(path: Path) -> dict:
    try:
        item = json.loads(path.read_text(encoding="utf-8"))
    except (ValueError, OSError) as exc:
        raise ValueError(f"Cannot read project metadata at {path}: {exc}") from exc
    if not isinstance(item, dict):
        raise ValueError(f"Project metadata must be a JSON object: {path}")
    return item


def _read_manifest_with_sha(path: Path) -> tuple[dict, str]:
    """Read an object passport and authenticate the exact parsed JSON bytes.

    This is a concurrency revision token, not an original-FITS checksum
    or scientific acceptance receipt. Never follow an external symlink.
    """
    if path.is_symlink():
        raise ValueError("Unsafe symlinked object project manifest")
    try:
        payload = path.read_bytes()
        record = json.loads(payload.decode("utf-8"))
    except (OSError, ValueError, UnicodeError) as exc:
        raise ValueError(f"Cannot read project metadata at {path}: {exc}") from exc
    if not isinstance(record, dict):
        raise ValueError(f"Project metadata must be a JSON object: {path}")
    return record, hashlib.sha256(payload).hexdigest()


def _field_from_payload(data: object) -> SkyField | None:
    if data is None:
        return None
    if not isinstance(data, dict):
        raise ValueError("Invalid stored ICRS coordinates")
    return SkyField(float(data["ra_deg"]), float(data["dec_deg"]),
                    float(data["width_arcmin"]))


def _field_dict(value: SkyField | None) -> dict | None:
    if value is None:
        return None
    return {"ra_deg": value.ra_deg, "dec_deg": value.dec_deg,
            "width_arcmin": value.width_arcmin}


def _project_id(raw: str) -> str:
    try:
        result = uuid.UUID(raw)
    except (ValueError, TypeError, AttributeError) as exc:
        raise ValueError("Invalid immutable project UUID") from exc
    if result.version != 4:
        raise ValueError("Project ID must be UUIDv4")
    return str(result)


@dataclass
class ObjectProject:
    project_id: str
    display_name: str
    root: Path
    field: SkyField | None = None
    aliases: list[str] = dc_field(default_factory=list)
    selected_dataset: tuple[str, str, str, str] | None = None
    data_root: Path | None = None
    target_source: str = "user"
    created_at: str = dc_field(default_factory=_now)
    last_opened_at: str = dc_field(default_factory=_now)
    last_science: dict = dc_field(default_factory=dict)
    catalog_match: dict | None = None
    last_dataset_audit: dict = dc_field(default_factory=dict)
    target_radius_arcsec: float = 60.
    sky_inner_arcsec: float = 90.
    sky_outer_arcsec: float = 180.
    # Private JSON revision from an actually read/published passport.
    # These fields are NEVER serialized to the on-disk science schema.
    manifest_sha256: str | None = dc_field(
        default=None, repr=False, compare=False)
    manifest_observed: bool = dc_field(
        default=False, repr=False, compare=False)

    @property
    def status(self) -> str:
        if self.field is None:
            return "NEEDS_TARGET_POSITION"
        if self.data_root is None:
            return "TARGET_DEFINED_NO_FITS"
        return "LOCAL_DATA_NOT_SCIENTIFICALLY_ADOPTED"


def _validate(project: ObjectProject) -> None:
    project.project_id = _project_id(project.project_id)
    if not isinstance(project.display_name, str) or not project.display_name.strip():
        raise ValueError("Project name cannot be empty")
    if len(project.display_name) > 140 or any(
            ord(c) < 32 or c in r'\/:*?"<>|' for c in project.display_name):
        raise ValueError("Project name contains unsupported characters")
    project.display_name = project.display_name.strip()
    if project.selected_dataset is not None and (
            len(project.selected_dataset) != 4 or
            any(not isinstance(x, str) for x in project.selected_dataset)):
        raise ValueError("Invalid selected survey dataset")
    if project.catalog_match is not None and not isinstance(
            project.catalog_match, dict):
        raise ValueError("Invalid external catalog association")
    if not isinstance(project.last_dataset_audit, dict):
        raise ValueError("Invalid dataset audit index")
    if not all(math.isfinite(v) for v in (
            project.target_radius_arcsec, project.sky_inner_arcsec,
            project.sky_outer_arcsec)) or not (
            0 < project.target_radius_arcsec <= project.sky_inner_arcsec
            < project.sky_outer_arcsec <= 1200):
        raise ValueError("Invalid project radii: require 0 < target <= inner < outer <= 1200")


def project_path(root: str | Path) -> Path:
    return Path(root).expanduser().resolve() / PROJECT_FILENAME


def save_object_project(project: ObjectProject) -> Path:
    """Publish only a revision-observed object passport, never stale data.

    A passport owns target coordinates, aliases, dataset selection and QA
    references. The stable OS lock guards read/check/atomic-replace, while
    the private full-file SHA catches a stale independent GUI instance.
    """
    _validate(project)
    project.root = project.root.expanduser().resolve()
    if not project.root.is_dir():
        raise FileNotFoundError(f"Project directory missing: {project.root}")
    path = project_path(project.root)
    with exclusive_metadata_lock(path.with_name(path.name + ".lock")):
        if path.is_symlink():
            raise ValueError("Unsafe symlinked object project manifest")
        if path.exists():
            current, current_sha = _read_manifest_with_sha(path)
            if (current.get("schema") != SCHEMA or
                    _project_id(current.get("project_id")) != project.project_id):
                raise ValueError(
                    "Refusing to overwrite a different project identity")
        else:
            current_sha = None
        if not project.manifest_observed:
            if current_sha is not None:
                raise RuntimeError(
                    "Object project passport was never opened at this revision; "
                    "reopen the galaxy project before saving")
        elif current_sha != project.manifest_sha256:
            raise RuntimeError(
                "Object project passport changed in another process; "
                "reopen the galaxy project before saving")
        root = project.root
        # Store data_root relative only if within project root. External paths
        # remain absolute and must be relinked if moved.
        data_root = project.data_root.resolve() if project.data_root else None
        if data_root is None:
            data_link = None
        elif data_root.is_relative_to(root):
            data_link = {"relative": data_root.relative_to(root).as_posix()}
        else:
            data_link = {"absolute": str(data_root)}
        _atomic_json(path, {
            "schema": SCHEMA,
            "project_id": project.project_id,
            "display_name": project.display_name,
            "aliases": list(project.aliases),
            "field": _field_dict(project.field),
            "target_source": project.target_source,
            "selected_dataset": (
                list(project.selected_dataset)
                if project.selected_dataset else None),
            "data_root": data_link,
            "created_at": project.created_at,
            "last_opened_at": project.last_opened_at,
            "last_science": project.last_science,
            "catalog_match": project.catalog_match,
            "last_dataset_audit": project.last_dataset_audit,
            "target_radius_arcsec": project.target_radius_arcsec,
            "sky_inner_arcsec": project.sky_inner_arcsec,
            "sky_outer_arcsec": project.sky_outer_arcsec,
            "status": project.status,
            "science_status": "DIAGNOSTIC_NOT_ADOPTED",
        })
        project.manifest_sha256 = hashlib.sha256(path.read_bytes()).hexdigest()
        project.manifest_observed = True
    return path

def read_object_project(root: str | Path) -> ObjectProject:
    folder = Path(root).expanduser().resolve()
    path = project_path(folder)
    if path.is_symlink():
        raise ValueError("Unsafe symlinked object project manifest")
    info, manifest_sha256 = _read_manifest_with_sha(path)
    if info.get("schema") != SCHEMA:
        raise ValueError(f"Unsupported project schema in {path}")
    link = info.get("data_root")
    if link is None:
        data_root = None
    elif isinstance(link, dict) and isinstance(link.get("relative"), str):
        raw_relative = link["relative"]
        value = Path(raw_relative)
        # save_object_project() legitimately writes {"relative": "."}
        # when the original FITS inventory is rooted at the project folder.
        # Path(".").parts is empty, unlike an unsafe empty string.
        if (not raw_relative or value.is_absolute()
                or ".." in value.parts
                or (not value.parts and raw_relative != ".")):
            raise ValueError("Unsafe relative data directory in project manifest")
        data_root = folder / value
        # A lexical relative path may traverse via a symlink/junction
        # even without '..'. This format is for *portable in-project*
        # data only. External data remains supported by the existing
        # explicit absolute path field; never silently change provenance.
        try:
            inside = data_root.resolve().is_relative_to(folder)
        except (OSError, RuntimeError) as exc:
            raise ValueError(
                "Unsafe relative data directory in project manifest") from exc
        if not inside:
            raise ValueError(
                "Unsafe relative data directory escapes project root")
    elif isinstance(link, dict) and isinstance(link.get("absolute"), str):
        data_root = Path(link["absolute"])
        if not data_root.is_absolute():
            raise ValueError("External data directory must be absolute")
    else:
        raise ValueError("Invalid data location in project manifest")
    selected = info.get("selected_dataset")
    project = ObjectProject(
        project_id=_project_id(info.get("project_id")),
        display_name=info.get("display_name"), root=folder,
        field=_field_from_payload(info.get("field")),
        aliases=list(info.get("aliases") or []),
        selected_dataset=tuple(selected) if selected else None,
        data_root=data_root, target_source=str(info.get("target_source", "user")),
        created_at=str(info.get("created_at", "")),
        last_opened_at=str(info.get("last_opened_at", "")),
        last_science=info.get("last_science") or {},
        catalog_match=info.get("catalog_match"),
        last_dataset_audit=info.get("last_dataset_audit") or {},
        target_radius_arcsec=float(info.get("target_radius_arcsec", 60.)),
        sky_inner_arcsec=float(info.get("sky_inner_arcsec", 90.)),
        sky_outer_arcsec=float(info.get("sky_outer_arcsec", 180.)),
        manifest_sha256=manifest_sha256, manifest_observed=True)
    _validate(project)
    return project


def _registry_payload() -> dict:
    """Read optional dismissal state without changing the existing schema."""
    path = library_index_path()
    if path.is_symlink():
        raise ValueError("Unsafe symlinked project library index")
    if not path.is_file():
        return {"schema": REGISTRY_SCHEMA, "projects": [],
                "dismissed_legacy_roots": []}
    raw = _read_json(path)
    if raw.get("schema") != REGISTRY_SCHEMA or not isinstance(
            raw.get("projects"), list):
        raise ValueError("Unrecognized project registry schema; no overwrite")
    hidden = raw.get("dismissed_legacy_roots", [])
    if not isinstance(hidden, list) or any(
            not isinstance(root, str) for root in hidden):
        raise ValueError("Invalid dismissed legacy project paths; no overwrite")
    return {**raw, "dismissed_legacy_roots": hidden}


def _registry_rows() -> list[dict]:
    return _registry_payload()["projects"]


def forget_recent_project(project_id: str, root: str | Path) -> None:
    """Remove exactly the selected library entry; NEVER unlink project files.

    Older versions can inject a legacy session pointer into recent_projects
    without registering it. Such suggestions are suppressed by a registry
    tombstone instead of deleting their original last-session pointer.
    """
    with _registry_transaction():
        selected_root = str(Path(root).expanduser().resolve())
        if not selected_root:
            raise ValueError("A project root is required")
        legacy_id = "legacy:" + selected_root
        if project_id != legacy_id:
            project_id = _project_id(project_id)
        payload = _registry_payload()
        rows = payload["projects"]
        matches = [
            row for row in rows
            if isinstance(row, dict) and row.get("project_id") == project_id
            and isinstance(row.get("root"), str)
            and str(Path(row["root"]).expanduser().resolve()) == selected_root
        ]
        if project_id != legacy_id and not matches:
            raise ValueError(
                "Selected project entry is no longer in the library; refresh first")
        if project_id == legacy_id:
            # Only allow dismissing the exact synthetic entry produced by
            # the legacy pointer, not arbitrary user-provided legacy: strings.
            old = last_project_directory()
            if old is None:
                raise ValueError("Legacy suggestion no longer exists")
            original = old.parent if old.name == "data" else old
            if str(original.resolve()) != selected_root:
                raise ValueError("Selected legacy suggestion has changed")
        kept = [
            row for row in rows
            if not (isinstance(row, dict)
                    and row.get("project_id") == project_id
                    and isinstance(row.get("root"), str)
                    and str(Path(row["root"]).expanduser().resolve())
                    == selected_root)
        ]
        # A now-missing registered path may also be the old application's
        # last-session pointer. Keep that suggestion from reappearing.
        hidden = list(payload["dismissed_legacy_roots"])
        if selected_root not in hidden:
            hidden.append(selected_root)
        _atomic_json(library_index_path(), {
            **payload, "projects": kept,
            "dismissed_legacy_roots": hidden,
        })


def recent_projects() -> list[dict]:
    """Registry plus a non-mutating legacy last-project migration suggestion."""
    registry = _registry_payload()
    rows = registry["projects"]
    hidden = set(registry["dismissed_legacy_roots"])
    items = [
        {**r, "available": project_path(r["root"]).is_file()}
        for r in rows if isinstance(r, dict) and
        isinstance(r.get("root"), str) and isinstance(r.get("project_id"), str)
    ]
    old = last_project_directory()
    if old is not None:
        original = old.parent if old.name == "data" else old
        old_key = str(original.resolve())
        if (old_key not in hidden and
                not project_path(original).is_file() and
                not any(Path(r["root"]).resolve() == original.resolve()
                        for r in items)):
            items.append({
                "project_id": "legacy:" + old_key,
                "display_name": original.name, "root": old_key,
                "last_opened_at": "", "available": True,
                "legacy": True,
            })
    return sorted(
        items, key=lambda x: (x.get("last_opened_at", ""), x.get("legacy", False)),
        reverse=True)


def register_object_project(project: ObjectProject) -> None:
    """Serialize project-library updates across Desidera processes."""
    with _registry_transaction():
        _register_object_project_unlocked(project)


def _register_object_project_unlocked(project: ObjectProject) -> None:
    """Only call while already holding the same stable library lock."""
    _validate(project)
    # Use ONE consistent registry snapshot for the whole update.
    # The previous double read could combine a stale 'projects' list
    # with newer dismissal metadata, silently discarding entries
    # registered between those reads. This function runs under the
    # same OS-level advisory lock throughout read-modify-atomic-replace.
    payload = _registry_payload()
    rows = payload["projects"]
    canonical = str(project.root.resolve())
    updated = [x for x in rows if x.get("project_id") != project.project_id and
               Path(x.get("root", "")).resolve() != project.root.resolve()]
    updated.append({
        "project_id": project.project_id, "display_name": project.display_name,
        "root": canonical, "last_opened_at": project.last_opened_at})
    # Preserve legacy dismissals when later imports/touches rewrite the index.
    _atomic_json(library_index_path(), {
        **payload, "projects": updated})


def _touch_object_project_unlocked(project: ObjectProject) -> None:
    """Update passport then registry while holding the library lock.

    Lock ordering is ALWAYS library-index -> passport when both are
    required. A crash between these atomic files leaves the immutable
    passport available for explicit re-import with the same project UUID.
    """
    project.last_opened_at = _now()
    save_object_project(project)
    _register_object_project_unlocked(project)


def touch_object_project(project: ObjectProject) -> None:
    """Keep passport and library publication in one cooperative transaction."""
    with _registry_transaction():
        _touch_object_project_unlocked(project)


def create_object_project(name: str, directory: str | Path, *,
                          field: SkyField | None = None) -> ObjectProject:
    with _registry_transaction():
        root = Path(directory).expanduser().resolve()
        if root.exists() and any(root.iterdir()):
            raise FileExistsError(
                "Project folder is not empty; use Import existing project instead")
        if any(item["display_name"].casefold() == name.strip().casefold()
               and item.get("available") for item in recent_projects()):
            raise ValueError("A project with this name is already in the library")
        project = ObjectProject(
            project_id=str(uuid.uuid4()), display_name=name, root=root,
            field=field)
        # Validate before creating a user-visible directory: an invalid
        # name/radius must not leave a confusing empty pseudo-project.
        _validate(project)
        root.mkdir(parents=True, exist_ok=True)
        save_object_project(project)
        _register_object_project_unlocked(project)
        return project


def import_object_project(directory: str | Path, *,
                          fallback_field: SkyField | None = None,
                          display_name: str | None = None) -> ObjectProject:
    """Adopt or recover one legacy folder without changing original FITS.

    Legacy FITS discovery may require directory traversal, so it occurs
    OUTSIDE the user-library lock. Re-check the passport INSIDE the lock
    before generating any new UUID; a concurrent importer must reuse
    the winner's actual immutable project identity.
    """
    root = Path(directory).expanduser().resolve()
    with _registry_transaction():
        if project_path(root).is_file():
            project = read_object_project(root)
            _touch_object_project_unlocked(project)
            return project

    prior = open_local_project(root, fallback_field=fallback_field)
    with _registry_transaction():
        # Both a parent folder and its data/ child can refer to the same
        # object. A passport may have appeared during the FITS inventory.
        # In particular this repairs a crash between passport publication
        # and project-library registration without ever assigning a new UUID.
        if project_path(prior.root).is_file():
            project = read_object_project(prior.root)
            _touch_object_project_unlocked(project)
            return project

        # Legacy imports use the folder name only as a descriptive default,
        # not an official or scientifically confirmed galaxy identity.
        project = ObjectProject(
            project_id=str(uuid.uuid4()),
            display_name=display_name or prior.root.name, root=prior.root,
            data_root=prior.data_root, field=prior.field,
            target_source=prior.field_source, selected_dataset=prior.selected,
            last_science=prior.last_science,
            target_radius_arcsec=prior.target_radius_arcsec,
            sky_inner_arcsec=prior.sky_inner_arcsec,
            sky_outer_arcsec=prior.sky_outer_arcsec)
        save_object_project(project)
        _register_object_project_unlocked(project)
        return project


def relink_object_project(project_id: str, new_directory: str | Path) -> ObjectProject:
    """Change only the library location after verifying the project's UUID."""
    path = Path(new_directory).expanduser().resolve()
    candidate = read_object_project(path)
    if candidate.project_id != _project_id(project_id):
        raise ValueError("Relink refused: selected directory is another project")
    register_object_project(candidate)
    return candidate
