"""Read-only snapshots for the object-centric Galaxy Projects card gallery.

This adapter never changes project manifests, touches the last-opened timestamp,
scans FITS directories, starts science jobs, or queries remote catalogues.
A saved data path is not evidence that its FITS products passed integrity/QA.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable

from lsbpipe.surveys.local_project import open_local_project
from lsbpipe.surveys.galaxy_passport_readiness import (
    PassportReadiness, inspect_passport_readiness,
)
from lsbpipe.surveys.galaxy_passport_audit import (
    SavedAuditEvidence, inspect_saved_dataset_audit,
)
from .card_previews import inspect_cached_preview
from .remote_card_previews import inspect_cached_remote_preview
from lsbpipe.surveys.object_projects import read_object_project
from lsbpipe.surveys.knowledge_summary import (
    CatalogueIdentity, load_galaxy_knowledge, NO_CONFIRMATION,
)


@dataclass(frozen=True)
class GalaxyCardSnapshot:
    project_id: str
    display_name: str
    root: Path
    available: bool
    legacy: bool
    coordinates: str
    dataset: str
    fits_status: str
    last_opened: str
    note: str = ""
    aliases: tuple[str, ...] = ()
    thumbnail: Path | None = None
    preview_status: str = "NOT_GENERATED"
    preview_source: str = "NONE"
    local_thumbnail: Path | None = None
    remote_preview_status: str = "REMOTE_RGB_NOT_CACHED"
    knowledge_status: str = NO_CONFIRMATION
    catalogue_identity: CatalogueIdentity | None = None
    catalogue_identities: tuple[CatalogueIdentity, ...] = ()
    provider_offset_arcsec: float | None = None
    ned_status: str = "NOT_QUERIED"
    knowledge_error: str = ""
    readiness: PassportReadiness | None = None
    audit: SavedAuditEvidence | None = None

    def matches(self, search: str) -> bool:
        query = search.casefold().strip()
        return not query or any(
            query in value.casefold()
            for value in (
                self.display_name, str(self.root), self.project_id,
                self.coordinates, self.dataset, *self.aliases,
                *(value.main_id for value in self.catalogue_identities),
                self.catalogue_identity.main_id
                if self.catalogue_identity is not None else "",
            )
        )


def _coordinates(field: object) -> str:
    if field is None:
        return "Coordinates not specified"
    return (
        f"RA {field.ra_deg:.6f}°  ·  "
        f"Dec {field.dec_deg:+.6f}° (ICRS)"
    )


def _thumbnail(root: Path) -> Path | None:
    """Only a dedicated, explicitly generated preview; never decode FITS here."""
    candidate = root / ".desidera" / "card-preview.png"
    return candidate if candidate.is_file() else None


def snapshot_from_record(
    record: dict,
    *,
    read_project: Callable = read_object_project,
    read_legacy: Callable = open_local_project,
) -> GalaxyCardSnapshot:
    """Build a defensively labelled card from the existing recent-project row."""
    root = Path(record["root"])
    available = bool(record.get("available", False))
    legacy = bool(record.get("legacy", False))
    common = dict(
        project_id=str(record["project_id"]),
        display_name=str(record.get("display_name") or root.name or "Unnamed"),
        root=root,
        available=available,
        legacy=legacy,
        last_opened=str(record.get("last_opened_at") or ""),
    )
    if not available:
        return GalaxyCardSnapshot(
            **common, coordinates="Coordinates unavailable",
            dataset="Not inspected", fits_status="Location missing",
            note="Relink the folder in Galaxy Projects.",
        )
    try:
        if legacy:
            project = read_legacy(root)
            return GalaxyCardSnapshot(
                **common, coordinates=_coordinates(project.field),
                dataset="Older local project", fits_status="Not verified",
                note="Import to create an object project.",
                thumbnail=_thumbnail(root),
            )
        project = read_project(root)
        readiness = inspect_passport_readiness(project)
        audit = inspect_saved_dataset_audit(project)
        dataset = (
            " / ".join(str(part) for part in project.selected_dataset)
            if project.selected_dataset else "Not selected"
        )
        if project.data_root is None:
            fits_status = "No local FITS link"
        elif not project.data_root.is_dir():
            fits_status = "Saved FITS path unavailable"
        else:
            fits_status = "FITS linked · integrity not checked"
        # Knowledge is an independent read-only projection, not project
        # metadata. A corrupt receipt is shown as an evidence error and must
        # never be reinterpreted as an accepted SIMBAD association.
        knowledge_status = NO_CONFIRMATION
        catalogue_identity = None
        catalogue_identities = ()
        provider_offset_arcsec = None
        ned_status = "NOT_QUERIED"
        knowledge_error = ""
        try:
            knowledge = load_galaxy_knowledge(project.root, project=project)
            knowledge_status = knowledge.status
            catalogue_identities = knowledge.identities
            catalogue_identity = (
                knowledge.identities[0] if knowledge.identities else None)
            provider_offset_arcsec = knowledge.provider_offset_arcsec
            ned_status = knowledge.ned_status
        except (OSError, ValueError, KeyError, TypeError, AttributeError) as exc:
            knowledge_status = "EVIDENCE_ERROR"
            knowledge_error = f"{type(exc).__name__}: {exc}"
        try:
            thumbnail, preview_state = inspect_cached_preview(project)
        except (OSError, ValueError, KeyError, AttributeError, TypeError):
            thumbnail, preview_state = _thumbnail(root), "LEGACY_PREVIEW_UNVERIFIED"
        try:
            color, remote_state = inspect_cached_remote_preview(project)
        except (OSError, ValueError, KeyError, AttributeError, TypeError):
            color, remote_state = None, "INVALID_REMOTE_RGB_CACHE"
        # RGB is a separate decorative representation. If unavailable,
        # retain the existing local FITS thumbnail without network access.
        if color is not None:
            display_thumbnail = color
            display_status = remote_state
            display_source = "REMOTE_LEGACY_RGB"
        else:
            display_thumbnail = thumbnail
            display_status = preview_state
            display_source = "LOCAL_FITS_GRAYSCALE" if thumbnail else "NONE"
        return GalaxyCardSnapshot(
            **{**common, "display_name": str(project.display_name)},
            coordinates=_coordinates(project.field),
            dataset=dataset,
            fits_status=fits_status,
            aliases=tuple(project.aliases),
            thumbnail=display_thumbnail,
            preview_status=display_status,
            preview_source=display_source,
            local_thumbnail=thumbnail,
            remote_preview_status=remote_state,
            knowledge_status=knowledge_status,
            catalogue_identity=catalogue_identity,
            catalogue_identities=catalogue_identities,
            provider_offset_arcsec=provider_offset_arcsec,
            ned_status=ned_status,
            knowledge_error=knowledge_error,
            readiness=readiness,
            audit=audit,
        )
    except (OSError, ValueError, KeyError, TypeError, AttributeError) as exc:
        return GalaxyCardSnapshot(
            **common, coordinates="Metadata unavailable",
            dataset="Not inspected", fits_status="Not verified",
            note=f"Project metadata needs attention: {exc}",
        )


def snapshots_from_records(
    rows: Iterable[dict],
    *,
    read_project: Callable = read_object_project,
    read_legacy: Callable = open_local_project,
) -> list[GalaxyCardSnapshot]:
    """Preserve the registry's recency order, including missing/legacy entries."""
    return [
        snapshot_from_record(row, read_project=read_project, read_legacy=read_legacy)
        for row in rows
    ]
