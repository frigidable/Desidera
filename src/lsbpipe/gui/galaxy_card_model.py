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
from lsbpipe.surveys.object_projects import read_object_project


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

    def matches(self, search: str) -> bool:
        query = search.casefold().strip()
        return not query or any(
            query in value.casefold()
            for value in (
                self.display_name, str(self.root), self.project_id,
                self.coordinates, self.dataset, *self.aliases,
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
        return GalaxyCardSnapshot(
            **{**common, "display_name": str(project.display_name)},
            coordinates=_coordinates(project.field),
            dataset=dataset,
            fits_status=fits_status,
            aliases=tuple(project.aliases),
            thumbnail=_thumbnail(root),
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
