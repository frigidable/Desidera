"""Offline, source-neutral catalogue identity read model for Galaxy Projects.

Identity assertions are separate from the immutable object-project UUID and
accepted ICRS target. Only explicit user-confirmed catalogue associations are
shown as identities; unconfirmed cone-search neighbours are never aliases.

SIMBAD and NED have separate verified evidence receipts. Neither provider's
coordinates, aliases or measurements automatically overwrite the project.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
import math
from pathlib import Path

from .galaxy_knowledge import (
    _digest, latest_simbad_confirmation, latest_simbad_snapshot,
    read_simbad_snapshot,
)
from .object_projects import ObjectProject, read_object_project
from .ned_knowledge import latest_ned_confirmation, latest_ned_snapshot, read_ned_snapshot

NO_CONFIRMATION = "NO_CONFIRMED_IDENTITY"
CANDIDATES_ONLY = "UNVERIFIED_CANDIDATES_ONLY"
USER_CONFIRMED = "USER_CONFIRMED_CATALOGUE_IDENTITY"


@dataclass(frozen=True)
class CatalogueIdentity:
    provider: str
    main_id: str
    ra_deg: float
    dec_deg: float
    project_offset_arcsec: float | None
    confirmed_at_utc: str
    confirmation_sha256: str
    query_snapshot_sha256: str
    coordinate_bibcode: str
    matched_id: str
    endpoint_host: str | None
    status: str = USER_CONFIRMED
    object_type: str = ""
    redshift: float | None = None
    redshift_refcode: str = ""


@dataclass(frozen=True)
class GalaxyKnowledgeSummary:
    project_id: str
    status: str
    identities: tuple[CatalogueIdentity, ...] = ()
    cached_candidate_count: int = 0
    ned_status: str = "NOT_QUERIED"
    provider_offset_arcsec: float | None = None

    @property
    def confirmed_ids(self) -> tuple[str, ...]:
        return tuple(value.main_id for value in self.identities)


def _distance_arcsec(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
    """Stable great-circle distance for close matches and RA wraparound."""
    a, b = math.radians(dec1), math.radians(dec2)
    dra = math.radians((ra2 - ra1 + 180) % 360 - 180)
    h = math.sin((b - a) / 2) ** 2 + (
        math.cos(a) * math.cos(b) * math.sin(dra / 2) ** 2
    )
    return math.degrees(2 * math.asin(math.sqrt(min(1., max(0., h))))) * 3600


def _include_ned(project: ObjectProject, summary: GalaxyKnowledgeSummary):
    """Append independently confirmed NED identity without reconciling claims."""
    query = latest_ned_snapshot(project.root)
    receipt = latest_ned_confirmation(project.root)
    if receipt is None:
        return GalaxyKnowledgeSummary(
            project_id=summary.project_id, status=summary.status,
            identities=summary.identities,
            cached_candidate_count=summary.cached_candidate_count,
            ned_status=query[1]["status"] if query else "NOT_QUERIED",
        )

    original = read_ned_snapshot(project.root, receipt["snapshot_sha256"])
    candidate = receipt["candidate"]
    if original.get("candidate") != candidate:
        raise ValueError("NED confirmation differs from its query evidence")
    name = candidate.get("main_id")
    if not isinstance(name, str) or not name.strip() or len(name) > 256:
        raise ValueError("Invalid confirmed NED name")
    ra, dec = float(candidate["ra_deg"]), float(candidate["dec_deg"])
    if not (math.isfinite(ra) and math.isfinite(dec) and
            0 <= ra < 360 and -90 <= dec <= 90):
        raise ValueError("Invalid confirmed NED ICRS position")
    redshift = candidate.get("redshift") or {}
    z = redshift.get("value")
    if z is not None:
        z = float(z)
        if not math.isfinite(z):
            raise ValueError("NED redshift is non-finite")
    project_offset = (
        _distance_arcsec(project.field.ra_deg, project.field.dec_deg, ra, dec)
        if project.field else None)
    identity = CatalogueIdentity(
        provider="NED", main_id=name, ra_deg=ra, dec_deg=dec,
        project_offset_arcsec=project_offset,
        confirmed_at_utc=receipt["confirmed_at_utc"],
        confirmation_sha256=_digest(receipt),
        query_snapshot_sha256=receipt["snapshot_sha256"],
        coordinate_bibcode=candidate.get("position_refcode") or "",
        matched_id="", endpoint_host="ned.ipac.caltech.edu",
        object_type=candidate.get("object_type") or "",
        redshift=z, redshift_refcode=redshift.get("refcode") or "",
    )
    separation = None
    if summary.identities:
        simbad = next((x for x in summary.identities
                       if x.provider == "SIMBAD"), None)
        if simbad:
            separation = _distance_arcsec(
                simbad.ra_deg, simbad.dec_deg, ra, dec)
    return GalaxyKnowledgeSummary(
        project_id=summary.project_id, status=USER_CONFIRMED,
        identities=summary.identities + (identity,),
        cached_candidate_count=summary.cached_candidate_count,
        ned_status="USER_CONFIRMED_CATALOGUE_IDENTITY",
        provider_offset_arcsec=separation,
    )


def load_galaxy_knowledge(
    root: str | Path, *, project: ObjectProject | None = None,
) -> GalaxyKnowledgeSummary:
    """Read-only, offline, fail-closed projection from immutable evidence.

    Raises ValueError when a saved receipt is corrupt or incorrectly bound to
    the project. Callers must surface that error, not silently downgrade it to
    'no confirmed identity'. Never downloads, rehashes FITS, or writes files.
    """
    if project is None:
        project = read_object_project(root)
    if Path(root).resolve() != project.root.resolve():
        raise ValueError("Knowledge path differs from the requested galaxy project")
    latest_query = latest_simbad_snapshot(project.root)
    count = len(latest_query[1]["candidates"]) if latest_query else 0
    receipt = latest_simbad_confirmation(project.root)
    if receipt is None:
        return _include_ned(project, GalaxyKnowledgeSummary(
            project_id=project.project_id,
            status=CANDIDATES_ONLY if latest_query else NO_CONFIRMATION,
            cached_candidate_count=count,
        ))

    if (receipt.get("status") != USER_CONFIRMED or
            receipt.get("coordinate_adoption") is not False or
            receipt.get("science_status") != "NOT_SCIENTIFICALLY_ADOPTED"):
        raise ValueError("SIMBAD confirmation has an invalid scientific status")
    name = receipt["candidate"]["main_id"]
    if not isinstance(name, str) or not name.strip() or len(name) > 256:
        raise ValueError("Invalid SIMBAD confirmed identifier")
    try:
        ra = float(receipt["candidate"]["ra_deg"])
        dec = float(receipt["candidate"]["dec_deg"])
    except (ValueError, TypeError, KeyError) as exc:
        raise ValueError("Invalid SIMBAD confirmed coordinate") from exc
    if not (math.isfinite(ra) and math.isfinite(dec) and
            0 <= ra < 360 and -90 <= dec <= 90):
        raise ValueError("Invalid SIMBAD confirmed ICRS position")
    original_snapshot = read_simbad_snapshot(
        project.root, receipt["snapshot_sha256"])
    if receipt["candidate"] not in original_snapshot["candidates"]:
        raise ValueError("SIMBAD identity does not match saved query snapshot")
    endpoint = original_snapshot.get("endpoint_evidence") or {}
    if not isinstance(endpoint, dict):
        raise ValueError("Invalid SIMBAD endpoint provenance")
    host = endpoint.get("endpoint_host")
    if host is not None and not isinstance(host, str):
        raise ValueError("Invalid SIMBAD endpoint name")
    timestamp = receipt.get("confirmed_at_utc")
    if not isinstance(timestamp, str):
        raise ValueError("Missing SIMBAD confirmation timestamp")
    try:
        datetime.fromisoformat(timestamp.replace("Z", "+00:00"))
    except ValueError as exc:
        raise ValueError("Invalid SIMBAD confirmation timestamp") from exc
    offset = None
    if project.field is not None:
        offset = _distance_arcsec(
            project.field.ra_deg, project.field.dec_deg, ra, dec)
    identity = CatalogueIdentity(
        provider="SIMBAD", main_id=name,
        ra_deg=ra, dec_deg=dec, project_offset_arcsec=offset,
        confirmed_at_utc=timestamp,
        confirmation_sha256=_digest(receipt),
        query_snapshot_sha256=receipt["snapshot_sha256"],
        coordinate_bibcode=str(receipt["candidate"].get(
            "coordinate_bibcode") or ""),
        matched_id=str(receipt["candidate"].get("matched_id") or ""),
        endpoint_host=host,
    )
    return _include_ned(project, GalaxyKnowledgeSummary(
        project_id=project.project_id, status=USER_CONFIRMED,
        identities=(identity,), cached_candidate_count=count,
    ))
