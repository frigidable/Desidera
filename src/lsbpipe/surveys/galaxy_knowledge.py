"""P0.7a: conservative SIMBAD identity discovery and immutable local evidence.

Remote candidates are never accepted as authoritative, and queries never
change galaxy coordinates, aliases, original FITS or science measurements.
"""
from __future__ import annotations

from datetime import datetime, timezone
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import re

import numpy as np

from .object_projects import read_object_project

SNAPSHOT_SCHEMA = "desidera.knowledge.simbad-query/1"
CONFIRMATION_SCHEMA = "desidera.knowledge.simbad-confirmation/1"
SHA = re.compile(r"^[0-9a-f]{64}$")


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _canonical(payload: dict) -> bytes:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False, allow_nan=False).encode("utf-8")


def _digest(payload: dict) -> str:
    return sha256(_canonical(payload)).hexdigest()


def _text(value: object) -> str:
    if value is None or np.ma.is_masked(value):
        return ""
    if isinstance(value, bytes):
        value = value.decode("utf-8", "replace")
    return str(value).strip()


def _optional_text(row, key: str) -> str:
    return _text(row[key]) if key in row.colnames else ""


def _candidate_rows(table, origin: str, centre) -> list[dict]:
    if table is None:
        return []
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    found = []
    for row in table:
        name = _optional_text(row, "main_id")
        if not name or len(name) > 256:
            continue
        try:
            ra, dec = float(row["ra"]), float(row["dec"])
        except (KeyError, ValueError, TypeError):
            continue
        if not (math.isfinite(ra) and math.isfinite(dec) and
                0 <= ra < 360 and -90 <= dec <= 90):
            continue
        separation = (float(centre.separation(SkyCoord(
            ra=ra * u.deg, dec=dec * u.deg, frame="icrs")).arcsec)
            if centre is not None else None)
        found.append({
            "main_id": name, "ra_deg": ra, "dec_deg": dec,
            "separation_arcsec": separation,
            "coordinate_bibcode": _optional_text(row, "coo_bibcode"),
            "matched_id": _optional_text(row, "matched_id"),
            "discovered_via": [origin],
        })
    return found


def query_simbad_candidates(*, name: str, ra_deg: float | None,
                            dec_deg: float | None,
                            radius_arcsec: float = 60.,
                            limit: int = 30, client=None) -> dict:
    """Manual external lookup only; preserve multiple candidates for review.

    Inject a fake client for offline tests. No remote service is contacted when
    reading stored snapshots or opening a project.
    """
    if not isinstance(name, str) or not name.strip() or len(name) > 140:
        raise ValueError("Galaxy search name must be 1–140 characters")
    if (ra_deg is None) != (dec_deg is None):
        raise ValueError("RA and Dec must be provided together")
    if not isinstance(limit, int) or isinstance(limit, bool) or not 1 <= limit <= 50:
        raise ValueError("SIMBAD candidate limit must be 1–50")
    radius_arcsec = float(radius_arcsec)
    if not math.isfinite(radius_arcsec) or not 1 <= radius_arcsec <= 300:
        raise ValueError("SIMBAD search radius must be between 1 and 300 arcsec")
    centre = None
    if ra_deg is not None:
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        ra_deg, dec_deg = float(ra_deg), float(dec_deg)
        if not (math.isfinite(ra_deg) and math.isfinite(dec_deg) and
                0 <= ra_deg < 360 and -90 <= dec_deg <= 90):
            raise ValueError("Invalid project ICRS coordinates")
        centre = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs")
    if client is None:
        try:
            from astroquery.simbad import Simbad
        except ImportError as exc:
            raise RuntimeError(
                'SIMBAD requires optional dependency: pip install -e ".[knowledge]"'
            ) from exc
        client = Simbad()
        client.TIMEOUT = 15
    client.ROW_LIMIT = limit + 1
    # Deliberately preserve alternatives; no closest-object auto-matching.
    by_id: dict[str, dict] = {}
    for origin, table in (
        ("name", client.query_object(name.strip(), wildcard=False)),
        ("cone", client.query_region(centre, radius=radius_arcsec * __import__(
            "astropy.units", fromlist=["arcsec"]).arcsec) if centre is not None else None),
    ):
        for candidate in _candidate_rows(table, origin, centre):
            key = candidate["main_id"].casefold()
            if key not in by_id:
                by_id[key] = candidate
            else:
                previous = by_id[key]
                previous["discovered_via"] = sorted(set(
                    previous["discovered_via"] + candidate["discovered_via"]))
    candidates = sorted(by_id.values(), key=lambda row: (
        "name" not in row["discovered_via"],
        row["separation_arcsec"] if row["separation_arcsec"] is not None
        else float("inf"), row["main_id"].casefold()))
    return {
        "schema": SNAPSHOT_SCHEMA, "source": "SIMBAD",
        "fetched_at_utc": _now(),
        "query": {"name": name.strip(), "ra_deg": ra_deg,
                  "dec_deg": dec_deg, "radius_arcsec": radius_arcsec,
                  "limit": limit},
        "candidates": candidates[:limit],
        "truncated": len(candidates) > limit,
        "status": "CANDIDATES_UNVERIFIED",
        "science_status": "NOT_SCIENTIFICALLY_ADOPTED",
    }


def _folder(root: str | Path) -> tuple[Path, str]:
    project = read_object_project(root)
    folder = project.root / "knowledge" / "simbad"
    if folder.is_symlink():
        raise ValueError("Refusing linked SIMBAD evidence directory")
    return folder, project.project_id


def _save_once(path: Path, payload: dict) -> str:
    """Content-addressed, create-only receipt; existing bytes never replaced."""
    raw = _canonical(payload) + b"\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_symlink():
        raise ValueError("Refusing symlink evidence")
    try:
        with path.open("xb") as output:
            output.write(raw)
            output.flush()
            os.fsync(output.fileno())
    except FileExistsError:
        if path.read_bytes() != raw:
            raise ValueError("Existing evidence content differs from immutable record")
    return path.name


def _read_receipt(path: Path, schema: str) -> dict:
    if path.is_symlink() or not path.is_file():
        raise ValueError("Missing or linked knowledge receipt")
    label = ("snapshot-" if schema == SNAPSHOT_SCHEMA else "confirmation-")
    raw_digest = path.stem.removeprefix(label)
    if not path.stem.startswith(label) or not SHA.fullmatch(raw_digest):
        raise ValueError("Incorrect knowledge receipt filename")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if (not isinstance(payload, dict) or payload.get("schema") != schema or
            _digest(payload) != raw_digest):
        raise ValueError("Corrupt or modified knowledge receipt")
    return payload


def store_simbad_snapshot(root: str | Path, result: dict) -> str:
    folder, project_id = _folder(root)
    if (not isinstance(result, dict) or result.get("schema") != SNAPSHOT_SCHEMA or
            result.get("source") != "SIMBAD" or
            result.get("science_status") != "NOT_SCIENTIFICALLY_ADOPTED" or
            result.get("status") != "CANDIDATES_UNVERIFIED" or
            not isinstance(result.get("candidates"), list) or
            len(result["candidates"]) > 50):
        raise ValueError("Invalid, untrusted SIMBAD query payload")
    payload = {**result, "project_id": project_id}
    digest = _digest(payload)
    _save_once(folder / f"snapshot-{digest}.json", payload)
    return digest


def read_simbad_snapshot(root: str | Path, digest: str) -> dict:
    folder, project_id = _folder(root)
    if not isinstance(digest, str) or not SHA.fullmatch(digest):
        raise ValueError("Invalid snapshot digest")
    result = _read_receipt(folder / f"snapshot-{digest}.json", SNAPSHOT_SCHEMA)
    if result.get("project_id") != project_id:
        raise ValueError("SIMBAD snapshot belongs to another galaxy project")
    return result


def latest_simbad_snapshot(root: str | Path) -> tuple[str, dict] | None:
    folder, _ = _folder(root)
    if not folder.is_dir():
        return None
    paths = list(folder.glob("snapshot-*.json"))
    if not paths:
        return None
    items = [(p.stem.removeprefix("snapshot-"),
              read_simbad_snapshot(root, p.stem.removeprefix("snapshot-")))
             for p in paths]
    return max(items, key=lambda entry: entry[1].get("fetched_at_utc", ""))


def confirm_simbad_identity(root: str | Path, snapshot_digest: str,
                            candidate_main_id: str) -> str:
    """Explicit user confirmation; never alter authoritative project target."""
    snapshot = read_simbad_snapshot(root, snapshot_digest)
    if not isinstance(candidate_main_id, str):
        raise ValueError("Candidate name is required")
    candidates = [c for c in snapshot["candidates"]
                  if c["main_id"] == candidate_main_id]
    if len(candidates) != 1:
        raise ValueError("Selected identity is not unique in snapshot")
    folder, project_id = _folder(root)
    payload = {
        "schema": CONFIRMATION_SCHEMA, "project_id": project_id,
        "snapshot_sha256": snapshot_digest,
        "confirmed_at_utc": _now(),
        "candidate": candidates[0],
        "status": "USER_CONFIRMED_CATALOGUE_IDENTITY",
        "coordinate_adoption": False,
        "science_status": "NOT_SCIENTIFICALLY_ADOPTED",
    }
    digest = _digest(payload)
    _save_once(folder / f"confirmation-{digest}.json", payload)
    return digest


def latest_simbad_confirmation(root: str | Path) -> dict | None:
    folder, project_id = _folder(root)
    if not folder.is_dir():
        return None
    results = []
    for path in folder.glob("confirmation-*.json"):
        receipt = _read_receipt(path, CONFIRMATION_SCHEMA)
        if receipt.get("project_id") != project_id:
            raise ValueError("SIMBAD confirmation belongs to another galaxy")
        snapshot = read_simbad_snapshot(root, receipt["snapshot_sha256"])
        if receipt.get("candidate") not in snapshot["candidates"]:
            raise ValueError("SIMBAD confirmation disagrees with query evidence")
        results.append(receipt)
    return max(results, key=lambda x: x["confirmed_at_utc"]) if results else None
