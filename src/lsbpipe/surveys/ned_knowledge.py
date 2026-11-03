"""NED ObjectLookup adapter: independent name resolution and evidence receipts.

The documented NED /srs/ObjectLookup JSON endpoint is deliberately kept
separate from SIMBAD. No automatic catalogue-to-project association, physical
parameter adoption, background processing or FITS operations.
"""
from __future__ import annotations

from datetime import datetime, timezone
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import re
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from .object_projects import read_object_project

API = "https://ned.ipac.caltech.edu/srs/ObjectLookup"
SCHEMA = "desidera.knowledge.ned-objectlookup/1"
CONFIRM_SCHEMA = "desidera.knowledge.ned-confirmation/1"
SHA = re.compile(r"^[a-f0-9]{64}$")
MAX_RESPONSE_BYTES = 262144
STATUSES = {
    0: "NAME_NOT_INTERPRETED",
    1: "AMBIGUOUS_NAME",
    2: "NAME_NOT_IN_DATABASE",
    3: "KNOWN_OBJECT",
}


def _now():
    return datetime.now(timezone.utc).isoformat()


def _canonical(x):
    return json.dumps(x, ensure_ascii=False, sort_keys=True,
                      separators=(",", ":"), allow_nan=False).encode("utf-8")


def _hash(x):
    return sha256(_canonical(x)).hexdigest()


def _number(x):
    if x is None or isinstance(x, bool):
        return None
    try:
        value = float(x)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def _optional_str(x):
    return str(x).strip() if isinstance(x, str) else ""


def _remote_lookup(name: str) -> dict:
    """One bounded HTTPS POST, opt-in. No undocumented legacy endpoints."""
    payload = {"name": {"v": name}, "aliases": {"v": True}}
    encoded = urlencode({"json": json.dumps(payload, separators=(",", ":"))})
    req = Request(
        API, data=encoded.encode("utf-8"), method="POST",
        headers={"Accept": "application/json",
                 "Content-Type": "application/x-www-form-urlencoded",
                 "User-Agent": "DESIDERA/0.1 scientific catalogue client"},
    )
    with urlopen(req, timeout=15) as response:
        if response.geturl().split("/")[2].lower() != "ned.ipac.caltech.edu":
            raise ValueError("NED redirected outside its authoritative host")
        size = response.read(MAX_RESPONSE_BYTES + 1)
    if len(size) > MAX_RESPONSE_BYTES:
        raise ValueError("NED lookup exceeded the allowed response size")
    return json.loads(size.decode("utf-8"))


def query_ned_identity(name: str, *, transport=None) -> dict:
    """Manual lookup; preserve ambiguity, source metadata, and typed claims."""
    if not isinstance(name, str) or not name.strip() or len(name.strip()) > 140:
        raise ValueError("NED search name must have 1–140 characters")
    source = transport or _remote_lookup
    raw = source(name.strip())
    if not isinstance(raw, dict):
        raise ValueError("NED returned a non-object JSON response")
    version = raw.get("Version")
    if (not isinstance(version, str) or not version.split(".")[0].isdigit()
            or int(version.split(".")[0]) != 2):
        raise ValueError("Unsupported NED ObjectLookup response major version")
    if raw.get("StatusCode") != 100:
        raise ValueError("NED ObjectLookup reported a non-success service status")
    result_code = raw.get("ResultCode")
    if type(result_code) is not int or result_code not in STATUSES:
        raise ValueError("NED ObjectLookup returned unknown result code")

    candidate = None
    aliases = []
    if result_code == 3:
        preferred = raw.get("Preferred")
        if not isinstance(preferred, dict):
            raise ValueError("Known NED object has no Preferred metadata")
        official_name = _optional_str(preferred.get("Name"))
        position = preferred.get("Position")
        if not official_name or not isinstance(position, dict):
            raise ValueError("NED known object lacks name/position metadata")
        ra = _number(position.get("RA"))
        dec = _number(position.get("Dec"))
        if ra is None or dec is None or not (0 <= ra < 360 and -90 <= dec <= 90):
            raise ValueError("NED returned invalid or missing ICRS coordinates")
        typ = preferred.get("ObjType") or {}
        z = preferred.get("Redshift") or {}
        if not isinstance(typ, dict) or not isinstance(z, dict):
            raise ValueError("Malformed NED object-type/redshift metadata")
        z_value = _number(z.get("Value"))
        z_error = _number(z.get("Uncertainty"))
        if z_error is not None and z_error < 0:
            raise ValueError("Negative NED redshift uncertainty")
        candidate = {
            "main_id": official_name[:256], "ra_deg": ra,
            "dec_deg": dec,
            "object_type": _optional_str(typ.get("Value")),
            "position_refcode": _optional_str(position.get("RefCode")),
            "redshift": {
                "value": z_value, "uncertainty": z_error,
                "refcode": _optional_str(z.get("RefCode")),
                "quality_flag": _optional_str(z.get("QualityFlag")),
                "frame": "heliocentric", "adopted": False,
            },
        }
        candidate["main_id"] = official_name
        if len(official_name) > 256:
            raise ValueError("NED identifier too long")
        aliases_raw = raw.get("Aliases") or []
        if not isinstance(aliases_raw, list):
            raise ValueError("Malformed NED aliases")
        aliases = [a.strip() for a in aliases_raw if isinstance(a, str)
                   and a.strip() and len(a) <= 256][:200]
    alternatives = []
    if result_code == 1:
        interpreted = raw.get("Interpreted") or {}
        if not isinstance(interpreted, dict):
            raise ValueError("Malformed ambiguous NED interpretation")
        possible = interpreted.get("Aliases") or []
        if not isinstance(possible, list):
            raise ValueError("Malformed ambiguous NED candidates")
        alternatives = [x.strip() for x in possible
                        if isinstance(x, str) and x.strip() and len(x) <= 256][:100]

    return {
        "schema": SCHEMA, "provider": "NED", "api": API,
        "fetched_at_utc": _now(), "query_name": name.strip(),
        "api_version": version, "result_code": result_code,
        "status": STATUSES[result_code], "candidate": candidate,
        "aliases_unadopted": aliases,
        "ambiguous_names": alternatives,
        "source_response": raw,
        "science_status": "NOT_SCIENTIFICALLY_ADOPTED",
    }


def _folder(root):
    p = read_object_project(root)
    folder = p.root / "knowledge" / "ned"
    if folder.is_symlink():
        raise ValueError("Refusing symlinked NED evidence directory")
    return folder, p.project_id


def _save_new(path, payload):
    raw = _canonical(payload) + b"\n"
    if path.is_symlink():
        raise ValueError("Refusing symlinked NED receipt")
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with path.open("xb") as fp:
            fp.write(raw)
            fp.flush()
            os.fsync(fp.fileno())
    except FileExistsError:
        if path.is_symlink() or path.read_bytes() != raw:
            raise ValueError("Existing immutable NED receipt differs")


def _read_file(path, schema, prefix):
    if path.is_symlink() or not path.is_file():
        raise ValueError("NED evidence receipt missing or symlinked")
    digest = path.stem.removeprefix(prefix)
    if not path.stem.startswith(prefix) or not SHA.fullmatch(digest):
        raise ValueError("Invalid NED evidence filename")
    try:
        record = json.loads(path.read_bytes())
    except (OSError, UnicodeError, ValueError) as exc:
        raise ValueError("Unreadable NED evidence receipt") from exc
    if (not isinstance(record, dict) or record.get("schema") != schema
            or _hash(record) != digest):
        raise ValueError("Corrupted or altered NED evidence receipt")
    return record


def save_ned_snapshot(root, result):
    folder, project_id = _folder(root)
    if (not isinstance(result, dict) or result.get("schema") != SCHEMA
            or result.get("provider") != "NED" or result.get("api") != API
            or result.get("status") not in STATUSES.values()
            or result.get("science_status") != "NOT_SCIENTIFICALLY_ADOPTED"):
        raise ValueError("Invalid NED query result")
    record = {**result, "project_id": project_id}
    digest = _hash(record)
    _save_new(folder / f"snapshot-{digest}.json", record)
    return digest


def read_ned_snapshot(root, digest):
    folder, project_id = _folder(root)
    if not isinstance(digest, str) or not SHA.fullmatch(digest):
        raise ValueError("Invalid NED snapshot digest")
    record = _read_file(folder / f"snapshot-{digest}.json", SCHEMA, "snapshot-")
    if record.get("project_id") != project_id:
        raise ValueError("NED snapshot belongs to another galaxy project")
    if record.get("science_status") != "NOT_SCIENTIFICALLY_ADOPTED":
        raise ValueError("NED snapshot scientific status altered")
    return record


def latest_ned_snapshot(root):
    folder, _ = _folder(root)
    if not folder.is_dir():
        return None
    paths = list(folder.glob("snapshot-*.json"))
    snapshots = [(p.stem.removeprefix("snapshot-"),
                  read_ned_snapshot(root, p.stem.removeprefix("snapshot-")))
                 for p in paths]
    return max(snapshots, key=lambda x: x[1]["fetched_at_utc"]) if snapshots else None


def confirm_ned_identity(root, snapshot_digest):
    snapshot = read_ned_snapshot(root, snapshot_digest)
    if snapshot.get("status") != "KNOWN_OBJECT" or not snapshot.get("candidate"):
        raise ValueError("NED identity may be confirmed only for known named objects")
    folder, project_id = _folder(root)
    payload = {
        "schema": CONFIRM_SCHEMA, "project_id": project_id,
        "snapshot_sha256": snapshot_digest,
        "confirmed_at_utc": _now(),
        "candidate": snapshot["candidate"],
        "status": "USER_CONFIRMED_CATALOGUE_IDENTITY",
        "coordinate_adoption": False,
        "science_status": "NOT_SCIENTIFICALLY_ADOPTED",
    }
    digest = _hash(payload)
    _save_new(folder / f"confirmation-{digest}.json", payload)
    return digest


def latest_ned_confirmation(root):
    folder, project_id = _folder(root)
    if not folder.is_dir():
        return None
    confirmations = []
    for p in folder.glob("confirmation-*.json"):
        receipt = _read_file(p, CONFIRM_SCHEMA, "confirmation-")
        if receipt.get("project_id") != project_id:
            raise ValueError("NED confirmation belongs to another galaxy project")
        original = read_ned_snapshot(root, receipt["snapshot_sha256"])
        if (original["status"] != "KNOWN_OBJECT"
                or receipt.get("candidate") != original.get("candidate")
                or receipt.get("status") != "USER_CONFIRMED_CATALOGUE_IDENTITY"
                or receipt.get("coordinate_adoption") is not False
                or receipt.get("science_status") != "NOT_SCIENTIFICALLY_ADOPTED"):
            raise ValueError("NED confirmation disagrees with source evidence")
        confirmations.append(receipt)
    return max(confirmations, key=lambda x: x["confirmed_at_utc"]) if confirmations else None
