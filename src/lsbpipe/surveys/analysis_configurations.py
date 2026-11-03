"""P0.8a: immutable, content-verified analysis SETUP snapshots.

A configuration captures project parameters and the *recorded* dataset
fingerprints. It neither rehashes FITS nor certifies scientific readiness.
"""
from __future__ import annotations

from datetime import datetime, timezone
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import re

from .object_projects import read_object_project
from .dataset_registry import read_dataset_index, CHECKSUM_VERIFIED

SCHEMA = "desidera.analysis-configuration/1"
SHA = re.compile(r"^[0-9a-f]{64}$")


def _canonical(value: dict) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False, allow_nan=False).encode("utf-8")


def _hash(value: dict) -> str:
    return sha256(_canonical(value)).hexdigest()


def _folder(root: str | Path) -> tuple[Path, str]:
    p = read_object_project(root)
    directory = p.root / "analyses" / "configurations"
    if directory.is_symlink():
        raise ValueError("Refusing symlinked analysis configuration directory")
    return directory, p.project_id


def capture_analysis_setup(root: str | Path) -> dict:
    """Get a deterministic, read-only snapshot of currently recorded choices."""
    p = read_object_project(root)
    registry = read_dataset_index(p.root)
    selected = tuple(p.selected_dataset) if p.selected_dataset else None
    source_by_id = {s["source_id"]: s for s in registry["sources"]}
    selected_assets = []
    if selected is not None:
        release, region, brick, band = selected
        for a in registry["assets"]:
            s = source_by_id[a["source_id"]]
            if (s.get("release"), s.get("region"), a.get("brick")) != (
                    release, region, brick):
                continue
            if a.get("band") not in (band, None):
                continue
            selected_assets.append({
                "asset_id": a["asset_id"], "source_id": a["source_id"],
                "provider": s.get("provider"), "release": s.get("release"),
                "region": s.get("region"), "brick": a.get("brick"),
                "band": a.get("band"), "product_kind": a.get("product_kind"),
                "product_family": a.get("product_family"),
                "sha256_recorded": a.get("sha256"),
                "checksum_status_at_freeze": a.get("checksum_status"),
                "size_bytes_recorded": a.get("size_bytes"),
                "fits_status_at_freeze": a.get("fits_status"),
            })
    selected_assets.sort(key=lambda a: (str(a["source_id"]), str(a["asset_id"])))
    coord = (None if p.field is None else {
        "ra_deg": p.field.ra_deg,
        "dec_deg": p.field.dec_deg,
        "width_arcmin": p.field.width_arcmin,
    })
    complete_hashes = bool(selected_assets) and all(
        a["checksum_status_at_freeze"] == CHECKSUM_VERIFIED and
        isinstance(a["sha256_recorded"], str) and
        SHA.fullmatch(a["sha256_recorded"])
        for a in selected_assets)
    if selected is None:
        evidence = "SETUP_ONLY_NO_DATASET"
    elif not selected_assets:
        evidence = "NO_REGISTERED_ORIGINAL_INPUTS"
    elif complete_hashes:
        evidence = "RECORDED_SHA256_NOT_REVALIDATED"
    else:
        evidence = "PARTIAL_OR_UNVERIFIED_INPUT_EVIDENCE"
    return {
        "target_icrs": coord,
        "selected_dataset": list(selected) if selected else None,
        "parameters": {
            "target_radius_arcsec": p.target_radius_arcsec,
            "sky_inner_arcsec": p.sky_inner_arcsec,
            "sky_outer_arcsec": p.sky_outer_arcsec,
        },
        "selected_registered_assets": selected_assets,
        "input_evidence_state": evidence,
    }


def freeze_analysis_configuration(root: str | Path, label: str) -> tuple[str, dict]:
    if (not isinstance(label, str) or not label.strip() or
            len(label.strip()) > 100 or any(ord(x) < 32 for x in label)):
        raise ValueError("Configuration label must contain 1–100 printable characters")
    directory, project_id = _folder(root)
    snapshot = capture_analysis_setup(root)
    payload = {
        "schema": SCHEMA, "project_id": project_id,
        "label": label.strip(),
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "setup": snapshot, "setup_sha256": _hash(snapshot),
        "code_revision": None,
        "status": "FROZEN_NOT_RUN",
        "science_status": "NOT_SCIENTIFICALLY_ADOPTED",
        "input_fits_rehashed_during_freeze": False,
    }
    digest = _hash(payload)
    path = directory / f"config-{digest}.json"
    if path.is_symlink():
        raise ValueError("Refusing symlinked analysis configuration")
    directory.mkdir(parents=True, exist_ok=True)
    raw = _canonical(payload) + b"\n"
    # Exclusive create: neither old configurations nor user science files
    # can be replaced. A crash-damaged partial receipt fails later validation.
    try:
        with path.open("xb") as output:
            output.write(raw)
            output.flush()
            os.fsync(output.fileno())
    except FileExistsError:
        if path.read_bytes() != raw:
            raise ValueError("Existing immutable configuration was altered")
    return digest, payload


def read_analysis_configuration(root: str | Path, digest: str) -> dict:
    directory, project_id = _folder(root)
    if not isinstance(digest, str) or not SHA.fullmatch(digest):
        raise ValueError("Invalid configuration identifier")
    path = directory / f"config-{digest}.json"
    if not path.is_file() or path.is_symlink():
        raise ValueError("Configuration is missing or linked")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ValueError(f"Cannot read immutable configuration: {exc}") from exc
    if (not isinstance(payload, dict) or payload.get("schema") != SCHEMA or
            payload.get("project_id") != project_id or
            payload.get("status") != "FROZEN_NOT_RUN" or
            payload.get("science_status") != "NOT_SCIENTIFICALLY_ADOPTED" or
            payload.get("input_fits_rehashed_during_freeze") is not False or
            not isinstance(payload.get("setup"), dict) or
            payload.get("setup_sha256") != _hash(payload["setup"]) or
            _hash(payload) != digest):
        raise ValueError("Immutable configuration is corrupted or belongs to another project")
    return payload


def list_analysis_configurations(root: str | Path) -> list[dict]:
    directory, _ = _folder(root)
    if not directory.is_dir():
        return []
    items = []
    for path in directory.glob("config-*.json"):
        digest = path.stem.removeprefix("config-")
        item = read_analysis_configuration(root, digest)
        items.append({"digest": digest, "label": item["label"],
                      "created_at_utc": item["created_at_utc"],
                      "input_evidence_state":
                          item["setup"]["input_evidence_state"]})
    return sorted(items, key=lambda r: (r["created_at_utc"], r["digest"]),
                  reverse=True)


def compare_current_setup(root: str | Path, digest: str) -> dict:
    old = read_analysis_configuration(root, digest)
    current = capture_analysis_setup(root)
    state = ("SETUP_METADATA_MATCHES_INPUTS_NOT_REHASHED"
             if _hash(current) == old["setup_sha256"] else
             "PROJECT_OR_DATASET_METADATA_CHANGED")
    return {
        "configuration_sha256": digest, "state": state,
        "input_fits_rehashed": False, "science_validity_certified": False,
        "current_setup_sha256": _hash(current),
        "stored_setup_sha256": old["setup_sha256"],
    }
