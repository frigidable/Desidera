"""Scientific execution fingerprints and explicit validity states (P0.3).

Pure Python, no Qt, no network, no astronomy imports. A recorded result
is never implicitly accepted as a measurement and a file's *historical*
SHA-256 is never misrepresented as proof that today's file is unchanged.
"""
from __future__ import annotations

from collections import defaultdict
from hashlib import sha256
from pathlib import Path
from uuid import UUID
import json
import math
import re

CONTRACT_SCHEMA = "lsbpipe.science-run-contract/1"
SCIENTIFIC_ACTIONS = frozenset(("square", "census", "masks", "figure", "compare"))
APERTURE_ACTIONS = frozenset(("census", "masks", "figure"))
VALID_HASH = re.compile(r"^[0-9a-f]{64}$")
RUN_ID = re.compile(r"^[0-9]{8}T[0-9]{12}Z-[0-9a-f]{10}$")


def canonical_digest(obj: object) -> str:
    """Deterministic across dict insertion order and platform path separators."""
    raw = json.dumps(obj, sort_keys=True, separators=(",", ":"),
                     allow_nan=False, ensure_ascii=False)
    return sha256(raw.encode("utf-8")).hexdigest()


def _uuid(value: str | None) -> str | None:
    if value is None:
        return None
    try:
        ident = UUID(value)
    except (ValueError, TypeError, AttributeError) as exc:
        raise ValueError("Invalid project UUID in scientific context") from exc
    if ident.version != 4:
        raise ValueError("Expected a UUIDv4 project identity")
    return str(ident)


def science_parameters(context: dict) -> dict:
    """Only algorithmically relevant inputs, NOT transient GUI file paths.

    The target/sky rings do not affect square coverage or paired reductions;
    changing their UI controls must not falsely mark those products stale.
    """
    action = context["action"]
    if action not in SCIENTIFIC_ACTIONS:
        raise ValueError(f"Unsupported persisted science action: {action}")
    values = {
        "action": action,
        "release": context["release"], "region": context["region"],
        "brick": context["brick"], "band": context["band"],
        "ra_deg": float(context["ra_deg"]),
        "dec_deg": float(context["dec_deg"]),
        "width_arcmin": float(context["width_arcmin"]),
    }
    if action in APERTURE_ACTIONS:
        values.update({
            "target_radius_arcsec": float(context["target_radius_arcsec"]),
            "sky_inner_arcsec": float(context["sky_inner_arcsec"]),
            "sky_outer_arcsec": float(context["sky_outer_arcsec"]),
        })
    if not all(math.isfinite(x) for x in values.values()
               if isinstance(x, float)):
        raise ValueError("Non-finite science parameters")
    return values


def dataset_identities(context: dict) -> list[str]:
    stem = f"{context['release']}/{context['brick']}/{context['band']}"
    if context["action"] == "compare":
        return [f"{context['release']}/north/{context['brick']}/{context['band']}",
                f"{context['release']}/south/{context['brick']}/{context['band']}"]
    return [f"{context['release']}/{context['region']}/{context['brick']}/{context['band']}"]


def input_content_digest(inputs: dict) -> str:
    """Hash content IDs only: moving a folder does not invalidate its data."""
    snapshot = {}
    if not isinstance(inputs, dict) or not inputs:
        raise ValueError("Scientific run has no original-FITS input snapshot")
    for label, quartet in sorted(inputs.items()):
        if not isinstance(quartet, dict) or not quartet:
            raise ValueError("Missing input products")
        snapshot[label] = {}
        for product, value in sorted(quartet.items()):
            digest = value.get("sha256_local")
            if not isinstance(digest, str) or not VALID_HASH.fullmatch(digest):
                raise ValueError("Input content digest must be a local SHA-256")
            snapshot[label][product] = digest
    return canonical_digest(snapshot)


def make_run_contract(*, context: dict, inputs: dict,
                      project_id: str | None = None,
                      dependencies: tuple[str, ...] = (),
                      code_revision: str | None = None) -> dict:
    identity = _uuid(project_id)
    deps = list(dependencies)
    if any(not isinstance(x, str) or not RUN_ID.fullmatch(x) for x in deps):
        raise ValueError("Dependencies must be immutable run IDs")
    if len(deps) != len(set(deps)):
        raise ValueError("Repeated scientific run dependency")
    params = science_parameters(context)
    parameter_hash = canonical_digest(params)
    input_hash = input_content_digest(inputs)
    return {
        "schema": CONTRACT_SCHEMA,
        "project_id": identity,
        "dataset_ids": dataset_identities(context),
        "parameters": params,
        "parameter_sha256": parameter_hash,
        "input_set_sha256": input_hash,
        "dependencies": deps,
        "code_revision": code_revision,
        "signature_sha256": canonical_digest({
            "project_id": identity, "parameters": parameter_hash,
            "inputs": input_hash, "code_revision": code_revision,
            "dependencies": deps,
        }),
        "status": "DIAGNOSTIC_NOT_ADOPTED",
    }


def check_run_receipt(manifest_path: str | Path, *,
                      project_id: str | None = None,
                      current_context: dict | None = None,
                      current_inputs: dict | None = None,
                      verify_artifacts: bool = True) -> dict:
    """Read-only check; never silently rehash all original FITS by default.

    CURRENT_INPUTS_AND_PARAMETERS requires the caller to pass a separately
    acquired current input snapshot, not merely the historical manifest.
    """
    path = Path(manifest_path)
    result = {"manifest_path": str(path), "state": "INVALID_RECEIPT"}
    if path.is_symlink() or not path.is_file():
        result["reason"] = "Receipt is absent or a symlink"
        return result
    try:
        record = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(record, dict) or record.get("schema") != (
                "lsbpipe.gui-science-qa/1"):
            raise ValueError("Unexpected QA receipt schema")
        result["run_id"] = record["run_id"]
        if record["science_status"] != "DIAGNOSTIC_NOT_ADOPTED":
            raise ValueError("Unexpected scientific adoption status")
        artifacts = record.get("artifacts", {})
        if not isinstance(artifacts, dict) or not artifacts:
            raise ValueError("No completed output artifacts")
        if verify_artifacts:
            for name, meta in artifacts.items():
                if (not isinstance(name, str) or
                        Path(name).name != name or name in ("", ".", "..")):
                    raise ValueError("Unsafe artifact filename")
                artifact = path.parent / name
                if artifact.is_symlink() or not artifact.is_file():
                    result.update(state="ARTIFACT_MISSING", reason=name)
                    return result
                digest = sha256()
                with artifact.open("rb") as source:
                    for block in iter(lambda: source.read(1024 * 1024), b""):
                        digest.update(block)
                if (digest.hexdigest() != meta["sha256_local"] or
                        artifact.stat().st_size != meta["bytes"]):
                    result.update(state="ARTIFACT_MODIFIED", reason=name)
                    return result
        contract = record.get("run_contract")
        if not contract:
            result.update(state="LEGACY_UNVERIFIED",
                          reason="Completed before P0.3 fingerprint contract")
            return result
        if contract.get("schema") != CONTRACT_SCHEMA:
            raise ValueError("Unexpected run contract schema")
        rebuilt = make_run_contract(
            context=record["context"],
            inputs=record["input_original_fits"],
            project_id=contract["project_id"],
            dependencies=tuple(contract["dependencies"]),
            code_revision=contract.get("code_revision"))
        if rebuilt != contract:
            raise ValueError("Scientific run contract fingerprint altered")
        expected = _uuid(project_id)
        if expected is not None and expected != contract["project_id"]:
            result.update(state="PROJECT_MISMATCH",
                          reason="Different immutable object project")
            return result
        if current_context is not None:
            if canonical_digest(science_parameters(current_context)) != (
                    contract["parameter_sha256"]):
                result.update(state="STALE_PARAMETERS",
                              reason="Algorithmic inputs differ")
                return result
        if current_inputs is not None:
            if input_content_digest(current_inputs) != (
                    contract["input_set_sha256"]):
                result.update(state="STALE_INPUTS",
                              reason="Current FITS content differs")
                return result
        if current_context is not None and current_inputs is not None:
            result.update(state="CURRENT_INPUTS_AND_PARAMETERS",
                          reason="Explicit parameter and content fingerprints match")
        else:
            result.update(state="RECORDED_NOT_REVALIDATED",
                          reason="Present-day FITS/parameters not both verified")
        return result
    except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError) as exc:
        result["reason"] = str(exc)
        return result


def discover_run_history(project_root: str | Path, *,
                         data_root: str | Path | None = None,
                         project_id: str | None = None) -> list[dict]:
    """Read completed receipts from project and optional older QA directory.

    Stale .partial directories are ignored; no automatic deletion or replay.
    """
    root = Path(project_root).expanduser().resolve()
    locations = [root / "qa-runs"]
    if data_root is not None:
        parent = Path(data_root).expanduser().resolve().parent / "qa-runs"
        if parent != locations[0]:
            locations.append(parent)
    found = []
    visited = set()
    for location in locations:
        if location.is_symlink() or not location.is_dir():
            continue
        for directory in location.iterdir():
            if not directory.is_dir() or directory.is_symlink() or (
                    directory.name.startswith(".")):
                continue
            manifest = directory / "manifest.json"
            if not manifest.is_file():
                continue
            check = check_run_receipt(manifest, project_id=project_id,
                                      verify_artifacts=True)
            if check.get("state") == "INVALID_RECEIPT":
                continue
            unique = (str(manifest.resolve()), check.get("run_id"))
            if unique in visited:
                continue
            visited.add(unique)
            found.append(check)
    return sorted(found, key=lambda r: r.get("run_id", ""), reverse=True)


def dependency_order(graph: dict[str, list[str]]) -> list[str]:
    """Topological order; unknown dependencies and cycles fail closed."""
    visited, active, ordered = set(), set(), []
    def visit(key):
        if key in active:
            raise ValueError("Scientific run dependency cycle")
        if key in visited:
            return
        if key not in graph:
            raise ValueError(f"Missing scientific dependency: {key}")
        active.add(key)
        for parent in graph[key]:
            visit(parent)
        active.remove(key)
        visited.add(key)
        ordered.append(key)
    for key in sorted(graph):
        visit(key)
    return ordered
