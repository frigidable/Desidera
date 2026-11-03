"""P0 foundation: multiple immutable observational assets per galaxy project.

Sidecar registry, not a replacement for lsbpipe-object.json or legacy session
files. Discovering/importing is offline and does NOT claim FITS/WCS validity,
exposure independence, photometric suitability, or scientific adoption.
"""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import uuid

from .dataset_integrity import discover_dataset_inventory
from .local_project import _atomic_json
from .object_projects import read_object_project

SCHEMA = "lsbpipe.dataset-index/1"
CONFIG_SCHEMA = "lsbpipe.analysis-configuration/1"
PRODUCT_KINDS = frozenset(("image", "invvar", "nexp", "maskbits"))
CHECKSUM_PENDING = "NOT_HASHED"
CHECKSUM_VERIFIED = "FILE_SHA256_VERIFIED"
CHECKSUM_CHANGED = "CHANGED_SINCE_CHECKSUM"
FITS_UNVERIFIED = "NOT_VALIDATED"
INDEPENDENCE_UNKNOWN = "UNKNOWN"


def _now():
    return datetime.now(timezone.utc).isoformat()


def registry_path(project_root: str | Path) -> Path:
    """New optional sidecar; never overwrite historical project metadata."""
    return Path(project_root).expanduser().resolve() / "datasets" / "index.json"


def _project(root: str | Path):
    return read_object_project(root)


def _path_ref(path: Path, project_root: Path) -> dict:
    path = path.expanduser().absolute()
    if not path.is_absolute():
        raise ValueError("Data asset must have an absolute filesystem path")
    if path.is_symlink():
        raise ValueError("Refusing to index symlinked original FITS")
    project_root = project_root.resolve()
    try:
        relative = path.relative_to(project_root)
        return {"relative": relative.as_posix()}
    except ValueError:
        return {"absolute": str(path)}


def resolve_asset_path(project_root: str | Path, record: dict) -> Path:
    """Resolve indexed file without allowing traversal outside project root."""
    root = Path(project_root).expanduser().resolve()
    ref = record.get("path_ref")
    if not isinstance(ref, dict):
        raise ValueError("Missing data asset path reference")
    if set(ref) == {"relative"}:
        raw = ref["relative"]
        if not isinstance(raw, str) or not raw:
            raise ValueError("Invalid indexed relative path")
        candidate = Path(raw)
        if candidate.is_absolute() or ".." in candidate.parts:
            raise ValueError("Unsafe indexed relative data path")
        path = root / candidate
        if not path.resolve().is_relative_to(root):
            raise ValueError("Indexed file escapes galaxy project directory")
    elif set(ref) == {"absolute"}:
        raw = ref["absolute"]
        if not isinstance(raw, str) or not Path(raw).is_absolute():
            raise ValueError("Invalid absolute external FITS path")
        path = Path(raw)
    else:
        raise ValueError("Unrecognized data asset path reference")
    if path.is_symlink():
        raise ValueError("Refusing symlinked original FITS")
    return path


def _empty(project_id: str) -> dict:
    return {
        "schema": SCHEMA,
        "project_id": project_id,
        "created_at": _now(),
        "sources": [],
        "assets": [],
        "collections": [],
        "analysis_configurations": [],
        "runs": [],
        "scientific_status": "NOT_SCIENTIFICALLY_ADOPTED",
    }


def read_dataset_index(project_root: str | Path) -> dict:
    """Read without creating files; a legacy project with no index is valid."""
    project = _project(project_root)
    path = registry_path(project.root)
    if not path.is_file():
        return _empty(project.project_id)
    try:
        record = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ValueError(f"Unreadable dataset index: {exc}") from exc
    if not isinstance(record, dict) or record.get("schema") != SCHEMA:
        raise ValueError("Unknown dataset index schema; no migration or overwrite")
    if record.get("project_id") != project.project_id:
        raise ValueError("Dataset index refers to a different galaxy project")
    for key in ("sources", "assets", "collections", "analysis_configurations", "runs"):
        if not isinstance(record.get(key), list):
            raise ValueError(f"Invalid dataset index collection: {key}")
    source_ids = [s.get("source_id") for s in record["sources"]]
    asset_ids = [a.get("asset_id") for a in record["assets"]]
    if (None in source_ids or None in asset_ids or
            len(source_ids) != len(set(source_ids)) or
            len(asset_ids) != len(set(asset_ids))):
        raise ValueError("Duplicate or missing immutable source/asset ID")
    if any(a.get("source_id") not in source_ids for a in record["assets"]):
        raise ValueError("Dataset asset references an unknown data source")
    return record


def _save(project_root: str | Path, index: dict) -> Path:
    # Re-validate identity before an atomic sidecar write.
    project = _project(project_root)
    if index.get("schema") != SCHEMA or index.get("project_id") != project.project_id:
        raise ValueError("Cannot write data assets to another galaxy project")
    path = registry_path(project.root)
    if path.exists():
        current = read_dataset_index(project.root)
        if current["project_id"] != index["project_id"]:
            raise ValueError("Existing registry belongs to a different project")
    _atomic_json(path, index)
    return path


def import_legacy_fits(project_root: str | Path, data_root: str | Path) -> dict:
    """Register arbitrarily many on-disk Legacy coadd products, no download.

    The existing P0.2 *discovery* function enumerates every matched science
    product; full FITS/WCS audit and file SHA-256 are separate explicit steps.
    Preserve the historical project selection and every source file in place.
    """
    project = _project(project_root)
    root = Path(data_root).expanduser().resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"Original data directory not found: {root}")
    # Legacy projects store original downloads either at the object root
    # or under its data/ child. Accept either folder selection without moving
    # files, while storing the exact resolved original asset location.
    if not any((root / release).is_dir()
               for release in ("dr9", "dr10", "dr11")) and (
                   root / "data").is_dir():
        root = (root / "data").resolve()
    inventory = discover_dataset_inventory(root)
    registry = read_dataset_index(project.root)
    source_by_key = {s.get("import_key"): s for s in registry["sources"]}
    asset_by_key = {a.get("import_key"): a for a in registry["assets"]}
    sources_added = assets_added = assets_changed = 0
    for entry in inventory:
        release, region, brick = (entry["release"], entry["region"], entry["brick"])
        # Distinguish multiple local reductions of the *same* release/region
        # by exact import root, rather than conflating their provenance.
        skey = f"legacy_surveys/{release}/{region}/{root}"
        source = source_by_key.get(skey)
        if source is None:
            source = {
                "source_id": str(uuid.uuid4()),
                "import_key": skey,
                "provider": "legacy_surveys",
                "release": release,
                "region": region,
                "data_root": _path_ref(root, project.root),
                "observation_lineage": INDEPENDENCE_UNKNOWN,
                "calibration_status": "UNASSESSED",
                "created_at": _now(),
            }
            registry["sources"].append(source)
            source_by_key[skey] = source
            sources_added += 1
        for kind, product in entry["products"].items():
            if kind not in PRODUCT_KINDS:
                continue
            relative = product["relative_path"]
            file = root / relative
            # Skip nonexistent products; incomplete quartets are represented
            # through the inventory receipt, not as invented FITS files.
            if product["presence"] == "MISSING":
                continue
            if file.is_symlink():
                raise ValueError(f"Original FITS symlink not allowed: {file}")
            stat = file.stat()
            akey = f"{skey}/{relative}"
            existing = asset_by_key.get(akey)
            if existing is not None:
                if (existing.get("size_bytes") != stat.st_size or
                        existing.get("mtime_ns") != stat.st_mtime_ns):
                    existing["presence"] = "PRESENT_UNVERIFIED"
                    existing["checksum_status"] = (
                        CHECKSUM_CHANGED if existing.get("sha256") else
                        CHECKSUM_PENDING)
                    existing["fits_status"] = FITS_UNVERIFIED
                    existing["size_bytes"] = stat.st_size
                    existing["mtime_ns"] = stat.st_mtime_ns
                    assets_changed += 1
                continue
            asset = {
                "asset_id": str(uuid.uuid4()),
                "source_id": source["source_id"],
                "import_key": akey,
                "brick": brick,
                "band": None if kind == "maskbits" else entry["band"],
                "product_kind": kind,
                "product_family": "resampled_coadd",
                "path_ref": _path_ref(file, project.root),
                "size_bytes": stat.st_size,
                "mtime_ns": stat.st_mtime_ns,
                "presence": ("EMPTY_FILE" if stat.st_size == 0 else
                             "PRESENT_UNVERIFIED"),
                "sha256": None,
                "checksum_status": CHECKSUM_PENDING,
                "fits_status": FITS_UNVERIFIED,
                "exposure_lineage": INDEPENDENCE_UNKNOWN,
                "created_at": _now(),
            }
            registry["assets"].append(asset)
            asset_by_key[akey] = asset
            assets_added += 1
    _save(project.root, registry)
    return {
        "sources_added": sources_added,
        "assets_added": assets_added,
        "assets_changed": assets_changed,
        "sources_total": len(registry["sources"]),
        "assets_total": len(registry["assets"]),
        "science_status": "NOT_SCIENTIFICALLY_ADOPTED",
        "location": str(registry_path(project.root)),
    }


def register_external_fits(
    project_root: str | Path, file_path: str | Path, *,
    provider: str, band: str, description: str = "",
) -> dict:
    """Register a user-selected external FITS; no invented survey metadata."""
    if not provider.strip() or not band.strip():
        raise ValueError("External provider and passband must be supplied")
    project = _project(project_root)
    file = Path(file_path).expanduser().absolute()
    if not file.is_file() or file.is_symlink():
        raise FileNotFoundError("External input must be an existing ordinary FITS")
    if not file.name.lower().endswith((".fits", ".fit", ".fits.fz", ".fts")):
        raise ValueError("Unsupported external FITS filename")
    registry = read_dataset_index(project.root)
    file_ref = _path_ref(file, project.root)
    for existing in registry["assets"]:
        if (existing.get("path_ref") == file_ref and
                existing.get("band") == band.strip() and
                any(source["source_id"] == existing["source_id"] and
                    source.get("provider") == provider.strip()
                    for source in registry["sources"])):
            # Reattaching the same original is idempotent; keep provenance
            # and SHA-256 receipts instead of creating duplicate assets.
            return {"source_id": existing["source_id"],
                    "asset_id": existing["asset_id"]}
    source = {
        "source_id": str(uuid.uuid4()),
        "import_key": None,
        "provider": provider.strip(), "release": None, "region": None,
        "data_root": None, "observation_lineage": INDEPENDENCE_UNKNOWN,
        "calibration_status": "UNASSESSED", "description": description,
        "created_at": _now(),
    }
    stat = file.stat()
    asset = {
        "asset_id": str(uuid.uuid4()),
        "source_id": source["source_id"], "import_key": None,
        "brick": None, "band": band.strip(),
        "product_kind": "image", "product_family": "external_original_unknown",
        "path_ref": _path_ref(file, project.root),
        "size_bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns,
        "presence": "PRESENT_UNVERIFIED",
        "sha256": None, "checksum_status": CHECKSUM_PENDING,
        "fits_status": FITS_UNVERIFIED,
        "exposure_lineage": INDEPENDENCE_UNKNOWN, "created_at": _now(),
    }
    registry["sources"].append(source)
    registry["assets"].append(asset)
    _save(project.root, registry)
    return {"source_id": source["source_id"], "asset_id": asset["asset_id"]}


def verify_asset_checksum(project_root: str | Path, asset_id: str) -> dict:
    """Read-only file SHA-256 receipt; not scientific FITS/WCS validation.

    Previously verified content that changes is refused; re-import as a new
    version, instead of silently replacing a reproducibility dependency.
    """
    root = _project(project_root).root
    registry = read_dataset_index(root)
    asset = next((a for a in registry["assets"]
                  if a["asset_id"] == asset_id), None)
    if asset is None:
        raise KeyError(f"Unknown asset ID: {asset_id}")
    path = resolve_asset_path(root, asset)
    if not path.is_file() or path.is_symlink():
        raise FileNotFoundError(f"Indexed FITS missing: {path}")
    before = path.stat()
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    after = path.stat()
    if (before.st_size != after.st_size or
            before.st_mtime_ns != after.st_mtime_ns):
        raise RuntimeError("FITS changed during SHA256 verification")
    sha = digest.hexdigest()
    previous = asset.get("sha256")
    if previous is not None and previous != sha:
        raise RuntimeError("Indexed FITS content changed; import a new asset version")
    asset.update({
        "sha256": sha, "size_bytes": after.st_size,
        "mtime_ns": after.st_mtime_ns,
        "checksum_status": CHECKSUM_VERIFIED,
        "presence": "PRESENT_UNVERIFIED",  # FITS validation remains separate
    })
    _save(root, registry)
    return {"asset_id": asset_id, "sha256": sha,
            "checksum_status": CHECKSUM_VERIFIED,
            "fits_status": asset["fits_status"]}


def science_cache_key(
    *, stage: str, algorithm_version: str, parameters: dict,
    selected_assets: list[dict], dependency_keys: list[str] | None = None,
) -> str:
    """Content-addressed cache key; no reuse on unknown source checksums.

    The caller must bind every stage-relevant processing parameter explicitly.
    This is *not* a pipeline scheduler or a scientific validity certificate.
    """
    if not stage or not algorithm_version or not isinstance(parameters, dict):
        raise ValueError("Stage, algorithm version and parameter map are required")
    if not selected_assets:
        raise ValueError("Cannot cache a science result with no data assets")
    signatures = []
    for asset in selected_assets:
        sha = asset.get("sha256")
        if (asset.get("checksum_status") != CHECKSUM_VERIFIED or
                not isinstance(sha, str) or len(sha) != 64):
            raise ValueError("Unverified FITS fingerprint; refuse cache reuse")
        signatures.append((asset["asset_id"], sha))
    if len(signatures) != len(set(x[0] for x in signatures)):
        raise ValueError("Duplicate asset in cache dependencies")
    payload = {
        "stage": stage, "algorithm_version": algorithm_version,
        "parameters": parameters,
        "assets": sorted(signatures),
        "dependencies": sorted(dependency_keys or []),
    }
    serialized = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(serialized.encode("utf-8")).hexdigest()


def compatibility_policy(sources: list[dict], assets: list[dict],
                         *, mode: str) -> dict:
    """Fail closed for unknown exposure lineage or calibration compatibility.

    Planning policy only; does not itself build a mosaic or coadd.
    """
    if mode not in ("primary", "compare", "mosaic", "joint"):
        raise ValueError("Unknown dataset composition mode")
    if not sources or not assets:
        return {"allowed": False, "reason": "No source assets selected"}
    source_ids = {s["source_id"] for s in sources}
    if any(a["source_id"] not in source_ids for a in assets):
        return {"allowed": False, "reason": "Unknown source in selection"}
    if mode in ("primary", "compare"):
        return {"allowed": True,
                "reason": "Separate data selections; no pixel coaddition"}
    if mode == "joint":
        return {"allowed": False, "reason":
                "Joint coaddition is disabled until independently certified "
                "exposure lineage, calibration, PSF and noise covariance"}
    # Mosaic geometry is allowed only *as a plan* for single-release
    # same-region, same-band brick cutouts. This grants no photometric coadd.
    keys = {(s.get("provider"), s.get("release"), s.get("region"))
            for s in sources}
    bands = {a.get("band") for a in assets if a.get("band") is not None}
    if len(keys) != 1 or len(bands) != 1 or None in bands:
        return {"allowed": False, "reason":
                "Mosaic requires one provider/release/region and one band"}
    if any(a.get("product_family") != "resampled_coadd" for a in assets):
        return {"allowed": False, "reason": "Mosaic needs validated coadd assets"}
    return {"allowed": True, "reason":
            "Geometry plan only; overlap pixels are not independent. "
            "Choose unique coverage; audit PSF, sky and correlated noise "
            "before any science mosaic"}
