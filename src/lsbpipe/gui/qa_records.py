"""Immutable GUI science context and write-once original-FITS QA records.

Run directories are siblings of the original FITS data root, never inside
the release/brick tree. Every analysis produces a distinct, atomic record.
SHA-256 values are *local* file fingerprints, not archive checksums.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from hashlib import sha256
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from contextlib import contextmanager
import json
import os
import secrets
import shutil
import subprocess
import re
from functools import lru_cache

from .run_contracts import make_run_contract, check_run_receipt
from lsbpipe.surveys.legacy import BrickHit, SkyField, product_for


def _path_identity(value: str | Path) -> str:
    if not str(value).strip():
        return ""
    return os.path.normcase(str(Path(value).expanduser().resolve(strict=False)))


@dataclass(frozen=True)
class ScienceContext:
    """Frozen job parameters; can be compared after a background result."""

    action: str
    release: str
    region: str
    brick: str
    band: str
    ra_deg: float
    dec_deg: float
    width_arcmin: float
    target_radius_arcsec: float
    sky_inner_arcsec: float
    sky_outer_arcsec: float
    directory: str
    north_directory: str
    south_directory: str
    project_id: str | None = None

    @classmethod
    def capture(cls, *, action: str, field: SkyField, hit: BrickHit,
                band: str, directory: str | Path,
                north_directory: str | Path = "",
                south_directory: str | Path = "",
                target_radius_arcsec: float = 35.,
                sky_inner_arcsec: float = 65.,
                sky_outer_arcsec: float = 110.,
                project_id: str | None = None):
        return cls(
            action=action, release=hit.release, region=hit.region,
            brick=hit.name, band=band, ra_deg=float(field.ra_deg),
            dec_deg=float(field.dec_deg), width_arcmin=float(field.width_arcmin),
            target_radius_arcsec=float(target_radius_arcsec),
            sky_inner_arcsec=float(sky_inner_arcsec),
            sky_outer_arcsec=float(sky_outer_arcsec),
            directory=_path_identity(directory),
            north_directory=_path_identity(north_directory),
            south_directory=_path_identity(south_directory),
            project_id=project_id,
        )

    def to_dict(self) -> dict:
        return asdict(self)


def original_quartet(data_root: Path, hit: BrickHit, band: str) -> dict[str, Path]:
    return {
        kind: Path(data_root) / product_for(
            hit, kind, None if kind == "maskbits" else band).relative_path
        for kind in ("image", "invvar", "nexp", "maskbits")
    }


def _hash_file(path: Path) -> dict:
    before = path.stat()
    digest = sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (
            after.st_size, after.st_mtime_ns):
        raise RuntimeError(f"Input changed while hashing: {path}")
    return {
        "path": str(path.resolve()),
        "bytes": after.st_size,
        "mtime_ns": after.st_mtime_ns,
        "sha256_local": digest.hexdigest(),
    }


def snapshot_inputs(quartets: dict[str, dict[str, Path]]) -> dict:
    """Record authoritative local paths and nonauthoritative local hashes."""
    result = {}
    for label, quartet in sorted(quartets.items()):
        result[label] = {
            name: _hash_file(path)
            for name, path in sorted(quartet.items())
        }
    return result


def assert_inputs_unchanged(snapshot: dict):
    """Reject changed files before publishing any completed QA record."""
    for quartet in snapshot.values():
        for info in quartet.values():
            current = Path(info["path"]).stat()
            if ((current.st_size, current.st_mtime_ns) != (
                    info["bytes"], info["mtime_ns"]) or
                    _hash_file(Path(info["path"]))["sha256_local"]
                    != info["sha256_local"]):
                raise RuntimeError(
                    "An original FITS changed during science QA: "
                    + info["path"])


def _version() -> str:
    try:
        return version("lsbpipe-prototype")
    except PackageNotFoundError:
        return "unknown-local-source"


@lru_cache(maxsize=1)
def _git_revision() -> str | None:
    """Current source revision when available; no network or shell execution."""
    revision = os.environ.get("LSBPIPE_BUILD_COMMIT") or os.environ.get("GITHUB_SHA")
    if revision and re.fullmatch(r"[0-9a-fA-F]{40}", revision):
        return revision.lower()
    try:
        checked = subprocess.run(
            ["git", "rev-parse", "--verify", "HEAD"],
            cwd=Path(__file__).resolve().parent, shell=False,
            capture_output=True, text=True, timeout=3, check=False)
    except (OSError, subprocess.TimeoutExpired):
        return None
    value = checked.stdout.strip().lower()
    return value if (checked.returncode == 0 and
                     re.fullmatch(r"[0-9a-f]{40}", value)) else None


@contextmanager
def new_qa_run(*, data_root: Path, context: ScienceContext,
               inputs: dict, project_root: Path | None = None,
               dependencies: tuple[str, ...] = ()):
    """Stage artifacts privately and publish the *whole* run on success.

    On errors, remove staging but preserve every earlier result. An
    incomplete run never appears as a completed QA record.
    """
    output_root = (
        Path(project_root).expanduser().resolve()
        if project_root is not None else Path(data_root).resolve().parent
    ) / "qa-runs"
    output_root.mkdir(parents=True, exist_ok=True)

    def check_dependencies() -> None:
        for ident in dependencies:
            candidate = output_root / ident / "manifest.json"
            result = check_run_receipt(
                candidate, project_id=context.project_id)
            if result.get("state") != "RECORDED_NOT_REVALIDATED":
                raise ValueError(
                    "Run dependency missing, modified or wrong project: "
                    f"{ident} ({result.get('state')})")

    # Dependencies are immutable, verified artifact references, not filenames
    # or user-provided paths. Validation is repeated before publishing.
    if dependencies:
        check_dependencies()
    stamp = datetime.now(timezone.utc)
    run_id = f"{stamp:%Y%m%dT%H%M%S%fZ}-{secrets.token_hex(5)}"
    staging = output_root / ("." + run_id + ".partial")
    committed = output_root / run_id
    staging.mkdir(exist_ok=False)
    try:
        record = {
            "schema": "lsbpipe.gui-science-qa/1",
            "run_id": run_id,
            "created_at_utc": stamp.isoformat(),
            "science_status": "DIAGNOSTIC_NOT_ADOPTED",
            "lsbpipe_version": _version(),
            "context": context.to_dict(),
            "run_contract": make_run_contract(
                context=context.to_dict(), inputs=inputs,
                project_id=context.project_id,
                dependencies=dependencies, code_revision=_git_revision()),
            "input_fingerprint_type": (
                "local SHA-256; no authoritative remote digest supplied"),
            "input_original_fits": inputs,
            "artifacts": {},
        }
        yield staging, record
        assert_inputs_unchanged(inputs)
        if dependencies:
            check_dependencies()
        record["artifacts"] = {
            file.name: {
                "sha256_local": _hash_file(file)["sha256_local"],
                "bytes": file.stat().st_size,
            }
            for file in staging.iterdir() if file.is_file()
        }
        (staging / "manifest.json").write_text(
            json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False),
            encoding="utf-8")
        if committed.exists():
            raise FileExistsError(f"QA record already exists: {committed}")
        staging.rename(committed)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
