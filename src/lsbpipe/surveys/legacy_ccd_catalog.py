"""Offline, fail-closed Legacy Survey per-brick CCD table importer.

DR9/DR10/DR11 officially document coadd/<AAA>/<brick>/
legacysurvey-<brick>-ccds.fits, including CAMERA, EXPNUM, CCDNAME,
FILTER, and CCD/brick overlap extents. This adapter does not trust a
filename alone, does not fetch data, and never certifies image independence.
"""
from __future__ import annotations

import hashlib
from pathlib import Path
import re

from astropy.io import fits
import numpy as np

from .dataset_registry import (
    CHECKSUM_VERIFIED, read_dataset_index, resolve_asset_path,
)
from .object_projects import read_object_project
from .observation_lineage import (
    compare_coadd_lineage, read_observation_lineage, register_ccd_assertions,
)

_FILE_PATTERN = re.compile(r"legacysurvey-(\d{4}[pm]\d{3})-ccds\.fits(?:\.gz)?$")
_MANDATORY = frozenset((
    "camera", "expnum", "ccdname", "filter", "brick_x0", "brick_x1",
    "brick_y0", "brick_y1",
))
_MAX_ROWS = 200_000
# A documented namespace for all camera names, not an assumption that
# all cameras share a physical observatory or FITS detector naming system.
_ARCHIVE_NAMESPACE = "legacy_surveys"
_ALLOWED_CAMERAS = frozenset(("decam", "90prime", "mosaic"))


def _as_text(value, column: str) -> str:
    if np.ma.is_masked(value) or value is None:
        raise ValueError(f"{column}: null or masked cell")
    if isinstance(value, (bytes, np.bytes_)):
        try:
            value = value.decode("ascii")
        except UnicodeDecodeError as exc:
            raise ValueError(f"{column}: non-ASCII identifier") from exc
    if not isinstance(value, str):
        raise ValueError(f"{column}: string expected")
    text = value.strip()
    if not text or not text.isascii():
        raise ValueError(f"{column}: blank/non-ASCII value")
    return text


def _positive_expnum(value) -> str:
    if np.ma.is_masked(value) or isinstance(value, (bool, np.bool_)):
        raise ValueError("EXPNUM must be a positive integer")
    if not isinstance(value, (int, np.integer)):
        raise ValueError("EXPNUM must be an integer FITS column")
    if int(value) <= 0:
        raise ValueError("EXPNUM must be a positive integer")
    # The official catalogue uses an integer, not a zero-padded string.
    # Never invent leading zeros that cannot be recovered from FITS.
    return str(int(value))


def _bounded_int(value, name: str) -> int:
    if np.ma.is_masked(value) or isinstance(value, (bool, np.bool_)) or (
            not isinstance(value, (int, np.integer))):
        raise ValueError(f"{name}: integer coordinate required")
    result = int(value)
    if result < 0 or result > 100_000:
        raise ValueError(f"{name}: invalid detector/brick coordinate")
    return result


def _path_reference(file: Path, project_root: Path) -> dict:
    try:
        return {"relative": file.relative_to(project_root).as_posix()}
    except ValueError:
        return {"absolute": str(file)}


def _read_table(file: Path, band: str) -> tuple[list[dict], int, dict]:
    """Validate an unambiguous BINTABLE and every row before registration."""
    with fits.open(file, mode="readonly", memmap=False,
                   checksum=True, do_not_scale_image_data=True) as hdul:
        hdul.verify("exception")
        tables = [(i, hdu) for i, hdu in enumerate(hdul)
                  if isinstance(hdu, fits.BinTableHDU)]
        if len(tables) != 1:
            raise ValueError(
                "Expected exactly one CCD binary table; no HDU guessing")
        number, table = tables[0]
        names = table.columns.names
        if names is None:
            raise ValueError("Missing CCD table columns")
        columns = {str(n).lower(): n for n in names}
        if len(columns) != len(names) or not _MANDATORY.issubset(columns):
            raise ValueError(
                "Malformed CCD table: missing or duplicate mandatory columns "
                + ", ".join(sorted(_MANDATORY - set(columns))))
        # No assumption that nonempty file, identical filenames or a
        # generic CCD catalogue proves a complete list of source exposures.
        if len(table.data) == 0 or len(table.data) > _MAX_ROWS:
            raise ValueError("Empty or excessively large CCD table")
        rows = []
        seen = set()
        band_counts = {}
        for index, record in enumerate(table.data):
            filter_name = _as_text(record[columns["filter"]], "FILTER").lower()
            if filter_name not in ("g", "r", "i", "z"):
                raise ValueError(f"Unsupported FILTER in CCD row {index}")
            band_counts[filter_name] = band_counts.get(filter_name, 0) + 1
            camera = _as_text(record[columns["camera"]], "CAMERA").lower()
            if camera not in _ALLOWED_CAMERAS:
                raise ValueError(f"Unrecognized CAMERA in row {index}")
            expnum = _positive_expnum(record[columns["expnum"]])
            ccdname = _as_text(record[columns["ccdname"]], "CCDNAME")
            bx0 = _bounded_int(record[columns["brick_x0"]], "BRICK_X0")
            bx1 = _bounded_int(record[columns["brick_x1"]], "BRICK_X1")
            by0 = _bounded_int(record[columns["brick_y0"]], "BRICK_Y0")
            by1 = _bounded_int(record[columns["brick_y1"]], "BRICK_Y1")
            if bx1 <= bx0 or by1 <= by0:
                raise ValueError(
                    f"CCD row {index} has nonpositive brick overlap")
            key = (camera, expnum, ccdname.upper(), filter_name)
            if key in seen:
                raise ValueError(
                    f"Duplicate physical CCD in catalogue, row {index}")
            seen.add(key)
            if filter_name == band:
                rows.append({
                    "archive": _ARCHIVE_NAMESPACE,
                    "instrument": camera,
                    "exposure_id": expnum,
                    "detector_id": ccdname,
                    "band": filter_name,
                    "record_locator": f"hdu:{number}:row:{index}",
                })
        if not rows:
            raise ValueError(
                f"No CCD catalogue rows in requested band {band}")
        checksum_states = []
        for hdu in hdul:
            if "CHECKSUM" in hdu.header:
                if hdu.verify_checksum() != 1:
                    raise ValueError("Invalid FITS CHECKSUM")
                checksum_states.append("CHECKSUM_VALID")
            if "DATASUM" in hdu.header:
                if hdu.verify_datasum() != 1:
                    raise ValueError("Invalid FITS DATASUM")
                checksum_states.append("DATASUM_VALID")
        extra = {
            "bintable_hdu": number, "catalogue_rows": len(table.data),
            "filter_counts": band_counts,
            "fits_integrity": (
                "FITS_CHECKSUM_VERIFIED" if checksum_states
                else "FITS_NO_EMBEDDED_CHECKSUM"),
        }
    return rows, len(rows), extra


def inspect_brick_ccd_catalog(
    project_root: str | Path, *, coadd_asset_id: str,
    catalogue_file: str | Path,
) -> dict:
    """Preflight an explicitly selected original coadd + local CCD file.

    Inspect never persists lineage. Its result is an untrusted catalogue
    assertion receipt even when the physical source file SHA-256 is computed.
    """
    project = read_object_project(project_root)
    index = read_dataset_index(project.root)
    asset = next((a for a in index["assets"]
                  if a["asset_id"] == coadd_asset_id), None)
    if (asset is None or asset.get("product_kind") != "image" or
            asset.get("product_family") != "resampled_coadd"):
        raise ValueError("Selected asset must be indexed original coadd IMAGE")
    if asset.get("checksum_status") != CHECKSUM_VERIFIED or (
            not isinstance(asset.get("sha256"), str)
            or len(asset["sha256"]) != 64):
        raise ValueError("Original coadd SHA-256 must be verified first")
    source = next((s for s in index["sources"]
                   if s["source_id"] == asset["source_id"]), None)
    if source is None or source.get("provider") != "legacy_surveys" or (
            source.get("release") not in ("dr9", "dr10", "dr11")) or (
            source.get("region") not in ("north", "south")):
        raise ValueError("Unsupported coadd source provenance")
    original = resolve_asset_path(project.root, asset)
    if not original.is_file() or original.is_symlink():
        raise FileNotFoundError("Indexed source coadd no longer exists")
    original_stat = original.stat()
    if (original_stat.st_size != asset["size_bytes"]
            or original_stat.st_mtime_ns != asset["mtime_ns"]):
        raise ValueError("Indexed coadd stat changed after checksum")
    # The existing checksum is a receipt, not evidence of FITS/WCS validity.
    file = Path(catalogue_file).expanduser().absolute()
    if file.is_symlink() or not file.is_file():
        raise ValueError("Catalogue must be an existing nonsymlink local file")
    name = _FILE_PATTERN.fullmatch(file.name)
    if not name or name.group(1) != asset.get("brick"):
        raise ValueError("Catalogue filename brick differs from indexed coadd")
    if asset.get("band") not in ("g", "r", "i", "z"):
        raise ValueError("Unsupported source band")
    before = file.stat()
    digest = hashlib.sha256()
    with file.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(block)
    sha = digest.hexdigest()
    rows, selected_count, extra = _read_table(file, asset["band"])
    after = file.stat()
    if ((before.st_size, before.st_mtime_ns) !=
            (after.st_size, after.st_mtime_ns)):
        raise RuntimeError("CCD catalogue changed during inspection")
    release = source["release"]
    region = source["region"]
    brick = asset["brick"]
    # Exact input digest participates in catalogue provenance identity.
    # The one per-brick CCD table is deliberately reused for g/r/i/z;
    # receipts are scoped to the exact registered image asset and band.
    ref = (f"legacy_brick_ccds:{release}:{region}:{brick}:"
           f"{asset['band']}:asset:{coadd_asset_id}:sha256:{sha}")
    return {
        "catalogue_reference": ref,
        "catalogue_path_ref": _path_reference(file, project.root),
        "catalogue_sha256": sha,
        "catalogue_size_bytes": after.st_size,
        "catalogue_mtime_ns": after.st_mtime_ns,
        "coadd_asset_id": coadd_asset_id,
        "coadd_sha256": asset["sha256"],
        "source_id": asset["source_id"],
        "release": release, "region": region, "brick": brick,
        "band": asset["band"],
        "records": rows, "selected_ccd_count": selected_count,
        **extra,
        "status": "CATALOGUE_ASSERTION_NOT_VALIDATED",
        "independence_certified": False,
        "joint_scientific_combination_allowed": False,
    }


def import_brick_ccd_catalog(
    project_root: str | Path, *, coadd_asset_id: str,
    catalogue_file: str | Path,
) -> dict:
    """Preflight, then persist all CCD identities/receipts in one atomic write.

    Repeated imports of identical bytes are idempotent. Different catalogue
    content for the same exact source asset/brick requires explicit review.
    """
    preflight = inspect_brick_ccd_catalog(
        project_root, coadd_asset_id=coadd_asset_id,
        catalogue_file=catalogue_file)
    # Do not persist a report if the catalogue was replaced after preflight.
    file = Path(catalogue_file).expanduser().absolute()
    current = hashlib.sha256()
    with file.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            current.update(chunk)
    if current.hexdigest() != preflight["catalogue_sha256"]:
        raise RuntimeError("CCD catalogue content changed after inspection")
    # Another process may also have changed/reindexed the source image.
    source_asset = next((a for a in read_dataset_index(project_root)["assets"]
                         if a["asset_id"] == coadd_asset_id), None)
    if (source_asset is None or source_asset.get("sha256") !=
            preflight["coadd_sha256"] or
            source_asset.get("checksum_status") != CHECKSUM_VERIFIED):
        raise RuntimeError("Coadd asset provenance changed during CCD import")
    previous = read_observation_lineage(project_root)
    for receipt in previous.get("catalogue_receipts", []):
        if (receipt["coadd_asset_id"] == coadd_asset_id and
                receipt["catalogue_sha256"] != preflight["catalogue_sha256"]):
            raise ValueError(
                "CCD catalogue changed; do not overwrite existing "
                "lineage without explicit version review")
    record = {
        key: value for key, value in preflight.items()
        if key not in ("records",)
    }
    result = register_ccd_assertions(
        project_root, coadd_asset_id=coadd_asset_id,
        catalogue_reference=preflight["catalogue_reference"],
        records=preflight["records"], catalogue_receipt=record)
    return {
        **result,
        "catalogue_reference": preflight["catalogue_reference"],
        "catalogue_sha256": preflight["catalogue_sha256"],
        "catalogue_rows": preflight["catalogue_rows"],
        "selected_ccd_count": preflight["selected_ccd_count"],
        "fits_integrity": preflight["fits_integrity"],
        "science_status": "NOT_SCIENTIFICALLY_ADOPTED",
    }
