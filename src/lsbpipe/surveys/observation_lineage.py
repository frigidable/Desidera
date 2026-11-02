"""P0.6: conservative observational identity and coadd CCD-lineage assertions.

This is a metadata-only, opt-in sidecar. It does not download or open FITS,
adopt scientific inputs, certify source independence, or compute photometry.
Catalogue-derived relationships remain unverified until separately audited.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re

from .dataset_registry import read_dataset_index
from .local_project import _atomic_json
from .object_projects import read_object_project

SCHEMA = "desidera.observation-lineage/1"
UNVERIFIED = "CATALOGUE_ASSERTION_NOT_VALIDATED"
_TOKEN = re.compile(r"^[A-Za-z0-9_.:+-]{1,128}$")


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _canonical_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      allow_nan=False, ensure_ascii=True)


def _key_part(value: str, field: str, *, case: str = "") -> str:
    if not isinstance(value, str):
        raise ValueError(f"{field} must be a string; preserve leading zeros")
    value = value.strip()
    if not _TOKEN.fullmatch(value):
        raise ValueError(f"Invalid {field}: expected a nonempty identifier")
    if case == "lower":
        return value.lower()
    if case == "upper":
        return value.upper()
    return value


@dataclass(frozen=True)
class CCDIdentity:
    """Release-independent *catalogue-asserted* CCD identity.

    Avoid normalizing numeric exposure identifiers: '00123' and '123'
    cannot be equated without an authoritative crosswalk. The archive
    namespace prevents collisions across observatories/instruments.
    """
    archive: str
    instrument: str
    exposure_id: str
    detector_id: str

    def __post_init__(self) -> None:
        for name, case in (("archive", "lower"), ("instrument", "lower"),
                           ("exposure_id", ""), ("detector_id", "upper")):
            object.__setattr__(
                self, name, _key_part(getattr(self, name), name, case=case))

    def as_dict(self) -> dict:
        return {
            "archive": self.archive, "instrument": self.instrument,
            "exposure_id": self.exposure_id, "detector_id": self.detector_id,
        }

    @property
    def ccd_id(self) -> str:
        digest = hashlib.sha256(
            _canonical_json(self.as_dict()).encode("ascii")).hexdigest()
        return "ccd-" + digest

    @property
    def exposure_id_stable(self) -> str:
        without_detector = {
            "archive": self.archive, "instrument": self.instrument,
            "exposure_id": self.exposure_id,
        }
        digest = hashlib.sha256(
            _canonical_json(without_detector).encode("ascii")).hexdigest()
        return "exp-" + digest


def lineage_path(project_root: str | Path) -> Path:
    """Separate optional metadata; the existing dataset index is unchanged."""
    project = read_object_project(project_root)
    return Path(project.root) / "datasets" / "observation-lineage.json"


def _empty(project_id: str) -> dict:
    return {
        "schema": SCHEMA,
        "project_id": project_id,
        "created_at": _now(),
        "observations": [],
        "contributions": [],
        "catalogue_receipts": [],
        "scientific_status": "NOT_SCIENTIFICALLY_ADOPTED",
    }


def _validated_observations(record: dict, index: dict) -> None:
    observations = record.get("observations")
    contributions = record.get("contributions")
    if not isinstance(observations, list) or not isinstance(contributions, list):
        raise ValueError("Lineage requires observation/contribution lists")
    if record.get("scientific_status") != "NOT_SCIENTIFICALLY_ADOPTED":
        raise ValueError("Lineage sidecar cannot adopt scientific results")
    identities = {}
    for item in observations:
        if not isinstance(item, dict) or not isinstance(item.get("identity"), dict):
            raise ValueError("Malformed observation identity")
        raw = item["identity"]
        if set(raw) != {"archive", "instrument", "exposure_id", "detector_id"}:
            raise ValueError("Unexpected observation identity fields")
        identity = CCDIdentity(**raw)
        if identity.as_dict() != raw or item.get("observation_id") != identity.ccd_id:
            raise ValueError("Observation ID or canonical identity mismatch")
        if item.get("exposure_id") != identity.exposure_id_stable:
            raise ValueError("Exposure group ID mismatch")
        if identity.ccd_id in identities:
            raise ValueError("Duplicate observation identity")
        identities[identity.ccd_id] = identity

    assets = {a["asset_id"]: a for a in index["assets"]}
    sources = {source["source_id"]: source for source in index["sources"]}
    receipts = record.get("catalogue_receipts", [])
    if not isinstance(receipts, list):
        raise ValueError("Invalid CCD catalogue receipt collection")
    by_reference = {}
    for receipt in receipts:
        if not isinstance(receipt, dict):
            raise ValueError("Malformed CCD catalogue receipt")
        ref = receipt.get("catalogue_reference")
        sha = receipt.get("catalogue_sha256")
        coadd = assets.get(receipt.get("coadd_asset_id"))
        source = sources.get(receipt.get("source_id"))
        if (not isinstance(ref, str) or not ref.startswith(
                "legacy_brick_ccds:") or not isinstance(sha, str) or
                len(sha) != 64 or any(ch not in "0123456789abcdef"
                                     for ch in sha)):
            raise ValueError("Invalid CCD catalogue digest/reference")
        if ref in by_reference:
            raise ValueError("Duplicate CCD catalogue reference")
        if (coadd is None or coadd.get("product_kind") != "image" or
                coadd.get("product_family") != "resampled_coadd" or
                source is None or source["source_id"] != coadd["source_id"]):
            raise ValueError("CCD catalogue refers to wrong coadd/source")
        if (receipt.get("coadd_sha256") != coadd.get("sha256") or
                receipt.get("band") != coadd.get("band") or
                receipt.get("brick") != coadd.get("brick") or
                receipt.get("release") != source.get("release") or
                receipt.get("region") != source.get("region")):
            raise ValueError("Catalogue receipt disagrees with indexed coadd")
        expected = (
            f"legacy_brick_ccds:{source['release']}:"
            f"{source['region']}:{coadd['brick']}:sha256:{sha}")
        if ref != expected:
            raise ValueError("CCD catalogue reference does not match receipt")
        if (not isinstance(receipt.get("catalogue_path_ref"), dict) or
                not isinstance(receipt.get("catalogue_size_bytes"), int) or
                receipt["catalogue_size_bytes"] <= 0 or
                not isinstance(receipt.get("catalogue_rows"), int) or
                receipt["catalogue_rows"] <= 0 or
                not isinstance(receipt.get("selected_ccd_count"), int) or
                receipt["selected_ccd_count"] <= 0 or
                receipt.get("status") != UNVERIFIED or
                receipt.get("independence_certified") is not False or
                receipt.get("joint_scientific_combination_allowed") is not False):
            raise ValueError("Incomplete or invalid catalogue import receipt")
        by_reference[ref] = receipt
    keys = set()
    locators = {}
    for item in contributions:
        if not isinstance(item, dict):
            raise ValueError("Malformed lineage contribution")
        aid, oid = item.get("asset_id"), item.get("observation_id")
        asset = assets.get(aid)
        if asset is None or asset.get("product_kind") != "image" or (
                asset.get("product_family") != "resampled_coadd"):
            raise ValueError("Lineage references unknown/non-coadd IMAGE asset")
        if oid not in identities or item.get("evidence_state") != UNVERIFIED:
            raise ValueError("Unknown observation or unverified evidence state changed")
        if item.get("source_id") != asset["source_id"]:
            raise ValueError("Lineage source does not match indexed asset")
        if item.get("band") != asset.get("band"):
            raise ValueError("Lineage band does not match indexed asset")
        ref = item.get("catalogue_reference")
        locator = item.get("record_locator")
        if not isinstance(ref, str) or not ref.strip() or (
                not isinstance(locator, str) or not locator.strip()):
            raise ValueError("Missing catalogue provenance")
        if ref.startswith("legacy_brick_ccds:"):
            receipt = by_reference.get(ref)
            if receipt is None or receipt["coadd_asset_id"] != aid:
                raise ValueError("CCD assertion missing source-file receipt")
        key = (aid, oid, ref, locator)
        if key in keys:
            raise ValueError("Duplicate contribution assertion")
        keys.add(key)
        claim = (aid, ref, locator)
        if claim in locators and locators[claim] != oid:
            raise ValueError("Conflicting CCD identity at one catalogue locator")
        locators[claim] = oid


def read_observation_lineage(project_root: str | Path) -> dict:
    """Read-only; absence is UNKNOWN, never evidence of no observations."""
    project = read_object_project(project_root)
    path = lineage_path(project.root)
    if not path.is_file():
        return _empty(project.project_id)
    try:
        record = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ValueError(f"Invalid observational lineage sidecar: {exc}") from exc
    if not isinstance(record, dict) or record.get("schema") != SCHEMA:
        raise ValueError("Unsupported lineage schema; no automatic migration")
    if record.get("project_id") != project.project_id:
        raise ValueError("Lineage sidecar belongs to a different galaxy")
    _validated_observations(record, read_dataset_index(project.root))
    return record


def register_ccd_assertions(
    project_root: str | Path, *, coadd_asset_id: str,
    catalogue_reference: str, records: list[dict],
    catalogue_receipt: dict | None = None,
) -> dict:
    """Index explicit catalogue statements about CCDs contributing to a coadd.

    This is NOT a parser or validation of a Legacy brick-CCDs FITS table.
    Callers must provide each exact detector identity and catalogue row key.
    Commit is all-or-nothing; no original data or old registry is mutated.
    Re-importing an identical set is idempotent, preserving UUIDs/timestamps.
    """
    if not isinstance(catalogue_reference, str) or not catalogue_reference.strip():
        raise ValueError("A catalogue reference is required")
    if not isinstance(records, list) or not records:
        raise ValueError("At least one catalogue record is required")
    project = read_object_project(project_root)
    index = read_dataset_index(project.root)
    assets = {a["asset_id"]: a for a in index["assets"]}
    asset = assets.get(coadd_asset_id)
    if not asset or asset.get("product_kind") != "image" or (
            asset.get("product_family") != "resampled_coadd"):
        raise ValueError("CCD lineage requires an indexed native coadd IMAGE")
    state = read_observation_lineage(project.root)
    by_ccd = {o["observation_id"]: o for o in state["observations"]}
    by_key = {(c["asset_id"], c["observation_id"],
               c["catalogue_reference"], c["record_locator"])
              for c in state["contributions"]}
    by_locator = {(c["asset_id"], c["catalogue_reference"],
                   c["record_locator"]): c["observation_id"]
                  for c in state["contributions"]}
    observations_added = contributions_added = 0
    receipt_added = False
    if catalogue_receipt is not None:
        if (not isinstance(catalogue_receipt, dict) or
                catalogue_receipt.get("catalogue_reference") !=
                catalogue_reference or
                catalogue_receipt.get("coadd_asset_id") != coadd_asset_id):
            raise ValueError("Catalogue receipt does not match import target")
        current_receipts = state.setdefault("catalogue_receipts", [])
        matches = [r for r in current_receipts
                   if r.get("catalogue_reference") == catalogue_reference]
        if matches:
            if matches != [catalogue_receipt]:
                raise ValueError("Conflicting catalogue receipt for same digest")
        else:
            # Reject a new revision of the same catalogue for this asset.
            # Version supersession must be an explicit separate operation.
            if any(r.get("coadd_asset_id") == coadd_asset_id and
                   r.get("catalogue_sha256") !=
                   catalogue_receipt.get("catalogue_sha256")
                   for r in current_receipts):
                raise ValueError("Conflicting catalogue revision for coadd")
            current_receipts.append(dict(catalogue_receipt))
            receipt_added = True
    for row in records:
        if not isinstance(row, dict):
            raise ValueError("CCD catalogue rows must be maps")
        required = {"archive", "instrument", "exposure_id",
                    "detector_id", "band", "record_locator"}
        if set(row) != required:
            raise ValueError("CCD row requires exact identity/band/locator fields")
        ident = CCDIdentity(row["archive"], row["instrument"],
                            row["exposure_id"], row["detector_id"])
        if row["band"] != asset.get("band"):
            raise ValueError("CCD band conflicts with indexed coadd band")
        locator = row["record_locator"]
        if not isinstance(locator, str) or not locator.strip():
            raise ValueError("Catalogue row locator is required")
        loc = (coadd_asset_id, catalogue_reference, locator)
        if loc in by_locator and by_locator[loc] != ident.ccd_id:
            raise ValueError("Conflicting CCD identity at catalogue locator")
        by_locator[loc] = ident.ccd_id
        if ident.ccd_id not in by_ccd:
            obs = {
                "observation_id": ident.ccd_id,
                "exposure_id": ident.exposure_id_stable,
                "identity": ident.as_dict(),
                "created_at": _now(),
            }
            state["observations"].append(obs)
            by_ccd[ident.ccd_id] = obs
            observations_added += 1
        key = (coadd_asset_id, ident.ccd_id, catalogue_reference, locator)
        if key not in by_key:
            state["contributions"].append({
                "asset_id": coadd_asset_id, "source_id": asset["source_id"],
                "observation_id": ident.ccd_id, "band": asset["band"],
                "catalogue_reference": catalogue_reference,
                "record_locator": locator, "evidence_state": UNVERIFIED,
                "created_at": _now(),
            })
            by_key.add(key)
            contributions_added += 1
    _validated_observations(state, index)
    if observations_added or contributions_added or receipt_added:
        _atomic_json(lineage_path(project.root), state)
    return {
        "observations_added": observations_added,
        "contributions_added": contributions_added,
        "catalogue_receipt_added": receipt_added,
        "total_observations": len(state["observations"]),
        "total_contributions": len(state["contributions"]),
        "evidence_state": UNVERIFIED,
        "science_status": state["scientific_status"],
    }


def compare_coadd_lineage(project_root: str | Path,
                          first_asset_id: str, second_asset_id: str) -> dict:
    """Conservative overlap signal, NEVER a certificate of independence."""
    index = read_dataset_index(project_root)
    images = {a["asset_id"]: a for a in index["assets"] if (
        a.get("product_family") == "resampled_coadd" and
        a.get("product_kind") == "image")}
    if first_asset_id not in images or second_asset_id not in images:
        raise ValueError("Only indexed original coadd IMAGE assets can be compared")
    state = read_observation_lineage(project_root)
    ccd_to_exposure = {o["observation_id"]: o["exposure_id"]
                       for o in state["observations"]}
    a = {c["observation_id"] for c in state["contributions"]
         if c["asset_id"] == first_asset_id}
    b = {c["observation_id"] for c in state["contributions"]
         if c["asset_id"] == second_asset_id}
    same_ccd = sorted(a & b)
    same_exposure = sorted({ccd_to_exposure[x] for x in a} &
                           {ccd_to_exposure[x] for x in b})
    if first_asset_id == second_asset_id:
        status = "SAME_ASSET"
    elif same_ccd:
        status = "SHARED_CCD_IDS_ASSERTED"
    elif same_exposure:
        status = "SHARED_EXPOSURE_IDS_ASSERTED"
    elif not a or not b:
        status = "UNKNOWN_LINEAGE_NOT_INDEXED"
    else:
        status = "NO_SHARED_REGISTERED_IDS_NOT_CERTIFIED"
    return {
        "status": status,
        "first_ccd_count": len(a), "second_ccd_count": len(b),
        "shared_ccd_ids": same_ccd,
        "shared_exposure_ids": same_exposure,
        "catalogue_evidence": UNVERIFIED,
        "independence_certified": False,
        "joint_scientific_combination_allowed": False,
    }
