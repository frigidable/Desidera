"""Sky-coordinate vector masks: non-destructive multi-band manual-edit drafts.

The same ICRS region is rasterized independently on each native WCS.
No native MASKBITS, invalid/zero-IVAR, saturation, or no-coverage pixel
can be unmasked by a manual edit. No photometric result is created.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
import math
from typing import Literal

import numpy as np

SCHEMA = "desidera.sky-mask-draft/1"
DRAFT_STATUS = "MASK_DRAFT_NOT_ADOPTED"
PURPOSES = frozenset(("exclude", "restore_auto", "protect_target", "sky_train", "sky_holdout"))
BANDS = frozenset(("g", "r", "i", "z"))


def _angle(value: float) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError("Non-finite sky coordinate")
    return number


def _position(pair: tuple[float, float]) -> tuple[float, float]:
    if not isinstance(pair, (tuple, list)) or len(pair) != 2:
        raise ValueError("A sky position needs RA and Dec")
    ra, dec = map(_angle, pair)
    if not 0 <= ra < 360 or not -90 <= dec <= 90:
        raise ValueError("Position must be valid ICRS degrees")
    return ra, dec


@dataclass(frozen=True)
class SkyMaskShape:
    """A sky region independent of CCD/brick/filter image coordinates.

    Polygon vertices are geodesic-small-field approximations: rasterization
    uses a local sky-offset plane; only regions <= 0.5 degree are supported
    to limit projection errors pending a spherical polygon backend.
    """
    kind: Literal["circle", "ellipse", "polygon"]
    center: tuple[float, float] | None = None
    radius_arcsec: float | None = None
    a_arcsec: float | None = None
    q: float | None = None
    pa_deg: float | None = None
    vertices: tuple[tuple[float, float], ...] = ()

    def __post_init__(self):
        if self.kind not in ("circle", "ellipse", "polygon"):
            raise ValueError("Unsupported mask shape")
        if self.kind == "polygon":
            if self.center is not None or self.radius_arcsec is not None or self.a_arcsec is not None:
                raise ValueError("Polygon shape cannot have a circle/ellipse center")
            if not 3 <= len(self.vertices) <= 256:
                raise ValueError("Polygon needs 3–256 vertices")
            points = tuple(_position(p) for p in self.vertices)
            if len(set(points)) < 3:
                raise ValueError("Polygon must have at least three distinct vertices")
            object.__setattr__(self, "vertices", points)
            if self.q is not None or self.pa_deg is not None:
                raise ValueError("Polygon has no axial ratio/position angle")
        else:
            if self.vertices:
                raise ValueError("Circle/ellipse must not have polygon vertices")
            object.__setattr__(self, "center", _position(self.center))
            if self.kind == "circle":
                radius = _angle(self.radius_arcsec)
                if not 0 < radius <= 1800:
                    raise ValueError("Circle radius must be (0,1800] arcsec")
                if any(v is not None for v in (self.a_arcsec, self.q, self.pa_deg)):
                    raise ValueError("Circle may only specify radius")
            else:
                a, q, pa = _angle(self.a_arcsec), _angle(self.q), _angle(self.pa_deg)
                if not 0 < a <= 1800 or not 0 < q <= 1:
                    raise ValueError("Ellipse semimajor/q outside allowed range")
                object.__setattr__(self, "pa_deg", pa % 180.0)
                if self.radius_arcsec is not None:
                    raise ValueError("Ellipse may not specify circle radius")

    def to_dict(self) -> dict:
        return {"kind": self.kind, "center": list(self.center) if self.center else None,
                "radius_arcsec": self.radius_arcsec, "a_arcsec": self.a_arcsec,
                "q": self.q, "pa_deg": self.pa_deg,
                "vertices": [list(p) for p in self.vertices]}

    @classmethod
    def from_dict(cls, info: dict) -> "SkyMaskShape":
        return cls(kind=info["kind"], center=info.get("center"),
                   radius_arcsec=info.get("radius_arcsec"), a_arcsec=info.get("a_arcsec"),
                   q=info.get("q"), pa_deg=info.get("pa_deg"),
                   vertices=tuple(tuple(p) for p in info.get("vertices", ())))


@dataclass(frozen=True)
class SkyMaskEdit:
    """One revisioned researcher action; default applies across optical bands."""
    edit_id: str
    purpose: str
    shape: SkyMaskShape
    band: str | None = None
    comment: str = ""

    def __post_init__(self):
        if (not isinstance(self.edit_id, str) or not self.edit_id or
                len(self.edit_id) > 96 or not all(c.isalnum() or c in "_-" for c in self.edit_id)):
            raise ValueError("Unsafe mask edit identifier")
        if self.purpose not in PURPOSES:
            raise ValueError("Unsupported edit purpose")
        if self.band is not None and self.band not in BANDS:
            raise ValueError("Mask scope must be all bands or one valid optical band")
        if not isinstance(self.shape, SkyMaskShape):
            raise TypeError("Mask edit requires a sky shape")
        if not isinstance(self.comment, str) or len(self.comment) > 512:
            raise ValueError("Mask comment must be <=512 characters")

    def applies_to(self, band: str) -> bool:
        if band not in BANDS:
            raise ValueError("Unsupported optical band")
        return self.band is None or self.band == band

    def to_dict(self) -> dict:
        return {"edit_id": self.edit_id, "purpose": self.purpose,
                "shape": self.shape.to_dict(), "band": self.band, "comment": self.comment}

    @classmethod
    def from_dict(cls, info: dict) -> "SkyMaskEdit":
        return cls(edit_id=info["edit_id"], purpose=info["purpose"],
                   shape=SkyMaskShape.from_dict(info["shape"]),
                   band=info.get("band"), comment=info.get("comment", ""))


@dataclass(frozen=True)
class SkyMaskDraft:
    project_id: str
    edits: tuple[SkyMaskEdit, ...] = ()
    revision: int = 0

    def __post_init__(self):
        if not isinstance(self.project_id, str) or not self.project_id.strip():
            raise ValueError("Project identity required")
        if type(self.revision) is not int or self.revision < 0:
            raise ValueError("Nonnegative revision required")
        edits = tuple(self.edits)
        if any(not isinstance(edit, SkyMaskEdit) for edit in edits):
            raise ValueError("Invalid sky edit entry")
        if len({e.edit_id for e in edits}) != len(edits):
            raise ValueError("Duplicate mask edit IDs")
        object.__setattr__(self, "edits", edits)

    def with_edit(self, edit: SkyMaskEdit) -> "SkyMaskDraft":
        if any(e.edit_id == edit.edit_id for e in self.edits):
            raise ValueError("Duplicate edit, remove explicitly first")
        return SkyMaskDraft(self.project_id, self.edits + (edit,), self.revision + 1)

    def without_edit(self, edit_id: str) -> "SkyMaskDraft":
        retained = tuple(e for e in self.edits if e.edit_id != edit_id)
        if len(retained) == len(self.edits):
            raise KeyError(edit_id)
        return SkyMaskDraft(self.project_id, retained, self.revision + 1)

    def to_dict(self) -> dict:
        return {"schema": SCHEMA, "project_id": self.project_id,
                "revision": self.revision, "status": DRAFT_STATUS,
                "edits": [e.to_dict() for e in self.edits]}

    @classmethod
    def from_dict(cls, record: dict) -> "SkyMaskDraft":
        if record.get("schema") != SCHEMA or record.get("status") != DRAFT_STATUS:
            raise ValueError("Invalid/unsupported sky mask draft")
        return cls(record["project_id"],
                   tuple(SkyMaskEdit.from_dict(x) for x in record["edits"]),
                   record["revision"])

    @property
    def digest(self) -> str:
        payload = json.dumps(self.to_dict(), sort_keys=True, allow_nan=False,
                             separators=(",", ":")).encode("utf-8")
        return sha256(payload).hexdigest()


def _polygon_contains(x: np.ndarray, y: np.ndarray, vx: np.ndarray, vy: np.ndarray) -> np.ndarray:
    """Ray-crossing test in the local sky-offset plane."""
    inside = np.zeros_like(x, dtype=bool)
    j = len(vx) - 1
    for i in range(len(vx)):
        yi, yj = vy[i], vy[j]
        cross = ((yi > y) != (yj > y))
        xp = (vx[j] - vx[i]) * (y - yi) / (yj - yi if yi != yj else 1.) + vx[i]
        inside ^= cross & (x < xp)
        j = i
    return inside


def rasterize_sky_shape(shape: SkyMaskShape, wcs, pixel_shape: tuple[int, int],
                        *, subpixels: int = 1, row_chunk: int = 128) -> np.ndarray:
    """Rasterize ICRS geometry onto original FITS (y,x) grid.

    Returns fractional pixel membership; subpixels=1 tests pixel centers.
    No image reprojecting, no mask rewriting, no photometry. The small-field
    polygon/ellipse mapping uses astropy spherical_offsets_to in the tangent
    directions and is not a general spherical polygon construction.
    """
    from astropy.coordinates import SkyCoord
    from astropy import units as u

    if (len(pixel_shape) != 2 or
            any(not isinstance(v, (int, np.integer)) or v <= 0 for v in pixel_shape)):
        raise ValueError("Positive 2D pixel shape required")
    if type(subpixels) is not int or not 1 <= subpixels <= 16:
        raise ValueError("Subpixel factor must be 1–16")
    if not isinstance(row_chunk, int) or row_chunk <= 0:
        raise ValueError("Positive row chunk required")
    ny, nx = map(int, pixel_shape)
    celestial = wcs.celestial
    if shape.kind == "polygon":
        vertices = SkyCoord(
            [v[0] for v in shape.vertices]*u.deg,
            [v[1] for v in shape.vertices]*u.deg, frame="icrs")
        # Vector mean of unit vectors avoids RA-wrap discontinuity.
        vec = vertices.cartesian.xyz.value.mean(axis=1)
        if not np.isfinite(vec).all() or np.linalg.norm(vec) < 1e-9:
            raise ValueError("Polygon has no valid local projection center")
        center = SkyCoord(x=vec[0], y=vec[1], z=vec[2],
                          representation_type="cartesian", frame="icrs")
        de, dn = center.spherical_offsets_to(vertices)
        vx, vy = de.to_value(u.arcsec), dn.to_value(u.arcsec)
        if max(float(np.ptp(vx)), float(np.ptp(vy))) > 3600:
            raise ValueError("Polygon exceeds safe small-field extent")
        if abs(float(np.dot(vx, np.roll(vy, -1)) -
                     np.dot(vy, np.roll(vx, -1)))) < 1e-8:
            raise ValueError("Degenerate polygon")
    else:
        center = SkyCoord(shape.center[0]*u.deg, shape.center[1]*u.deg, frame="icrs")
    out = np.zeros((ny, nx), dtype=np.float64)
    offsets = (np.arange(subpixels) + .5) / subpixels - .5
    for lo in range(0, ny, row_chunk):
        hi = min(ny, lo + row_chunk)
        rows, cols = np.mgrid[lo:hi, 0:nx]
        coverage = np.zeros((hi-lo, nx), dtype=float)
        for oy in offsets:
            for ox in offsets:
                ra, dec = celestial.all_pix2world(cols.astype(float)+ox,
                                                   rows.astype(float)+oy, 0)
                good = np.isfinite(ra) & np.isfinite(dec) & (abs(dec) <= 90)
                if not np.any(good):
                    continue
                locations = SkyCoord(np.mod(ra[good], 360)*u.deg,
                                     dec[good]*u.deg, frame="icrs")
                if shape.kind == "circle":
                    accepted = center.separation(locations).to_value(u.arcsec) < shape.radius_arcsec
                else:
                    east, north = center.spherical_offsets_to(locations)
                    east, north = east.to_value(u.arcsec), north.to_value(u.arcsec)
                    if shape.kind == "ellipse":
                        pa = np.deg2rad(shape.pa_deg)
                        major = east*np.sin(pa) + north*np.cos(pa)
                        minor = east*np.cos(pa) - north*np.sin(pa)
                        accepted = (major / shape.a_arcsec)**2 + (
                            minor / (shape.a_arcsec*shape.q))**2 < 1.
                    else:
                        accepted = _polygon_contains(east, north, vx, vy)
                chunk = np.zeros((hi-lo, nx), dtype=float)
                chunk[good] = accepted
                coverage += chunk
        out[lo:hi] = coverage / subpixels**2
    return out


def evaluate_mask_draft(draft: SkyMaskDraft, wcs, pixel_shape: tuple[int, int],
                        band: str, *, native_hard: np.ndarray,
                        automatic_contaminants: np.ndarray | None = None,
                        subpixels: int = 1) -> dict:
    """Preview only; a restore operation cannot change native-hard exclusions.

    Manual exclusions override restorations. Sky train/holdout are annotations,
    never reasons to 'restore' a native bad pixel. Thresholding fractional
    membership is explicit: >0 = any-overlap conservative rule.
    """
    if band not in BANDS:
        raise ValueError("Unknown band")
    hard = np.asarray(native_hard)
    if hard.shape != tuple(pixel_shape) or hard.dtype.kind != "b":
        raise ValueError("Native hard mask requires same-grid booleans")
    if automatic_contaminants is None:
        auto = np.zeros(pixel_shape, bool)
    else:
        auto = np.asarray(automatic_contaminants)
        if auto.shape != tuple(pixel_shape) or auto.dtype.kind != "b":
            raise ValueError("Automatic contaminant mask requires same-grid booleans")
    layers = {key: np.zeros(pixel_shape, dtype=bool) for key in PURPOSES}
    for edit in draft.edits:
        if not edit.applies_to(band):
            continue
        covered = rasterize_sky_shape(
            edit.shape, wcs, pixel_shape, subpixels=subpixels) > 0
        layers[edit.purpose] |= covered
    if np.any(layers["sky_train"] & layers["sky_holdout"]):
        raise ValueError("Sky training and validation regions must not overlap")
    editable_auto = auto & ~layers["restore_auto"]
    excluded = hard | editable_auto | layers["exclude"]
    return {
        "status": DRAFT_STATUS,
        "project_id": draft.project_id,
        "revision": draft.revision,
        "draft_sha256": draft.digest,
        "band": band,
        "native_hard": hard.copy(),
        "original_auto": auto.copy(),
        "effective_auto": editable_auto,
        "manual_exclude": layers["exclude"],
        "manual_restore_auto": layers["restore_auto"],
        "target_protect": layers["protect_target"],
        "sky_train": layers["sky_train"] & ~excluded & ~layers["protect_target"],
        "sky_holdout": layers["sky_holdout"] & ~excluded & ~layers["protect_target"],
        "excluded": excluded,
        "usable": ~excluded,
    }
