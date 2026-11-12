"""P0.5 strict, immutable inputs for future Legacy Survey surface photometry.

This module contains no sky fitting, aperture integration, PSF convolution,
photometric calibration correction, mosaicking or SBF estimator. Values
remain in native coadd units with explicit uncertainty limitations.
"""
from __future__ import annotations

from dataclasses import InitVar, dataclass, field
from hashlib import sha256
import copy
import json
import math
import re
from types import MappingProxyType
from typing import Mapping

import numpy as np

from lsbpipe.survey import legacy_validity, maskbits_to_uint
from lsbpipe.surveys.legacy import BRICK_RE, regions_for_release

CONTRACT_SCHEMA = "lsbpipe.photometry-input/1"
BANDS = frozenset(("g", "r", "i", "z"))
BACKGROUND_STATES = frozenset(("NOT_ESTIMATED", "ESTIMATED"))
PSF_STATES = frozenset(("UNKNOWN", "EMPIRICAL_KERNEL"))
COVARIANCE_STATE = "DIAGONAL_IVAR_CORRELATIONS_UNKNOWN"
AB_FLUX_ZEROPOINT_NMGY = 22.5


def _immutable(data, dtype, *, ndim: int = 2):
    array = np.asarray(data, dtype=dtype)
    if array.ndim != ndim or not array.size:
        raise ValueError("Expected a nonempty 2D raster")
    # A read-only numpy view of a writable buffer could be made mutable
    # again via .setflags(write=True). Own the bytes in a read-only buffer.
    frozen = np.frombuffer(
        np.ascontiguousarray(array).tobytes(), dtype=np.dtype(dtype))
    return frozen.reshape(array.shape)


def _copy_metadata(value):
    """Defensive recursive immutability while retaining JSON-compatible types."""
    encoded = json.dumps(value, sort_keys=True, allow_nan=False)
    def freeze(obj):
        if isinstance(obj, dict):
            return MappingProxyType({str(k): freeze(v) for k, v in obj.items()})
        if isinstance(obj, list):
            return tuple(freeze(v) for v in obj)
        return obj
    return freeze(json.loads(encoded))


def _area_from_wcs(wcs):
    """Local projection-plane area, not exact spherical pixel solid angle."""
    if not getattr(wcs, "has_celestial", False):
        raise ValueError("Image lacks celestial WCS")
    celestial = wcs.celestial
    if celestial.pixel_n_dim != 2 or celestial.world_n_dim != 2:
        raise ValueError("Expected exactly two celestial WCS axes")
    if getattr(celestial, "has_distortion", False):
        raise ValueError(
            "Distorted WCS requires an explicit per-pixel area model; "
            "a central affine area alone is not sufficient")
    matrix = np.asarray(celestial.pixel_scale_matrix, dtype=np.float64)
    if matrix.shape != (2, 2) or not np.all(np.isfinite(matrix)):
        raise ValueError("Invalid WCS pixel scale matrix")
    det = float(np.linalg.det(matrix))
    if not math.isfinite(det) or abs(det) < 1e-20:
        raise ValueError("Singular or numerically unstable celestial WCS")
    # Astropy's celestial WCS normalizes its celestial pixel scale to deg.
    area = abs(det) * (3600.0 ** 2)
    if not math.isfinite(area) or area <= 0:
        raise ValueError("Non-positive local projection-plane pixel area")
    return area


def _check_wcs_samples(wcs, shape):
    ny, nx = shape
    x = np.asarray([0, nx - 1, (nx - 1)/2, 0, nx - 1], dtype=float)
    y = np.asarray([0, 0, (ny - 1)/2, ny - 1, ny - 1], dtype=float)
    ra, dec = wcs.all_pix2world(x, y, 0)
    if not np.all(np.isfinite(ra)) or not np.all(np.isfinite(dec)):
        raise ValueError("WCS gives nonfinite positions inside raster")
    rx, ry = wcs.all_world2pix(ra, dec, 0)
    error = np.hypot(np.asarray(rx) - x, np.asarray(ry) - y)
    if not np.all(np.isfinite(error)) or np.max(error) > 1e-3:
        raise ValueError("WCS inverse roundtrip error exceeds 0.001 pixel")


@dataclass(frozen=True)
class MaskPolicy:
    """Separate hard exclusions from diagnostic Legacy flags."""
    reject_satur: bool = True
    reject_allmask: bool = True
    reject_nprimary: bool = False
    reject_bright: bool = False
    reject_galaxy: bool = False

    def __post_init__(self):
        if any(type(getattr(self, name)) is not bool for name in (
                "reject_satur", "reject_allmask", "reject_nprimary",
                "reject_bright", "reject_galaxy")):
            raise ValueError("Mask flags must be explicit booleans")

    def to_dict(self):
        return {key: getattr(self, key) for key in (
            "reject_satur", "reject_allmask", "reject_nprimary",
            "reject_bright", "reject_galaxy")}


@dataclass(frozen=True)
class PSFDescription:
    """Unknown unless an explicitly supplied normalized empirical kernel exists."""
    state: str = "UNKNOWN"
    kernel: np.ndarray | None = field(default=None, repr=False, compare=False)
    provenance: Mapping = field(default_factory=dict)

    def __post_init__(self):
        if self.state not in PSF_STATES:
            raise ValueError("Unsupported PSF state")
        if self.state == "UNKNOWN" and self.kernel is not None:
            raise ValueError("Unknown PSF must not contain an unverified kernel")
        if self.state == "EMPIRICAL_KERNEL":
            if self.kernel is None:
                raise ValueError("Empirical PSF requires a kernel")
            arr = _immutable(self.kernel, np.float64)
            if min(arr.shape) < 3 or any(d % 2 != 1 for d in arr.shape):
                raise ValueError("Empirical PSF must be odd-sized and at least 3x3")
            if (not np.all(np.isfinite(arr)) or np.any(arr < 0) or
                    not np.isclose(arr.sum(), 1.0, rtol=1e-6, atol=1e-9)):
                raise ValueError("Empirical PSF must be finite, nonnegative and normalized")
            object.__setattr__(self, "kernel", arr)
        object.__setattr__(self, "provenance", _copy_metadata(dict(self.provenance)))


@dataclass(frozen=True)
class SkyDescription:
    """Residual background is unknown until modeled; zero is not assumed."""
    state: str = "NOT_ESTIMATED"
    background_nMgy_per_pixel: np.ndarray | None = field(
        default=None, repr=False, compare=False)
    provenance: Mapping = field(default_factory=dict)

    def __post_init__(self):
        if self.state not in BACKGROUND_STATES:
            raise ValueError("Unsupported background state")
        if self.state == "NOT_ESTIMATED" and self.background_nMgy_per_pixel is not None:
            raise ValueError("A non-estimated background cannot carry a numerical model")
        if self.state == "ESTIMATED":
            if self.background_nMgy_per_pixel is None:
                raise ValueError("Estimated sky requires an explicit raster")
            image = _immutable(self.background_nMgy_per_pixel, np.float64)
            if not np.all(np.isfinite(image)):
                raise ValueError("Sky model must be finite")
            object.__setattr__(self, "background_nMgy_per_pixel", image)
        object.__setattr__(self, "provenance", _copy_metadata(dict(self.provenance)))


@dataclass(frozen=True)
class NativeBandFrame:
    """One native brick, one band, no implicit reprojection or resampling.

    IMAGE is integrated flux per native pixel [nMgy/pixel].
    IVAR describes per-pixel inverse variance [nMgy^-2]; resampling
    covariance is UNKNOWN, not the zero off-diagonal matrix.
    NEXP counts contributing exposures; it is not a noise estimator.
    """
    image_nMgy: np.ndarray = field(repr=False, compare=False)
    invvar_nMgy2: np.ndarray = field(repr=False, compare=False)
    nexp: np.ndarray = field(repr=False, compare=False)
    maskbits: np.ndarray = field(repr=False, compare=False)
    native_wcs: InitVar[object]
    release: str = "dr9"
    region: str = "south"
    brick: str = ""
    band: str = "r"
    policy: MaskPolicy = field(default_factory=MaskPolicy)
    source_evidence: Mapping = field(default_factory=dict, repr=False, compare=False)
    extra_exclusion: np.ndarray | None = field(
        default=None, repr=False, compare=False)
    # Populated only by __post_init__, not caller-controlled.
    _wcs: object = field(init=False, repr=False, compare=False)
    _valid: np.ndarray = field(init=False, repr=False, compare=False)
    _reasons: Mapping = field(init=False, repr=False, compare=False)
    _area: float = field(init=False, repr=False, compare=False)
    _native_input_digest: str | None = field(
        init=False, default=None, repr=False, compare=False)

    def __post_init__(self, native_wcs):
        if self.region not in regions_for_release(self.release):
            raise ValueError("Invalid release/region for native coadd")
        if not isinstance(self.brick, str) or not BRICK_RE.fullmatch(self.brick):
            raise ValueError("Expected 8-character Legacy brick identifier")
        if self.band not in BANDS or (self.region == "north" and self.band == "i"):
            raise ValueError("Optical band unavailable for region")
        if not isinstance(self.policy, MaskPolicy):
            raise ValueError("Explicit MaskPolicy required")
        image = _immutable(self.image_nMgy, np.float64)
        ivar = _immutable(self.invvar_nMgy2, np.float64)
        exposures = _immutable(self.nexp, np.float64)
        bits = _immutable(maskbits_to_uint(self.maskbits), np.uint32)
        if not (image.shape == ivar.shape == exposures.shape == bits.shape):
            raise ValueError("Four native FITS layers require identical 2D shape")
        if not np.all(np.isfinite(exposures) | np.isnan(exposures)):
            raise ValueError("NEXP must contain finite values or NaN")
        area = _area_from_wcs(native_wcs)
        # Only assign a *deep copy* to avoid modifying the input WCS object;
        # use .wcs to obtain another independent copy on every read.
        copied = copy.deepcopy(native_wcs.celestial)
        _check_wcs_samples(copied, image.shape)
        base_valid, reasons = legacy_validity(
            image, ivar, bits, self.band,
            reject_satur=self.policy.reject_satur,
            reject_allmask=self.policy.reject_allmask)
        valid = np.array(base_valid, dtype=bool, copy=True)
        nexp_bad = ~np.isfinite(exposures) | (exposures <= 0)
        reasons = dict(reasons)
        reasons["BAD_NEXP"] = nexp_bad
        valid &= ~nexp_bad
        for label, enabled in (
                ("NPRIMARY", self.policy.reject_nprimary),
                ("BRIGHT", self.policy.reject_bright),
                ("GALAXY", self.policy.reject_galaxy)):
            if enabled:
                valid &= ~reasons[label]
        if self.extra_exclusion is not None:
            extra = _immutable(self.extra_exclusion, np.bool_)
            if extra.shape != image.shape:
                raise ValueError("Additional exclusion mask shape mismatch")
            reasons["USER_EXCLUSION"] = extra
            valid &= ~extra
            object.__setattr__(self, "extra_exclusion", extra)
        object.__setattr__(self, "image_nMgy", image)
        object.__setattr__(self, "invvar_nMgy2", ivar)
        object.__setattr__(self, "nexp", exposures)
        object.__setattr__(self, "maskbits", bits)
        object.__setattr__(self, "_wcs", copied)
        object.__setattr__(self, "_valid", _immutable(valid, np.bool_))
        object.__setattr__(self, "_reasons", MappingProxyType({
            k: _immutable(v, np.bool_) for k, v in reasons.items()}))
        object.__setattr__(self, "_area", area)
        object.__setattr__(self, "source_evidence",
                           _copy_metadata(dict(self.source_evidence)))

    @property
    def dataset_id(self) -> str:
        return f"{self.release}/{self.region}/{self.brick}/{self.band}"

    @property
    def shape(self):
        return self.image_nMgy.shape

    @property
    def wcs(self):
        return copy.deepcopy(self._wcs)

    @property
    def valid(self):
        return self._valid

    @property
    def mask_reasons(self):
        return self._reasons

    @property
    def local_pixel_area_arcsec2(self):
        return self._area

    @property
    def native_input_sha256(self) -> str:
        """SHA-256 of immutable native IMAGE/IVAR/NEXP/MASKBITS and validity.

        This is a local in-memory array identity, NOT an external FITS
        authenticity or WCS certificate. Compute once only when requested
        and hash contiguous buffers without temporary full-array copies.
        """
        cached = self._native_input_digest
        if cached is not None:
            return cached
        digest = sha256(b"desidera-native-input-arrays-v1\n")
        digest.update(self.dataset_id.encode("utf-8"))
        for name in ("image_nMgy", "invvar_nMgy2", "nexp",
                     "maskbits", "_valid"):
            array = getattr(self, name)
            digest.update((f"\n{name}:{array.dtype.str}:"
                           f"{array.shape[0]}x{array.shape[1]}\n").encode("ascii"))
            digest.update(memoryview(array).cast("B"))
        value = digest.hexdigest()
        object.__setattr__(self, "_native_input_digest", value)
        return value

    @property
    def science_state(self):
        return "NATIVE_VALIDATED_NOT_PHOTOMETRY_READY"

    def summary(self) -> dict:
        return {
            "schema": CONTRACT_SCHEMA,
            "dataset_id": self.dataset_id,
            "shape": list(self.shape),
            "units": {
                "image": "nanomaggies/pixel",
                "invvar": "1/nanomaggies^2",
                "nexp": "exposures",
                "maskbits": "bit field",
            },
            "wcs": {
                "pixel_origin": 0,
                "area_arcsec2_local_projection_plane": self._area,
                "area_policy": "LOCAL_WCS_PROJECTION_APPROXIMATE",
                "reprojected": False,
            },
            "mask_policy": self.policy.to_dict(),
            "usable_pixels": int(np.count_nonzero(self.valid)),
            "all_pixels": int(self.valid.size),
            "covariance": COVARIANCE_STATE,
            "upstream_background_state": "UPSTREAM_SKY_SUBTRACTED_RESIDUAL_UNKNOWN",
            "science_state": self.science_state,
        }


@dataclass(frozen=True)
class PhotometryPreparation:
    """Preparation metadata, not permission to run any measurement yet."""
    frame: NativeBandFrame
    sky: SkyDescription = field(default_factory=SkyDescription)
    psf: PSFDescription = field(default_factory=PSFDescription)
    neighbor_model_nMgy_per_pixel: np.ndarray | None = field(
        default=None, repr=False, compare=False)

    def __post_init__(self):
        if not isinstance(self.frame, NativeBandFrame):
            raise ValueError("Photometry preparation requires a native validated frame")
        if not isinstance(self.sky, SkyDescription) or not isinstance(self.psf, PSFDescription):
            raise ValueError("Explicit sky and PSF states required")
        if self.sky.background_nMgy_per_pixel is not None:
            if self.sky.background_nMgy_per_pixel.shape != self.frame.shape:
                raise ValueError("Sky/background raster does not match native image")
        if self.neighbor_model_nMgy_per_pixel is not None:
            neighbors = _immutable(self.neighbor_model_nMgy_per_pixel, np.float64)
            if neighbors.shape != self.frame.shape or not np.all(np.isfinite(neighbors)):
                raise ValueError("Invalid neighbor model raster")
            object.__setattr__(self, "neighbor_model_nMgy_per_pixel", neighbors)

    @property
    def pending(self) -> tuple[str, ...]:
        issues = []
        if self.sky.state != "ESTIMATED":
            issues.append("RESIDUAL_SKY_NOT_ESTIMATED")
        if self.psf.state != "EMPIRICAL_KERNEL":
            issues.append("PSF_NOT_MEASURED")
        issues.append("COADD_NOISE_CORRELATION_NOT_CALIBRATED")
        issues.append("PHOTOMETRIC_ALGORITHM_NOT_IMPLEMENTED")
        return tuple(issues)


@dataclass(frozen=True)
class NativeFrameCollection:
    """Multiple originals stay independent; no implicit multi-brick mosaic."""
    frames: tuple[NativeBandFrame, ...]

    def __post_init__(self):
        frames = tuple(self.frames)
        if not frames or any(not isinstance(f, NativeBandFrame) for f in frames):
            raise ValueError("Collection requires one or more native band frames")
        identities = [f.dataset_id for f in frames]
        if len(set(identities)) != len(identities):
            raise ValueError("Duplicate native dataset in collection")
        object.__setattr__(self, "frames", frames)

    @property
    def processing_groups(self) -> dict[tuple[str, str, str], tuple[str, ...]]:
        groups: dict[tuple[str, str, str], list[str]] = {}
        for frame in self.frames:
            groups.setdefault((frame.release, frame.region, frame.band), []).append(
                frame.dataset_id)
        return {key: tuple(sorted(ids)) for key, ids in sorted(groups.items())}

    @property
    def mosaic_state(self):
        return "NATIVE_UNMOSAICKED_NO_REPROJECTION"

    def require_single_native(self) -> NativeBandFrame:
        if len(self.frames) != 1:
            raise ValueError("Multiple native frames require an explicit mosaic policy")
        return self.frames[0]
