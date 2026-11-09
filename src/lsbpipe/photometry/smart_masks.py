"""Researcher-seeded stellar-contaminant mask candidates on native Legacy FITS.

A researcher selects the source in ICRS and supplies a search radius. This
module proposes a *band-specific* connected isophotal region on validated
original native pixels. No source FITS are altered; proposals require review
and never imply that stellar PSF wings, bleed trails or blended targets
have been adequately masked.

A detection threshold on a smoothed coadd is a *contrast heuristic*, not
a calibrated independent-pixel significance or an optimal noise boundary.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math
import re

import numpy as np
from scipy.ndimage import binary_dilation, gaussian_filter, label, maximum_filter

from lsbpipe.photometry.contracts import NativeBandFrame
from lsbpipe.scene import LayerKind, MaskLayer, grid_fingerprint
from lsbpipe.sky_geometry import ellipse_from_wcs

SMART_MASK_STATUS = "SMART_MASK_CANDIDATE_NOT_ADOPTED"
_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,95}$")


@dataclass(frozen=True)
class StarMaskSeed:
    """User's shared ICRS click/circle, independent of optical passband."""
    source_id: str
    ra_deg: float
    dec_deg: float
    search_radius_arcsec: float

    def __post_init__(self):
        if not isinstance(self.source_id, str) or not _ID_RE.fullmatch(self.source_id):
            raise ValueError("Source identifier must be a safe nonempty label")
        numbers = (self.ra_deg, self.dec_deg, self.search_radius_arcsec)
        if not all(isinstance(x, (int, float)) and math.isfinite(x) for x in numbers):
            raise ValueError("Seed ICRS position and angular radius must be finite")
        if not 0 <= self.ra_deg < 360 or not -90 <= self.dec_deg <= 90:
            raise ValueError("Invalid ICRS star position")
        if not 0 < self.search_radius_arcsec <= 120:
            raise ValueError("Smart-mask pilot supports 0-120 arcsec search regions")


def _freeze_boolean(array: np.ndarray) -> np.ndarray:
    raw = np.asarray(array)
    if raw.dtype.kind != "b" or raw.ndim != 2:
        raise ValueError("Candidate mask must be a 2D Boolean raster")
    return np.frombuffer(raw.tobytes(order="C"), dtype=np.bool_).reshape(raw.shape)


@dataclass(frozen=True)
class SmartStarMaskCandidate:
    """Non-adopted, same-native-WCS raster proposal for a single optical band."""
    source_id: str
    dataset_id: str
    band: str
    grid_id: str
    mask: np.ndarray = field(repr=False, compare=False)
    peak_x_pix: int
    peak_y_pix: int
    n_pixels: int
    threshold_rms: float
    search_radius_arcsec: float
    flags: tuple[str, ...]
    status: str = SMART_MASK_STATUS

    def __post_init__(self):
        if self.status != SMART_MASK_STATUS or self.band not in ("g", "r", "i", "z"):
            raise ValueError("A proposed star mask is not an adopted measurement")
        object.__setattr__(self, "mask", _freeze_boolean(self.mask))
        if int(np.count_nonzero(self.mask)) != self.n_pixels or self.n_pixels < 1:
            raise ValueError("Mask pixel count disagrees with immutable candidate")

    def summary(self) -> dict:
        return {
            "source_id": self.source_id, "dataset_id": self.dataset_id,
            "band": self.band, "grid_id": self.grid_id,
            "peak_fits_pixel_xy": [self.peak_x_pix, self.peak_y_pix],
            "n_masked_valid_pixels": self.n_pixels,
            "threshold_rms": self.threshold_rms,
            "user_search_radius_arcsec": self.search_radius_arcsec,
            "flags": list(self.flags), "status": self.status,
            "mask_is_validated_science": False,
        }

    def to_reviewed_scene_layer(self, *, name: str, confirmed: bool) -> MaskLayer:
        """User confirmation required; native hard-mask rules still apply."""
        if confirmed is not True:
            raise ValueError("Researcher must explicitly confirm this candidate")
        return MaskLayer(
            name=name, kind=LayerKind.MANUAL_EXCLUDE,
            data=self.mask, grid_id=self.grid_id, band=self.band,
            source=f"Researcher-confirmed provisional smart mask: {self.source_id}")


def _region_array(value, shape: tuple[int, int], *, positive: bool, label_: str):
    data = np.asarray(value, dtype=float)
    if data.ndim == 0:
        if not math.isfinite(float(data)) or (positive and data <= 0):
            raise ValueError(f"{label_} must be finite and positive when required")
    elif data.shape != shape:
        raise ValueError(f"{label_} must be scalar or native-frame-shaped")
    return data


def propose_seeded_star_mask(
    frame: NativeBandFrame,
    seed: StarMaskSeed,
    *,
    residual_sky_nMgy_per_pixel: np.ndarray | float,
    local_rms_nMgy_per_pixel: np.ndarray | float,
    target_protection: np.ndarray | None = None,
    threshold_rms: float = 2.5,
    smoothing_sigma_pix: float = 1.0,
    min_pixels: int = 5,
    peak_search_fraction: float = 0.3,
    dilation_pixels: int = 0,
) -> SmartStarMaskCandidate:
    """Find a user-anchored connected star candidate, never automatically adopt.

    FITS are original; residual sky and local RMS must be *supplied*, not
    invented from INVVAR or a default zero. The algorithm only uses a
    bounded patch near the seed, and reprojects the same sky anchor for
    each requested band. A local affine WCS approximation limits use to
    small angular regions and must be independently checked on real FITS.

    This exploratory mode is intended for isolated ordinary stars.
    Bright/saturated stars, PSF wings and blended sources require separate
    PSF-aware models or explicit human-drawn masks.
    """
    if not isinstance(frame, NativeBandFrame) or not isinstance(seed, StarMaskSeed):
        raise TypeError("Requires a validated native photometry frame and sky seed")
    if (not math.isfinite(threshold_rms) or not 0 < threshold_rms <= 30
            or not math.isfinite(smoothing_sigma_pix)
            or not 0 <= smoothing_sigma_pix <= 10
            or type(min_pixels) is not int or min_pixels < 3
            or not math.isfinite(peak_search_fraction)
            or not 0 < peak_search_fraction <= 0.5
            or type(dilation_pixels) is not int or not 0 <= dilation_pixels <= 25):
        raise ValueError("Invalid contrast, smoothing, search or dilation parameters")

    shape = frame.shape
    sky = _region_array(residual_sky_nMgy_per_pixel, shape,
                        positive=False, label_="Residual sky")
    rms = _region_array(local_rms_nMgy_per_pixel, shape,
                        positive=True, label_="Local RMS")
    valid = frame.valid
    if target_protection is not None:
        protect = np.asarray(target_protection)
        if protect.dtype.kind != "b" or protect.shape != shape:
            raise ValueError("Target protection must be native-grid Boolean")
    else:
        protect = None

    # Use the native WCS per band; do not transport pixel coordinates
    # across the Legacy g/r/z/i arrays.
    wcs = frame.wcs
    geometry = ellipse_from_wcs(
        wcs, seed.ra_deg, seed.dec_deg, sky_q=1., sky_pa_deg=0.)
    center = geometry.pixel
    cx, cy = center.x_pix, center.y_pix
    ny, nx = shape
    if not (0 <= cx < nx and 0 <= cy < ny):
        raise ValueError("Clicked ICRS seed lies outside this native image")

    a_pix = seed.search_radius_arcsec * geometry.pixel_a_per_sky_arcsec
    if a_pix < 2:
        raise ValueError("Search radius spans fewer than two pixels")
    # Bounding square of the local projected sky ellipse.
    extent = int(math.ceil(a_pix)) + max(3, dilation_pixels+2)
    left, right = max(0, int(math.floor(cx))-extent), min(nx, int(math.ceil(cx))+extent+1)
    bottom, top = max(0, int(math.floor(cy))-extent), min(ny, int(math.ceil(cy))+extent+1)
    raw = np.asarray(frame.image_nMgy[bottom:top, left:right], dtype=float)
    sky_roi = np.broadcast_to(sky, shape)[bottom:top, left:right]
    rms_roi = np.broadcast_to(rms, shape)[bottom:top, left:right]
    valid_roi = valid[bottom:top, left:right]
    if (np.any(~np.isfinite(sky_roi[valid_roi])) or
            np.any(~np.isfinite(rms_roi[valid_roi])) or
            np.any(rms_roi[valid_roi] <= 0)):
        raise ValueError("Sky and local RMS must be finite and RMS positive on valid pixels")

    yy, xx = np.mgrid[bottom:top, left:right]
    radius_pix = center.elliptical_radius(xx, yy)
    within = radius_pix <= a_pix
    good = valid_roi & within & np.isfinite(raw) & np.isfinite(sky_roi)
    if not np.any(good):
        raise ValueError("No usable science pixels near the clicked source")
    contrast = raw-sky_roi
    if smoothing_sigma_pix > 0:
        weights = good.astype(float)
        numerator = gaussian_filter(np.where(good, contrast, 0.), smoothing_sigma_pix)
        denominator = gaussian_filter(weights, smoothing_sigma_pix)
        filtered = np.divide(numerator, denominator,
                             out=np.zeros_like(numerator), where=denominator > .65)
        supported = good & (denominator > .65)
    else:
        filtered = np.where(good, contrast, 0.)
        supported = good
    statistic = np.divide(filtered, rms_roi,
                          out=np.zeros_like(filtered), where=rms_roi > 0)
    # The user's anchor may be slightly offset. Search nearby, never
    # jump to a brighter galaxy/star elsewhere in the field.
    close = supported & (radius_pix <= a_pix * peak_search_fraction)
    if not np.any(close):
        raise ValueError("No usable seed pixels within allowed center refinement")
    selection = np.where(close, statistic, -np.inf)
    peak_row, peak_col = np.unravel_index(int(np.argmax(selection)), selection.shape)
    if statistic[peak_row, peak_col] <= threshold_rms:
        raise ValueError("No connected source above selected contrast threshold")
    detected = supported & (statistic >= threshold_rms)
    components, _ = label(detected, structure=np.ones((3, 3), dtype=np.uint8))
    component_id = int(components[peak_row, peak_col])
    if component_id == 0:
        raise ValueError("No source component anchored to selected peak")
    component = components == component_id
    n_raw = int(np.count_nonzero(component))
    if n_raw < min_pixels:
        raise ValueError("Seeded source is smaller than minimum connected area")
    if n_raw > 0.7 * np.count_nonzero(within):
        # The click may actually be on the galaxy or a blended region.
        raise ValueError("Seeded component fills most of the user search region")

    warnings = ["SMOOTHED_RMS_CONTRAST_NOT_CALIBRATED_SIGNIFICANCE",
                "EXTENDED_PSF_WINGS_NOT_VALIDATED"]
    if np.any(component & (radius_pix >= .92*a_pix)):
        warnings.append("COMPONENT_TOUCHES_USER_SEARCH_BOUNDARY")
    if np.any(within & ~valid_roi):
        warnings.append("INVALID_NATIVE_PIXELS_NEAR_SOURCE")
    if np.any(component & (np.hypot(xx-(left+peak_col),
                                  yy-(bottom+peak_row)) >
                          .75*a_pix)):
        warnings.append("POSSIBLE_BLEND_OR_ASYMMETRIC_SOURCE")

    # Check for competing peaks within one connected component; a
    # saddle may indicate two blended objects instead of one star.
    prominent = (component &
                 (statistic >= max(threshold_rms, .7*statistic[peak_row, peak_col])) &
                 (statistic == maximum_filter(statistic, size=5)))
    if int(np.count_nonzero(prominent)) > 1:
        warnings.append("MULTIPLE_LOCAL_PEAKS_REVIEW")

    if dilation_pixels:
        component = binary_dilation(component, iterations=dilation_pixels) & within
        component &= valid_roi
        warnings.append("USER_REQUESTED_MORPHOLOGICAL_DILATION")
    if protect is not None and np.any(component & protect[bottom:top, left:right]):
        warnings.append("TARGET_PROTECTION_OVERLAP_REVIEW_REQUIRED")

    output = np.zeros(shape, dtype=bool)
    output[bottom:top, left:right] = component
    return SmartStarMaskCandidate(
        source_id=seed.source_id, dataset_id=frame.dataset_id, band=frame.band,
        grid_id=grid_fingerprint(wcs, shape), mask=output,
        peak_x_pix=left+int(peak_col), peak_y_pix=bottom+int(peak_row),
        n_pixels=int(output.sum()), threshold_rms=float(threshold_rms),
        search_radius_arcsec=float(seed.search_radius_arcsec),
        flags=tuple(warnings))
