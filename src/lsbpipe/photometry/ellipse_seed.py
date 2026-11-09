"""Preliminary geometry seeds from masked, background-subtracted isophotal regions.

This is a *starting guess*, not an isophote fit, sky estimate, Sérsic
measurement, or final photometry. The supplied residual sky and RMS must
be established upstream; detection thresholds are NOT calibrated sigma
significances after smoothing/correlated noise.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy.ndimage import gaussian_filter, label


SEED_STATUS = "ISOPHOTAL_SEED_NOT_ADOPTED"


@dataclass(frozen=True)
class EllipseSeed:
    x_pix: float
    y_pix: float
    sma_pix: float
    q: float
    theta_native_deg: float
    component_pixels: int
    centroid_offset_pix: float
    state: str = SEED_STATUS
    flags: tuple[str, ...] = ()

    def __post_init__(self):
        numbers = (self.x_pix, self.y_pix, self.sma_pix,
                   self.q, self.theta_native_deg, self.centroid_offset_pix)
        if not all(math.isfinite(v) for v in numbers):
            raise ValueError("Non-finite ellipse seed")
        if not self.sma_pix > 0 or not 0 < self.q <= 1:
            raise ValueError("Invalid ellipse geometry")
        if self.component_pixels < 1 or self.state != SEED_STATUS:
            raise ValueError("Invalid preliminary seed evidence")


def propose_isophotal_seed(image: np.ndarray, rms: np.ndarray | float,
                           valid: np.ndarray, center_xy: tuple[float, float],
                           *, threshold_rms: float = 1.8,
                           smooth_sigma_pix: float = 2.,
                           search_radius_pix: float = 15.,
                           min_pixels: int = 20,
                           center_mode: str = "fixed") -> EllipseSeed:
    """Connected isophotal component + second moments, conditional on a supplied sky.

    The user/astrometric target center anchors selection of the component
    to avoid jumping to the brightest contaminant elsewhere in the image.
    The result is explicitly diagnostic. It is not Jedrzejewski iteration.
    """
    image = np.asarray(image, dtype=float)
    valid = np.asarray(valid)
    if image.ndim != 2 or image.size == 0 or valid.shape != image.shape or valid.dtype.kind != "b":
        raise ValueError("Require same-grid 2D image and Boolean validity mask")
    sigma = np.asarray(rms, dtype=float)
    if sigma.ndim == 0:
        if not np.isfinite(sigma) or sigma <= 0:
            raise ValueError("RMS must be finite and positive")
    elif sigma.shape != image.shape or np.any(~np.isfinite(sigma[valid]) | (sigma[valid] <= 0)):
        raise ValueError("Valid pixels need positive RMS on same grid")
    if (not np.isfinite(image[valid]).all() or
            len(center_xy) != 2 or not all(math.isfinite(float(v)) for v in center_xy)):
        raise ValueError("Invalid source image/center")
    x0, y0 = (float(v) for v in center_xy)
    ny, nx = image.shape
    if not (-.5 <= x0 < nx-.5 and -.5 <= y0 < ny-.5):
        raise ValueError("Initial ICRS-projected center lies outside the image")
    if (not math.isfinite(threshold_rms) or threshold_rms <= 0 or
            not math.isfinite(smooth_sigma_pix) or smooth_sigma_pix < 0 or
            not math.isfinite(search_radius_pix) or search_radius_pix < 0):
        raise ValueError("Invalid threshold/smoothing/search settings")
    if type(min_pixels) is not int or min_pixels < 6:
        raise ValueError("At least six connected pixels required")
    if center_mode not in ("fixed", "moments"):
        raise ValueError("Unknown center mode")

    # Local normalized convolution avoids replacing invalid pixels by a
    # meaningful zero in the smoothing stage. Empty neighborhoods stay invalid.
    good = valid & np.isfinite(image)
    weights = good.astype(float)
    if smooth_sigma_pix:
        numerator = gaussian_filter(np.where(good, image, 0.), smooth_sigma_pix)
        norm = gaussian_filter(weights, smooth_sigma_pix)
        filtered = np.divide(numerator, norm,
                             out=np.zeros_like(numerator), where=norm > 0.5)
        supported = good & (norm > 0.5)
    else:
        filtered = np.where(good, image, 0.)
        supported = good
    if sigma.ndim == 0:
        contrast = filtered / float(sigma)
    else:
        contrast = np.divide(filtered, sigma, out=np.zeros_like(filtered),
                             where=np.asarray(sigma) > 0)
    candidates = supported & (contrast > threshold_rms)
    components, ncomp = label(candidates, structure=np.ones((3, 3), dtype=np.int8))
    if ncomp == 0:
        raise ValueError("No isophotal component above chosen contrast threshold")

    # Consider only components with an actual pixel within the search radius.
    yy, xx = np.indices(image.shape)
    distances = np.hypot(xx - x0, yy - y0)
    vicinity = distances <= search_radius_pix
    viable = []
    for idx in range(1, ncomp + 1):
        region = components == idx
        npx = int(region.sum())
        if npx < min_pixels or not np.any(region & vicinity):
            continue
        rmin = float(np.min(distances[region]))
        viable.append((rmin, -npx, idx))
    if not viable:
        raise ValueError("No sufficiently extended component near target center")
    viable.sort()
    region = components == viable[0][2]
    yy0, xx0 = np.nonzero(region)
    # Positive contrast weights reflect the observed isophotal light,
    # without claiming statistical weighting or a validated flux estimate.
    signal = np.clip(contrast[region] - threshold_rms, 0., None)
    if not np.isfinite(signal).all() or signal.sum() <= 0:
        raise ValueError("No usable positive isophotal contrast")
    cx = float(np.average(xx0, weights=signal))
    cy = float(np.average(yy0, weights=signal))
    dx, dy = xx0 - cx, yy0 - cy
    cov = np.cov(np.vstack((dx, dy)), aweights=signal, bias=True)
    vals, vecs = np.linalg.eigh(cov)
    if vals[0] <= 0 or not np.isfinite(vals).all():
        raise ValueError("Degenerate preliminary isophote")
    major = vecs[:, 1]
    q = float(np.sqrt(vals[0] / vals[1]))
    if q < 0.08:
        raise ValueError("Initial component is too elongated/degenerate")
    theta = float(np.degrees(np.arctan2(major[1], major[0])) % 180.)
    r_area = math.sqrt(len(xx0) / (math.pi*q))
    offset = math.hypot(cx-x0, cy-y0)
    flags = []
    if offset > max(2.5, search_radius_pix/2):
        flags.append("CENTROID_OFFSET_REVIEW")
    if q > .94:
        flags.append("PA_UNDEFINED_NEAR_CIRCULAR")
    if (xx0.min() == 0 or yy0.min() == 0 or
            xx0.max() == nx-1 or yy0.max() == ny-1):
        flags.append("ISOPHOTE_TRUNCATED_AT_IMAGE_EDGE")
    # The threshold is defined in upstream RMS units, but the smoothing
    # introduces different noise statistics and correlations.
    flags.append("SMOOTHED_CONTRAST_NOT_CALIBRATED_SIGNIFICANCE")
    return EllipseSeed(
        x_pix=cx if center_mode == "moments" else x0,
        y_pix=cy if center_mode == "moments" else y0,
        sma_pix=r_area, q=q, theta_native_deg=theta,
        component_pixels=len(xx0), centroid_offset_pix=offset,
        flags=tuple(flags))


def seed_as_sky_geometry(seed: EllipseSeed, wcs) -> dict:
    """Convert native-pixel ellipse by local WCS to ICRS center, PA and axis ratio.

    Uses tangent-plane differential mapping, so a large/strongly distorted
    field needs a separate nonlinear footprint check before science use.
    """
    from astropy.coordinates import SkyCoord
    from astropy import units as u

    center = SkyCoord(*[x*u.deg for x in wcs.celestial.all_pix2world(
        seed.x_pix, seed.y_pix, 0)], frame="icrs")
    offsets = []
    for xx, yy in ((1., 0.), (0., 1.)):
        lon, lat = wcs.celestial.all_pix2world(
            seed.x_pix+xx, seed.y_pix+yy, 0)
        point = SkyCoord(lon*u.deg, lat*u.deg, frame="icrs")
        de, dn = center.spherical_offsets_to(point)
        offsets.append((de.to_value(u.arcsec), dn.to_value(u.arcsec)))
    jac = np.array(offsets, dtype=float).T
    if (not np.isfinite(jac).all() or
            abs(np.linalg.det(jac)) < 1e-10 or np.linalg.cond(jac) > 1e8):
        raise ValueError("Degenerate local celestial WCS")
    angle = math.radians(seed.theta_native_deg)
    pix_basis = np.array([[math.cos(angle), -math.sin(angle)],
                          [math.sin(angle), math.cos(angle)]])
    mapped = jac @ pix_basis @ np.diag([seed.sma_pix, seed.sma_pix*seed.q])
    vectors, singular, _ = np.linalg.svd(mapped)
    pa = float(np.degrees(np.arctan2(vectors[0, 0], vectors[1, 0])) % 180.)
    return {
        "status": seed.state,
        "ra_deg": float(center.ra.deg),
        "dec_deg": float(center.dec.deg),
        "sma_arcsec": float(singular[0]),
        "q_sky": float(singular[1]/singular[0]),
        "pa_north_through_east_deg": pa,
        "source": "LOCAL_WCS_ISOPHOTE_MOMENT_SEED",
        "flags": seed.flags + ("LOCAL_AFFINE_WCS_APPROXIMATION",),
    }
