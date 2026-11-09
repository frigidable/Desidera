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


@dataclass(frozen=True)
class LocalIsophotalProposal:
    """Display-only connected contour and moment ellipse, not science isophotes."""
    seed: EllipseSeed
    contour_xy: tuple[tuple[float, float], ...]
    local_sky: float
    raw_sky_rms: float
    smoothed_sky_rms: float
    sky_pixels: int
    threshold_contrast: float
    smoothing_pix: float
    flags: tuple[str, ...]
    # Sorted flat indices on the *native* pixel grid for reproducible
    # multi-threshold intersection-over-union. Not a science segmentation.
    component_indices: np.ndarray | None = None
    state: str = "LOCAL_ISOPHOTE_PROPOSAL_NOT_ADOPTED"


def _clipped_location(values):
    """Robust location and scale on finite data, no astrophysical sky model."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size < 200:
        raise ValueError("Insufficient valid local background samples")
    for _ in range(4):
        med = float(np.median(values))
        mad = 1.4826 * float(np.median(np.abs(values - med)))
        if mad <= 0 or not math.isfinite(mad):
            raise ValueError("Cannot estimate positive local background dispersion")
        kept = values[np.abs(values - med) <= 3.5 * mad]
        if kept.size < 200:
            raise ValueError("Local background has too few uncontaminated pixels")
        if kept.size == values.size:
            break
        values = kept
    return float(np.median(values)), max(
        1.4826 * float(np.median(np.abs(values - np.median(values)))), 1e-12)


def _outer_boundary(region, xoff, yoff, *, max_vertices=192):
    """Trace the exterior grid-edge ring of a 4-connected pixel component.

    No skimage/matplotlib dependency. Inner holes are deliberately ignored:
    hard-invalid FITS pixels remain invalid independently of region geometry.
    Coordinates are FITS pixel centres, not row-inverted Qt coordinates.
    """
    pixels = np.argwhere(region)
    if len(pixels) < 6:
        raise ValueError("Isophotal component is too small for a contour")
    ny, nx = region.shape
    edges = {}
    for y, x in pixels:
        y, x = int(y), int(x)
        if y == 0 or not region[y-1, x]:
            edges[(x, y)] = (x+1, y)
        if x == nx-1 or not region[y, x+1]:
            edges[(x+1, y)] = (x+1, y+1)
        if y == ny-1 or not region[y+1, x]:
            edges[(x+1, y+1)] = (x, y+1)
        if x == 0 or not region[y, x-1]:
            edges[(x, y+1)] = (x, y)
    rings = []
    while edges:
        first = next(iter(edges))
        cursor = first
        ring = [first]
        for _ in range(len(edges)+1):
            successor = edges.pop(cursor, None)
            if successor is None:
                raise ValueError("Disconnected isophotal pixel boundary")
            if successor == first:
                break
            ring.append(successor)
            cursor = successor
        else:
            raise ValueError("Isophotal boundary could not be closed")
        if len(ring) >= 8:
            vertices = np.asarray(ring, dtype=float)
            area = abs(float(np.dot(vertices[:, 0], np.roll(vertices[:, 1], 1))
                             - np.dot(vertices[:, 1], np.roll(vertices[:, 0], 1))))/2
            rings.append((area, vertices))
    if not rings:
        raise ValueError("No closed isophotal outline")
    exterior = max(rings, key=lambda row: row[0])[1]
    # Reduce stair-step noise and bound immutable sky-mask snapshots to 256
    # vertices. Equal arc-length spacing preserves the contour's spatial
    # character better than taking every Nth grid edge.
    from scipy.ndimage import gaussian_filter1d
    exterior = gaussian_filter1d(exterior, 1.5, axis=0, mode="wrap")
    closed = np.vstack([exterior, exterior[0]])
    distances = np.hypot(np.diff(closed[:, 0]), np.diff(closed[:, 1]))
    perimeter = float(np.sum(distances))
    if perimeter < 8.:
        raise ValueError("Isophotal contour has insufficient perimeter")
    count = min(max_vertices, max(24, int(math.ceil(perimeter / 2.5))))
    lengths = np.concatenate([[0.], np.cumsum(distances)])
    target = np.linspace(0., perimeter, count, endpoint=False)
    xs = np.interp(target, lengths, closed[:, 0]) + xoff - .5
    ys = np.interp(target, lengths, closed[:, 1]) + yoff - .5
    return tuple((float(x), float(y)) for x, y in zip(xs, ys))


def propose_local_isophote(image: np.ndarray, valid: np.ndarray,
                           center_xy: tuple[float, float], *,
                           threshold_contrast: float = 2.2,
                           smooth_sigma_pix: float = 4.,
                           max_radius_pix: float = 380.,
                           bridge_guard_pix: int = 0) -> LocalIsophotalProposal:
    """Propose a TARGET-ANCHORED local threshold isophote, never an adopted one.

    Estimates a constant provisional sky from a distant annulus of the
    original band, measures correlated *smoothed* background fluctuations
    empirically in that same annulus, then finds a 4-connected region near
    the selected ICRS target. This is an exploratory morphological contour,
    NOT a calibrated surface-brightness isophote or reliable noise sigma.
    No source FITS are written and no photometric sums are produced.
    """
    from scipy.ndimage import (
        binary_dilation, binary_opening, gaussian_filter, label,
    )

    data = np.asarray(image, dtype=float)
    accepted = np.asarray(valid)
    if (data.ndim != 2 or data.size == 0 or
            accepted.shape != data.shape or accepted.dtype.kind != "b"):
        raise ValueError("Require native 2D IMAGE and same-grid Boolean validity")
    if len(center_xy) != 2 or not all(math.isfinite(float(p)) for p in center_xy):
        raise ValueError("Invalid proposed galaxy center")
    x0, y0 = map(float, center_xy)
    ny, nx = data.shape
    if not (0 <= x0 < nx and 0 <= y0 < ny):
        raise ValueError("Galaxy centre outside this native FITS cutout")
    if (not 1. <= threshold_contrast <= 8. or
            not 1. <= smooth_sigma_pix <= 12. or
            not 60. <= max_radius_pix <= 750. or
            type(bridge_guard_pix) is not int or
            not 0 <= bridge_guard_pix <= 8):
        raise ValueError("Isophote contrast/smoothing/bridge settings out of range")
    good = accepted & np.isfinite(data)
    yy, xx = np.ogrid[:ny, :nx]
    distance_sq = (xx-x0)**2+(yy-y0)**2
    radius = float(min(nx, ny))
    outer = .47*radius
    inner = .32*radius
    sky_mask = good & (distance_sq >= inner**2) & (distance_sq <= outer**2)
    sky, raw_rms = _clipped_location(data[sky_mask])
    numerator = gaussian_filter(np.where(good, data-sky, 0.), smooth_sigma_pix)
    coverage = gaussian_filter(good.astype(float), smooth_sigma_pix)
    smoothed = np.divide(numerator, coverage,
                         out=np.zeros_like(numerator), where=coverage > .8)
    supported = good & (coverage > .8)
    _, smooth_rms = _clipped_location(smoothed[sky_mask & supported])
    if smooth_rms < raw_rms * .015:
        raise ValueError("Smoothed background is degenerate; no reliable outline")

    # Bounded connected-component discovery avoids scanning millions of
    # isolated stellar speckles. The aperture is centred on the target,
    # not the brightest source in the image.
    box = min(float(max_radius_pix), .42*radius)
    x1, x2 = max(0, int(x0-box)), min(nx, int(x0+box)+2)
    y1, y2 = max(0, int(y0-box)), min(ny, int(y0+box)+2)
    strength = smoothed[y1:y2, x1:x2] / smooth_rms
    positive = supported[y1:y2, x1:x2] & (strength >= threshold_contrast)
    if bridge_guard_pix:
        # A morphology diagnostic: remove links thinner than the selected
        # native-pixel opening radius, never claiming real tidal features
        # are unphysical. The user can turn this guard OFF to compare.
        positive = binary_opening(positive, iterations=bridge_guard_pix)
    connected, count = label(positive)  # 4-connected: diagonal specks do not bridge
    if not count:
        raise ValueError(
            "No connected galaxy-like emission above local contrast. "
            "Lower the contrast threshold or select the diffuse light centre.")
    cy, cx = np.ogrid[y1:y2, x1:x2]
    close = (cx-x0)**2+(cy-y0)**2 <= max(15.,3.*smooth_sigma_pix)**2
    nearby = connected[close]
    nearby = nearby[nearby > 0]
    if nearby.size == 0:
        raise ValueError(
            "No connected emission near selected galaxy centre. "
            "Pick the diffuse component manually.")
    labels, counts = np.unique(nearby, return_counts=True)
    sizes = np.bincount(connected.ravel())
    ranked = sorted(
        (int(i) for i in labels if sizes[i] >= 60),
        key=lambda i: (-int(counts[np.where(labels == i)[0][0]]),
                       -int(sizes[i])))
    if not ranked:
        raise ValueError("No extended galaxy-like component near the target")
    region = connected == ranked[0]
    if (region[0].any() or region[-1].any() or
            region[:, 0].any() or region[:, -1].any()):
        raise ValueError(
            "Proposed region reaches search boundary: enlarge the view "
            "or raise contrast; no clipped contour is offered.")
    pixel_y, pixel_x = np.nonzero(region)
    extent = len(pixel_x)
    signal = np.maximum(strength[region] - threshold_contrast, 0.)
    # Bound bright knots/stellar peaks so they cannot dominate all moments.
    weights = np.minimum(signal, max(float(np.quantile(signal, .85)), 1e-6))
    if not np.isfinite(weights).all() or weights.sum() <= 0:
        raise ValueError("Insufficient diffuse contrast for ellipse moments")
    bary_x = float(np.average(pixel_x+x1, weights=weights))
    bary_y = float(np.average(pixel_y+y1, weights=weights))
    dx, dy = pixel_x+x1-bary_x, pixel_y+y1-bary_y
    covariance = np.cov(np.vstack((dx, dy)), aweights=weights, bias=True)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    if eigenvalues[0] <= 0 or not np.isfinite(eigenvalues).all():
        raise ValueError("Detected region has degenerate geometry")
    q = float(np.sqrt(eigenvalues[0]/eigenvalues[1]))
    if q < .08:
        raise ValueError("Detected region is too thin for a trustworthy ellipse")
    theta = float(np.degrees(np.arctan2(
        eigenvectors[1, 1], eigenvectors[0, 1])) % 180.)
    centroid_offset = float(math.hypot(bary_x-x0,bary_y-y0))
    seed = EllipseSeed(
        x_pix=x0, y_pix=y0, sma_pix=float(math.sqrt(extent/(math.pi*q))),
        q=q, theta_native_deg=theta, component_pixels=extent,
        centroid_offset_pix=centroid_offset,
        flags=("PROVISIONAL_LOCAL_SKY", "SMOOTHED_NOISE_IS_EMPIRICAL_NOT_SIGNIFICANCE"))
    # A small outward display margin is explicitly part of the proposal.
    padded = binary_dilation(region, iterations=2)
    if (padded[0].any() or padded[-1].any() or
            padded[:,0].any() or padded[:,-1].any()):
        raise ValueError("Padded contour reaches search boundary")
    contour = _outer_boundary(padded, x1, y1)
    flags = ["BACKGROUND_CONSTANT_NOT_VALIDATED",
             "NO_VALIDATED_PSFWINGS_OR_MODEL_SUBTRACTION",
             "CONTOUR_PADDING_2_NATIVE_PIXELS"]
    if bridge_guard_pix:
        flags.append(f"THIN_BRIDGE_OPENING_{bridge_guard_pix}PX")
    if centroid_offset > max(12., seed.sma_pix*.30):
        flags.append("CENTROID_OFFSET_REVIEW")
    if q > .94:
        flags.append("PA_UNDEFINED_NEAR_CIRCULAR")
    if (x1 == 0 or y1 == 0 or x2 == nx or y2 == ny):
        flags.append("SEARCH_APERTURE_NEAR_FITS_EDGE")
    component_indices = np.sort(
        (pixel_y.astype(np.int64)+y1)*nx+(pixel_x.astype(np.int64)+x1))
    component_indices.setflags(write=False)
    return LocalIsophotalProposal(
        seed=seed, contour_xy=contour, local_sky=sky,
        raw_sky_rms=raw_rms, smoothed_sky_rms=smooth_rms,
        sky_pixels=int(np.count_nonzero(sky_mask)),
        threshold_contrast=threshold_contrast,
        smoothing_pix=smooth_sigma_pix, flags=tuple(flags),
        component_indices=component_indices)


def assess_local_isophote_sensitivity(
        image: np.ndarray, valid: np.ndarray,
        center_xy: tuple[float, float], reference: LocalIsophotalProposal,
        *, threshold_contrast: float, smooth_sigma_pix: float,
        bridge_guard_pix: int = 0) -> dict:
    """Evaluate multi-threshold/multi-scale robustness without an S/N claim.

    This is an exploratory QA gate. The local sky remains provisional and
    correlated noise is NOT calibrated. Failed nearby thresholds count as
    low confidence, rather than being removed from the denominator.
    """
    variants = (
        (max(1., threshold_contrast-.5), smooth_sigma_pix),
        (min(8., threshold_contrast+.65), smooth_sigma_pix),
        (threshold_contrast, min(12., smooth_sigma_pix*1.35)),
    )
    trials = [{
        "threshold": float(threshold_contrast),
        "smooth_pix": float(smooth_sigma_pix),
        "connected_pixels": reference.seed.component_pixels,
        "q": reference.seed.q,
        "centroid_offset_pix": reference.seed.centroid_offset_pix,
        "overlap_to_reference": 1.,
        "state": "DETECTED",
    }]
    failed = 0
    overlaps = []
    reference_pixels = reference.component_indices
    for contrast, smooth in variants:
        try:
            model = propose_local_isophote(
                image, valid, center_xy,
                threshold_contrast=contrast,
                smooth_sigma_pix=smooth,
                bridge_guard_pix=bridge_guard_pix)
            overlap = None
            if (reference_pixels is not None and
                    model.component_indices is not None):
                shared = np.intersect1d(
                    reference_pixels, model.component_indices,
                    assume_unique=True).size
                union = (reference_pixels.size
                         + model.component_indices.size-shared)
                overlap = float(shared/union) if union > 0 else 0.
                overlaps.append(overlap)
            trials.append({
                "threshold": float(contrast),
                "smooth_pix": float(smooth),
                "connected_pixels": model.seed.component_pixels,
                "q": model.seed.q,
                "centroid_offset_pix": model.seed.centroid_offset_pix,
                "overlap_to_reference": overlap,
                "state": "DETECTED",
            })
        except ValueError as exc:
            failed += 1
            trials.append({
                "threshold": float(contrast),
                "smooth_pix": float(smooth),
                "state": "NO_STABLE_COMPONENT",
                "reason": str(exc)[:240],
            })
    detected = [trial for trial in trials if trial["state"]=="DETECTED"]
    areas = [trial["connected_pixels"] for trial in detected]
    ratios = [trial["q"] for trial in detected]
    offsets = [trial["centroid_offset_pix"] for trial in detected]
    area_ratio = float(max(areas)/min(areas))
    q_span = float(max(ratios)-min(ratios))
    max_offset = float(max(offsets))
    flags = []
    if failed:
        flags.append("CONTOUR_FAILS_NEARBY_THRESHOLDS")
    if area_ratio > 2.2:
        flags.append("CONTOUR_AREA_UNSTABLE")
    if q_span > .25:
        flags.append("ELLIPTICITY_UNSTABLE")
    if overlaps and min(overlaps) < .50:
        flags.append("LOW_SPATIAL_OVERLAP_BETWEEN_ISOPHOTES")
    if max_offset > max(16., .25*reference.seed.sma_pix):
        flags.append("UNANCHORED_CENTROID_DRIFT")
    # No unreviewed automatic adoption even in the best case; this only
    # distinguishes a candidate worth reviewing from an unstable one.
    confidence = ("REVIEWABLE_PRELIMINARY" if not flags
                  else "LOW_CONFIDENCE_REQUIRES_EXPLICIT_OVERRIDE")
    return {
        "confidence": confidence,
        "flags": tuple(flags),
        "trial_count": len(trials),
        "failed_trials": failed,
        "area_ratio": area_ratio,
        "q_spread": q_span,
        "minimum_spatial_iou": min(overlaps) if overlaps else None,
        "max_centroid_offset_pix": max_offset,
        "trials": tuple(trials),
        "sensitivity_status": "MULTISCALE_DIAGNOSTIC_NOT_SCIENCE_VALIDATED",
    }
