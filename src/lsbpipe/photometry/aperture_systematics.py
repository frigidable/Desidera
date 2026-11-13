"""Stage 04.3: read-only conditional systematics and blank-sky diagnostics.

All results describe the CURRENT native coadd, draft masks and fixed aperture
hypotheses. They are not adopted sky corrections, missing-light reconstructions,
or calibrated uncertainties of a total galaxy magnitude.
"""
from __future__ import annotations

from dataclasses import asdict
from hashlib import sha256
import json
import math
from bisect import bisect_right

import numpy as np

from .aperture_flux import (
    ApertureGeometryICRS, ApertureProvenance, ApertureSkyScenario,
    measure_native_aperture, project_icrs_ellipse_to_native,
    _snapshot_native_boolean_mask,
)
from .contracts import NativeBandFrame
from lsbpipe.scene import grid_fingerprint


_SCHEMA = "desidera.photometry.aperture-systematics/1"
_BLANK_SCHEMA = "desidera.photometry.blank-annulus-covariance/3"


def _scene_masks(frame, source_mask, manual_exclude):
    """Snapshot GUI masks once; nested radial probes reuse immutable buffers."""
    arrays = []
    for name, raw in (("source", source_mask), ("manual", manual_exclude)):
        array = np.asarray(raw)
        if array.dtype.kind != "b" or array.shape != frame.shape:
            raise ValueError(name + " mask must be native-grid Boolean")
        arrays.append(_snapshot_native_boolean_mask(
            name, array, frame.shape))
    return arrays[0], arrays[1]


def _ring_aperture(x, y, projected, radii, index):
    """Exact annulus (including 0..r0), translated without resampling."""
    from photutils.aperture import EllipticalAnnulus, EllipticalAperture

    scale = radii[index] / radii[-1]
    if index == 0:
        return EllipticalAperture(
            (x, y), projected.a_pix * scale,
            projected.b_pix * scale,
            theta=projected.theta_rad_pixel_ccw)
    return EllipticalAnnulus(
        (x, y),
        a_in=projected.a_pix * radii[index - 1] / radii[-1],
        a_out=projected.a_pix * scale,
        b_in=projected.b_pix * radii[index - 1] / radii[-1],
        b_out=projected.b_pix * scale,
        theta=projected.theta_rad_pixel_ccw)


def _bounded_quadrant_blank_candidates(shape, margin, step, *, seed=20261112):
    """Predeclare a staggered four-quadrant lattice with O(rows) memory.

    A fully enumerated `(x, y)` list can be millions of entries on a
    large native survey image, defeating the accepted/evaluation budgets.
    Instead, each image row contributes at most TWO contiguous
    arithmetic-progression blocks. A seeded coprime stride visits each
    block's flattened quadrant index exactly once without constructing
    any list of pixel centres. No fluxes or masks enter candidate order.

    Returns (candidate_counts_by_quadrant, candidate_iterator).
    """
    ny, nx = shape
    if not (isinstance(ny, (int, np.integer)) and
            isinstance(nx, (int, np.integer)) and ny > 0 and nx > 0
            and isinstance(margin, int) and margin >= 0
            and isinstance(step, int) and step > 0):
        raise ValueError("Invalid bounded blank candidate geometry")
    blocks = [[], [], [], []]
    totals = [0, 0, 0, 0]
    for row, y in enumerate(range(margin, ny - margin, step)):
        x_first = margin + (step // 2 if row % 2 else 0)
        n_cols = len(range(x_first, nx - margin, step))
        if n_cols == 0:
            continue
        # The two halves are defined solely in original native pixels.
        n_left = max(0, min(
            n_cols, math.ceil((nx / 2. - x_first) / step)))
        choices = (
            (0, x_first, n_left),
            (1, x_first + n_left * step, n_cols - n_left),
        )
        for side, first, count in choices:
            if count == 0:
                continue
            quadrant = side + 2 * int(y >= ny / 2.)
            totals[quadrant] += count
            # prefix end, original y row, first x column
            blocks[quadrant].append(
                (totals[quadrant], int(y), int(first)))
    # Seeded permutation of a quadrant's flattenable integer indices,
    # no permutation array and no sampling with replacement.
    rng = np.random.default_rng(seed)
    offsets, strides, ends = [], [], []
    for quadrant, count in enumerate(totals):
        if count == 0:
            offsets.append(0)
            strides.append(1)
        else:
            offsets.append(int(rng.integers(0, count)))
            stride = int(rng.integers(1, count + 1))
            while math.gcd(stride, count) != 1:
                stride = 1 if stride == count else stride + 1
            strides.append(stride)
        ends.append([record[0] for record in blocks[quadrant]])

    def iterator():
        for i in range(max(totals, default=0)):
            for quadrant, count in enumerate(totals):
                if i >= count:
                    continue
                flat_index = (offsets[quadrant] +
                              i * strides[quadrant]) % count
                block_index = bisect_right(ends[quadrant], flat_index)
                preceding = (ends[quadrant][block_index - 1]
                             if block_index else 0)
                _, y, x_first = blocks[quadrant][block_index]
                yield (quadrant,
                       x_first + step * (flat_index - preceding),
                       y)
    return totals, iterator()



def _translated_exact_aperture_overlap_slices(mask, shape, dx, dy):
    """Reindex an already exact Photutils mask by INTEGER native pixels.

    Fractional pixel overlap is translation invariant on a rectangular
    native pixel lattice. Shifting the *mask bbox* (not its data, WCS,
    science image or fractional weights) exactly preserves that operator.
    Mirrors BoundingBox.get_overlap_slices, including cutout edges.
    """
    if type(dx) is not int or type(dy) is not int:
        raise TypeError("Exact aperture template shifts must be integer pixels")
    ny, nx = shape
    x0, x1 = mask.bbox.ixmin + dx, mask.bbox.ixmax + dx
    y0, y1 = mask.bbox.iymin + dy, mask.bbox.iymax + dy
    if x0 >= nx or x1 <= 0 or y0 >= ny or y1 <= 0:
        return None, None
    xa, xb = max(0, x0), min(nx, x1)
    ya, yb = max(0, y0), min(ny, y1)
    return ((slice(ya, yb), slice(xa, xb)),
            (slice(ya-y0, yb-y0), slice(xa-x0, xb-x0)))


def estimate_blank_annulus_covariance(
        frame: NativeBandFrame, outer: ApertureGeometryICRS,
        radii_arcsec, *, effective_source_mask, manual_exclude,
        min_placements=12, max_placements=None,
        full_target_exclusion_outer: ApertureGeometryICRS | None = None,
        sky: ApertureSkyScenario | None = None):
    """Empirical covariance of signed fluxes in non-overlapping blank rings.

    Candidates have the same native-pixel shape as the galaxy rings; they
    cannot overlap each other or the target aperture. A placement with any
    masked/invalid geometric area is rejected. No assumption of pixel-white
    noise, no sky smoothing, and no inverse covariance are introduced.

    The target exclusion extends only to its specified outer ellipse;
    undetected galaxy outskirts, background gradients and residual
    contaminants remain possible biases. Hence even successful sampling
    remains DIAGNOSTIC and cannot certify a science covariance.
    """
    if not isinstance(frame, NativeBandFrame) or not isinstance(
            outer, ApertureGeometryICRS):
        raise TypeError("Validated native frame and ICRS geometry required")
    sky = ApertureSkyScenario() if sky is None else sky
    if not isinstance(sky, ApertureSkyScenario):
        raise TypeError("Versioned conditional blank-sky scenario required")
    radii = tuple(float(v) for v in radii_arcsec)
    if (not 2 <= len(radii) <= 40 or
            any(not math.isfinite(v) or v <= 0 for v in radii) or
            any(b <= a for a, b in zip(radii[:-1], radii[1:])) or
            not math.isclose(radii[-1], outer.semimajor_arcsec,
                             rel_tol=1e-9, abs_tol=1e-10)):
        raise ValueError("Increasing radii must end at the outer ellipse")
    if (isinstance(min_placements, bool) or
            not isinstance(min_placements, int) or
            min_placements < 3 or min_placements > 256 or
            (max_placements is not None and (
                isinstance(max_placements, bool) or
                not isinstance(max_placements, int) or
                max_placements < min_placements or max_placements > 256))):
        raise ValueError("Invalid bounded blank-aperture sampling limits")

    source, manual = _scene_masks(
        frame, effective_source_mask, manual_exclude)
    invalid = ~frame.valid | source | manual
    projected = project_icrs_ellipse_to_native(outer, frame.wcs)
    # Exact fractional-pixel aperture weights are sensitive to the
    # subpixel phase of the aperture centre. Integer-centred blank
    # samples and a fractionally centred target integrate DIFFERENT
    # linear combinations of the original pixel noise, even for
    # identical axis lengths/angle and a uniform pixel variance.
    # Translate the deterministic native lattice by the target's
    # *fractional* pixel phase (never resample pixel values). This
    # preserves candidate separations and the original sky masks.
    phase_x = float(projected.x_center - math.floor(projected.x_center))
    phase_y = float(projected.y_center - math.floor(projected.y_center))
    # A smaller-scale blank probe must STILL avoid the full galaxy aperture.
    # The centre-distance test is deliberately conservative for ellipses.
    if full_target_exclusion_outer is None:
        full_radius_pix = projected.a_pix
    else:
        if not isinstance(full_target_exclusion_outer, ApertureGeometryICRS):
            raise TypeError("Full target exclusion must be ICRS ellipse")
        full_projection = project_icrs_ellipse_to_native(
            full_target_exclusion_outer, frame.wcs)
        if (math.hypot(
                full_projection.x_center - projected.x_center,
                full_projection.y_center - projected.y_center) > 1e-3):
            raise ValueError("Full target exclusion must share project anchor")
        full_radius_pix = full_projection.a_pix
        if full_radius_pix < projected.a_pix:
            raise ValueError("Full target exclusion cannot be smaller than blank rings")
    # The blank sampling distribution can have a very different
    # diagonal inverse variance from the actual galaxy aperture,
    # especially across NEXP boundaries. Compare the *same fractional
    # pixel operators* before attempting any covariance interpretation.
    # Hard and source/manual masks remain effective at the target.
    target_ring_formal_variances = []
    target_ring_usable_areas = []
    # Store only the used native pixels in each exact ring. Even rings
    # disjoint on the sky share boundary PIXELS after fractional overlap;
    # hence their formal white-noise cross-covariance is generally nonzero.
    # A sparse Gram matrix distinguishes this known integration effect
    # from the unknown correlated noise of the resampled survey coadd.
    target_ring_pixel_indices = []
    target_ring_whitened_weights = []
    # The same fractional-pixel operator is used at EVERY blank location,
    # translated by an integer lattice step after M23's phase matching.
    # Keep the full exact Photutils masks, including their bbox, once:
    # recomputing n_rings × n_candidates expensive exact ellipse-pixel
    # intersections has no scientific benefit and can subtly change
    # boundary rounding from one absolute image position to another.
    exact_ring_templates = []
    for i in range(len(radii)):
        target_mask = _ring_aperture(
            projected.x_center, projected.y_center,
            projected, radii, i).to_mask(method="exact")
        exact_ring_templates.append(target_mask)
        tlarge, tsmall = target_mask.get_overlap_slices(frame.shape)
        if tlarge is None:
            target_ring_formal_variances.append(None)
            target_ring_usable_areas.append(0.)
            target_ring_pixel_indices.append(np.empty(0, dtype=np.int64))
            target_ring_whitened_weights.append(np.empty(0, dtype=float))
            continue
        tw = np.asarray(target_mask.data[tsmall], dtype=float)
        if (not np.isfinite(tw).all() or
                np.any((tw < -1.e-10) | (tw > 1.+1.e-10))):
            raise ValueError("Invalid target-ring exact weights in blank QA")
        tw = np.clip(tw, 0., 1.)
        target_usable = (tw > 0.) & ~invalid[tlarge]
        target_ring_usable_areas.append(float(np.sum(
            tw[target_usable], dtype=np.float64)))
        if np.any(target_usable):
            local_y, local_x = np.nonzero(target_usable)
            pixel_indices = (
                (local_y.astype(np.int64) + tlarge[0].start) *
                frame.shape[1] + local_x + tlarge[1].start)
            whitened_weights = (
                tw[target_usable] /
                np.sqrt(frame.invvar_nMgy2[tlarge][target_usable]))
            target_ring_pixel_indices.append(pixel_indices)
            target_ring_whitened_weights.append(whitened_weights)
            target_ring_formal_variances.append(float(np.sum(
                whitened_weights**2, dtype=np.float64)))
        else:
            target_ring_pixel_indices.append(np.empty(0, dtype=np.int64))
            target_ring_whitened_weights.append(np.empty(0, dtype=float))
            target_ring_formal_variances.append(None)

    # C_ij = sum_p w_ip * w_jp / INVVAR_p on ACTUAL target pixels.
    # A single pixel straddling two exact annuli contributes to BOTH,
    # even for statistically independent native pixel errors. A sparse
    # (n_ring x n_pixel) operator avoids a dense image covariance or
    # allocating one full native image per ring.
    from scipy.sparse import csr_matrix

    ring_sizes = [len(indices) for indices in target_ring_pixel_indices]
    operator = csr_matrix((
        np.concatenate(target_ring_whitened_weights),
        (np.repeat(np.arange(len(radii)), ring_sizes),
         np.concatenate(target_ring_pixel_indices))),
        shape=(len(radii), frame.shape[0] * frame.shape[1]))
    target_diagonal_covariance = (operator @ operator.T).toarray()
    for i, variance in enumerate(target_ring_formal_variances):
        if variance is not None and not np.isclose(
                target_diagonal_covariance[i, i], variance,
                rtol=1.e-11, atol=1.e-10):
            raise ValueError("Inconsistent exact target-ring diagonal covariance")
    measured_ring_indices = [
        i for i, variance in enumerate(target_ring_formal_variances)
        if variance is not None]
    off_diagonal_pairs = [
        (i, j) for i in measured_ring_indices
        for j in measured_ring_indices if i < j and
        target_diagonal_covariance[i, j] > 0.]
    max_boundary_correlation = max((
        float(target_diagonal_covariance[i, j] /
              math.sqrt(target_diagonal_covariance[i, i] *
                        target_diagonal_covariance[j, j]))
        for i, j in off_diagonal_pairs), default=0.)
    target_diagonal_covariance_json = [
        [float(target_diagonal_covariance[i, j])
         if i in measured_ring_indices and j in measured_ring_indices
         else None for j in range(len(radii))]
        for i in range(len(radii))]
    # Cumulative aperture weights sum the underlying exact ring weights.
    # Omitting their shared boundary-pixel covariance would understate
    # even the purely DIAGONAL-INVVAR cumulative measurement variance.
    transfer = np.tril(np.ones((len(radii), len(radii))))
    target_cumulative_diagonal_covariance = (
        transfer @ target_diagonal_covariance @ transfer.T)
    independent_ring_variance_sum = np.cumsum(
        np.diag(target_diagonal_covariance))
    cumulative_valid = [
        any(v is not None for v in target_ring_formal_variances[:i+1])
        for i in range(len(radii))]
    cumulative_formal_variances = [
        float(target_cumulative_diagonal_covariance[i, i])
        if cumulative_valid[i] else None for i in range(len(radii))]
    cumulative_variance_increase = [
        float(target_cumulative_diagonal_covariance[i, i] /
              independent_ring_variance_sum[i] - 1.)
        if cumulative_valid[i] and independent_ring_variance_sum[i] > 0.
        else None for i in range(len(radii))]
    exclusion_distance_pix = full_radius_pix + projected.a_pix + 3.
    margin = int(math.ceil(projected.a_pix + 3))
    step = max(2 * margin + 2, 8)
    ny, nx = frame.shape
    centres, samples, blank_ring_formal_variances = [], [], []
    accepted_quadrants = []
    examined, rejected = 0, 0
    required = max(min_placements, 2 * len(radii) + 2)
    # A heavily masked large CCD must not trigger an unbounded search.
    # For long radial vectors, a hardcoded 48-position default can be
    # LESS than the required 2*n_rings+2. Auto-scale the default
    # sampling budget, but never silently override an explicit user cap.
    requested_max_placements = max_placements
    if max_placements is None:
        max_placements = min(256, max(48, required))
    candidate_budget = max(256, 32 * max_placements)
    # Enumerate at most two integer-progression blocks per native row.
    # Never allocate all possible (x, y) positions: an image can have
    # millions of lattice nodes even when at most 256 blanks are used.
    available_by_quadrant, candidate_positions = (
        _bounded_quadrant_blank_candidates(frame.shape, margin, step))
    accepted_by_quadrant = [0, 0, 0, 0]
    translated_ring_evaluations = 0
    target_integer_x = math.floor(projected.x_center)
    target_integer_y = math.floor(projected.y_center)
    for quadrant, lattice_x, lattice_y in candidate_positions:
        if len(centres) >= max_placements or examined >= candidate_budget:
            break
        examined += 1
        x = float(lattice_x + phase_x)
        y = float(lattice_y + phase_y)
        # At least the full outer aperture separates each sky placement
        # from the putative target. Extended LSB light is NOT excluded.
        if math.hypot(
                x - projected.x_center,
                y - projected.y_center) <= exclusion_distance_pix:
            rejected += 1
            continue
        vector = []
        formal_vector = []
        # x - projected.x_center and y - projected.y_center are exact
        # integer translations by CONSTRUCTION; never round a measured
        # floating-point WCS displacement to decide how to shift.
        shift_x = int(lattice_x - target_integer_x)
        shift_y = int(lattice_y - target_integer_y)
        for i, mask in enumerate(exact_ring_templates):
            large, small = _translated_exact_aperture_overlap_slices(
                mask, frame.shape, shift_x, shift_y)
            if large is None:
                vector = []
                break
            translated_ring_evaluations += 1
            weights = np.asarray(mask.data[small], dtype=float)
            if (not np.isfinite(weights).all() or
                    np.any((weights < -1e-10) | (weights > 1 + 1e-10))):
                raise ValueError("Invalid blank-aperture exact weights")
            weights = np.clip(weights, 0., 1.)
            # This checks all partially covered boundary pixels too.
            # Reject every positive-area contact with a bad pixel. A
            # tiny tolerated masked fraction is not an unmasked blank.
            # Crucially, even a zero fractional weight times NaN is NaN
            # in NumPy: never sum the entire rectangular mask bounding
            # box, which may contain invalid pixels OUTSIDE the ring.
            support = weights > 0.
            if (not np.any(support) or
                    np.any(invalid[large] & support)):
                vector = []
                break
            pixel_values = frame.image_nMgy[large][support]
            if sky.method != "zero_additional_correction":
                # The target aperture uses the declared sky plane.
                # A blank-ring covariance measured on UNcorrected
                # coadd pixels cannot represent residuals of that
                # same conditional photometric measurement.
                ry, rx = np.nonzero(support)
                xref, yref = sky.reference_xy_pix
                pixel_values = pixel_values - (
                    sky.b0_nMgy_per_pixel +
                    sky.bx_nMgy_per_pixel2 *
                    (rx+large[1].start-xref) +
                    sky.by_nMgy_per_pixel2 *
                    (ry+large[0].start-yref))
            vector.append(float(np.sum(
                weights[support] * pixel_values,
                dtype=np.float64)))
            diagonal_variance = float(np.sum(
                weights[support]**2 /
                frame.invvar_nMgy2[large][support], dtype=np.float64))
            if not np.isfinite(diagonal_variance):
                raise ValueError("Nonfinite blank native diagonal variance")
            formal_vector.append(diagonal_variance)
        if len(vector) != len(radii):
            rejected += 1
            continue
        centres.append([float(x), float(y)])
        # Quadrant-stratified lattice labels refer to the INTEGER
        # candidate node. A fractional phase shift can cross a quadrant
        # boundary: use the actual placed centre for the scientific
        # spatial-mean contrast (but preserve nominal lattice counts).
        actual_quadrant = int(x >= nx/2.) + 2*int(y >= ny/2.)
        accepted_quadrants.append(actual_quadrant)
        accepted_by_quadrant[quadrant] += 1
        samples.append(vector)
        blank_ring_formal_variances.append(formal_vector)

    if blank_ring_formal_variances:
        blank_formal_array = np.asarray(
            blank_ring_formal_variances, dtype=float)
        median_blank_formal = np.median(
            blank_formal_array, axis=0).tolist()
        target_to_blank_formal = [
            (float(t/b) if t is not None and b > 0. else None)
            for t, b in zip(target_ring_formal_variances,
                            median_blank_formal)]
        blank_formal_min = np.min(
            blank_formal_array, axis=0).tolist()
        blank_formal_max = np.max(
            blank_formal_array, axis=0).tolist()
    else:
        median_blank_formal = [None]*len(radii)
        target_to_blank_formal = [None]*len(radii)
        blank_formal_min = [None]*len(radii)
        blank_formal_max = [None]*len(radii)

    physical_quadrant_counts = [
        int(accepted_quadrants.count(q)) for q in range(4)]
    report = {
        "schema": _BLANK_SCHEMA,
        "dataset_id": frame.dataset_id,
        "radii_arcsec": list(radii),
        "conditional_sky_scenario": asdict(sky),
        "blank_signal_uses_same_frozen_sky_plane_as_target": True,
        "sky_plane_not_calibrated_or_adopted": True,
        "placement_policy": (
            "STRATIFIED_QUADRANT_BOUNDED_MEMORY_NATIVE_RINGS_V3"),
        "candidate_order_independent_of_image_brightness": True,
        "candidate_order_random_seed": 20261112,
        "blank_centre_phase_matches_target_native_xy": True,
        "target_native_subpixel_phase_xy": [phase_x, phase_y],
        "blank_lattice_translation_is_integer_pixels_after_phase_match": True,
        "blank_ring_operator_reuse_policy":
            "ONE_EXACT_NATIVE_RING_TEMPLATE_PER_RADIUS_INTEGER_BBOX_TRANSLATION_V1",
        "n_exact_ring_templates_rasterized": len(exact_ring_templates),
        "n_integer_translated_ring_evaluations": translated_ring_evaluations,
        "blank_ring_weights_not_interpolated_or_approximated": True,
        "no_science_covariance_calibration_from_phase_matching": True,
        "native_diagonal_variance_comparability": {
            "schema": "desidera.photometry.blank-vs-target-diagonal-variance/1",
            "target_partial_ring_usable_area_native_pixel2":
                target_ring_usable_areas,
            "target_partial_ring_formal_variance_nMgy2":
                target_ring_formal_variances,
            "median_full_blank_ring_formal_variance_nMgy2":
                median_blank_formal,
            "min_full_blank_ring_formal_variance_nMgy2":
                blank_formal_min,
            "max_full_blank_ring_formal_variance_nMgy2":
                blank_formal_max,
            "target_to_median_blank_diagonal_variance_ratio":
                target_to_blank_formal,
            "blank_measurements_use_complete_valid_rings": True,
            "target_measurements_keep_original_source_and_manual_exclusions":
                True,
            "native_ivar_is_diagonal_only_not_coadd_covariance": True,
            "target_exact_ring_diagonal_pixel_covariance": {
                "schema": "desidera.photometry.exact-target-ring-diagonal-pixel-covariance/1",
                "covariance_nMgy2": target_diagonal_covariance_json,
                "n_nonzero_off_diagonal_pairs": len(off_diagonal_pairs),
                "max_off_diagonal_correlation": max_boundary_correlation,
                "cumulative_diagonal_pixel_formal_variance_nMgy2":
                    cumulative_formal_variances,
                "cumulative_variance_increase_vs_independent_ring_sum":
                    cumulative_variance_increase,
                "cumulative_covariance_uses_exact_ring_transfer": True,
                "null_row_and_column_means_no_usable_target_pixels": True,
                "exact_fractional_boundaries_can_share_native_pixels": True,
                "off_diagonal_terms_are_present_even_for_white_pixel_noise": True,
                "excludes_resampling_noise_correlations_and_sky_uncertainty": True,
                "no_science_covariance_calibration": True,
            },
            "ratio_not_a_science_covariance_rescaling_factor": True,
            "no_automatic_target_noise_correction": True,
        },
        "candidate_storage_policy": "ROW_BLOCKS_AND_COPRIME_INDEX_STRIDE_NO_CENTER_ARRAY",
        "candidate_grid_nodes_by_quadrant": available_by_quadrant,
        # A phase shift can cross the half-image boundary: spatial
        # coverage QA must refer to ACTUALLY PLACED native centres,
        # not the quadrant originally assigned to each integer lattice node.
        "accepted_centres_by_quadrant": physical_quadrant_counts,
        "accepted_lattice_nodes_by_quadrant": accepted_by_quadrant,
        "accepted_centres_are_actual_phase_shifted_native_positions": True,
        "n_quadrants_with_accepted_blanks": sum(
            n > 0 for n in physical_quadrant_counts),
        "n_examined": int(examined),
        "n_rejected": int(rejected),
        "n_accepted": len(centres),
        "n_required": int(required),
        "n_max_placements": int(max_placements),
        "placement_cap_explicit": requested_max_placements is not None,
        "n_candidate_budget": int(candidate_budget),
        "target_exclusion_radius_native_pix": float(full_radius_pix),
        "blank_outer_radius_native_pix": float(projected.a_pix),
        "target_plus_blank_exclusion_distance_native_pix": float(
            exclusion_distance_pix),
        "candidate_search_budget_exhausted": examined >= candidate_budget,
        "centres_native_xy": centres,
        "annulus_mean_blank_flux_nMgy": None,
        "quadrant_mean_flux_scatter_decomposition": {
            "schema": "desidera.photometry.blank-spatial-mean-scatter/1",
            "status": "INSUFFICIENT_VALID_BLANK_COVARIANCE",
            "counts_by_native_quadrant": physical_quadrant_counts,
            "lattice_stratification_quadrant_counts":
                accepted_by_quadrant.copy(),
            "quadrants_assigned_by_actual_phase_shifted_centres": True,
            "per_quadrant_mean_signed_ring_flux_nMgy": {},
            "within_quadrant_flux_sum_squares_nMgy2": None,
            "between_quadrant_mean_flux_sum_squares_nMgy2": None,
            "total_flux_sum_squares_nMgy2": None,
            "between_quadrant_mean_fraction_of_total_sum_squares": None,
            "same_coadd_quadrants_are_not_independent_exposures": True,
            "between_quadrant_contrast_is_not_proven_sky_gradient": True,
            "not_a_noise_covariance_debiasing_or_science_gate": True,
        },
        "annulus_covariance_nMgy2": None,
        "annulus_sigma_empirical_nMgy": None,
        "cumulative_covariance_nMgy2": None,
        "cumulative_sigma_empirical_nMgy": None,
        "status": "INSUFFICIENT_INDEPENDENT_BLANK_POSITIONS",
        "covariance_calibrated_for_science": False,
        "warnings": [
            "TARGET_LIGHT_OUTSIDE_TESTED_ELLIPSE_NOT_EXCLUDED",
            "UNDETECTED_SOURCES_OR_BACKGROUND_GRADIENTS_MAY_REMAIN",
            "BLANK_POSITIONS_FROM_ONE_COADD_NOT_AN_INDEPENDENT_EXPOSURE_TEST",
            "NONOVERLAPPING_PLACEMENTS_MAY_STILL_SHARE_LARGE_SCALE_NOISE",
        ],
    }
    # Quasi-independent positions across a single coadd are not
    # automatically representative. Capture that failure explicitly
    # even if the nominal sample covariance matrix is finite.
    if len(centres) > 0:
        if report["n_quadrants_with_accepted_blanks"] < 3:
            report["warnings"].append(
                "BLANK_APERTURE_SPATIAL_QUADRANT_COVERAGE_LIMITED")
        if max(physical_quadrant_counts) / len(centres) > .5:
            report["warnings"].append(
                "BLANK_APERTURE_SAMPLES_DOMINATED_BY_ONE_QUADRANT")
    if required > max_placements:
        report["status"] = "EXPLICIT_PLACEMENT_CAP_BELOW_REQUIRED_COVARIANCE_SAMPLE"
        report["warnings"].append(
            "REQUESTED_CAP_CANNOT_MEET_MINIMUM_SAMPLES_FOR_RADIAL_DIMENSION")
    if len(centres) >= required:
        observations = np.asarray(samples, dtype=float)
        covariance = np.atleast_2d(np.cov(observations, rowvar=False, ddof=1))
        diagonal = np.diag(covariance)
        if (np.isfinite(covariance).all() and
                np.all(diagonal >= 0) and np.any(diagonal > 0)):
            report["status"] = "EMPIRICAL_BLANK_COVARIANCE_DIAGNOSTIC_ONLY"
            report["annulus_mean_blank_flux_nMgy"] = (
                np.mean(observations, axis=0).tolist())
            report["annulus_covariance_nMgy2"] = covariance.tolist()
            report["annulus_sigma_empirical_nMgy"] = np.sqrt(
                diagonal).tolist()
            # The sample covariance conflates within-location noise
            # with coherent differences in means across quadrants.
            # Decompose the exact SAME signed blank observations by
            # their original, geometrically predeclared native quadrant,
            # without treating each quadrant as an independent exposure.
            labels = np.asarray(accepted_quadrants, dtype=int)
            n_quadrants = int(np.count_nonzero(
                physical_quadrant_counts))
            qa = report["quadrant_mean_flux_scatter_decomposition"]
            if n_quadrants >= 2:
                global_mean = np.mean(observations, axis=0)
                total_ss = np.sum(
                    (observations-global_mean)**2, axis=0)
                within_ss = np.zeros(len(radii), dtype=float)
                between_ss = np.zeros(len(radii), dtype=float)
                per_quadrant = {}
                for q in range(4):
                    part = observations[labels == q]
                    if not len(part):
                        continue
                    local_mean = np.mean(part, axis=0)
                    per_quadrant[str(q)] = local_mean.tolist()
                    within_ss += np.sum((part-local_mean)**2, axis=0)
                    between_ss += len(part)*(local_mean-global_mean)**2
                if not np.allclose(
                        total_ss, within_ss+between_ss,
                        rtol=1.e-10, atol=1.e-8):
                    raise ValueError(
                        "Inconsistent native blank spatial scatter decomposition")
                fractions = [
                    float(b/t) if t > 0. else None
                    for b,t in zip(between_ss,total_ss)]
                qa.update({
                    "status": "SPATIAL_MEAN_CONTRAST_DESCRIPTIVE_ONLY",
                    "per_quadrant_mean_signed_ring_flux_nMgy": per_quadrant,
                    "within_quadrant_flux_sum_squares_nMgy2":
                        within_ss.tolist(),
                    "between_quadrant_mean_flux_sum_squares_nMgy2":
                        between_ss.tolist(),
                    "total_flux_sum_squares_nMgy2": total_ss.tolist(),
                    "between_quadrant_mean_fraction_of_total_sum_squares":
                        fractions,
                    "exact_sum_of_squares_identity_verified": True,
                })
            else:
                qa["status"] = (
                    "ONE_SPATIAL_QUADRANT_NO_BETWEEN_QUADRANT_CONTRAST")
            # Nested apertures share their inner rings: C_CoG = T C_ring T^T.
            # No independent-error quadrature and no covariance inversion.
            transfer = np.tril(np.ones((len(radii), len(radii))))
            cumulative = transfer @ covariance @ transfer.T
            report["cumulative_covariance_nMgy2"] = cumulative.tolist()
            report["cumulative_sigma_empirical_nMgy"] = np.sqrt(
                np.maximum(0., np.diag(cumulative))).tolist()
        else:
            report["status"] = "DEGENERATE_OR_NONFINITE_BLANK_VARIANCE"
    json.dumps(report, allow_nan=False)
    return report


def _smaller_blank_scales(frame, full_outer, full_radii, source, manual, *,
                          full_scale_available, sky=None):
    """Fixed prefix tests, never represent them as full-size covariance.

    A 105-arcsec blank ellipse cannot be tiled in a 4-arcmin cutout. Smaller
    apertures can still constrain native-coadd fluctuations at shorter scales
    IF they remain outside the *full* target footprint and are clean.
    """
    result = {
        "status": "NOT_REQUIRED_FULL_SCALE_DIAGNOSTIC_AVAILABLE"
        if full_scale_available else "SMALLER_SCALES_NOT_YET_ESTABLISHED",
        "full_target_semimajor_arcsec": full_outer.semimajor_arcsec,
        "full_scale_covariance_available": bool(full_scale_available),
        "probes": [],
        "first_available_probe_outer_arcsec": None,
        "must_not_replace_full_scale_covariance": True,
        "no_target_mask_or_fits_modified": True,
        "covariance_calibrated_for_science": False,
        "adopted_for_science": False,
    }
    if full_scale_available or len(full_radii) < 3:
        return result
    # Predeclared, decreasing radii: no data-dependent sampling selection.
    indices = []
    for value in (
            len(full_radii) - 4, len(full_radii) // 2 - 1, 2, 1):
        if 1 <= value < len(full_radii) - 1 and value not in indices:
            indices.append(value)
    indices.sort(reverse=True)
    for index in indices:
        radius = float(full_radii[index])
        small_outer = ApertureGeometryICRS(
            ra_deg=full_outer.ra_deg, dec_deg=full_outer.dec_deg,
            semimajor_arcsec=radius, axis_ratio=full_outer.axis_ratio,
            pa_deg_east_of_north=full_outer.pa_deg_east_of_north)
        probe = estimate_blank_annulus_covariance(
            frame, small_outer, full_radii[:index+1],
            effective_source_mask=source, manual_exclude=manual,
            full_target_exclusion_outer=full_outer, sky=sky)
        result["probes"].append(probe)
        if probe["annulus_covariance_nMgy2"] is not None:
            result["status"] = "SMALL_SCALE_COVARIANCE_DIAGNOSTIC_ONLY"
            result["first_available_probe_outer_arcsec"] = radius
            break
    if result["first_available_probe_outer_arcsec"] is None:
        result["status"] = "NO_CLEAN_SMALL_SCALE_COVARIANCE_AVAILABLE"
    return result



def _verify_frozen_growth_curves(frame, anchor_icrs, provenance,
                                 source, manual, curves, sky):
    """Fail closed on stale Stage 04.2 inputs, including same-area mask edits.

    Comparing only valid pixel area cannot identify changed source pixels,
    image values, inverse variances or a different conditional sky model.
    Verify each frozen receipt and recompute its *outermost* direct aperture
    against the original native frame and the exact current mask/sky recipe.
    This is bounded by the number of declared geometry hypotheses; no
    photometry, FITS or user draft is written.
    """
    if not isinstance(curves, (tuple, list)) or not curves:
        raise ValueError("At least one frozen Stage 04.2 curve required")
    expected_grid = grid_fingerprint(frame.wcs, frame.shape)
    expected_native_input = frame.native_input_sha256
    reference_radii = tuple(curves[0]["radii_arcsec"])
    anchor = tuple(map(float, anchor_icrs))
    for curve in curves:
        if not isinstance(curve, dict):
            raise ValueError("Stage 04.2 curve must be a frozen QA mapping")
        receipt = sha256(json.dumps(
            {key: value for key, value in curve.items()
             if key != "receipt_sha256"},
            sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode("utf-8")).hexdigest()
        if curve.get("receipt_sha256") != receipt:
            raise ValueError("Stage 04.2 growth receipt integrity mismatch")
        expected = {
            "project_id": provenance.project_id,
            "dataset_id": frame.dataset_id,
            "grid_fingerprint": expected_grid,
            "native_input_sha256": expected_native_input,
            "draft_sha256": provenance.scene_draft_sha256,
            "catalogue_sha256": provenance.catalogue_sha256,
            "source_policy_id": provenance.source_policy_id,
            "sky_method": sky.method,
            "sky_evidence_id": sky.evidence_id,
        }
        for key, value in expected.items():
            if curve.get(key) != value:
                raise ValueError(
                    "Stage 04.2 lineage or sky scenario mismatch: " + key)
        if (tuple(curve.get("fixed_project_anchor_icrs", ())) != anchor or
                tuple(curve.get("radii_arcsec", ())) != reference_radii):
            raise ValueError("Stage 04.2 anchor or radial grid mismatch")
        final = curve["cumulative"][-1]
        hypothesis = curve["geometry_hypothesis"]
        outer = ApertureGeometryICRS(
            ra_deg=anchor[0], dec_deg=anchor[1],
            semimajor_arcsec=float(reference_radii[-1]),
            axis_ratio=float(hypothesis["q"]),
            pa_deg_east_of_north=float(hypothesis["pa_deg"]))
        fresh = measure_native_aperture(
            frame, outer, provenance,
            effective_source_mask=source, manual_exclude=manual, sky=sky)
        if fresh.recipe_sha256 != final.get("recipe_sha256"):
            raise ValueError(
                "Stage 04.2 aperture recipe mismatch: mask pixels, "
                "geometry, calibration lineage or sky changed")
        for name, recorded, current in (
            ("valid area", final["area_native_pixel2"]["valid_native_pixel2"],
             fresh.area_native_pixel2["valid_native_pixel2"]),
            ("signed flux", final["signed_valid_flux_nMgy"],
             fresh.signed_valid_flux_nMgy),
            ("formal variance", final["formal_diagonal_sigma_nMgy"],
             fresh.formal_diagonal_sigma_nMgy),
        ):
            if recorded is None and current is None:
                continue
            if (recorded is None or current is None or
                    not math.isclose(float(recorded), float(current),
                                     rel_tol=1e-11, abs_tol=1e-8)):
                raise ValueError(
                    "Stage 04.2 original native pixels or weights changed: "
                    + name)

def diagnose_aperture_systematics(
        frame: NativeBandFrame, *, anchor_icrs, provenance: ApertureProvenance,
        effective_source_mask, manual_exclude, curves,
        sky: ApertureSkyScenario | None = None):
    """Conditional sensitivity at the largest common, observed aperture.

    - sky response = exact derivative of observed valid-area flux to a
      residual constant or tilted native-pixel sky at the target centre;
    - source-mask counterfactual = recompute with ONLY provisional catalogue
      exclusions disabled; never restores hard-invalid or manual-excluded data;
    - geometry span = difference of distinct observed regions, NOT 1-sigma;
    - blank covariance = non-overlapping empty-field diagnostic, if feasible.
    """
    if not isinstance(frame, NativeBandFrame) or not isinstance(
            provenance, ApertureProvenance):
        raise TypeError("Validated frame and project provenance required")
    source, manual = _scene_masks(
        frame, effective_source_mask, manual_exclude)
    sky = ApertureSkyScenario() if sky is None else sky
    if not isinstance(sky, ApertureSkyScenario):
        raise TypeError("Versioned frozen sky scenario required")
    _verify_frozen_growth_curves(
        frame, anchor_icrs, provenance, source, manual, curves, sky)
    first = curves[0]
    last = first["cumulative"][-1]
    radius = float(first["radii_arcsec"][-1])
    hypothesis = first["geometry_hypothesis"]
    outer = ApertureGeometryICRS(
        ra_deg=float(anchor_icrs[0]), dec_deg=float(anchor_icrs[1]),
        semimajor_arcsec=radius, axis_ratio=float(hypothesis["q"]),
        pa_deg_east_of_north=float(hypothesis["pa_deg"]))
    projected = project_icrs_ellipse_to_native(outer, frame.wcs)

    # Apply the same exact native-pixel integration operator to EVERY
    # frozen radius, not only the outermost aperture.  The response is
    # conditional on the fixed masks, WCS and sky pivot, and its spatial
    # moments are derivatives, not an assumed uncertainty.
    from photutils.aperture import EllipticalAperture
    xref, yref = sky.reference_xy_pix
    radial_rows = []
    for radial_index, radial_radius in enumerate(first["radii_arcsec"]):
        radial_geom = ApertureGeometryICRS(
            ra_deg=outer.ra_deg, dec_deg=outer.dec_deg,
            semimajor_arcsec=float(radial_radius),
            axis_ratio=outer.axis_ratio,
            pa_deg_east_of_north=outer.pa_deg_east_of_north)
        radial_projected = project_icrs_ellipse_to_native(
            radial_geom, frame.wcs)
        aperture = EllipticalAperture(
            (radial_projected.x_center, radial_projected.y_center),
            radial_projected.a_pix, radial_projected.b_pix,
            theta=radial_projected.theta_rad_pixel_ccw)
        mask = aperture.to_mask(method="exact")
        large, small = mask.get_overlap_slices(frame.shape)
        area_i = moment_x = moment_y = 0.
        direct_flux = direct_sigma = None
        source_only_area = source_only_flux = 0.
        if large is not None:
            weights = np.asarray(mask.data[small], dtype=float)
            if (not np.isfinite(weights).all() or
                    np.any((weights < -1e-10) |
                           (weights > 1.+1e-10))):
                raise ValueError(
                    "Invalid Stage 04.3 exact aperture overlap fractions")
            weights = np.clip(weights, 0., 1.)
            usable = frame.valid[large] & ~source[large] & ~manual[large]
            yy, xx = np.nonzero(usable & (weights > 0.))
            if xx.size:
                wt = weights[yy, xx]
                area_i = float(np.sum(wt, dtype=np.float64))
                moment_x = float(np.sum(
                    wt*(xx+large[1].start-xref), dtype=np.float64))
                moment_y = float(np.sum(
                    wt*(yy+large[0].start-yref), dtype=np.float64))
                pixel_values = frame.image_nMgy[large][yy, xx]
                inverse_variance = frame.invvar_nMgy2[large][yy, xx]
                fitted_sky = (
                    sky.b0_nMgy_per_pixel +
                    sky.bx_nMgy_per_pixel2*(xx+large[1].start-xref) +
                    sky.by_nMgy_per_pixel2*(yy+large[0].start-yref))
                direct_flux = float(np.sum(
                    wt*(pixel_values-fitted_sky), dtype=np.float64))
                direct_sigma = float(np.sqrt(np.sum(
                    wt*wt/inverse_variance, dtype=np.float64)))
            # In the alternative treatment, ONLY provisional catalogue
            # exclusions are removed. Never unmask hard-invalid or
            # manual-excluded pixels, nor mistake their flux for proven
            # target light. Compute the exact native-pixel flux increment
            # directly; this avoids another full aperture run per radius.
            extra = (frame.valid[large] & source[large] &
                     ~manual[large] & (weights > 0.))
            ey, ex = np.nonzero(extra)
            if ex.size:
                ew = weights[ey, ex]
                source_only_area = float(np.sum(ew, dtype=np.float64))
                sky_at_pixels = (
                    sky.b0_nMgy_per_pixel +
                    sky.bx_nMgy_per_pixel2 *
                    (ex+large[1].start-xref) +
                    sky.by_nMgy_per_pixel2 *
                    (ey+large[0].start-yref))
                source_only_flux = float(np.sum(
                    ew*(frame.image_nMgy[large][ey, ex] -
                        sky_at_pixels), dtype=np.float64))
        recorded_area = first["cumulative"][radial_index][
            "area_native_pixel2"]["valid_native_pixel2"]
        if not math.isclose(
                area_i, float(recorded_area), rel_tol=1e-10,
                abs_tol=1e-6):
            raise ValueError(
                "Stage 04.3 radial masks disagree with frozen Stage 04.2")
        # Stage 04.2's receipt authenticates a *self-consistent*
        # dictionary, not the truth of each saved radial measurement.
        # Stage 04.3 already has exact original weights at all radii;
        # independently recompute the flux and diagonal variance here.
        # The prior verification only recomputed the outermost aperture,
        # allowing internally re-receipted stale inner flux/sigma through.
        stored = first["cumulative"][radial_index]
        for field, original, measured in (
                ("signed flux", stored["signed_valid_flux_nMgy"],
                 direct_flux),
                ("formal diagonal sigma",
                 stored["formal_diagonal_sigma_nMgy"], direct_sigma)):
            if original is None and measured is None:
                continue
            if (original is None or measured is None or
                    not math.isfinite(float(original)) or
                    not math.isfinite(float(measured)) or
                    not math.isclose(
                        float(original), float(measured),
                        rel_tol=1.e-11, abs_tol=1.e-8)):
                raise ValueError(
                    "Stage 04.3 stale frozen Stage 04.2 radial " + field)
        radial_rows.append({
            "radius_arcsec": float(radial_radius),
            "valid_native_pixel2": area_i,
            "frozen_valid_flux_nMgy": first["cumulative"][radial_index][
                "signed_valid_flux_nMgy"],
            "derivative_vector_for_b0_bx_by": [
                -area_i, -moment_x, -moment_y],
            "provisional_source_exclusion_area_native_pixel2":
                source_only_area,
            "flux_added_if_only_source_mask_disabled_nMgy":
                source_only_flux,
            # An entirely source-masked baseline has no measured flux.
            # Still, after relaxing ONLY that mask, the exposed original
            # pixels can be integrated exactly; keep the baseline null.
            # This is a conditional CONTAMINANT+TARGET flux, not repaired
            # galaxy photometry or a missing-light estimate.
            "conditional_relaxed_flux_nMgy": (
                float((first["cumulative"][radial_index][
                    "signed_valid_flux_nMgy"] or 0.) + source_only_flux)
                if (first["cumulative"][radial_index][
                    "signed_valid_flux_nMgy"] is not None or
                    source_only_area > 0.) else None),
            "relaxed_flux_measured_only_after_unmasking_source": bool(
                first["cumulative"][radial_index][
                    "signed_valid_flux_nMgy"] is None and
                source_only_area > 0.),
            "status": ("CONDITIONAL_EXACT_SKY_RESPONSE"
                       if area_i > 0. else
                       "NO_VALID_NATIVE_PIXELS_NO_SKY_INFERENCE"),
        })

    # Disjoint annular derivatives must equal differences of the same
    # nested cumulative-aperture derivatives.  This introduces no
    # assumption of independence between adjacent or cumulative rings.
    radial_rings = []
    for i in range(1, len(radial_rows)):
        inner, outer_row = radial_rows[i-1], radial_rows[i]
        vector = [float(a-b) for a, b in zip(
            outer_row["derivative_vector_for_b0_bx_by"],
            inner["derivative_vector_for_b0_bx_by"])]
        radial_rings.append({
            "inner_radius_arcsec": inner["radius_arcsec"],
            "outer_radius_arcsec": outer_row["radius_arcsec"],
            "valid_native_pixel2": float(
                outer_row["valid_native_pixel2"] -
                inner["valid_native_pixel2"]),
            "derivative_vector_for_b0_bx_by": vector,
            "source_only_mask_area_native_pixel2": float(
                outer_row["provisional_source_exclusion_area_native_pixel2"] -
                inner["provisional_source_exclusion_area_native_pixel2"]),
            "source_only_flux_change_nMgy": float(
                outer_row["flux_added_if_only_source_mask_disabled_nMgy"] -
                inner["flux_added_if_only_source_mask_disabled_nMgy"]),
            "status": ("CONDITIONAL_EXACT_RING_RESPONSE"
                       if outer_row["valid_native_pixel2"] >
                       inner["valid_native_pixel2"] + 1.e-10
                       else "NO_VALID_ANNULAR_AREA_NO_SKY_INFERENCE"),
        })
    area = radial_rows[-1]["valid_native_pixel2"]
    dx_sum_scenario = -radial_rows[-1][
        "derivative_vector_for_b0_bx_by"][1]
    dy_sum_scenario = -radial_rows[-1][
        "derivative_vector_for_b0_bx_by"][2]
    dx_sum = dx_sum_scenario - area*(projected.x_center-xref)
    dy_sum = dy_sum_scenario - area*(projected.y_center-yref)

    relaxed = measure_native_aperture(
        frame, outer, provenance,
        effective_source_mask=np.zeros(frame.shape, dtype=bool),
        manual_exclude=manual, sky=sky)
    original_flux = last["signed_valid_flux_nMgy"]
    relaxed_flux = relaxed.signed_valid_flux_nMgy
    # An original galaxy aperture can be entirely source-masked.  The
    # M30 counterfactual is then still observable if relaxing ONLY the
    # source mask reveals valid native pixels. Verify that result against
    # independent direct reintegration too; do not restrict this identity
    # to cases that already had a baseline observed-area flux.
    predicted_relaxed_flux = radial_rows[-1]["conditional_relaxed_flux_nMgy"]
    if (relaxed_flux is None) != (predicted_relaxed_flux is None):
        raise ValueError(
            "Stage 04.3 source-only radial mask counterfactual "
            "observability disagrees with independent native integration")
    if relaxed_flux is not None and not math.isclose(
            predicted_relaxed_flux, relaxed_flux,
            rel_tol=1.e-11, abs_tol=1.e-7):
        raise ValueError(
            "Stage 04.3 source-only radial mask counterfactual "
            "disagrees with independent native integration")
    candidate_fluxes = [
        float(curve["cumulative"][-1]["signed_valid_flux_nMgy"])
        for curve in curves
        if (math.isclose(float(curve["radii_arcsec"][-1]), radius,
                         rel_tol=1e-10, abs_tol=1e-10) and
            curve["cumulative"][-1]["signed_valid_flux_nMgy"] is not None)
    ]
    blanks = estimate_blank_annulus_covariance(
        frame, outer, first["radii_arcsec"],
        effective_source_mask=source, manual_exclude=manual,
        sky=sky)
    small_scale = _smaller_blank_scales(
        frame, outer, first["radii_arcsec"], source, manual,
        full_scale_available=(blanks["annulus_covariance_nMgy2"] is not None),
        sky=sky)
    report = {
        "schema": _SCHEMA,
        "project_id": provenance.project_id,
        "dataset_id": frame.dataset_id,
        "draft_sha256": provenance.scene_draft_sha256,
        "catalogue_sha256": provenance.catalogue_sha256,
        "source_policy_id": provenance.source_policy_id,
        "native_input_sha256": frame.native_input_sha256,
        "stage_04_2_receipts_and_outer_apertures_verified": len(curves),
        "conditional_sky_scenario": asdict(sky),
        "baseline_geometry": hypothesis,
        "reference_semimajor_arcsec": radius,
        "reference_valid_flux_nMgy": original_flux,
        "radial_sky_response": {
            "schema": "desidera.photometry.radial-sky-jacobian/1",
            "sky_reference_native_xy": [xref, yref],
            "sky_parameter_order": [
                "b0_nMgy_per_pixel",
                "bx_nMgy_per_pixel2",
                "by_nMgy_per_pixel2"],
            "cumulative": radial_rows,
            "disjoint_annuli": radial_rings,
            "all_frozen_radial_areas_independently_verified": True,
            "all_frozen_radial_fluxes_independently_verified": True,
            "all_frozen_radial_diagonal_sigmas_independently_verified": True,
            "no_unverified_internal_radial_samples": True,
            "sky_parameter_covariance": None,
            "confidence_interval_or_sigma_calibrated": False,
            "conditional_on_unchanged_scene_masks_geometry": True,
            "no_sky_or_psf_model_adopted": True,
        },
        "sky_response": {
            "reference_native_xy": [xref, yref],
            "aperture_center_native_xy": [
                projected.x_center, projected.y_center],
            "gradient_reference_matches_frozen_sky_scenario": True,
            "dflux_db0_native_pixel2": -area,
            "dflux_dbx_native_pixel3": -dx_sum_scenario,
            "dflux_dby_native_pixel3": -dy_sum_scenario,
            "dflux_dbx_at_aperture_center_native_pixel3": -dx_sum,
            "dflux_dby_at_aperture_center_native_pixel3": -dy_sum,
            "sky_level_uncertainty_assumed": None,
            "sky_plane_uncertainty_assumed": None,
            "sky_science_adopted": False,
        },
        "radial_catalogue_mask_counterfactual": {
            "schema": "desidera.photometry.radial-source-mask-ablation/1",
            "radii_arcsec": list(first["radii_arcsec"]),
            "cumulative": [{
                "radius_arcsec": row["radius_arcsec"],
                "valid_area_without_source_relaxation_native_pixel2":
                    row["valid_native_pixel2"],
                "area_newly_available_source_only_native_pixel2":
                    row["provisional_source_exclusion_area_native_pixel2"],
                "signed_flux_newly_available_nMgy":
                    row["flux_added_if_only_source_mask_disabled_nMgy"],
                "conditional_flux_with_source_mask_removed_nMgy":
                    row["conditional_relaxed_flux_nMgy"],
                "baseline_underconstrained":
                    row["frozen_valid_flux_nMgy"] is None,
                "relaxed_flux_measured_only_after_unmasking_source":
                    row["relaxed_flux_measured_only_after_unmasking_source"],
            } for row in radial_rows],
            "disjoint_annuli": [{
                "inner_radius_arcsec": row["inner_radius_arcsec"],
                "outer_radius_arcsec": row["outer_radius_arcsec"],
                "area_newly_available_source_only_native_pixel2":
                    row["source_only_mask_area_native_pixel2"],
                "signed_flux_newly_available_nMgy":
                    row["source_only_flux_change_nMgy"],
            } for row in radial_rings],
            "outer_native_reintegration_verified": True,
            "outer_reintegration_includes_null_baseline_case": True,
            "outer_reintegration_checks_null_observability_as_well_as_flux": True,
            "relaxed_flux_may_be_observable_with_null_baseline": True,
            "no_missing_galaxy_light_claim_for_relaxed_flux": True,
            "hard_and_manual_masks_unchanged": True,
            "positive_or_negative_flux_is_not_membership_proof": True,
            "missing_target_light_not_reconstructed": True,
            "no_psf_or_source_model_subtracted": True,
            "adopted_for_science": False,
        },
        "catalogue_mask_counterfactual": {
            "flux_with_provisional_source_mask_nMgy": original_flux,
            "flux_without_provisional_source_mask_nMgy": relaxed_flux,
            "conditional_delta_nMgy": (
                float(relaxed_flux - original_flux)
                if relaxed_flux is not None and original_flux is not None
                else None),
            "manual_and_hard_exclusions_retained": True,
            "delta_is_not_missing_galaxy_light": True,
        },
        "geometry_hypothesis_span": {
            "n_valid_hypotheses": len(candidate_fluxes),
            "min_flux_nMgy": min(candidate_fluxes) if candidate_fluxes else None,
            "max_flux_nMgy": max(candidate_fluxes) if candidate_fluxes else None,
            "span_nMgy": (
                max(candidate_fluxes) - min(candidate_fluxes)
                if candidate_fluxes else None),
            "not_a_statistical_uncertainty": True,
        },
        "blank_apertures": blanks,
        "smaller_scale_blank_diagnostics": small_scale,
        "adopted_for_science": False,
        "covariance_calibrated_for_science": False,
        "status": "CONDITIONAL_SYSTEMATICS_DIAGNOSTIC_ONLY",
        "warnings": [
            "SKY_RESPONSE_IS_A_DERIVATIVE_NOT_A_SKY_ERROR_BAR",
            "RELAXED_CATALOGUE_MASK_INCLUDES_CONTAMINANTS_AND_TARGET_MEMBERS",
            "GEOMETRY_SPAN_COMPARES_DIFFERENT_SKY_FOOTPRINTS",
            "NO_TOTAL_FLUX_OR_PSF_CORRECTION_INFERRED",
        ] + ([
            "BLANK_COVARIANCE_CONDITIONAL_ON_UNADOPTED_SKY_PLANE"
        ] if sky.method != "zero_additional_correction" else []),
    }
    json.dumps(report, allow_nan=False)
    return report
