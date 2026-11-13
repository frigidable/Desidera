"""Reproducible, offline computational-kernel benchmark for Desidera.

Measures ONLY selected optimized inner operations on synthetic native
arrays. It is not a full pipeline benchmark, astrophysical certification,
or a stable CI timing threshold. Never accesses or writes original FITS.

Usage:
    python scripts/benchmark_desidera_performance.py --size 512 --repeats 3
"""
from __future__ import annotations

import argparse
import json
import math
import statistics
import time

import numpy as np
from astropy.wcs import WCS

from lsbpipe.photometry.aperture_flux import (
    ApertureGeometryICRS, _mask_for, _snapshot_native_boolean_mask,
    project_icrs_ellipse_to_native,
)
from lsbpipe.photometry.aperture_systematics import (
    _ring_aperture, _translated_exact_aperture_overlap_slices,
)
from lsbpipe.photometry.contracts import NativeBandFrame


def _timed(func, repeats):
    durations = []
    for _ in range(repeats):
        start = time.perf_counter()
        func()
        durations.append((time.perf_counter() - start) * 1000.)
    return float(statistics.median(durations))


def _compare(old, new, repeats):
    slow = _timed(old, repeats)
    fast = _timed(new, repeats)
    return {
        "reference_median_ms": slow,
        "optimized_median_ms": fast,
        "observed_reference_over_optimized_ratio": (
            slow / fast if fast > 0 else None),
        "repeats": repeats,
        "wall_time_is_hardware_and_cache_dependent": True,
    }


def _wcs(size):
    w = WCS(naxis=2)
    w.wcs.crpix = [size / 2. + .5, size / 2. + .5]
    w.wcs.crval = [155.4, 18.]
    w.wcs.cd = np.diag([-.262, .262]) / 3600.
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    return w


def run(*, size, repeats):
    if size < 192 or size > 2048 or repeats < 1 or repeats > 10:
        raise ValueError("Synthetic benchmark size must be 192–2048 and repeats 1–10")
    shape = (size, size)
    rng = np.random.default_rng(20261113)
    image = rng.normal(.02, .1, shape).astype(np.float32)
    ivar = np.full(shape, 4., dtype=np.float32)
    nexp = np.full(shape, 2, dtype=np.int16)
    maskbits = np.zeros(shape, dtype=np.int16)
    wcs = _wcs(size)

    def make_frame():
        return NativeBandFrame(
            image_nMgy=image, invvar_nMgy2=ivar,
            nexp=nexp, maskbits=maskbits, native_wcs=wcs,
            release="dr11", region="south", brick="1554p180", band="r",
            source_evidence={"synthetic_performance_fixture": True})

    def reference_two_frame_materializations():
        first = make_frame()
        second = NativeBandFrame(
            image_nMgy=first.image_nMgy,
            invvar_nMgy2=first.invvar_nMgy2,
            nexp=first.nexp, maskbits=first.maskbits,
            native_wcs=first.wcs, release="dr11", region="south",
            brick="1554p180", band="r", policy=first.policy,
            source_evidence={"synthetic_performance_fixture": True})
        return second

    frame_result = _compare(
        reference_two_frame_materializations, make_frame, repeats)

    source = np.zeros(shape, bool)
    manual = np.zeros(shape, bool)
    source[size//3:size//3+17, size//2:size//2+19] = True
    manual[size//2:size//2+13, size//4:size//4+11] = True
    n_radii = 10

    def repeated_mask_copies():
        for _ in range(n_radii):
            _mask_for("source", source, shape)
            _mask_for("manual", manual, shape)

    def reused_frozen_masks():
        fixed_source = _snapshot_native_boolean_mask("source", source, shape)
        fixed_manual = _snapshot_native_boolean_mask("manual", manual, shape)
        for _ in range(n_radii):
            assert _mask_for("source", fixed_source, shape) is fixed_source
            assert _mask_for("manual", fixed_manual, shape) is fixed_manual

    mask_result = _compare(repeated_mask_copies, reused_frozen_masks, repeats)

    galaxy = ApertureGeometryICRS(
        ra_deg=155.4, dec_deg=18., semimajor_arcsec=5.,
        axis_ratio=.68, pa_deg_east_of_north=31.)
    projected = project_icrs_ellipse_to_native(galaxy, wcs)
    radii = (0.65, 1.1, 1.6, 2.2, 2.8, 3.5, 4.2, 5.)
    floor_x = math.floor(projected.x_center)
    floor_y = math.floor(projected.y_center)
    phase_x = projected.x_center-floor_x
    phase_y = projected.y_center-floor_y
    lattice = [
        (int((i % 5 + 1) * size/6),
         int((i // 5 + 1) * size/5))
        for i in range(20)]

    def recompute_exact_ring_weights():
        total = 0.
        for x, y in lattice:
            for i in range(len(radii)):
                mask = _ring_aperture(
                    x+phase_x, y+phase_y, projected, radii, i).to_mask(
                        method="exact")
                large, small = mask.get_overlap_slices(shape)
                if large is not None:
                    total += float(np.sum(mask.data[small]))
        return total

    def translate_exact_ring_weights():
        templates = [
            _ring_aperture(
                projected.x_center, projected.y_center,
                projected, radii, i).to_mask(method="exact")
            for i in range(len(radii))]
        total = 0.
        for x, y in lattice:
            dx, dy = x-floor_x, y-floor_y
            for mask in templates:
                large, small = _translated_exact_aperture_overlap_slices(
                    mask, shape, dx, dy)
                if large is not None:
                    total += float(np.sum(mask.data[small]))
        return total

    expected = recompute_exact_ring_weights()
    reused = translate_exact_ring_weights()
    if not np.isclose(expected, reused, rtol=1.e-10, atol=1.e-7):
        raise RuntimeError("Integer-translated exact ring operator changed numerical area")
    ring_result = _compare(
        recompute_exact_ring_weights, translate_exact_ring_weights, repeats)
    ring_result.update({
        "ring_count": len(radii), "blank_positions": len(lattice),
        "reference_exact_rasterizations": len(radii)*len(lattice),
        "optimized_exact_rasterizations": len(radii),
        "independent_exact_overlap_area_verified": True,
    })

    # Reference Stage03 implementation materialized two full float64
    # coordinate grids and recomputed the same elliptical radius once
    # for each factor. The optimized one uses broadcast 1-D offsets
    # and a single reusable ellipse-distance map.
    x_anchor = size/2. - .37
    y_anchor = size/2. + .21
    trial_scales = (14., 28., 56.)

    def dense_stage03_distance_trials():
        yy, xx = np.indices(shape, dtype=float)
        selection = np.hypot(xx-x_anchor, yy-y_anchor) < 70.
        output = []
        for h in trial_scales:
            ellipse = np.hypot(xx-x_anchor, (yy-y_anchor)/.65)
            signal = .005*np.exp(-ellipse/h)
            output.append((float(np.sum(signal)),
                           float(np.sum(signal[selection]))))
        return output

    def broadcast_stage03_distance_trials():
        dx = np.arange(size, dtype=float)-x_anchor
        dy = np.arange(size, dtype=float)-y_anchor
        selection = np.hypot(dx[None, :], dy[:, None]) < 70.
        ellipse = np.hypot(dx[None, :], dy[:, None]/.65)
        output = []
        for h in trial_scales:
            signal = .005*np.exp(-ellipse/h)
            output.append((float(np.sum(signal)),
                           float(np.sum(signal[selection]))))
        return output

    direct = dense_stage03_distance_trials()
    broadcast = broadcast_stage03_distance_trials()
    if not np.allclose(direct, broadcast, rtol=1.e-13, atol=1.e-11):
        raise RuntimeError(
            "Broadcast Stage03 ellipse-distance operator changed injected flux")
    injection_result = _compare(
        dense_stage03_distance_trials,
        broadcast_stage03_distance_trials, repeats)
    injection_result.update({
        "n_scale_length_trials": len(trial_scales),
        "full_float64_coordinate_grids_avoided": 2,
        "nominal_persistent_coordinate_grid_bytes_avoided":
            2*size*size*np.dtype(np.float64).itemsize,
        "repeated_elliptical_distance_maps_avoided":
            len(trial_scales)-1,
        "full_synthetic_flux_reference_verified": True,
        "peak_rss_not_measured": True,
    })

    return {
        "schema": "desidera.performance.synthetic-operator-benchmark/1",
        "size_native_pixels": [size, size],
        "n_native_pixels": size*size,
        "python_wall_clock_method": "time.perf_counter",
        "comparison": {
            "native_frame_two_materializations_vs_one": frame_result,
            "ten_radii_full_mask_copies_vs_immutable_snapshot_reuse":
                mask_result,
            "exact_blank_ring_rasterization_vs_integer_template_translation":
                ring_result,
            "sky_frozen_injection_dense_xy_vs_broadcast_reusable_distance":
                injection_result,
        },
        "interpretation": [
            "Four distinct computational kernels, not end-to-end FITS loading time",
            "No network, original FITS, catalogue fetch or GUI involved",
            "No scientific sky, PSF, photometric, SBF or covariance adoption",
            "Ratio depends on machine, NumPy, Photutils, OS cache and FITS size",
            "Do not compare absolute durations between different machines",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--size", type=int, default=512)
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    report = run(size=args.size, repeats=args.repeats)
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
