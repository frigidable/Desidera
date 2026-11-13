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

    # The native WCS/shape identity is unchanged during a whole
    # growth curve. Previously every cumulative aperture serialized
    # that same WCS header again; preserve the exact original SHA
    # while doing that expensive conversion only once per frame.
    from lsbpipe.scene import grid_fingerprint
    frozen_frame = make_frame()
    expected_grid = grid_fingerprint(frozen_frame.wcs, frozen_frame.shape)
    if expected_grid != frozen_frame.native_grid_fingerprint:
        raise RuntimeError("Memoized native grid fingerprint changed WCS identity")

    def repeated_grid_header_serialization():
        return [grid_fingerprint(frozen_frame.wcs, frozen_frame.shape)
                for _ in range(n_radii)]

    def cached_native_grid_header_identity():
        return [frozen_frame.native_grid_fingerprint
                for _ in range(n_radii)]

    grid_identity_result = _compare(
        repeated_grid_header_serialization,
        cached_native_grid_header_identity, repeats)
    grid_identity_result.update({
        "wcs_header_serializations_reference": n_radii,
        "wcs_header_serializations_optimized": 1,
        "exact_header_and_shape_identity_verified": True,
        "scientific_wcs_affine_projection_checks_still_required": True,
    })

    # The six predeclared same-project Stage03 sky-grid variants often
    # change only tile boundaries, not original target/source veto
    # geometry. Reusing baseline Boolean EDT guards is safe only
    # after an exact raster comparison (the implementation does this).
    from lsbpipe.background import _automatic_scene_sky_pixels

    sky_target = np.zeros(shape, bool)
    yy1, xx1 = np.ogrid[:size, :size]
    sky_target[(xx1-size/2.)**2+(yy1-size/2.)**2 < (size/12.)**2] = True
    sky_source = np.zeros(shape, bool)
    sky_source[size//4:size//4+9, size//5:size//5+11] = True
    sky_source[size//2:size//2+8, 4*size//5:4*size//5+10] = True
    sky_dataset = {
        "pixel_shape": shape, "wcs": wcs,
        "native_hard": np.zeros(shape, bool),
        "data": {"image": image, "invvar": ivar},
    }
    sky_preview = {
        "manual_exclude": np.zeros(shape, bool),
        "target_protect": np.zeros(shape, bool),
        "isophote_candidate": np.zeros(shape, bool),
    }
    sky_contamination = {
        "auto_mask": sky_source,
        "protected": sky_target,
        "source_inventory": {
            "entries": [], "stored_count": 0,
            "total_relevant_sources": 0, "truncated": 0,
            "screening_evidence_incomplete": False,
        },
    }
    variants = (
        {}, {"tile_edge_native_pix": 72},
        {"tile_edge_native_pix": 96,
         "grid_offset_xy_native_pix": (48, 48)},
        {"tile_edge_native_pix": 128},
        {"target_guard_factor": 1.5},
        {"source_guard_factor": 1.5},
        {"sample_policy": "SPATIAL_8X8_PROPORTIONAL_V2"},
    )

    def sky_grid_trials(cache):
        results = []
        for options in variants:
            train, hold, qa = _automatic_scene_sky_pixels(
                sky_dataset, sky_preview, sky_contamination,
                _guard_cache=cache, **options)
            results.append((int(np.count_nonzero(train)),
                            int(np.count_nonzero(hold)),
                            qa["target_outskirts_guard_native_pix"],
                            qa["source_neighbour_guard_native_pix"],
                            tuple(qa["training_quadrants"]),
                            tuple(qa["holdout_quadrants"])))
        return results

    if sky_grid_trials(None) != sky_grid_trials({}):
        raise RuntimeError(
            "Sky cached guards changed deterministic selected-tile QA")
    sky_guard_result = _compare(
        lambda: sky_grid_trials(None),
        lambda: sky_grid_trials({}), repeats)
    sky_guard_result.update({
        "sky_grid_trials": len(variants),
        "reference_distance_transforms_maximum": 2*len(variants),
        "optimized_distance_transforms_maximum": 4,
        "exact_native_source_target_geometry_verified": True,
        "cache_lifetime": "ONE_STAGE03_CALL_NO_GLOBAL_STATE",
        "cache_memory_cost_not_peak_rss_measured": True,
    })

    # Stage04 pixel-centre diagnostic sectors formerly constructed
    # two full float64 native coordinate grids and an additional
    # 2×N dense matrix for tangent-frame multiplication. Preserve
    # the exact eight-sector labels with 1-D native broadcasting.
    from lsbpipe.photometry.aperture_flux import (
        _sector_ids_from_broadcast_native_centres)
    inv_tangent = np.array([[.91, .04], [-.075, 1.07]])
    x_offset, y_offset = -size*.47-.163, -size*.52+.217

    def dense_native_sectors():
        ygrid, xgrid = np.indices(shape, dtype=float)
        xgrid += x_offset
        ygrid += y_offset
        tangent = inv_tangent @ np.stack((xgrid.ravel(), ygrid.ravel()))
        angles = np.mod(
            np.arctan2(tangent[0], tangent[1]), 2*np.pi)
        return (np.floor(angles/(np.pi/4)).astype(int) % 8).reshape(shape)

    def broadcast_native_sectors():
        return _sector_ids_from_broadcast_native_centres(
            shape, x_offset, y_offset, inv_tangent)

    if not np.array_equal(dense_native_sectors(), broadcast_native_sectors()):
        raise RuntimeError("Broadcast Stage04 sectors changed native PA labels")
    sector_result = _compare(
        dense_native_sectors, broadcast_native_sectors, repeats)
    sector_result.update({
        "full_coordinate_grids_not_constructed": 2,
        "dense_tangent_input_stack_not_constructed": True,
        "eight_native_pa_sectors_identical": True,
        "pixel_center_qa_not_exact_angular_sector_photometry": True,
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
            "native_wcs_grid_header_serialization_vs_per_frame_fingerprint":
                grid_identity_result,
            "sky_six_fixed_trials_full_edt_vs_exact_guard_reuse":
                sky_guard_result,
            "native_pa_sector_dense_xy_vs_broadcast_tangent":
                sector_result,
        },
        "interpretation": [
            "Seven distinct computational kernels, not end-to-end FITS loading time",
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
