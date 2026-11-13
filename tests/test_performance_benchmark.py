"""Non-gating smoke test of deterministic synthetic performance measurements.

Performance ratios are reported, not thresholded: CI runners differ.
Scientific numerical identities and immutable input constraints are
tested separately in the photometry regression modules.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

pytest.importorskip("astropy")
pytest.importorskip("photutils")


def test_performance_benchmark_runs_without_original_fits_or_thresholds():
    script = Path(__file__).resolve().parents[1] / (
        "scripts/benchmark_desidera_performance.py")
    spec = importlib.util.spec_from_file_location(
        "_desidera_performance_benchmark", script)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    report = module.run(size=192, repeats=1)
    assert report["schema"] == (
        "desidera.performance.synthetic-operator-benchmark/1")
    assert report["n_native_pixels"] == 192**2
    assert set(report["comparison"]) == {
        "native_frame_two_materializations_vs_one",
        "ten_radii_full_mask_copies_vs_immutable_snapshot_reuse",
        "exact_blank_ring_rasterization_vs_integer_template_translation",
        "sky_frozen_injection_dense_xy_vs_broadcast_reusable_distance",
        "native_wcs_grid_header_serialization_vs_per_frame_fingerprint",
        "sky_six_fixed_trials_full_edt_vs_exact_guard_reuse",
    }
    for row in report["comparison"].values():
        assert row["reference_median_ms"] >= 0.
        assert row["optimized_median_ms"] >= 0.
        assert row["wall_time_is_hardware_and_cache_dependent"]
    ring = report["comparison"][
        "exact_blank_ring_rasterization_vs_integer_template_translation"]
    assert ring["reference_exact_rasterizations"] == 160
    assert ring["optimized_exact_rasterizations"] == 8
    assert ring["independent_exact_overlap_area_verified"]
    injection = report["comparison"][
        "sky_frozen_injection_dense_xy_vs_broadcast_reusable_distance"]
    assert injection["full_synthetic_flux_reference_verified"]
    assert injection["nominal_persistent_coordinate_grid_bytes_avoided"] == (
        2*192*192*8)
    assert injection["full_float64_coordinate_grids_avoided"] == 2
    assert injection["repeated_elliptical_distance_maps_avoided"] == 2
    assert injection["peak_rss_not_measured"]
    wcs_qa = report["comparison"][
        "native_wcs_grid_header_serialization_vs_per_frame_fingerprint"]
    assert wcs_qa["wcs_header_serializations_reference"] == 10
    assert wcs_qa["wcs_header_serializations_optimized"] == 1
    assert wcs_qa["exact_header_and_shape_identity_verified"]
    guards = report["comparison"][
        "sky_six_fixed_trials_full_edt_vs_exact_guard_reuse"]
    assert guards["sky_grid_trials"] == 7
    assert guards["reference_distance_transforms_maximum"] == 14
    assert guards["optimized_distance_transforms_maximum"] == 4
    assert guards["exact_native_source_target_geometry_verified"]
    assert guards["cache_lifetime"] == "ONE_STAGE03_CALL_NO_GLOBAL_STATE"
    json.dumps(report, allow_nan=False)
    with pytest.raises(ValueError, match="Synthetic benchmark size"):
        module.run(size=32, repeats=1)
