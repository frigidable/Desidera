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
    json.dumps(report, allow_nan=False)
    with pytest.raises(ValueError, match="Synthetic benchmark size"):
        module.run(size=32, repeats=1)
