"""One-click project QA must distinguish missing, running, current and stale."""
import json

import numpy as np
from types import SimpleNamespace

import pytest

pytest.importorskip("pytestqt")
pytest.importorskip("PySide6")
from PySide6.QtWidgets import QApplication

from lsbpipe.gui.photometry_workspace import PhotometryWorkspace
from lsbpipe.gui.project_qa_bundle import (
    build_project_diagnostic_bundle, serialize_project_diagnostic_bundle,
    _compact_ablation, _compact_apertures, _guard_raw_arrays)
from lsbpipe.photometry.sky_masks import SkyMaskDraft


def _bound_workspace(qtbot, tmp_path):
    widget = PhotometryWorkspace()
    qtbot.addWidget(widget)
    widget._project_root = tmp_path
    widget._bound_project_id = "one-project-uuid"
    widget._bound_source_identity = ("local-root", ("dr11", "south", "1554p180"))
    scene = widget.scene_editor
    scene._project_root = tmp_path
    scene._project_id = widget._bound_project_id
    scene._project_target = (155.438012, 18.088025)
    scene._draft = SkyMaskDraft(widget._bound_project_id)
    return widget


def test_copy_qa_requires_project_without_mutating_clipboard(qtbot):
    widget = PhotometryWorkspace()
    qtbot.addWidget(widget)
    before = "sentinel: leave clipboard untouched"
    QApplication.clipboard().setText(before)
    widget.copy_project_qa_btn.click()
    assert QApplication.clipboard().text() == before
    assert "Open a saved galaxy project" in widget.project_qa_status.text()
    assert widget.stage_tabs.count() == 4


def test_single_copy_includes_all_stages_and_reproducibility(qtbot, tmp_path):
    widget = _bound_workspace(qtbot, tmp_path)
    before = sorted(tmp_path.iterdir())
    assert not widget.sky_stage.run_btn.isEnabled()
    assert not widget.aperture_stage.run_btn.isEnabled()
    widget.copy_project_qa_btn.click()
    bundle = json.loads(QApplication.clipboard().text())
    assert bundle["schema"] == "desidera.gui.project-diagnostic-bundle/1"
    assert bundle["project"]["project_id"] == "one-project-uuid"
    assert bundle["environment"]["python"]
    assert bundle["science_hold"] is True
    assert bundle["scientific_adoption"] is False
    assert bundle["workers_started_by_copy"] is False
    assert list(bundle["stages"]) == [
        "aperture", "background", "psf0", "scene_masks"]
    assert bundle["stages"]["scene_masks"]["state"] == "NOT_LOADED"
    assert bundle["stages"]["background"]["state"] == "BLOCKED_BY_SCENE_INPUTS"
    assert bundle["stages"]["aperture"]["state"] == "BLOCKED_BY_SCENE_INPUTS"
    assert all(v["state"] == "NOT_RUN"
               for v in bundle["stages"]["psf0"]["bands"].values())
    assert sorted(tmp_path.iterdir()) == before
    assert widget.sky_stage._future is None
    assert widget.aperture_stage._future is None
    assert "Copied all-stage project QA" in widget.project_qa_status.text()


def test_stale_background_is_explicitly_withheld(qtbot, tmp_path):
    widget = _bound_workspace(qtbot, tmp_path)
    widget.sky_stage._result = {"old_flux_that_must_not_leak": 999}
    widget.sky_stage._result_identity = ("wrong",)
    report = build_project_diagnostic_bundle(widget)
    sky = report["stages"]["background"]
    assert sky["state"] == "STALE_RESULT_WITHHELD"
    assert sky["current_result_provenance_matches"] is False
    assert "result" not in sky
    assert "999" not in json.dumps(sky)
    assert widget.sky_stage._result is not None  # copy never clears state


def test_per_stage_collection_failure_does_not_hide_other_stages(
        qtbot, tmp_path, monkeypatch):
    widget = _bound_workspace(qtbot, tmp_path)
    def fail():
        raise RuntimeError("controlled source registry failure")
    monkeypatch.setattr(widget.scene_editor, "_complete_scene_qa_snapshot", fail)
    bundle = json.loads(serialize_project_diagnostic_bundle(widget))
    scene = bundle["stages"]["scene_masks"]
    assert scene["state"] == "COLLECTION_FAILED"
    assert scene["error_type"] == "RuntimeError"
    assert "controlled source registry failure" in scene["error_message"]
    assert bundle["stages"]["background"]["state"] == "BLOCKED_BY_SCENE_INPUTS"
    assert bundle["scientific_adoption"] is False


def test_ablation_keeps_every_scenario_and_radial_bin_without_huge_octants():
    bins = [
        dict(a_mid_arcsec=1., outer_arcsec=2., pixels=30, valid_pixels=30,
             coverage=1., sectors_with_data=8, sectors_positive=8,
             median_smoothed_rms=2., supported=True,
             sector_contrast_by_octant=list(range(8))),
        dict(a_mid_arcsec=2., outer_arcsec=3., pixels=40, valid_pixels=40,
             coverage=1., sectors_with_data=8, sectors_positive=4,
             median_smoothed_rms=.1, supported=False,
             sector_contrast_by_octant=list(range(8))),
        dict(a_mid_arcsec=3., outer_arcsec=4., pixels=50, valid_pixels=25,
             coverage=.5, sectors_with_data=7, sectors_positive=2,
             median_smoothed_rms=-.1, supported=False,
             sector_contrast_by_octant=list(range(8))),
    ]
    source = {"scenarios": [
        {"policy": "none", "component_pixels": 50, "sga_radial_bins": bins},
        {"policy": "both", "component_pixels": 22, "sga_radial_bins": bins},
    ]}
    out = _compact_ablation(source)
    assert len(out["scenarios"]) == 2
    assert "sga_radial_bins" in source["scenarios"][0]
    assert "sga_radial_bins" not in out["scenarios"][0]
    trace = out["scenarios"][0]["sga_radial_trace"]
    assert len(trace["rows"]) == 3
    assert trace["all_radial_rows_included"]
    assert trace["rows"][2][4] == .5
    assert trace["rows"][2][7] == -.1
    assert trace["rows"][1][-1] is False
    assert len(trace["octant_contrasts_at_boundaries"]) == 3


def test_aperture_bundle_keeps_signed_fluxes_mask_areas_and_04_3():
    data = {
        "status": "DIAGNOSTIC", "project_id": "p", "band": "i",
        "systematics_04_3": {"sky_response": {"dflux_db0_native_pixel2": -99.}},
        "curves": [{
            "geometry_hypothesis": {"name": "anchor_circle"},
            "max_abs_annulus_closure_nMgy": 2e-14,
            "cumulative": [{
                "project_id": "p", "recipe_sha256": "test",
                "geometry_icrs": {"semimajor_arcsec": 10.},
                "area_native_pixel2": {"valid_native_pixel2": 99},
                "signed_valid_flux_nMgy": -3.,
                "formal_diagonal_sigma_nMgy": 2.,
                "sector_qa": [{"sector_id_pa_45deg": 1, "valid_fraction": .8}],
            }],
            "independently_integrated_annuli": [
                {"signed_flux_nMgy": -3., "ab_surface_brightness_mag_arcsec2": None}],
        }],
    }
    compact = _compact_apertures(data)
    c = compact["curves"][0]
    assert c["cumulative"][0]["signed_valid_flux_nMgy"] == -3.
    assert "project_id" not in c["cumulative"][0]
    assert c["cumulative_sector_qa"][0][0]["valid_fraction"] == .8
    assert c["independently_integrated_annuli"][0]["signed_flux_nMgy"] == -3.
    assert compact["systematics_04_3"]["sky_response"][
        "dflux_db0_native_pixel2"] == -99.
    json.dumps(compact, allow_nan=False)


def test_loaded_synthetic_scene_collects_current_sky_and_aperture_then_vetoes_stale(
        qtbot, tmp_path):
    from test_aperture_flux import frame

    widget = _bound_workspace(qtbot, tmp_path)
    scene = widget.scene_editor
    scene.set_band("r")
    native = frame(value=.003)
    source = "DR11 south · 1554p180 · r"
    scene._source_key = (
        "fake-data-root", widget._bound_project_id,
        ("dr11", "south", "1554p180", "r"), "native")
    scene._loaded = {
        "project_id": widget._bound_project_id,
        "status": "NATIVE_SCENE_CUTOUT_NOT_SCIENCE_ADOPTED",
        "source": source, "band": "r", "pixel_shape": native.shape,
        "display_width_arcmin": 4., "project_field_width_arcmin": 4.,
        "hdu_indices": {"image": 0},
        "native_hard": np.zeros(native.shape, dtype=bool),
        "wcs": native.wcs,
    }
    scene._contamination = {
        "status": "CATALOGUE_CONTAMINATION_MASK_NOT_ADOPTED",
        "project_id": widget._bound_project_id,
        "band": "r", "dataset_source": source,
        "draft_digest": scene._draft.digest,
        "catalog_sha256": "a" * 64,
        "auto_mask": np.zeros(native.shape, dtype=bool),
        "mask": np.zeros(native.shape, dtype=bool),
    }
    sky_identity = widget.sky_stage._snapshot()[0]
    widget.sky_stage._result_identity = sky_identity
    widget.sky_stage._result = {
        "project_id": widget._bound_project_id, "band": "r",
        "dataset_source": source, "draft_sha256": scene._draft.digest,
        "catalog_sha256": "a" * 64,
        "status": "SKY_DIAGNOSTIC_ONLY",
        "methods": {"zero": {"status": "BASELINE"}},
        "sky_grid_sensitivity_qa": {"trials": [
            {"name": "test", "status": "UNADOPTED"}]},
    }
    aperture_identity = widget.aperture_stage._snapshot()[0]
    widget.aperture_stage._result_identity = aperture_identity
    widget.aperture_stage._result = {
        "project_id": widget._bound_project_id, "band": "r",
        "dataset_source": source, "scene_draft_sha256": scene._draft.digest,
        "catalog_sha256": "a" * 64,
        "status": "DIRECT_FINITE_APERTURE_DIAGNOSTICS_NOT_ADOPTED",
        "systematics_04_3": {
            "blank_apertures": {"annulus_covariance_nMgy2": None}},
        "curves": [{
            "geometry_hypothesis": {"name": "anchor_circle"},
            "cumulative": [{
                "geometry_icrs": {"semimajor_arcsec": 105.},
                "signed_valid_flux_nMgy": 2.5,
                "formal_diagonal_sigma_nMgy": .5,
                "area_native_pixel2": {"valid_native_pixel2": 50.},
                "sector_qa": [],
            }],
            "independently_integrated_annuli": [{
                "signed_flux_nMgy": -1.2, "valid_fraction_of_geometric": .7}],
        }],
    }
    report = json.loads(serialize_project_diagnostic_bundle(widget))
    assert report["stages"]["scene_masks"]["state"] == "CURRENT_LOADED"
    assert report["stages"]["scene_masks"]["native_fits"]["band"] == "r"
    assert report["stages"]["background"]["state"] == "CURRENT_DIAGNOSTIC"
    assert report["stages"]["aperture"]["state"] == "CURRENT_DIAGNOSTIC"
    assert report["stages"]["aperture"]["result"]["curves"][0][
        "independently_integrated_annuli"][0]["signed_flux_nMgy"] == -1.2
    assert "BLANK_ANNULUS_COVARIANCE_UNAVAILABLE" in report[
        "diagnostic_findings"]
    assert "SIGNED_NEGATIVE_ANNULAR_FLUX_REQUIRES_SKY_QA" in report[
        "diagnostic_findings"]
    assert "APERTURE_VALID_AREA_BELOW_90_PERCENT" not in report[
        "diagnostic_findings"]

    # An asynchronous band change must never let the old Stage 03/04
    # results pass as evidence for the newly selected native scene.
    scene._band = "i"
    stale = json.loads(serialize_project_diagnostic_bundle(widget))
    assert stale["stages"]["scene_masks"]["state"] == "STALE_SCENE_INPUTS"
    for name in ("background", "aperture"):
        assert stale["stages"][name]["state"] == "STALE_RESULT_WITHHELD"
        assert "result" not in stale["stages"][name]
    assert "SCENE_BAND_MISMATCH" in stale["diagnostic_findings"]
    assert widget.sky_stage._result is not None
    assert widget.aperture_stage._result is not None


def test_nonfinite_debug_values_and_pixel_grids_never_break_qa_copy():
    incoming = {
        "sky_nonfinite": float("nan"),
        "annulus_nonfinite": np.float64("-inf"),
        "raw_fits_pixel_grid": np.ones((400, 400)),
        "small_array_with_nan": np.array([1., np.nan]),
        "signed_negative": -8.5,
    }
    result = _guard_raw_arrays(incoming)
    assert result["sky_nonfinite"]["diagnostic_nonfinite_float"] == "nan"
    assert result["annulus_nonfinite"]["diagnostic_nonfinite_float"] == "-inf"
    assert result["raw_fits_pixel_grid"]["omitted_raw_array"]
    assert result["raw_fits_pixel_grid"]["shape"] == [400, 400]
    assert result["small_array_with_nan"][1][
        "diagnostic_nonfinite_float"] == "nan"
    assert result["signed_negative"] == -8.5
    json.dumps(result, allow_nan=False)
