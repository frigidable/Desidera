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
    _compact_ablation, _compact_apertures, _guard_raw_arrays,
    _foreground_aperture_overlap_qa, _diagnostic_findings,
    _conditional_sky_aperture_response_qa, _cross_stage_source_lineage_qa)
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
        "aperture", "background", "conditional_sky_aperture_response",
        "cross_stage_source_lineage", "foreground_aperture_overlap",
        "psf0", "scene_masks"]
    assert bundle["stages"]["cross_stage_source_lineage"]["state"] == (
        "SCENE_SOURCE_NOT_CURRENT")
    assert bundle["stages"]["conditional_sky_aperture_response"]["state"] == (
        "BACKGROUND_OR_APERTURE_NOT_CURRENT")
    assert bundle["stages"]["foreground_aperture_overlap"]["state"] == (
        "BLOCKED_BY_MISSING_OR_STALE_STAGE")
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


def test_gaia_protected_foreground_spike_is_explicit_not_subtracted():
    """An r-band Leo-P-like Gaia PSF star sits in the 18--25.2 arcsec ring."""
    inventory = {
        "status": "SOURCE_EVIDENCE_REVIEW_ONLY",
        "dataset_id": "dr11/south/1554p180/r",
        "entries": [
            {"objid": 2035,
             "gaia_dr3_source_id": 623919492384922624,
             "tractor_type": "PSF", "offset_arcsec": 3.0,
             "ra_deg": 155.4, "dec_deg": 18. + 19.870608956040662 / 3600.,
             "parallax_significance": 30.6958, "parallax_mas": 1.08808,
             "gaia_g_mag": 15.0252, "band_model_flux_nMgy": 1078.746,
             "centre_in_target_protection": True,
             "screening_disposition": "PROTECTED_OVERLAP",
             "diagnostic_inner_gaia_mask_preview_active": False,
             "evidence_flags": ["BRIGHT_STAR_WINGS_REVIEW"]},
            {"objid": 42, "gaia_dr3_source_id": None, "tractor_type": "PSF",
             "offset_arcsec": 20., "parallax_significance": 0.,
             "band_model_flux_nMgy": 3000.},
            {"objid": 43, "gaia_dr3_source_id": 12345, "tractor_type": "PSF",
             "offset_arcsec": 25., "parallax_mas": .15,
             "parallax_significance": 2.,
             "band_model_flux_nMgy": 4000.},
            {"objid": 44, "gaia_dr3_source_id": 23456, "tractor_type": "PSF",
             "offset_arcsec": 21., "parallax_mas": -1.5,
             "parallax_significance": 15.,
             "band_model_flux_nMgy": 9000.},
        ],
    }
    circle = {
        "geometry_hypothesis": {"name": "anchor_circle"},
        "radii_arcsec": [7.2, 12., 18., 25.2, 33.6],
        "cumulative": [
            {"signed_valid_flux_nMgy": 31.},
            {"signed_valid_flux_nMgy": 58.},
            {"signed_valid_flux_nMgy": 100.},
            {"signed_valid_flux_nMgy": 1226.},
            {"signed_valid_flux_nMgy": 1382.},
        ],
        "independently_integrated_annuli": [
            {"signed_flux_nMgy": 27.},
            {"signed_flux_nMgy": 42.},
            {"signed_flux_nMgy": 1126.},
            {"signed_flux_nMgy": 156.},
        ],
    }
    stages = {
        "scene_masks": {"state": "CURRENT_LOADED",
                        "source_inventory": inventory,
                        "contaminant_screening": {}},
        "background": {"state": "CURRENT_DIAGNOSTIC"},
        "aperture": {"state": "CURRENT_DIAGNOSTIC",
                     "result": {"dataset_id": "dr11/south/1554p180/r",
                                "target_anchor_icrs": [155.4, 18.],
                                "curves": [circle], "systematics_04_3": {}}},
    }
    diag = _foreground_aperture_overlap_qa(stages)
    assert diag["state"] == "FOREGROUND_ASTROMETRY_APERTURE_OVERLAP_QA_ONLY"
    assert diag["n_qualifying_sources_within_outer_aperture"] == 1
    star = diag["sources_ranked_by_tractor_model_flux"][0]
    assert star["objid"] == 2035
    assert star["distance_to_project_anchor_arcsec"] == (
        pytest.approx(19.870608956040662, abs=1e-7))
    # Deliberately contradictory source-inventory offset is from another
    # reference centre and MUST NOT determine the Stage04 ring.
    assert star["source_inventory_offset_arcsec_potentially_different_anchor"] == 3.0
    assert star["ring_inner_arcsec"] == 18.
    assert star["ring_outer_arcsec"] == 25.2
    assert star["ring_measured_signed_flux_nMgy"] == 1126.
    assert star["model_to_signed_ring_flux_ratio_not_a_fraction"] == (
        pytest.approx(1078.746 / 1126.))
    assert star["source_centre_protected"]
    assert diag["science_adopted"] is False
    assert diag["source_flux_removed"] is False
    stages["foreground_aperture_overlap"] = diag
    findings = _diagnostic_findings(stages, {})
    assert "GAIA_ASTROMETRIC_FOREGROUND_CANDIDATES_INSIDE_APERTURE" in findings
    assert "PROTECTED_BRIGHT_GAIA_SOURCE_MAY_DOMINATE_APERTURE_RING" in findings
    # No reference to the user supplied dataset is needed at runtime.
    json.dumps(diag, allow_nan=False)


def test_foreground_aperture_diagnostic_fails_closed_on_stale_band():
    base = {"scene_masks": {
                "state": "CURRENT_LOADED",
                "source_inventory": {
                    "status": "SOURCE_EVIDENCE_REVIEW_ONLY",
                    "dataset_id": "dr11/south/1554p180/r", "entries": []}},
            "aperture": {"state": "CURRENT_DIAGNOSTIC",
                         "result": {"dataset_id": "dr11/south/1554p180/i",
                                    "curves": []}}}
    r = _foreground_aperture_overlap_qa(base)
    assert r["state"] == "SOURCE_INVENTORY_NOT_CURRENT"
    assert r["source_flux_removed"] is False


def test_unadopted_plane_sky_flux_response_uses_correct_reference_origin():
    hashes={"image":"1"*64,"invvar":"2"*64,
            "nexp":"3"*64,"maskbits":"4"*64}
    key=["dr11","south","1554p180","r"]
    fixed_center=[155.4,18.]
    stages = {
        "scene_masks":{"state":"CURRENT_LOADED","native_fits":{
            "source_quartet_key":key,
            "grid_fingerprint":"wcs-native-grid",
            "original_fits_sha256_local":hashes,
            "loaded_target_icrs":fixed_center,
            "image_unit_evidence":"FITS_IMAGE_BUNIT_NANOMAGGIES_CONFIRMED"}},
        "background": {"state": "CURRENT_DIAGNOSTIC", "result": {
            "band": "r", "project_id": "p", "draft_sha256": "a",
            "catalog_sha256": "b", "dataset_source": "coadd-r",
            "source_quartet_key":key,
            "grid_fingerprint":"wcs-native-grid",
            "original_fits_sha256_local":hashes,
            "image_unit_evidence":"FITS_IMAGE_BUNIT_NANOMAGGIES_CONFIRMED",
            "loaded_project_target_icrs":fixed_center,
            "source_fits_hashes_verified_at_scene_load":True,
            "methods": {
                "constant": {"fit_status": "FITTED",
                             "coefficients_nMgy": [0.002, 0., 0.],
                             "reference_xy_pix": [100., 200.]},
                "tilted_plane": {"fit_status": "FITTED",
                                 "coefficients_nMgy": [.002, .0001, -.0002],
                                 "reference_xy_pix": [100., 200.]},
            }}},
        "aperture": {"state": "CURRENT_DIAGNOSTIC", "result": {
            "band": "r", "project_id": "p", "scene_draft_sha256": "a",
            "catalog_sha256": "b", "dataset_source": "coadd-r",
            "dataset_id": "dr11/south/test/r",
            "source_quartet_key":key,
            "grid_fingerprint":"wcs-native-grid",
            "original_fits_sha256_local":hashes,
            "image_unit_evidence":"FITS_IMAGE_BUNIT_NANOMAGGIES_CONFIRMED",
            "loaded_project_target_icrs":fixed_center,
            "target_anchor_icrs":fixed_center,
            "source_fits_hashes_verified":True,
            "systematics_04_3": {
                "reference_semimajor_arcsec": 100.,
                "sky_response": {
                    "reference_native_xy": [105., 195.],
                    "dflux_db0_native_pixel2": -1000.,
                    "dflux_dbx_native_pixel3": -500.,
                    "dflux_dby_native_pixel3": 250.}}}}
    }
    response = _conditional_sky_aperture_response_qa(stages)
    assert response["state"] == "FITTED_SKY_COUNTERFACTUALS_AVAILABLE"
    assert response["methods"]["constant"][
        "conditional_flux_change_if_subtracted_nMgy"] == pytest.approx(-2.)
    at_target = .002 + .0001*5 + (-.0002)*(-5)
    expected = -1000*at_target + (-500)*.0001 + 250*(-.0002)
    assert response["methods"]["tilted_plane"][
        "conditional_flux_change_if_subtracted_nMgy"] == pytest.approx(
            expected)
    assert response["sky_uncertainty_inferred"] is False
    assert response["correction_applied"] is False
    assert response["science_adopted"] is False
    stages["conditional_sky_aperture_response"] = response
    assert "APERTURE_SENSITIVE_TO_UNADOPTED_SKY_FITS" in (
        _diagnostic_findings(stages, {}))
    stages["background"]["result"]["draft_sha256"] = "other"
    assert _conditional_sky_aperture_response_qa(stages)["state"] == (
        "SCENE_OR_BAND_PROVENANCE_MISMATCH")
    stages["background"]["result"]["draft_sha256"]="a"
    ap=stages["aperture"]["result"]
    ap["original_fits_sha256_local"]=dict(hashes,image="f"*64)
    mismatch=_conditional_sky_aperture_response_qa(stages)
    assert mismatch["state"]=="SCENE_OR_BAND_PROVENANCE_MISMATCH"
    assert mismatch["mismatch_key"]=="original_fits_sha256_local"
    ap["original_fits_sha256_local"]=hashes
    ap["grid_fingerprint"]="different-wcs"
    assert _conditional_sky_aperture_response_qa(stages)["mismatch_key"] == (
        "grid_fingerprint")
    ap["grid_fingerprint"]="wcs-native-grid"
    stages["scene_masks"]["native_fits"]["grid_fingerprint"]="other"
    assert _conditional_sky_aperture_response_qa(stages)["mismatch_key"] == (
        "full_native_source_or_target")
    stages["scene_masks"]["native_fits"]["grid_fingerprint"]="wcs-native-grid"
    ap["target_anchor_icrs"]=[155.401,18.]
    assert _conditional_sky_aperture_response_qa(stages)["mismatch_key"] == (
        "full_native_source_or_target")
    ap["target_anchor_icrs"]=fixed_center
    stages["scene_masks"]["native_fits"]["original_fits_sha256_local"]=None
    assert _conditional_sky_aperture_response_qa(stages)["state"] == (
        "SCENE_OR_BAND_PROVENANCE_MISMATCH")


def test_sky_flux_counterfactual_refuses_different_mask_raster_recipe():
    hashes={k:str(i)*64 for i,k in enumerate(
        ("image","invvar","nexp","maskbits"),1)}
    key=["dr11","south","1554p180","r"]
    anchor=[155.4,18.]
    policy="PROVISIONAL_SCENE_EFFECTIVE_AUTO_PLUS_EXCLUDE_V2_TAN_WCS"
    sky={"band":"r","project_id":"p","draft_sha256":"a",
         "catalog_sha256":"b","dataset_source":"coadd-r",
         "source_quartet_key":key,"grid_fingerprint":"wcs",
         "original_fits_sha256_local":hashes,
         "image_unit_evidence":"FITS_IMAGE_BUNIT_NANOMAGGIES_CONFIRMED",
         "loaded_project_target_icrs":anchor,
         "source_fits_hashes_verified_at_scene_load":True,
         "source_mask_policy_id":policy,
         "methods":{}}
    ap={"band":"r","project_id":"p","scene_draft_sha256":"a",
        "catalog_sha256":"b","dataset_source":"coadd-r",
        "source_quartet_key":key,"grid_fingerprint":"wcs",
        "original_fits_sha256_local":hashes,
        "image_unit_evidence":"FITS_IMAGE_BUNIT_NANOMAGGIES_CONFIRMED",
        "loaded_project_target_icrs":anchor,
        "target_anchor_icrs":anchor,
        "source_fits_hashes_verified":True,
        "source_policy_id":policy}
    stage={
        "scene_masks":{"state":"CURRENT_LOADED",
                       "native_fits":{
                           "source_quartet_key":key,
                           "grid_fingerprint":"wcs",
                           "original_fits_sha256_local":hashes,
                           "loaded_target_icrs":anchor,
                           "image_unit_evidence":
                               "FITS_IMAGE_BUNIT_NANOMAGGIES_CONFIRMED",
                           "original_fits_provenance_status":
                               "FULL_ORIGINAL_FITS_SHA256_VERIFIED_AT_SCENE_LOAD"},
                       "contaminant_screening":{
                           "source_mask_policy_id":policy}},
        "background":{"state":"CURRENT_DIAGNOSTIC","result":sky},
        "aperture":{"state":"CURRENT_DIAGNOSTIC","result":ap}}
    assert _cross_stage_source_lineage_qa(stage)["state"]==(
        "CURRENT_AVAILABLE_STAGES_MATCH_ORIGINAL_SCENE")
    assert _conditional_sky_aperture_response_qa(stage)["state"]==(
        "NO_EXACT_APERTURE_RESPONSE")
    ap["source_policy_id"]="PROVISIONAL_SCENE_EFFECTIVE_AUTO_PLUS_EXCLUDE_V1"
    mismatch=_cross_stage_source_lineage_qa(stage)
    assert mismatch["state"]=="CROSS_STAGE_SOURCE_IDENTITY_MISMATCH"
    assert not mismatch["stage_comparisons"]["aperture"][
        "source_mask_geometry_policy_matches"]
    resp=_conditional_sky_aperture_response_qa(stage)
    assert resp["state"]=="SCENE_OR_BAND_PROVENANCE_MISMATCH"
    assert resp["mismatch_key"]=="source_mask_geometry_policy"
    ap["source_policy_id"]=policy
    source_sha="a"*64
    manual_sha="b"*64
    screening=stage["scene_masks"]["contaminant_screening"]
    screening["effective_source_mask_sha256"]=source_sha
    screening["manual_exclude_mask_sha256"]=manual_sha
    for result in (sky,ap):
        result["effective_source_mask_sha256"]=source_sha
        result["manual_exclude_mask_sha256"]=manual_sha
    assert _cross_stage_source_lineage_qa(stage)["state"]==(
        "CURRENT_AVAILABLE_STAGES_MATCH_ORIGINAL_SCENE")
    assert _conditional_sky_aperture_response_qa(stage)["state"]==(
        "NO_EXACT_APERTURE_RESPONSE")
    ap["effective_source_mask_sha256"]="e"*64
    mismatch=_cross_stage_source_lineage_qa(stage)
    assert mismatch["state"]=="CROSS_STAGE_SOURCE_IDENTITY_MISMATCH"
    assert mismatch["stage_comparisons"]["aperture"][
        "native_mask_fingerprint_matches"]["effective_source_mask_sha256"] is False
    stage["cross_stage_source_lineage"]=mismatch
    findings=_diagnostic_findings(stage,{})
    assert "CROSS_STAGE_SOURCE_OR_MASK_LINEAGE_MISMATCH" in findings
    assert "CROSS_STAGE_NATIVE_MASK_PIXEL_FINGERPRINT_MISMATCH" in findings
    assert "CROSS_STAGE_SOURCE_FITS_PROVENANCE_MISMATCH" not in findings
    withheld=_conditional_sky_aperture_response_qa(stage)
    assert withheld["state"]=="SCENE_OR_BAND_PROVENANCE_MISMATCH"
    assert withheld["mismatch_key"]=="effective_source_mask_sha256"
    ap["effective_source_mask_sha256"]=source_sha
    sky["manual_exclude_mask_sha256"]="0"*64
    assert _conditional_sky_aperture_response_qa(stage)["mismatch_key"]==(
        "manual_exclude_mask_sha256")


def test_psf_band_navigation_never_discards_loaded_scene(qtbot, tmp_path):
    widget = _bound_workspace(qtbot, tmp_path)
    scene = widget.scene_editor
    # The minimal project fixture has no physical FITS or PSF inputs,
    # so initialize the fixed-band selectors without opening any files.
    for band, panel in widget.psf0_panels.items():
        panel.band.clear()
        panel.band.addItem(band.upper(), band)
    scene.set_band("r")
    scene._loaded = {"band": "r", "source": "native-r"}
    # The initially selected PSF tab is G; Scene native remains R.
    assert widget.band_tabs.currentWidget().band.currentData() == "g"
    widget.stage_tabs.setCurrentWidget(scene)
    assert scene._band == "r"
    assert scene._loaded["source"] == "native-r"
    widget.stage_tabs.setCurrentWidget(widget.stage_tabs.widget(0))
    widget.band_tabs.setCurrentIndex(2)
    assert widget.band_tabs.currentWidget().band.currentData() == "i"
    widget.stage_tabs.setCurrentWidget(scene)
    assert scene._band == "r"
    assert scene._loaded["source"] == "native-r"
    # In an empty Scene the PSF selector is still allowed to initialize
    # its band, retaining the original first-visit convenience.
    scene._loaded = None
    widget.stage_tabs.setCurrentWidget(widget.stage_tabs.widget(0))
    widget.stage_tabs.setCurrentWidget(scene)
    assert scene._band == "i"


def test_cross_stage_source_lineage_qa_rejects_foreign_fits_and_units():
    local_hashes={
        "image":"1"*64, "invvar":"2"*64,
        "nexp":"3"*64, "maskbits":"4"*64}
    key=["dr11","south","1554p180","r"]
    scene={"state":"CURRENT_LOADED","native_fits":{
        "source_quartet_key":key,
        "grid_fingerprint":"native-wcs-001",
        "loaded_target_icrs":[155.438012,18.088025],
        "original_fits_provenance_status":
            "FULL_ORIGINAL_FITS_SHA256_VERIFIED_AT_SCENE_LOAD",
        "original_fits_sha256_local":local_hashes,
        "image_unit_evidence":
            "FITS_IMAGE_BUNIT_NANOMAGGIES_CONFIRMED"}}
    target={"source_quartet_key":list(key),
            "grid_fingerprint":"native-wcs-001",
            "loaded_project_target_icrs":[155.438012,18.088025],
            "original_fits_sha256_local":dict(local_hashes),
            "image_unit_evidence":
                "FITS_IMAGE_BUNIT_NANOMAGGIES_CONFIRMED"}
    stages={"scene_masks":scene,
            "background":{"state":"CURRENT_DIAGNOSTIC",
                          "result":dict(target)},
            "aperture":{"state":"CURRENT_DIAGNOSTIC",
                        "result":dict(target,
                            target_anchor_icrs=[155.438012,18.088025])}}
    result=_cross_stage_source_lineage_qa(stages)
    assert result["state"]=="CURRENT_AVAILABLE_STAGES_MATCH_ORIGINAL_SCENE"
    assert all(row["state"]=="MATCH" for row in
               result["stage_comparisons"].values())
    assert result["science_adopted"] is False
    stages["aperture"]["result"]=dict(
        target,target_anchor_icrs=[155.438012,18.09])
    wrong_center=_cross_stage_source_lineage_qa(stages)
    assert wrong_center["state"]=="CROSS_STAGE_SOURCE_IDENTITY_MISMATCH"
    assert wrong_center["stage_comparisons"]["aperture"][
        "project_target_icrs_matches"] is False
    # A stage cannot be relabelled as the original simply because its
    # target display name, brick text or band agrees.
    stages["aperture"]["result"]=dict(target,
        original_fits_sha256_local=dict(local_hashes,image="f"*64))
    mismatch=_cross_stage_source_lineage_qa(stages)
    assert mismatch["state"]=="CROSS_STAGE_SOURCE_IDENTITY_MISMATCH"
    assert mismatch["inconsistent_stages"]==["aperture"]
    stages["cross_stage_source_lineage"]=mismatch
    assert "CROSS_STAGE_SOURCE_FITS_PROVENANCE_MISMATCH" in (
        _diagnostic_findings(stages, {}))
    stages["aperture"]["result"]=dict(target,image_unit_evidence="ADU")
    assert _cross_stage_source_lineage_qa(stages)["state"] == (
        "CROSS_STAGE_SOURCE_IDENTITY_MISMATCH")
    stages["aperture"]["result"]=dict(
        target,grid_fingerprint="different-wcs-grid")
    bad_grid=_cross_stage_source_lineage_qa(stages)
    assert bad_grid["state"]=="CROSS_STAGE_SOURCE_IDENTITY_MISMATCH"
    assert bad_grid["stage_comparisons"]["aperture"][
        "native_grid_fingerprint_matches"] is False
    stages["background"]["state"]="NOT_RUN"
    stages["aperture"]["state"]="NOT_RUN"
    assert _cross_stage_source_lineage_qa(stages)["state"] == (
        "NO_CURRENT_DEPENDENT_STAGE_TO_COMPARE")
    scene["native_fits"]["original_fits_sha256_local"]=None
    assert _cross_stage_source_lineage_qa(stages)["state"] == (
        "SCENE_SOURCE_NOT_FULLY_VERIFIED")
    scene["native_fits"]["original_fits_sha256_local"]=dict(
        local_hashes,image="non-sha-hex")
    assert _cross_stage_source_lineage_qa(stages)["state"] == (
        "SCENE_SOURCE_NOT_FULLY_VERIFIED")
