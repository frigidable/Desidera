"""Synthetic validation of residual-background fitting and spatial holdouts."""
import numpy as np
import pytest

from lsbpipe.scene import MaskLayer, LayerKind, SceneState
from lsbpipe.background import (
    BackgroundSettings, FitStatus, HoldoutStatus, fit_residual_background)


def make_scene(shape, train, hold, extra=()):
    grid = 'synthetic-grid-A'
    layers = [
        MaskLayer('training', LayerKind.SKY_TRAIN, train, grid, 'r'),
        MaskLayer('holdout', LayerKind.SKY_HOLDOUT, hold, grid, 'r')]
    for i, (kind, mask) in enumerate(extra):
        layers.append(MaskLayer(f'extra_{i}', kind, mask, grid, 'r'))
    return SceneState('galaxy', grid, 'r', shape, tuple(layers))


def regions(shape=(48, 52)):
    yy, xx = np.indices(shape)
    train = (yy % 5 != 0) & (xx % 6 != 0)
    hold = ~train
    target = (xx - 25)**2 + (yy - 23)**2 < 7**2
    train &= ~target
    hold &= ~target
    return train, hold


def fit(image, scene, method='constant', settings=None, ivar=None):
    return fit_residual_background(
        image, np.ones_like(image) if ivar is None else ivar,
        scene, method=method, grid_id='synthetic-grid-A', settings=settings)


def test_robust_constant_rejects_outliers_without_mutating_data():
    rng = np.random.default_rng(123)
    shape = (48, 52)
    train, hold = regions(shape)
    image = .008 + rng.normal(0., .001, shape)
    image[train & (np.indices(shape)[1] % 29 == 0)] += .5
    original = image.copy()
    result = fit(image, make_scene(shape, train, hold))
    assert result.fit_status == FitStatus.FITTED
    assert result.holdout_status == HoldoutStatus.MEASURED_UNCALIBRATED
    assert result.coefficients_nMgy[0] == pytest.approx(.008, abs=.001)
    assert result.coefficients_nMgy[1:] == (0., 0.)
    assert result.training_pixels > 100 and result.holdout_pixels > 100
    assert len(result.tiles) == 4
    assert np.array_equal(image, original)


def test_tilted_plane_restores_exact_gradient_and_outperforms_constant():
    shape = (55, 53)
    yy, xx = np.indices(shape)
    train, hold = regions(shape)
    image = .012 + .0008*(xx-26) - .0004*(yy-27)
    scene = make_scene(shape, train, hold)
    plane = fit(image, scene, 'tilted_plane')
    const = fit(image, scene, 'constant')
    assert plane.fit_status == FitStatus.FITTED
    assert plane.coefficients_nMgy == pytest.approx((.012, .0008, -.0004), abs=1e-10)
    assert plane.holdout_rms_residual_nMgy < 1e-10
    assert const.holdout_rms_residual_nMgy > .005
    assert np.max(np.abs(plane.image(shape) - image)) < 1e-10


def test_user_defined_holdout_threshold_flags_unfit_model():
    shape = (48, 52)
    yy, xx = np.indices(shape)
    train, hold = regions(shape)
    image = .02 + .002*(xx-25)
    scene = make_scene(shape, train, hold)
    settings = BackgroundSettings(max_abs_holdout_tile_median_nMgy=.003)
    assert fit(image, scene, 'constant', settings).holdout_status == (
        HoldoutStatus.REJECTED_BY_USER_THRESHOLD)
    assert fit(image, scene, 'tilted_plane', settings).holdout_status == (
        HoldoutStatus.WITHIN_USER_THRESHOLD)


def test_insufficient_and_collinear_sky_are_flagged_without_model():
    shape = (30, 30)
    train = np.zeros(shape, bool)
    train[3, 3] = True
    hold = np.zeros(shape, bool)
    hold[19:23, 2:27] = True
    result = fit(np.ones(shape), make_scene(shape, train, hold))
    assert result.fit_status == FitStatus.INSUFFICIENT_SKY
    with pytest.raises(ValueError, match='unsuccessful'):
        result.image(shape)
    train = np.zeros(shape, bool)
    train[3, 2:28] = True
    result = fit(np.ones(shape), make_scene(shape, train, hold), 'tilted_plane')
    assert result.fit_status == FitStatus.ILL_CONDITIONED
    assert result.design_condition > BackgroundSettings().max_design_condition


def test_target_protection_and_neighbor_are_excluded_from_sky_samples():
    shape = (48, 52)
    train, hold = regions(shape)
    target = np.zeros(shape, bool)
    target[:3, :] = True
    neighbor = np.zeros(shape, bool)
    neighbor[40:42, 10:13] = True
    scene = make_scene(shape, train, hold, [
        (LayerKind.TARGET_PROTECT, target),
        (LayerKind.NEIGHBOR, neighbor)])
    image = np.full(shape, .01)
    image[target | neighbor] = 500.
    result = fit(image, scene)
    assert result.coefficients_nMgy[0] == pytest.approx(.01, abs=1e-10)
    assert result.training_pixels < np.count_nonzero(train)
    with pytest.raises(ValueError, match='grid_id'):
        fit_residual_background(image, np.ones_like(image), scene,
                                method='constant', grid_id='wrong')


def test_single_quadrant_holdout_does_not_pass_spatial_qa():
    shape = (30, 30)
    train = np.zeros(shape, bool)
    train[:15, :] = True
    hold = np.zeros(shape, bool)
    hold[2:8, 2:9] = True
    train[hold] = False
    result = fit(np.ones(shape), make_scene(shape, train, hold))
    assert result.fit_status == FitStatus.FITTED
    assert result.holdout_status == HoldoutStatus.INSUFFICIENT_SPATIAL_COVERAGE


def test_no_valid_holdout_pixels_are_reported_not_adopted():
    shape = (48, 52)
    train, hold = regions(shape)
    ivar = np.ones(shape)
    ivar[hold] = 0.
    result = fit(np.ones(shape), make_scene(shape, train, hold), ivar=ivar)
    assert result.fit_status == FitStatus.FITTED
    assert result.holdout_pixels == 0
    assert result.holdout_status == HoldoutStatus.INSUFFICIENT_HOLDOUT
    assert result.holdout_median_residual_nMgy is None


def test_missing_explicit_sky_regions_fails_closed():
    scene = SceneState('galaxy', 'synthetic-grid-A', 'r', (10, 11))
    with pytest.raises(ValueError, match='Explicit sky_train'):
        fit(np.ones((10, 11)), scene)


def test_invalid_numerical_config_is_rejected():
    for options in ({'huber_t': 0}, {'max_design_condition': -1},
                    {'min_holdout_quadrants': 5},
                    {'max_abs_holdout_tile_median_nMgy': np.nan}):
        with pytest.raises(ValueError):
            BackgroundSettings(**options)



def _stage03_project_fixture(*, with_regions=True):
    from astropy.wcs import WCS
    from lsbpipe.photometry.sky_masks import (
        SkyMaskDraft, SkyMaskEdit, SkyMaskShape)
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [60.5, 60.5]
    wcs.wcs.crval = [155., 18.]
    wcs.wcs.cd = np.array([[-1., 0.], [0., 1.]]) / 3600.
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    shape = (120, 120)
    yy, xx = np.indices(shape, dtype=float)
    sky = .012 + .00045*(xx-59.5) - .00032*(yy-59.5)
    image = sky.copy()
    ivar = np.ones(shape)*1e6
    hard = np.zeros(shape, bool)
    hard[14:17, 13:16] = True
    auto = np.zeros(shape, bool)
    auto[100:111,100:111] = True
    protected = (xx-59.5)**2+(yy-59.5)**2 < 15.**2
    draft = SkyMaskDraft("stage03-project")
    if with_regions:
        train = ((14,14),(105,14),(14,105),(105,105),(60,60))
        hold = ((15,45),(105,45),(15,77),(105,77))
        for i, (purpose, locations) in enumerate((
                ("sky_train",train),("sky_holdout",hold))):
            for j,(x,y) in enumerate(locations):
                ra,dec = wcs.all_pix2world(float(x),float(y),0)
                shape_edit = SkyMaskShape(
                    "circle", (float(ra),float(dec)), radius_arcsec=8.)
                draft = draft.with_edit(SkyMaskEdit(
                    f"{purpose}_{i}_{j}",purpose,shape_edit,band="r"))
    dataset = {
        "project_id":"stage03-project",
        "band":"r","source":"dr11/south/test-brick/r",
        "pixel_shape":shape,"wcs":wcs,"native_hard":hard,
        "data":{"image":image,"invvar":ivar},
    }
    contamination = {
        "status":"CATALOGUE_CONTAMINATION_MASK_NOT_ADOPTED",
        "project_id":"stage03-project","band":"r",
        "dataset_source":dataset["source"],"draft_digest":draft.digest,
        "catalog_sha256":"a"*64,
        "auto_mask":auto,"protected":protected,"hard":hard.copy(),
        "gaia_inner_preview_enabled":True,
        "red_extended_preview_enabled":False,
        "source_inventory":{"entries":[]},
        "anchor_icrs":tuple(map(float,wcs.all_pix2world(59.5,59.5,0))),
    }
    return dataset, draft, contamination


def test_project_stage03_two_robust_sky_models_share_holdout_and_never_adopt():
    from lsbpipe.background import diagnose_current_project_residual_sky
    dataset,draft,contamination=_stage03_project_fixture()
    originals = (dataset["data"]["image"].copy(),
                 dataset["native_hard"].copy(),
                 contamination["auto_mask"].copy(),
                 contamination["protected"].copy())
    out=diagnose_current_project_residual_sky(dataset,draft,contamination)
    assert out["status"]=="AUTOMATIC_SKY_CANDIDATES_MEASURED_NOT_ADOPTED"
    assert out["schema"]=="desidera.sky.project-native-residual-diagnostic/2"
    assert out["automatic_sky_selection_qa"][
        "sky_regions_drawn_by_user_required"] is False
    assert out["automatic_sky_selection_qa"][
        "geometry_sufficient_for_diagnostic"] is True
    assert len(out["automatic_sky_selection_qa"]["holdout_quadrants"])>=3
    assert out["provisional_automatic_model"]["method"]=="tilted_plane"
    baseline = out["zero_additional_correction_baseline"]
    assert baseline["method"] == "zero_additional_correction_DR11_coadd"
    assert baseline["same_original_holdout_as_fitted_candidates"] is True
    assert baseline["n_holdout_samples"] == out["sky_holdout_effective_pixels"]
    assert baseline["n_training_samples"] == out["sky_train_effective_pixels"]
    assert baseline["holdout_rms_nMgy_per_pixel"] > 0.
    assert out["model_improvement_vs_zero_baseline"]["tilted_plane"][
        "rms_improvement_fraction"] > 0.99
    assert out["recommended_science_action"] == (
        "NO_ADDITIONAL_SKY_SUBTRACTION_UNTIL_VALIDATED")
    assert out["automatic_sky_background_correction_authorized"] is False
    assert out["provisional_automatic_model"][
        "choice_only_between_two_fitted_nonzero_candidates"] is True
    qa = out["automatic_sky_spatial_qa"]
    assert qa["schema"] == "desidera.sky.automatic-tile-spatial-qa/1"
    assert qa["no_tile_reselection_after_inspecting_residuals"] is True
    assert qa["spatial_support_complete_for_science"] is False
    hold_tiles = qa["methods"]["tilted_plane"]["holdout"]["tiles"]
    assert len(hold_tiles) == out["automatic_sky_selection_qa"]["holdout_tiles"]
    assert sum(t["n_samples"] for t in hold_tiles) == (
        out["sky_holdout_effective_pixels"])
    assert max(abs(t["median_residual_nMgy_per_pixel"])
               for t in hold_tiles) < 1.e-9
    assert out["band"]=="r"
    assert out["project_id"]=="stage03-project"
    assert out["sky_train_effective_pixels"]>100
    assert out["sky_holdout_effective_pixels"]>100
    assert out["automatic_sky_selection_qa"]["training_tiles"]>=8
    assert out["automatic_sky_selection_qa"]["holdout_tiles"]>=4
    assert out["preview_gaia_inner_mask_used"] is True
    assert out["protected_target_pixels"]>100
    assert out["methods"]["constant"]["holdout_status"]==(
        "MEASURED_UNCALIBRATED")
    plane=out["methods"]["tilted_plane"]
    assert plane["fit_status"]=="FITTED"
    assert plane["holdout_status"]=="MEASURED_UNCALIBRATED"
    assert plane["coefficients_nMgy"]==pytest.approx(
        (.012,.00045,-.00032),abs=1.e-9)
    assert plane["holdout_rms_residual_nMgy"] < 1.e-9
    assert out["methods"]["constant"]["holdout_rms_residual_nMgy"] > .01
    assert out["method_comparison"][
        "both_methods_share_exact_same_sky_train_and_holdout"] is True
    assert out["science_status"]=="UNADOPTED_DIAGNOSTIC_ONLY"
    assert out["sky_model_adopted"] is False
    assert out["psf_wings_validated"] is False
    assert out["noise_covariance_validated"] is False
    assert out["additional_bricks_or_catalogues_used"] is False
    assert out["background_subtracted_image_created"] is False
    assert np.array_equal(dataset["data"]["image"],originals[0])
    assert np.array_equal(dataset["native_hard"],originals[1])
    assert np.array_equal(contamination["auto_mask"],originals[2])
    assert np.array_equal(contamination["protected"],originals[3])


def test_project_stage03_no_manual_sky_regions_still_builds_automatic_sky():
    from lsbpipe.background import diagnose_current_project_residual_sky
    dataset,draft,contamination=_stage03_project_fixture(with_regions=False)
    out=diagnose_current_project_residual_sky(dataset,draft,contamination)
    assert out["status"]=="AUTOMATIC_SKY_CANDIDATES_MEASURED_NOT_ADOPTED"
    assert out["methods"]["tilted_plane"]["fit_status"]=="FITTED"
    assert out["automatic_sky_selection_qa"][
        "manual_sky_annotations_ignored_not_deleted"] is True
    assert out["sky_train_effective_pixels"]>100
    assert out["sky_holdout_effective_pixels"]>100
    assert out["sky_model_adopted"] is False


def test_project_stage03_rejects_stale_draft_and_bad_native_grid():
    from lsbpipe.background import diagnose_current_project_residual_sky
    dataset,draft,contamination=_stage03_project_fixture()
    wrong=dict(contamination,draft_digest="0"*64)
    with pytest.raises(ValueError, match="Stale"):
        diagnose_current_project_residual_sky(dataset,draft,wrong)
    wrong=dict(contamination,band="g")
    with pytest.raises(ValueError, match="Stale"):
        diagnose_current_project_residual_sky(dataset,draft,wrong)
    wrong=dict(contamination,hard=np.ones(dataset["pixel_shape"],bool))
    with pytest.raises(ValueError, match="native-hard"):
        diagnose_current_project_residual_sky(dataset,draft,wrong)
    wrong=dict(contamination,auto_mask=np.zeros((4,4),bool))
    with pytest.raises(ValueError, match="grid"):
        diagnose_current_project_residual_sky(dataset,draft,wrong)


def test_project_stage03_source_protection_hides_mistaken_sky_annotation():
    from lsbpipe.background import diagnose_current_project_residual_sky
    dataset,draft,contamination=_stage03_project_fixture()
    # One 8-arcsec sky circle is deliberately inside the protected galaxy.
    # It must contribute *zero* pixels to either numerical sky fit.
    result=diagnose_current_project_residual_sky(dataset,draft,contamination)
    # Removing the central training annotation must leave exactly the
    # same effective pixel footprints and both models unchanged.
    central=next(e.edit_id for e in draft.edits
                 if e.purpose=="sky_train" and
                 tuple(round(float(x)) for x in dataset["wcs"].all_world2pix(
                     *e.shape.center,0))==(60,60))
    clean=draft.without_edit(central)
    aligned=dict(contamination,draft_digest=clean.digest)
    other=diagnose_current_project_residual_sky(dataset,clean,aligned)
    assert result["sky_train_effective_pixels"]==other[
        "sky_train_effective_pixels"]
    assert result["methods"]["tilted_plane"]["coefficients_nMgy"] == (
        other["methods"]["tilted_plane"]["coefficients_nMgy"])
    assert result["method_comparison"] == other["method_comparison"]
    assert result["automatic_sky_selection_qa"] == other[
        "automatic_sky_selection_qa"]



def test_project_stage03_fails_closed_when_no_safe_automatic_sky_tiles():
    from lsbpipe.background import diagnose_current_project_residual_sky
    dataset,draft,contamination=_stage03_project_fixture(with_regions=False)
    bad=dict(contamination,auto_mask=np.ones(dataset["pixel_shape"],dtype=bool))
    out=diagnose_current_project_residual_sky(dataset,draft,bad)
    assert out["status"]=="INSUFFICIENT_AUTOMATIC_SCENE_SKY_GEOMETRY"
    assert out["zero_additional_correction_baseline"] is None
    assert out["automatic_sky_background_correction_authorized"] is False
    assert out["methods"]=={}
    assert out["provisional_automatic_model"] is None
    assert out["sky_model_adopted"] is False
    assert out["sky_train_effective_pixels"]==0


def test_project_stage03_requires_existing_bright_star_inventory():
    from lsbpipe.background import diagnose_current_project_residual_sky
    dataset,draft,contamination=_stage03_project_fixture(with_regions=False)
    bad=dict(contamination)
    bad.pop("source_inventory")
    out=diagnose_current_project_residual_sky(dataset,draft,bad)
    assert out["status"]=="SOURCE_INVENTORY_REQUIRED_FOR_AUTOMATIC_SKY"
    assert out["methods"]=={}
    assert out["provisional_automatic_model"] is None
    assert "BRIGHT_GAIA_SOURCE_INVENTORY_ABSENT" in out["warnings"]


def test_project_stage03_excludes_galaxy_leakage_and_existing_bright_gaia():
    from lsbpipe.background import diagnose_current_project_residual_sky
    dataset,draft,contamination=_stage03_project_fixture(with_regions=False)
    image=dataset["data"]["image"]
    yy,xx=np.indices(image.shape)
    target=(xx-59.5)**2+(yy-59.5)**2<=22.**2
    original=image.copy()
    image[target]+=50.0
    auto=contamination.copy()
    auto["source_inventory"]={"entries":[
        {"gaia_g_mag":14.5,"x":13.,"y":107.,"objid":25}]}
    out=diagnose_current_project_residual_sky(dataset,draft,auto)
    assert out["automatic_sky_selection_qa"][
        "n_existing_bright_gaia_sources_vetoed"]==1
    assert out["automatic_sky_selection_qa"]["bright_star_veto_pixels"]>100
    assert out["sky_model_adopted"] is False
    if out["methods"]:
        assert out["methods"]["tilted_plane"]["coefficients_nMgy"]==(
            pytest.approx((.012,.00045,-.00032),abs=1e-8))
    assert np.all(image[target]>original[target])
    assert np.array_equal(dataset["native_hard"],
                          contamination["hard"])


def test_project_stage03_automatic_selection_independent_of_image_brightness():
    from lsbpipe.background import diagnose_current_project_residual_sky
    dataset,draft,contamination=_stage03_project_fixture(with_regions=False)
    first=diagnose_current_project_residual_sky(dataset,draft,contamination)
    bright=dataset["data"]["image"].copy()
    yy,xx=np.indices(bright.shape)
    bright[30:90,30:90]+=200.0
    second=diagnose_current_project_residual_sky(
        dict(dataset,data=dict(dataset["data"],image=bright)),
        draft,contamination)
    assert first["automatic_sky_selection_qa"]==second[
        "automatic_sky_selection_qa"]
    assert first["sky_train_effective_pixels"]==second[
        "sky_train_effective_pixels"]
    assert first["sky_holdout_effective_pixels"]==second[
        "sky_holdout_effective_pixels"]
    assert first["automatic_sky_selection_qa"][
        "selection_predeclared_geometry_only_not_flux_optimized"] is True



def test_project_stage03_tile_qa_detects_geometric_cluster_and_rms_imbalance():
    from types import SimpleNamespace
    from lsbpipe.background import _audit_automatic_sky_tile_residuals
    img = np.zeros((96, 96), dtype=float)
    train = np.zeros(img.shape, dtype=bool)
    holdout = np.zeros(img.shape, dtype=bool)

    def make(tiles, mask):
        result = []
        for x, y in tiles:
            mask[y*32+5, x*32+5] = True
            result.append({
                "tile_index": [x, y],
                "quadrant": 0 if y == 0 else 2,
                "n_sampled_pixels": 1,
            })
        return result

    train_rows = make([(0,0),(1,0),(2,0)], train)
    holdout_rows = make([(0,1),(1,1),(2,1)], holdout)
    selection = {
        "tile_edge_native_pix": 32,
        "training_tile_evidence": train_rows,
        "holdout_tile_evidence": holdout_rows,
        "training_quadrants": [0],
        "holdout_quadrants": [2],
    }
    candidate = SimpleNamespace(
        fit_status=FitStatus.FITTED,
        training_rms_residual_nMgy=.05,
        holdout_rms_residual_nMgy=.005,
        evaluate_at=lambda x, y: np.zeros_like(x, dtype=float),
    )
    qa = _audit_automatic_sky_tile_residuals(
        img, train, holdout, selection, {"constant":candidate})
    assert qa["geometry"]["training"]["dominant_single_row_fraction"] == 1.
    assert qa["geometry"]["holdout"]["dominant_single_row_fraction"] == 1.
    assert qa["geometry"]["training"]["missing_quadrants"] == [1,2,3]
    assert qa["methods"]["constant"]["training_to_holdout_pixel_rms_ratio"] == 10.
    assert "TRAINING_TILES_CLUSTERED_IN_SINGLE_ROW" in qa["warning_flags"]
    assert "HOLDOUT_TILES_CLUSTERED_IN_SINGLE_ROW" in qa["warning_flags"]
    assert "CONSTANT_TRAINING_RMS_OVER_TWICE_HOLDOUT" in qa["warning_flags"]
    assert "TRAINING_QUADRANTS_NOT_ALL_REPRESENTED" in qa["warning_flags"]
    assert all(not row["n_partial_edge_tiles"]
               for row in qa["geometry"].values())
    assert len(qa["methods"]["constant"]["holdout"]["tiles"]) == 3
    assert qa["science_holds_even_without_warning"] is True
    assert not np.any(img)


def test_project_stage03_tile_qa_rejects_mismatched_sample_count():
    from types import SimpleNamespace
    from lsbpipe.background import _audit_automatic_sky_tile_residuals
    image = np.zeros((64,64))
    mask = np.zeros(image.shape, dtype=bool)
    mask[3,3] = True
    records = [{
        "tile_index":[0,0], "quadrant":0,"n_sampled_pixels":2,
    }]
    selection = {
        "tile_edge_native_pix":32,
        "training_tile_evidence":records,
        "holdout_tile_evidence":[],
        "training_quadrants":[0],
        "holdout_quadrants":[],
    }
    candidate = SimpleNamespace(
        fit_status=FitStatus.FITTED,
        training_rms_residual_nMgy=0.,
        holdout_rms_residual_nMgy=0.,
        evaluate_at=lambda x, y: np.zeros_like(x, dtype=float),
    )
    with pytest.raises(ValueError, match="tile geometry/provenance"):
        _audit_automatic_sky_tile_residuals(
            image, mask, np.zeros_like(mask), selection,
            {"constant":candidate})



def test_project_stage03_zero_baseline_is_exact_same_heldout_native_pixels():
    """The model-free DR11 coadd is always compared fairly with fitted sky."""
    from lsbpipe.background import (
        _audit_zero_sky_baseline, _automatic_scene_sky_pixels)
    from lsbpipe.photometry.sky_masks import evaluate_mask_draft
    dataset,draft,contamination=_stage03_project_fixture(with_regions=False)
    image = dataset["data"]["image"]
    masks = evaluate_mask_draft(
        draft,dataset["wcs"],dataset["pixel_shape"],"r",
        native_hard=dataset["native_hard"],
        automatic_contaminants=contamination["auto_mask"])
    train,held,selection = _automatic_scene_sky_pixels(
        dataset,masks,contamination)
    original = image.copy()
    baseline = _audit_zero_sky_baseline(image,train,held,selection)
    assert baseline["holdout_rms_nMgy_per_pixel"] == pytest.approx(
        np.sqrt(np.mean(image[held]**2)),abs=1e-13)
    assert baseline["training_rms_nMgy_per_pixel"] == pytest.approx(
        np.sqrt(np.mean(image[train]**2)),abs=1e-13)
    assert baseline["holdout_median_nMgy_per_pixel"] == pytest.approx(
        np.median(image[held]))
    assert sum(t["n_samples"] for t in
               baseline["holdout_tile_medians_nMgy_per_pixel"]) == held.sum()
    assert not baseline["background_subtracted_from_image"]
    assert not baseline["validates_original_sky_level"]
    assert np.array_equal(original,image)


def test_project_stage03_fit_preference_is_not_sky_subtraction_authorization():
    from lsbpipe.background import diagnose_current_project_residual_sky
    dataset,draft,contamination=_stage03_project_fixture(with_regions=False)
    # This sky image has a perfectly real non-zero plane. Even if the
    # numerical plane fits it exactly, this is not independent evidence
    # that the already sky-subtracted DESI coadd may be altered in science.
    out=diagnose_current_project_residual_sky(dataset,draft,contamination)
    assert out["provisional_automatic_model"]["method"] == "tilted_plane"
    assert out["provisional_automatic_model"][
        "scientifically_recommended_for_subtraction"] is False
    assert out["automatic_sky_background_correction_authorized"] is False
    assert out["sky_model_adopted"] is False
    assert out["background_subtracted_image_created"] is False
    assert out["zero_additional_correction_baseline"] is not None
    assert out["model_improvement_vs_zero_baseline"]["tilted_plane"][
        "improvement_is_not_statistical_significance"] is True
