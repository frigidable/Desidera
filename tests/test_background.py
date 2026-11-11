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
    assert out["status"]=="TWO_SKY_MODELS_FITTED_HOLDOUT_UNCALIBRATED"
    assert out["band"]=="r"
    assert out["project_id"]=="stage03-project"
    assert out["sky_train_effective_pixels"]>100
    assert out["sky_holdout_effective_pixels"]>100
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


def test_project_stage03_missing_explicit_regions_fails_closed():
    from lsbpipe.background import diagnose_current_project_residual_sky
    dataset,draft,contamination=_stage03_project_fixture(with_regions=False)
    out=diagnose_current_project_residual_sky(dataset,draft,contamination)
    assert out["status"]=="EXPLICIT_DISJOINT_SKY_REGIONS_REQUIRED"
    assert out["methods"]=={}
    assert out["method_comparison"] is None
    assert out["sky_train_effective_pixels"]==0
    assert out["sky_holdout_effective_pixels"]==0
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
                 tuple(round(x) for x in dataset["wcs"].all_world2pix(
                     *e.shape.center,0))==(60,60))
    clean=draft.without_edit(central)
    aligned=dict(contamination,draft_digest=clean.digest)
    other=diagnose_current_project_residual_sky(dataset,clean,aligned)
    assert result["sky_train_effective_pixels"]==other[
        "sky_train_effective_pixels"]
    assert result["methods"]["tilted_plane"]["coefficients_nMgy"] == (
        other["methods"]["tilted_plane"]["coefficients_nMgy"])
    assert result["method_comparison"] == other["method_comparison"]
