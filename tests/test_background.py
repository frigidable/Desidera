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
        "source_inventory":{
            "entries": [], "stored_count": 0,
            "total_relevant_sources": 0, "truncated": 0,
            "screening_evidence_incomplete": False,
        },
        "anchor_icrs":tuple(map(float,wcs.all_pix2world(59.5,59.5,0))),
    }
    return dataset, draft, contamination


def test_anisotropic_wcs_gaia_sky_veto_uses_true_sky_circle():
    """One angular exclusion must not become a wrong circular pixel mask."""
    from lsbpipe.background import (
        _sky_circle_native_patch, _automatic_scene_sky_pixels)
    from lsbpipe.photometry.sky_masks import evaluate_mask_draft
    jac=np.diag([-1.,.5])
    ys,xs,patch=_sky_circle_native_patch(
        (145,145),70.,70.,18.,jac)
    def included(x,y):
        # The native patch is intentionally bounded. Pixels outside its
        # conservative angular support are correctly excluded as False,
        # not indexed into a non-existent local array.
        if not (xs.start <= x < xs.stop and ys.start <= y < ys.stop):
            return False
        return bool(patch[y-ys.start,x-xs.start])
    assert included(70,70)
    assert included(70,102)   # dy=32 native px -> 16 arcsec
    assert not included(95,70) # dx=25 native px -> 25 arcsec
    assert not included(70,110) # dy=40 native px -> 20 arcsec
    with pytest.raises(ValueError,match="Degenerate celestial WCS"):
        _sky_circle_native_patch((145,145),70.,70.,18.,
                                 np.zeros((2,2)))

    dataset,draft,contamination=_stage03_project_fixture(
        with_regions=False)
    dataset["wcs"].wcs.cd=np.diag([-1.,.5])/3600.
    preview=evaluate_mask_draft(
        draft,dataset["wcs"],dataset["pixel_shape"],"r",
        native_hard=dataset["native_hard"],
        automatic_contaminants=contamination["auto_mask"])
    _,_,baseline=_automatic_scene_sky_pixels(
        dataset,preview,contamination)
    contamination["source_inventory"]["entries"]=[
        {"gaia_g_mag":17.5,"x":25.,"y":25.}]
    contamination["source_inventory"]["stored_count"] = 1
    contamination["source_inventory"]["total_relevant_sources"] = 1
    train,held,qa=_automatic_scene_sky_pixels(
        dataset,preview,contamination)
    assert qa["n_existing_bright_gaia_sources_vetoed"]==1
    assert qa["bright_star_veto_pixels"]>0
    assert baseline["bright_star_veto_pixels"]==0
    assert qa["bright_gaia_veto_uses_local_tangent_plane_wcs"] is True
    assert qa["guard_minimum_native_angular_scale_arcsec_per_pixel"] == (
        pytest.approx(.5))
    assert qa["native_pixel_scale_arcsec"] == (
        pytest.approx(np.sqrt(.5)))
    assert not np.any(train&held)
    assert not qa["galaxy_free_certified"]


def test_automatic_sky_refuses_empty_target_protection():
    """No target mask must never allow the galaxy to become sky training."""
    from lsbpipe.background import diagnose_current_project_residual_sky
    dataset,draft,contamination=_stage03_project_fixture(
        with_regions=False)
    mask_before=contamination["protected"].copy()
    image_before=dataset["data"]["image"].copy()
    contamination["protected"]=np.zeros(dataset["pixel_shape"],bool)
    with pytest.raises(ValueError,match="No protected target region"):
        diagnose_current_project_residual_sky(
            dataset,draft,contamination)
    assert np.array_equal(mask_before,mask_before.copy())
    assert np.array_equal(image_before,dataset["data"]["image"])
    assert not contamination["protected"].any()


def test_stage03_exact_mask_receipts_are_verified_before_sky_samples():
    from lsbpipe.background import diagnose_current_project_residual_sky
    from lsbpipe.photometry.sky_masks import native_boolean_mask_sha256
    dataset,draft,contamination=_stage03_project_fixture(with_regions=False)
    original_image=dataset["data"]["image"].copy()
    source_hash=native_boolean_mask_sha256(contamination["auto_mask"])
    manual_hash=native_boolean_mask_sha256(
        np.zeros(dataset["pixel_shape"],dtype=bool))
    contamination["effective_source_mask_sha256"]=source_hash
    contamination["manual_exclude_mask_sha256"]=manual_hash
    good=diagnose_current_project_residual_sky(
        dataset,draft,contamination)
    assert good["effective_source_mask_sha256"]==source_hash
    assert good["manual_exclude_mask_sha256"]==manual_hash
    assert good["sky_model_adopted"] is False
    contamination["manual_exclude_mask_sha256"]="f"*64
    with pytest.raises(ValueError,match="mask geometry fingerprint mismatch"):
        diagnose_current_project_residual_sky(dataset,draft,contamination)
    contamination["manual_exclude_mask_sha256"]=manual_hash
    changed=contamination["auto_mask"].copy()
    changed[5,5] = ~changed[5,5]
    contamination["auto_mask"]=changed
    with pytest.raises(ValueError,match="mask geometry fingerprint mismatch"):
        diagnose_current_project_residual_sky(dataset,draft,contamination)
    np.testing.assert_array_equal(dataset["data"]["image"],original_image)


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
    hull = out["target_training_geometry_qa"]
    assert hull["schema"] == "desidera.sky.target-training-convex-support/2"
    assert hull["n_hull_candidate_row_extrema"] <= 2*dataset["pixel_shape"][0]
    assert hull["exact_convex_hull_from_native_row_extrema"]
    assert hull["n_provisional_target_protected_pixels"] > 0
    assert hull["n_provisional_target_pixels_outside_training_hull"] == 0
    assert hull["provisional_target_footprint_fully_inside_training_hull"] is True
    assert hull["target_inside_training_pixel_convex_hull"] is True
    assert hull["status"] == "TARGET_INTERPOLATION_GEOMETRY_ONLY"
    assert hull["max_signed_supporting_line_distance_native_pix"] < 0
    assert hull["geometry_only_independent_of_image_flux"]
    assert hull["inside_hull_does_not_prove_unbiased_sky"]
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
    native_key=["dr11","south","1554p180","r"]
    dataset["source_quartet_key"]=native_key
    wrong=dict(contamination,source_quartet_key=[
        "dr11","south","1555p180","r"])
    with pytest.raises(ValueError, match="quartet identity mismatch"):
        diagnose_current_project_residual_sky(dataset,draft,wrong)
    correct=dict(contamination,source_quartet_key=native_key)
    output=diagnose_current_project_residual_sky(dataset,draft,correct)
    assert output["source_quartet_key"] == native_key
    hashes={key:str(i)*64 for i,key in enumerate(
        ("image","invvar","nexp","maskbits"),1)}
    dataset["original_fits_sha256_local"]=hashes
    dataset["original_fits_provenance_status"]=(
        "FULL_ORIGINAL_FITS_SHA256_VERIFIED_AT_SCENE_LOAD")
    dataset["image_unit_evidence"]="FITS_IMAGE_BUNIT_NANOMAGGIES_CONFIRMED"
    verified=diagnose_current_project_residual_sky(dataset,draft,correct)
    assert verified["source_fits_hashes_verified_at_scene_load"] is True
    assert verified["original_fits_sha256_local"] == hashes
    assert verified["sky_model_adopted"] is False
    dataset["image_unit_evidence"]=None
    with pytest.raises(ValueError,match="Scene FITS provenance"):
        diagnose_current_project_residual_sky(dataset,draft,correct)


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
    auto["source_inventory"]={
        "entries":[{"gaia_g_mag":14.5,"x":13.,"y":107.,"objid":25}],
        "stored_count":1,"total_relevant_sources":1,
        "truncated":0,"screening_evidence_incomplete":False,
    }
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
                "tile_bounds_native_xyxy": [x*32, y*32, (x+1)*32, (y+1)*32],
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
        "tile_index":[0,0],
        "tile_bounds_native_xyxy":[0,0,32,32],
        "quadrant":0,"n_sampled_pixels":2,
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



def test_shifted_sky_grid_spatial_qa_uses_true_native_pixel_bounds():
    """QA must not confuse grid indices with unshifted pixel origins."""
    from types import SimpleNamespace
    from lsbpipe.background import _audit_automatic_sky_tile_residuals

    image = np.zeros((96, 96), dtype=float)
    train = np.zeros_like(image, dtype=bool)
    holdout = np.zeros_like(image, dtype=bool)
    train[21, 25] = True
    holdout[21, 57] = True
    image[21, 25] = .125
    image[21, 57] = -.075
    # A tile grid offset (8, 12): grid index [1,1] begins at (24,20),
    # not at the naive (32,32).
    selection = {
        "tile_edge_native_pix": 32,
        "training_tile_evidence": [
            {"tile_index":[1,1],
             "tile_bounds_native_xyxy":[24,20,56,52],
             "n_sampled_pixels":1,
             "quadrant":0}],
        "holdout_tile_evidence": [
            {"tile_index":[2,1],
             "tile_bounds_native_xyxy":[56,20,88,52],
             "n_sampled_pixels":1,
             "quadrant":1}],
        "training_quadrants": [0],
        "holdout_quadrants": [1],
    }
    candidate = SimpleNamespace(
        fit_status=FitStatus.FITTED,
        training_rms_residual_nMgy=.125,
        holdout_rms_residual_nMgy=.075,
        evaluate_at=lambda x,y: np.zeros_like(x, dtype=float),
    )
    qa = _audit_automatic_sky_tile_residuals(
        image, train, holdout, selection, {"constant":candidate})
    assert qa["methods"]["constant"]["training"]["tiles"][0][
        "median_residual_nMgy_per_pixel"] == pytest.approx(.125)
    assert qa["methods"]["constant"]["holdout"]["tiles"][0][
        "median_residual_nMgy_per_pixel"] == pytest.approx(-.075)
    assert qa["geometry"]["training"]["n_partial_edge_tiles"] == 0
    assert qa["geometry"]["holdout"]["n_partial_edge_tiles"] == 0


def test_zero_sky_baseline_uses_shifted_native_tile_bounds_exactly():
    """A shifted holdout must not be reassigned to unshifted grid pixels."""
    from lsbpipe.background import _audit_zero_sky_baseline
    image=np.zeros((96,96),dtype=float)
    train=np.zeros((96,96),dtype=bool)
    heldout=np.zeros((96,96),dtype=bool)
    train[21,25]=True
    heldout[21,57]=True
    image[21,25]=.125
    image[21,57]=-.075
    selection={
        "tile_edge_native_pix":32,
        "holdout_tile_evidence":[{
            "tile_index":[2,1],
            "tile_bounds_native_xyxy":[56,20,88,52],
            "n_sampled_pixels":1,
            "quadrant":1}],
    }
    result=_audit_zero_sky_baseline(image,train,heldout,selection)
    assert result["holdout_tile_medians_nMgy_per_pixel"][0][
        "median_nMgy_per_pixel"]==pytest.approx(-.075)
    assert result["holdout_rms_nMgy_per_pixel"]==pytest.approx(.075)
    assert result["training_rms_nMgy_per_pixel"]==pytest.approx(.125)
    assert result["n_holdout_samples"]==1
    assert result["validates_original_sky_level"] is False
    selection["holdout_tile_evidence"][0][
        "tile_bounds_native_xyxy"]=[64,32,96,64]
    with pytest.raises(ValueError,match="sample lineage differs"):
        _audit_zero_sky_baseline(image,train,heldout,selection)


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



def test_project_stage03_fixed_grid_sensitivity_and_frozen_aperture_injections():
    """Automatic checks use the same Scene and never adopt sky or alter pixels."""
    from lsbpipe.background import diagnose_current_project_residual_sky
    dataset, draft, contamination = _stage03_project_fixture(with_regions=False)
    original_image = dataset["data"]["image"].copy()
    original_invvar = dataset["data"]["invvar"].copy()
    original_hard = dataset["native_hard"].copy()
    original_auto = contamination["auto_mask"].copy()
    report = diagnose_current_project_residual_sky(
        dataset, draft, contamination)
    assert report["status"] == "AUTOMATIC_SKY_CANDIDATES_MEASURED_NOT_ADOPTED"
    grid = report["sky_grid_sensitivity_qa"]
    assert grid["schema"] == "desidera.sky.internal-grid-sensitivity/1"
    assert grid["n_predeclared_variants"] == 6
    assert len(grid["trials"]) == 6
    assert grid["original_tile_selection_unchanged"] is True
    assert grid["no_science_adoption"] is True
    assert {r["name"] for r in grid["trials"]} == {
        "tile_72_origin", "tile_96_shift48", "tile_128_origin",
        "expanded_target_1p5", "expanded_sources_1p5",
        "spatial_cells_8x8_v2"}
    assert all(r["parameters"]["target_guard_factor"] >= 1.
               and r["parameters"]["source_guard_factor"] >= 1.
               for r in grid["trials"])
    assert all(r["no_reselection_based_on_pixel_brightness"]
               for r in grid["trials"])
    assert grid["summary"]["constant"]["not_an_uncertainty_interval"] is True
    assert grid["summary"]["tilted_plane"]["not_an_uncertainty_interval"] is True

    injection = report["sky_injection_recovery_qa"]
    assert injection["schema"] == (
        "desidera.sky.frozen-mask-exponential-injection/1")
    assert injection["status"] == (
        "SAME_COADD_FROZEN_MASK_INJECTION_DIAGNOSTIC_ONLY")
    assert len(injection["scenarios"]) == 3
    assert injection["not_full_galaxy_detections_or_photometric_pipeline"]
    assert injection["does_not_test_historic_DR11_sky_subtraction"] is True
    assert injection["no_noise_or_psf_realization_synthesized"] is True
    assert injection["original_arrays_mutated"] is False
    measured_medians = []
    for scenario in injection["scenarios"]:
        assert scenario["injected_flux_in_unmasked_aperture_nMgy"] > 0.
        leakage = scenario["injected_signal_on_frozen_sky_samples"]
        assert leakage["n_training_pixels"] == (
            report["sky_train_effective_pixels"])
        assert leakage["n_holdout_pixels"] == (
            report["sky_holdout_effective_pixels"])
        assert leakage["training_median_nMgy_per_pixel"] > 0.
        assert leakage["holdout_median_nMgy_per_pixel"] > 0.
        assert leakage["training_max_nMgy_per_pixel"] <= .005
        assert leakage["training_90th_percentile_nMgy_per_pixel"] <= .005
        assert leakage["training_signal_is_not_sky_or_noise"] is True
        measured_medians.append(leakage["training_median_nMgy_per_pixel"])
        assert scenario["zero_additional_correction_recovery_fraction"] == 1.
        assert scenario["unmasked_aperture_fraction"] <= 1.
        for method in ("constant", "tilted_plane"):
            candidate = scenario["models"][method]
            assert candidate["status"] == "INJECTED_FIT_MEASURED_UNCALIBRATED"
            assert np.isfinite(candidate[
                "injection_flux_recovery_fraction_fixed_aperture"])
            assert candidate["absolute_galaxy_flux_bias_identified"] is False
    assert measured_medians == sorted(measured_medians)
    assert len(set(measured_medians)) == 3
    assert report["automatic_sky_background_correction_authorized"] is False
    assert report["sky_model_adopted"] is False
    assert report["background_subtracted_image_created"] is False
    assert np.array_equal(dataset["data"]["image"], original_image)
    assert np.array_equal(dataset["data"]["invvar"], original_invvar)
    assert np.array_equal(dataset["native_hard"], original_hard)
    assert np.array_equal(contamination["auto_mask"], original_auto)
    assert (report["sky_grid_sensitivity_qa"] ==
            diagnose_current_project_residual_sky(
                dataset, draft, contamination)["sky_grid_sensitivity_qa"])


def test_project_stage03_sensitivity_selector_rejects_smaller_source_guards():
    from lsbpipe.background import _automatic_scene_sky_pixels
    from lsbpipe.photometry.sky_masks import evaluate_mask_draft
    dataset, draft, contamination = _stage03_project_fixture(with_regions=False)
    preview = evaluate_mask_draft(
        draft, dataset["wcs"], dataset["pixel_shape"], "r",
        native_hard=dataset["native_hard"],
        automatic_contaminants=contamination["auto_mask"])
    for options in (
            {"source_guard_factor":.5},
            {"target_guard_factor":.9},
            {"grid_offset_xy_native_pix":(96,0),
             "tile_edge_native_pix":96},
            {"tile_edge_native_pix":0}):
        with pytest.raises(ValueError):
            _automatic_scene_sky_pixels(
                dataset, preview, contamination, **options)
    a,b,base = _automatic_scene_sky_pixels(
        dataset,preview,contamination)
    x,y,expanded = _automatic_scene_sky_pixels(
        dataset,preview,contamination,target_guard_factor=1.5)
    assert expanded["target_outskirts_guard_native_pix"] >= (
        base["target_outskirts_guard_native_pix"])
    assert expanded["geometry_eligible_pixels_before_tile_sampling"] <= (
        base["geometry_eligible_pixels_before_tile_sampling"])
    # Downsampling to 384 pixels per tile can change membership even
    # when the enlarged exclusion removes no additional eligible tile.
    # The footprint must shrink, not the precise sparse subsample.
    assert not np.any(x & y)
    assert not np.any(a & b)


def test_project_stage03_frozen_mask_injection_fails_without_target_geometry():
    from lsbpipe.background import _audit_current_project_sky_injection_recovery
    image = np.zeros((80, 80), dtype=float)
    valid = np.ones_like(image)
    exclude = np.zeros(image.shape, dtype=bool)
    out = _audit_current_project_sky_injection_recovery(
        image, valid, None, {},
        anchor_native_xy=(40.,40.), protected=exclude,
        source_mask=exclude, manual_mask=exclude, hard=exclude,
        settings=None)
    assert out["status"] == "NOT_MEASURED_NO_TARGET_PROTECTION"
    assert out["science_adopted"] is False
    assert out["scenarios"] == []


def test_robust_sky_irls_iteration_cap_is_not_misreported_as_fitted():
    """One nonconverged Huber update cannot enter the sky comparison."""
    rng = np.random.default_rng(5031)
    shape = (48, 52)
    train, hold = regions(shape)
    img = .010 + rng.normal(0., .002, shape)
    yy, xx = np.indices(shape)
    img[train & (xx % 9 == 0)] += .15
    scene = make_scene(shape, train, hold)
    capped = BackgroundSettings(max_iterations=1, convergence_relative=1e-15)
    candidate = fit(img, scene, "constant", settings=capped)
    assert candidate.fit_status == FitStatus.NONCONVERGED
    assert candidate.holdout_status == HoldoutStatus.NOT_EVALUATED
    assert candidate.iterations == 1
    assert candidate.coefficients_nMgy is None
    with pytest.raises(ValueError, match="unsuccessful"):
        candidate.image(shape)
    # The same predeclared pixels may converge with a sufficient budget.
    resolved = fit(img, scene, "constant")
    assert resolved.fit_status == FitStatus.FITTED
    assert resolved.coefficients_nMgy is not None


def test_stage03_fails_closed_on_truncated_or_unverified_star_inventory():
    """An incomplete bright-star veto cannot produce a fitted sky plane."""
    from lsbpipe.background import diagnose_current_project_residual_sky
    dataset, draft, contamination = _stage03_project_fixture(
        with_regions=False)
    initial = diagnose_current_project_residual_sky(
        dataset, draft, contamination)
    assert initial["status"] == (
        "AUTOMATIC_SKY_CANDIDATES_MEASURED_NOT_ADOPTED")
    cases = [
        dict(truncated=2,total_relevant_sources=2),
        dict(screening_evidence_incomplete=True),
        dict(stored_count=5,total_relevant_sources=5),
        dict(total_relevant_sources=1),
        dict(truncated=None),
    ]
    for updates in cases:
        bad = dict(contamination)
        bad["source_inventory"] = dict(
            contamination["source_inventory"], **updates)
        report = diagnose_current_project_residual_sky(
            dataset, draft, bad)
        assert report["status"] == (
            "SOURCE_INVENTORY_INCOMPLETE_FOR_AUTOMATIC_SKY")
        assert report["methods"] == {}
        assert report["provisional_automatic_model"] is None
        assert report["sky_model_adopted"] is False
        assert report["automatic_sky_selection_qa"][
            "bright_gaia_inventory_complete"] is False


def test_residual_sky_settings_reject_boolean_fractional_and_string_thresholds():
    """Malformed sky QA controls must fail before optimizer/holdout gates."""
    invalid = (
        {"min_training_pixels": True},
        {"min_training_pixels": 12.5},
        {"min_holdout_pixels": False},
        {"min_holdout_pixels": 10.25},
        {"min_holdout_quadrants": True},
        {"min_holdout_quadrants": 2.5},
        {"max_iterations": True},
        {"max_iterations": 2.25},
        {"max_design_condition": True},
        {"max_design_condition": "10000"},
        {"huber_t": False},
        {"convergence_relative": "1e-10"},
        {"max_abs_holdout_tile_median_nMgy": True},
        {"max_abs_holdout_tile_median_nMgy": "0.01"},
    )
    for arguments in invalid:
        with pytest.raises(ValueError):
            BackgroundSettings(**arguments)
    settings = BackgroundSettings(
        min_training_pixels=np.int64(15),
        min_holdout_quadrants=np.int64(3),
        huber_t=np.float32(1.345),
        max_abs_holdout_tile_median_nMgy=np.float64(.01))
    assert type(settings.min_training_pixels) is int
    assert type(settings.min_holdout_quadrants) is int
    assert type(settings.huber_t) is float
    assert type(settings.max_abs_holdout_tile_median_nMgy) is float


def test_predeclared_spatial_sky_tile_sampling_balances_native_subcells():
    """Compare the unchanged row-major sampler with a 2-D fixed alternative."""
    from lsbpipe.background import _sample_native_sky_tile_pixels
    yy, xx = np.indices((80, 80))
    yy, xx = yy.ravel(), xx.ravel()
    kwargs = dict(footprint_shape=(80, 80),
                  tile_origin=(37, -14))
    a_y, a_x = _sample_native_sky_tile_pixels(
        yy, xx, policy="SPATIAL_8X8_PROPORTIONAL_V2", **kwargs)
    b_y, b_x = _sample_native_sky_tile_pixels(
        yy, xx, policy="SPATIAL_8X8_PROPORTIONAL_V2", **kwargs)
    assert len(a_x) == 384
    assert np.array_equal(a_x, b_x)
    assert np.array_equal(a_y, b_y)
    assert len(set(zip(a_y, a_x))) == 384
    cells = (np.minimum(7, a_y*8//80)*8 +
             np.minimum(7, a_x*8//80))
    assert np.bincount(cells, minlength=64).tolist() == [6]*64
    legacy_y, legacy_x = _sample_native_sky_tile_pixels(
        yy, xx, policy="ROW_MAJOR_UNIFORM_V1", **kwargs)
    expected = np.linspace(0, len(xx)-1, 384, dtype=np.int64)
    np.testing.assert_array_equal(legacy_y, yy[expected])
    np.testing.assert_array_equal(legacy_x, xx[expected])
    assert not np.array_equal(a_x, legacy_x)

    # An unevenly masked tile must still sample exactly 384 distinct
    # eligible positions without fabricating sky pixels or changing data.
    valid = ((yy > 20) & (xx < 53)) | ((yy < 16) & (xx > 60))
    uneven_y, uneven_x = _sample_native_sky_tile_pixels(
        yy[valid], xx[valid],
        policy="SPATIAL_8X8_PROPORTIONAL_V2", **kwargs)
    assert len(uneven_x) == 384
    assert len(set(zip(uneven_y, uneven_x))) == 384
    assert all(valid.reshape(80, 80)[y, x]
               for y, x in zip(uneven_y, uneven_x))
    with pytest.raises(ValueError, match="sampling policy"):
        _sample_native_sky_tile_pixels(
            yy, xx, policy="UNKNOWN_SKY_SAMPLER", **kwargs)


def test_stage03_spatial_sampler_is_alternative_not_default_selection():
    from lsbpipe.background import _automatic_scene_sky_pixels
    from lsbpipe.photometry.sky_masks import evaluate_mask_draft
    dataset, draft, contamination = _stage03_project_fixture(
        with_regions=False)
    preview = evaluate_mask_draft(
        draft, dataset["wcs"], dataset["pixel_shape"], "r",
        native_hard=dataset["native_hard"],
        automatic_contaminants=contamination["auto_mask"])
    _, _, baseline = _automatic_scene_sky_pixels(
        dataset, preview, contamination)
    _, _, alternative = _automatic_scene_sky_pixels(
        dataset, preview, contamination,
        sample_policy="SPATIAL_8X8_PROPORTIONAL_V2")
    assert baseline["sky_tile_subsampling_policy"] == (
        "ROW_MAJOR_UNIFORM_V1")
    assert alternative["sky_tile_subsampling_policy"] == (
        "SPATIAL_8X8_PROPORTIONAL_V2")
    assert baseline["sky_tile_subsampling_seed"] is None
    assert alternative["sky_tile_subsampling_seed"] == 20261113
    assert not alternative["sky_tile_subsampling_uses_original_flux"]
    assert alternative["sky_tile_stratification_cells_per_axis"] == 8
    assert baseline["training_tiles"] == alternative["training_tiles"]
    assert baseline["holdout_tiles"] == alternative["holdout_tiles"]
    with pytest.raises(ValueError, match="sampling policy"):
        _automatic_scene_sky_pixels(
            dataset, preview, contamination, sample_policy="BAD")


def test_stage03_incomplete_gaia_rows_fail_closed_without_silent_star_skip():
    """A numerically complete SC1 inventory must contain usable veto evidence."""
    from lsbpipe.background import diagnose_current_project_residual_sky
    dataset, draft, contamination = _stage03_project_fixture(
        with_regions=False)
    invalid_entries = (
        [None],
        [["not", "a", "source"]],
        [{"x": "invalid", "y": 20., "gaia_g_mag": 13.}],
        [{"x": np.nan, "y": 20., "gaia_g_mag": 13.}],
        [{"x": 20., "y": 20., "gaia_g_mag": True}],
        [{"x": 20., "y": 20., "gaia_g_mag": "13.0"}],
        [{"x": 20., "y": 20., "gaia_g_mag": None,
          "gaia_dr3_source_id": 123456789}],
    )
    for entries in invalid_entries:
        bad = dict(contamination)
        bad["source_inventory"] = {
            "entries": entries, "stored_count": 1,
            "total_relevant_sources": 1, "truncated": 0,
            "screening_evidence_incomplete": False,
        }
        report = diagnose_current_project_residual_sky(
            dataset, draft, bad)
        assert report["status"] == (
            "SOURCE_INVENTORY_INCOMPLETE_FOR_AUTOMATIC_SKY")
        assert report["methods"] == {}
        assert report["sky_model_adopted"] is False
        evidence = report["automatic_sky_selection_qa"]
        assert not evidence["bright_gaia_inventory_complete"]
        assert (evidence["source_inventory_invalid_or_unverifiable_rows"] +
                evidence["gaia_reference_sources_without_g_magnitude"]) >= 1


def test_stage03_zero_sky_reference_is_ranked_without_authorizing_subtraction():
    """The current sky-proposal heuristic must not hide a better zero baseline."""
    from lsbpipe.background import (
        diagnose_current_project_residual_sky,
        _automatic_scene_sky_pixels)
    from lsbpipe.photometry.sky_masks import evaluate_mask_draft
    dataset, draft, contamination = _stage03_project_fixture(
        with_regions=False)
    preview = evaluate_mask_draft(
        draft, dataset["wcs"], dataset["pixel_shape"], "r",
        native_hard=dataset["native_hard"],
        automatic_contaminants=contamination["auto_mask"])
    train, hold, _ = _automatic_scene_sky_pixels(
        dataset, preview, contamination)
    # This artificial field has a biased training footprint but an
    # unmodified, already perfectly sky-subtracted held-out footprint.
    # A fitted nonzero plane would *worsen* the held-out sky relative
    # to doing nothing; all masks and tile choices remain unchanged.
    image = np.zeros(dataset["pixel_shape"], dtype=float)
    image[train] = .025
    assert np.count_nonzero(hold) > 0
    dataset["data"]["image"] = image
    output = diagnose_current_project_residual_sky(
        dataset, draft, contamination)
    ranked = output["zero_inclusive_same_coadd_model_comparison"]
    assert ranked["lowest_rms_method_descriptive_only"] == (
        "zero_additional_correction")
    assert ranked["model_holdout_rms_nMgy_per_pixel"][
        "zero_additional_correction"] == 0.
    assert ranked["model_holdout_rms_nMgy_per_pixel"][
        "constant"] > 0.
    assert not ranked["science_sky_model_selected"]
    assert ranked["no_sky_subtraction_authorized"] is True
    assert output["automatic_sky_background_correction_authorized"] is False
    assert output["sky_model_adopted"] is False


def test_stage03_sky_veto_for_ambiguous_target_source_preserves_flux_mask():
    from lsbpipe.background import _automatic_scene_sky_pixels
    from lsbpipe.photometry.sky_masks import (
        evaluate_mask_draft, native_boolean_mask_sha256)
    dataset, draft, contamination = _stage03_project_fixture(with_regions=False)
    preview = evaluate_mask_draft(
        draft, dataset["wcs"], dataset["pixel_shape"], "r",
        native_hard=dataset["native_hard"],
        automatic_contaminants=contamination["auto_mask"])
    base_train, base_hold, before = _automatic_scene_sky_pixels(
        dataset, preview, contamination)
    yy, xx = np.indices(dataset["pixel_shape"])
    ambiguous = ((xx-94.)**2+(yy-59.5)**2 < 25.**2)
    amended = dict(contamination)
    amended["unresolved_target_overlap_sky_veto_mask"] = ambiguous
    amended["unresolved_target_overlap_sky_veto_sha256"] = (
        native_boolean_mask_sha256(ambiguous))
    train, hold, after = _automatic_scene_sky_pixels(
        dataset, preview, amended)
    assert after["unresolved_target_overlap_sky_veto_present"]
    assert after["unresolved_target_overlap_sky_veto_native_pixels"] > 0
    assert after["geometry_eligible_pixels_before_tile_sampling"] < (
        before["geometry_eligible_pixels_before_tile_sampling"])
    assert not np.any((train | hold) & ambiguous)
    assert np.array_equal(contamination["auto_mask"], amended["auto_mask"])
    assert after["unresolved_overlap_source_mask_and_manual_mask_unchanged"]
    bad = dict(amended)
    bad["unresolved_target_overlap_sky_veto_sha256"] = "0"*64
    with pytest.raises(ValueError, match="provenance mismatch"):
        _automatic_scene_sky_pixels(dataset, preview, bad)
    bad = dict(amended)
    bad["unresolved_target_overlap_sky_veto_mask"] = ambiguous.astype(np.uint8)
    with pytest.raises(ValueError, match="native Boolean"):
        _automatic_scene_sky_pixels(dataset, preview, bad)
    assert np.any(base_train | base_hold)


def test_stage03_ivar_weighted_heldout_rms_matches_direct_pixel_operator():
    """The optional noise-weighted comparison uses unchanged fixed pixels."""
    import json
    from lsbpipe.background import (
        diagnose_current_project_residual_sky,
        _automatic_scene_sky_pixels)
    from lsbpipe.photometry.sky_masks import evaluate_mask_draft
    dataset, draft, contamination = _stage03_project_fixture(with_regions=False)
    ivar = np.array(dataset["data"]["invvar"], copy=True)
    ivar[:, :ivar.shape[1]//2] *= .02
    dataset["data"]["invvar"] = ivar
    preview = evaluate_mask_draft(
        draft, dataset["wcs"], dataset["pixel_shape"], "r",
        native_hard=dataset["native_hard"],
        automatic_contaminants=contamination["auto_mask"])
    _, held, _ = _automatic_scene_sky_pixels(
        dataset, preview, contamination)
    report = diagnose_current_project_residual_sky(
        dataset, draft, contamination)
    assert report["status"] == "AUTOMATIC_SKY_CANDIDATES_MEASURED_NOT_ADOPTED"
    comparison = report["zero_inclusive_same_coadd_model_comparison"]
    assert comparison["weighted_qa_uses_exact_same_native_holdout_pixels"]
    assert comparison["weighted_qa_uses_original_diagonal_ivar_not_full_covariance"]
    assert comparison["weighted_qa_not_chi_squared_or_statistical_significance"]
    y, x = np.nonzero(held)
    w = ivar[y, x]
    samples = np.asarray(dataset["data"]["image"])[y, x]
    for method in ("zero_additional_correction", "constant", "tilted_plane"):
        residual = samples.copy()
        if method != "zero_additional_correction":
            row = report["methods"][method]
            b0, bx, by = row["coefficients_nMgy"]
            xref, yref = row["reference_xy_pix"]
            residual -= b0 + bx*(x-xref) + by*(y-yref)
        expected = float(np.sqrt(
            np.sum(w*residual**2) / np.sum(w)))
        assert comparison[
            "model_holdout_ivar_weighted_rms_nMgy_per_pixel"][method] == (
                pytest.approx(expected, rel=1e-12, abs=1e-12))
    assert report["sky_model_adopted"] is False
    assert report["automatic_sky_background_correction_authorized"] is False
    json.dumps(comparison, allow_nan=False)


def test_sky_target_convex_support_detects_extrapolation_independently_of_flux():
    """Three quadrants or a good fitted plane cannot certify spatial support."""
    from lsbpipe.background import _audit_sky_training_target_geometry

    train = np.zeros((60, 60), bool)
    train[10:16, 10:16] = True
    train[10:16, 30:36] = True
    train[30:36, 10:16] = True
    train[30:36, 30:36] = True
    original = train.copy()
    inside = _audit_sky_training_target_geometry(train, (22.5, 22.5))
    outside = _audit_sky_training_target_geometry(train, (51., 22.5))
    assert inside["status"] == "TARGET_INTERPOLATION_GEOMETRY_ONLY"
    assert inside["target_inside_training_pixel_convex_hull"] is True
    assert inside["max_signed_supporting_line_distance_native_pix"] < 0.
    assert outside["status"] == "TARGET_EXTRAPOLATION_GEOMETRY_ONLY"
    assert outside["target_inside_training_pixel_convex_hull"] is False
    assert outside["max_signed_supporting_line_distance_native_pix"] > 0.
    assert (inside["n_original_training_pixels"] ==
            outside["n_original_training_pixels"] == np.count_nonzero(train))
    assert outside["inside_hull_does_not_prove_unbiased_sky"]
    assert outside["no_sky_or_psf_correction_authorized"]
    np.testing.assert_array_equal(train, original)

    collinear = np.zeros((60, 60), bool)
    collinear[10, 10] = collinear[20, 20] = collinear[30, 30] = True
    invalid = _audit_sky_training_target_geometry(collinear, (20., 20.))
    assert invalid["status"] == "NO_NONDEGENERATE_SPATIAL_TRAINING_HULL"
    assert invalid["target_inside_training_pixel_convex_hull"] is None
    assert invalid["max_signed_supporting_line_distance_native_pix"] is None
    with pytest.raises(ValueError, match="Boolean"):
        _audit_sky_training_target_geometry(
            np.zeros((60, 60), dtype=float), (20., 20.))
    with pytest.raises(ValueError, match="finite native"):
        _audit_sky_training_target_geometry(train, (True, 22.))


def test_project_sky_report_flags_target_extrapolation_without_adopting(monkeypatch):
    """Keep numeric sky candidates but expose extrapolation of the target."""
    from lsbpipe.background import diagnose_current_project_residual_sky
    dataset, draft, contamination = _stage03_project_fixture(
        with_regions=False)
    # The actual Scene protector remains central, while this intentionally
    # inconsistent target-anchor alternative illustrates a displaced
    # galaxy against unchanged sky training/heldout samples.
    xfar, yfar = -10., 59.5
    ra, dec = dataset["wcs"].all_pix2world(xfar, yfar, 0)
    contamination["anchor_icrs"] = (float(ra), float(dec))
    out = diagnose_current_project_residual_sky(
        dataset, draft, contamination)
    assert out["status"] == "AUTOMATIC_SKY_CANDIDATES_MEASURED_NOT_ADOPTED"
    qa = out["target_training_geometry_qa"]
    assert qa["status"] == "TARGET_EXTRAPOLATION_GEOMETRY_ONLY"
    assert qa["target_inside_training_pixel_convex_hull"] is False
    assert "SKY_AT_TARGET_EXTRAPOLATED_OUTSIDE_TRAINING_CONVEX_HULL" in (
        out["warnings"])
    assert out["sky_model_adopted"] is False
    assert out["automatic_sky_background_correction_authorized"] is False



def test_sky_convex_hull_row_extrema_equal_full_pixel_hull_and_bounded_memory():
    """Every on-row interior point is a convex combination of its endpoints."""
    from scipy.spatial import ConvexHull
    from lsbpipe.background import _audit_sky_training_target_geometry

    rng = np.random.default_rng(309)
    train = rng.random((173, 251)) < .75
    yy, xx = np.nonzero(train)
    full = ConvexHull(np.column_stack((xx, yy)).astype(float))
    assert len(xx) > 30000
    for x, y in ((2., 3.), (120.5, 94.25), (250., 172.),
                 (-4., 88.), (275., 90.), (120., 210.)):
        result = _audit_sky_training_target_geometry(train, (x, y))
        direct = float(np.max(
            full.equations[:, :2] @ np.array([x, y])
            + full.equations[:, 2]))
        assert result["max_signed_supporting_line_distance_native_pix"] == (
            pytest.approx(direct, abs=1e-9))
        assert result["n_convex_hull_vertices"] == len(full.vertices)
        assert result["n_hull_candidate_row_extrema"] <= 2*train.shape[0]
        assert (result["n_hull_candidate_row_extrema"] <
                result["n_original_training_pixels"]/60)
        assert result["exact_convex_hull_from_native_row_extrema"]


def test_sky_protected_footprint_can_extrapolate_with_center_inside_hull():
    """Measuring the centre alone cannot certify a long/asymmetric galaxy."""
    from lsbpipe.background import _audit_sky_training_target_geometry

    train = np.zeros((80, 90), dtype=bool)
    for xs, ys in ((slice(9, 16), slice(10, 17)),
                   (slice(42, 49), slice(10, 17)),
                   (slice(9, 16), slice(45, 52)),
                   (slice(42, 49), slice(45, 52))):
        train[ys, xs] = True
    protected = np.zeros_like(train)
    protected[25:33, 25:40] = True       # Bracketed by sky train
    protected[25:33, 65:75] = True       # Target's outer asymmetry
    original = train.copy(), protected.copy()
    result = _audit_sky_training_target_geometry(
        train, (28.5, 29.), target_protection=protected)
    assert result["status"] == "TARGET_INTERPOLATION_GEOMETRY_ONLY"
    assert result["target_inside_training_pixel_convex_hull"] is True
    assert result["provisional_target_footprint_fully_inside_training_hull"] is False
    assert result["n_provisional_target_protected_pixels"] == 200
    assert result["n_provisional_target_pixels_outside_training_hull"] == 80
    assert result["fraction_provisional_target_pixels_outside_training_hull"] == (
        pytest.approx(.4))
    assert result["max_provisional_target_extrapolation_distance_native_pix"] > 0
    assert result["provisional_target_protection_not_a_measured_isophote"]
    assert result["protected_pixel_coverage_does_not_certify_actual_galaxy_extent"]
    np.testing.assert_array_equal(train, original[0])
    np.testing.assert_array_equal(protected, original[1])
    with pytest.raises(ValueError, match="same-grid Boolean"):
        _audit_sky_training_target_geometry(
            train, (28.5, 29.), target_protection=protected.astype(int))
    with pytest.raises(ValueError, match="same-grid Boolean"):
        _audit_sky_training_target_geometry(
            train, (28.5, 29.), target_protection=protected[:, :-1])


def test_sky_protected_footprint_null_for_unmeasured_or_collinear_training():
    from lsbpipe.background import _audit_sky_training_target_geometry

    train = np.zeros((60, 60), bool)
    train[5, 5] = train[10, 10] = train[15, 15] = True
    target = np.zeros_like(train)
    target[15:25, 25:35] = True
    report = _audit_sky_training_target_geometry(
        train, (27., 20.), target_protection=target)
    assert report["status"] == "NO_NONDEGENERATE_SPATIAL_TRAINING_HULL"
    assert report["n_provisional_target_protected_pixels"] == 100
    assert report["n_provisional_target_pixels_outside_training_hull"] is None
    assert report["fraction_provisional_target_pixels_outside_training_hull"] is None
    assert report["provisional_target_footprint_fully_inside_training_hull"] is None


def test_weighted_sky_lstsq_uses_its_own_singular_values_without_second_svd(
        monkeypatch):
    """LAPACK least squares returns the same conditioning evidence as SVD."""
    import lsbpipe.background as module

    rng = np.random.default_rng(20261113)
    x = np.linspace(-1., 1., 180)
    design = np.column_stack((np.ones_like(x), x, x*x))
    weights = rng.uniform(.05, 50., len(x))
    truth = np.array([.012, -.0008, .002])
    signal = design @ truth + rng.normal(0., 1.e-6, len(x))
    sqrtw = np.sqrt(weights / np.median(weights))
    original_singular = np.linalg.svd(
        design * sqrtw[:, None], compute_uv=False)
    expected_condition = original_singular[0]/original_singular[-1]

    def svd_must_not_be_called(*args, **kwargs):
        raise AssertionError("duplicate weighted design SVD")
    monkeypatch.setattr(module.np.linalg, "svd", svd_must_not_be_called)
    fitted, condition = module._weighted_lstsq(
        design, signal, weights, max_condition=1.e8)
    assert fitted == pytest.approx(
        np.linalg.lstsq(
            design*sqrtw[:, None], signal*sqrtw, rcond=None)[0],
        rel=1.e-11, abs=1.e-11)
    assert condition == pytest.approx(expected_condition, rel=1.e-11)
    rejected, poor_condition = module._weighted_lstsq(
        design, signal, weights, max_condition=1.01)
    assert rejected is None
    assert poor_condition == pytest.approx(expected_condition, rel=1.e-11)

    collinear = np.column_stack((np.ones_like(x), x, x))
    invalid, condition = module._weighted_lstsq(
        collinear, signal, weights, max_condition=1.e12)
    assert invalid is None
    assert condition > 1.e12



def test_frozen_sky_injection_broadcast_grid_matches_independent_dense_reference():
    """One reusable ellipse-distance grid is numerically identical to dense xy."""
    from lsbpipe.background import diagnose_current_project_residual_sky

    dataset, draft, contamination = _stage03_project_fixture(
        with_regions=False)
    frozen_image = dataset["data"]["image"].copy()
    report = diagnose_current_project_residual_sky(
        dataset, draft, contamination)
    qa = report["sky_injection_recovery_qa"]
    assert qa["status"] == (
        "SAME_COADD_FROZEN_MASK_INJECTION_DIAGNOSTIC_ONLY")
    shape = dataset["pixel_shape"]
    yy, xx = np.indices(shape, dtype=float)
    ax, ay = map(float, dataset["wcs"].all_world2pix(
        *contamination["anchor_icrs"], 0))
    target_radius = np.sqrt(
        np.count_nonzero(contamination["protected"])/np.pi)
    usable = (
        np.hypot(xx-ax, yy-ay) <= 2.*target_radius) & (
        ~dataset["native_hard"]) & (
        ~contamination["auto_mask"])
    for factor, scenario in zip((.5, 1., 2.), qa["scenarios"]):
        h = factor*target_radius
        independent_signal = .005*np.exp(
            -np.hypot(xx-ax, (yy-ay)/.65)/h)
        assert scenario["scale_length_native_pix"] == pytest.approx(h)
        assert scenario["injected_total_image_flux_nMgy"] == pytest.approx(
            float(independent_signal.sum()), rel=1e-12)
        assert scenario["injected_flux_in_unmasked_aperture_nMgy"] == (
            pytest.approx(float(independent_signal[usable].sum()),
                          rel=1e-12))
        assert scenario["unmasked_aperture_pixels"] == np.count_nonzero(usable)
    np.testing.assert_array_equal(dataset["data"]["image"], frozen_image)
    assert report["sky_model_adopted"] is False
    assert qa["original_arrays_mutated"] is False



def test_stage03_guard_cache_reuses_exact_baseline_and_rejects_mutated_masks(
        monkeypatch):
    """Six fixed grid trials need not repeat the same two full-grid EDTs."""
    import scipy.ndimage as ndimage
    from lsbpipe.background import _automatic_scene_sky_pixels
    from lsbpipe.photometry.sky_masks import evaluate_mask_draft

    dataset, draft, contamination = _stage03_project_fixture(
        with_regions=False)
    preview = evaluate_mask_draft(
        draft, dataset["wcs"], dataset["pixel_shape"], "r",
        native_hard=dataset["native_hard"],
        automatic_contaminants=contamination["auto_mask"])
    original_edt = ndimage.distance_transform_edt
    edt_calls = []
    def track_edt(*args, **kwargs):
        edt_calls.append(np.shape(args[0]))
        return original_edt(*args, **kwargs)
    monkeypatch.setattr(ndimage, "distance_transform_edt", track_edt)
    cache = {}
    variants = (
        {},
        {"tile_edge_native_pix": 72},
        {"tile_edge_native_pix": 96,
         "grid_offset_xy_native_pix": (48, 48)},
        {"tile_edge_native_pix": 128},
        {"target_guard_factor": 1.5},
        {"source_guard_factor": 1.5},
        {"sample_policy": "SPATIAL_8X8_PROPORTIONAL_V2"},
    )
    cached_outputs = [
        _automatic_scene_sky_pixels(
            dataset, preview, contamination,
            _guard_cache=cache, **opts)
        for opts in variants]
    assert len(edt_calls) == 4  # 2 baseline + 1 per expanded guard
    assert set(cache) == {"target", "catalogue"}
    assert all(not item[1].flags.writeable and
               not item[2].flags.writeable
               for item in cache.values())
    for opts, (train, hold, report) in zip(variants, cached_outputs):
        reference = _automatic_scene_sky_pixels(
            dataset, preview, contamination, **opts)
        np.testing.assert_array_equal(train, reference[0])
        np.testing.assert_array_equal(hold, reference[1])
        assert report == reference[2]
    assert len(edt_calls) == 4 + 2*len(variants)

    # Cached radii are not evidence of identical masks. An edit to a
    # source exclusion must trigger re-evaluation of the altered mask.
    before = len(edt_calls)
    preview["manual_exclude"][37:43, 40:46] = True
    changed = _automatic_scene_sky_pixels(
        dataset, preview, contamination, _guard_cache=cache)
    assert len(edt_calls) == before + 1
    direct = _automatic_scene_sky_pixels(
        dataset, preview, contamination)
    np.testing.assert_array_equal(changed[0], direct[0])
    np.testing.assert_array_equal(changed[1], direct[1])
    assert changed[2] == direct[2]

    # The protected target is a separate cache input with its own
    # exact-pixel invalidation, independent of centre/area estimates.
    contamination["protected"][60:63, 83:86] = True
    before = len(edt_calls)
    changed = _automatic_scene_sky_pixels(
        dataset, preview, contamination, _guard_cache=cache)
    assert len(edt_calls) == before + 2  # changed target AND source
    direct = _automatic_scene_sky_pixels(
        dataset, preview, contamination)
    np.testing.assert_array_equal(changed[0], direct[0])
    np.testing.assert_array_equal(changed[1], direct[1])
    assert changed[2] == direct[2]
