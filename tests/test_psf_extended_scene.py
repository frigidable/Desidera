"""P15 morphology must be real DR11 catalog metadata, never guessed."""
import numpy as np
import pytest

pytest.importorskip("astropy")
from astropy.wcs import WCS

from lsbpipe.photometry.psf_extended_scene import (
    audit_tractor_shape, local_tangent_jacobian,
    convolved_native_sersic, pilot_extended_scene,
    morphology_stress_variants, audit_approximate_sersic_sampling,
    compare_analytic_stamp_to_external_tractor,
    _morphology_stress_p15d,
)
from lsbpipe.photometry.psf_joint_scene import _shifted_native_kernel


def _wcs():
    w=WCS(naxis=2)
    w.wcs.crpix=[100.5,100.5]
    w.wcs.crval=[155.43,18.0]
    w.wcs.cdelt=[-.262/3600.,.262/3600.]
    w.wcs.ctype=["RA---TAN","DEC--TAN"]
    return w


def _kernel():
    yy,xx=np.indices((81,81),dtype=float)
    kernel=np.exp(-.5*((xx-40)**2+(yy-40)**2)/3.**2)
    return kernel/kernel.sum()


def _item(typ="DEV",**fields):
    row={"tractor_objid":1941,"tractor_type":typ,
         "shape_r_arcsec":1.6,"shape_e1":0.25,"shape_e2":0.,
         "flux_nMgy":5., "x_pix":118.6,"y_pix":100.4}
    row.update(fields)
    return row


def test_dr11_triple_of_ellipse_models_and_wcs_convention():
    w=_wcs()
    jac=local_tangent_jacobian(w,(100.,100.))
    assert jac[0,0]==pytest.approx(-.262,rel=2e-6)
    assert jac[1,1]==pytest.approx(.262,rel=2e-6)
    assert abs(jac[0,1])<1e-5
    assert abs(jac[1,0])<1e-5
    for typ,n in (("REX",1.),("EXP",1.),("DEV",4.)):
        morphology=audit_tractor_shape(_item(typ))
        assert morphology["status"]=="VALID_APPROXIMATE_EXTENDED_TEMPLATE"
        assert morphology["sersic_n"]==n
        assert morphology["shape_r_arcsec"]==1.6
        assert morphology["not_exact_tractor_renderer"]
    rex=audit_tractor_shape(_item("REX",shape_e1=.7,shape_e2=.4))
    assert rex["axis_ratio"]==1.
    assert rex["position_angle_east_of_north_deg"] is None
    exp=audit_tractor_shape(_item("EXP"))
    assert exp["axis_ratio"]==pytest.approx(.6)
    assert exp["position_angle_east_of_north_deg"]==pytest.approx(0.)
    east=audit_tractor_shape(_item("EXP",shape_e1=-.25))
    assert east["position_angle_east_of_north_deg"]==pytest.approx(90.)


def test_missing_invalid_and_unsupported_shapes_fail_closed():
    assert audit_tractor_shape(_item(shape_r_arcsec=None))["status"]==(
        "MISSING_OR_INVALID_HALF_LIGHT_RADIUS_ARCSEC")
    assert audit_tractor_shape(_item(shape_r_arcsec=-1.))["status"]==(
        "MISSING_OR_INVALID_HALF_LIGHT_RADIUS_ARCSEC")
    assert audit_tractor_shape(_item("DEV",shape_e1=None))["status"]==(
        "MISSING_TRACTOR_ELLIPTICITY")
    assert audit_tractor_shape(_item("EXP",shape_e1=1.))["status"]==(
        "INVALID_TRACTOR_ELLIPTICITY")
    assert audit_tractor_shape(_item("SER",sersic_n=2.))["status"]==(
        "UNSUPPORTED_TRACTOR_PROFILE_TYPE")
    assert audit_tractor_shape(_item("DUP"))["science_adopted"] is False


def test_convolution_preserves_orientation_and_finite_nonnegative_pixels():
    kernel=_kernel()
    w=_wcs()
    source=(100.,100.)
    north=convolved_native_sersic(
        kernel=kernel,wcs=w,source_xy=source,origin_xy=(20.,20.),
        width=161,morphology=audit_tractor_shape(_item(
            "EXP",shape_r_arcsec=2.,shape_e1=.4,shape_e2=0.)))
    east=convolved_native_sersic(
        kernel=kernel,wcs=w,source_xy=source,origin_xy=(20.,20.),
        width=161,morphology=audit_tractor_shape(_item(
            "EXP",shape_r_arcsec=2.,shape_e1=-.4,shape_e2=0.)))
    assert np.isfinite(north).all()
    assert np.min(north)>=0.
    assert north.shape==(161,161)
    # At the same pixel radius, the major axis has higher flux.
    assert north[93,80] > north[80,93]
    assert east[80,93] > east[93,80]
    round_rex=convolved_native_sersic(
        kernel=kernel,wcs=w,source_xy=source,origin_xy=(20.,20.),
        width=161,morphology=audit_tractor_shape(_item(
            "REX",shape_r_arcsec=2.,shape_e1=.4,shape_e2=.2)))
    assert round_rex[93,80]==pytest.approx(round_rex[80,93],rel=.02)
    dev=convolved_native_sersic(
        kernel=kernel,wcs=w,source_xy=source,origin_xy=(20.,20.),
        width=161,morphology=audit_tractor_shape(_item("DEV")))
    assert np.isfinite(dev).all()
    assert dev.max()>0.


def _truth_scene(include_extended=True,missing=False):
    kernel=_kernel()
    w=_wcs()
    yy,xx=np.indices((201,201),dtype=float)
    star=(100.25,99.65)
    source=_item("DEV",shape_r_arcsec=1.6,
                 x_pix=118.6,y_pix=100.4)
    image=.004+.000015*(xx-star[0])-.000008*(yy-star[1])
    image += 120.*_shifted_native_kernel(kernel,xx,yy,*star)
    if include_extended:
        gal=convolved_native_sersic(
            kernel=kernel,wcs=w,source_xy=(118.6,100.4),
            origin_xy=(20.,20.),width=161,
            morphology=audit_tractor_shape(source))
        # Insert the exact native test model on the same scene as pilot.
        image[20:181,20:181] += 7.*gal
    if missing:
        source["shape_r_arcsec"]=None
    return image,kernel,w,star,source


def test_extended_pilot_recovers_synthetic_model_and_does_not_mutate_pixels():
    image,k,w,star,gal=_truth_scene()
    original=image.copy()
    out=pilot_extended_scene(
        image,np.ones_like(image),np.ones_like(image,dtype=bool),
        wcs=w,kernel=k,target_xy=star,point_sources=[],
        extended_sources=[gal],guessed_fwhm_pix=6.57)
    assert out["status"]=="PARTIAL_APPROXIMATE_EXTENDED_SCENE_DIAGNOSTIC"
    assert out["n_extended_profiles_rendered"]==1
    assert out["n_extended_catalog_unmodeled"]==0
    assert out["fitted_extended_source_objids"]==[1941]
    assert out["n_extended_templates_with_positive_fitted_nnls_amplitude"]==1
    fit=out["fitted_extended_components"][0]
    assert fit["tractor_objid"]==1941
    assert fit["tractor_type"]=="DEV"
    assert not fit["nnls_amplitude_exactly_zero"]
    assert fit["nnls_fitted_peak_nMgy_per_native_pixel"]>0
    assert fit["fitted_model_sum_inside_4_guessed_fwhm_nMgy"]>0
    assert fit["finite_footprint_model_not_catalog_total_flux"]
    assert out["weighted_rss_with_extended"] < 1e-10
    assert out["weighted_rss_without_extended"] > 0.
    assert out["full_scene_fractional_rss_change"]>.99
    ring=out["predeclared_2_to_4_fwhm_ring_prediction"]
    assert ring["status"]=="PREDECLARED_RING_OUT_OF_SAMPLE_PIXEL_PREDICTION"
    assert ring["n_reserved_native_pixels"] > 0
    assert ring["n_fitting_native_pixels"] + ring["n_reserved_native_pixels"] == (
        out["n_valid_original_native_pixels"])
    assert ring["weighted_reserved_ring_rss_without_extended"] > 0
    assert ring["weighted_reserved_ring_rss_with_extended"] < 1e-10
    assert ring["reserved_ring_rss_ratio_extended_over_point_only"] < 1e-10
    assert not ring["extended_prediction_worse_than_point_only"]
    assert ring["signed_predictive_ring_rss_improvement_fraction"] > .99
    assert ring["both_fits_exclude_all_reserved_ring_pixels"]
    assert ring["pixel_resampling_covariance_across_ring_boundary_unmodeled"]
    assert ring["psf_itself_not_independently_trained_for_every_star"]
    assert ring["no_model_selection_or_science_adoption"]
    assert ring["ring_component_fits_use_only_nonreserved_original_pixels"]
    assert ring["n_extended_components_positive_with_ring_excluded"]==1
    assert ring["n_extended_components_full_positive_but_ring_excluded_zero"]==0
    assert ring["n_extended_components_full_zero_but_ring_excluded_positive"]==0
    assert len(ring["ring_excluded_extended_components"])==1
    split=ring["ring_excluded_extended_components"][0]
    assert split["tractor_objid"]==1941
    assert split["tractor_type"]=="DEV"
    assert split["design_matrix_column"]==1
    assert not split["ring_excluded_nnls_amplitude_exactly_zero"]
    assert split["ring_excluded_nnls_peak_amplitude_nMgy_per_pixel"]>0
    assert split["ring_predicted_component_sum_on_reserved_native_pixels_nMgy"]>0
    assert split["ring_excluded_to_full_peak_amplitude_ratio"]==pytest.approx(
        1.,rel=1e-7)
    assert split["split_fit_amplitudes_are_conditional_not_source_detection"]
    assert not split["full_scene_active_ring_excluded_inactive"]
    assert not split["full_scene_inactive_ring_excluded_active"]
    assert ring["ring_split_diagnostics_are_not_independent_star_holdouts"]
    assert out["central_aperture_comparison"]["four_guessed_fwhm"][
        "same_original_native_pixels_globally_fit_models"]
    assert out["no_coadd_covariance_or_parameter_uncertainty_model"]
    assert out["profile_amplitudes_not_catalog_total_fluxes"]
    p15e=out["p15e_analytic_subpixel_convergence"]
    assert p15e["status"]==(
        "THREE_GRID_ANALYTIC_CONVERGENCE_MEASURED_NOT_TRACTOR")
    assert p15e["selected_tractor_objid"]==1941
    assert p15e["evaluated_intrinsic_subpixel_grids"]==[3,5,7]
    assert p15e["n_original_valid_scene_pixels"]==161**2
    assert p15e["three_to_five"]["shape_l1_unit_flux_difference"]>0.
    assert p15e["five_to_seven"]["shape_l1_unit_flux_difference"]>0.
    assert p15e["not_a_comparison_against_tractor_mixture_of_gaussians"]
    assert p15e["no_automatic_threshold_or_best_grid_selection"]
    assert not p15e["science_adopted"]
    assert out["p15e_external_tractor_reference"]["status"]==(
        "ORIGINAL_MATCHED_TRACTOR_REFERENCE_NOT_SUPPLIED")
    assert out["science_adopted"] is False
    np.testing.assert_array_equal(image,original)


def test_unknown_original_shape_not_silently_fit_as_psf():
    image,k,w,star,source=_truth_scene(include_extended=False,missing=True)
    out=pilot_extended_scene(
        image,np.ones_like(image),np.ones_like(image,dtype=bool),
        wcs=w,kernel=k,target_xy=star,point_sources=[],
        extended_sources=[source],guessed_fwhm_pix=6.57)
    assert out["status"]=="NO_VALID_EXTENDED_SOURCE_PROFILE_TO_TEST"
    assert out["n_extended_catalog_unmodeled"]==1
    assert out["extended_source_morphology_audit"][0][
        "morphology"]["status"]=="MISSING_OR_INVALID_HALF_LIGHT_RADIUS_ARCSEC"
    assert not out["science_adopted"]


def test_scene_skips_unsupported_types_without_declaring_completeness():
    image,k,w,star,dev=_truth_scene()
    unknown=_item("SER",tractor_objid=2033,shape_r_arcsec=1.,
                  sersic_n=2.6,x_pix=140.,y_pix=120.)
    out=pilot_extended_scene(
        image,np.ones_like(image),np.ones_like(image,dtype=bool),
        wcs=w,kernel=k,target_xy=star,point_sources=[],
        extended_sources=[dev,unknown],guessed_fwhm_pix=6.57)
    assert out["status"]=="PARTIAL_APPROXIMATE_EXTENDED_SCENE_DIAGNOSTIC"
    assert out["n_extended_profiles_rendered"]==1
    assert out["n_extended_catalog_unmodeled"]==1
    assert not out["all_catalog_extended_profiles_accounted_for"]
    assert out["omitted_catalog_sources_prevent_complete_deblending"]



def test_reserved_ring_no_valid_pixels_reports_failure_without_claiming_cv():
    image,k,w,star,dev=_truth_scene()
    yy,xx=np.indices(image.shape,dtype=float)
    dist=np.hypot(xx-star[0],yy-star[1])
    valid=np.ones_like(image,dtype=bool)
    valid[(dist>2*6.57)&(dist<=4*6.57)]=False
    out=pilot_extended_scene(
        image,np.ones_like(image),valid,wcs=w,kernel=k,
        target_xy=star,point_sources=[],extended_sources=[dev],
        guessed_fwhm_pix=6.57)
    assert out["status"]=="PARTIAL_APPROXIMATE_EXTENDED_SCENE_DIAGNOSTIC"
    ring=out["predeclared_2_to_4_fwhm_ring_prediction"]
    assert ring["status"]=="INSUFFICIENT_RING_OR_TRAINING_PIXELS"
    assert ring["n_reserved_native_pixels"]==0
    assert ring["no_model_selection_or_science_adoption"]


def test_reserved_ring_partition_is_deterministic_under_preserved_native_inputs():
    image,k,w,star,gal=_truth_scene()
    kwargs=dict(wcs=w,kernel=k,target_xy=star,
                point_sources=[],extended_sources=[gal],
                guessed_fwhm_pix=6.57)
    a=pilot_extended_scene(
        image,np.ones_like(image),np.ones_like(image,dtype=bool),
        **kwargs)["predeclared_2_to_4_fwhm_ring_prediction"]
    b=pilot_extended_scene(
        image,np.ones_like(image),np.ones_like(image,dtype=bool),
        **kwargs)["predeclared_2_to_4_fwhm_ring_prediction"]
    assert a==b
    assert a["status"]=="PREDECLARED_RING_OUT_OF_SAMPLE_PIXEL_PREDICTION"



def test_injected_wrong_catalog_shape_worsens_reserved_ring_prediction():
    """An intentionally wrong DEV shape cannot mimic known exact truth."""
    image,k,w,star,truth=_truth_scene()
    valid=np.ones_like(image,dtype=bool)
    ivar=np.ones_like(image)
    good=pilot_extended_scene(
        image,ivar,valid,wcs=w,kernel=k,
        target_xy=star,point_sources=[],extended_sources=[truth],
        guessed_fwhm_pix=6.57)
    wrong=dict(truth,shape_r_arcsec=4.5,shape_e1=-.35,shape_e2=.15)
    bad=pilot_extended_scene(
        image,ivar,valid,wcs=w,kernel=k,
        target_xy=star,point_sources=[],extended_sources=[wrong],
        guessed_fwhm_pix=6.57)
    heldout_good=good["predeclared_2_to_4_fwhm_ring_prediction"]
    heldout_bad=bad["predeclared_2_to_4_fwhm_ring_prediction"]
    assert heldout_good["status"]==(
        "PREDECLARED_RING_OUT_OF_SAMPLE_PIXEL_PREDICTION")
    assert heldout_bad["status"]==(
        "PREDECLARED_RING_OUT_OF_SAMPLE_PIXEL_PREDICTION")
    assert heldout_good["weighted_reserved_ring_rss_with_extended"]<1e-10
    assert heldout_bad["weighted_reserved_ring_rss_with_extended"]>1e-9
    assert heldout_bad["n_reserved_native_pixels"]==(
        heldout_good["n_reserved_native_pixels"])
    assert bad["science_adopted"] is False


def test_large_shape_valid_catalog_fails_before_expensive_rendering():
    image,k,w,star,dev=_truth_scene(include_extended=False)
    candidates=[
        dict(dev,tractor_objid=3000+j,x_pix=110.+j*.1)
        for j in range(29)
    ]
    result=pilot_extended_scene(
        image,np.ones_like(image),np.ones_like(image,dtype=bool),
        wcs=w,kernel=k,target_xy=star,point_sources=[],
        extended_sources=candidates,guessed_fwhm_pix=6.57)
    assert result["status"]==(
        "POTENTIAL_EXTENDED_SOURCE_BUDGET_EXCEEDED_BEFORE_RENDER")
    assert result["n_potential_shape_valid_extended_sources"]==29
    assert result["source_budget"]==28
    assert result["no_partial_catalog_source_truncation"]
    assert not result["science_adopted"]



def test_p15d_morphology_cases_are_predeclared_and_keep_tractor_conventions():
    original=_item("DEV",shape_e1=.25,shape_e2=0.)
    copy=dict(original)
    cases=morphology_stress_variants(original)
    assert len(cases)==6
    assert all(x["status"]=="VALID_CONTROLLED_VARIANT" for x in cases)
    by_case={x["case"]:x for x in cases}
    smaller=by_case["half_light_radius_minus_10_percent"]["catalog_row"]
    larger=by_case["half_light_radius_plus_10_percent"]["catalog_row"]
    assert smaller["shape_r_arcsec"]==pytest.approx(1.44)
    assert larger["shape_r_arcsec"]==pytest.approx(1.76)
    e_minus=by_case["ellipticity_norm_minus_0p05"]["catalog_row"]
    e_plus=by_case["ellipticity_norm_plus_0p05"]["catalog_row"]
    assert np.hypot(e_minus["shape_e1"],e_minus["shape_e2"])==pytest.approx(.2)
    assert np.hypot(e_plus["shape_e1"],e_plus["shape_e2"])==pytest.approx(.3)
    pa_minus=by_case["position_angle_minus_10_deg"]["catalog_row"]
    pa_plus=by_case["position_angle_plus_10_deg"]["catalog_row"]
    assert audit_tractor_shape(pa_minus)[
        "position_angle_east_of_north_deg"]==pytest.approx(170.)
    assert audit_tractor_shape(pa_plus)[
        "position_angle_east_of_north_deg"]==pytest.approx(10.)
    assert np.hypot(pa_minus["shape_e1"],pa_minus["shape_e2"])==pytest.approx(.25)
    assert np.hypot(pa_plus["shape_e1"],pa_plus["shape_e2"])==pytest.approx(.25)
    assert original==copy  # input Tractor metadata is strictly read-only


def test_p15d_rex_has_radius_only_and_round_exp_has_undefined_pa():
    rex=_item("REX",shape_e1=None,shape_e2=None)
    cases=morphology_stress_variants(rex)
    assert len(cases)==6
    assert sum(x["status"]=="VALID_CONTROLLED_VARIANT" for x in cases)==2
    assert sum(x["status"]=="NOT_A_REX_PROFILE_PARAMETER" for x in cases)==4
    circle=morphology_stress_variants(_item("EXP",shape_e1=0.,shape_e2=0.))
    assert sum(x["status"]=="VALID_CONTROLLED_VARIANT" for x in circle)==2
    assert sum(x["status"]=="UNDEFINED_CIRCULAR_REFERENCE_PA" for x in circle)==4


def test_p15d_ellipticity_boundary_excluded_not_clamped():
    cases=morphology_stress_variants(_item(
        "DEV",shape_e1=.98,shape_e2=0.))
    ups=[x for x in cases if x["case"]=="ellipticity_norm_plus_0p05"]
    assert ups[0]["status"]=="OUTSIDE_PHYSICAL_ELLIPTICITY_DOMAIN"
    assert "catalog_row" not in ups[0]
    assert not morphology_stress_variants(_item("DEV",shape_r_arcsec=None))


def test_p15d_same_original_pixels_stress_response_on_exact_injected_dev():
    image,k,w,star,dev=_truth_scene()
    # Catalog inverse variances are separate input evidence, not the
    # physical tolerances chosen for our predeclared stress experiment.
    dev=dict(dev,shape_r_ivar_arcsec_m2=25.,
             shape_e1_ivar=100.,shape_e2_ivar=225.)
    original=image.copy()
    result=pilot_extended_scene(
        image,np.ones_like(image),np.ones_like(image,dtype=bool),
        wcs=w,kernel=k,target_xy=star,point_sources=[],
        extended_sources=[dev],guessed_fwhm_pix=6.57)
    stress=result["p15d_morphology_controlled_stress"]
    assert stress["status"]=="PREDECLARED_MORPHOLOGY_STRESS_DIAGNOSTIC"
    assert stress["schema"]=="desidera.psf0-p15d-tractor-morphology-stress/1"
    assert stress["selected_tractor_objid"]==1941
    assert stress["n_rendered_extended_sources_in_original_fit"]==1
    assert stress["source_selection_policy"]==(
        "NEAREST_RENDERABLE_EXTENDED_CENTROID_PRIOR_TO_FITTING")
    assert stress["n_successful_stress_cases"]==6
    assert stress["n_skipped_or_failed_cases"]==0
    assert stress["n_fitting_original_native_pixels"]+(
        stress["n_reserved_original_native_pixels"])==(
            result["n_valid_original_native_pixels"])
    assert stress["perturbations_are_not_catalog_sigma_or_uncertainty"]
    ivars=stress["nominal_catalog_formal_shape_ivars_not_covariance_or_stress_scale"]
    assert ivars["shape_r_ivar_arcsec_m2"]==25.
    assert ivars["shape_e1_ivar"]==100.
    assert ivars["shape_e2_ivar"]==225.
    assert stress["radius_fractional_stress"]==pytest.approx(.1)
    assert stress["ellipticity_norm_absolute_stress"]==pytest.approx(.05)
    assert stress["position_angle_degree_stress"]==pytest.approx(10.)
    assert stress["stress_envelope_is_not_a_confidence_interval"]
    assert stress["reserved_2_to_4_guessed_fwhm_never_used_in_fits"]
    assert stress["science_adopted"] is False
    assert stress["nominal_ring_rss"]<1e-10
    assert stress["nominal_sky_at_star_nMgy_per_native_pixel"]==pytest.approx(
        .004,abs=1e-8)
    for case in stress["cases"]:
        assert case["status"]=="CONTROLLED_PERTURBATION_FIT_SUCCESS"
        assert case["reserved_ring_weighted_rss"]>stress["nominal_ring_rss"]
        assert case["delta_reserved_ring_weighted_rss"]>0.
        assert case["ring_score_is_descriptive_not_noise_covariance_likelihood"]
    envelope=stress["descriptive_min_max_stress_response_not_error_bars"]
    assert envelope["delta_reserved_ring_weighted_rss"]["max"]>0.
    np.testing.assert_array_equal(image,original)


def test_p15d_preselected_nearest_source_has_stable_id_not_residual_ranking():
    image,k,w,star,truth=_truth_scene()
    far=_item(
        "EXP",tractor_objid=2033,x_pix=140.,y_pix=120.,flux_nMgy=3.)
    out=pilot_extended_scene(
        image,np.ones_like(image),np.ones_like(image,dtype=bool),
        wcs=w,kernel=k,target_xy=star,point_sources=[],
        extended_sources=[far,truth],guessed_fwhm_pix=6.57)
    stress=out["p15d_morphology_controlled_stress"]
    assert stress["selected_tractor_objid"]==1941
    assert stress["n_rendered_extended_sources_in_original_fit"]==2
    assert stress["only_one_extended_source_perturbed_at_a_time"]
    assert stress["all_other_rendered_catalog_components_frozen"]
    assert not stress["science_adopted"]



def test_p15e_same_native_morphology_grid_convergence_is_not_tractor_reference():
    kernel=_kernel()
    w=_wcs()
    morphology=audit_tractor_shape(_item("EXP",shape_r_arcsec=1.5))
    m=np.ones((81,81),dtype=bool)
    result=audit_approximate_sersic_sampling(
        kernel=kernel,wcs=w,source_xy=(100.33,99.66),
        origin_xy=(60.,60.),width=81,morphology=morphology,valid=m)
    assert result["status"]==(
        "THREE_GRID_ANALYTIC_CONVERGENCE_MEASURED_NOT_TRACTOR")
    assert result["evaluated_intrinsic_subpixel_grids"]==[3,5,7]
    assert result["same_catalog_morphology_and_native_wcs_all_grids"]
    assert result["same_scene_geometry_and_original_valid_pixel_mask"]
    for key in ("three_to_five","five_to_seven","three_to_seven"):
        diff=result[key]
        assert np.isfinite(diff["shape_l1_unit_flux_difference"])
        assert diff["shape_l1_unit_flux_difference"] >= 0.
        assert diff["relative_patch_normalization_not_source_total_flux"]
    again=audit_approximate_sersic_sampling(
        kernel=kernel,wcs=w,source_xy=(100.33,99.66),
        origin_xy=(60.,60.),width=81,morphology=morphology,valid=m)
    assert result==again
    assert result["not_a_comparison_against_tractor_mixture_of_gaussians"]
    assert not result["science_adopted"]


def test_p15e_malformed_reference_or_missing_provenance_cannot_claim_validation():
    kernel=_kernel()
    truth=np.ones((81,81),dtype=float)
    mask=np.ones_like(truth,dtype=bool)
    good=dict(
        renderer="dstndstn/tractor",
        renderer_version_or_commit="external-commit-mock",
        tractor_source_type="DEV",
        source_identifier="test-obj-1941",
        native_wcs_fingerprint="wcs-mock",
        input_psf_fingerprint="psf-mock",
        reference_generation_recipe="test fixture, not real tractor",
    )
    missing=compare_analytic_stamp_to_external_tractor(
        truth,truth,mask,provenance={})
    assert missing["status"]=="INCOMPLETE_EXTERNAL_REFERENCE_PROVENANCE"
    assert "native_wcs_fingerprint" in missing["missing_fields"]
    other=dict(good,renderer="unrelated-stamp-generator")
    refused=compare_analytic_stamp_to_external_tractor(
        truth,truth,mask,provenance=other)
    assert refused["status"]=="NOT_AN_IDENTIFIED_TRACTOR_RENDERER"
    mismatch=compare_analytic_stamp_to_external_tractor(
        truth,truth[:60],mask,provenance=good)
    assert mismatch["status"]=="INCOMPARABLE_EXTERNAL_REFERENCE_ARRAYS"
    # Even a pixel-perfect manually-supplied mock cannot be reported
    # as independent external Tractor validation by this API.
    numerical=compare_analytic_stamp_to_external_tractor(
        truth,truth,mask,provenance=good)
    assert numerical["status"]==(
        "EXTERNAL_REFERENCE_ARRAYS_COMPARED_NOT_SCIENCE_CERTIFIED")
    assert numerical["common_native_pixel_shape_comparison"][
        "shape_l1_unit_flux_difference"]==pytest.approx(0.)
    assert numerical["external_reference_identity_self_reported_not_verified"]
    assert numerical["still_requires_independent_psf_and_wcs_lineage_verification"]
    assert numerical["comparison_does_not_certify_tractor_reproduction"]
    assert numerical["science_adopted"] is False


def test_p15e_grid_validation_fails_closed_and_does_not_invent_reference():
    with pytest.raises(ValueError,match="UNSUPPORTED_DECLARED_SUBPIXEL_GRID"):
        convolved_native_sersic(
            kernel=_kernel(),wcs=_wcs(),
            source_xy=(100.,100.),origin_xy=(60.,60.),width=81,
            morphology=audit_tractor_shape(_item("DEV")),
            subpixel_grid=4)
    no_mask=audit_approximate_sersic_sampling(
        kernel=_kernel(),wcs=_wcs(),source_xy=(100.,100.),
        origin_xy=(60.,60.),width=81,
        morphology=audit_tractor_shape(_item("DEV")),
        valid=np.zeros((81,81),dtype=bool))
    assert no_mask["status"]=="INVALID_OR_EMPTY_NATIVE_REFERENCE_MASK"
    assert no_mask["science_adopted"] is False


def test_p15e_masked_stamps_are_never_compared_on_different_pixel_sets():
    scene=np.arange(81*81,dtype=float).reshape(81,81)+1.
    valid=np.ones_like(scene,dtype=bool)
    valid[:,40:]=False
    provided=dict(
        renderer="dstndstn/tractor",
        renderer_version_or_commit="test-only-not-real",
        tractor_source_type="EXP",source_identifier="mock",
        native_wcs_fingerprint="native-identity-test",
        input_psf_fingerprint="test-frozen-psf",
        reference_generation_recipe="synthetic fixture",
    )
    candidate=scene.copy()
    candidate[:,40:]*=1000.
    report=compare_analytic_stamp_to_external_tractor(
        candidate,scene,valid,provenance=provided)
    assert report["status"]=="EXTERNAL_REFERENCE_ARRAYS_COMPARED_NOT_SCIENCE_CERTIFIED"
    assert report["common_native_pixel_shape_comparison"][
        "shape_l1_unit_flux_difference"]==pytest.approx(0.)
    assert not report["science_adopted"]



def test_p15f_inactive_nominal_selected_source_stress_does_not_claim_robustness():
    # Construct an exact baseline star+plane signal with *negative*
    # extended-source-shaped mismatch. Its nonnegative NNLS coefficient
    # is exactly at the inactive boundary; geometry stresses then lack
    # evidence about the catalogued galaxy.
    yy,xx=np.indices((81,81),dtype=float)
    star=np.exp(-.5*((xx-40.)**2+(yy-40.)**2)/3.**2)
    neighbor=np.exp(-.5*((xx-59.)**2+(yy-40.)**2)/4.**2)
    P=np.column_stack((star.ravel(),neighbor.ravel()))
    B=np.column_stack((
        np.ones(81*81),
        ((xx-40.)/40.).ravel(),
        ((yy-40.)/40.).ravel()))
    z=(.2+25*star-3*neighbor).ravel()
    result=_morphology_stress_p15d(
        z=z,B=B,P=P,weights=np.ones_like(z),
        native_radius=np.hypot(xx-40.,yy-40.).ravel(),
        mask=np.ones((81,81),dtype=bool),
        kernel=_kernel(),wcs=_wcs(),
        scene_origin_xy=(0,0),scene_width_pix=81,
        guessed_fwhm_pix=6.57,
        source_row=_item("DEV",x_pix=59.,y_pix=40.),
        template_column=1,n_rendered_extended_sources=1)
    assert result["status"]==(
        "PREDECLARED_SELECTED_EXTENDED_SOURCE_NNLS_INACTIVE")
    assert result["selected_extended_template_peak_amplitude_nMgy_per_pixel"]==0
    assert not result["selected_extended_template_active_in_ring_excluded_nnls"]
    assert result["zero_sensitivity_from_inactive_source_not_evidence_of_stability"]
    assert result["n_successful_stress_cases"]==0
    assert result["n_skipped_or_failed_cases"]==6
    assert all(case["status"]==(
        "NOT_EVALUATED_SELECTED_EXTENDED_SOURCE_INACTIVE")
        for case in result["cases"])
    assert not result["science_adopted"]


def test_p15f_increased_grid_step_is_explicit_not_hidden():
    comparison=audit_approximate_sersic_sampling(
        kernel=_kernel(),wcs=_wcs(),source_xy=(40.06,40.01),
        origin_xy=(0,0),width=81,
        morphology=audit_tractor_shape(
            _item("DEV",shape_r_arcsec=.4)),
        valid=np.ones((81,81),dtype=bool))
    assert comparison["status"]==(
        "THREE_GRID_ANALYTIC_CONVERGENCE_MEASURED_NOT_TRACTOR")
    a=comparison["three_to_five"]["shape_l1_unit_flux_difference"]
    b=comparison["five_to_seven"]["shape_l1_unit_flux_difference"]
    assert comparison["l1_step_5_to_7_larger_than_3_to_5"]==(b>a)
    assert comparison["l1_step_ratio_5_to_7_over_3_to_5"]==pytest.approx(
        b/a)
    assert comparison["decreasing_l1_step_not_proof_of_convergence"]
    assert not comparison["science_adopted"]


def test_p15f_zero_or_worse_reserved_ring_score_reports_explicit_ratio():
    from lsbpipe.photometry.psf_extended_scene import (
        _predeclared_ring_prediction)
    yy,xx=np.indices((81,81),dtype=float)
    radii=np.hypot(xx-40.,yy-40.).ravel()
    B=np.column_stack((np.ones(81*81),
        ((xx-40.)/40.).ravel(),((yy-40.)/40.).ravel()))
    star=np.exp(-.5*((xx-40.)**2+(yy-40.)**2)/3.**2).ravel()
    ext=np.exp(-.5*((xx-61.)**2+(yy-40.)**2)/5.**2).ravel()
    P=np.column_stack((star,ext))
    # Native profile differs between the fitting region and ring.
    # This guarantees no unsupported universal "improved" claim.
    # The source is present only on the fit pixels, but its assumed
    # extrapolation has the wrong sign in reserved pixels.
    z=.2+20*star+5*ext
    ring=(radii>2*6.57)&(radii<=4*6.57)
    z[ring]-=10*ext[ring]
    a=_predeclared_ring_prediction(
        z,B,P,np.ones_like(z),radii,
        guessed_fwhm_pix=6.57,n_base_columns=1)
    assert a["status"]=="PREDECLARED_RING_OUT_OF_SAMPLE_PIXEL_PREDICTION"
    assert a["reserved_ring_rss_ratio_extended_over_point_only"]==pytest.approx(
        a["weighted_reserved_ring_rss_with_extended"] /
        a["weighted_reserved_ring_rss_without_extended"])
    assert a["extended_prediction_worse_than_point_only"]
    assert a["reserved_ring_rss_ratio_extended_over_point_only"]>1.
    assert a["worsened_ring_is_explicit_falsification_not_a_p_value"]



def test_p15g_inactive_ring_component_is_separate_from_full_native_fit():
    # Extended profile appears only in the reserved ring. A global
    # regression can fit it, but a fit that does not see the ring
    # must not claim that it measured its amplitude.
    from lsbpipe.photometry.psf_extended_scene import (
        _predeclared_ring_prediction, _projected_nnls)
    yy,xx=np.indices((81,81),dtype=float)
    radius=np.hypot(xx-40.,yy-40.).ravel()
    B=np.column_stack((np.ones(81*81),
                       ((xx-40.)/40.).ravel(),
                       ((yy-40.)/40.).ravel()))
    star=np.exp(-.5*((xx-40.)**2+(yy-40.)**2)/3.**2).ravel()
    galaxy=np.exp(-.5*((xx-60.)**2+(yy-40.)**2)/5.**2).ravel()
    P=np.column_stack((star,galaxy))
    reserved=(radius>2*6.57)&(radius<=4*6.57)
    z=.1+20*star
    z[reserved]+=4*galaxy[reserved]
    weights=np.ones_like(z)
    global_amplitudes,_,_,_=_projected_nnls(z,B,P,weights)
    out=_predeclared_ring_prediction(
        z,B,P,weights,radius,guessed_fwhm_pix=6.57,
        n_base_columns=1)
    assert out["status"]=="PREDECLARED_RING_OUT_OF_SAMPLE_PIXEL_PREDICTION"
    assert global_amplitudes[1]>0
    # Floating arithmetic may make an inactive boundary tiny
    # rather than identically zero, so test magnitude rather than
    # imposing an unphysical diagnostic threshold.
    assert out["ring_excluded_extended_components"][0][
        "ring_excluded_nnls_peak_amplitude_nMgy_per_pixel"]<1e-8
    assert out["ring_excluded_extended_components"][0][
        "ring_predicted_component_sum_on_reserved_native_pixels_nMgy"]<1e-7
    assert out["ring_component_fits_use_only_nonreserved_original_pixels"]
    assert out["n_reserved_native_pixels"]>0
    assert not out["no_model_selection_or_science_adoption"] is False
