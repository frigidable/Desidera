"""P15 morphology must be real DR11 catalog metadata, never guessed."""
import numpy as np
import pytest

pytest.importorskip("astropy")
from astropy.wcs import WCS

from lsbpipe.photometry.psf_extended_scene import (
    audit_tractor_shape, local_tangent_jacobian,
    convolved_native_sersic, pilot_extended_scene,
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
    assert ring["signed_predictive_ring_rss_improvement_fraction"] > .99
    assert ring["both_fits_exclude_all_reserved_ring_pixels"]
    assert ring["pixel_resampling_covariance_across_ring_boundary_unmodeled"]
    assert ring["psf_itself_not_independently_trained_for_every_star"]
    assert ring["no_model_selection_or_science_adoption"]
    assert out["central_aperture_comparison"]["four_guessed_fwhm"][
        "same_original_native_pixels_globally_fit_models"]
    assert out["no_coadd_covariance_or_parameter_uncertainty_model"]
    assert out["profile_amplitudes_not_catalog_total_fluxes"]
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
