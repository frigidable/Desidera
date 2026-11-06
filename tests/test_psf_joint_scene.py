"""Joint native-scene fitting must remain read-only, honest and non-adopting."""
from types import SimpleNamespace

import numpy as np
import pytest

from lsbpipe.photometry.psf_joint_scene import (
    _shifted_native_kernel, fit_fixed_psf_scene,
    audit_fixed_psf_joint_scenes, MAX_NUISANCE_SOURCES,
)


def _synthetic_field(width=201, sigma=3., *, center=(100.25, 99.65),
                     neighbor=(118.6, 100.4), source_flux=120.,
                     neighbor_flux=15., include_neighbor=True):
    ky, kx = np.indices((81, 81), dtype=float)
    kernel = np.exp(-.5*((kx-40)**2+(ky-40)**2)/sigma**2)
    kernel /= np.sum(kernel)
    yy, xx = np.indices((width, width), dtype=float)
    image = (.004 + .000015*(xx-center[0])
             - .000008*(yy-center[1]))
    image += source_flux * _shifted_native_kernel(
        kernel, xx, yy, center[0], center[1])
    if include_neighbor:
        image += neighbor_flux * _shifted_native_kernel(
            kernel, xx, yy, neighbor[0], neighbor[1])
    return image, kernel


def _source_neighbor(x=118.6, y=100.4):
    return {"x_pix":x, "y_pix":y, "tractor_objid":33,
            "tractor_type":"PSF", "flux_nMgy":15.}


def test_native_joint_scene_recovers_known_star_sky_neighbor_without_changing_input():
    image, kernel = _synthetic_field()
    original = image.copy()
    ivar = np.full_like(image, 1500.)
    mask = np.ones_like(image, dtype=bool)
    result = fit_fixed_psf_scene(
        image, ivar, mask, kernel=kernel,
        target_xy=(100.25,99.65),
        neighbors=[_source_neighbor()], guessed_fwhm_pix=6.57,
        original_sky_nMgy_pix=.001)
    assert result["status"] == "FIXED_PSF_JOINT_SCENE_DIAGNOSTIC"
    assert result["scene_width_pix"] == 161
    assert result["n_original_native_pixels"] == 161**2
    assert result["target_fitted_model_support_flux_nMgy"] == pytest.approx(
        120.,rel=2e-3)
    assert result["neighbors"][0]["fitted_flux_within_model_support_nMgy"] == pytest.approx(
        15.,rel=2e-3)
    assert result["fitted_sky_at_star_nMgy_per_native_pixel"] == pytest.approx(
        .004,abs=1e-8)
    assert result["fitted_sky_gradient_x_nMgy_per_native_pixel_per_pix"] == pytest.approx(
        .000015,abs=1e-9)
    assert result["fitted_sky_gradient_y_nMgy_per_native_pixel_per_pix"] == pytest.approx(
        -.000008,abs=1e-9)
    assert result["joint_minus_annulus_sky_nMgy_per_native_pixel"] == pytest.approx(
        .003,abs=1e-8)
    assert result["weighted_squared_residual_star_neighbors_plane"] < 1e-10
    assert result["weighted_squared_residual_star_plane_only"] > 0.
    assert result["weighted_residual_improvement_fraction"] > .99
    assert result["neighbor_to_central_core_model_flux_ratio"] is not None
    assert result["neighbor_to_central_4fwhm_model_flux_ratio"] > 0.
    assert result["correlated_noise_and_fit_parameter_covariance_unmodeled"]
    assert result["kernel_support_beyond_observed_ePSF_unmodeled"]
    assert not result["science_adopted"]
    np.testing.assert_array_equal(image,original)
    np.testing.assert_array_equal(ivar,1500.)


def test_joint_scene_empty_neighbors_recovers_same_sky_plane():
    image, kernel = _synthetic_field(include_neighbor=False)
    result = fit_fixed_psf_scene(
        image, np.ones_like(image), np.ones_like(image,dtype=bool),
        kernel=kernel, target_xy=(100.25,99.65),
        neighbors=[],guessed_fwhm_pix=6.57,
        original_sky_nMgy_pix=0.)
    assert result["status"] == "FIXED_PSF_JOINT_SCENE_DIAGNOSTIC"
    assert result["n_catalog_point_neighbors"] == 0
    assert result["weighted_squared_residual_star_neighbors_plane"] == pytest.approx(
        result["weighted_squared_residual_star_plane_only"], abs=1e-10)


def test_joint_scene_refuses_masked_central_pixel_and_stamp_edge():
    image, kernel = _synthetic_field()
    valid = np.ones_like(image,dtype=bool)
    # Rounded central pixel must not be reconstructed or inpainted.
    valid[100,100] = False
    common=dict(kernel=kernel, neighbors=[], guessed_fwhm_pix=6.57,
                original_sky_nMgy_pix=0.)
    result=fit_fixed_psf_scene(
        image,np.ones_like(image),valid,
        target_xy=(100.25,99.65),**common)
    assert result["status"]=="INVALID_CENTRAL_NATIVE_PIXEL"
    result=fit_fixed_psf_scene(
        image,np.ones_like(image),np.ones_like(valid),
        target_xy=(60.,60.),**common)
    assert result["status"]=="FULL_EXPANDED_NATIVE_STAMP_OUTSIDE_IMAGE"


def test_joint_scene_refuses_more_than_bounded_source_budget_without_partial_fit():
    image,kernel=_synthetic_field()
    neighbors=[_source_neighbor(x=115.+i*.01) for i in range(
        MAX_NUISANCE_SOURCES+1)]
    result=fit_fixed_psf_scene(
        image,np.ones_like(image),np.ones_like(image,dtype=bool),
        target_xy=(100.25,99.65),kernel=kernel,
        neighbors=neighbors, guessed_fwhm_pix=6.57,
        original_sky_nMgy_pix=.004)
    assert result["status"]=="TOO_MANY_PSF_NEIGHBORS_FOR_BOUNDED_FIT"
    assert result["science_adopted"] is False


def test_joint_scene_catalog_audit_keeps_original_holdout_and_non_psf_unmodeled():
    image,kernel=_synthetic_field()
    frame=SimpleNamespace(
        image_nMgy=image, invvar_nMgy2=np.ones_like(image),
        valid=np.ones_like(image,dtype=bool), mask_reasons={},
        band="g",dataset_id="dr11/south/abc/g",shape=image.shape)
    cat=SimpleNamespace(
        source_sha256="a"*64,
        dataset_id=frame.dataset_id,match_radius_pix=2.5,
        rows=(
            {"ref_cat":"G3","ref_id":10,"objid":14,"type":"PSF","flux":120.},
            {"ref_cat":"","ref_id":0,"objid":33,"type":"PSF","flux":15.},
            {"ref_cat":"","ref_id":0,"objid":34,"type":"DEV","flux":.1},
        ),
        pixel_positions=np.array(
            [[100.25,99.65],[118.6,100.4],[140.,110.]])
    )
    star=SimpleNamespace(
        star_id="star-0032",x_pix=100.25,y_pix=99.65,
        gaia_source_id=10,cutout_sky_nMgy=.001,role="holdout")
    result=audit_fixed_psf_joint_scenes(
        frame,kernel,[star],cat,
        protected_mask=np.zeros_like(image,dtype=bool),
        guessed_fwhm_pix=6.57)
    assert result["status"]=="FIXED_PSF_ORIGINAL_NATIVE_SCENE_AUDIT_NO_SCIENCE_ADOPTION"
    assert result["n_completed_fixed_holdouts"]==1
    row=result["per_star"][0]
    assert row["source_id"]=="star-0032"
    assert row["role"]=="holdout"
    assert row["n_unmodeled_extended_or_nonpositive_catalog_sources"]==1
    assert row["catalog_scene_fully_modelled"] is False
    assert row["known_unmodeled_catalog_sources_prevent_science_certification"] is True
    assert result["n_completed_scenes_with_unmodeled_catalog_sources"]==1
    assert row["catalog_point_neighbor_count_full_support"]==1
    assert row["neighbors"][0]["tractor_objid"]==33
    assert result["nuisance_parameter_fit_to_holdouts_not_blind_validation"]
    assert result["science_adopted"] is False


def test_joint_scene_gaia_mismatch_cannot_relabel_another_catalog_source():
    image,kernel=_synthetic_field()
    frame=SimpleNamespace(
        image_nMgy=image, invvar_nMgy2=np.ones_like(image),
        valid=np.ones_like(image,dtype=bool), mask_reasons={},
        band="g",dataset_id="dr11/south/abc/g",shape=image.shape)
    cat=SimpleNamespace(source_sha256="a"*64,
        dataset_id=frame.dataset_id,match_radius_pix=2.5,
        rows=({"ref_cat":"G3","ref_id":99,"objid":1,"type":"PSF",
               "flux":100.},),
        pixel_positions=np.array([[100.25,99.65]]))
    star=SimpleNamespace(
        star_id="star-0032",x_pix=100.25,y_pix=99.65,
        gaia_source_id=10,cutout_sky_nMgy=.001,role="holdout")
    report=audit_fixed_psf_joint_scenes(
        frame,kernel,[star],cat,
        protected_mask=np.zeros_like(image,dtype=bool),
        guessed_fwhm_pix=6.57)
    assert report["n_completed_joint_scene_fits"]==0
    assert report["per_star"][0]["status"]=="CENTRAL_GAIA_ID_MISMATCH"


def test_joint_scene_does_not_pretend_fixed_kernel_has_extended_wings():
    image,kernel=_synthetic_field()
    # The actual source has a known additional broad halo; the fitted
    # kernel lacks it. The diagnostic must never report wide-wing recovery.
    yy,xx=np.indices(image.shape,dtype=float)
    halo=np.exp(-.5*(
        (xx-100.25)**2+(yy-99.65)**2)/35.**2)
    original=image+halo*.0003
    result=fit_fixed_psf_scene(
        original,np.ones_like(original),np.ones_like(original,dtype=bool),
        kernel=kernel,target_xy=(100.25,99.65),
        neighbors=[_source_neighbor()],guessed_fwhm_pix=6.57,
        original_sky_nMgy_pix=.004)
    assert result["status"]=="FIXED_PSF_JOINT_SCENE_DIAGNOSTIC"
    assert result["kernel_support_beyond_observed_ePSF_unmodeled"]
    assert result["science_adopted"] is False


def test_joint_scene_gui_summary_has_no_false_science_adoption():
    pytest.importorskip("PySide6")
    from lsbpipe.gui.psf0_panel import PSF0Panel
    metrics={
        "fixed_psf_joint_native_scene_qa": {
            "status": "FIXED_PSF_ORIGINAL_NATIVE_SCENE_AUDIT_NO_SCIENCE_ADOPTION",
            "n_completed_joint_scene_fits": 1,
            "n_preselected_stars": 2,
            "n_completed_fixed_holdouts": 1,
            "n_holdouts_expected": 1,
            "per_star": [
                {"source_id":"star-0032", "role":"holdout",
                 "status":"FIXED_PSF_JOINT_SCENE_DIAGNOSTIC",
                 "n_catalog_point_neighbors":1,
                 "neighbor_to_central_4fwhm_model_flux_ratio": .025,
                 "joint_minus_annulus_sky_nMgy_per_native_pixel":.0004,
                 "n_unmodeled_extended_or_nonpositive_catalog_sources":2},
                {"source_id":"star-0011", "role":"train",
                 "status":"FULL_EXPANDED_NATIVE_STAMP_OUTSIDE_IMAGE"},
            ]
        }
    }
    summary=PSF0Panel._calibration_diagnostic_summary(metrics)
    assert "1/2 stars" in summary
    assert "1/1 reserved Gaia holdouts" in summary
    assert "star-0032" in summary
    assert "2.5%" in summary
    assert "0.0004 nMgy/pix" in summary
    assert "2 unmodeled catalog sources" in summary
    assert "No science adoption" in summary


def test_joint_scene_rejects_fractional_footprint_without_silent_truncation():
    image,kernel=_synthetic_field()
    with pytest.raises(ValueError,match="footprint"):
        fit_fixed_psf_scene(
            image,np.ones_like(image),np.ones_like(image,dtype=bool),
            kernel=kernel,target_xy=(100.25,99.65),
            neighbors=[],guessed_fwhm_pix=6.57,
            original_sky_nMgy_pix=.004,scene_size_pix=160.5)
