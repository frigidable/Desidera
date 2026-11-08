"""P16B wide-radius native source data are not a calibrated wide PSF."""
from types import SimpleNamespace

import numpy as np
import pytest

from lsbpipe.photometry.psf_wing_profile import (
    audit_native_broad_radius_star_profiles, SCHEMA,
)

FWHM=5.136902957677135


def _frame(*, with_neighbor=True):
    shape=(401,401)
    yy,xx=np.indices(shape,dtype=float)
    r=np.hypot(xx-200.,yy-200.)
    original=(.01 + .001*(xx-200.)/100
              + 10*np.exp(-.5*(r/2.8)**2)
              + .03/(1+(r/14)**2))
    if with_neighbor:
        original=original+60*np.exp(
            -.5*((xx-225.)**2+(yy-200.)**2)/2.5**2)
    return SimpleNamespace(
        image_nMgy=original,
        invvar_nMgy2=np.ones(shape)*100,
        valid=np.ones(shape,dtype=bool),
        mask_reasons={},
        band="g",
        shape=shape,
    )


def _star(x=200.,y=200.,role="holdout"):
    return SimpleNamespace(star_id="star-0032",x_pix=x,y_pix=y,
                           role=role,peak_nMgy=10.)


def _geometry():
    return {
        "status":"FINITE_NATIVE_CORE_GEOMETRY_MEASURED",
        "measured_core_fwhm_axis_geometric_mean_native_pix":FWHM,
    }


def test_p16b_read_only_native_profile_masks_catalog_neighbor_and_sky():
    frame=_frame()
    original=frame.image_nMgy.copy()
    matcher=SimpleNamespace(pixel_positions=np.array([
        [200.,200.],[225.,200.]]))
    result=audit_native_broad_radius_star_profiles(
        frame,[_star()],matcher,_geometry())
    assert result["schema"]==SCHEMA
    assert result["status"]=="ORIGINAL_NATIVE_OUTER_PROFILE_FEASIBILITY_ONLY"
    assert result["n_selected_stars"]==1
    assert result["n_scenes_with_two_outer_sky_hypotheses"]==1
    assert result["n_fixed_holdouts_in_original_selected_inventory"]==1
    assert result["n_fixed_holdouts_with_outer_profile"]==1
    star=result["per_star"][0]
    assert star["status"]=="MASKED_OUTER_PROFILE_DESCRIPTIVE_SKY_DEGENERATE"
    assert star["n_catalog_neighbors_whose_cores_were_masked"]==1
    assert star["catalog_mask_only_covers_neighbor_cores_not_their_wings"]
    assert star["sky_hypotheses_from_disjoint_original_native_annuli"]
    assert star["conditional_sky_inner_ring"]["radius_inner_in_measured_fwhm"]==22.
    assert star["conditional_sky_inner_ring"]["radius_outer_in_measured_fwhm"]==26.
    assert star["conditional_sky_outer_ring"]["radius_inner_in_measured_fwhm"]==26.
    assert star["conditional_sky_outer_ring"]["radius_outer_in_measured_fwhm"]==30.
    assert star["sky_inner_minus_outer_nMgy_per_native_pixel"]>0
    assert len(star["radial_profile"])==9
    assert all(x["n_unmasked_native_pixels"]>0
               for x in star["radial_profile"])
    assert all(0<x["unmasked_coverage_fraction"]<=1
               for x in star["radial_profile"])
    assert star["radial_profile"][0]["unmasked_coverage_fraction"]<1
    assert star["bright_stars_independent_of_existing_psf_training_not_identified"]
    assert star["stellar_wing_detection_not_claimed"]
    assert result["no_original_pixels_or_psf_kernel_modified"]
    assert not result["science_adopted"]
    np.testing.assert_array_equal(original,frame.image_nMgy)


def test_p16b_neighbor_core_mask_removes_large_catalog_point_contamination():
    mask_catalog=SimpleNamespace(
        pixel_positions=np.array([[200.,200.],[225.,200.]]))
    a=audit_native_broad_radius_star_profiles(
        _frame(with_neighbor=True),[_star()],mask_catalog,_geometry())
    b=audit_native_broad_radius_star_profiles(
        _frame(with_neighbor=False),[_star()],mask_catalog,_geometry())
    assert a["per_star"][0]["status"]==b["per_star"][0]["status"]
    for first,second in zip(
            a["per_star"][0]["radial_profile"],
            b["per_star"][0]["radial_profile"]):
        assert first["median_original_nMgy_per_native_pixel"]==(
            pytest.approx(
                second["median_original_nMgy_per_native_pixel"],
                abs=1e-6))
    assert a["per_star"][0]["catalog_mask_only_covers_neighbor_cores_not_their_wings"]


def test_p16b_missing_catalog_fails_closed_not_a_neighbor_free_measurement():
    result=audit_native_broad_radius_star_profiles(
        _frame(),[_star()],None,_geometry())
    assert result["status"]=="NO_VERIFIED_NATIVE_CATALOG_FOR_NEIGHBOR_MASK"
    assert result["per_star"]==[]
    assert result["science_adopted"] is False


def test_p16b_edge_star_does_not_clip_its_sky_apertures():
    matcher=SimpleNamespace(
        pixel_positions=np.array([[120.,200.]]))
    result=audit_native_broad_radius_star_profiles(
        _frame(),[_star(x=120.)],matcher,_geometry())
    assert result["n_scenes_not_measurable"]==1
    assert result["per_star"][0]["status"]==(
        "FULL_OUTER_NATIVE_CUTOUT_OUTSIDE_IMAGE")
    assert "radial_profile" not in result["per_star"][0]


def test_p16b_ambiguous_catalog_central_counterpart_fails_closed():
    matcher=SimpleNamespace(pixel_positions=np.array([
        [200.,200.],[201.,199.],[250.,250.]]))
    result=audit_native_broad_radius_star_profiles(
        _frame(),[_star()],matcher,_geometry())
    assert result["per_star"][0]["status"]=="CATALOG_CENTRAL_MATCH_NOT_UNIQUE"
    assert result["per_star"][0]["n_possible_catalog_centers"]==2
    assert result["n_scenes_with_two_outer_sky_hypotheses"]==0


def test_p16b_no_measured_fwhm_fails_closed():
    matcher=SimpleNamespace(
        pixel_positions=np.array([[200.,200.]]))
    result=audit_native_broad_radius_star_profiles(
        _frame(),[_star()],matcher,{"status":"UNBRACKETED_CORE_HALF_MAXIMUM"})
    assert result["status"]=="NO_MEASURED_NATIVE_CORE_FWHM"
    assert not result["science_adopted"]


def test_p16b_saturated_or_bad_original_native_pixels_are_excluded():
    frame=_frame(with_neighbor=False)
    yy,xx=np.indices(frame.shape)
    bad=(xx>200)&(xx<230)&(yy>190)&(yy<210)
    frame.mask_reasons["SATUR_G"]=bad
    matcher=SimpleNamespace(pixel_positions=np.array([[200.,200.]]))
    result=audit_native_broad_radius_star_profiles(
        frame,[_star()],matcher,_geometry())
    star=result["per_star"][0]
    assert star["status"]=="MASKED_OUTER_PROFILE_DESCRIPTIVE_SKY_DEGENERATE"
    first=star["radial_profile"][0]
    assert first["unmasked_coverage_fraction"]<1
    assert star["original_image_and_inverse_variance_unchanged"]


def test_p16b_repeated_input_is_deterministic():
    frame=_frame(with_neighbor=True)
    matcher=SimpleNamespace(pixel_positions=np.array([
        [200.,200.],[225.,200.]]))
    a=audit_native_broad_radius_star_profiles(
        frame,[_star()],matcher,_geometry())
    b=audit_native_broad_radius_star_profiles(
        frame,[_star()],matcher,_geometry())
    assert a==b



def test_p16b_gui_reports_original_profile_without_wing_adoption():
    pytest.importorskip("PySide6")
    from lsbpipe.gui.psf0_panel import PSF0Panel
    matcher=SimpleNamespace(pixel_positions=np.array([
        [200.,200.],[225.,200.]]))
    q=audit_native_broad_radius_star_profiles(
        _frame(),[_star()],matcher,_geometry())
    summary=PSF0Panel._calibration_diagnostic_summary({
        "original_native_outer_stellar_profile_qa":q})
    assert "P16B original-pixel broad-radius" in summary
    assert "fixed holdouts 1/1" in summary
    assert "Two disjoint outer-sky hypotheses" in summary
    assert "P16B sky-annulus median difference / stellar peak=" in summary
    assert "NOT a calibrated wide PSF or sky" in summary
