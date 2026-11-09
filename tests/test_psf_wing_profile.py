"""P16B wide-radius native source data are not a calibrated wide PSF."""
from copy import deepcopy
from types import SimpleNamespace

import numpy as np
import pytest

from lsbpipe.photometry.psf_wing_profile import (
    audit_native_broad_radius_star_profiles,
    audit_disjoint_bright_gaia_wing_probes,
    PROSPECTIVE_SCHEMA, SCHEMA, P16D_SCHEMA,
    audit_fixed_bright_gaia_neighbor_sky_stress,
    audit_extended_bright_gaia_mask_experiment,
    _p17_adaptive_neighbor_radii,
    P17_SCHEMA, P17_COHORT_LIMIT, P17_MASK_CHOICES,
    P18_SCHEMA, audit_catalog_shape_neighbor_and_sky_sensitivity,
    _p18_neighbor_footprint, _p18_shape_mask,
    _p18_sky_sector_deletion,
)
from lsbpipe.photometry.psf_core_wings_validation import (
    finalize_single_brick_empirical_psf, SINGLE_BRICK_CLOSURE_SCHEMA,
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



def test_p16b_existing_galaxy_target_mask_blocks_outer_sky_pixels():
    frame=_frame(with_neighbor=False)
    yy,xx=np.indices(frame.shape)
    r=np.hypot(xx-200.,yy-200.)
    protected=(r>22.*FWHM)&(r<26.*FWHM)&(xx>200)
    matcher=SimpleNamespace(
        pixel_positions=np.array([[200.,200.]]))
    result=audit_native_broad_radius_star_profiles(
        frame,[_star()],matcher,_geometry(),protected_mask=protected)
    assert result["status"]=="ORIGINAL_NATIVE_OUTER_PROFILE_FEASIBILITY_ONLY"
    assert result["original_target_protection_mask_applied"]
    star=result["per_star"][0]
    assert star["original_target_protection_mask_applied"]
    inner=star["conditional_sky_inner_ring"]
    outer=star["conditional_sky_outer_ring"]
    assert inner["unmasked_coverage_fraction"]<.6
    assert outer["unmasked_coverage_fraction"]==pytest.approx(1.)
    assert star["stellar_wing_detection_not_claimed"]


def test_p16b_mismatched_target_protection_mask_refuses_to_measure():
    matcher=SimpleNamespace(
        pixel_positions=np.array([[200.,200.]]))
    result=audit_native_broad_radius_star_profiles(
        _frame(),[_star()],matcher,_geometry(),
        protected_mask=np.zeros((3,3),dtype=bool))
    assert result["status"]=="INVALID_TARGET_PROTECTION_MASK_SHAPE"
    assert not result["science_adopted"]



def _catalog_row(i, flux, **changes):
    row = {
        "ref_cat":"G3", "ref_id":i, "type":"PSF",
        "fitbits":0, "maskbits":0, "gaia_duplicated_source":False,
        "flux":float(flux), "flux_ivar":100., "objid":1000+i,
    }
    row.update(changes)
    return row


def _prospective_catalog():
    return SimpleNamespace(
        source_sha256="catalog-deterministic-test-hash",
        pixel_positions=np.array([
            [200.,200.], [225.,200.], [200.,225.],
            [225.,225.], [230.,230.],
        ]),
        rows=(
            _catalog_row(11,300.),
            _catalog_row(12,120.),
            _catalog_row(13,130.),
            _catalog_row(14,1000.,maskbits=4),
            _catalog_row(15,200.,gaia_duplicated_source=True),
        ),
    )


def test_p16c_preregistered_bright_catalog_cohort_excludes_train_holdout():
    frame=_frame()
    matcher=_prospective_catalog()
    selected=_star(role="train")
    selected.gaia_source_id=11
    original=frame.image_nMgy.copy()
    result=audit_disjoint_bright_gaia_wing_probes(
        frame,[selected],matcher,_geometry())
    assert result["schema"]==PROSPECTIVE_SCHEMA
    assert result["status"]==(
        "DISJOINT_BRIGHT_GAIA_NATIVE_PROFILE_FEASIBILITY_ONLY")
    assert result["n_native_catalog_entries_examined"]==5
    assert result["n_catalog_probes_fixed_before_original_pixel_reads"]==2
    assert result["n_profiles_with_both_conditional_outer_skies"]==2
    assert result["predeclared_catalog_rejection_counts"]["EXISTING_TRAIN_OR_HOLDOUT_STAR"]==1
    assert result["predeclared_catalog_rejection_counts"]["UNSAFE_CATALOG_QUALITY_FLAGS"]==2
    assert [(r["gaia_source_id"],r["catalog_flux_rank_in_fixed_new_cohort"])
            for r in result["per_star"]]==[(13,1),(12,2)]
    assert all(row["not_frozen_psf_training_or_holdout"]
               for row in result["per_star"])
    assert all(row["native_peak_from_postselection_3x3_original_pixels"]
               for row in result["per_star"])
    assert all(row["science_adopted"] is False
               for row in result["per_star"])
    assert result["no_profile_based_star_replacement"]
    assert result["science_adopted"] is False
    np.testing.assert_array_equal(frame.image_nMgy,original)


def test_p16c_all_qualified_stars_are_fixed_before_pixel_read_no_backfill(
        monkeypatch):
    import lsbpipe.photometry.psf_wing_profile as mod
    monkeypatch.setattr(mod,"PROSPECTIVE_LIMIT",1)
    frame=_frame()
    matcher=_prospective_catalog()
    selected=_star(role="holdout")
    selected.gaia_source_id=11
    # The highest-ranked NEW source is Gaia 13 at (200,225).
    # Removing its observed core pixels must NEVER select Gaia 12
    # instead, since the latter would be a posthoc sample replacement.
    frame.valid[224:227,199:202]=False
    result=mod.audit_disjoint_bright_gaia_wing_probes(
        frame,[selected],matcher,_geometry())
    assert result["n_catalog_probes_fixed_before_original_pixel_reads"]==1
    assert result["per_star"][0]["gaia_source_id"]==13
    assert result["per_star"][0]["status"]==(
        "NO_UNMASKED_ORIGINAL_NATIVE_CORE_PEAK")
    assert result["n_probe_profiles_not_measurable"]==1
    assert result["no_profile_based_star_replacement"]


def test_p16c_new_catalog_cohort_not_fabricated_without_native_rows():
    matcher=SimpleNamespace(pixel_positions=np.array([[200.,200.]]))
    result=audit_disjoint_bright_gaia_wing_probes(
        _frame(),[_star()],matcher,_geometry())
    assert result["status"]=="NO_FULL_NATIVE_GAIA_TRACTOR_CATALOG_METADATA"
    assert result["per_star"]==[]


def test_p16c_prospective_star_cohort_changes_on_catalog_flux_not_image():
    matcher=_prospective_catalog()
    a=audit_disjoint_bright_gaia_wing_probes(
        _frame(with_neighbor=True),[_star()],matcher,_geometry())
    b=audit_disjoint_bright_gaia_wing_probes(
        _frame(with_neighbor=False),[_star()],matcher,_geometry())
    assert [s["gaia_source_id"] for s in a["per_star"]]==(
           [s["gaia_source_id"] for s in b["per_star"]])
    assert [s["catalog_flux_rank_in_fixed_new_cohort"] for s in a["per_star"]]==(
           [s["catalog_flux_rank_in_fixed_new_cohort"] for s in b["per_star"]])
    assert a["catalog_rank_uses_native_tractor_band_flux_descending"]
    assert not a["catalog_selection_uses_image_or_wing_residuals"]


def test_p16c_invalid_protected_target_shape_refuses_prospective_sample():
    result=audit_disjoint_bright_gaia_wing_probes(
        _frame(),[_star()],_prospective_catalog(),_geometry(),
        protected_mask=np.ones((1,2),dtype=bool))
    assert result["status"]=="INVALID_TARGET_PROTECTION_MASK_SHAPE"
    assert not result["science_adopted"]


def test_p16c_catalog_neighbor_overlap_is_excluded_before_profiling():
    catalog=SimpleNamespace(
        source_sha256="catalog",
        pixel_positions=np.array([[200.,200.],[210.,200.]]),
        rows=(_catalog_row(11,300.),_catalog_row(12,120.)))
    result=audit_disjoint_bright_gaia_wing_probes(
        _frame(),[],catalog,_geometry())
    assert result["n_catalog_probes_fixed_before_original_pixel_reads"]==0
    assert result["predeclared_catalog_rejection_counts"]["NEIGHBOR_CATALOG_CORE_OVERLAP"]==2
    assert not result["science_adopted"]



def test_p16c_gui_prospective_cohort_is_descriptive_only():
    pytest.importorskip("PySide6")
    from lsbpipe.gui.psf0_panel import PSF0Panel
    sample=_star(role="train")
    sample.gaia_source_id=11
    qa=audit_disjoint_bright_gaia_wing_probes(
        _frame(),[sample],_prospective_catalog(),_geometry())
    summary=PSF0Panel._calibration_diagnostic_summary({
        "disjoint_bright_gaia_outer_profile_qa":qa})
    assert "P16C disjoint Gaia/Tractor bright-star candidates" in summary
    assert "2/2 original native profiles" in summary
    assert "Distinct from kernel train/holdout identities" in summary
    assert "NOT independently calibrated PSF wings or background" in summary
    assert "wing-cat-0001" in summary



def test_p16c_duplicate_gaia_ref_id_is_not_an_independent_star():
    catalog=SimpleNamespace(
        source_sha256="catalog",
        pixel_positions=np.array([
            [200.,200.],[230.,200.],[200.,230.]]),
        rows=(
            _catalog_row(11,300.),
            _catalog_row(11,400.),
            _catalog_row(12,200.),
        ))
    result=audit_disjoint_bright_gaia_wing_probes(
        _frame(),[],catalog,_geometry())
    assert [v["gaia_source_id"] for v in result["per_star"]]==[12]
    assert result["predeclared_catalog_rejection_counts"][
        "DUPLICATE_GAIA_ID_IN_NATIVE_TRACTOR_CATALOG"]==2
    assert result["n_catalog_probes_fixed_before_original_pixel_reads"]==1
    assert result["science_adopted"] is False



def test_p16d_fixed_star_identities_baseline_and_read_only_stress():
    frame = _frame()
    matcher = _prospective_catalog()
    selected = _star(role="train")
    selected.gaia_source_id = 11
    original = frame.image_nMgy.copy()
    invvar = frame.invvar_nMgy2.copy()
    p16c = audit_disjoint_bright_gaia_wing_probes(
        frame, [selected], matcher, _geometry())
    result = audit_fixed_bright_gaia_neighbor_sky_stress(
        frame, matcher, p16c, _geometry())
    assert result["schema"] == P16D_SCHEMA
    assert result["status"] == "FIXED_COHORT_SKY_NEIGHBOR_STRESS_QA_ONLY"
    assert result["n_fixed_p16c_star_identities"] == 2
    assert result["n_baseline_p16c_reproduced"] == 2
    assert result["n_stress_profiles_with_both_sky_and_plane"] == 2
    assert [v["gaia_source_id"] for v in result["per_star"]] == [
        v["gaia_source_id"] for v in p16c["per_star"]]
    assert [v["catalog_flux_rank_in_fixed_new_cohort"]
            for v in result["per_star"]] == [1, 2]
    assert all(v["status"] == "DESCRIPTIVE_SKY_NEIGHBOR_STRESS_ONLY"
               for v in result["per_star"])
    assert all(v["reproduces_original_p16c_three_fwhm_sky_and_profile"]
               for v in result["per_star"])
    assert all(v["bright_stars_independent_of_existing_psf_training_not_identified"]
               is False for v in p16c["per_star"])
    for star in result["per_star"]:
        baseline, stress = star["mask_variants"]
        assert [v["mask_radius_in_measured_fwhm"]
                for v in star["mask_variants"]] == [3., 6.]
        assert len(baseline["radial_profile"]) == 6
        assert len(stress["radial_profile"]) == 6
        assert baseline["conditional_outer_sky_plane_coefficients_nMgy_per_pix"]
        assert stress["conditional_outer_sky_plane_coefficients_nMgy_per_pix"]
        for smaller, larger in zip(baseline["radial_profile"],
                                   stress["radial_profile"]):
            assert larger["n_unmasked_native_pixels"] <= (
                smaller["n_unmasked_native_pixels"])
            assert len(larger["sky_b_subtracted_octant_medians_over_peak"]) == 8
            assert "six_minus_three_fwhm_mask_residual_over_peak" in larger
    assert result["fixed_p16c_cohort_never_reranked_or_replaced"]
    assert not result["science_adopted"]
    np.testing.assert_array_equal(original, frame.image_nMgy)
    np.testing.assert_array_equal(invvar, frame.invvar_nMgy2)


def test_p16d_requires_preexisting_matching_catalog_and_fixed_p16c_qa():
    frame = _frame()
    matcher = _prospective_catalog()
    report = audit_fixed_bright_gaia_neighbor_sky_stress(
        frame, matcher, {}, _geometry())
    assert report["status"] == "NO_FIXED_P16C_COHORT"
    p16c = audit_disjoint_bright_gaia_wing_probes(
        frame, [], matcher, _geometry())
    wrong = SimpleNamespace(
        pixel_positions=matcher.pixel_positions, source_sha256="other-catalog")
    bad = audit_fixed_bright_gaia_neighbor_sky_stress(
        frame, wrong, p16c, _geometry())
    assert bad["status"] == "P16C_CATALOG_LINEAGE_MISMATCH"
    assert bad["per_star"] == []
    assert not bad["science_adopted"]


def test_p16d_selected_but_failed_source_is_not_replaced():
    frame = _frame()
    matcher = _prospective_catalog()
    selected = _star(role="train")
    selected.gaia_source_id = 11
    frame.valid[224:227, 199:202] = False
    p16c = audit_disjoint_bright_gaia_wing_probes(
        frame, [selected], matcher, _geometry())
    result = audit_fixed_bright_gaia_neighbor_sky_stress(
        frame, matcher, p16c, _geometry())
    assert len(result["per_star"]) == 2
    assert result["per_star"][0]["gaia_source_id"] == 13
    assert result["per_star"][0]["status"] == "P16C_PROFILE_NOT_AVAILABLE"
    assert result["per_star"][1]["gaia_source_id"] == 12
    assert result["n_baseline_p16c_reproduced"] == 1
    assert not result["science_adopted"]


def test_p16d_repeated_input_is_deterministic():
    frame = _frame()
    matcher = _prospective_catalog()
    p16c = audit_disjoint_bright_gaia_wing_probes(
        frame, [], matcher, _geometry())
    a = audit_fixed_bright_gaia_neighbor_sky_stress(
        frame, matcher, p16c, _geometry())
    b = audit_fixed_bright_gaia_neighbor_sky_stress(
        frame, matcher, p16c, _geometry())
    assert a == b
    assert a["n_baseline_p16c_reproduced"] == 3


def test_p16d_synthetic_linear_sky_gradient_is_recorded_not_adopted():
    frame = _frame(with_neighbor=False)
    matcher = _prospective_catalog()
    p16c = audit_disjoint_bright_gaia_wing_probes(
        frame, [], matcher, _geometry())
    result = audit_fixed_bright_gaia_neighbor_sky_stress(
        frame, matcher, p16c, _geometry())
    assert result["n_baseline_p16c_reproduced"] == 3
    assert result["n_stress_profiles_with_both_sky_and_plane"] == 3
    assert any(abs(v["mask_variants"][0][
        "conditional_outer_sky_plane_coefficients_nMgy_per_pix"][1]) > 0
        for v in result["per_star"])
    assert result["plane_is_conditional_and_can_absorb_stellar_wings"]
    assert result["no_covariance_calibration_or_statistical_detection_test"]
    assert not result["science_adopted"]


def test_p16d_gui_summary_does_not_claim_wing_detection():
    pytest.importorskip("PySide6")
    from lsbpipe.gui.psf0_panel import PSF0Panel
    frame = _frame()
    matcher = _prospective_catalog()
    p16c = audit_disjoint_bright_gaia_wing_probes(
        frame, [], matcher, _geometry())
    stress = audit_fixed_bright_gaia_neighbor_sky_stress(
        frame, matcher, p16c, _geometry())
    summary = PSF0Panel._calibration_diagnostic_summary({
        "fixed_bright_gaia_neighbor_sky_stress_qa": stress})
    assert "P16D fixed P16C bright Gaia cohort" in summary
    assert "NOT calibrated wide-angle PSF" in summary



def _p17_synthetic_field():
    """48 clean Gaia positions, 1 catalog galaxy neighbor; original pixels."""
    shape = (840, 880)
    yy, xx = np.indices(shape, dtype=float)
    data = .001 + .000006 * xx + .000004 * yy
    positions = []
    rows = []
    for index in range(48):
        x = 180 + 70 * (index % 8)
        y = 180 + 70 * (index // 8)
        data[y, x] += 5. + (48 - index) / 48.
        positions.append([float(x), float(y)])
        rows.append(_catalog_row(10000 + index, 500 - 3 * index))
    # This non-Gaia neighbor is >3 FWHM away: does not change 48-star
    # selection, but must receive the conservative adaptive mask.
    positions.append([215., 180.])
    rows.append(_catalog_row(999, 100., ref_cat="", type="EXP"))
    frame = SimpleNamespace(
        image_nMgy=data,
        invvar_nMgy2=np.ones(shape, dtype=float) * 100.,
        valid=np.ones(shape, dtype=bool),
        mask_reasons={},
        band="g",
        shape=shape,
    )
    matcher = SimpleNamespace(
        source_sha256="fixed-p17-test-catalog-hash",
        pixel_positions=np.asarray(positions, dtype=float),
        rows=tuple(rows),
    )
    return frame, matcher


def test_p17_expands_frozen_legacy_cohort_and_preserves_holdout():
    frame, matcher = _p17_synthetic_field()
    original = frame.image_nMgy.copy()
    invvar = frame.invvar_nMgy2.copy()
    p16c = audit_disjoint_bright_gaia_wing_probes(
        frame, [], matcher, _geometry())
    assert len(p16c["per_star"]) == 16
    p17 = audit_extended_bright_gaia_mask_experiment(
        frame, [], matcher, _geometry(), p16c)
    assert p17["schema"] == P17_SCHEMA
    assert p17["status"] == "FIXED_EXPANDED_COHORT_MASK_FEASIBILITY_ONLY"
    assert p17["n_catalog_only_eligible_before_limit"] == 48
    assert p17["n_catalog_selected"] == P17_COHORT_LIMIT
    assert p17["n_legacy_p16c_stars"] == 16
    assert p17["n_new_exploratory"] == 24
    assert p17["n_new_identity_holdout"] == 8
    assert p17["n_legacy_baseline_reproduced"] == 16
    assert p17["n_multimask_stars_measured"] == 48
    assert [v["gaia_source_id"] for v in p17["per_star"][:16]] == [
        v["gaia_source_id"] for v in p16c["per_star"]]
    assert [v["catalog_rank"] for v in p17["per_star"]
            if v["cohort_role"] == "NEW_CATALOG_HOLDOUT"] == [
                20, 24, 28, 32, 36, 40, 44, 48]
    assert [v["mask_variants"][0]["strategy"] for v in p17["per_star"]] == [
        "fixed_3_fwhm"] * 48
    assert all([v["strategy"] for v in item["mask_variants"]] ==
               list(P17_MASK_CHOICES) for item in p17["per_star"])
    assert all(item["status"] == "MULTIMASK_DESCRIPTIVE_ONLY"
               for item in p17["per_star"])
    for item in p17["per_star"]:
        a, b, adaptive, c = item["mask_variants"]
        for small, medium, mix, large in zip(
                a["radial_profile"], b["radial_profile"],
                adaptive["radial_profile"], c["radial_profile"]):
            assert (large["n_unmasked_native_pixels"] <=
                    mix["n_unmasked_native_pixels"] <=
                    small["n_unmasked_native_pixels"])
            assert (large["n_unmasked_native_pixels"] <=
                    medium["n_unmasked_native_pixels"] <=
                    small["n_unmasked_native_pixels"])
            assert len(mix["sky_b_subtracted_octant_medians_over_peak"]) == 8
    assert not p17["science_adopted"]
    assert p17["new_control_stars_not_independent_of_same_coadd"]
    np.testing.assert_array_equal(original, frame.image_nMgy)
    np.testing.assert_array_equal(invvar, frame.invvar_nMgy2)


def test_p17_catalog_flux_adaptive_mask_tiers_fixed_without_pixels():
    rows = (
        {"type": "PSF", "flux": .5},          # .005 -> 3
        {"type": "PSF", "flux": 5.},          # .05 -> 4.5
        {"type": "PSF", "flux": 20.},         # .20 -> 6
        {"type": "EXP", "flux": .1},          # morphology -> 6
        {"type": "PSF", "flux": None},        # absent -> 6
        {"type": "PSF", "flux": float("nan")},# nan -> 6
        {"type": "PSF", "flux": 1.},          # exactly .01 -> 4.5
        {"type": "PSF", "flux": 10.},         # exactly .1 -> 6
    )
    result = _p17_adaptive_neighbor_radii(rows, 100., 2.)
    np.testing.assert_array_equal(result,
                                  np.array([6., 9., 12., 12., 12.,
                                            12., 9., 12.]))


def test_p17_invalid_catalog_or_incomplete_p16c_fails_closed():
    frame = _frame()
    matcher = _prospective_catalog()
    result = audit_extended_bright_gaia_mask_experiment(
        frame, [], matcher, _geometry(), {})
    assert result["status"] == "NO_COMPLETE_FROZEN_16_STAR_P16C_COHORT"
    p16c = audit_disjoint_bright_gaia_wing_probes(
        frame, [], matcher, _geometry())
    result = audit_extended_bright_gaia_mask_experiment(
        frame, [], matcher, _geometry(), p16c)
    assert result["status"] == "NO_COMPLETE_FROZEN_16_STAR_P16C_COHORT"
    assert result["per_star"] == []
    assert not result["science_adopted"]


def test_p17_gui_summary_marks_masks_as_uncalibrated():
    pytest.importorskip("PySide6")
    from lsbpipe.gui.psf0_panel import PSF0Panel
    summary = PSF0Panel._calibration_diagnostic_summary({
        "extended_bright_gaia_mask_experiment_qa": {
            "schema": P17_SCHEMA,
            "status": "FIXED_EXPANDED_COHORT_MASK_FEASIBILITY_ONLY",
            "n_catalog_selected": 48,
            "n_legacy_p16c_stars": 16,
            "n_new_exploratory": 24,
            "n_new_identity_holdout": 8,
            "n_legacy_baseline_reproduced": 16,
            "n_multimask_stars_measured": 48,
        },
    })
    assert "P17 fixed catalog-only Leo P cohort" in summary
    assert "24 new exploratory" in summary
    assert "8 new identity holdout" in summary
    assert "NOT calibrated PSF wings" in summary



def test_p18_tractor_shape_ellipticity_and_missing_shape_fallback():
    fwhm, scale = 2., .5
    point = _p18_neighbor_footprint({
        "type": "PSF", "flux": 1.
    }, 100., fwhm, scale)
    assert point == (3., 3., 0., "CATALOG_PSF")
    extended = _p18_neighbor_footprint({
        "type": "EXP", "flux": 1., "shape_r_arcsec": 2.,
        "shape_e1": .6, "shape_e2": 0.
    }, 100., fwhm, scale)
    assert extended[0] > extended[1] > point[0]
    assert extended[2] == 0.
    assert extended[3] == "EXTENDED_SHAPE_ELLIPSE"
    missing = _p18_neighbor_footprint({
        "type": "DEV", "flux": .01, "shape_r_arcsec": None,
        "shape_e1": 0., "shape_e2": 0.
    }, 100., fwhm, scale)
    assert missing[-1] == "MISSING_SHAPE_FALLBACK"
    assert missing[:2] == (6., 6.)
    circular = _p18_neighbor_footprint({
        "type": "REX", "flux": 1., "shape_r_arcsec": 2.,
        "shape_e1": None, "shape_e2": None
    }, 100., fwhm, scale)
    assert circular[0] == circular[1]
    assert circular[-1] == "EXTENDED_SHAPE_ELLIPSE"
    bad = _p18_neighbor_footprint({
        "type": "EXP", "flux": 1., "shape_r_arcsec": 2.,
        "shape_e1": 1.2, "shape_e2": 0.
    }, 100., fwhm, scale)
    assert bad[-1] == "INVALID_ELLIPTICITY_FALLBACK"


def test_p18_rotated_ellipse_masks_major_axis_not_all_octants():
    catalog = [
        {"type": "PSF", "flux": 100.},
        {"type": "EXP", "flux": 1., "shape_r_arcsec": 2.,
         "shape_e1": .6, "shape_e2": 0.},
    ]
    keep, flags, n = _p18_shape_mask(
        (101,101), 0, 0,
        np.array([[50.,50.], [70.,50.]]), catalog,
        0, 100., 2., np.array([[.5,0.],[0.,.5]]))
    assert n == 1
    assert flags == {"EXTENDED_SHAPE_ELLIPSE": 1}
    assert not keep[50,80]
    assert keep[60,70]
    # 90-degree rotation of ellipticity must rotate the footprint.
    catalog[1]["shape_e1"] = -.6
    rotated, _, _ = _p18_shape_mask(
        (101,101), 0, 0,
        np.array([[50.,50.], [70.,50.]]), catalog,
        0, 100., 2., np.array([[.5,0.],[0.,.5]]))
    assert rotated[50,80]
    assert not rotated[60,70]


def test_p18_sky_sector_deletion_is_sensitivity_not_errorbar():
    yy, xx = np.indices((241,241), dtype=float)
    radius = np.hypot(xx-120,yy-120)
    octants = (np.floor(
        (np.arctan2(yy-120,xx-120)%(2*np.pi))*8/(2*np.pi)
    ).astype(int))
    pixels = .01 + .00001*xx + .00002*yy
    full, replicates, sensitivity = _p18_sky_sector_deletion(
        pixels, np.ones_like(pixels,dtype=bool), radius, octants,
        3., 10.)
    assert np.isfinite(full)
    assert len(replicates) == 8
    assert sensitivity > 0
    assert len(set(replicates)) > 1


def test_p18_fixed_cohort_shape_footprints_and_sky_are_read_only():
    from astropy.wcs import WCS
    frame, matcher = _p17_synthetic_field()
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [440.,420.]
    wcs.wcs.crval = [155.,18.]
    wcs.wcs.cdelt = [-.262/3600., .262/3600.]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    frame.wcs = wcs
    # Original rows 48 Gaia stars. One non-Gaia extended
    # neighbor with genuine SHAPE_R/E1/E2, not an invented point PSF.
    matcher.rows[48]["shape_r_arcsec"] = 1.2
    matcher.rows[48]["shape_e1"] = .2
    matcher.rows[48]["shape_e2"] = 0.
    original = frame.image_nMgy.copy()
    inverse = frame.invvar_nMgy2.copy()
    valid = frame.valid.copy()
    p16c = audit_disjoint_bright_gaia_wing_probes(
        frame, [], matcher, _geometry())
    p17 = audit_extended_bright_gaia_mask_experiment(
        frame, [], matcher, _geometry(), p16c)
    p18 = audit_catalog_shape_neighbor_and_sky_sensitivity(
        frame, matcher, _geometry(), p17)
    assert p18["schema"] == P18_SCHEMA
    assert p18["status"] == (
        "FIXED_P17_COHORT_MORPHOLOGY_SKY_DIAGNOSTIC_ONLY")
    assert p18["n_fixed_star_identities"] == 48
    assert p18["n_new_identity_holdout"] == 8
    assert p18["n_profiles_with_both_sky_and_plane"] == 48
    assert [r["gaia_source_id"] for r in p18["per_star"]] == [
        r["gaia_source_id"] for r in p17["per_star"]]
    assert all(v["status"] ==
               "MORPHOLOGY_MASK_AND_SKY_DESCRIPTIVE_ONLY"
               for v in p18["per_star"])
    assert all(len(v["radial_profile"]) == 6
               for v in p18["per_star"])
    assert all(len(v["sky_b_leave_one_octant_out_medians_nMgy"]) == 8
               for v in p18["per_star"])
    assert any(v["catalog_footprint_class_counts"].get(
        "EXTENDED_SHAPE_ELLIPSE",0) > 0 for v in p18["per_star"])
    assert not p18["science_adopted"]
    np.testing.assert_array_equal(original,frame.image_nMgy)
    np.testing.assert_array_equal(inverse,frame.invvar_nMgy2)
    np.testing.assert_array_equal(valid,frame.valid)


def test_p18_no_wcs_or_catalog_lineage_fails_closed():
    frame, matcher = _p17_synthetic_field()
    p16c = audit_disjoint_bright_gaia_wing_probes(
        frame, [], matcher, _geometry())
    p17 = audit_extended_bright_gaia_mask_experiment(
        frame, [], matcher, _geometry(), p16c)
    altered = SimpleNamespace(
        source_sha256="another-catalog-hash",
        rows=matcher.rows, pixel_positions=matcher.pixel_positions)
    bad = audit_catalog_shape_neighbor_and_sky_sensitivity(
        frame, altered, _geometry(), p17)
    assert bad["status"] == "CATALOG_LINEAGE_MISMATCH"
    assert bad["per_star"] == []
    without_wcs = audit_catalog_shape_neighbor_and_sky_sensitivity(
        frame, matcher, _geometry(), p17)
    assert without_wcs["n_profiles_with_both_sky_and_plane"] == 0
    assert all(v["status"] == "MISSING_OR_UNSAFE_WCS_TANGENT_MATRIX"
               for v in without_wcs["per_star"])
    assert not without_wcs["science_adopted"]


def test_p18_gui_displays_only_descriptive_status():
    pytest.importorskip("PySide6")
    from lsbpipe.gui.psf0_panel import PSF0Panel
    text = PSF0Panel._calibration_diagnostic_summary({
        "catalog_shape_neighbor_sky_qa": {
            "schema": P18_SCHEMA,
            "status": "FIXED_P17_COHORT_MORPHOLOGY_SKY_DIAGNOSTIC_ONLY",
            "n_fixed_star_identities":48,
            "n_profiles_with_both_sky_and_plane":44,
        },
    })
    assert "P18 catalog morphology / sky sensitivity" in text
    assert "48 fixed P17 identities" in text
    assert "NO calibrated sky uncertainty" in text



def _single_brick_closure_synthetic_evidence():
    from astropy.wcs import WCS
    frame, matcher = _p17_synthetic_field()
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [440.,420.]
    wcs.wcs.crval = [155.,18.]
    wcs.wcs.cdelt = [-.262/3600., .262/3600.]
    wcs.wcs.ctype = ["RA---TAN","DEC--TAN"]
    frame.wcs = wcs
    p16c = audit_disjoint_bright_gaia_wing_probes(
        frame, [], matcher, _geometry())
    p17 = audit_extended_bright_gaia_mask_experiment(
        frame, [], matcher, _geometry(), p16c)
    p18 = audit_catalog_shape_neighbor_and_sky_sensitivity(
        frame, matcher, _geometry(), p17)
    core = {
        "status": "SCIENCE_GATE_HOLD",
        "n_train_stars": 10,
        "n_fixed_holdout_stars": 3,
        "core_geometry": dict(
            _geometry(), fully_contained_radial_radius_limit_native_pix=40.),
    }
    return frame, p17, p18, core


def test_single_brick_closure_uses_only_frozen_matched_qa_and_keeps_hold():
    frame, p17, p18, core = _single_brick_closure_synthetic_evidence()
    snapshot = frame.image_nMgy.copy()
    args = dict(
        dataset_id="dr11/south/1554p180/g",
        image_shape=frame.shape,
        core_validation=core,
        p17_qa=p17, p18_qa=p18,
        verified_native_sources=True,
    )
    finished = finalize_single_brick_empirical_psf(**args)
    assert finished["schema"] == SINGLE_BRICK_CLOSURE_SCHEMA
    assert finished["status"] == (
        "BRICK_ONLY_EMPIRICAL_CORE_AND_WING_CONSTRAINTS_COMPLETE")
    assert finished["engineering_stage_finished"]
    assert not finished["science_adopted"]
    assert not finished["scientifically_validated_extended_wings"]
    assert not finished["wide_wing_photometric_correction_permitted"]
    assert finished["no_other_brick_required_or_used"]
    assert finished["n_frozen_stars"] == 48
    assert finished["n_original_kernel_training_stars"] == 10
    assert finished["n_original_kernel_holdout_stars"] == 3
    assert finished["outer_diagnostic_radii_fwhm"] == [4., 22.]
    assert finished["finite_core_kernel_full_radial_support_native_pix"] == 40.
    for label, expected in zip((
            "P16C_LEGACY_EXPLORATORY", "NEW_CATALOG_EXPLORATORY",
            "NEW_CATALOG_HOLDOUT"), (16, 24, 8)):
        assert finished["cohorts"][label]["n_stars"] == expected
        assert len(finished["cohorts"][label]["radial_profiles"]) == 6
        assert finished["cohorts"][label][
            "sky_sensitivity"]["eight_deletions_not_independent_error_estimates"]
    assert sum(x["n_fixed_stars"] for x in
               finished["same_brick_catalog_pixel_quadrants"]) == 48
    assert sum(x["n_legacy_bright_stars"] for x in
               finished["same_brick_catalog_pixel_quadrants"]) == 16
    np.testing.assert_array_equal(snapshot, frame.image_nMgy)
    # The diagnostic must follow IDs, never ordering of a JSON list.
    reverse = dict(args, p18_qa=dict(p18, per_star=p18["per_star"][::-1]))
    shuffled = finalize_single_brick_empirical_psf(**reverse)
    assert shuffled["cohorts"] == finished["cohorts"]
    assert shuffled["same_brick_catalog_pixel_quadrants"] == (
        finished["same_brick_catalog_pixel_quadrants"])


def test_single_brick_closure_fails_closed_on_provenance_and_identity():
    frame, p17, p18, core = _single_brick_closure_synthetic_evidence()
    args = dict(
        dataset_id="dr11/south/1554p180/g",
        image_shape=frame.shape,
        core_validation=core,
        p17_qa=p17, p18_qa=p18,
        verified_native_sources=True,
    )
    assert finalize_single_brick_empirical_psf(**dict(
        args, verified_native_sources=False))["status"] == (
            "ORIGINAL_NATIVE_SOURCE_HASHES_UNVERIFIED")
    bad_core = dict(core, core_geometry={"status":"UNMEASURED"})
    assert finalize_single_brick_empirical_psf(**dict(
        args, core_validation=bad_core))["status"] == (
            "NO_MEASURED_FINITE_CORE")
    altered = dict(p18, native_catalog_sha256="different")
    assert finalize_single_brick_empirical_psf(**dict(
        args, p18_qa=altered))["status"] == (
            "NATIVE_TRACTOR_CATALOG_HASH_MISMATCH")
    missing = dict(p18, per_star=p18["per_star"][:-1])
    assert finalize_single_brick_empirical_psf(**dict(
        args, p18_qa=missing))["status"] == "INCOMPLETE_48_STAR_COHORT"
    altered_id = deepcopy(p18)
    altered_id["per_star"][4]["gaia_source_id"] = -8
    assert finalize_single_brick_empirical_psf(**dict(
        args, p18_qa=altered_id))["status"] == (
            "P18_P17_STAR_IDENTITY_MISMATCH")
    assert not finalize_single_brick_empirical_psf(**dict(
        args, p18_qa=altered_id))["engineering_stage_finished"]


def test_single_brick_closure_sky_available_does_not_hide_empty_ring():
    frame, p17, p18, core = _single_brick_closure_synthetic_evidence()
    # Preserve 48 successful sky solutions, but one mask destroys an
    # internal stellar annulus. No fake zero-valued measurement.
    altered = deepcopy(p18)
    ring = altered["per_star"][19]["radial_profile"][0]
    ring["sky_b_subtracted_median_over_peak"] = None
    ring["unmasked_fraction"] = 0.
    ring["n_octants_with_at_least_25_pixels"] = 0
    result = finalize_single_brick_empirical_psf(
        dataset_id="dr11/south/1554p180/g",
        image_shape=frame.shape,
        core_validation=core,
        p17_qa=p17, p18_qa=altered,
        verified_native_sources=True)
    assert result["engineering_stage_finished"]
    assert result["n_p18_radial_annuli_with_no_pixels"] == 1
    row = result["cohorts"]["all_48"]["radial_profiles"][0]
    assert row["n_nonnull_p18_radial_residuals"] == 47
    assert row["n_stars"] == 48
    assert not result["scientifically_validated_extended_wings"]


def test_single_brick_closure_gui_never_implies_science_adoption():
    pytest.importorskip("PySide6")
    from lsbpipe.gui.psf0_panel import PSF0Panel
    summary = PSF0Panel._calibration_diagnostic_summary({
        "single_brick_empirical_psf_closure_qa": {
            "schema": SINGLE_BRICK_CLOSURE_SCHEMA,
            "status": "BRICK_ONLY_EMPIRICAL_CORE_AND_WING_CONSTRAINTS_COMPLETE",
            "engineering_stage_finished": True,
            "n_frozen_stars": 48,
            "finite_core_measured_fwhm_native_pix": 5.13690,
        },
    })
    assert "single-brick characterization CLOSED" in summary
    assert "NO validated outer-wing kernel" in summary
    assert "THIS native brick" in summary
