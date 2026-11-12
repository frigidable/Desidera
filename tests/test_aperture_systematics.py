"""Stage 04.3 scientific diagnostics: empirical blanks and exact sensitivity.

Tests deliberately do NOT certify a science total, a sky uncertainty, or
restoration of galaxy light hidden beneath provisional source masks.
"""
import json

import numpy as np
import pytest

pytest.importorskip("astropy")
pytest.importorskip("photutils")

from test_aperture_flux import frame, provenance, geometry, wcs
from lsbpipe.photometry.aperture_flux import (
    ApertureSkyScenario, measure_native_aperture,
)
from lsbpipe.photometry.aperture_growth import (
    GrowthHypothesis, measure_growth_diagnostic,
)
from lsbpipe.photometry.aperture_systematics import (
    diagnose_aperture_systematics, estimate_blank_annulus_covariance,
)


def _curve(f, source, manual, *, q=1., pa=0., radii=(1., 2.5, 4.), sky=None):
    return measure_growth_diagnostic(
        f, anchor_icrs=(155.4, 18.), provenance=provenance(),
        effective_source_mask=source, manual_exclude=manual,
        radii_arcsec=radii,
        hypothesis=GrowthHypothesis("fixed_test", q, pa, "SYNTHETIC"),
        sky=sky)


def test_blank_annulus_covariance_is_deterministic_signed_and_psd():
    rng = np.random.default_rng(20261111)
    img = rng.normal(0., .5, (321, 321))
    f = frame(image=img, shape=img.shape)
    original = f.image_nMgy.copy()
    none = np.zeros(f.shape, bool)
    outer = geometry(a=1.6, q=.65, pa=29.)
    kw = dict(
        frame=f, outer=outer, radii_arcsec=(.5, 1., 1.6),
        effective_source_mask=none, manual_exclude=none)
    result = estimate_blank_annulus_covariance(**kw)
    assert result["status"] == "EMPIRICAL_BLANK_COVARIANCE_DIAGNOSTIC_ONLY"
    assert result["n_accepted"] >= result["n_required"] >= 12
    assert result["n_accepted"] <= 48
    assert result["covariance_calibrated_for_science"] is False
    matrix = np.asarray(result["annulus_covariance_nMgy2"])
    assert matrix.shape == (3, 3)
    np.testing.assert_allclose(matrix, matrix.T, atol=1e-10)
    transfer = np.tril(np.ones((3, 3)))
    cumulative = np.asarray(result["cumulative_covariance_nMgy2"])
    np.testing.assert_allclose(cumulative, transfer @ matrix @ transfer.T)
    np.testing.assert_allclose(
        np.square(result["cumulative_sigma_empirical_nMgy"]),
        np.diag(cumulative), rtol=1e-10, atol=1e-10)
    assert np.linalg.eigvalsh(matrix).min() >= -1e-8
    assert np.all(np.diag(matrix) > 0)
    assert any(abs(v) > 0 for v in result["annulus_mean_blank_flux_nMgy"])
    assert result == estimate_blank_annulus_covariance(**kw)
    assert np.array_equal(f.image_nMgy, original)
    json.dumps(result, allow_nan=False)


def test_bounded_memory_stratified_candidates_large_grid_does_not_allocate_centres():
    """Tens of millions of candidate nodes require only row-block metadata."""
    from itertools import islice
    from lsbpipe.photometry.aperture_systematics import (
        _bounded_quadrant_blank_candidates)
    shape=(40000,40000)
    margin=6
    step=12
    counts,seq=_bounded_quadrant_blank_candidates(
        shape,margin,step)
    assert len(counts)==4
    assert sum(counts)>10_000_000
    first=list(islice(seq,80))
    assert len(first)==80
    assert {v[0] for v in first}=={0,1,2,3}
    assert len({(x,y) for _,x,y in first})==80
    assert all(margin<=x<shape[1]-margin and
               margin<=y<shape[0]-margin
               for _,x,y in first)
    assert first==list(islice(
        _bounded_quadrant_blank_candidates(shape,margin,step)[1],80))
    # The prior implementation allocated tens of millions of (x,y)
    # tuples before considering any science mask. The V3 generator
    # only stores a pair of row blocks per native grid row.
    assert type(seq).__name__=="generator"
    assert seq.__sizeof__()<4096


def test_blank_sampling_is_spatially_stratified_not_first_clean_rows():
    """A capped clean field must sample distinct quadrants, not one edge."""
    rng=np.random.default_rng(29)
    values=rng.normal(0.,.3,(321,321))
    f=frame(image=values,shape=values.shape)
    none=np.zeros(f.shape,dtype=bool)
    result=estimate_blank_annulus_covariance(
        f,geometry(a=1.6,q=.8,pa=25.),(.5,1.,1.6),
        effective_source_mask=none,manual_exclude=none,
        max_placements=48)
    counts=result["accepted_centres_by_quadrant"]
    assert result["n_accepted"]==48
    assert sum(counts)==48
    assert result["n_quadrants_with_accepted_blanks"]==4
    assert min(counts)>=8
    assert max(counts)-min(counts)<=4
    assert result["candidate_order_independent_of_image_brightness"] is True
    assert result["placement_policy"]==(
        "STRATIFIED_QUADRANT_BOUNDED_MEMORY_NATIVE_RINGS_V3")
    assert result==estimate_blank_annulus_covariance(
        f,geometry(a=1.6,q=.8,pa=25.),(.5,1.,1.6),
        effective_source_mask=none,manual_exclude=none,
        max_placements=48)
    # The centre selection is entirely geometric and must not change
    # when the coadd pixel fluxes change while validity stays fixed.
    alternate=frame(image=-values,shape=values.shape)
    other=estimate_blank_annulus_covariance(
        alternate,geometry(a=1.6,q=.8,pa=25.),(.5,1.,1.6),
        effective_source_mask=none,manual_exclude=none,
        max_placements=48)
    assert other["centres_native_xy"]==result["centres_native_xy"]
    assert result["covariance_calibrated_for_science"] is False
    coords=np.asarray(result["centres_native_xy"],dtype=float)
    # All accepted footprints are separated by the original 2r+6 lattice.
    d2=((coords[:,None,:]-coords[None,:,:])**2).sum(axis=2)
    np.fill_diagonal(d2,np.inf)
    assert np.min(d2)>=(2*1.6)**2


def test_blank_covariance_flags_available_positions_in_one_quadrant():
    """Nominally numerous blanks can still give poor sky coverage."""
    pixwcs=wcs()
    pixwcs.wcs.crpix=[161.,161.]
    data=np.random.default_rng(91).normal(0.,.2,(321,321))
    f=frame(image=data,shape=data.shape,w=pixwcs)
    source=np.ones(f.shape,dtype=bool)
    source[:145,:145]=False
    result=estimate_blank_annulus_covariance(
        f,geometry(a=1.6),(.5,1.,1.6),
        effective_source_mask=source,
        manual_exclude=np.zeros(f.shape,dtype=bool),
        max_placements=48)
    assert result["n_accepted"] >= 12
    assert result["accepted_centres_by_quadrant"][0] == (
        result["n_accepted"])
    assert result["n_quadrants_with_accepted_blanks"] == 1
    assert "BLANK_APERTURE_SPATIAL_QUADRANT_COVERAGE_LIMITED" in (
        result["warnings"])
    assert "BLANK_APERTURE_SAMPLES_DOMINATED_BY_ONE_QUADRANT" in (
        result["warnings"])
    assert result["covariance_calibrated_for_science"] is False


def test_blank_apertures_ignore_nan_outside_exact_ring_support():
    """A masked NaN in a rectangular bbox with zero ring weight is harmless.

    Previously `weights * full_bbox_image` turned 0*NaN into NaN,
    poisoning otherwise valid blank covariance realizations.
    """
    from itertools import islice
    from lsbpipe.photometry.aperture_systematics import (
        _bounded_quadrant_blank_candidates)
    from lsbpipe.photometry.aperture_flux import project_icrs_ellipse_to_native
    original=np.random.default_rng(122).normal(0.,.3,(321,321))
    f0=frame(image=original,shape=original.shape)
    outer=geometry(a=1.6,q=1.,pa=0.)
    projected=project_icrs_ellipse_to_native(outer,f0.wcs)
    margin=int(np.ceil(projected.a_pix+3))
    step=max(2*margin+2,8)
    _,centres=_bounded_quadrant_blank_candidates(f0.shape,margin,step)
    with_nans=original.copy()
    for quadrant,x,y in islice(centres,48):
        assert x+6 < with_nans.shape[1] and y+6 < with_nans.shape[0]
        with_nans[y+6,x+6]=np.nan
    f=frame(image=with_nans,shape=with_nans.shape)
    none=np.zeros(f.shape,bool)
    result=estimate_blank_annulus_covariance(
        f,outer,(.5,1.,1.6),
        effective_source_mask=none,manual_exclude=none)
    assert result["n_accepted"]>=result["n_required"]
    assert result["status"]=="EMPIRICAL_BLANK_COVARIANCE_DIAGNOSTIC_ONLY"
    assert np.isfinite(result["annulus_covariance_nMgy2"]).all()
    assert result["covariance_calibrated_for_science"] is False
    assert np.isnan(with_nans).any()
    np.testing.assert_array_equal(f.image_nMgy,with_nans)


def test_missing_independent_blanks_refuses_covariance_instead_of_zero():
    f = frame(value=.1)
    source = np.ones(f.shape, dtype=bool)
    manual = np.zeros(f.shape, dtype=bool)
    kw = dict(
        frame=f, outer=geometry(a=6., q=1.),
        radii_arcsec=(1., 3., 6.),
        effective_source_mask=source, manual_exclude=manual)
    result = estimate_blank_annulus_covariance(**kw)
    assert result["status"] == "INSUFFICIENT_INDEPENDENT_BLANK_POSITIONS"
    assert result["n_accepted"] == 0
    assert result["annulus_covariance_nMgy2"] is None
    assert result["annulus_sigma_empirical_nMgy"] is None
    assert result["cumulative_covariance_nMgy2"] is None
    with pytest.raises(ValueError, match="Boolean"):
        estimate_blank_annulus_covariance(
            **dict(kw, manual_exclude=np.zeros((4, 4), bool)))


def test_constant_sky_response_and_mask_counterfactual_are_exact():
    f = frame(value=2.)
    source = np.zeros(f.shape, dtype=bool)
    source[56:63, 57:64] = True
    manual = np.zeros(f.shape, dtype=bool)
    manual[59:61, 60:62] = True
    radii = (1., 2.5, 4.)
    c1 = _curve(f, source, manual, radii=radii)
    c2 = _curve(f, source, manual, q=.7, pa=45., radii=radii)
    result = diagnose_aperture_systematics(
        f, anchor_icrs=(155.4, 18.), provenance=provenance(),
        effective_source_mask=source, manual_exclude=manual,
        curves=[c1, c2])
    derivative = result["sky_response"]["dflux_db0_native_pixel2"]
    assert derivative == pytest.approx(
        -c1["cumulative"][-1]["area_native_pixel2"][
            "valid_native_pixel2"], rel=1e-10)
    baseline = result["reference_valid_flux_nMgy"]
    b0 = .0013
    sky = ApertureSkyScenario(
        method="constant", b0_nMgy_per_pixel=b0,
        evidence_id="SYNTHETIC_FIXED_B0")
    corrected = measure_native_aperture(
        f, geometry(a=4., q=1., pa=0.), provenance(),
        effective_source_mask=source, manual_exclude=manual, sky=sky)
    assert corrected.signed_valid_flux_nMgy - baseline == (
        pytest.approx(b0 * derivative, abs=1e-8))
    mask = result["catalogue_mask_counterfactual"]
    unrestricted = measure_native_aperture(
        f, geometry(a=4., q=1., pa=0.), provenance(),
        effective_source_mask=np.zeros(f.shape, bool),
        manual_exclude=manual)
    assert mask["flux_without_provisional_source_mask_nMgy"] == (
        pytest.approx(unrestricted.signed_valid_flux_nMgy))
    assert mask["conditional_delta_nMgy"] > 0
    assert mask["delta_is_not_missing_galaxy_light"]
    assert result["geometry_hypothesis_span"]["n_valid_hypotheses"] == 2
    assert result["geometry_hypothesis_span"][
        "not_a_statistical_uncertainty"]
    assert not result["covariance_calibrated_for_science"]
    assert not result["adopted_for_science"]
    assert np.all(f.image_nMgy == 2.)
    json.dumps(result, allow_nan=False)


def test_tilted_sky_derivatives_use_target_center_native_coordinates():
    f = frame(value=.5)
    source = np.zeros(f.shape, dtype=bool)
    source[60:65, 65:69] = True  # asymmetry makes sky slopes observable
    manual = np.zeros(f.shape, dtype=bool)
    c1 = _curve(f, source, manual, q=.7, pa=27.)
    output = diagnose_aperture_systematics(
        f, anchor_icrs=(155.4, 18.), provenance=provenance(),
        effective_source_mask=source, manual_exclude=manual, curves=[c1])
    response = output["sky_response"]
    b0, bx, by = .002, .0004, -.0006
    scenario = ApertureSkyScenario(
        method="tilted_plane", b0_nMgy_per_pixel=b0,
        bx_nMgy_per_pixel2=bx, by_nMgy_per_pixel2=by,
        reference_xy_pix=tuple(response["reference_native_xy"]),
        evidence_id="SYNTHETIC_EXACT_PLANE")
    conditional = measure_native_aperture(
        f, geometry(a=4., q=.7, pa=27.), provenance(),
        effective_source_mask=source, manual_exclude=manual, sky=scenario)
    predicted = (
        b0 * response["dflux_db0_native_pixel2"] +
        bx * response["dflux_dbx_native_pixel3"] +
        by * response["dflux_dby_native_pixel3"])
    assert (conditional.signed_valid_flux_nMgy -
            output["reference_valid_flux_nMgy"]) == pytest.approx(
                predicted, abs=1e-8)
    assert response["sky_level_uncertainty_assumed"] is None
    assert response["sky_plane_uncertainty_assumed"] is None


def test_all_source_masked_is_underconstrained_not_a_flux_repair():
    f = frame(value=.8)
    source = np.ones(f.shape, dtype=bool)
    manual = np.zeros(f.shape, dtype=bool)
    c = _curve(f, source, manual)
    report = diagnose_aperture_systematics(
        f, anchor_icrs=(155.4, 18.), provenance=provenance(),
        effective_source_mask=source, manual_exclude=manual, curves=[c])
    assert report["reference_valid_flux_nMgy"] is None
    assert report["sky_response"]["dflux_db0_native_pixel2"] == 0.
    assert report["catalogue_mask_counterfactual"][
        "conditional_delta_nMgy"] is None
    assert report["blank_apertures"]["annulus_covariance_nMgy2"] is None
    assert report["adopted_for_science"] is False


def test_invalid_sampling_limits_and_radial_grid_fail_closed():
    f = frame()
    source = np.zeros(f.shape, bool)
    kw = dict(frame=f, outer=geometry(a=4.),
              radii_arcsec=(1., 2., 4.),
              effective_source_mask=source, manual_exclude=source)
    with pytest.raises(ValueError):
        estimate_blank_annulus_covariance(**dict(kw, max_placements=2))
    with pytest.raises(ValueError):
        estimate_blank_annulus_covariance(
            **dict(kw, radii_arcsec=(4., 2., 1.)))
    with pytest.raises(ValueError):
        estimate_blank_annulus_covariance(
            **dict(kw, radii_arcsec=(1., 2., 3.)))


def test_blank_covariance_sample_budget_tracks_radial_dimension():
    """The default cannot require more independent blanks than it searches."""
    f=frame()
    source=np.ones(f.shape,dtype=bool)
    radii=tuple(float(x) for x in np.linspace(.1,1.6,31))
    kw=dict(frame=f,outer=geometry(a=1.6),
            radii_arcsec=radii,effective_source_mask=source,
            manual_exclude=np.zeros(f.shape,dtype=bool))
    default=estimate_blank_annulus_covariance(**kw)
    assert default["n_required"]==64
    assert default["n_max_placements"]>=64
    assert not default["placement_cap_explicit"]
    assert default["status"]=="INSUFFICIENT_INDEPENDENT_BLANK_POSITIONS"
    assert default["annulus_covariance_nMgy2"] is None
    explicit=estimate_blank_annulus_covariance(**kw,max_placements=48)
    assert explicit["n_required"]==64
    assert explicit["n_max_placements"]==48
    assert explicit["placement_cap_explicit"] is True
    assert explicit["status"]==(
        "EXPLICIT_PLACEMENT_CAP_BELOW_REQUIRED_COVARIANCE_SAMPLE")
    assert explicit["annulus_covariance_nMgy2"] is None
    assert explicit["covariance_calibrated_for_science"] is False
    json.dumps(explicit,allow_nan=False)


def test_small_scale_blanks_keep_full_target_exclusion_and_never_replace_outer():
    """Smaller clean positions remain outside the FULL target aperture."""
    rng = np.random.default_rng(2035)
    w = wcs()
    w.wcs.crpix = [161., 161.]
    image = rng.normal(0., .03, (321, 321))
    f = frame(image=image, shape=image.shape, w=w)
    source = np.zeros(f.shape, bool)
    radii = (1., 2., 3., 4., 7., 25.)
    curve = _curve(f, source, source, radii=radii)
    result = diagnose_aperture_systematics(
        f, anchor_icrs=(155.4, 18.), provenance=provenance(),
        effective_source_mask=source, manual_exclude=source, curves=[curve])
    full = result["blank_apertures"]
    small = result["smaller_scale_blank_diagnostics"]
    assert full["annulus_covariance_nMgy2"] is None
    assert small["status"] == "SMALL_SCALE_COVARIANCE_DIAGNOSTIC_ONLY"
    assert small["must_not_replace_full_scale_covariance"]
    assert small["full_target_semimajor_arcsec"] == 25.
    assert small["first_available_probe_outer_arcsec"] in (2., 3., 4., 7.)
    probe = small["probes"][-1]
    assert probe["annulus_covariance_nMgy2"] is not None
    assert probe["n_accepted"] >= probe["n_required"]
    assert probe["target_exclusion_radius_native_pix"] > (
        probe["blank_outer_radius_native_pix"])
    from lsbpipe.photometry.aperture_flux import (
        project_icrs_ellipse_to_native)
    outer = project_icrs_ellipse_to_native(
        geometry(a=25., q=1., pa=0.), w)
    for x, y in probe["centres_native_xy"]:
        assert np.hypot(x-outer.x_center, y-outer.y_center) > (
            probe["target_plus_blank_exclusion_distance_native_pix"])
    assert result["covariance_calibrated_for_science"] is False
    assert small["adopted_for_science"] is False
    json.dumps(result, allow_nan=False)


def test_blank_subscale_refuses_mismatched_or_smaller_target_exclusion():
    f = frame()
    source = np.zeros(f.shape, bool)
    normal = dict(
        frame=f, outer=geometry(a=2., q=1.),
        radii_arcsec=(1., 2.), effective_source_mask=source,
        manual_exclude=source)
    with pytest.raises(ValueError, match="share project anchor"):
        estimate_blank_annulus_covariance(
            **normal, full_target_exclusion_outer=geometry(
                ra=155.401, a=5., q=1.))
    with pytest.raises(ValueError, match="cannot be smaller"):
        estimate_blank_annulus_covariance(
            **normal, full_target_exclusion_outer=geometry(a=1., q=1.))


def test_nonzero_sky_is_held_fixed_in_catalogue_mask_counterfactual():
    """Turning off a provisional mask must never also turn off sky subtraction."""
    f = frame(value=2.)
    source = np.zeros(f.shape, dtype=bool)
    source[57:62, 57:62] = True
    manual = np.zeros(f.shape, dtype=bool)
    scenario = ApertureSkyScenario(
        method="tilted_plane", b0_nMgy_per_pixel=.15,
        bx_nMgy_per_pixel2=.001, by_nMgy_per_pixel2=-.002,
        reference_xy_pix=(60., 60.), evidence_id="KNOWN_SYNTHETIC_PLANE")
    curve = _curve(f, source, manual, sky=scenario)
    report = diagnose_aperture_systematics(
        f, anchor_icrs=(155.4, 18.), provenance=provenance(),
        effective_source_mask=source, manual_exclude=manual,
        curves=[curve], sky=scenario)
    direct = measure_native_aperture(
        f, geometry(a=4., q=1., pa=0.), provenance(),
        effective_source_mask=np.zeros(f.shape, dtype=bool),
        manual_exclude=manual, sky=scenario)
    masked = curve["cumulative"][-1]["signed_valid_flux_nMgy"]
    counter = report["catalogue_mask_counterfactual"]
    assert counter["flux_with_provisional_source_mask_nMgy"] == pytest.approx(masked)
    assert counter["flux_without_provisional_source_mask_nMgy"] == pytest.approx(
        direct.signed_valid_flux_nMgy)
    assert counter["conditional_delta_nMgy"] == pytest.approx(
        direct.signed_valid_flux_nMgy - masked)
    assert report["conditional_sky_scenario"]["method"] == "tilted_plane"
    assert report["stage_04_2_receipts_and_outer_apertures_verified"] == 1
    assert "BLANK_COVARIANCE_USES_UNCORRECTED_NATIVE_COADD" in report["warnings"]
    with pytest.raises(ValueError, match="sky scenario|recipe mismatch"):
        diagnose_aperture_systematics(
            f, anchor_icrs=(155.4, 18.), provenance=provenance(),
            effective_source_mask=source, manual_exclude=manual,
            curves=[curve])


def test_stage_04_3_rejects_same_area_modified_source_pixels_and_stale_image():
    """Equal mask area and matching WCS do not certify matching input pixels."""
    f = frame(value=2.)
    source = np.zeros(f.shape, dtype=bool)
    source[56:60, 57:61] = True
    manual = np.zeros(f.shape, dtype=bool)
    curve = _curve(f, source, manual)
    modified_source = source.copy()
    modified_source[56:60, 57:61] = False
    modified_source[60:64, 60:64] = True
    assert np.count_nonzero(source) == np.count_nonzero(modified_source)
    args = dict(anchor_icrs=(155.4, 18.), provenance=provenance(),
                manual_exclude=manual, curves=[curve])
    with pytest.raises(ValueError, match="recipe mismatch"):
        diagnose_aperture_systematics(
            f, effective_source_mask=modified_source, **args)
    altered = frame(value=3.)
    with pytest.raises(ValueError, match="native pixels or weights changed"):
        diagnose_aperture_systematics(
            altered, effective_source_mask=source, **args)


def test_stage_04_3_rejects_altered_receipts_and_mixed_growth_lineages():
    f = frame(value=.9)
    source = np.zeros(f.shape, bool)
    manual = np.zeros(f.shape, bool)
    curve = _curve(f, source, manual)
    args = dict(frame=f, anchor_icrs=(155.4, 18.),
                provenance=provenance(), effective_source_mask=source,
                manual_exclude=manual)
    corrupted = dict(curve)
    corrupted["source_policy_id"] = "DIFFERENT_MASK_POLICY"
    with pytest.raises(ValueError, match="receipt integrity"):
        diagnose_aperture_systematics(curves=[corrupted], **args)
    other = _curve(f, source, manual, radii=(1., 3., 4.))
    with pytest.raises(ValueError, match="radial grid"):
        diagnose_aperture_systematics(curves=[curve, other], **args)
    good = diagnose_aperture_systematics(curves=[curve], **args)
    assert good["stage_04_2_receipts_and_outer_apertures_verified"] == 1
    assert good["conditional_sky_scenario"]["method"] == "zero_additional_correction"
