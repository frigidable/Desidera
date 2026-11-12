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


def test_tilted_sky_derivatives_use_declared_frozen_sky_reference():
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
    assert "BLANK_COVARIANCE_CONDITIONAL_ON_UNADOPTED_SKY_PLANE" in report["warnings"]
    assert report["blank_apertures"]["conditional_sky_scenario"] == (
        report["conditional_sky_scenario"])
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
    with pytest.raises(ValueError, match="native pixels or weights changed|native_input_sha256"):
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


def test_compensated_native_pixel_permutations_fail_closed():
    """Matching whole-aperture sums cannot prove unchanged IMAGE or INVVAR."""
    image = np.ones((121, 121), dtype=float)
    image[60, 60], image[60, 61] = 2., 4.
    base = frame(image=image)
    zeros = np.zeros(base.shape, dtype=bool)
    curve = _curve(base, zeros, zeros)
    switched = image.copy()
    switched[60, 60], switched[60, 61] = 4., 2.
    changed = frame(image=switched)
    args = dict(anchor_icrs=(155.4, 18.), provenance=provenance(),
                effective_source_mask=zeros, manual_exclude=zeros,
                curves=[curve])
    from lsbpipe.photometry.aperture_flux import measure_native_aperture
    def outer(f):
        return measure_native_aperture(
            f, geometry(a=4., q=1., pa=0.), provenance(),
            effective_source_mask=zeros, manual_exclude=zeros)
    first, second = outer(base), outer(changed)
    assert first.recipe_sha256 == second.recipe_sha256
    assert first.signed_valid_flux_nMgy == pytest.approx(
        second.signed_valid_flux_nMgy)
    assert first.formal_diagonal_sigma_nMgy == pytest.approx(
        second.formal_diagonal_sigma_nMgy)
    assert base.native_input_sha256 != changed.native_input_sha256
    with pytest.raises(ValueError, match="native_input_sha256"):
        diagnose_aperture_systematics(changed, **args)

    initial_ivar = np.full(base.shape, 4.)
    initial_ivar[60, 60], initial_ivar[60, 61] = 2., 6.
    variance_frame = frame(image=image, ivar=initial_ivar)
    variance_curve = _curve(variance_frame, zeros, zeros)
    changed_ivar = initial_ivar.copy()
    changed_ivar[60, 60], changed_ivar[60, 61] = 6., 2.
    changed_variance = frame(image=image, ivar=changed_ivar)
    a, b = outer(variance_frame), outer(changed_variance)
    assert a.formal_diagonal_sigma_nMgy == pytest.approx(
        b.formal_diagonal_sigma_nMgy)
    assert a.recipe_sha256 == b.recipe_sha256
    with pytest.raises(ValueError, match="native_input_sha256"):
        diagnose_aperture_systematics(
            changed_variance, **dict(args, curves=[variance_curve]))
    assert variance_frame.native_input_sha256 == variance_frame.native_input_sha256


@pytest.mark.parametrize("bad_weight", (-0.2, 1.2, float("nan")))
def test_stage_04_3_rejects_invalid_exact_weights_in_sky_derivative(
        monkeypatch, bad_weight):
    """04.3 must reject bad overlap fractions even after valid 04.2 receipts."""
    import photutils.aperture as ap
    baseline_class = ap.EllipticalAperture
    f = frame(value=1.)
    zeros = np.zeros(f.shape, dtype=bool)
    curve = _curve(f, zeros, zeros)
    calls = [0]

    class CorruptedSecondAperture(baseline_class):
        def to_mask(self, method="exact"):
            result = super().to_mask(method=method)
            calls[0] += 1
            # First use is the independent fresh Stage 04.2 aperture
            # verification; the second use is the Stage 04.3 derivative.
            if calls[0] == 2:
                result.data[result.data.shape[0]//2,
                            result.data.shape[1]//2] = bad_weight
            return result

    monkeypatch.setattr(ap, "EllipticalAperture", CorruptedSecondAperture)
    with pytest.raises(ValueError, match="Invalid Stage 04.3 exact"):
        diagnose_aperture_systematics(
            f, anchor_icrs=(155.4, 18.),
            provenance=provenance(), effective_source_mask=zeros,
            manual_exclude=zeros, curves=[curve])
    assert calls[0] == 2


def test_stage_04_3_sky_gradient_derivatives_reproduce_direct_flux_for_offcenter_reference():
    """An off-centre sky-plane pivot must not be replaced by the target centre."""
    f = frame(value=.6)
    source = np.zeros(f.shape, dtype=bool)
    source[59:64, 64:70] = True
    manual = np.zeros(f.shape, dtype=bool)
    base = ApertureSkyScenario(
        method="tilted_plane", b0_nMgy_per_pixel=.03,
        bx_nMgy_per_pixel2=.0004, by_nMgy_per_pixel2=-.0009,
        reference_xy_pix=(7.5, 92.25),
        evidence_id="SYNTHETIC_FIXED_REFERENCE")
    curve = _curve(f, source, manual, q=.7, pa=27., sky=base)
    report = diagnose_aperture_systematics(
        f, anchor_icrs=(155.4, 18.), provenance=provenance(),
        effective_source_mask=source, manual_exclude=manual,
        curves=[curve], sky=base)
    response = report["sky_response"]
    assert response["reference_native_xy"] == [7.5, 92.25]
    assert response["gradient_reference_matches_frozen_sky_scenario"] is True
    assert response["aperture_center_native_xy"] != [7.5, 92.25]
    db0, dbx, dby = .002, .0003, -.0002
    changed = ApertureSkyScenario(
        method="tilted_plane", b0_nMgy_per_pixel=base.b0_nMgy_per_pixel+db0,
        bx_nMgy_per_pixel2=base.bx_nMgy_per_pixel2+dbx,
        by_nMgy_per_pixel2=base.by_nMgy_per_pixel2+dby,
        reference_xy_pix=base.reference_xy_pix,
        evidence_id=base.evidence_id)
    measured = measure_native_aperture(
        f, geometry(a=4., q=.7, pa=27.), provenance(),
        effective_source_mask=source, manual_exclude=manual, sky=changed)
    observed_change = (measured.signed_valid_flux_nMgy -
                       report["reference_valid_flux_nMgy"])
    predicted_change = (
        db0*response["dflux_db0_native_pixel2"] +
        dbx*response["dflux_dbx_native_pixel3"] +
        dby*response["dflux_dby_native_pixel3"])
    assert observed_change == pytest.approx(predicted_change, abs=1e-8)
    # A target-centred derivative is a legitimate alternative
    # parameterization, but is NOT the actual stored scenario derivative.
    assert abs(response["dflux_dbx_native_pixel3"] -
               response["dflux_dbx_at_aperture_center_native_pixel3"]) > 1.
    assert abs(response["dflux_dby_native_pixel3"] -
               response["dflux_dby_at_aperture_center_native_pixel3"]) > 1.
    json.dumps(report, allow_nan=False)


def test_radial_sky_jacobian_predicts_every_native_growth_aperture_and_ring():
    """Each saved curve point has the SAME frozen sky reference as direct flux."""
    f = frame(value=1.7)
    source = np.zeros(f.shape, dtype=bool)
    source[58:65, 61:73] = True
    manual = np.zeros(f.shape, dtype=bool)
    manual[62:70, 49:55] = True
    radii = (.7, 1.5, 2.8, 4.)
    base = ApertureSkyScenario(
        method="tilted_plane", b0_nMgy_per_pixel=.05,
        bx_nMgy_per_pixel2=.0013,
        by_nMgy_per_pixel2=-.0008,
        reference_xy_pix=(9.25, 104.5),
        evidence_id="FIXED_OFFCENTRE_SKY")
    curve = _curve(f, source, manual, q=.73, pa=38.,
                   radii=radii, sky=base)
    result = diagnose_aperture_systematics(
        f, anchor_icrs=(155.4, 18.), provenance=provenance(),
        effective_source_mask=source, manual_exclude=manual,
        curves=[curve], sky=base)
    radial = result["radial_sky_response"]
    assert radial["sky_reference_native_xy"] == [9.25, 104.5]
    assert radial["sky_parameter_covariance"] is None
    assert radial["confidence_interval_or_sigma_calibrated"] is False
    assert radial["all_frozen_radial_areas_independently_verified"]
    assert len(radial["cumulative"]) == len(radii)
    assert len(radial["disjoint_annuli"]) == len(radii)-1
    db0, dbx, dby = (.011, .0007, -.0004)
    changed = ApertureSkyScenario(
        method="tilted_plane", b0_nMgy_per_pixel=.05+db0,
        bx_nMgy_per_pixel2=.0013+dbx,
        by_nMgy_per_pixel2=-.0008+dby,
        reference_xy_pix=base.reference_xy_pix,
        evidence_id=base.evidence_id)
    perturbed_fluxes = []
    for i, radius in enumerate(radii):
        row = radial["cumulative"][i]
        direct = measure_native_aperture(
            f, geometry(a=radius, q=.73, pa=38.),
            provenance(), effective_source_mask=source,
            manual_exclude=manual, sky=changed)
        expected = sum(delta*derivative for delta, derivative in zip(
            (db0, dbx, dby), row["derivative_vector_for_b0_bx_by"]))
        observed = direct.signed_valid_flux_nMgy - (
            curve["cumulative"][i]["signed_valid_flux_nMgy"])
        assert observed == pytest.approx(expected, abs=2.e-8)
        assert row["valid_native_pixel2"] == pytest.approx(
            curve["cumulative"][i]["area_native_pixel2"][
                "valid_native_pixel2"], abs=1.e-6)
        perturbed_fluxes.append(direct.signed_valid_flux_nMgy)
    for i, ring in enumerate(radial["disjoint_annuli"], start=1):
        change = (
            perturbed_fluxes[i] - perturbed_fluxes[i-1] -
            curve["cumulative"][i]["signed_valid_flux_nMgy"] +
            curve["cumulative"][i-1]["signed_valid_flux_nMgy"])
        expected = sum(delta*derivative for delta, derivative in zip(
            (db0, dbx, dby), ring["derivative_vector_for_b0_bx_by"]))
        assert change == pytest.approx(expected, abs=2.e-8)
        assert ring["valid_native_pixel2"] == pytest.approx(
            curve["independently_integrated_annuli"][i-1][
                "valid_area_pixel2"], abs=1.e-6)
    assert (radial["cumulative"][-1][
            "derivative_vector_for_b0_bx_by"][1] ==
            pytest.approx(result["sky_response"][
                "dflux_dbx_native_pixel3"]))
    json.dumps(result, allow_nan=False)


def test_radial_sky_jacobian_reports_missing_area_without_false_errorbars():
    f = frame(value=0.5)
    source = np.ones(f.shape, dtype=bool)
    zeros = np.zeros(f.shape, dtype=bool)
    curve = _curve(f, source, zeros, radii=(1., 2., 4.))
    result = diagnose_aperture_systematics(
        f, anchor_icrs=(155.4, 18.), provenance=provenance(),
        effective_source_mask=source, manual_exclude=zeros,
        curves=[curve])
    radial = result["radial_sky_response"]
    assert all(row["status"] ==
               "NO_VALID_NATIVE_PIXELS_NO_SKY_INFERENCE"
               for row in radial["cumulative"])
    assert all(row["derivative_vector_for_b0_bx_by"] ==
               [0., 0., 0.] for row in radial["cumulative"])
    assert radial["sky_parameter_covariance"] is None
    assert not radial["confidence_interval_or_sigma_calibrated"]


def test_radial_source_only_mask_ablation_matches_exact_native_reintegration():
    """A galaxy knot hidden by an automatic mask is measurable sensitivity,
    not proof that the knot belongs to the galaxy or should be restored.
    """
    f = frame(value=1.8)
    source = np.zeros(f.shape, dtype=bool)
    source[55:65, 62:70] = True
    manual = np.zeros(f.shape, dtype=bool)
    manual[57:61, 63:67] = True
    original_ivar = np.full(f.shape, 4.)
    original_ivar[58:60, 68:70] = 0.
    f = frame(value=1.8, ivar=original_ivar)
    radii = (1., 2., 3.5, 5.)
    sky = ApertureSkyScenario(
        method="tilted_plane",
        b0_nMgy_per_pixel=.02,
        bx_nMgy_per_pixel2=.0005,
        by_nMgy_per_pixel2=-.0007,
        reference_xy_pix=(19., 91.),
        evidence_id="FROZEN_MASK_ALTERNATIVE")
    curve = _curve(f, source, manual, q=.75, pa=32.,
                   radii=radii, sky=sky)
    report = diagnose_aperture_systematics(
        f, anchor_icrs=(155.4, 18.), provenance=provenance(),
        effective_source_mask=source, manual_exclude=manual,
        curves=[curve], sky=sky)
    radial = report["radial_catalogue_mask_counterfactual"]
    assert radial["outer_native_reintegration_verified"]
    assert radial["missing_target_light_not_reconstructed"]
    assert radial["hard_and_manual_masks_unchanged"]
    assert radial["adopted_for_science"] is False
    deltas = []
    for i, radius in enumerate(radii):
        row = radial["cumulative"][i]
        relaxed = measure_native_aperture(
            f, geometry(a=radius, q=.75, pa=32.), provenance(),
            effective_source_mask=np.zeros(f.shape, dtype=bool),
            manual_exclude=manual, sky=sky)
        baseline = curve["cumulative"][i]["signed_valid_flux_nMgy"]
        delta = relaxed.signed_valid_flux_nMgy - baseline
        assert row["signed_flux_newly_available_nMgy"] == pytest.approx(
            delta, abs=1.e-7)
        assert row["conditional_flux_with_source_mask_removed_nMgy"] == (
            pytest.approx(relaxed.signed_valid_flux_nMgy, abs=1.e-7))
        assert row["area_newly_available_source_only_native_pixel2"] >= 0.
        deltas.append(delta)
    assert deltas[-1] > 0.
    for i, ring in enumerate(radial["disjoint_annuli"], start=1):
        assert ring["signed_flux_newly_available_nMgy"] == pytest.approx(
            deltas[i] - deltas[i-1], abs=1.e-7)
    json.dumps(report, allow_nan=False)


def test_radial_source_mask_ablation_preserves_signed_negative_contributions():
    image = np.ones((121, 121))
    image[60:63, 62:66] = -2.
    f = frame(image=image)
    source = np.zeros(f.shape, dtype=bool)
    source[60:63, 62:66] = True
    manual = np.zeros(f.shape, dtype=bool)
    curve = _curve(f, source, manual, q=1., pa=0.,
                   radii=(1., 2.5, 4.))
    out = diagnose_aperture_systematics(
        f, anchor_icrs=(155.4, 18.), provenance=provenance(),
        effective_source_mask=source, manual_exclude=manual,
        curves=[curve])
    radial = out["radial_catalogue_mask_counterfactual"]
    assert radial["cumulative"][-1][
        "signed_flux_newly_available_nMgy"] < 0.
    assert out["catalogue_mask_counterfactual"][
        "conditional_delta_nMgy"] < 0.
    assert out["adopted_for_science"] is False


def test_blank_covariance_respects_exact_same_conditional_plane_as_target():
    """A known tilted plane must be subtracted at EACH blank native pixel."""
    from lsbpipe.photometry.aperture_systematics import (
        estimate_blank_annulus_covariance)
    w = wcs()
    w.wcs.crpix = [161., 161.]
    yy, xx = np.indices((321, 321))
    rng = np.random.default_rng(290187)
    noise = rng.normal(0., .007, (321, 321))
    b0, bx, by = .13, .004, -.002
    origin = (23.5, 197.25)
    image = (b0 + bx*(xx-origin[0]) +
             by*(yy-origin[1]) + noise)
    f = frame(image=image, shape=image.shape, w=w)
    zeros = np.zeros(f.shape, dtype=bool)
    kw = dict(frame=f, outer=geometry(a=1.6, q=.77, pa=31.),
              radii_arcsec=(.4, 1., 1.6),
              effective_source_mask=zeros, manual_exclude=zeros)
    raw = estimate_blank_annulus_covariance(**kw)
    sky = ApertureSkyScenario(
        method="tilted_plane",
        b0_nMgy_per_pixel=b0,
        bx_nMgy_per_pixel2=bx,
        by_nMgy_per_pixel2=by,
        reference_xy_pix=origin,
        evidence_id="KNOWN_NATIVE_PLANE")
    residual = estimate_blank_annulus_covariance(**kw, sky=sky)
    assert raw["status"] == residual["status"] == (
        "EMPIRICAL_BLANK_COVARIANCE_DIAGNOSTIC_ONLY")
    assert residual["schema"] == (
        "desidera.photometry.blank-annulus-covariance/3")
    assert residual["centres_native_xy"] == raw["centres_native_xy"]
    assert residual["conditional_sky_scenario"]["reference_xy_pix"] == origin or (
        residual["conditional_sky_scenario"]["reference_xy_pix"] == list(origin))
    assert residual["blank_signal_uses_same_frozen_sky_plane_as_target"]
    assert residual["sky_plane_not_calibrated_or_adopted"] is True
    # The large plane gradient contributes covariance across blank
    # locations in the raw image; removing the known plane greatly
    # lowers the conditional residual variance.
    raw_diag = np.diag(raw["annulus_covariance_nMgy2"])
    residual_diag = np.diag(residual["annulus_covariance_nMgy2"])
    assert np.all(residual_diag < raw_diag*.1)
    assert np.max(np.abs(residual["annulus_mean_blank_flux_nMgy"])) < .03
    assert residual["covariance_calibrated_for_science"] is False
    json.dumps(residual, allow_nan=False)


def test_smaller_scale_blank_probe_uses_same_nonadopted_sky_scenario():
    w = wcs()
    w.wcs.crpix = [161., 161.]
    rng = np.random.default_rng(20261113)
    image = 1. + rng.normal(0., .002, (321, 321))
    f = frame(image=image, shape=(321, 321), w=w)
    masks = np.zeros(f.shape, dtype=bool)
    sky = ApertureSkyScenario(
        method="constant", b0_nMgy_per_pixel=1.,
        evidence_id="SYNTHETIC_CONDITIONAL")
    curve = _curve(f, masks, masks, radii=(1., 2., 3., 4., 7., 25.),
                   sky=sky)
    report = diagnose_aperture_systematics(
        f, anchor_icrs=(155.4, 18.), provenance=provenance(),
        effective_source_mask=masks, manual_exclude=masks,
        curves=[curve], sky=sky)
    assert report["blank_apertures"]["conditional_sky_scenario"][
        "b0_nMgy_per_pixel"] == 1.
    small = report["smaller_scale_blank_diagnostics"]
    assert small["status"] == "SMALL_SCALE_COVARIANCE_DIAGNOSTIC_ONLY"
    assert len(small["probes"]) >= 1
    assert small["probes"][-1]["annulus_covariance_nMgy2"] is not None
    assert all(p["conditional_sky_scenario"]["method"] == "constant"
               for p in small["probes"])
    assert all(p["blank_signal_uses_same_frozen_sky_plane_as_target"]
               for p in small["probes"])
    assert not report["covariance_calibrated_for_science"]


def test_stage043_rejects_re_receipted_inner_flux_and_variance_corruption():
    """A matching outer point and coherent receipt do not verify inner data."""
    from copy import deepcopy
    from hashlib import sha256
    f = frame(value=.31)
    source = np.zeros(f.shape, dtype=bool)
    source[63:66, 61:67] = True
    manual = np.zeros(f.shape, dtype=bool)
    curve = _curve(f, source, manual, q=.72, pa=31.,
                   radii=(.8, 1.7, 2.9, 4.))
    opts = dict(
        frame=f, anchor_icrs=(155.4, 18.), provenance=provenance(),
        effective_source_mask=source, manual_exclude=manual)
    clean = diagnose_aperture_systematics(**opts, curves=[curve])
    radial = clean["radial_sky_response"]
    assert radial["all_frozen_radial_fluxes_independently_verified"]
    assert radial["all_frozen_radial_diagonal_sigmas_independently_verified"]
    assert radial["no_unverified_internal_radial_samples"]
    for field, shift in (("signed_valid_flux_nMgy", .75),
                         ("formal_diagonal_sigma_nMgy", .15)):
        forged = deepcopy(curve)
        assert forged["cumulative"][1][field] is not None
        forged["cumulative"][1][field] += shift
        # A caller can recompute a receipt of an internally inconsistent
        # Stage 04.2 mapping. The frozen *outer* aperture is unchanged.
        forged["receipt_sha256"] = sha256(json.dumps(
            {k: v for k, v in forged.items() if k != "receipt_sha256"},
            sort_keys=True, separators=(",", ":"),
            allow_nan=False).encode()).hexdigest()
        with pytest.raises(ValueError, match="stale frozen Stage 04.2 radial"):
            diagnose_aperture_systematics(**opts, curves=[forged])
    # Underconstrained radial bins must preserve None rather than be
    # silently replaced by zero or a fake positive flux.
    all_masked = np.ones(f.shape, dtype=bool)
    none_curve = _curve(f, all_masked, manual,
                        radii=(1., 2.5, 4.))
    all_qa = diagnose_aperture_systematics(
        **dict(opts, effective_source_mask=all_masked),
        curves=[none_curve])
    assert all_qa["radial_sky_response"][
        "all_frozen_radial_fluxes_independently_verified"]
    assert all(r["frozen_valid_flux_nMgy"] is None for r in all_qa[
        "radial_sky_response"]["cumulative"])


def test_blank_native_centres_match_target_subpixel_phase_and_noise_weights():
    """Equal native-grid phase is required for comparable exact pixel weights."""
    from photutils.aperture import EllipticalAperture
    from lsbpipe.photometry.aperture_flux import project_icrs_ellipse_to_native
    pixel_wcs = wcs()
    pixel_wcs.wcs.crpix = [161., 161.]
    f = frame(
        image=np.random.default_rng(20261113).normal(0., .2, (321,321)),
        shape=(321,321), w=pixel_wcs)
    none = np.zeros(f.shape, dtype=bool)
    target = geometry(
        ra=155.400037, dec=18.000041, a=1.6, q=.67, pa=37.)
    projected = project_icrs_ellipse_to_native(target, f.wcs)
    result = estimate_blank_annulus_covariance(
        f, target, (.6, 1.1, 1.6),
        effective_source_mask=none, manual_exclude=none)
    assert result["status"] == "EMPIRICAL_BLANK_COVARIANCE_DIAGNOSTIC_ONLY"
    assert result["blank_centre_phase_matches_target_native_xy"]
    phase = result["target_native_subpixel_phase_xy"]
    assert 0. < phase[0] < 1. and 0. < phase[1] < 1.
    assert phase == pytest.approx([
        projected.x_center % 1., projected.y_center % 1.])
    target_ap = EllipticalAperture(
        (projected.x_center, projected.y_center),
        projected.a_pix, projected.b_pix,
        theta=projected.theta_rad_pixel_ccw)
    target_w = target_ap.to_mask(method="exact").data
    target_weight_square = float(np.sum(target_w**2))
    # The formal white-noise variance of the same-size aperture would
    # depend on sum(w^2), not just sum(w) or ellipse geometric area.
    for x, y in result["centres_native_xy"]:
        assert x % 1. == pytest.approx(phase[0], abs=1.e-10)
        assert y % 1. == pytest.approx(phase[1], abs=1.e-10)
        blank_ap = EllipticalAperture(
            (x, y), projected.a_pix, projected.b_pix,
            theta=projected.theta_rad_pixel_ccw)
        blank_w = blank_ap.to_mask(method="exact").data
        assert np.sum(blank_w**2) == pytest.approx(
            target_weight_square, abs=2.e-7)
    assert not result["covariance_calibrated_for_science"]
    assert result["no_science_covariance_calibration_from_phase_matching"]
    json.dumps(result, allow_nan=False)
