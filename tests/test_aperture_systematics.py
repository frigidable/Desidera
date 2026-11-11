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


def _curve(f, source, manual, *, q=1., pa=0., radii=(1., 2.5, 4.)):
    return measure_growth_diagnostic(
        f, anchor_icrs=(155.4, 18.), provenance=provenance(),
        effective_source_mask=source, manual_exclude=manual,
        radii_arcsec=radii,
        hypothesis=GrowthHypothesis("fixed_test", q, pa, "SYNTHETIC"))


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
