"""Stage 04.2 exact native finite curves of growth: fixed shapes, no total."""
import json

import numpy as np
import pytest

pytest.importorskip("astropy")
pytest.importorskip("photutils")

from test_aperture_flux import frame, geometry, provenance
from lsbpipe.photometry.aperture_growth import (
    GrowthHypothesis, geometry_hypotheses, predeclared_growth_radii,
    measure_growth_diagnostic, measure_native_annulus,
)
from lsbpipe.photometry.aperture_flux import (
    ApertureGeometryICRS, ApertureSkyScenario, measure_native_aperture,
)


def run_curve(f, *, source=None, manual=None, radii=(1.2,2.7,4.5,7.4),
              q=.72, pa=34., sky=None):
    if source is None:
        source=np.zeros(f.shape,bool)
    if manual is None:
        manual=np.zeros(f.shape,bool)
    return measure_growth_diagnostic(
        f, anchor_icrs=(155.4,18.),provenance=provenance(),
        effective_source_mask=source, manual_exclude=manual,
        radii_arcsec=tuple(radii),hypothesis=GrowthHypothesis(
            "test_ellipse",q,pa,"SYNTHETIC_TEST"),
        sky=sky)


def test_fixed_growth_constant_image_exact_ring_identity_and_area():
    f=frame(value=.35)
    report=run_curve(f)
    assert report["schema"]=="desidera.photometry.native-growth-diagnostic/1"
    assert report["science_status"]=="DIAGNOSTIC_MEASURED_PARTIAL_NOT_ADOPTED"
    assert not report["adopted_for_science"]
    assert not report["masked_galaxy_light_restored"]
    assert len(report["cumulative"])==4
    assert len(report["independently_integrated_annuli"])==3
    assert report["receipt_sha256"] and len(report["receipt_sha256"])==64
    for c in report["cumulative"]:
        assert c["sky_method"]=="zero_additional_correction"
        assert c["signed_valid_flux_nMgy"]==pytest.approx(
            .35*c["area_native_pixel2"]["valid_native_pixel2"],abs=1e-7)
    for ring in report["independently_integrated_annuli"]:
        assert ring["signed_flux_nMgy"]==pytest.approx(
            .35*ring["valid_area_pixel2"],abs=1e-7)
        assert ring["valid_fraction_of_geometric"]==pytest.approx(1,rel=1e-4)
        assert ring["formal_diagonal_sigma_nMgy"]>0.
    assert report["max_abs_annulus_closure_nMgy"]<1.e-6
    json.dumps(report,allow_nan=False)
    other=run_curve(f)
    assert other["receipt_sha256"]==report["receipt_sha256"]


def test_growth_never_forces_negative_outer_annulus_to_be_monotone():
    f=frame(value=0.)
    arr=np.full(f.shape,-.1,dtype=float)
    cy,cx=60,60
    yy,xx=np.indices(f.shape)
    arr[np.hypot(xx-cx,yy-cy)<15]=1.
    changed=frame(image=arr)
    out=run_curve(changed,radii=(1.,2.,3.,4.,6.),q=1.,pa=0.)
    fluxes=[c["signed_valid_flux_nMgy"] for c in out["cumulative"]]
    assert fluxes[1]>fluxes[0]
    assert fluxes[-1]<fluxes[-2]
    assert out["independently_integrated_annuli"][-1][
        "signed_flux_nMgy"]<0
    assert out["independently_integrated_annuli"][-1][
        "ab_surface_brightness_mag_arcsec2"] is None
    assert out["max_abs_annulus_closure_nMgy"]<1e-6


def test_annulus_masks_and_nonzero_conditional_sky_propagate_exactly():
    f=frame(value=1.)
    masks=np.zeros(f.shape,bool)
    masks[54:65,55:63]=True
    excluded=np.zeros(f.shape,bool)
    excluded[56:63,56:67]=True
    sky=ApertureSkyScenario(
        method="constant",b0_nMgy_per_pixel=.02,
        evidence_id="UNADOPTED_STAGE03_SENSITIVITY")
    r=run_curve(f,source=masks,manual=excluded,sky=sky)
    for annulus in r["independently_integrated_annuli"]:
        if annulus["signed_flux_nMgy"] is not None:
            assert annulus["signed_flux_nMgy"]==pytest.approx(
                .98*annulus["valid_area_pixel2"],abs=1e-6)
    assert r["cumulative"][-1]["area_native_pixel2"][
        "valid_fraction_of_ideal"]<1
    assert r["max_abs_annulus_closure_nMgy"]<1e-6
    assert r["cumulative"][-1]["sky_science_adopted"] is False


def test_frozen_anchor_radii_are_not_measured_isophotes():
    f=frame()
    radii,meta=predeclared_growth_radii(
        f,(155.4,18.),target_protection_native_pixels=650)
    assert len(radii)==10
    assert all(a<b for a,b in zip(radii,radii[1:]))
    assert meta["radii_are_not_measured_isophotes"]
    assert meta["target_protection_is_not_photometry_mask"]
    assert "PROVISIONAL_TARGET_PROTECTION" in meta["basis"]
    assert radii[-1] < 16.
    fallback,info=predeclared_growth_radii(f,(155.4,18.))
    assert info["basis"].startswith("NATIVE_FIELD_SCALE_ONLY")
    assert len(fallback)==10
    with pytest.raises(ValueError):
        predeclared_growth_radii(
            f,(155.4,18.),target_protection_native_pixels=0)


def test_elliptic_hypotheses_are_only_anchor_shape_proposals():
    fallback=geometry_hypotheses()
    assert len(fallback)==3
    assert fallback[0].q==1.
    assert fallback[0].pa_deg==0.
    assert len({x.name for x in fallback})==3
    sga=geometry_hypotheses(sga_q=.4,sga_pa_deg=210)
    assert sga[0].name=="anchor_circle"
    assert sga[1].name=="sga_q_pa_candidate"
    assert sga[1].pa_deg==30.
    assert len(sga)==2
    with pytest.raises(ValueError):
        GrowthHypothesis("bad",0.,90.,"SYNTHETIC")


def test_invalid_radial_grids_and_annular_geometry_fail_closed():
    f=frame()
    for radii in ((1.,),(3.,2.),(1.,1.),(0.,2.),
                  (float("nan"),3.)):
        with pytest.raises(ValueError):
            run_curve(f,radii=radii)
    g=geometry(a=5.,q=.7)
    mask=np.zeros(f.shape,bool)
    with pytest.raises(ValueError):
        measure_native_annulus(
            f,g,6.,provenance(),effective_source_mask=mask,
            manual_exclude=mask,sky=ApertureSkyScenario())
    with pytest.raises(ValueError):
        measure_native_annulus(
            f,g,1.,provenance(),
            effective_source_mask=np.zeros((4,4),bool),
            manual_exclude=mask,sky=ApertureSkyScenario())


def test_annulus_no_valid_pixels_does_not_invent_flux():
    f=frame()
    allmasked=np.ones(f.shape,bool)
    annulus=measure_native_annulus(
        f,geometry(a=6.,q=1.),2.,provenance(),
        effective_source_mask=allmasked,
        manual_exclude=np.zeros(f.shape,bool),
        sky=ApertureSkyScenario())
    assert annulus.no_valid_pixels is True
    assert annulus.signed_flux_nMgy is None
    assert annulus.ab_surface_brightness_mag_arcsec2 is None
    assert annulus.valid_area_pixel2==0



def test_entirely_masked_growth_sequence_is_not_claimed_measured():
    f=frame(value=5.)
    mask=np.ones(f.shape,dtype=bool)
    out=run_curve(f,source=mask,radii=(1.,2.,4.))
    assert out["science_status"]==(
        "UNDERCONSTRAINED_INSUFFICIENT_VALID_RADIAL_SAMPLES")
    assert out["n_measured_direct_apertures"]==0
    assert out["n_measured_independent_annuli"]==0
    assert out["max_abs_annulus_closure_nMgy"] is None
    assert out["adopted_for_science"] is False
    assert "NO_USABLE_RADIAL_CURVE_MEASUREMENT" in out["warning_flags"]
    assert all(r["signed_valid_flux_nMgy"] is None
               for r in out["cumulative"])
    json.dumps(out,allow_nan=False)


@pytest.mark.parametrize("invalid_weight", (-.2, 1.2, float("nan")))
def test_annulus_invalid_exact_pixel_weights_must_fail_closed(
        monkeypatch, invalid_weight):
    """Photutils output outside [0,1] must not silently alter signed flux."""
    import photutils.aperture as aperture
    original = aperture.EllipticalAnnulus
    class CorruptedAnnulus(original):
        def to_mask(self, method="exact"):
            result = super().to_mask(method=method)
            result.data[result.data.shape[0]//2,
                        result.data.shape[1]//2] = invalid_weight
            return result
    monkeypatch.setattr(aperture, "EllipticalAnnulus", CorruptedAnnulus)
    f = frame(value=.12)
    zeros = np.zeros(f.shape, dtype=bool)
    with pytest.raises(ValueError, match="Invalid exact annulus"):
        measure_native_annulus(
            f, geometry(a=6., q=.8, pa=33.), 2.0, provenance(),
            effective_source_mask=zeros, manual_exclude=zeros,
            sky=ApertureSkyScenario())


def test_offframe_growth_annulus_still_refuses_stale_or_nonboolean_scene_masks():
    """A no-overlap result must not skip input geometry/provenance checks."""
    f = frame()
    off = geometry(ra=155.45, dec=18., a=2., q=.8, pa=11.)
    good = np.zeros(f.shape, dtype=bool)
    empty = measure_native_annulus(
        f, off, 1., provenance(),
        effective_source_mask=good, manual_exclude=good,
        sky=ApertureSkyScenario())
    assert empty.no_valid_pixels is True
    assert empty.signed_flux_nMgy is None
    for bad in (np.zeros((2, 2), dtype=bool),
                np.zeros(f.shape, dtype=float)):
        with pytest.raises(ValueError, match="Annulus Scene masks"):
            measure_native_annulus(
                f, off, 1., provenance(),
                effective_source_mask=bad, manual_exclude=good,
                sky=ApertureSkyScenario())
    with pytest.raises(ValueError, match="Annulus Scene masks"):
        measure_native_annulus(
            f, off, 1., provenance(),
            effective_source_mask=good,
            manual_exclude=np.ones((1, 1), dtype=bool),
            sky=ApertureSkyScenario())


def test_growth_geometry_and_radii_reject_boolean_numeric_coercion():
    f=frame()
    for q, pa in ((True, 20.), (.7, False), ("0.7", 20.)):
        with pytest.raises(ValueError):
            GrowthHypothesis("invalid", q, pa, "SYNTHETIC")
    for inner in (True, "1.0", float("inf")):
        with pytest.raises(ValueError):
            measure_native_annulus(
                f, geometry(a=4.), inner, provenance(),
                effective_source_mask=np.zeros(f.shape, bool),
                manual_exclude=np.zeros(f.shape, bool),
                sky=ApertureSkyScenario())
    for radii in ((True,2.), ("1.0",2.), (1.,False)):
        with pytest.raises(ValueError):
            run_curve(f, radii=radii)
    accepted=GrowthHypothesis("numpy", np.float32(.8),
                              np.float64(42.), "SYNTHETIC")
    assert type(accepted.q) is float
    assert type(accepted.pa_deg) is float


def test_stage042_refuses_numerically_inconsistent_independent_annulus(
        monkeypatch):
    """A closure error is not merely descriptive when both operators exist."""
    from dataclasses import replace
    import lsbpipe.photometry.aperture_growth as module
    original = module.measure_native_annulus
    calls = []

    def corrupt_ring(*args, **kwargs):
        row = original(*args, **kwargs)
        calls.append(row)
        if row.signed_flux_nMgy is None:
            return row
        return replace(row, signed_flux_nMgy=row.signed_flux_nMgy+.01)
    monkeypatch.setattr(module, "measure_native_annulus", corrupt_ring)
    f=frame(value=.5)
    with pytest.raises(ValueError, match="independent annulus and cumulative flux closure"):
        run_curve(f, radii=(1., 2., 4.))
    assert calls and calls[0].signed_flux_nMgy is not None


def test_stage042_closure_gate_does_not_invent_flux_in_unobserved_ring():
    f=frame(value=.25)
    # If only the inner aperture survives, the larger masked ring
    # cannot acquire a synthetic zero-flux measurement to satisfy closure.
    yy,xx=np.indices(f.shape)
    source=np.hypot(xx-60, yy-60) > 5.
    result=run_curve(f, source=source, radii=(1., 2., 4., 6.),
                     q=1., pa=0.)
    assert result["all_comparable_annular_flux_closures_within_numerical_tolerance"]
    assert result["closure_is_not_a_noise_or_science_quality_test"]
    assert result["n_comparable_nested_radial_closures_verified"] >= 2
    assert all(
        row["signed_flux_nMgy"] is None or
        np.isfinite(row["signed_flux_nMgy"])
        for row in result["independently_integrated_annuli"])
    assert result["adopted_for_science"] is False


def test_growth_reuses_immutable_mask_snapshots_without_resnapshotting_each_radius(
        monkeypatch):
    """Stage04.2 must use exactly two owned byte-backed masks for all radii."""
    import lsbpipe.photometry.aperture_flux as aperture_flux
    original = aperture_flux._mask_for
    observed = []
    def inspect_mask(name, mask, shape):
        result = original(name, mask, shape)
        observed.append((name, mask, result))
        return result
    monkeypatch.setattr(aperture_flux, "_mask_for", inspect_mask)
    f = frame(value=.7, shape=(121, 121))
    source = np.zeros(f.shape, dtype=bool)
    manual = np.zeros(f.shape, dtype=bool)
    source[53:60, 58:64] = True
    manual[63:67, 65:69] = True
    radii = (1., 2., 3., 4., 5., 6., 7., 8.)
    result = run_curve(
        f, source=source, manual=manual, radii=radii)
    assert result["n_measured_direct_apertures"] == len(radii)
    assert len(observed) == 2*len(radii)
    for kind in ("effective source mask", "manual exclusion"):
        rows = [(inp, out) for key, inp, out in observed if key == kind]
        assert len(rows) == len(radii)
        assert len({id(inp) for inp, _ in rows}) == 1
        assert all(out is inp for inp, out in rows)
        assert all(not inp.flags.writeable for inp, _ in rows)
        for inp, _ in rows:
            with pytest.raises(ValueError):
                inp.setflags(write=True)
    assert all(
        row["signed_valid_flux_nMgy"] == pytest.approx(
            row["area_native_pixel2"]["valid_native_pixel2"]*.7, abs=1e-7)
        for row in result["cumulative"])


def test_byte_backed_mask_fast_path_rejects_mutable_and_readonly_views():
    """Readonly flags do not establish real immutability of a NumPy view."""
    from lsbpipe.photometry.aperture_flux import (
        _mask_for, _snapshot_native_boolean_mask)

    source = np.zeros((22, 23), bool)
    source[7, 9] = True
    verified = _snapshot_native_boolean_mask(
        "effective source mask", source, source.shape)
    assert _mask_for(
        "effective source mask", verified, source.shape) is verified
    with pytest.raises(ValueError):
        verified.setflags(write=True)

    mutable_view = source.view()
    mutable_view.flags.writeable = False
    detached = _mask_for(
        "effective source mask", mutable_view, source.shape)
    assert detached is not mutable_view
    source[7, 9] = False
    assert detached[7, 9] and verified[7, 9]
    source[:] = True
    assert not np.all(detached)
    assert not np.all(verified)

    with pytest.raises(ValueError, match="same-grid boolean"):
        _snapshot_native_boolean_mask(
            "manual exclusion", np.ones(source.shape, float), source.shape)



def test_growth_receipts_share_one_exact_immutable_native_grid_fingerprint(
        monkeypatch):
    """Twenty cumulative/ring computations cannot reserialize WCS metadata."""
    import lsbpipe.scene as scene

    original = scene.grid_fingerprint
    calls = []
    def counted(wcs, shape):
        calls.append(1)
        return original(wcs, shape)
    monkeypatch.setattr(scene, "grid_fingerprint", counted)
    f = frame(value=.12)
    output = run_curve(
        f, radii=(1., 2., 3., 4., 5., 6., 7., 8., 9., 10.))
    assert len(calls) == 1
    assert output["grid_fingerprint"] == original(f.wcs, f.shape)
    assert all(
        row["grid_fingerprint"] == output["grid_fingerprint"]
        for row in output["cumulative"])
    assert output["receipt_sha256"]
    assert not output["adopted_for_science"]
