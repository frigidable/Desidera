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
