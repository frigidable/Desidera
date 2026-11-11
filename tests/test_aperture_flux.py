"""Stage 04 finite signed exact-overlap aperture science-core, no adoption."""
import json
import math

import numpy as np
import pytest

pytest.importorskip("astropy")
pytest.importorskip("photutils")
from astropy.wcs import WCS

from lsbpipe.photometry.contracts import NativeBandFrame
from lsbpipe.photometry.aperture_flux import (
    ApertureGeometryICRS, ApertureProvenance,
    ApertureSkyScenario, measure_native_aperture,
    project_icrs_ellipse_to_native,
)


def wcs(*, ra=155.4, dec=18.0, scale=0.262, angle_deg=0.):
    w = WCS(naxis=2)
    w.wcs.crpix = [61.,61.]
    s = scale/3600.
    t = np.deg2rad(angle_deg)
    w.wcs.cd = s*np.array([
        [-np.cos(t),np.sin(t)],
        [np.sin(t),np.cos(t)]])
    w.wcs.crval = [ra,dec]
    w.wcs.ctype = ["RA---TAN","DEC--TAN"]
    return w


def frame(*, value=1., shape=(121,121), w=None,
          ivar=None, nexp=None, image=None):
    w = wcs() if w is None else w
    return NativeBandFrame(
        image_nMgy=(np.full(shape,value) if image is None else image),
        invvar_nMgy2=(np.full(shape,4.) if ivar is None else ivar),
        nexp=(np.full(shape,2.) if nexp is None else nexp),
        maskbits=np.zeros(shape,dtype=np.int16),
        native_wcs=w,release="dr11",region="south",
        brick="1554p180",band="r",
        source_evidence={"synthetic":True})


def provenance():
    return ApertureProvenance(
        project_id="leo-p-aperture-test",
        scene_draft_sha256="a"*64,
        catalogue_sha256="b"*64,
        source_policy_id="PROVISIONAL_SCENE_MASKS_V1")


def masks(f):
    return (np.zeros(f.shape,dtype=bool),
            np.zeros(f.shape,dtype=bool))


def geometry(*, ra=155.4,dec=18.0,a=7.5,q=.6,pa=25.):
    return ApertureGeometryICRS(
        ra_deg=ra,dec_deg=dec,semimajor_arcsec=a,
        axis_ratio=q,pa_deg_east_of_north=pa)


def measure(f, g=None, *, source=None, manual=None,sky=None):
    src,man=masks(f)
    return measure_native_aperture(
        f, geometry() if g is None else g,provenance(),
        effective_source_mask=src if source is None else source,
        manual_exclude=man if manual is None else manual,
        sky=sky)


def test_exact_fractional_ellipse_constant_image_flux_and_signed_calibration():
    f=frame(value=2.5)
    g=geometry()
    r=measure(f,g)
    area=r.area_native_pixel2
    assert r.schema=="desidera.photometry.native-aperture/1"
    assert r.dataset_id=="dr11/south/1554p180/r"
    assert r.science_state=="MEASURED_PARTIAL"
    assert area["inside_image_native_pixel2"]==pytest.approx(
        area["ideal_native_pixel2"],rel=1.e-5)
    assert area["valid_native_pixel2"]==pytest.approx(
        area["ideal_native_pixel2"],rel=1.e-5)
    assert r.signed_valid_flux_nMgy==pytest.approx(
        2.5*area["inside_image_native_pixel2"],abs=1.e-8)
    assert r.measured_mean_intensity_nMgy_per_arcsec2==pytest.approx(
        2.5/r.geometry_native["local_area_arcsec2_per_native_pixel"],
        rel=1.e-10)
    assert r.ab_magnitude_valid_area_only==pytest.approx(
        22.5-2.5*np.log10(r.signed_valid_flux_nMgy))
    assert r.formal_diagonal_sigma_nMgy > 0
    assert r.adopted_for_science is False
    assert r.covariance_calibrated is False
    assert r.psf_validated is False
    assert r.sky_science_adopted is False
    assert "COADD_PIXEL_COVARIANCE_UNKNOWN_FORMAL_IVAR_ONLY" in r.warnings
    assert "FULL_SOURCE_FITS_HASHES_NOT_VERIFIED" in r.warnings
    assert len(r.sector_qa)==8
    assert sum(s["valid_area_pixel2_by_centre_sector"] for s in r.sector_qa) == (
        pytest.approx(area["valid_native_pixel2"]))
    json.dumps(r.to_dict(),allow_nan=False)


def test_exact_fractional_weights_not_center_only():
    f=frame(value=1.)
    g=geometry(a=2.731,q=.73,pa=47.)
    r=measure(f,g)
    # A true ellipse has partial boundary pixels and a precise area.
    ideal=r.area_native_pixel2["ideal_native_pixel2"]
    assert r.area_native_pixel2["inside_image_native_pixel2"] == pytest.approx(
        ideal,rel=1e-4)
    assert abs(r.signed_valid_flux_nMgy-round(r.signed_valid_flux_nMgy)) > .01


def test_masks_accounting_no_double_count_and_no_auto_area_repair():
    f=frame(value=4.)
    src,man=masks(f)
    src[58:62,60:64]=True
    man[60:65,62:68]=True  # overlap must not count twice
    native=measure(f)
    r=measure(f,source=src,manual=man)
    a=r.area_by_reason_native_pixel2
    full=r.area_native_pixel2
    assert a["source_excluded_after_valid"]>0
    assert a["manual_excluded_after_source"]>0
    assert (full["valid_native_pixel2"]+a["frame_hard_or_invalid"]+
            a["source_excluded_after_valid"]+
            a["manual_excluded_after_source"])==pytest.approx(
                full["inside_image_native_pixel2"])
    assert r.signed_valid_flux_nMgy==pytest.approx(
        4.*full["valid_native_pixel2"],abs=1.e-8)
    assert r.signed_valid_flux_nMgy<native.signed_valid_flux_nMgy
    assert r.area_native_pixel2["valid_fraction_of_ideal"]<1
    assert "PARTIAL_VALID_APERTURE_NO_MASK_FLUX_RECONSTRUCTION" in r.warnings
    assert r.recipe_sha256 != native.recipe_sha256
    assert np.all(f.image_nMgy==4)
    assert not np.any(f.valid == False)


def test_hard_invalid_and_nexp_take_priority_over_empty_source_masks():
    im=np.ones((121,121))
    inv=np.full(im.shape,4.)
    exp=np.full(im.shape,2.)
    im[58,59]=np.nan
    inv[61,60]=0.
    exp[60,62]=0.
    f=frame(image=im,ivar=inv,nexp=exp)
    r=measure(f)
    assert r.area_by_reason_native_pixel2["frame_hard_or_invalid"]>=2
    assert np.isfinite(r.signed_valid_flux_nMgy)
    assert r.signed_valid_flux_nMgy==pytest.approx(
        r.area_native_pixel2["valid_native_pixel2"],abs=1.e-8)
    assert f.image_nMgy[58,59]!=f.image_nMgy[58,59]


def test_signed_zero_negative_flux_preserved_and_magnitude_undefined():
    f=frame(value=-.5)
    r=measure(f)
    assert r.signed_valid_flux_nMgy<0.
    assert r.ab_magnitude_valid_area_only is None
    assert r.measured_mean_intensity_nMgy_per_arcsec2<0
    zero=measure(frame(value=0.))
    assert zero.signed_valid_flux_nMgy==0.
    assert zero.ab_magnitude_valid_area_only is None


def test_exact_sky_constant_and_tilt_use_same_fractional_pixel_weights():
    f=frame(value=1.)
    original=measure(f)
    pedestal=measure(f,sky=ApertureSkyScenario(
        method="constant",b0_nMgy_per_pixel=.02,
        evidence_id="NONADOPTED"))
    area=original.area_native_pixel2["valid_native_pixel2"]
    assert pedestal.signed_valid_flux_nMgy==pytest.approx(
        original.signed_valid_flux_nMgy-.02*area,abs=1e-8)
    assert pedestal.sky_science_adopted is False
    assert "CONDITIONAL_NONZERO_SKY_SCENARIO_NOT_ADOPTED" in (
        pedestal.warnings)
    # Symmetric image centered on (x,y)=(60,60): first moment = 0.
    plane=measure(f,sky=ApertureSkyScenario(
        method="tilted_plane",
        b0_nMgy_per_pixel=.02,
        bx_nMgy_per_pixel2=.01,by_nMgy_per_pixel2=-.008,
        reference_xy_pix=(60.,60.),evidence_id="QA_ONLY"))
    assert plane.signed_valid_flux_nMgy==pytest.approx(
        pedestal.signed_valid_flux_nMgy,rel=1.e-7)
    assert original.recipe_sha256 != plane.recipe_sha256


def test_local_wcs_orientation_reflection_and_ra_wrap():
    g=geometry(ra=359.9999,dec=18.0,a=8.,q=.44,pa=0.)
    results=[]
    for angle in (0.,47.,-80.):
        w=wcs(ra=359.9999,angle_deg=angle)
        p=project_icrs_ellipse_to_native(g,w)
        assert p.a_pix>p.b_pix>0
        assert p.max_projection_error_pix<.15
        assert p.local_area_arcsec2_per_native_pixel==pytest.approx(
            .262**2,rel=1e-4)
        f=frame(w=w,value=1.)
        r=measure(f,g)
        assert r.area_native_pixel2["inside_image_native_pixel2"]==pytest.approx(
            r.area_native_pixel2["ideal_native_pixel2"],rel=1e-5)
        results.append(r.signed_valid_flux_nMgy)
    assert max(results)-min(results)<.05


def test_offframe_and_partially_truncated_apertures_have_honest_area():
    f=frame(value=1.)
    w=f.wcs
    ra,dec=w.all_pix2world(8.,60.,0)
    g=geometry(ra=float(ra),dec=float(dec),a=7.,q=1.)
    partial=measure(f,g)
    assert partial.area_native_pixel2["outside_image_native_pixel2"]>0
    assert partial.area_native_pixel2["valid_fraction_of_ideal"]<1
    assert "APERTURE_PARTLY_OUTSIDE_NATIVE_FITS" in partial.warnings
    assert partial.signed_valid_flux_nMgy==pytest.approx(
        partial.area_native_pixel2["valid_native_pixel2"],abs=1.e-8)
    ra,dec=w.all_pix2world(310.,330.,0)
    no_overlap=measure(f,geometry(
        ra=float(ra),dec=float(dec),a=3.,q=1.))
    assert no_overlap.signed_valid_flux_nMgy is None
    assert no_overlap.ab_magnitude_valid_area_only is None
    assert no_overlap.formal_diagonal_sigma_nMgy is None
    assert no_overlap.science_state=="UNDERCONSTRAINED_NO_VALID_PIXELS"
    assert no_overlap.area_native_pixel2["valid_native_pixel2"]==0.
    json.dumps(no_overlap.to_dict(),allow_nan=False)


def test_schema_rejects_bad_sky_geometry_and_mask_shapes():
    for kwargs in (
        {"semimajor_arcsec":0.},
        {"axis_ratio":0.},
        {"axis_ratio":1.01},
        {"ra_deg":360.},
        {"pa_deg_east_of_north":-5.},
        {"dec_deg":91.},
    ):
        args=dict(ra_deg=155.,dec_deg=18.,
                  semimajor_arcsec=8.,axis_ratio=.5,
                  pa_deg_east_of_north=20.)
        args.update(kwargs)
        with pytest.raises(ValueError):
            ApertureGeometryICRS(**args)
    with pytest.raises(ValueError):
        ApertureSkyScenario(
            method="zero_additional_correction",
            b0_nMgy_per_pixel=.1)
    with pytest.raises(ValueError):
        ApertureSkyScenario(
            method="constant",bx_nMgy_per_pixel2=.1)
    with pytest.raises(ValueError):
        ApertureProvenance(
            project_id="p",scene_draft_sha256="oops",
            catalogue_sha256="b"*64,source_policy_id="mask")
    f=frame()
    with pytest.raises(ValueError,match="same-grid boolean"):
        measure_native_aperture(
            f,geometry(),provenance(),
            effective_source_mask=np.zeros((3,3),bool),
            manual_exclude=np.zeros(f.shape,bool))
    with pytest.raises(ValueError,match="same-grid boolean"):
        measure_native_aperture(
            f,geometry(),provenance(),
            effective_source_mask=np.zeros(f.shape,int),
            manual_exclude=np.zeros(f.shape,bool))
