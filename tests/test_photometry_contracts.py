"""P0.5 acceptance contract on arrays and local original DESI Legacy FITS."""
from __future__ import annotations

import hashlib
import json
import numpy as np
import pytest

pytest.importorskip("astropy")
from astropy.io import fits
from astropy.wcs import WCS

from lsbpipe.photometry import (
    MaskPolicy, NativeBandFrame, NativeFrameCollection,
    PSFDescription, PhotometryPreparation, SkyDescription,
    load_native_legacy_quartet, CONTRACT_SCHEMA,
)
from lsbpipe.survey import LEGACY_BITS
from lsbpipe.surveys.legacy import BrickHit, product_for
from lsbpipe.surveys.task_control import OperationCancelled


def wcs(*, ra=40.096208, dec=1.22925):
    obj = WCS(naxis=2)
    obj.wcs.crpix = [6.5, 6.5]
    obj.wcs.cdelt = [-0.262/3600, 0.262/3600]
    obj.wcs.crval = [ra, dec]
    obj.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    return obj


def frame(*, region="south", band="r", brick="0401p012",
          policy=None, extra=None, bits=None, ivar=None, nexp=None):
    shape = (12, 12)
    return NativeBandFrame(
        image_nMgy=np.ones(shape, dtype=np.float32),
        invvar_nMgy2=(np.ones(shape, dtype=np.float32) * 4
                      if ivar is None else ivar),
        nexp=(np.ones(shape, dtype=np.int16) * 2
              if nexp is None else nexp),
        maskbits=(np.zeros(shape, dtype=np.int16) if bits is None else bits),
        native_wcs=wcs(), release="dr9", region=region,
        brick=brick, band=band,
        policy=policy if policy is not None else MaskPolicy(),
        extra_exclusion=extra,
        source_evidence={"synthetic": True, "source": ["array fixture"]})


def make_fits(root, *, band="r", region="south", brick="0401p012",
              image_bunit="nanomaggies", shifted=None, missing=None):
    hit = BrickHit("dr9", region, brick, 40, 40.2, 1, 1.4, {band: 2})
    originals = {}
    for kind in ("image", "invvar", "nexp", "maskbits"):
        if kind == missing:
            continue
        prod = product_for(hit, kind, None if kind == "maskbits" else band)
        dest = root / prod.relative_path
        dest.parent.mkdir(parents=True, exist_ok=True)
        hdr = wcs().to_header()
        if kind == shifted:
            hdr["CRVAL1"] += 0.001
        if kind == "image" and image_bunit is not None:
            hdr["BUNIT"] = image_bunit
        if kind == "invvar":
            hdr["BUNIT"] = "1/nanomaggies^2"
        ar = (np.zeros((12, 12), dtype=np.int16)
              if kind == "maskbits"
              else np.ones((12, 12),
                           dtype=np.int16 if kind == "nexp" else np.float32))
        if kind == "maskbits":
            # Optical mask is HDU1, WISE masks in later extensions.
            fits.HDUList([
                fits.PrimaryHDU(),
                fits.ImageHDU(data=ar, header=hdr),
                fits.ImageHDU(data=np.ones_like(ar) * 100),
            ]).writeto(dest)
        else:
            fits.PrimaryHDU(data=ar, header=hdr).writeto(dest)
        originals[kind] = dest
    return originals


def test_immutable_native_frame_is_not_writable_even_via_setflags():
    image = np.ones((12, 12))
    f = frame()
    assert f.image_nMgy.flags.writeable is False
    assert f.valid.flags.writeable is False
    for arr in (f.image_nMgy, f.invvar_nMgy2, f.maskbits,
                f.nexp, f.mask_reasons["BRIGHT"]):
        with pytest.raises(ValueError):
            arr.setflags(write=True)
    with pytest.raises(TypeError):
        f.source_evidence["synthetic"] = False
    image[:] = 12
    assert np.all(f.image_nMgy == 1)
    independent = f.wcs
    independent.wcs.crval[0] += 0.3
    assert f.wcs.wcs.crval[0] == pytest.approx(40.096208)


def test_pixel_area_local_not_equivalent_to_exact_solid_angle():
    f = frame()
    assert f.local_pixel_area_arcsec2 == pytest.approx(0.262**2)
    record = f.summary()
    assert record["schema"] == CONTRACT_SCHEMA
    assert record["wcs"]["area_policy"] == "LOCAL_WCS_PROJECTION_APPROXIMATE"
    assert record["wcs"]["reprojected"] is False
    assert record["covariance"] == "DIAGONAL_IVAR_CORRELATIONS_UNKNOWN"
    assert record["science_state"] == "NATIVE_VALIDATED_NOT_PHOTOMETRY_READY"
    json.dumps(record, allow_nan=False)


def test_mask_policy_distinguishes_default_hard_and_diagnostic_flags():
    mask = np.zeros((12, 12), dtype=np.int16)
    mask[0, 0] = 1 << LEGACY_BITS["SATUR_R"]
    mask[0, 1] = 1 << LEGACY_BITS["ALLMASK_R"]
    mask[0, 2] = 1 << LEGACY_BITS["BRIGHT"]
    mask[0, 3] = 1 << LEGACY_BITS["GALAXY"]
    mask[0, 4] = 1 << LEGACY_BITS["NPRIMARY"]
    f = frame(bits=mask)
    assert not f.valid[0, 0] and not f.valid[0, 1]
    assert f.valid[0, 2] and f.valid[0, 3] and f.valid[0, 4]
    assert f.mask_reasons["SATUR_R"][0, 0]
    assert f.mask_reasons["GALAXY"][0, 3]
    aggressive = frame(
        bits=mask, policy=MaskPolicy(
            reject_nprimary=True, reject_bright=True, reject_galaxy=True))
    assert not any(aggressive.valid[0, :5])
    permissive = frame(
        bits=mask, policy=MaskPolicy(
            reject_satur=False, reject_allmask=False))
    assert all(permissive.valid[0, :5])


def test_ivar_and_nexp_rejected_independently_without_resetting_flux():
    ivar = np.full((12, 12), 4.)
    nexp = np.full((12, 12), 2.)
    ivar[0, 0] = 0
    ivar[0, 1] = np.nan
    nexp[0, 2] = 0
    nexp[0, 3] = np.nan
    f = frame(ivar=ivar, nexp=nexp)
    assert not np.any(f.valid[0, :4])
    assert np.all(f.valid[0, 4:])
    assert f.mask_reasons["BAD_NEXP"][0, 3]
    assert f.image_nMgy[0, 0] == 1  # original data NOT zero-filled


def test_user_mask_is_separate_and_has_no_silent_shape_coercion():
    extra = np.zeros((12, 12), dtype=bool)
    extra[4, 6] = True
    f = frame(extra=extra)
    assert not f.valid[4, 6]
    assert f.mask_reasons["USER_EXCLUSION"][4, 6]
    with pytest.raises(ValueError, match="shape"):
        frame(extra=np.zeros((12, 11), dtype=bool))


def test_preparation_requires_explicit_sky_psf_and_correlations_remain_unknown():
    f = frame()
    unready = PhotometryPreparation(frame=f)
    assert "RESIDUAL_SKY_NOT_ESTIMATED" in unready.pending
    assert "PSF_NOT_MEASURED" in unready.pending
    with pytest.raises(ValueError, match="non-estimated"):
        SkyDescription(background_nMgy_per_pixel=np.zeros(f.shape))
    with pytest.raises(ValueError, match="kernel"):
        PSFDescription(state="EMPIRICAL_KERNEL")
    with pytest.raises(ValueError, match="normalized"):
        PSFDescription(
            state="EMPIRICAL_KERNEL", kernel=np.ones((3, 3)))
    psf = PSFDescription(
        state="EMPIRICAL_KERNEL",
        kernel=np.ones((3, 3)) / 9,
        provenance={"star_count": 24})
    sky = SkyDescription(
        state="ESTIMATED", background_nMgy_per_pixel=np.zeros(f.shape),
        provenance={"algorithm": "synthetic constant"})
    prepared = PhotometryPreparation(frame=f, sky=sky, psf=psf)
    assert "RESIDUAL_SKY_NOT_ESTIMATED" not in prepared.pending
    assert "PSF_NOT_MEASURED" not in prepared.pending
    assert "COADD_NOISE_CORRELATION_NOT_CALIBRATED" in prepared.pending
    assert "PHOTOMETRIC_ALGORITHM_NOT_IMPLEMENTED" in prepared.pending
    with pytest.raises(ValueError, match="match"):
        PhotometryPreparation(
            frame=f, sky=SkyDescription(
                state="ESTIMATED", background_nMgy_per_pixel=np.zeros((2, 2))))
    with pytest.raises(ValueError, match="odd-sized"):
        PSFDescription(state="EMPIRICAL_KERNEL",
                       kernel=np.ones((4, 4)) / 16)


def test_multiple_bricks_are_separate_and_north_south_not_summed():
    frames = (
        frame(brick="0401p012", region="south"),
        frame(brick="0402p012", region="south"),
        frame(brick="0401p012", region="north"),
    )
    coll = NativeFrameCollection(frames)
    assert len(coll.processing_groups[("dr9", "south", "r")]) == 2
    assert len(coll.processing_groups[("dr9", "north", "r")]) == 1
    assert coll.mosaic_state == "NATIVE_UNMOSAICKED_NO_REPROJECTION"
    with pytest.raises(ValueError, match="explicit mosaic"):
        coll.require_single_native()
    with pytest.raises(ValueError, match="Duplicate"):
        NativeFrameCollection((frames[0], frames[0]))
    assert NativeFrameCollection((frames[0],)).require_single_native() is frames[0]
    with pytest.raises(ValueError, match="unavailable"):
        frame(region="north", band="i")


def test_loader_reads_optical_mask_hdu1_and_never_modifies_fits(tmp_path):
    paths = make_fits(tmp_path)
    before = {p: hashlib.sha256(p.read_bytes()).hexdigest()
              for p in paths.values()}
    f = load_native_legacy_quartet(
        tmp_path, release="dr9", region="south",
        brick="0401p012", band="r")
    assert f.shape == (12, 12)
    assert f.dataset_id == "dr9/south/0401p012/r"
    assert f.source_evidence["hdu_indices"]["maskbits"] == 1
    assert f.source_evidence["image_unit_evidence"] == "IMAGE_FITS_BUNIT_CONFIRMED"
    assert len(f.source_evidence["local_hashes"]) == 4
    assert np.all(f.maskbits == 0), "optical HDU1, not the WISE HDU2"
    assert {p: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in paths.values()} == before


def test_loader_missing_bunit_requires_explicit_provenance(tmp_path):
    make_fits(tmp_path, image_bunit=None)
    kwargs = dict(
        release="dr9", region="south", brick="0401p012", band="r")
    with pytest.raises(ValueError, match="BUNIT missing"):
        load_native_legacy_quartet(tmp_path, **kwargs)
    accepted = load_native_legacy_quartet(
        tmp_path, **kwargs, require_explicit_bunit=False)
    assert accepted.source_evidence["image_unit_evidence"] == (
        "LEGACY_SURVEY_DOCUMENTATION_ASSUMPTION")


def test_loader_rejects_wrong_units_wcs_and_missing_fits(tmp_path):
    kwargs = dict(release="dr9", region="south",
                  brick="0401p012", band="r")
    make_fits(tmp_path, image_bunit="electrons")
    with pytest.raises(ValueError, match="BUNIT"):
        load_native_legacy_quartet(tmp_path, **kwargs)
    make_fits(tmp_path / "wrong-wcs", shifted="invvar")
    with pytest.raises(ValueError, match="WCS"):
        load_native_legacy_quartet(tmp_path / "wrong-wcs", **kwargs)
    make_fits(tmp_path / "missing", missing="nexp")
    with pytest.raises(FileNotFoundError, match="nexp"):
        load_native_legacy_quartet(tmp_path / "missing", **kwargs)


def test_loader_cooperative_cancel_preserves_original_files(tmp_path):
    paths = make_fits(tmp_path)
    before = {p: p.read_bytes() for p in paths.values()}
    count = [0]
    def stop():
        count[0] += 1
        if count[0] >= 4:
            raise OperationCancelled("photometry contract load cancelled")
    with pytest.raises(OperationCancelled):
        load_native_legacy_quartet(
            tmp_path, release="dr9", region="south",
            brick="0401p012", band="r", cancel_check=stop)
    assert all(p.read_bytes() == raw for p, raw in before.items())


def test_loader_rejects_path_traversal_before_touching_files(tmp_path):
    with pytest.raises(ValueError, match="brick"):
        load_native_legacy_quartet(
            tmp_path, release="dr9", region="south",
            brick="../../..", band="r")
    with pytest.raises(ValueError, match="brick"):
        load_native_legacy_quartet(
            tmp_path, release="dr9", region="south",
            brick="0401p0/2", band="r")


def test_signed_optical_maskbits_are_interpreted_as_uint16_bits():
    # Optical bit 15 (ALLMASK_I) appears negative when read as signed int16.
    bits = np.zeros((12, 12), dtype=np.int16)
    bits[0, 0] = np.array([32768], dtype=np.uint16).view(np.int16)[0]
    f = frame(bits=bits)
    assert f.maskbits[0, 0] == 32768
    assert not f.mask_reasons["ALLMASK_R"][0, 0]


def test_input_wcs_can_change_after_construction_without_mutating_frame():
    original = wcs()
    arr = np.ones((12, 12))
    f = NativeBandFrame(
        image_nMgy=arr, invvar_nMgy2=arr, nexp=arr,
        maskbits=np.zeros((12, 12), dtype=np.int16),
        native_wcs=original, brick="0401p012")
    original.wcs.crval[0] = 90.0
    assert f.wcs.wcs.crval[0] == pytest.approx(40.096208)
