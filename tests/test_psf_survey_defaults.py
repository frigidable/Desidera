"""Legacy native PSFSIZE is an OPTIONAL per-pixel seeing prior, not ePSF."""
from pathlib import Path

import numpy as np
import pytest
pytest.importorskip("astropy")
from astropy.io import fits
from astropy.wcs import WCS

from lsbpipe.photometry.contracts import NativeBandFrame
from lsbpipe.photometry.psf_survey_defaults import (
    suggest_detection_fwhm_from_psfsize,
)


def _make_frame(tmp_path, *, write_map=True, moved_wcs=False):
    shape = (61, 61)
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [31,31]
    wcs.wcs.crval = [155.438012,18.088025]
    wcs.wcs.cdelt = [-.262/3600,.262/3600]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    name = tmp_path/"legacysurvey-1554p180-image-g.fits.fz"
    name.write_bytes(b"Original image path fixture only")
    if write_map:
        header = wcs.to_header()
        if moved_wcs:
            header["CRVAL1"] += .002
        # FITS reader recognizes the content independently of extension.
        fits.PrimaryHDU(
            data=np.full(shape, 1.20, dtype=np.float32),
            header=header).writeto(
                tmp_path/"legacysurvey-1554p180-psfsize-g.fits.fz")
    frame = NativeBandFrame(
        image_nMgy=np.zeros(shape, dtype=float),
        invvar_nMgy2=np.ones(shape, dtype=float),
        nexp=np.ones(shape, dtype=np.int16),
        maskbits=np.zeros(shape, dtype=np.int16),
        native_wcs=wcs,
        release="dr11", region="south", brick="1554p180",
        band="g", source_evidence={
            "synthetic":True, "files":{"image":str(name)},
        },
    )
    return frame


def test_psfsize_map_proposes_native_fwhm_from_real_wcs(tmp_path):
    frame = _make_frame(tmp_path)
    suggestion = suggest_detection_fwhm_from_psfsize(
        frame, target_icrs_deg=(155.438012,18.088025))
    assert suggestion["available"] is True
    assert suggestion["band"] == "g"
    assert suggestion["psfsize_arcsec"] == pytest.approx(1.2,rel=1e-6)
    assert suggestion["image_pixel_scale_arcsec"] == pytest.approx(.262,rel=1e-5)
    assert suggestion["suggested_detection_fwhm_pix"] == pytest.approx(
        1.2/.262,rel=1e-5)
    assert len(suggestion["source_sha256_local"]) == 64
    assert suggestion["is_measured_empirical_psf"] is False
    assert suggestion["requires_researcher_review"] is True


def test_no_psfsize_map_never_invents_seeing(tmp_path):
    frame = _make_frame(tmp_path,write_map=False)
    suggestion = suggest_detection_fwhm_from_psfsize(
        frame, target_icrs_deg=(155.438012,18.088025))
    assert suggestion["available"] is False
    assert "ABSENT" in suggestion["source"]
    assert "suggested_detection_fwhm_pix" not in suggestion


def test_wrong_psfsize_wcs_must_fail_closed(tmp_path):
    frame = _make_frame(tmp_path,moved_wcs=True)
    with pytest.raises(ValueError, match="WCS"):
        suggest_detection_fwhm_from_psfsize(
            frame, target_icrs_deg=(155.438012,18.088025))
