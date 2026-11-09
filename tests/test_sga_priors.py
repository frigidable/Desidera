"""Reproducible SGA-2020/2025 local geometry-prior reader tests."""
from __future__ import annotations

import numpy as np
import pytest
from astropy.io import fits

from lsbpipe.photometry.sga_priors import read_local_sga_priors


def _write_2020(path, rows=None):
    if rows is None:
        rows=[(71,"Leo P",155.438012,18.088025,.62,"SB26",.70,122.,0),
              (72,"Other",155.438020,18.088030,.40,"LEDA",.9,12.,16),
              (74,"Distant",170.0,25.0,.5,"SB26",.8,50.,0)]
    cols = [
        fits.Column(name="SGA_ID",format="K",array=[r[0] for r in rows]),
        fits.Column(name="GALAXY",format="24A",array=[r[1] for r in rows]),
        fits.Column(name="RA",format="D",array=[r[2] for r in rows]),
        fits.Column(name="DEC",format="D",array=[r[3] for r in rows]),
        fits.Column(name="D26",format="E",array=[r[4] for r in rows]),
        fits.Column(name="D26_REF",format="4A",array=[r[5] for r in rows]),
        fits.Column(name="BA",format="E",array=[r[6] for r in rows]),
        fits.Column(name="PA",format="E",array=[r[7] for r in rows]),
        fits.Column(name="ELLIPSEBIT",format="J",array=[r[8] for r in rows]),
    ]
    fits.HDUList([fits.PrimaryHDU(),
                  fits.BinTableHDU.from_columns(cols,name="ELLIPSE")]).writeto(path)
    return path


def test_sga_2020_multiple_nearby_matches_require_explicit_selection(tmp_path):
    path=_write_2020(tmp_path/"sga2020.fits")
    result=read_local_sga_priors(path,(155.438012,18.088025))
    assert result.release=="SGA2020"
    assert result.status=="SGA_LOCAL_MATCHES_REQUIRE_RESEARCHER_SELECTION"
    assert len(result.matches)==2
    assert result.matches[0].sga_id==71
    assert result.matches[1].sga_id==72
    assert "D26_NOT_DIRECT_SB26" in result.matches[1].flags
    assert "SGA2020_ELLIPSE_NOT_FITTED" in result.matches[1].flags
    shape=result.matches[0].core_shape()
    assert shape.kind=="ellipse"
    assert shape.a_arcsec==pytest.approx(.62*30.,rel=1e-5)
    assert shape.q==pytest.approx(.7)
    assert shape.pa_deg==pytest.approx(122.)
    assert result.matches[0].source_sha256==result.source_sha256
    assert len(result.source_sha256)==64
    assert result.matches[0].status=="SGA_CORE_PRIOR_NOT_SCIENCE_ADOPTED"


def test_sga_2025_draft_quality_flags_and_sentinel_rows(tmp_path):
    path=tmp_path/"sga2025.fits"
    cols=[
        fits.Column(name="SGAID",format="K",array=[17,18]),
        fits.Column(name="GALAXY",format="24A",array=["Leo P","Bad shape"]),
        fits.Column(name="RA",format="D",array=[155.438,155.440]),
        fits.Column(name="DEC",format="D",array=[18.088,18.089]),
        fits.Column(name="D26",format="E",array=[.40,-1.]),
        fits.Column(name="D26_REF",format="3A",array=["mom","fix"]),
        fits.Column(name="BA",format="E",array=[.6,-1.]),
        fits.Column(name="PA",format="E",array=[55.,-1.]),
        fits.Column(name="ELLIPSEBIT",format="J",array=[16384+16+8,4]),
        fits.Column(name="ELLIPSEMODE",format="J",array=[1,1]),
    ]
    fits.HDUList([fits.PrimaryHDU(),fits.BinTableHDU.from_columns(
        cols,name="SGA2025")]).writeto(path)
    result=read_local_sga_priors(path,(155.438,18.088))
    assert len(result.matches)==2
    a,b=result.matches
    assert a.release=="SGA2025"
    assert "SGA2025_DRAFT_NOT_FINAL_RELEASE" in a.flags
    assert "SGA_FIT_GEOMETRY_FAILED" in a.flags
    assert "SGA_LARGE_CENTRE_SHIFT" in a.flags
    assert "SGA_BLENDED" in a.flags
    assert "SGA_FIXED_GEOMETRY" in a.flags
    assert b.core_shape() is None
    assert "NO_VALID_D26" in b.flags


def test_sga_empty_wrong_hdu_symlink_and_wrap_match(tmp_path):
    path=_write_2020(tmp_path/"sgawrap.fits",[
        (91,"RA wrap",359.9999,60.,.6,"SB26",.75,60.,0)])
    lookup=read_local_sga_priors(path,(.0001,60.))
    assert len(lookup.matches)==1
    assert lookup.matches[0].offset_arcsec<1.
    empty=read_local_sga_priors(path,(155.,18.))
    assert empty.matches==()
    invalid=tmp_path/"invalid.fits"
    fits.HDUList([fits.PrimaryHDU(),
                  fits.BinTableHDU.from_columns([
                      fits.Column(name="RA",format="D",array=[1.])],
                      name="NOT_SGA")]).writeto(invalid)
    with pytest.raises(ValueError,match="Neither SGA2025"):
        read_local_sga_priors(invalid,(155.,18.))
    shortcut=tmp_path/"link.fits"
    try:
        shortcut.symlink_to(path)
    except (OSError,NotImplementedError):
        pass
    else:
        with pytest.raises(FileNotFoundError):
            read_local_sga_priors(shortcut,(155.,18.))


def test_sga_malformed_geometry_does_not_become_protection(tmp_path):
    path=_write_2020(tmp_path/"badrow.fits",[
        (11,"uncertain",155.,18.,-1.,"LEDA",.2,50.,0),
        (12,"invalid q",155.01,18.,.5,"SB26",-1.,40.,0)])
    result=read_local_sga_priors(path,(155.,18.))
    assert all(match.core_shape() is None for match in result.matches)
