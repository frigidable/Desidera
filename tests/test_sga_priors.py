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


def test_official_viewer_sga2025_small_json_lookup_and_not_d26(monkeypatch):
    """Viewer geometry != merged scientific SGA2025 catalogue."""
    import json
    import io
    from urllib.parse import urlsplit, parse_qs
    import lsbpipe.photometry.sga_priors as module

    calls=[]
    reply=json.dumps({
        "rd":[[155.438,18.088],[155.441,18.089],[155.8,18.8]],
        "name":["SGA-731","SGA-732","SGA-900"],
        "radiusArcsec":[19.0,12.0,13.0],
        "abRatio":[.65,.9,.3],"posAngle":[70.,110.,10.],
        "pgc":[731,732,900],
        "groupname":["alpha","beta","far"],
        "type":["","",""],"redshift":[-1,-1,-1],
        "color":["#377eb8"]*3,
    }).encode()
    class Response(io.BytesIO):
        def __init__(self, url):
            super().__init__(reply)
            self.url=url
        def geturl(self):
            return self.url

    def fake_open(url,timeout):
        calls.append((url,timeout))
        return Response(url)
    monkeypatch.setattr(module,"urlopen",fake_open)
    result=module.fetch_viewer_sga2025_priors(
        (155.438,18.088),region="south",radius_arcmin=3.)
    assert len(calls)==1
    assert calls[0][0].startswith(
        "https://www.legacysurvey.org/viewer/sga2025-south/1/cat.json?")
    query=parse_qs(urlsplit(calls[0][0]).query)
    assert set(query)=={"ralo","rahi","declo","dechi"}
    assert len(result.matches)==2
    first=result.matches[0]
    assert first.sga_id==731
    assert first.source_path==calls[0][0]
    assert first.source_sha256==result.source_sha256
    assert first.d26_ref=="VIEWER_RADIUS_ARCSEC"
    assert "VIEWER_REF_ID_NOT_VERIFIED_SGAID" in first.flags
    assert "NO_SGA2025_D26_REF_OR_FIT_QUALITY" in first.flags
    assert first.core_shape().a_arcsec==pytest.approx(19.)
    assert first.ba==pytest.approx(.65)
    assert first.pa_deg==pytest.approx(70.)


def test_official_viewer_no_entries_and_ra_wrap(monkeypatch):
    import json
    import io
    from urllib.parse import urlsplit,parse_qs
    import lsbpipe.photometry.sga_priors as module

    empty=json.dumps({
       "rd":[],"name":[],"radiusArcsec":[],"abRatio":[],
       "posAngle":[],"pgc":[],"type":[],"redshift":[],
    }).encode()
    class Response(io.BytesIO):
        def __init__(self,payload,url):
            super().__init__(payload)
            self.url=url
        def geturl(self):return self.url

    captured=[]
    def fake_open(url,timeout):
        captured.append(url)
        return Response(empty,url)
    monkeypatch.setattr(module,"urlopen",fake_open)
    result=module.fetch_viewer_sga2025_priors(
        (359.999,60.),region="north")
    assert result.matches==()
    assert result.release=="SGA2025"
    assert "/sga2025-north/" in result.source_path
    q=parse_qs(urlsplit(captured[0]).query)
    assert float(q["ralo"][0]) > float(q["rahi"][0])


def test_official_viewer_fails_closed_on_redirect_or_malformed_catalog(monkeypatch):
    import io
    import json
    import lsbpipe.photometry.sga_priors as module
    class Response(io.BytesIO):
        def __init__(self,data,url):
            super().__init__(data)
            self.url=url
        def geturl(self):return self.url

    monkeypatch.setattr(module,"urlopen",
        lambda url,timeout:Response(b'{}',url))
    with pytest.raises(ValueError,match="schema mismatch"):
        module.fetch_viewer_sga2025_priors((155.,18.))
    monkeypatch.setattr(module,"urlopen",
        lambda url,timeout:Response(b'{}',"https://example.org/evil"))
    with pytest.raises(ValueError,match="Untrusted"):
        module.fetch_viewer_sga2025_priors((155.,18.))
    data=json.dumps({
        "rd":[[155.,18.]],"name":[],"radiusArcsec":[20.],
        "abRatio":[.8],"posAngle":[20.],"pgc":[123],
    }).encode()
    monkeypatch.setattr(module,"urlopen",
        lambda url,timeout:Response(data,url))
    with pytest.raises(ValueError,match="inconsistent"):
        module.fetch_viewer_sga2025_priors((155.,18.))
    with pytest.raises(ValueError,match="safe limits"):
        module.fetch_viewer_sga2025_priors((155.,89.))
