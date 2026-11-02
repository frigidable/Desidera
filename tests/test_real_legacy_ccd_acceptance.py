"""Offline controls. Real archive acceptance runs separately with network."""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path
import numpy as np
from astropy.io import fits

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from accept_legacy_ccd_archives import archive_url, audit_file, run  # noqa: E402


def fixture_catalogue(path: Path) -> Path:
    data = (("decam", 574299, "N10", "r", 0, 40, 2, 25),
            ("decam", 574300, "N11", "g", 0, 40, 2, 25))
    names = ("CAMERA", "EXPNUM", "CCDNAME", "FILTER",
             "BRICK_X0", "BRICK_X1", "BRICK_Y0", "BRICK_Y1")
    formats = ("8A", "K", "8A", "1A", "J", "J", "J", "J")
    cols = [fits.Column(name=n, format=f, array=np.array([row[i] for row in data]))
            for i, (n, f) in enumerate(zip(names, formats))]
    fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU.from_columns(cols)]).writeto(path)
    return path


def test_official_archive_url():
    url = archive_url("dr9", "north", "2100p542")
    assert url.endswith("/dr9/north/coadd/210/2100p542/legacysurvey-2100p542-ccds.fits")
    assert url.startswith("https://portal.nersc.gov/")


def test_independent_astropy_agreement(tmp_path):
    path = fixture_catalogue(tmp_path / "legacysurvey-1554p180-ccds.fits")
    result, rows = audit_file(path, release="dr9", region="south", brick="1554p180",
                              source_url=archive_url("dr9", "south", "1554p180"))
    assert result["row_count"] == 2
    assert result["bands"]["r"]["count"] == 1
    assert result["bands"]["g"]["count"] == 1
    assert result["sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert len(rows) == 2
    assert result["release_and_region_binding"] == "URL_PATH_ASSERTION_NOT_INDEPENDENTLY_CERTIFIED"
    assert result["science_status"] == "NOT_SCIENTIFICALLY_ADOPTED"


def test_missing_real_cases_fail_closed(tmp_path):
    cache = tmp_path / "cache"
    file = cache / "dr9" / "south" / "1554p180" / "legacysurvey-1554p180-ccds.fits"
    file.parent.mkdir(parents=True)
    fixture_catalogue(file)
    out = tmp_path / "report"
    assert run(cache, out, skip_download=True) is False
    report = (out / "legacy_ccd_real_archive_report.json").read_text()
    assert '"status": "INCOMPLETE_OR_FAILED"' in report
    assert '"certifies_scientific_independence": false' in report
    assert (out / "legacy_ccd_real_archive_report.md").exists()
