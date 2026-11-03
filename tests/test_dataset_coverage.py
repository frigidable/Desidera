"""P0.5b geometry: actual header WCS, never guessed brick boundaries."""
from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("astropy")
from astropy.io import fits
from astropy.wcs import WCS

from lsbpipe.surveys.object_projects import create_object_project
from lsbpipe.surveys.dataset_registry import (
    register_external_fits, import_legacy_fits, read_dataset_index,
)
from lsbpipe.surveys.dataset_coverage import survey_coverage_footprints
from lsbpipe.surveys.legacy import SkyField, BrickHit, product_for


def _original_wcs(path, *, dx=0., valid=True):
    w = WCS(naxis=2)
    w.wcs.crpix = [25., 25.]
    w.wcs.cdelt = [-.262/3600, .262/3600]
    w.wcs.crval = [40.096208 + dx/3600., 1.22925]
    w.wcs.ctype = (["RA---TAN", "DEC--TAN"] if valid
                   else ["LINEAR", "LINEAR"])
    path.parent.mkdir(parents=True, exist_ok=True)
    fits.PrimaryHDU(
        data=np.ones((50, 50), dtype=np.float32),
        header=w.to_header()).writeto(path)


@pytest.fixture
def project(tmp_path, monkeypatch):
    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path / "session"))
    return create_object_project(
        "UGC 2162", tmp_path / "UGC_2162",
        field=SkyField(40.096208, 1.22925, 5.))


def test_header_only_geometry_from_actual_fits_not_brick_name(project, tmp_path):
    p1 = tmp_path / "HSC_g.fits"
    p2 = tmp_path / "CFHT_r.fits"
    _original_wcs(p1, dx=0.)
    _original_wcs(p2, dx=20.)
    register_external_fits(project.root, p1, provider="HSC", band="g")
    register_external_fits(project.root, p2, provider="CFHT", band="r")
    before = (p1.read_bytes(), p2.read_bytes())
    coverage = survey_coverage_footprints(project.root)
    assert coverage["status"] == "GEOMETRY_ONLY_NOT_ADOPTED"
    assert len(coverage["footprints"]) == 2
    assert not coverage["rejected"]
    a, b = coverage["footprints"]
    assert len(a["points"]) == len(b["points"]) == 36
    assert a["geometry_status"] == "FITS_WCS_HEADER_ONLY"
    assert a["native_shape"] == [50, 50]
    assert (np.mean([x for x, y in b["points"]]) <
            np.mean([x for x, y in a["points"]]))  # RA increases left
    assert (p1.read_bytes(), p2.read_bytes()) == before


def test_compressed_fits_extension_has_real_header_wcs_geometry(
        project, tmp_path):
    """Legacy survey files are commonly compressed FITS, not plain .fits."""
    w = WCS(naxis=2)
    w.wcs.crpix = [35., 40.]
    w.wcs.cdelt = [-.262/3600., .262/3600.]
    w.wcs.crval = [40.096208, 1.22925]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    path = tmp_path / "compressed_image.fits.fz"
    image = np.ones((80, 70), dtype=np.int16)
    fits.HDUList([
        fits.PrimaryHDU(),
        fits.CompImageHDU(data=image, header=w.to_header(),
                          compression_type="RICE_1"),
    ]).writeto(path)
    register_external_fits(project.root, path, provider="Test DECam", band="r")
    response = survey_coverage_footprints(project.root)
    assert len(response["footprints"]) == 1
    assert not response["rejected"]
    assert response["footprints"][0]["native_shape"] == [80, 70]


def test_invalid_fits_and_no_wcs_are_reported_not_mapped(project, tmp_path):
    bad = tmp_path / "corrupt.fits"
    bad.write_bytes(b"NOT_A_FITS")
    unregistered = tmp_path / "unregistered.fits"
    _original_wcs(unregistered)
    register_external_fits(project.root, bad, provider="Unknown", band="r")
    report = survey_coverage_footprints(project.root)
    assert report["footprints"] == []
    assert len(report["rejected"]) == 1
    assert "Error" in report["rejected"][0]["reason"]
    # Never scan arbitrary nearby FITS or fabricate a brick footprint.
    assert unregistered.is_file()
    assert not report["omitted_count"]


def test_header_coverage_cap_reports_omissions_not_silent_skips(
        project, tmp_path):
    for n in range(4):
        f = tmp_path / f"image_{n}.fits"
        _original_wcs(f, dx=n*10.)
        register_external_fits(project.root, f, provider="HSC", band="g")
    receipt = survey_coverage_footprints(project.root, max_drawn=2)
    assert len(receipt["footprints"]) == 2
    assert receipt["omitted_count"] == 2
    assert len(read_dataset_index(project.root)["assets"]) == 4
    with pytest.raises(ValueError, match="Positive coverage"):
        survey_coverage_footprints(project.root, max_drawn=0)


def test_external_file_reattach_keeps_stable_identifiers(project, tmp_path):
    path = tmp_path / "same_original.fits"
    _original_wcs(path)
    first = register_external_fits(
        project.root, path, provider="HSC", band="r")
    second = register_external_fits(
        project.root, path, provider="HSC", band="r")
    assert first == second
    index = read_dataset_index(project.root)
    assert len(index["assets"]) == 1
    assert len(index["sources"]) == 1


def test_legacy_import_accepts_existing_object_root_and_nested_data(
        project):
    hit = BrickHit("dr9", "north", "0401p012", 0, 0, 0, 0, {"r": 1})
    for kind in ("image", "invvar", "nexp", "maskbits"):
        rel = product_for(
            hit, kind, None if kind == "maskbits" else "r").relative_path
        path = project.root / "data" / rel
        _original_wcs(path)
    first = import_legacy_fits(project.root, project.root)
    second = import_legacy_fits(project.root, project.root / "data")
    assert first["assets_added"] == 4
    assert second["assets_added"] == 0
    assert second["sources_total"] == 1
    assert all(a["fits_status"] == "NOT_VALIDATED" for a in
               read_dataset_index(project.root)["assets"])
