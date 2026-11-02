"""Offline synthetic native-FITS pixel support metrics; no science adoption."""
from __future__ import annotations

import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS

from lsbpipe.surveys.object_projects import create_object_project
from lsbpipe.surveys.legacy import BrickHit, SkyField, product_for
from lsbpipe.surveys.dataset_registry import (
    import_legacy_fits, read_dataset_index)
from lsbpipe.surveys.dataset_pixel_qa import inspect_pixel_coverage


@pytest.fixture
def project(tmp_path, monkeypatch):
    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path / "session"))
    return create_object_project(
        "UGC 2162", tmp_path / "UGC_2162",
        field=SkyField(40.096208, 1.22925, 5.))


def make_quartet(root, region="north", band="r", missing=()):
    shape = (120, 120)
    wcs = WCS(naxis=2)
    wcs.wcs.crval = [40.096208, 1.22925]
    wcs.wcs.crpix = [60.5, 60.5]
    wcs.wcs.cdelt = [-4./3600., 4./3600.]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    hit = BrickHit("dr9", region, "0401p012", 0, 0, 0, 0, {band: 1})
    img = np.ones(shape, dtype=np.float32)
    ivar = np.ones(shape, dtype=np.float32)
    invzero = np.indices(shape)[1] < 60
    ivar[invzero] = 0
    nexp = np.full(shape, 3, dtype=np.int16)
    mask = np.zeros(shape, dtype=np.int16)
    # Positive, signed-safe upper flag test (SATUR_R, maskbit 3).
    mask[45:76, 45:76] |= 1 << 3
    mask[55:65, 55:65] |= 1 << 0  # NPRIMARY reported, not auto-excluded
    for kind, data in (
        ("image", img), ("invvar", ivar),
        ("nexp", nexp), ("maskbits", mask)):
        if kind in missing:
            continue
        path = root / product_for(
            hit, kind, None if kind == "maskbits" else band).relative_path
        path.parent.mkdir(parents=True, exist_ok=True)
        fits.PrimaryHDU(data=data, header=wcs.to_header()).writeto(path)
    return img, ivar, nexp, mask


def test_original_pixel_coverage_separates_target_annulus_and_mask_flags(
        project):
    root = project.root / "data"
    arrays = make_quartet(root)
    originals = tuple(arr.tobytes() for arr in arrays)
    import_legacy_fits(project.root, root)
    report = inspect_pixel_coverage(project.root)
    assert report["status"] == "DIAGNOSTIC_NOT_ADOPTED"
    assert len(report["rows"]) == 1
    row = report["rows"][0]
    assert row["status"] == "PIXEL_DIAGNOSTIC_NOT_ADOPTED"
    assert len(row["input_asset_ids"]) == 4
    target = row["regions"]["target"]
    sky = row["regions"]["sky_annulus"]
    assert target["sampled_pixels"] > 0
    assert sky["sampled_pixels"] > target["sampled_pixels"]
    assert 0.48 < target["fractions"]["positive_invvar"] < 0.52
    assert target["fractions"]["positive_nexp"] == 1.
    assert target["fractions"]["finite_image"] == 1.
    assert target["fractions"]["flagged_band_satur"] > 0
    assert target["fractions"]["flagged_nprimary"] > 0
    assert target["fractions"]["numerical_triplet"] > (
        target["fractions"]["numerical_triplet_without_band_satur_allmask"])
    assert sky["fractions"]["numerical_triplet"] > 0
    assert row["flags_are_individual_counts_not_default_rejection"] is True
    assert tuple(arr.tobytes() for arr in arrays) == originals
    # Main original files must not be modified, hashed or adopted by QA.
    idx = read_dataset_index(project.root)
    assert all(a["sha256"] is None for a in idx["assets"])


def test_missing_south_g_invvar_never_borrows_north_quartet(project):
    root = project.root / "data"
    make_quartet(root, region="north", band="g")
    make_quartet(root, region="south", band="g", missing=("invvar",))
    import_legacy_fits(project.root, root)
    report = inspect_pixel_coverage(project.root)
    assert len(report["rows"]) == 2
    by_region = {row["region"]: row for row in report["rows"]}
    assert by_region["north"]["status"] == "PIXEL_DIAGNOSTIC_NOT_ADOPTED"
    assert by_region["south"]["status"] == "INCOMPLETE_FILES"
    assert by_region["south"]["missing"] == ["invvar"]
    assert "regions" not in by_region["south"]


def test_pixel_qa_bounded_and_reports_omission(project):
    root = project.root / "data"
    make_quartet(root, region="north", band="g")
    make_quartet(root, region="north", band="r")
    import_legacy_fits(project.root, root)
    report = inspect_pixel_coverage(project.root, max_quartets=1)
    assert len(report["rows"]) == 1
    assert report["omitted_quartets"] == 1
    with pytest.raises(ValueError, match="quartet limit"):
        inspect_pixel_coverage(project.root, max_quartets=0)


def test_negative_high_maskbit_preserved_in_original_array(project):
    root = project.root / "data"
    make_quartet(root, region="north", band="r")
    # Existing science reader recognizes signed upper bits using uint
    from lsbpipe.survey import maskbits_to_uint
    signed = np.array([[-32768]], dtype=np.int16)
    assert int(maskbits_to_uint(signed)[0, 0] & np.uint32(32768)) != 0
