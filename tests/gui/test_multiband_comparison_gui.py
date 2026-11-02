"""Offline regression for original-FITS optical-band comparisons."""
from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("pytestqt")
pytest.importorskip("PySide6.QtWebEngineWidgets")
pytest.importorskip("astropy")

from astropy.wcs import WCS
from PySide6.QtCore import Qt
from lsbpipe.gui.multiband_comparison import (
    MultibandComparison, load_multiband_originals, MAX_TOTAL_PIXELS,
)
from lsbpipe.surveys.legacy import SkyField, BrickHit, product_for
from scripts.offline_gui_fixture import offline_legacy_viewer


def dataset(band, *, crpix=(5., 5.)):
    field = SkyField(40.096208, 1.22925, 5.)
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = list(crpix)
    wcs.wcs.cdelt = [-0.262/3600., 0.262/3600.]
    wcs.wcs.crval = [field.ra_deg, field.dec_deg]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    image = np.arange(80, dtype=np.float32).reshape(8, 10)
    image.setflags(write=False)
    return {
        "status": "DIAGNOSTIC_NOT_ADOPTED", "band": band,
        "source": f"DR9 north · 0401p012 · {band}",
        "field": field, "wcs": wcs,
        "pixel_shape": (8, 10),
        "outlines": {"square": ((1., 1.), (8., 1.),
                                (8., 6.), (1., 6.), (1., 1.))},
        "data": {
            "image": image,
            "invvar": np.full((8, 10), 32., dtype=np.float32),
            "nexp": np.full((8, 10), 3, dtype=np.int16),
            "maskbits": np.zeros((8, 10), dtype=np.int16),
        }
    }


def test_multiband_native_views_missing_i_and_shared_sky(qtbot):
    widget = MultibandComparison()
    qtbot.addWidget(widget)
    widget.resize(1100, 750)
    widget.show()
    original_g = dataset("g", crpix=(5., 5.))
    original_r = dataset("r", crpix=(6., 4.))
    pristine = original_g["data"]["image"].copy()
    widget.set_data({
        "status": "DIAGNOSTIC_NOT_ADOPTED",
        "datasets": {"g": original_g, "r": original_r},
        "missing": {"i": "not present", "z": "not present"},
    })
    qtbot.waitUntil(
        lambda: widget._future is None
                and not widget.images["g"].pixmap().isNull()
                and not widget.images["r"].pixmap().isNull(),
        timeout=10000)
    assert widget.views["g"].isVisible()
    assert not widget.views["i"].isVisible()
    assert widget.placeholders["i"].isVisible()
    assert "unavailable" in widget.status.text()
    assert widget.link.isChecked()
    widget.views["g"].scale(1.5, 1.5)
    widget.views["g"]._relative_zoom = 1.5
    widget._sync(widget.views["g"])
    assert widget.views["r"].transform().m11() > 0
    # Different WCS reference pixels must map the same sky position
    # to different native pixel indices; they must never be assumed equal.
    gx, gy = original_g["wcs"].all_world2pix(40.096208, 1.22925, 0)
    rx, ry = original_r["wcs"].all_world2pix(40.096208, 1.22925, 0)
    assert gx != pytest.approx(rx)
    assert gy != pytest.approx(ry)
    widget._inspect("g", 3., 5.)
    assert "ICRS" in widget.readout.text()
    assert "nearest native pixels" in widget.readout.text()
    np.testing.assert_array_equal(original_g["data"]["image"], pristine)
    widget.shutdown_renderer()
    assert not widget._timer.isActive()


def test_multiband_rejects_cross_brick_and_science_adoption(qtbot):
    widget = MultibandComparison()
    qtbot.addWidget(widget)
    g = dataset("g")
    r = dataset("r")
    r["source"] = "DR9 south · 0401p012 · r"
    with pytest.raises(ValueError, match="different objects or bricks"):
        widget.set_data({
            "status": "DIAGNOSTIC_NOT_ADOPTED",
            "datasets": {"g": g, "r": r},
        })
    g["status"] = "SCIENCE_ADOPTED"
    with pytest.raises(ValueError, match="diagnostic"):
        widget.set_data({
            "status": "DIAGNOSTIC_NOT_ADOPTED",
            "datasets": {"g": g},
        })
    assert widget._future is None
    widget.shutdown_renderer()


def test_multiband_loader_missing_bands_are_not_downloaded(
        tmp_path, monkeypatch):
    from lsbpipe.gui import fits_raster
    field = SkyField(40.096208, 1.22925, 5.)
    hit = BrickHit("dr9", "north", "0401p012", 0, 0, 0, 0,
                   {"g": 1, "r": 1})
    for band in ("g", "r"):
        for kind in ("image", "invvar", "nexp", "maskbits"):
            product = product_for(
                hit, kind, None if kind == "maskbits" else band)
            target = tmp_path / "data" / product.relative_path
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(b"test bytes")
    calls = []

    def fake_prepare(*, directory, field, hit, band, **kwargs):
        calls.append(band)
        return dataset(band)

    monkeypatch.setattr(fits_raster, "prepare_fits_explorer", fake_prepare)
    data = load_multiband_originals(
        directory=tmp_path, field=field, hit=hit, radii=(30., 60., 90.))
    assert sorted(data["datasets"]) == ["g", "r"]
    assert sorted(data["missing"]) == ["i", "z"]
    assert calls == ["g", "r"]
    assert data["status"] == "DIAGNOSTIC_NOT_ADOPTED"


def test_bad_band_is_reported_without_discarding_other_fits(
        tmp_path, monkeypatch):
    from lsbpipe.gui import fits_raster
    field = SkyField(40.096208, 1.22925, 5.)
    hit = BrickHit("dr9", "north", "0401p012", 0, 0, 0, 0,
                   {"g": 1, "r": 1})
    for band in ("g", "r"):
        for kind in ("image", "invvar", "nexp", "maskbits"):
            product = product_for(
                hit, kind, None if kind == "maskbits" else band)
            path = tmp_path / "data" / product.relative_path
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b"original placeholder")

    def fake_loader(*, band, **kwargs):
        if band == "r":
            raise ValueError("Mismatched WCS in r quartet")
        return dataset(band)

    monkeypatch.setattr(fits_raster, "prepare_fits_explorer", fake_loader)
    result = load_multiband_originals(
        directory=tmp_path, field=field, hit=hit, radii=(30., 60., 90.))
    assert list(result["datasets"]) == ["g"]
    assert "Mismatched WCS" in result["missing"]["r"]
    assert "Band unavailable" in result["missing"]["i"]


def test_multiband_workspace_on_demand_and_context_invalidation(
        qtbot, monkeypatch):
    from lsbpipe.gui.survey_explorer import SurveyExplorer
    with offline_legacy_viewer():
        workspace = SurveyExplorer(restore_session=False, launch_remote=False)
    qtbot.addWidget(workspace)
    workspace.show()
    assert workspace.preview_tabs.count() == 2
    assert workspace.band_compare_btn.objectName() == "fits_open_band_comparison"
    assert "quality" in workspace.comparison_btn.text().lower()
    workspace._multiband_ready({
        "status": "DIAGNOSTIC_NOT_ADOPTED",
        "datasets": {"g": dataset("g"), "r": dataset("r")},
        "missing": {"i": "not present", "z": "not present"},
    })
    assert workspace.preview_tabs.count() == 3
    assert workspace.preview_tabs.currentWidget() is workspace.multiband_comparison
    workspace._invalidate_qa_figure()
    assert workspace.preview_tabs.count() == 2
    assert not workspace.multiband_comparison.datasets
    workspace.fits_explorer.shutdown_renderer()
    workspace.four_layer_comparison.shutdown_renderer()
    workspace.multiband_comparison.shutdown_renderer()
