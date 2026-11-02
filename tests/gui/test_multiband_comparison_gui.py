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
    MultibandComparison, load_multiband_originals,
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
    assert not widget.cells["i"].isVisible()  # grz uses available width
    assert "unavailable" in widget.status.text()
    assert widget.scaling.currentData() == "log"
    assert widget.stretch.currentData() == "faint"
    widget.show_missing.setChecked(True)
    assert widget.cells["i"].isVisible()
    assert widget.placeholders["i"].isVisible()
    widget.show_missing.setChecked(False)
    assert not widget.cells["i"].isVisible()
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


def test_three_band_grz_layout_uses_three_columns_not_empty_quadrant(qtbot):
    widget = MultibandComparison()
    qtbot.addWidget(widget)
    widget.resize(1300, 800)
    widget.show()
    widget.set_data({
        "status": "DIAGNOSTIC_NOT_ADOPTED",
        "datasets": {b: dataset(b) for b in ("g", "r", "z")},
        "missing": {"i": "Not available in DR9 north"},
    })
    qtbot.waitUntil(
        lambda: widget._future is None and
                not widget.images["r"].pixmap().isNull(), timeout=10000)
    assert not widget.cells["i"].isVisible()
    assert widget.band_grid.columnCount() >= 3
    assert [widget.band_grid.getItemPosition(
        widget.band_grid.indexOf(widget.cells[b]))[:2]
        for b in ("g", "r", "z")] == [(0, 0), (0, 1), (0, 2)]
    assert "unavailable: I" in widget.status.text()
    widget.shutdown_renderer()


def test_grz_fit_all_preserves_native_field_and_angular_scale(qtbot):
    """Regression: three narrow panels cannot crop source FITS at Fit all."""
    from PySide6.QtCore import Qt
    from astropy.wcs.utils import proj_plane_pixel_scales
    widget = MultibandComparison()
    qtbot.addWidget(widget)
    widget.resize(1350, 870)
    widget.show()
    sources = {}
    shape = (384, 384)
    for band, crpix in (
            ("g", (192., 192.)), ("r", (190., 194.)),
            ("z", (195., 189.))):
        d = dataset(band, crpix=crpix)
        d["pixel_shape"] = shape
        d["data"] = {
            "image": np.broadcast_to(np.array([[0.1]], dtype="f4"), shape),
            "invvar": np.broadcast_to(np.array([[25.]], dtype="f4"), shape),
            "nexp": np.broadcast_to(np.array([[5]], dtype="i2"), shape),
            "maskbits": np.broadcast_to(np.array([[0]], dtype="i2"), shape),
        }
        sources[band] = d
    widget.set_data({
        "status": "DIAGNOSTIC_NOT_ADOPTED",
        "datasets": sources,
        "missing": {"i": "DR9 North has no i band"},
    })
    qtbot.waitUntil(
        lambda: (widget._future is None and not widget._fit_queued
                 and all(not widget.images[b].pixmap().isNull()
                         for b in sources)),
        timeout=10000)
    widget.fit_views()
    angular_pixel_scales = []
    world_centers = []
    for band in sources:
        view = widget.views[band]
        assert view.horizontalScrollBarPolicy() == Qt.ScrollBarAlwaysOff
        assert view.verticalScrollBarPolicy() == Qt.ScrollBarAlwaysOff
        displayed = view.mapToScene(
            view.viewport().rect()).boundingRect()
        native = widget.scenes[band].sceneRect()
        # Virtual margins belong to the VIEW only; all original FITS pixel
        # and overlay scene coordinates remain in the native 384x384 grid.
        assert widget.views[band].sceneRect().width() > native.width()
        assert native.width() == shape[1]
        assert native.height() == shape[0]
        # Within a 1-pixel numerical allowance, every *native* cutout
        # must be visible when 'Fit all' is requested.
        assert displayed.adjusted(-1, -1, 1, 1).contains(native)
        assert view._relative_zoom == pytest.approx(1.)
        angle = float(np.mean(np.abs(
            proj_plane_pixel_scales(sources[band]["wcs"].celestial))))
        angular_pixel_scales.append(angle / view.transform().m11())
        center = view.mapToScene(view.viewport().rect().center())
        ny = sources[band]["pixel_shape"][0]
        ra, dec = sources[band]["wcs"].all_pix2world(
            center.x(), ny - 1. - center.y(), 0)
        world_centers.append((float(np.asarray(ra)),
                              float(np.asarray(dec))))
    assert max(angular_pixel_scales) - min(angular_pixel_scales) < 1e-10
    # Qt viewport centering uses integer scrollbar coordinates even when
    # scrollbars are visually hidden. Approximate WCS center agreement
    # should therefore be evaluated in physical *screen pixels*, not
    # at a sub-device-pixel angular tolerance (native WCS is untouched).
    for ra, dec in world_centers[1:]:
        separation_deg = np.hypot(
            (ra - world_centers[0][0]) *
            np.cos(np.deg2rad(world_centers[0][1])),
            dec - world_centers[0][1])
        assert separation_deg <= 2.0 * max(angular_pixel_scales)

    source = widget.views["g"]
    source.scale(1.22, 1.22)
    source._relative_zoom = 1.22
    widget._sync(source)
    assert widget.views["r"]._relative_zoom == pytest.approx(1.22)
    assert widget.views["z"]._relative_zoom == pytest.approx(1.22)
    widget.fit_views()
    assert all(widget.views[b]._relative_zoom == 1. for b in sources)
    widget.shutdown_renderer()


def test_ugc2162_sized_cutouts_fit_after_hidden_workspace_tab_is_shown(qtbot):
    """Reproduce real 1677x895 UGC 2162 screenshot, not tiny 384px Qt cases.

    All three 1833x1833 native IMAGE cutouts must fit narrow ~440px
    viewports after the previously hidden comparison tab becomes visible.
    """
    from lsbpipe.gui.survey_explorer import SurveyExplorer
    from scripts.offline_gui_fixture import offline_legacy_viewer
    with offline_legacy_viewer():
        workspace = SurveyExplorer(restore_session=False, launch_remote=False)
    qtbot.addWidget(workspace)
    workspace.resize(1677, 895)
    shape = (1833, 1833)
    datasets = {}
    for band in ("g", "r", "z"):
        d = dataset(band, crpix=(917., 917.))
        d["pixel_shape"] = shape
        source = np.array([[0.08]], dtype=np.float32)
        d["data"] = {
            "image": np.broadcast_to(source, shape),
            "invvar": np.broadcast_to(
                np.array([[40.]], dtype=np.float32), shape),
            "nexp": np.broadcast_to(
                np.array([[4]], dtype=np.int16), shape),
            "maskbits": np.broadcast_to(
                np.array([[0]], dtype=np.int16), shape),
        }
        datasets[band] = d
    # This matches production: set_data() is called while the comparison
    # tab has not yet been inserted/shown, and Qt allocates width later.
    workspace._multiband_ready({
        "status": "DIAGNOSTIC_NOT_ADOPTED",
        "datasets": datasets,
        "missing": {"i": "DR9 north has no i"},
    })
    workspace.show()
    qtbot.waitUntil(
        lambda: workspace.multiband_comparison._future is None and
                all(not workspace.multiband_comparison.images[b]
                        .pixmap().isNull() for b in datasets),
        timeout=25000)
    qtbot.wait(100)
    comparison = workspace.multiband_comparison
    assert workspace.preview_tabs.currentWidget() is comparison
    assert all(comparison.views[b].viewport().width() < 550
               for b in datasets)
    # Initial automatic Fit all, after Qt layout settles.
    for band in datasets:
        view = comparison.views[band]
        native = comparison.scenes[band].sceneRect()
        visible = view.mapToScene(view.viewport().rect()).boundingRect()
        assert visible.adjusted(-2, -2, 2, 2).contains(native), (
            f"{band}: cropped at startup: viewport={view.viewport().size()}, "
            f"visible={visible}, native={native}, "
            f"transform={view.transform().m11()}")
        assert (view.transform().m11() <=
                (view.viewport().width() - 4) / shape[1] + 0.001)
    # Repeat after explicit Fit all and a resize at the actual monitor aspect.
    comparison.fit_views()
    workspace.resize(1550, 850)
    qtbot.wait(100)
    comparison.fit_views()
    for band in datasets:
        view = comparison.views[band]
        visible = view.mapToScene(view.viewport().rect()).boundingRect()
        assert visible.adjusted(-2, -2, 2, 2).contains(
            comparison.scenes[band].sceneRect()), band
    workspace.close()
    workspace.fits_explorer.shutdown_renderer()
    workspace.four_layer_comparison.shutdown_renderer()
    workspace.multiband_comparison.shutdown_renderer()


def test_multiband_faint_log_matches_single_band_transfer_and_is_reversible(
        qtbot):
    from lsbpipe.gui.fits_raster import rgb_for_layer
    from lsbpipe.gui.multiband_preview import display_preview
    widget = MultibandComparison()
    qtbot.addWidget(widget)
    widget.resize(1200, 790)
    widget.show()
    d = dataset("r")
    yy, xx = np.mgrid[:128, :128]
    rng = np.random.default_rng(1909)
    background = rng.normal(-0.01, 0.018, (128, 128))
    diffuse = 0.1 * np.exp(-(
        ((xx - 64.) / 27.) ** 2 + ((yy - 64.) / 19.) ** 2))
    img = (background + diffuse).astype(np.float32)
    img[110, 5:120] = 3.0  # strong bright artifact in field
    img.setflags(write=False)
    d["data"]["image"] = img
    d["data"]["invvar"] = np.ones_like(img)
    d["data"]["nexp"] = np.ones_like(img, dtype=np.int16)
    d["data"]["maskbits"] = np.zeros_like(img, dtype=np.int16)
    d["pixel_shape"] = img.shape
    original = img.copy()

    rgb, factor, lo, hi, mode = widget._prepare_overview(
        d, "r", "faint", "log")
    assert factor == 1
    assert mode == "log" and 0 < lo < hi
    np.testing.assert_array_equal(
        rgb, rgb_for_layer(
            img, "image", "r", contrast="faint", scaling="log"))
    assert float(rgb[64, 64, 0]) > float(rgb[10, 10, 0])
    assert np.any(img <= 0)
    assert np.all(rgb[img <= 0, 0] == 0)

    widget.set_data({
        "status": "DIAGNOSTIC_NOT_ADOPTED",
        "datasets": {"r": d},
        "missing": {b: "Unavailable" for b in ("g", "i", "z")},
    })
    qtbot.waitUntil(
        lambda: widget._future is None and not
                widget.images["r"].pixmap().isNull(), timeout=10000)
    assert "log" in widget.preview_info["r"].text()
    assert "nMgy" in widget.preview_info["r"].text()
    widget.scaling.setCurrentIndex(1)  # signed asinh comparison
    qtbot.waitUntil(
        lambda: widget._future is None and "asinh" in
                widget.preview_info["r"].text(), timeout=10000)
    widget.scaling.setCurrentIndex(0)
    qtbot.waitUntil(
        lambda: widget._future is None and "log" in
                widget.preview_info["r"].text(), timeout=10000)
    np.testing.assert_array_equal(img, original)
    assert not img.flags.writeable
    widget.shutdown_renderer()


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


def test_bounded_preview_area_average_and_originals_are_immutable():
    from lsbpipe.gui.multiband_preview import display_preview
    image = np.asarray([
        [1., 2., 3., 4., 5.],
        [3., np.nan, 5., 6., 7.],
        [10., 12., 14., 16., 18.],
        [20., 22., 24., 26., 28.],
        [40., 42., 44., 46., 48.],
    ], dtype=np.float32)
    # max_side must be >=32. Test non-divisible dimensions using a 70x70
    # input, with exactly controlled 2x2 block averages.
    large = np.full((70, 70), 2., dtype=np.float32)
    large[0:2, 0:2] = [1., 3.], [5., 7.]
    large[0, 0] = np.nan
    original = large.copy()
    large.setflags(write=False)
    overview, factor = display_preview(large, max_side=40)
    assert factor == 2
    assert overview.shape == (35, 35)
    assert overview[0, 0] == pytest.approx(5.)
    assert overview[-1, -1] == pytest.approx(2.)
    np.testing.assert_array_equal(large, original)
    assert not large.flags.writeable
    same, no_bin = display_preview(image)
    assert same is image and no_bin == 1


def test_large_g_r_z_cutouts_above_old_eight_million_limit_display(
        qtbot):
    from lsbpipe.gui.multiband_preview import (
        MAX_PREVIEW_SIDE, native_dataset_bytes)
    widget = MultibandComparison()
    qtbot.addWidget(widget)
    widget.resize(1050, 760)
    widget.show()
    shape = (2048, 2048)
    sources = {}
    for band in ("g", "r", "z"):
        d = dataset(band)
        d["pixel_shape"] = shape
        # Nominal ~50 MiB of science arrays per band, despite tiny actual
        # allocation from read-only broadcast fixtures.
        d["data"] = {
            "image": np.broadcast_to(
                np.array([[3.]], dtype=np.float32), shape),
            "invvar": np.broadcast_to(
                np.array([[25.]], dtype=np.float32), shape),
            "nexp": np.broadcast_to(
                np.array([[5]], dtype=np.int16), shape),
            "maskbits": np.broadcast_to(
                np.array([[0]], dtype=np.int16), shape),
        }
        sources[band] = d
    assert sum(np.prod(d["pixel_shape"]) for d in sources.values()) > 8_000_000
    assert sum(native_dataset_bytes(d) for d in sources.values()) < 512 * 2**20
    widget.set_data({
        "status": "DIAGNOSTIC_NOT_ADOPTED",
        "datasets": sources,
        "missing": {"i": "Not available in DR9 north"},
    })
    qtbot.waitUntil(
        lambda: widget._future is None and all(
            not widget.images[band].pixmap().isNull()
            for band in sources), timeout=25000)
    for band in sources:
        source_ny, source_nx = sources[band]["pixel_shape"]
        item = widget.images[band]
        pix = item.pixmap()
        assert pix.width() <= MAX_PREVIEW_SIDE
        assert pix.height() <= MAX_PREVIEW_SIDE
        assert (pix.width() * item.transform().m11()
                == pytest.approx(source_nx))
        assert (pix.height() * item.transform().m22()
                == pytest.approx(source_ny))
        assert widget.scenes[band].sceneRect().width() == source_nx
        assert "bin" in widget.preview_info[band].text()
        assert widget.preview_factors[band] == 2
    # Native WCS/hovers still sample ORIGINAL data, not preview pixels.
    widget._inspect("r", 1150., 1200.)
    assert "r: 3 nMgy" in widget.readout.text()
    assert "nearest native pixels" in widget.readout.text()
    assert not sources["r"]["data"]["image"].flags.writeable
    widget.shutdown_renderer()


def test_multiband_rejects_excessive_native_array_memory(qtbot):
    from lsbpipe.gui.multiband_preview import MAX_NATIVE_BYTES
    widget = MultibandComparison()
    qtbot.addWidget(widget)
    d = dataset("g")
    shape = (10000, 10000)
    d["pixel_shape"] = shape
    d["data"] = {
        "image": np.broadcast_to(np.array([[1.]], dtype=np.float32), shape),
        "invvar": np.broadcast_to(np.array([[2.]], dtype=np.float32), shape),
        "nexp": np.broadcast_to(np.array([[3]], dtype=np.int16), shape),
        "maskbits": np.broadcast_to(np.array([[0]], dtype=np.int16), shape),
    }
    assert sum(a.nbytes for a in d["data"].values()) > MAX_NATIVE_BYTES
    with pytest.raises(ValueError, match="original FITS arrays|Original FITS"):
        widget.set_data({
            "status": "DIAGNOSTIC_NOT_ADOPTED", "datasets": {"g": d}})
    assert widget._future is None
    widget.shutdown_renderer()


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
