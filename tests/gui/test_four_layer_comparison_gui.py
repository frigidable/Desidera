"""Real Qt acceptance for FITS side controls and interactive 2x2 comparison."""
from __future__ import annotations
import numpy as np
import pytest

pytest.importorskip("pytestqt")
pytest.importorskip("PySide6.QtWebEngineWidgets")
pytest.importorskip("astropy")
from astropy.wcs import WCS
from PySide6.QtCore import Qt, QPoint
from PySide6.QtGui import QWheelEvent
from PySide6.QtWidgets import QWidget
from lsbpipe.gui.fits_explorer import FITSExplorer
from lsbpipe.gui.four_layer_comparison import FourLayerComparison
from lsbpipe.gui.survey_explorer import SurveyExplorer
from lsbpipe.surveys.legacy import SkyField
from scripts.offline_gui_fixture import offline_legacy_viewer


def data_fixture():
    field = SkyField(40.096208, 1.22925, 5.)
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [5., 5.]
    wcs.wcs.cdelt = [-0.262/3600., 0.262/3600.]
    wcs.wcs.crval = [field.ra_deg, field.dec_deg]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    image = np.arange(80, dtype=np.float32).reshape(8, 10)
    return {
        "status": "DIAGNOSTIC_NOT_ADOPTED",
        "source": "DR9 south · 0401p012 · r",
        "pixel_shape": image.shape, "band": "r", "wcs": wcs,
        "field": field,
        "hdu_indices": {"image": 1, "invvar": 1, "nexp": 1, "maskbits": 1},
        "data": {
            "image": image,
            "invvar": np.full_like(image, 32.),
            "nexp": np.full(image.shape, 3, dtype=np.int16),
            "maskbits": np.zeros(image.shape, dtype=np.int16),
        },
        "outlines": {"square": ((1, 1), (8, 1), (8, 6), (1, 6), (1, 1))},
    }


def test_fits_sidebar_occupies_left_column_and_diagnostics_hidden(qtbot):
    viewer = FITSExplorer()
    qtbot.addWidget(viewer)
    viewer.resize(1100, 700)
    viewer.show()
    assert viewer.workspace_splitter.count() == 2
    assert viewer.workspace_splitter.widget(0) is viewer.sidebar_scroll
    assert viewer.workspace_splitter.widget(1) is viewer.image_pane
    assert viewer.sidebar_scroll.width() < viewer.view.width()
    assert not viewer.detail_panel.isVisible()
    viewer.detail_toggle.click()
    assert viewer.detail_panel.isVisible()
    viewer.detail_toggle.click()
    assert not viewer.detail_panel.isVisible()
    viewer.set_dataset(data_fixture())
    qtbot.waitUntil(
        lambda: viewer._render_future is None and
                not viewer.base_item.pixmap().isNull(), timeout=8000)
    assert not viewer.base_item.pixmap().isNull()
    assert viewer.context.text().startswith("DR9 south")
    assert viewer.view.width() > viewer.sidebar_scroll.width()
    viewer.shutdown_renderer()


def test_four_layer_native_views_zoom_pan_hover_without_png(qtbot):
    compare = FourLayerComparison()
    qtbot.addWidget(compare)
    compare.resize(1050, 830)
    compare.show()
    dataset = data_fixture()
    initial = dataset["data"]["image"].copy()
    compare.set_dataset(dataset)
    qtbot.waitUntil(
        lambda: compare._future is None and all(
            not compare.images[k].pixmap().isNull()
            for k, _ in compare.LAYERS), timeout=10000)
    assert "not four optical filters" in compare.findChild(
        type(compare.status), "fits_four_layer_explanation").text()
    assert set(compare.views) == {"image", "invvar", "nexp", "maskbits"}
    assert all(compare.views[k].scene() is compare.scenes[k]
               for k in compare.views)
    compare._inspect("image", 3., 5.)
    assert "ICRS" in compare.readout.text()
    assert "IMAGE" in compare.readout.text()
    compare.views["image"].scale(1.5, 1.5)
    compare.views["image"]._relative_zoom = 1.5
    compare._sync_zoom(compare.views["image"])
    assert compare.views["invvar"].transform() == (
        compare.views["image"].transform())
    compare.fit_all.click()
    assert compare.views["image"]._relative_zoom == 1.
    np.testing.assert_array_equal(dataset["data"]["image"], initial)
    compare.shutdown_renderer()
    compare.shutdown_renderer()
    assert not compare._timer.isActive()


def test_four_layer_workspace_is_on_demand_and_closable(qtbot, monkeypatch):
    with offline_legacy_viewer():
        workspace = SurveyExplorer(restore_session=False, launch_remote=False)
    qtbot.addWidget(workspace)
    workspace.resize(1400, 820)
    workspace.show()
    data = data_fixture()
    monkeypatch.setattr(workspace, "_qa_context_current", lambda *_: True)
    workspace.fits_explorer.set_dataset(data)
    assert workspace.preview_tabs.count() == 2
    workspace.comparison_btn.click()
    assert workspace.preview_tabs.count() == 3
    assert workspace.preview_tabs.tabText(1) == "FITS Explorer"
    assert workspace.preview_tabs.tabText(2) == "Data quality · 4 maps"
    assert workspace.preview_tabs.currentWidget() is workspace.four_layer_comparison
    qtbot.waitUntil(
        lambda: workspace.four_layer_comparison._future is None and
                not workspace.four_layer_comparison.images["image"].pixmap().isNull(),
        timeout=10000)
    workspace.comparison_btn.click()
    assert workspace.preview_tabs.count() == 2
    assert workspace.preview_tabs.currentWidget() is workspace.fits_explorer
    workspace.fits_explorer.shutdown_renderer()
    workspace.four_layer_comparison.shutdown_renderer()


def test_four_layer_rejects_oversize_or_non_diagnostic_without_allocating(qtbot):
    widget = FourLayerComparison()
    qtbot.addWidget(widget)
    bad = data_fixture()
    bad["status"] = "SCIENCE_ADOPTED"
    with pytest.raises(ValueError, match="diagnostic"):
        widget.set_dataset(bad)
    bad = data_fixture()
    bad["pixel_shape"] = (3000, 3000)
    with pytest.raises(ValueError, match="safe viewer"):
        widget.set_dataset(bad)
    assert widget._future is None
    widget.shutdown_renderer()
