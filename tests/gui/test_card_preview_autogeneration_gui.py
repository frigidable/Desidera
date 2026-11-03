"""P10c production Qt: selected galaxy generates local preview once."""
import pytest

pytest.importorskip("pytestqt")
pytest.importorskip("PySide6.QtWebEngineWidgets")

import numpy as np

from scripts.offline_gui_fixture import offline_legacy_viewer
from lsbpipe.gui import project_hub as hub_module
from lsbpipe.surveys.legacy import SkyField, BrickHit, product_for
from lsbpipe.surveys.object_projects import create_object_project, save_object_project


def test_selected_leo_p_builds_cache_offline_and_reopens_without_fits(
        qtbot, tmp_path, monkeypatch):
    fits = pytest.importorskip("astropy.io.fits")
    WCS = pytest.importorskip("astropy.wcs").WCS
    p = create_object_project(
        "Leo P", tmp_path / "Leo_P",
        field=SkyField(155.438012, 18.088025, 2.))
    p.selected_dataset = ("dr11", "south", "1554p180", "r")
    p.data_root = p.root
    save_object_project(p)
    manifest = p.root / "lsbpipe-object.json"
    saved = manifest.read_bytes()

    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [150, 150]
    wcs.wcs.cdelt = [-1 / 3600, 1 / 3600]
    wcs.wcs.crval = [155.438012, 18.088025]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    yy, xx = np.mgrid[:300, :300]
    source = p.root / product_for(
        BrickHit("dr11", "south", "1554p180",
                 0, 0, 0, 0, {"r": 1}), "image", "r").relative_path
    source.parent.mkdir(parents=True)
    fits.HDUList([
        fits.PrimaryHDU(),
        fits.CompImageHDU(
            data=(np.exp(-((xx - 150)**2 + (yy - 150)**2)
                         / (2 * 20**2))).astype("f4"),
            header=wcs.to_header()),
    ]).writeto(source)
    before_fits = source.read_bytes()

    with offline_legacy_viewer():
        hub = hub_module.ProjectHub()
    qtbot.addWidget(hub)
    hub.show()
    qtbot.waitUntil(
        lambda: (hub._preview_future is None
                 and hub.galaxy_card.snapshot is not None
                 and hub.galaxy_card.snapshot.project_id == p.project_id
                 and hub.galaxy_card.snapshot.thumbnail is not None),
        timeout=15000)
    assert hub.galaxy_card.preview.pixmap() is not None
    assert not hub.galaxy_card.preview.pixmap().isNull()
    assert "display only" in hub.galaxy_card.preview.toolTip().lower()
    assert not hub.galaxy_card.preview_note.isVisible()
    assert hub.preview_btn.isEnabled()
    assert manifest.read_bytes() == saved
    assert source.read_bytes() == before_fits

    def forbidden(*args):
        pytest.fail("Reopening cached Galaxy Projects must not decode FITS")

    monkeypatch.setattr(hub_module, "build_local_preview", forbidden)
    with offline_legacy_viewer():
        reopened = hub_module.ProjectHub()
    qtbot.addWidget(reopened)
    reopened.show()
    assert reopened._preview_future is None
    assert reopened.galaxy_card.snapshot.thumbnail is not None
    assert reopened.galaxy_card.preview.pixmap() is not None
    assert not reopened.galaxy_card.preview.pixmap().isNull()
    assert source.read_bytes() == before_fits
