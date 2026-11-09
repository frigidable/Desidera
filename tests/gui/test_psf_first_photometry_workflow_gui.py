"""Offline Qt acceptance: operational PSF workflow has no fake stages."""
from __future__ import annotations

import pytest

pytest.importorskip("pytestqt")
pytest.importorskip("PySide6")

from PySide6.QtWidgets import QGroupBox
from lsbpipe.gui.photometry_workspace import PhotometryWorkspace


def test_photometry_uses_real_psf_workflow_without_unimplemented_stage_cards(qtbot):
    w = PhotometryWorkspace()
    qtbot.addWidget(w)
    w.resize(1366, 768)
    w.show()
    # The eight proposal-only stages were not executable. Keep the
    # scientific method and planned pipeline in docs, not fake GUI cards.
    assert not hasattr(w, "workflow_toggle")
    assert not hasattr(w, "workflow_details")
    assert not hasattr(w, "method_details_toggle")
    assert not any(
        w.findChild(QGroupBox, f"photometry_stage_{i:02d}") is not None
        for i in range(1, 9))
    assert w.band_selector.count() == 4
    assert tuple(w.psf0_panels) == ("g", "r", "i", "z")
    assert w.psf0_panels["g"].run_btn.text() == "Build empirical PSF"
    assert w.psf0_panels["g"].experiment.currentData() == (
        "outer_sky_isolation_footprint")
    assert "no adopted PSF, sky or galaxy photometry" in w.subtitle.text()
    assert w.table.rowCount() == 0
    assert w._project_root is None
