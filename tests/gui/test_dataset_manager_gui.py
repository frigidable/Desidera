"""Qt offline acceptance: project-bound Dataset Manager, no downloads."""
from __future__ import annotations

import pytest

pytest.importorskip("pytestqt")
pytest.importorskip("PySide6.QtWebEngineWidgets")

from lsbpipe.surveys.object_projects import create_object_project
from lsbpipe.surveys.dataset_registry import (
    register_external_fits, read_dataset_index)
from lsbpipe.surveys.legacy import SkyField
from lsbpipe.gui.dataset_manager import DatasetManager
from lsbpipe.gui.survey_explorer import SurveyExplorer
from scripts.offline_gui_fixture import offline_legacy_viewer


@pytest.fixture
def project(tmp_path, monkeypatch):
    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path / "session"))
    return create_object_project(
        "UGC 2162", tmp_path / "UGC_2162",
        field=SkyField(40.096208, 1.22925, 5.))


def test_workspace_navigation_to_dataset_manager_keeps_old_fits_tabs(
        qtbot, project):
    with offline_legacy_viewer():
        workspace = SurveyExplorer(restore_session=False, launch_remote=False)
    qtbot.addWidget(workspace)
    workspace.attach_object_project(project)
    workspace.show()
    assert workspace.dataset_manager_btn.isVisible()
    assert workspace.preview_tabs.count() == 2
    workspace.dataset_manager_btn.click()
    assert workspace.preview_tabs.currentWidget() is workspace.dataset_manager
    assert workspace.preview_tabs.count() == 3
    assert workspace.workflow_stage.text() == "02 / 04 · Data sources"
    assert not workspace.band_compare_btn.isVisible()
    assert not workspace.comparison_btn.isVisible()
    assert not workspace.data_toggle.isEnabled()
    workspace._toggle_data_drawer()
    assert not workspace.data_drawer.isVisible()
    assert workspace.dataset_manager.tree.topLevelItemCount() == 0
    assert "No indexed sources" in workspace.dataset_manager.source_summary.text()
    assert "not inspected" in workspace.dataset_manager.map_status.text().lower()
    assert not workspace.dataset_manager.map_btn.isEnabled()
    workspace.dataset_manager.back_btn.click()
    assert workspace.preview_tabs.currentWidget() is workspace.sky
    assert workspace.preview_tabs.indexOf(workspace.fits_explorer) == 1
    assert not (project.root / "datasets" / "index.json").exists()
    workspace._returning_to_hub = True
    workspace.close()


def test_dataset_manager_shows_actual_external_source_and_hash_action(
        qtbot, project, tmp_path):
    image = tmp_path / "hsc_external_i.fits"
    image.write_bytes(b"USER ORIGINAL BYTES")
    register_external_fits(project.root, image, provider="HSC", band="i")
    widget = DatasetManager()
    qtbot.addWidget(widget)
    widget.bind_project(project, read_dataset_index(project.root))
    widget.resize(1200, 700)
    widget.show()
    assert widget.tree.topLevelItemCount() == 1
    parent = widget.tree.topLevelItem(0)
    assert "HSC" in parent.text(0)
    assert parent.text(1) == "i"
    assert parent.childCount() == 1
    assert not widget.hash_btn.isEnabled()
    asset = parent.child(0)
    widget.tree.setCurrentItem(asset)
    assert widget.hash_btn.isEnabled()
    requested = []
    widget.checksumRequested.connect(requested.append)
    widget.hash_btn.click()
    assert requested == [asset.data(0, 0x0100)[1]]  # Qt.UserRole
    assert "not yet verified" in widget.details.text()
    assert widget.map_btn.isEnabled()
    assert image.read_bytes() == b"USER ORIGINAL BYTES"


def test_manager_real_wcs_footprint_display_is_read_only(
        qtbot, project):
    widget = DatasetManager()
    qtbot.addWidget(widget)
    widget.bind_project(project, read_dataset_index(project.root))
    original = widget.scene.itemsBoundingRect()
    assert original.width() >= 300  # target search square, arcseconds
    widget.set_coverage({
        "status": "GEOMETRY_ONLY_NOT_ADOPTED",
        "footprints": [{
            "asset_id": "asset-1", "source_id": "source-1",
            "provider": "Synthetic", "release": None, "region": None,
            "brick": None, "band": "r",
            "points": [[-100, -100], [100, -100],
                       [100, 100], [-100, 100]],
        }],
        "rejected": [],
        "omitted_count": 0,
    })
    assert "1 original IMAGE WCS headers" in widget.map_status.text()
    assert "geometry only" in widget.map_status.text()
    assert widget.scene.itemsBoundingRect().width() >= original.width()


def _legacy_index(project):
    """Screenshot-realistic registry: North complete, South missing g INVVAR."""
    import uuid
    index = {
        "project_id": project.project_id, "schema": "lsbpipe.dataset-index/1",
        "sources": [], "assets": [], "collections": [],
        "analysis_configurations": [], "runs": [],
    }
    for region in ("north", "south"):
        sid = f"source-{region}"
        index["sources"].append({
            "source_id": sid, "provider": "legacy_surveys",
            "release": "dr9", "region": region,
            "calibration_status": "UNASSESSED",
        })
        for band in ("g", "r", "z"):
            kinds = ["image", "invvar", "nexp"]
            if region == "south" and band == "g":
                kinds.remove("invvar")
            for kind in kinds:
                index["assets"].append({
                    "asset_id": f"{region}-{band}-{kind}",
                    "source_id": sid,
                    "brick": "0401p012", "band": band,
                    "product_kind": kind,
                    "checksum_status": "NOT_HASHED",
                    "fits_status": "NOT_VALIDATED", "sha256": None,
                })
        index["assets"].append({
            "asset_id": f"{region}-maskbits", "source_id": sid,
            "brick": "0401p012", "band": None, "product_kind": "maskbits",
            "checksum_status": "NOT_HASHED",
            "fits_status": "NOT_VALIDATED", "sha256": None,
        })
    return index


def test_ugc2162_incomplete_quartet_explicit_and_not_science_approved(
        qtbot, project):
    widget = DatasetManager()
    qtbot.addWidget(widget)
    widget.bind_project(project, _legacy_index(project))
    widget.show()
    assert widget.tree.topLevelItemCount() == 2
    north = widget.tree.topLevelItem(0)
    south = widget.tree.topLevelItem(1)
    assert "3/3 quartets indexed" == north.text(2)
    assert "2/3 quartets indexed" == south.text(2)
    assert "1 incomplete" in widget.source_summary.text()
    missing = [south.child(i) for i in range(south.childCount())
               if "missing invvar" in south.child(i).text(0)]
    assert len(missing) == 1
    assert missing[0].text(1) == "g"
    widget.tree.setCurrentItem(missing[0])
    assert not widget.hash_btn.isEnabled()
    assert "not registered locally" in widget.details.text()


def test_map_merges_only_identical_geometry_within_each_source(
        qtbot, project):
    from PySide6.QtWidgets import QGraphicsPolygonItem
    from PySide6.QtCore import Qt
    widget = DatasetManager()
    qtbot.addWidget(widget)
    index = _legacy_index(project)
    widget.bind_project(project, index)
    widget.resize(1280, 750)
    widget.show()
    polygon = [[-230, -230], [230, -230],
               [230, 230], [-230, 230]]
    # Six original headers, two source footprints; no filled scene overlay.
    rows = []
    for source in index["sources"]:
        for band in ("g", "r", "z"):
            rows.append({
                "asset_id": f"{source['region']}-{band}-image",
                "source_id": source["source_id"],
                "release": "dr9", "region": source["region"],
                "brick": "0401p012", "band": band,
                "points": polygon,
            })
    widget.set_coverage({
        "status": "GEOMETRY_ONLY_NOT_ADOPTED",
        "footprints": rows, "rejected": [], "omitted_count": 0,
    })
    paths = [item for item in widget.scene.items()
             if isinstance(item, QGraphicsPolygonItem)]
    assert len(paths) == 2
    assert all(p.brush().style() == Qt.NoBrush for p in paths)
    assert "6 original IMAGE WCS headers" in widget.map_status.text()
    assert "2 distinct source outlines" in widget.map_status.text()
    assert "bands g, r, z" in paths[0].toolTip()
    # Default is target-centered, not forced to show a giant 8' frame.
    target_scale = widget.map_view.transform().m11()
    widget.fit_all_btn.click()
    full_scale = widget.map_view.transform().m11()
    assert full_scale < target_scale
    widget.fit_target_btn.click()
    assert widget.map_view.transform().m11() > full_scale


def test_sha_verification_preserves_selected_asset_and_refreshes_details(
        qtbot, project):
    widget = DatasetManager()
    qtbot.addWidget(widget)
    original = _legacy_index(project)
    widget.bind_project(project, original)
    widget.show()
    north = widget.tree.topLevelItem(0)
    asset = north.child(0)
    aid = asset.data(0, 0x0100)[1]
    widget.tree.setCurrentItem(asset)
    widget.set_busy(True)
    changed = _legacy_index(project)
    record = next(a for a in changed["assets"] if a["asset_id"] == aid)
    record["sha256"] = "a" * 64
    record["checksum_status"] = "FILE_SHA256_VERIFIED"
    widget.set_index(changed)
    widget.set_busy(False)
    assert widget.tree.currentItem().data(0, 0x0100)[1] == aid
    assert "a" * 64 in widget.details.text()
    assert "Calculating" not in widget.details.text()
    assert widget.hash_btn.isEnabled()
