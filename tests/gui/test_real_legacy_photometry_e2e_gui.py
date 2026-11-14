"""Opt-in native DESI Legacy brick acceptance: NEVER a synthetic science PASS.

Usage requires an existing on-disk archive quartet and same-brick Tractor
catalogue. All original files are opened read-only. The project, mask drafts,
GUI receipts and optional screenshots are written to pytest's temp directory.
The normal CI suite SKIPS this test unless explicitly given a real data root
and an independently established brick name. A locally claimed archive
source is not remotely authenticated by this test.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("pytestqt")
pytest.importorskip("PySide6")
pytest.importorskip("astropy")
pytest.importorskip("photutils")

from PySide6.QtWidgets import QApplication

from lsbpipe.gui.project_hub import ProjectHub
from lsbpipe.surveys.dataset_integrity import (
    _assert_safe_original_paths, _expected_paths, _sha256_file, verify_dataset,
)
from lsbpipe.surveys.legacy import BrickHit, SkyField, product_for
from lsbpipe.surveys.object_projects import create_object_project, save_object_project
from lsbpipe.photometry.mask_drafts import read_active_mask_draft


def _archive_test_input():
    root = os.environ.get("DESIDERA_REAL_FITS_ROOT", "").strip()
    if not root and os.environ.get("DESIDERA_REAL_BRICK", "").strip():
        pytest.fail("DESIDERA_REAL_FITS_ROOT is missing for the requested archive brick")
    if not root:
        pytest.skip(
            "Authentic Legacy FITS not supplied: set DESIDERA_REAL_FITS_ROOT "
            "and DESIDERA_REAL_BRICK for explicit archival acceptance")
    brick = os.environ.get("DESIDERA_REAL_BRICK", "").strip()
    if not re.fullmatch(r"[0-9]{4}[pm][0-9]{3}", brick):
        pytest.fail("DESIDERA_REAL_BRICK must be an exact verified 8-character brick")
    release = os.environ.get("DESIDERA_REAL_RELEASE", "dr11").strip()
    region = os.environ.get("DESIDERA_REAL_REGION", "south").strip()
    band = os.environ.get("DESIDERA_REAL_BAND", "r").strip()
    if release not in ("dr9", "dr10", "dr11"):
        pytest.fail("Unsupported source release")
    if region not in ("north", "south") or (release == "dr10" and region != "south"):
        pytest.fail("Unsupported source region for release")
    if band not in ("g", "r", "i", "z") or (region == "north" and band == "i"):
        pytest.fail("Unsupported optical band for source region")
    data_root = Path(root).expanduser().resolve()
    if not data_root.is_dir():
        pytest.fail(f"Original archive root is missing: {data_root}")
    return data_root, (release, region, brick, band)


def _archive_paths(root, key):
    paths = _expected_paths(root, key)
    release, region, brick, band = key
    hit = BrickHit(release, region, brick, 0., 0., 0., 0., {band: 1})
    paths["tractor"] = root / product_for(hit, "tractor").relative_path
    for kind, path in paths.items():
        if path.is_symlink() or not path.is_file() or not path.stat().st_size:
            pytest.fail(f"Missing/unsafe original archive input ({kind}): {path}")
    _assert_safe_original_paths(root, paths)
    return paths


def _source_hashes(paths, root):
    return {name: _sha256_file(path, data_root=root)["sha256_local"]
            for name, path in paths.items()}



def test_real_archive_test_is_explicitly_opt_in(monkeypatch):
    monkeypatch.delenv("DESIDERA_REAL_FITS_ROOT", raising=False)
    monkeypatch.delenv("DESIDERA_REAL_BRICK", raising=False)
    with pytest.raises(pytest.skip.Exception, match="not supplied"):
        _archive_test_input()


def test_real_archive_partial_or_invalid_input_fails_closed(monkeypatch, tmp_path):
    monkeypatch.delenv("DESIDERA_REAL_FITS_ROOT", raising=False)
    monkeypatch.setenv("DESIDERA_REAL_BRICK", "0401p012")
    with pytest.raises(pytest.fail.Exception, match="ROOT is missing"):
        _archive_test_input()
    monkeypatch.setenv("DESIDERA_REAL_FITS_ROOT", str(tmp_path))
    monkeypatch.setenv("DESIDERA_REAL_BRICK", "not_a_brick")
    with pytest.raises(pytest.fail.Exception, match="verified 8-character brick"):
        _archive_test_input()
    monkeypatch.setenv("DESIDERA_REAL_BRICK", "0401p012")
    monkeypatch.setenv("DESIDERA_REAL_BAND", "")
    with pytest.raises(pytest.fail.Exception, match="Unsupported optical band"):
        _archive_test_input()


def test_real_archive_missing_original_fits_is_failure_not_skip(
        monkeypatch, tmp_path):
    monkeypatch.setenv("DESIDERA_REAL_FITS_ROOT", str(tmp_path))
    monkeypatch.setenv("DESIDERA_REAL_BRICK", "0401p012")
    monkeypatch.setenv("DESIDERA_REAL_BAND", "r")
    root, key = _archive_test_input()
    with pytest.raises(pytest.fail.Exception, match="Missing/unsafe original"):
        _archive_paths(root, key)


def test_supplied_archive_brick_full_gui_route(qtbot, tmp_path, monkeypatch):
    root, key = _archive_test_input()
    paths = _archive_paths(root, key)
    original_hashes = _source_hashes(paths, root)
    release, region, brick, band = key

    # This gate checks full, native-brick source products rather than Viewer
    # cutouts or synthetic 384-pixel stand-ins. A valid FITS is NOT proof of
    # authentic archive provenance or scientifically adopted photometry.
    entry = {
        "release": release, "region": region, "brick": brick, "band": band,
        "status": "PRESENT_UNVERIFIED",
        "products": {kind: {"relative_path": str(p.relative_to(root))}
                     for kind, p in paths.items() if kind != "tractor"},
    }
    audit = verify_dataset(root, entry, hash_inputs=False)
    assert audit["status"] == "FITS_WCS_VALIDATED_DIAGNOSTIC", audit
    assert audit["shape"] == [3600, 3600], (
        "Expected full 3600x3600 Legacy brick; cutouts cannot certify "
        "original-brick GUI acceptance", audit["shape"])

    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path / "session"))
    monkeypatch.setenv("DESIDERA_DISABLE_WEBENGINE", "1")
    monkeypatch.setenv("DESIDERA_DISABLE_AUTOMATIC_RGB", "1")
    project = create_object_project(
        "UGC 2162 archival acceptance", tmp_path / "object_project",
        field=SkyField(40.096208, 1.229250, 8.0))
    project.selected_dataset = key
    project.data_root = root
    project.target_radius_arcsec = 30.0
    save_object_project(project)
    before_passport = json.loads(
        (project.root / "lsbpipe-object.json").read_text(encoding="utf-8"))

    # Exercise actual Project Hub -> Galaxy Workspace -> Photometry controls.
    hub = ProjectHub()
    qtbot.addWidget(hub)
    hub.show()
    hub.filter.setText("UGC 2162 archival acceptance")
    assert hub.projects.currentItem() is not None
    hub.galaxy_card.open_button.click()
    assert hub.workspace is not None
    explorer = hub.workspace
    assert explorer._object_project.project_id == project.project_id
    explorer.photometry_btn.click()
    widget = explorer.photometry_workspace
    assert explorer.preview_tabs.currentWidget() is widget
    scene, sky, aperture = (
        widget.scene_editor, widget.sky_stage, widget.aperture_stage)
    assert scene._loaded is None
    widget.restore_navigation(stage="scene", band=band)
    assert scene.load_btn.isEnabled()
    scene.width.setValue(8.0)
    scene.load_btn.click()
    qtbot.waitUntil(lambda: scene._future is None, timeout=120000)
    assert scene._loaded is not None, scene.status.text()
    native = scene._loaded
    assert native["original_fits_provenance_status"] == (
        "FULL_ORIGINAL_FITS_SHA256_VERIFIED_AT_SCENE_LOAD")
    assert native["source_quartet_key"] == list(key)
    assert native["original_fits_sha256_local"] == {
        k: original_hashes[k] for k in ("image", "invvar", "nexp", "maskbits")}
    assert all(not native["data"][k].flags.writeable for k
               in ("image", "invvar", "nexp", "maskbits"))
    ny, nx = native["pixel_shape"]
    assert ny > 300 and nx > 300
    midpoint = native["data"]["image"][ny // 2, nx // 2]
    assert np.isfinite(midpoint), "Target center lacks valid IMAGE pixels"

    # One physical mouse-drawn exclusion, saved only to the TEMP project.
    scene.scene_task_buttons["draw"].click()
    scene.mode.setCurrentIndex(scene.mode.findData("circle"))
    scene.purpose.setCurrentIndex(scene.purpose.findData("exclude"))
    scene.scope.setCurrentIndex(scene.scope.findData("band"))
    x, y = nx * .70, ny * .68
    scene._on_click(x + .5, ny - y - .5)
    scene._on_click(x + 7.5, ny - y - .5)
    assert scene.draft is not None and len(scene.draft.edits) == 1
    scene.save_btn.click()
    stored, digest = read_active_mask_draft(project.root)
    assert stored is not None and digest == scene.draft.digest

    # The original, *unmodified* same-brick Tractor BINTABLE must be used.
    scene.scene_task_buttons["isophote"].click()
    scene.iso_prepare.click()
    qtbot.waitUntil(lambda: scene._future is None, timeout=120000)
    assert scene._contamination is not None, scene.iso_mask_status.text()
    contamination = scene._contamination
    assert contamination["catalog_sha256"] == original_hashes["tractor"]
    assert contamination["draft_digest"] == digest
    assert contamination["source_quartet_key"] == list(key)
    assert contamination["stats"]["catalog_rows"] > 0

    widget.restore_navigation(stage="sky", band=band)
    assert sky.run_btn.isEnabled(), sky.status.text()
    sky.run_btn.click()
    qtbot.waitUntil(lambda: sky._future is None, timeout=180000)
    assert sky._result is not None, sky.status.text()
    sky_result = sky._result
    assert sky_result["catalog_sha256"] == original_hashes["tractor"]
    assert sky_result["science_status"] == "UNADOPTED_DIAGNOSTIC_ONLY"
    assert sky_result.get("sky_model_adopted") is False

    widget.restore_navigation(stage="aperture", band=band)
    assert aperture.run_btn.isEnabled(), aperture.status.text()
    aperture.run_btn.click()
    qtbot.waitUntil(lambda: aperture._future is None, timeout=180000)
    assert aperture._result is not None, aperture.status.text()
    aperture_result = aperture._result
    assert aperture_result["source_quartet_key"] == list(key)
    assert aperture_result["source_fits_hashes_verified"]
    assert aperture_result["original_fits_sha256_local"] == {
        k: original_hashes[k] for k in ("image", "invvar", "nexp", "maskbits")}
    assert not aperture_result["adopted_for_science"]
    assert not aperture_result["masked_light_reconstructed"]

    # Editing the mask invalidates prior sky and aperture evidence.
    widget.restore_navigation(stage="scene", band=band)
    scene.scene_task_buttons["draw"].click()
    scene._on_click(x + .5, ny - (y - 25.) - .5)
    scene._on_click(x + 7.5, ny - (y - 25.) - .5)
    assert len(scene.draft.edits) == 2
    sky.refresh_context()
    aperture.refresh_context()
    assert sky._result is None and aperture._result is None
    scene.save_btn.click()
    _, edited_digest = read_active_mask_draft(project.root)
    assert edited_digest != digest

    after_passport = json.loads(
        (project.root / "lsbpipe-object.json").read_text(encoding="utf-8"))
    assert {k: v for k, v in after_passport.items()
            if k != "last_opened_at"} == {
                k: v for k, v in before_passport.items()
                if k != "last_opened_at"}
    assert _source_hashes(paths, root) == original_hashes

    # All receipt/capture bytes go OUTSIDE the archive, in an explicit
    # test-only output path. Never embed large source FITS in the repository.
    output = os.environ.get("DESIDERA_REAL_E2E_REPORT_DIR")
    if output:
        report_dir = Path(output).expanduser().resolve()
        if report_dir == root or root in report_dir.parents:
            pytest.fail("Report output cannot be inside the original FITS root")
        report_dir.mkdir(parents=True, exist_ok=True)
        report = {
            "classification": "USER_SUPPLIED_NATIVE_BRICK_NOT_REMOTE_AUTHENTICATED",
            "target": "UGC 2162",
            "source_key": key,
            "files_sha256_local": original_hashes,
            "full_brick_qa": {
                "status": audit["status"], "shape": audit["shape"],
                "fractions": audit["fractions"],
                "unit_evidence": audit["unit_evidence"],
            },
            "scene_status": native["original_fits_provenance_status"],
            "tractor_rows": contamination["stats"]["catalog_rows"],
            "sky_status": sky_result["status"],
            "aperture_status": aperture_result["status"],
            "mask_draft_invalidation_verified": True,
            "source_files_unchanged": True,
            "science_adopted": False,
        }
        (report_dir / "archive_candidate_e2e_receipt.json").write_text(
            json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
        widget.restore_navigation(stage="scene", band=band)
        QApplication.processEvents()
        assert explorer.grab().save(str(report_dir / "archive_candidate_scene.png"))
