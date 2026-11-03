"""Offline Qt screenshot baseline for the *existing* Survey Explorer UI.

Captures geometry and layout only. The remote Legacy Viewer is deliberately
NOT loaded and is never represented as real survey imagery.
Requires: pip install -e ".[gui]" (Windows CI sets QT_QPA_PLATFORM=offscreen).
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from offline_gui_fixture import (offline_legacy_viewer, apply_offline_state,
                                 make_synthetic_fits_view)

# Intentional: this is a deterministic noninteractive CI capture, NOT a
# replacement for real Windows/on-device Qt WebEngine validation.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("QTWEBENGINE_CHROMIUM_FLAGS", "--disable-gpu")

from PySide6.QtCore import QPoint
from PySide6.QtGui import QFont, QFontInfo, QFontDatabase
from PySide6.QtWidgets import QApplication, QGroupBox, QTableWidget

from lsbpipe.gui.survey_explorer import SurveyExplorer
from lsbpipe.gui.project_hub import ProjectHub, CreateProjectDialog
from lsbpipe.surveys.object_projects import create_object_project
from lsbpipe.surveys.legacy import SkyField


SIZES = ((1366, 768), (1600, 900), (1920, 1080))


def _bounds(widget, root):
    pos = widget.mapTo(root, QPoint(0, 0))
    rect = widget.rect()
    return {
        "x": pos.x(), "y": pos.y(),
        "width": rect.width(), "height": rect.height(),
        "visible": bool(widget.isVisible()),
    }


def _capture(window, app, outdir: Path, width: int, height: int, state: str):
    window.resize(width, height)
    window.show()
    for _ in range(8):
        app.processEvents()

    screenshot = outdir / f"survey_{state}_{width}x{height}.png"
    pixmap = window.grab()
    if pixmap.isNull() or not pixmap.save(str(screenshot), "PNG"):
        raise RuntimeError(f"Unable to capture Qt screenshot: {screenshot}")
    if (pixmap.width(), pixmap.height()) != (width, height):
        raise AssertionError(
            f"Window expanded beyond requested size for {state}: "
            f"requested {width}x{height}, captured "
            f"{pixmap.width()}x{pixmap.height()}")
    if screenshot.stat().st_size < 1024:
        raise RuntimeError(f"Screenshot appears truncated: {screenshot}")

    splitter = window.centralWidget()
    if splitter is None:
        raise RuntimeError("Survey Explorer has no central QWidget")
    widgets = {
        "search_button": window.search_btn,
        "target_form": window.findChild(QGroupBox),
        "brick_table": window.bricks,
        "file_table": window.files,
        "coverage_button": window.qa_btn,
        "download_button": window.download_btn,
        "web_viewer_placeholder": window.sky,
        "science_qa_tab": window.science_qa,
        "science_qa_result": window.science_qa.result,
    }
    viewer_share = _bounds(splitter.widget(1), window)["width"] / width
    if viewer_share < 0.55:
        raise AssertionError(
            f"{state} viewer insufficient at {width}x{height}: "
            f"{viewer_share:.1%} (minimum 55%)")
    visible_table_height = 0
    if state == "S2":  # explicit advanced file table
        table_bounds = _bounds(window.files, window)
        viewport_bounds = _bounds(window.data_scroll.viewport(), window)
        visible_table_height = max(
            0,
            min(table_bounds["y"] + table_bounds["height"],
                viewport_bounds["y"] + viewport_bounds["height"]) -
            max(table_bounds["y"], viewport_bounds["y"]))
        if visible_table_height < 105:
            raise AssertionError(
                f"{state} file-table viewport too small at "
                f"{width}x{height}: {visible_table_height}px visible "
                "(need table header and several files)")
        primary = _bounds(window.download_btn, window)
        if (not primary["visible"] or primary["y"] < 0 or
                primary["y"] + primary["height"] > height):
            raise AssertionError(
                f"{state} Download selected action is not reachable "
                f"within the {width}x{height} window: {primary}")
    if state in ("S4", "S6"):
        if window.inspector_tabs.currentWidget() is not window.science_qa:
            raise AssertionError("Science QA tab not selected")
        if not window.science_qa.result_summary.isVisible():
            raise AssertionError("Current QA outcome summary is not reachable")
        # Extra diagnostics and full receipts are intentionally behind
        # explicit disclosure. S4/S6 capture the expanded scientific mode,
        # while S5 captures the uncluttered ordinary FITS inspection.
        if state in ("S4", "S6"):
            window.science_qa.advanced_toggle.setChecked(True)
            window.science_qa.result_toggle.setChecked(True)
            for _ in range(3):
                app.processEvents()
            if not window.science_qa.result.isVisible():
                raise AssertionError("Expanded QA report was not reachable")
            expected = "compare" if state == "S6" else "square"
            if not window.science_qa.buttons[expected].isVisible():
                raise AssertionError(
                    f"{state} expanded QA action {expected!r} unreachable")
            if state == "S6" and (
                    not window.science_qa.north_directory.isVisible() or
                    not window.science_qa.south_directory.isVisible()):
                raise AssertionError(
                    "Separate original North/South inputs are unreachable")
    if state == "S5":
        if window.data_drawer.isVisible():
            raise AssertionError(
                "FITS viewer should not be squeezed by archive controls")
        if window.preview_tabs.currentWidget() is not window.fits_explorer:
            raise AssertionError("S5 native FITS Explorer is not visible")
        if not window.fits_explorer.layer.isVisible():
            raise AssertionError("S5 image-layer control is clipped or hidden")
        if window.fits_explorer.base_item.pixmap().isNull():
            raise AssertionError("S5 synthetic image not rendered")
    report = {
        "state": state + "_offline_synthetic",
        "viewer_width_fraction": round(viewer_share, 4),
        "visible_file_table_height_px": visible_table_height,
        "capture_font_family": QFontInfo(app.font()).family(),
        "font_catalog_families": QFontDatabase.families()[:40],
        "size_requested_logical": [width, height],
        "screenshot": screenshot.name,
        "screenshot_px": [pixmap.width(), pixmap.height()],
        "screenshot_bytes": screenshot.stat().st_size,
        "github_sha": os.getenv("GITHUB_SHA"),
        "network_used": False,
        "viewer_content": "BLANK_OFFLINE_PLACEHOLDER_NOT_SURVEY_IMAGERY",
        "scientific_status": "NO_SCIENCE_RESULT",
        "science_workbench_visible": (
            window.inspector_tabs.currentWidget() is window.science_qa),
        "splitter": _bounds(splitter, window),
        "regions": {
            "control_panel": _bounds(splitter.widget(0), window),
            "viewer": _bounds(splitter.widget(1), window),
        },
        "widgets": {
            name: _bounds(widget, window)
            for name, widget in widgets.items() if widget is not None
        },
        "group_boxes": [
            {"title": box.title(), **_bounds(box, window)}
            for box in window.findChildren(QGroupBox)
        ],
        "tables": [
            {"columns": table.columnCount(), "rows": table.rowCount(),
             **_bounds(table, window)}
            for table in window.findChildren(QTableWidget)
        ],
    }
    (outdir / f"survey_{state}_{width}x{height}.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="gui-artifacts")
    args = parser.parse_args()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    app = QApplication.instance() or QApplication([])
    # The Windows offscreen QPA backend may report an EMPTY font database.
    # Register system fonts already installed on the CI machine (do not
    # redistribute fonts). This must precede any widget construction.
    windows_font_dir = Path(os.environ.get("WINDIR", "C:/Windows")) / "Fonts"
    loaded_fonts = []
    for filename in ("segoeui.ttf", "arial.ttf", "tahoma.ttf"):
        font_path = windows_font_dir / filename
        if font_path.is_file():
            font_id = QFontDatabase.addApplicationFont(str(font_path))
            if font_id >= 0:
                loaded_fonts.extend(QFontDatabase.applicationFontFamilies(font_id))
    fonts = set(QFontDatabase.families())
    preferred = next((name for name in (
        "Segoe UI", "Arial", "Tahoma", "DejaVu Sans")
        if name in fonts), app.font().family())
    app.setFont(QFont(preferred, 10))
    print(f"Qt capture font database: {len(fonts)} families; "
          f"registered: {loaded_fonts}; selected: {preferred!r}")
    # No Viewer navigation or archive activity from this controlled
    # capture process. The application source itself is NOT modified.
    reports = []
    # Isolate launcher registry even when this script runs outside pytest:
    # GUI evidence must NEVER alter the real user's recent-project index.
    os.environ["LSBPIPE_SESSION_DIR"] = str(outdir.resolve() / ".session")
    from lsbpipe.surveys.object_projects import library_index_path
    fixture_root = outdir.resolve() / "test-project"
    with offline_legacy_viewer():
        hub = ProjectHub()
        try:
            for w, h in SIZES:
                hub.resize(w, h)
                hub.show()
                app.processEvents()
                path = outdir / f"hub_empty_{w}x{h}.png"
                if not hub.grab().save(str(path), "PNG"):
                    raise RuntimeError(f"Could not capture {path}")
                reports.append({
                    "state": "H0_empty_library_offline", "screenshot": path.name,
                    "screenshot_bytes": path.stat().st_size,
                    "network_used": False, "scientific_status": "NO_SCIENCE_RESULT"})
            create_object_project(
                "UGC 2162 · synthetic library entry", fixture_root,
                field=SkyField(40.096208, 1.229250, 5.))
            hub.refresh()
            for w, h in SIZES:
                hub.resize(w, h)
                hub.show()
                app.processEvents()
                path = outdir / f"hub_recent_{w}x{h}.png"
                if not hub.grab().save(str(path), "PNG"):
                    raise RuntimeError(f"Could not capture {path}")
                reports.append({
                    "state": "H1_recent_project_offline",
                    "screenshot": path.name,
                    "screenshot_bytes": path.stat().st_size,
                    "network_used": False, "scientific_status": "NO_SCIENCE_RESULT"})
        finally:
            hub.close()
        # Capture the real single-page creation form with synthetic input.
        # Do not press Create: screenshot capture may never mutate science data.
        wizard = CreateProjectDialog(
            default_root=outdir.resolve() / "fixture-new-projects")
        try:
            wizard.name.setText("Synthetic LSB galaxy")
            wizard.ra.setText("40.096208")
            wizard.dec.setText("+1.229250")
            wizard.show()
            for _ in range(4):
                app.processEvents()
            path = outdir / "hub_new_project_form.png"
            shot = wizard.grab()
            if shot.isNull() or not shot.save(str(path), "PNG"):
                raise RuntimeError(f"Could not capture {path}")
            reports.append({
                "state": "H2_new_project_form_offline",
                "screenshot": path.name,
                "screenshot_px": [shot.width(), shot.height()],
                "screenshot_bytes": path.stat().st_size,
                "network_used": False,
                "scientific_status": "NO_SCIENCE_RESULT",
                "project_created": False,
            })
        finally:
            wizard.close()
        # P0.5b Data Manager: synthetic native WCS headers only.
        # Never show this screenshot as evidence of real Legacy/HSC coverage.
        from astropy.io import fits
        from astropy.wcs import WCS
        import numpy as np
        from lsbpipe.gui.dataset_manager import DatasetManager
        from lsbpipe.surveys.dataset_registry import (
            register_external_fits, read_dataset_index)
        from lsbpipe.surveys.dataset_coverage import survey_coverage_footprints
        from lsbpipe.surveys.object_projects import read_object_project
        synthetic_project = read_object_project(fixture_root)
        for n, (provider, band, dx) in enumerate((
                ("SYNTHETIC HSC", "g", -45.),
                ("SYNTHETIC CFHT", "r", 50.))):
            wcs = WCS(naxis=2)
            wcs.wcs.crpix = [36., 36.]
            wcs.wcs.cdelt = [-2.5/3600, 2.5/3600]
            wcs.wcs.crval = [40.096208+dx/3600, 1.22925]
            wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
            file = outdir / "synthetic-original-headers" / f"{n}-{band}.fits"
            file.parent.mkdir(exist_ok=True, parents=True)
            fits.PrimaryHDU(
                data=np.ones((72, 72), dtype=np.float32),
                header=wcs.to_header()).writeto(file, overwrite=True)
            register_external_fits(
                synthetic_project.root, file,
                provider=provider, band=band)
        manager = DatasetManager()
        try:
            manager.bind_project(
                synthetic_project, read_dataset_index(synthetic_project.root))
            manager.set_coverage(
                survey_coverage_footprints(synthetic_project.root))
            for width, height in ((1366, 768), (1600, 900)):
                manager.resize(width, height)
                manager.show()
                for _ in range(8):
                    app.processEvents()
                path = outdir / f"dataset_manager_SYNTHETIC_{width}x{height}.png"
                shot = manager.grab()
                if shot.isNull() or not shot.save(str(path), "PNG"):
                    raise RuntimeError(f"Unable to save {path}")
                reports.append({
                    "state": "DM1_synthetic_external_wcs_headers",
                    "screenshot": path.name,
                    "screenshot_px": [shot.width(), shot.height()],
                    "screenshot_bytes": path.stat().st_size,
                    "network_used": False,
                    "viewer_content": "SYNTHETIC_NOT_SURVEY_COVERAGE",
                    "scientific_status": "NO_SCIENCE_RESULT",
                })
        finally:
            manager.close()
        assert library_index_path().is_file()
    # Fresh widget instance for every state: no inadvertent selection/resize
    # carryover between screenshots, but all components are real Qt widgets.
    with offline_legacy_viewer():
        for state in ("S0", "S1", "S2", "S3", "S4", "S5", "S6"):
            window = SurveyExplorer()
            try:
                # Exercise the populated-state transition AFTER Qt has shown
                # the window, like a real user search (not before layout).
                window.show()
                for _ in range(4):
                    app.processEvents()
                apply_offline_state(window, state)
                for _ in range(4):
                    app.processEvents()
                for w, h in SIZES:
                    reports.append(_capture(
                        window, app, outdir, w, h, state))
            finally:
                window.close()
                app.processEvents()
    # Synthetic three-band comparison at a realistic desktop viewport.
    # This fixture contains NO DESI pixels; it is strictly geometry evidence.
    with offline_legacy_viewer():
        window = SurveyExplorer(restore_session=False, launch_remote=False)
        try:
            images = {}
            for band in ("g", "r", "z"):
                fake = make_synthetic_fits_view()
                fake["band"] = band
                fake["source"] = f"DR9 north · 1544p465 · {band}"
                fake["data"]["image"] = (
                    fake["data"]["image"] * {"g": 0.75, "r": 1., "z": 1.3}[band])
                images[band] = fake
            window._multiband_ready({
                "status": "DIAGNOSTIC_NOT_ADOPTED",
                "datasets": images,
                "missing": {"i": "Synthetic missing band"},
            })
            window.preview_label.setText(
                "SYNTHETIC VISUAL FIXTURE · NO ORIGINAL SURVEY FITS")
            window.resize(1600, 900)
            window.show()
            import time
            for _ in range(500):
                app.processEvents()
                if (window.multiband_comparison._future is None and
                        all(not window.multiband_comparison.images[b]
                                .pixmap().isNull() for b in images)):
                    break
                time.sleep(0.01)
            if window.multiband_comparison._future is not None:
                raise AssertionError("Synthetic multi-band previews did not render")
            path = outdir / "multiband_synthetic_1600x900.png"
            shot = window.grab()
            if not shot.save(str(path), "PNG"):
                raise RuntimeError(f"Could not capture {path}")
            reports.append({
                "state": "S7_multiband_synthetic_not_original",
                "screenshot": path.name,
                "screenshot_px": [shot.width(), shot.height()],
                "screenshot_bytes": path.stat().st_size,
                "viewer_content": "SYNTHETIC_NOT_SURVEY_FITS",
                "network_used": False,
                "scientific_status": "NO_SCIENCE_RESULT",
            })
        finally:
            window.close()
            app.processEvents()
    index = {
        "evidence_level": "E1_CAPTURED_NOT_VISUALLY_REVIEWED",
        "screens": reports,
        "limitations": [
            "WebEngine intentionally replaced by a blank local offline placeholder",
            "All catalog and product S1–S6 data are synthetic fixtures",
            "S4–S6 are mock UI content; no original FITS are opened",
            "S5 FITS raster is artificial and cannot support scientific conclusions",
            "S7 band comparison is synthetic and contains no original survey pixels",
            "Offscreen Qt capture cannot verify real monitor / DPI rendering",
            "This job does not certify accessibility or scientific validity",
        ],
    }
    (outdir / "manifest.json").write_text(
        json.dumps(index, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"Offline GUI capture: {len(reports)} screenshots in {outdir}")
    for report in reports:
        print(report["screenshot"], report["screenshot_bytes"], "bytes")
    # Explicit Qt destruction matters on Windows offscreen runner. Returning
    # while still retaining closed WebEngine/graphics widgets can provoke
    # interpreter-shutdown errors AFTER all 31 files have been written.
    import gc
    for widget in (window, manager, wizard, hub):
        widget.close()
        widget.deleteLater()
    app.closeAllWindows()
    app.processEvents()
    del window, manager, wizard, hub
    gc.collect()
    app.processEvents()
    print("Qt capture cleanup complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
