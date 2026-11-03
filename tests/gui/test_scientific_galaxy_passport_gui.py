"""P10g Qt acceptance for meaningful, evidence-bound galaxy passports."""
import pytest

pytest.importorskip("pytestqt")
pytest.importorskip("PySide6.QtWebEngineWidgets")

from PySide6.QtCore import Qt

from scripts.offline_gui_fixture import offline_legacy_viewer
from lsbpipe.gui.project_hub import ProjectHub
from lsbpipe.gui.galaxy_card_model import snapshot_from_record
from lsbpipe.surveys.legacy import SkyField
from lsbpipe.surveys.object_projects import (
    create_object_project, read_object_project, save_object_project,
)
from lsbpipe.surveys.galaxy_knowledge import (
    store_simbad_snapshot, confirm_simbad_identity,
)
from lsbpipe.surveys.ned_knowledge import (
    query_ned_identity, save_ned_snapshot, confirm_ned_identity,
)


def _hub(qtbot):
    with offline_legacy_viewer():
        window = ProjectHub()
    qtbot.addWidget(window)
    window.resize(1300, 800)
    window.show()
    return window


def _record(project):
    return {"project_id": project.project_id, "display_name": project.display_name,
            "root": str(project.root), "available": True, "legacy": False,
            "last_opened_at": ""}


def test_no_dataset_is_normal_state_not_error_message(qtbot, tmp_path):
    project = create_object_project(
        "dw 1343+58", tmp_path / "DW1343",
        field=SkyField(205.779167, 58.227778, 4))
    original = (project.root / "lsbpipe-object.json").read_bytes()
    hub = _hub(qtbot)
    card = hub.galaxy_card
    assert card.snapshot.project_id == project.project_id
    assert card.snapshot.readiness.code == "DATASET_NOT_SELECTED"
    assert "No survey dataset selected" in card.fits.text()
    assert card.open_button.text() == "Open project"
    assert "find data" in card.open_button.toolTip()
    assert "Why identify this galaxy?" in card.catalogue_help.text()
    assert card.catalogue_help.isVisible()
    assert not card.literature.isVisible()
    assert card.preview.text() == "No image yet"
    assert "Traceback" not in card.preview.text()
    assert not hub.local_preview_action.isEnabled()
    assert "Select survey data" in hub.projects.item(0).text()

    # A stale previous task failure must not fill the image square if there
    # is no original IMAGE source to preview.
    hub._preview_errors[project.project_id] = (
        "FileNotFoundError: No original Legacy IMAGE FITS")
    hub._selection_changed()
    assert card.preview.text() == "No image yet"
    assert project.project_id not in hub._preview_errors
    hub._regenerate_preview()
    assert "FileNotFoundError" not in card.preview.text()
    assert (project.root / "lsbpipe-object.json").read_bytes() == original


def _ned_response():
    return {
        "Version": "2.1", "StatusCode": 100, "ResultCode": 3,
        "Supplied": "Leo P", "Interpreted": {"Name": "Leo P"},
        "Aliases": ["Leo P", "AGC 208583"],
        "Preferred": {
            "Name": "PGC1 5065058 NED001",
            "Position": {"RA": 155.43708330, "Dec": 18.08944444,
                         "RefCode": "2018ApJ...861..49H"},
            "ObjType": {"Value": "G", "RefCode": None},
            "Redshift": {"Value": 0.0008806092,
                         "Uncertainty": None,
                         "RefCode": "2026ApJ...997..271M",
                         "QualityFlag": None},
        },
    }


def _simbad_response():
    return {
        "schema": "desidera.knowledge.simbad-query/1",
        "source": "SIMBAD",
        "fetched_at_utc": "2026-11-03T12:00:00+00:00",
        "query": {"name": "Leo P"},
        "candidates": [
            {"main_id": "NAME Leo P", "ra_deg": 155.438011,
             "dec_deg": 18.088024, "separation_arcsec": 0.005,
             "matched_id": "Leo P", "coordinate_bibcode": "",
             "discovered_via": ["name"]}
        ],
        "truncated": False, "status": "CANDIDATES_UNVERIFIED",
        "science_status": "NOT_SCIENTIFICALLY_ADOPTED",
    }


def test_confirmed_catalogues_supply_scoped_literature_not_adopted_science(
        qtbot, tmp_path):
    project = create_object_project(
        "Leo P", tmp_path / "LeoP",
        field=SkyField(155.438012, 18.088025, 4))
    project.selected_dataset = ("dr11", "south", "1554p180", "r")
    save_object_project(project)
    saved = (project.root / "lsbpipe-object.json").read_bytes()

    simbad = store_simbad_snapshot(project.root, _simbad_response())
    confirm_simbad_identity(project.root, simbad, "NAME Leo P")
    result = query_ned_identity("Leo P", transport=lambda _: _ned_response())
    sha = save_ned_snapshot(project.root, result)
    confirm_ned_identity(project.root, sha)
    hub = _hub(qtbot)
    card = hub.galaxy_card

    assert "NAME Leo P" in card.catalogue.text()
    assert "PGC1 5065058 NED001" in card.catalogue.text()
    assert card.catalogue.text().count('<a href=') == 2
    assert "Type (NED)" in card.literature.text()
    assert "G · galaxy" in card.literature.text()
    assert "0.0008806092" in card.literature.text()
    assert "2026ApJ...997..271M" not in card.literature.text()
    assert "2026ApJ...997..271M" in card.catalogue_provenance.text()
    assert "not adopted" in card.literature.text()
    assert "6.02" in card.cross_catalogue.text()
    assert card.literature.isVisible()
    assert not card.catalogue_help.isVisible()
    assert "Survey selected" in card.fits.text()
    assert card.open_button.text() == "Open project"
    assert "get original FITS" in card.open_button.toolTip()
    assert "No image yet" in card.preview.text()
    assert card.details_panel.isVisible() is False
    assert not card.cross_catalogue.isVisible()
    assert not card.observation_explanation.isVisible()
    card.details_toggle.click()
    assert card.cross_catalogue.isVisible()
    assert card.observation_explanation.isVisible()
    assert "Separation from project centre" in (
        card.catalogue_provenance.text())
    assert read_object_project(project.root).field == project.field
    assert (project.root / "lsbpipe-object.json").read_bytes() == saved


def test_stored_but_unconfirmed_candidates_never_become_literature(qtbot, tmp_path):
    project = create_object_project(
        "Leo P", tmp_path / "LeoP",
        field=SkyField(155.438012, 18.088025, 4))
    store_simbad_snapshot(project.root, _simbad_response())
    hub = _hub(qtbot)
    card = hub.galaxy_card
    assert "no identity confirmed" in card.catalogue.text()
    assert "Why identify" in card.catalogue_help.text()
    assert not card.literature.isVisible()
    assert "<a" not in card.catalogue.text()
