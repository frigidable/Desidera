"""P10a: production Galaxy Projects shows saved SIMBAD identity offline."""
import pytest

pytest.importorskip("pytestqt")
pytest.importorskip("PySide6.QtWebEngineWidgets")

from PySide6.QtCore import Qt

from scripts.offline_gui_fixture import offline_legacy_viewer
from lsbpipe.gui.project_hub import ProjectHub
from lsbpipe.surveys.object_projects import create_object_project, read_object_project
from lsbpipe.surveys.legacy import SkyField
from lsbpipe.surveys.galaxy_knowledge import (
    store_simbad_snapshot, confirm_simbad_identity,
)


def _evidence():
    return {
        "schema": "desidera.knowledge.simbad-query/1",
        "source": "SIMBAD",
        "fetched_at_utc": "2026-11-03T12:00:00+00:00",
        "query": {"name": "Leo P"},
        "candidates": [
            {"main_id": "NAME Leo P", "ra_deg": 155.438011,
             "dec_deg": 18.088024, "separation_arcsec": 0.,
             "coordinate_bibcode": "", "matched_id": "Leo P",
             "discovered_via": ["name"]},
            {"main_id": "Nearby field star", "ra_deg": 155.44,
             "dec_deg": 18.089, "separation_arcsec": 3.,
             "coordinate_bibcode": "", "matched_id": "",
             "discovered_via": ["cone"]},
        ],
        "truncated": False, "status": "CANDIDATES_UNVERIFIED",
        "science_status": "NOT_SCIENTIFICALLY_ADOPTED",
        "endpoint_evidence": {
            "provider": "SIMBAD", "endpoint_host": "simbad.harvard.edu",
            "endpoint_policy": "auto",
            "attempted_hosts": ["simbad.cds.unistra.fr", "simbad.harvard.edu"],
            "failed_attempts": [{"host": "simbad.cds.unistra.fr",
                                 "error_type": "TimeoutError"}],
            "fallback_used": True,
        },
    }


def test_passport_identity_search_and_offline_reopen(qtbot, tmp_path):
    p = create_object_project(
        "Leo P", tmp_path / "Leo_P",
        field=SkyField(155.438012, 18.088025, 4.))
    original = p.root / "original-image-r.fits.fz"
    original.write_bytes(b"NO FITS READ OR WRITE")
    manifest = p.root / "lsbpipe-object.json"
    before = manifest.read_bytes()
    with offline_legacy_viewer():
        hub = ProjectHub()
    qtbot.addWidget(hub)
    hub.show()
    assert "Not yet linked to SIMBAD or NED" in hub.galaxy_card.catalogue.text()
    assert "Why identify this galaxy?" in hub.galaxy_card.catalogue_help.text()

    d = store_simbad_snapshot(p.root, _evidence())
    assert "No user-confirmed" in hub.galaxy_card.catalogue.text()
    confirm_simbad_identity(p.root, d, "NAME Leo P")
    hub.refresh()
    assert "NAME Leo P" in hub.galaxy_card.catalogue.text()
    assert 'href="https://simbad.harvard.edu/simbad/sim-id?Ident=NAME+Leo+P"' in (
        hub.galaxy_card.catalogue.text())
    assert hub.galaxy_card.catalogue.openExternalLinks()
    assert "user-confirmed" in hub.galaxy_card.catalogue.toolTip()
    assert "Project coordinates unchanged" in (
        hub.galaxy_card.catalogue_provenance.text())
    assert "Separation from project centre" in (
        hub.galaxy_card.catalogue_provenance.text())
    assert "simbad.harvard.edu" not in hub.galaxy_card.name.text()
    assert read_object_project(p.root).aliases == []

    hub.filter.setText("NAME Leo P")
    assert hub.projects.count() == 1
    assert not hub.projects.item(0).isHidden()
    assert hub.galaxy_card.snapshot.project_id == p.project_id
    hub.filter.setText("Nearby field star")
    assert hub.projects.item(0).isHidden()
    assert hub.galaxy_card.snapshot is None

    # Closing and reopening the application reads only local evidence.
    with offline_legacy_viewer():
        reopened = ProjectHub()
    qtbot.addWidget(reopened)
    assert "NAME Leo P" in reopened.galaxy_card.catalogue.text()
    assert manifest.read_bytes() == before
    assert original.read_bytes() == b"NO FITS READ OR WRITE"


def test_corrupt_knowledge_displays_an_error_not_a_confirmed_match(
        qtbot, tmp_path):
    p = create_object_project(
        "Leo P", tmp_path / "Leo_P",
        field=SkyField(155.438012, 18.088025, 4.))
    d = store_simbad_snapshot(p.root, _evidence())
    confirm_simbad_identity(p.root, d, "NAME Leo P")
    bad = p.root / "knowledge" / "simbad" / f"snapshot-{d}.json"
    bad.write_text(bad.read_text().replace("NAME Leo P", "FALSIFIED"), "utf-8")

    with offline_legacy_viewer():
        hub = ProjectHub()
    qtbot.addWidget(hub)
    hub.show()
    assert "requires inspection" in hub.galaxy_card.catalogue.text()
    assert "<a " not in hub.galaxy_card.catalogue.text()
    assert "Corrupt" in hub.galaxy_card.catalogue_provenance.text()
    hub.filter.setText("NAME Leo P")
    assert hub.galaxy_card.snapshot is None
