"""P10f integration: independent SIMBAD/NED links in one evidence-safe passport."""
import pytest

pytest.importorskip("pytestqt")
pytest.importorskip("PySide6.QtWebEngineWidgets")

from scripts.offline_gui_fixture import offline_legacy_viewer
from lsbpipe.gui.project_hub import ProjectHub
from lsbpipe.surveys.galaxy_knowledge import (
    store_simbad_snapshot, confirm_simbad_identity,
)
from lsbpipe.surveys.ned_knowledge import (
    save_ned_snapshot, confirm_ned_identity,
)
from lsbpipe.surveys.object_projects import create_object_project, read_object_project
from lsbpipe.surveys.legacy import SkyField


def _simbad_snapshot():
    return {
        "schema": "desidera.knowledge.simbad-query/1",
        "source": "SIMBAD",
        "fetched_at_utc": "2026-11-03T12:00:00+00:00",
        "query": {"name": "Leo P"},
        "candidates": [
            {"main_id": "NAME Leo P", "ra_deg": 155.438011,
             "dec_deg": 18.088024, "separation_arcsec": 0.005,
             "matched_id": "Leo P", "coordinate_bibcode": "",
             "discovered_via": ["name"]},
            {"main_id": "unconfirmed nearby source", "ra_deg": 155.4383,
             "dec_deg": 18.089, "separation_arcsec": 4.2,
             "matched_id": "", "coordinate_bibcode": "",
             "discovered_via": ["cone"]},
        ],
        "truncated": False, "status": "CANDIDATES_UNVERIFIED",
        "science_status": "NOT_SCIENTIFICALLY_ADOPTED",
        "endpoint_evidence": {
            "provider": "SIMBAD", "endpoint_host": "simbad.harvard.edu",
            "endpoint_policy": "auto", "attempted_hosts": [
                "simbad.cds.unistra.fr", "simbad.harvard.edu"],
            "failed_attempts": [], "fallback_used": True,
        },
    }


def _ned_snapshot():
    return {
        "schema": "desidera.knowledge.ned-objectlookup/1",
        "provider": "NED", "api": "https://ned.ipac.caltech.edu/srs/ObjectLookup",
        "fetched_at_utc": "2026-11-03T14:00:00+00:00",
        "query_name": "Leo P", "api_version": "2.1",
        "result_code": 3, "status": "KNOWN_OBJECT",
        "candidate": {
            "main_id": "PGC1 5065058 NED001", "ra_deg": 155.4370833,
            "dec_deg": 18.0894444, "object_type": "G",
            "position_refcode": "2018ApJ...861..49H",
            "redshift": {"value": 0.0008806092, "uncertainty": None,
                         "refcode": "2026ApJ...997..271M",
                         "quality_flag": "", "frame": "heliocentric",
                         "adopted": False},
        },
        "aliases_unadopted": ["Leo P", "AGC 208583"],
        "ambiguous_names": [],
        "source_response": {"Version": "2.1", "StatusCode": 100,
                            "ResultCode": 3},
        "science_status": "NOT_SCIENTIFICALLY_ADOPTED",
    }


def test_both_confirmed_names_openable_without_changing_project_or_science(
        qtbot, tmp_path):
    p = create_object_project(
        "Leo P", tmp_path / "Leo_P",
        field=SkyField(155.438012, 18.088025, 4.))
    manifest = p.root / "lsbpipe-object.json"
    before = manifest.read_bytes()
    d = store_simbad_snapshot(p.root, _simbad_snapshot())
    # Search candidates are not catalogue associations.
    with offline_legacy_viewer():
        hub = ProjectHub()
    qtbot.addWidget(hub)
    hub.show()
    assert "<a " not in hub.galaxy_card.catalogue.text()

    confirm_simbad_identity(p.root, d, "NAME Leo P")
    ned_digest = save_ned_snapshot(p.root, _ned_snapshot())
    confirm_ned_identity(p.root, ned_digest)
    hub.refresh()
    displayed = hub.galaxy_card.catalogue.text()
    assert displayed.count("<a href=") == 2
    assert 'simbad.harvard.edu/simbad/sim-id?Ident=NAME+Leo+P' in displayed
    assert 'ned.ipac.caltech.edu/byname?objname=PGC1+5065058+NED001' in displayed
    assert "unconfirmed nearby source" not in displayed
    assert hub.galaxy_card.catalogue.openExternalLinks()
    assert hub.galaxy_card.snapshot.provider_offset_arcsec is not None
    assert "SIMBAD–NED positions differ" in (
        hub.galaxy_card.catalogue_provenance.text())
    assert read_object_project(p.root).field == p.field
    assert manifest.read_bytes() == before

    # Read-only/offline reopen must retain same links.
    with offline_legacy_viewer():
        reopened = ProjectHub()
    qtbot.addWidget(reopened)
    assert reopened.galaxy_card.catalogue.text() == displayed
    assert manifest.read_bytes() == before


def test_tampered_ned_receipt_never_produces_trusted_href(qtbot, tmp_path):
    p = create_object_project(
        "Leo P", tmp_path / "Leo_P",
        field=SkyField(155.438012, 18.088025, 4.))
    digest = save_ned_snapshot(p.root, _ned_snapshot())
    confirm_ned_identity(p.root, digest)
    original = p.root / "knowledge" / "ned" / f"snapshot-{digest}.json"
    original.write_text(original.read_text().replace(
        "PGC1 5065058 NED001", "FAKE OBJECT"), encoding="utf-8")
    with offline_legacy_viewer():
        hub = ProjectHub()
    qtbot.addWidget(hub)
    assert hub.galaxy_card.snapshot.knowledge_status == "EVIDENCE_ERROR"
    assert "<a " not in hub.galaxy_card.catalogue.text()
    assert hub.galaxy_card.details_panel.isVisible() is False  # not shown until hub.show()
    hub.show()
    assert hub.galaxy_card.details_toggle.isChecked()
