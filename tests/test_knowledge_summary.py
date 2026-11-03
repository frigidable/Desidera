"""P10a: SIMBAD confirmed identity becomes read-only galaxy passport data.

No network, data downloads, FITS mutation, or scientific adoption.
"""
from pathlib import Path

import pytest

from lsbpipe.surveys.legacy import SkyField
from lsbpipe.surveys.object_projects import (
    create_object_project, read_object_project, save_object_project,
)
from lsbpipe.surveys.galaxy_knowledge import (
    store_simbad_snapshot, confirm_simbad_identity,
)
from lsbpipe.surveys.knowledge_summary import (
    CatalogueIdentity, load_galaxy_knowledge, _distance_arcsec,
    NO_CONFIRMATION, CANDIDATES_ONLY, USER_CONFIRMED,
)
from lsbpipe.gui.galaxy_card_model import snapshot_from_record


@pytest.fixture
def leo(tmp_path):
    return create_object_project(
        "Leo P", tmp_path / "Leo_P",
        field=SkyField(155.438012, 18.088025, 4.0),
    )


def _snapshot():
    return {
        "schema": "desidera.knowledge.simbad-query/1",
        "source": "SIMBAD",
        "fetched_at_utc": "2026-11-03T10:00:00+00:00",
        "query": {"name": "Leo P", "ra_deg": 155.438012,
                  "dec_deg": 18.088025, "radius_arcsec": 60},
        "candidates": [
            {"main_id": "NAME Leo P", "ra_deg": 155.438011,
             "dec_deg": 18.088024, "separation_arcsec": 0.005,
             "coordinate_bibcode": "2020A&A...000A..01Z",
             "matched_id": "Leo P", "discovered_via": ["name"]},
            {"main_id": "Foreground star", "ra_deg": 155.4382,
             "dec_deg": 18.0891, "separation_arcsec": 4.0,
             "coordinate_bibcode": "", "matched_id": "",
             "discovered_via": ["cone"]},
        ],
        "truncated": False,
        "endpoint_evidence": {
            "provider": "SIMBAD", "endpoint_host": "simbad.harvard.edu",
            "endpoint_policy": "auto", "attempted_hosts": [
                "simbad.cds.unistra.fr", "simbad.harvard.edu"],
            "failed_attempts": [{
                "host": "simbad.cds.unistra.fr",
                "error_type": "TimeoutError"}],
            "fallback_used": True,
        },
        "status": "CANDIDATES_UNVERIFIED",
        "science_status": "NOT_SCIENTIFICALLY_ADOPTED",
    }


def _record(p):
    return {"project_id": p.project_id, "root": str(p.root),
            "display_name": p.display_name, "last_opened_at": "",
            "available": True, "legacy": False}


def test_no_catalogue_is_not_an_error_and_open_is_read_only(leo):
    manifest = (leo.root / "lsbpipe-object.json").read_bytes()
    summary = load_galaxy_knowledge(leo.root)
    assert summary.status == NO_CONFIRMATION
    assert summary.identities == ()
    card = snapshot_from_record(_record(leo))
    assert card.catalogue_identity is None
    assert card.knowledge_status == NO_CONFIRMATION
    assert (leo.root / "lsbpipe-object.json").read_bytes() == manifest
    assert not (leo.root / "knowledge").exists()


def test_unconfirmed_simbad_neighbours_never_become_project_aliases(leo):
    digest = store_simbad_snapshot(leo.root, _snapshot())
    assert len(digest) == 64
    summary = load_galaxy_knowledge(leo.root)
    assert summary.status == CANDIDATES_ONLY
    assert summary.cached_candidate_count == 2
    card = snapshot_from_record(_record(leo))
    assert not card.matches("Foreground star")
    assert not card.matches("NAME Leo P")
    assert read_object_project(leo.root).aliases == ()


def test_confirmed_identity_becomes_provenance_linked_card_without_mutation(leo):
    original = leo.root / "image-r.fits.fz"
    original.write_bytes(b"ORIGINAL FITS MUST NOT CHANGE")
    manifest = leo.root / "lsbpipe-object.json"
    before = manifest.read_bytes()
    digest = store_simbad_snapshot(leo.root, _snapshot())
    confirmation_sha = confirm_simbad_identity(leo.root, digest, "NAME Leo P")

    summary = load_galaxy_knowledge(leo.root)
    assert summary.status == USER_CONFIRMED
    assert summary.confirmed_ids == ("NAME Leo P",)
    assert len(summary.identities) == 1
    identity = summary.identities[0]
    assert identity.provider == "SIMBAD"
    assert identity.confirmation_sha256 == confirmation_sha
    assert identity.query_snapshot_sha256 == digest
    assert identity.endpoint_host == "simbad.harvard.edu"
    assert identity.coordinate_bibcode == "2020A&A...000A..01Z"
    assert identity.project_offset_arcsec is not None
    assert identity.project_offset_arcsec < 0.01

    card = snapshot_from_record(_record(leo))
    assert card.knowledge_status == USER_CONFIRMED
    assert card.catalogue_identity == identity
    assert card.matches("name leo p")
    assert not card.matches("Foreground star")
    assert "NAME Leo P" not in read_object_project(leo.root).aliases
    assert manifest.read_bytes() == before
    assert original.read_bytes() == b"ORIGINAL FITS MUST NOT CHANGE"


def test_project_coordinates_remain_authoritative_even_after_user_move(leo):
    digest = store_simbad_snapshot(leo.root, _snapshot())
    confirm_simbad_identity(leo.root, digest, "NAME Leo P")
    leo.field = SkyField(155.440000, 18.088025, 4.)
    save_object_project(leo)
    summary = load_galaxy_knowledge(leo.root)
    identity = summary.identities[0]
    assert 6 < identity.project_offset_arcsec < 8
    assert identity.ra_deg == 155.438011
    assert read_object_project(leo.root).field.ra_deg == 155.440000


def test_corrupted_saved_evidence_is_shown_as_error_not_unconfirmed(leo):
    digest = store_simbad_snapshot(leo.root, _snapshot())
    confirm_simbad_identity(leo.root, digest, "NAME Leo P")
    path = (leo.root / "knowledge" / "simbad" /
            f"snapshot-{digest}.json")
    path.write_text(path.read_text().replace(
        "NAME Leo P", "WRONG NAME"), encoding="utf-8")
    with pytest.raises(ValueError, match="Corrupt"):
        load_galaxy_knowledge(leo.root)
    card = snapshot_from_record(_record(leo))
    assert card.catalogue_identity is None
    assert card.knowledge_status == "EVIDENCE_ERROR"
    assert "Corrupt" in card.knowledge_error


def test_cross_project_receipt_rejected(leo, tmp_path):
    other = create_object_project("Other", tmp_path / "other")
    digest = store_simbad_snapshot(leo.root, _snapshot())
    confirmation = confirm_simbad_identity(leo.root, digest, "NAME Leo P")
    source = leo.root / "knowledge" / "simbad"
    dest = other.root / "knowledge" / "simbad"
    dest.mkdir(parents=True)
    for filename in (f"snapshot-{digest}.json",
                     f"confirmation-{confirmation}.json"):
        (dest / filename).write_bytes((source / filename).read_bytes())
    with pytest.raises(ValueError, match="another galaxy"):
        load_galaxy_knowledge(other.root)
    card = snapshot_from_record(_record(other))
    assert card.catalogue_identity is None
    assert card.knowledge_status == "EVIDENCE_ERROR"


def test_separation_wraps_ra_and_is_zero_for_identical_positions():
    assert _distance_arcsec(359.999, 0, 0.001, 0) == pytest.approx(7.2)
    assert _distance_arcsec(155.438, 18.088, 155.438, 18.088) == 0
