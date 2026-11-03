"""P10b offline NED ObjectLookup and independent catalogue evidence contracts."""
import json
import math

import pytest

from lsbpipe.surveys.legacy import SkyField
from lsbpipe.surveys.object_projects import create_object_project, read_object_project
from lsbpipe.surveys.galaxy_knowledge import store_simbad_snapshot, confirm_simbad_identity
from lsbpipe.surveys.ned_knowledge import (
    API, query_ned_identity, save_ned_snapshot, read_ned_snapshot,
    latest_ned_snapshot, confirm_ned_identity, latest_ned_confirmation,
)
from lsbpipe.surveys.knowledge_summary import load_galaxy_knowledge
from lsbpipe.gui.galaxy_card_model import snapshot_from_record


@pytest.fixture
def leo(tmp_path):
    return create_object_project(
        "Leo P", tmp_path / "Leo_P",
        field=SkyField(155.438012, 18.088025, 4.))


def _known(name="Leo P", ra=155.4382, dec=18.0880):
    return {
        "Version": "2.0", "StatusCode": 100, "ResultCode": 3,
        "Supplied": name, "Interpreted": {"Name": name},
        "Aliases": ["Leo P", "AGC 208583"],
        "Preferred": {
            "Name": name,
            "Position": {"RA": ra, "Dec": dec,
                         "RefCode": "2022ApJ...000A..12B"},
            "ObjType": {"Value": "G", "RefCode": None},
            "Redshift": {"Value": 0.0009, "Uncertainty": 0.0001,
                         "RefCode": "2021AJ....000...1A",
                         "QualityFlag": None},
        },
    }


def _record(project):
    return {"root": str(project.root), "project_id": project.project_id,
            "display_name": project.display_name, "available": True,
            "legacy": False, "last_opened_at": ""}


def test_manual_ned_response_normalization_is_not_science_adoption():
    seen = []
    def fake(name):
        seen.append(name)
        return _known(name)
    got = query_ned_identity("Leo P", transport=fake)
    assert seen == ["Leo P"]
    assert got["api"] == API
    assert got["status"] == "KNOWN_OBJECT"
    assert got["candidate"]["object_type"] == "G"
    assert got["candidate"]["redshift"]["value"] == .0009
    assert got["candidate"]["redshift"]["frame"] == "heliocentric"
    assert got["candidate"]["redshift"]["adopted"] is False
    assert got["aliases_unadopted"] == ["Leo P", "AGC 208583"]
    assert got["source_response"]["Preferred"]["Name"] == "Leo P"
    assert got["science_status"] == "NOT_SCIENTIFICALLY_ADOPTED"


@pytest.mark.parametrize("code,expected", [
    (0, "NAME_NOT_INTERPRETED"),
    (1, "AMBIGUOUS_NAME"),
    (2, "NAME_NOT_IN_DATABASE"),
])
def test_unresolved_names_are_not_confirmable(leo, code, expected):
    raw = {"Version": "2.0", "StatusCode": 100, "ResultCode": code}
    if code == 1:
        raw["Interpreted"] = {"Aliases": ["LEO I", "LEO II"]}
    result = query_ned_identity("Leo", transport=lambda _: raw)
    assert result["status"] == expected
    assert result["candidate"] is None
    digest = save_ned_snapshot(leo.root, result)
    with pytest.raises(ValueError, match="only for known"):
        confirm_ned_identity(leo.root, digest)
    assert latest_ned_confirmation(leo.root) is None


@pytest.mark.parametrize("bad", [
    {"Version": "3.0", "StatusCode": 100, "ResultCode": 3},
    {"Version": "2.0", "StatusCode": 500, "ResultCode": 3},
    {"Version": "2.0", "StatusCode": 100, "ResultCode": 404},
    {"Version": "2.0", "StatusCode": 100, "ResultCode": 3,
     "Preferred": {"Name": "Bad", "Position": {"RA": float("nan"), "Dec": 1}}},
    {"Version": "2.0", "StatusCode": 100, "ResultCode": 3,
     "Preferred": {"Name": "Bad", "Position": {"RA": 360, "Dec": 1}}},
])
def test_malformed_ned_response_fails_closed(bad):
    with pytest.raises(ValueError):
        query_ned_identity("Leo P", transport=lambda _: bad)


def test_ned_consent_and_simbad_both_preserved_and_not_reconciled(leo):
    image = leo.root / "original-image-r.fits.fz"
    image.write_bytes(b"SCIENCE ORIGINAL")
    manifest = leo.root / "lsbpipe-object.json"
    before = manifest.read_bytes()
    lookup = query_ned_identity("Leo P", transport=lambda _: _known())
    d = save_ned_snapshot(leo.root, lookup)
    assert read_ned_snapshot(leo.root, d)["candidate"]["main_id"] == "Leo P"
    assert latest_ned_snapshot(leo.root)[0] == d
    assert load_galaxy_knowledge(leo.root).identities == ()
    nd = confirm_ned_identity(leo.root, d)
    assert len(nd) == 64
    assert latest_ned_confirmation(leo.root)["candidate"]["main_id"] == "Leo P"
    summary = load_galaxy_knowledge(leo.root)
    assert [x.provider for x in summary.identities] == ["NED"]
    assert summary.provider_offset_arcsec is None
    assert summary.confirmed_ids == ("Leo P",)
    assert summary.identities[0].redshift == pytest.approx(.0009)

    simbad = {
        "schema": "desidera.knowledge.simbad-query/1", "source": "SIMBAD",
        "fetched_at_utc": "2026-11-03T12:00:00+00:00",
        "query": {"name": "Leo P"},
        "candidates": [{"main_id": "NAME Leo P", "ra_deg": 155.438011,
                        "dec_deg": 18.088024, "separation_arcsec": .004,
                        "coordinate_bibcode": "", "matched_id": "",
                        "discovered_via": ["name"]}],
        "truncated": False, "status": "CANDIDATES_UNVERIFIED",
        "science_status": "NOT_SCIENTIFICALLY_ADOPTED",
    }
    sd = store_simbad_snapshot(leo.root, simbad)
    confirm_simbad_identity(leo.root, sd, "NAME Leo P")
    both = load_galaxy_knowledge(leo.root)
    assert [x.provider for x in both.identities] == ["SIMBAD", "NED"]
    assert both.provider_offset_arcsec > 0
    assert set(both.confirmed_ids) == {"NAME Leo P", "Leo P"}
    card = snapshot_from_record(_record(leo))
    assert card.matches("NAME Leo P")
    assert card.matches("Leo P")
    assert len(card.catalogue_identities) == 2
    assert read_object_project(leo.root).aliases == []
    assert manifest.read_bytes() == before
    assert image.read_bytes() == b"SCIENCE ORIGINAL"


def test_corrupt_ned_evidence_never_downgrades_to_no_match(leo):
    result = query_ned_identity("Leo P", transport=lambda _: _known())
    digest = save_ned_snapshot(leo.root, result)
    confirm_ned_identity(leo.root, digest)
    path = leo.root / "knowledge" / "ned" / f"snapshot-{digest}.json"
    path.write_text(path.read_text().replace("Leo P", "WRONG OBJECT"))
    with pytest.raises(ValueError, match="Corrupted"):
        load_galaxy_knowledge(leo.root)
    card = snapshot_from_record(_record(leo))
    assert card.knowledge_status == "EVIDENCE_ERROR"
    assert card.catalogue_identity is None


def test_cross_project_ned_receipt_refused(leo, tmp_path):
    other = create_object_project("Other", tmp_path / "other")
    result = query_ned_identity("Leo P", transport=lambda _: _known())
    digest = save_ned_snapshot(leo.root, result)
    other_dir = other.root / "knowledge" / "ned"
    other_dir.mkdir(parents=True)
    source = leo.root / "knowledge" / "ned" / f"snapshot-{digest}.json"
    (other_dir / source.name).write_bytes(source.read_bytes())
    with pytest.raises(ValueError, match="another galaxy"):
        read_ned_snapshot(other.root, digest)
    assert snapshot_from_record(_record(other)).knowledge_status == "EVIDENCE_ERROR"
