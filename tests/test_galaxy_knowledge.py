"""Offline P0.7a SIMBAD identity and immutable evidence contracts."""
import json

import pytest
from astropy.table import Table

from lsbpipe.surveys.legacy import SkyField
from lsbpipe.surveys.object_projects import create_object_project
from lsbpipe.surveys.galaxy_knowledge import (
    query_simbad_candidates, store_simbad_snapshot, read_simbad_snapshot,
    latest_simbad_snapshot, confirm_simbad_identity,
    latest_simbad_confirmation,
)


class FakeSimbad:
    ROW_LIMIT = None

    def query_object(self, name, wildcard=False):
        assert name == "UGC 2162"
        assert not wildcard
        return Table(rows=[
            ("UGC 2162", 40.096208, 1.22925, "2011A&A...000A..01X"),
        ], names=("main_id", "ra", "dec", "coo_bibcode"))

    def query_region(self, centre, radius):
        assert radius.to_value("arcsec") == 60.
        return Table(rows=[
            ("UGC 2162", 40.096208, 1.22925, "2011A&A...000A..01X"),
            ("Background star", 40.0963, 1.2294, ""),
        ], names=("main_id", "ra", "dec", "coo_bibcode"))


@pytest.fixture
def project(tmp_path):
    return create_object_project("UGC 2162", tmp_path / "UGC_2162",
                                 field=SkyField(40.096208, 1.22925, 5.))


def test_fake_simbad_preserves_ambiguous_candidates(project):
    value = query_simbad_candidates(
        name=project.display_name, ra_deg=project.field.ra_deg,
        dec_deg=project.field.dec_deg, client=FakeSimbad())
    assert len(value["candidates"]) == 2
    assert value["candidates"][0]["main_id"] == "UGC 2162"
    assert set(value["candidates"][0]["discovered_via"]) == {"name", "cone"}
    assert value["status"] == "CANDIDATES_UNVERIFIED"
    assert latest_simbad_confirmation(project.root) is None


def test_snapshot_and_explicit_confirm_never_mutate_originals(project):
    original = project.root / "original-r.fits.fz"
    original.write_bytes(b"EXACT ORIGINAL OBSERVATION")
    manifest = project.root / "lsbpipe-object.json"
    before = manifest.read_bytes()
    query = query_simbad_candidates(
        name=project.display_name, ra_deg=project.field.ra_deg,
        dec_deg=project.field.dec_deg, client=FakeSimbad())
    digest = store_simbad_snapshot(project.root, query)
    assert store_simbad_snapshot(project.root, query) == digest
    restored = read_simbad_snapshot(project.root, digest)
    assert restored["project_id"] == project.project_id
    assert latest_simbad_snapshot(project.root)[0] == digest
    with pytest.raises(ValueError, match="not unique"):
        confirm_simbad_identity(project.root, digest, "NOT A CANDIDATE")
    confirm_simbad_identity(project.root, digest, "UGC 2162")
    latest = latest_simbad_confirmation(project.root)
    assert latest["candidate"]["main_id"] == "UGC 2162"
    assert latest["coordinate_adoption"] is False
    assert latest["science_status"] == "NOT_SCIENTIFICALLY_ADOPTED"
    assert manifest.read_bytes() == before
    assert original.read_bytes() == b"EXACT ORIGINAL OBSERVATION"


def test_corrupt_snapshot_fails_closed(project):
    payload = query_simbad_candidates(
        name=project.display_name, ra_deg=40.096208, dec_deg=1.22925,
        client=FakeSimbad())
    digest = store_simbad_snapshot(project.root, payload)
    file = project.root / "knowledge" / "simbad" / f"snapshot-{digest}.json"
    data = json.loads(file.read_text())
    data["candidates"][0]["main_id"] = "WRONG GALAXY"
    file.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="Corrupt or modified"):
        read_simbad_snapshot(project.root, digest)
    with pytest.raises(ValueError):
        latest_simbad_snapshot(project.root)


@pytest.mark.parametrize("ra,dec,radius", [
    (float("nan"), 1.2, 60), (40., None, 60), (40., 1.2, 0),
    (360., 1.2, 60), (40., -92., 60),
])
def test_invalid_query_never_contacts_simbad(ra, dec, radius):
    with pytest.raises(ValueError):
        query_simbad_candidates(name="UGC 2162", ra_deg=ra, dec_deg=dec,
                                radius_arcsec=radius, client=FakeSimbad())


def test_cross_project_knowledge_receipt_rejected(project, tmp_path):
    p2 = create_object_project("Galaxy B", tmp_path / "Galaxy_B")
    query = query_simbad_candidates(
        name="UGC 2162", ra_deg=40.096208, dec_deg=1.22925,
        client=FakeSimbad())
    digest = store_simbad_snapshot(project.root, query)
    directory = p2.root / "knowledge" / "simbad"
    directory.mkdir(parents=True)
    (directory / f"snapshot-{digest}.json").write_bytes(
        (project.root / "knowledge" / "simbad" /
         f"snapshot-{digest}.json").read_bytes())
    with pytest.raises(ValueError, match="another galaxy"):
        read_simbad_snapshot(p2.root, digest)
