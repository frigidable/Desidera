"""Offline Galaxy Card model contracts: no downloads, no manifest mutations."""
from pathlib import Path
from types import SimpleNamespace

from lsbpipe.gui.galaxy_card_model import (
    snapshot_from_record, snapshots_from_records,
)


def record(root: Path, *, name="Leo P", available=True, legacy=False):
    return {
        "project_id": "id-leo-p",
        "display_name": name,
        "root": str(root),
        "available": available,
        "legacy": legacy,
        "last_opened_at": "2026-11-02T10:00:00Z",
    }


def project(root, **overrides):
    data = dict(
        display_name="Leo P",
        field=SimpleNamespace(ra_deg=155.438012, dec_deg=18.088025),
        selected_dataset=("dr11", "south", "1554p180", "r"),
        data_root=None, aliases=["AGC 208583"],
    )
    data.update(overrides)
    return SimpleNamespace(**data)


def should_not_read(_):
    raise AssertionError("Missing records must not trigger file reads")


def test_card_uses_project_manifest_not_catalogue_guess(tmp_path):
    read = lambda _: project(tmp_path)
    card = snapshot_from_record(record(tmp_path), read_project=read)
    assert card.display_name == "Leo P"
    assert "155.438012" in card.coordinates
    assert "dr11 / south" in card.dataset
    assert card.fits_status == "No local FITS link"
    assert card.matches("agc 208583")
    assert card.matches("155.438")
    assert not card.matches("IC 10")
    assert card.thumbnail is None


def test_linked_directory_does_not_claim_integrity_or_scientific_qa(tmp_path):
    data_root = tmp_path / "data"
    data_root.mkdir()
    card = snapshot_from_record(
        record(tmp_path), read_project=lambda _: project(tmp_path, data_root=data_root)
    )
    assert "integrity not checked" in card.fits_status.lower()
    assert "valid" not in card.fits_status.lower()
    assert not list(data_root.iterdir())


def test_missing_path_keeps_project_identity_and_skips_file_access(tmp_path):
    card = snapshot_from_record(
        record(tmp_path / "moved", available=False),
        read_project=should_not_read,
        read_legacy=should_not_read,
    )
    assert card.project_id == "id-leo-p"
    assert card.fits_status == "Location missing"
    assert "relink" in card.note.lower()


def test_corrupt_metadata_has_honest_nonfatal_state(tmp_path):
    def fail(_):
        raise ValueError("Invalid manifest")
    card = snapshot_from_record(record(tmp_path), read_project=fail)
    assert card.coordinates == "Metadata unavailable"
    assert "Invalid manifest" in card.note
    assert card.fits_status == "Not verified"


def test_legacy_project_is_labelled_separately(tmp_path):
    older = SimpleNamespace(
        field=SimpleNamespace(ra_deg=40.096208, dec_deg=1.229250)
    )
    card = snapshot_from_record(
        record(tmp_path, legacy=True), read_legacy=lambda _: older,
        read_project=should_not_read,
    )
    assert card.legacy
    assert card.dataset == "Older local project"
    assert "import" in card.note.lower()


def test_snapshots_preserve_registry_order(tmp_path):
    rows = [
        record(tmp_path / "a", name="A", available=False),
        {**record(tmp_path / "b", name="B", available=False),
         "project_id": "id-b"},
    ]
    cards = snapshots_from_records(
        rows, read_project=should_not_read, read_legacy=should_not_read
    )
    assert [card.display_name for card in cards] == ["A", "B"]


def test_only_explicit_png_preview_is_considered(tmp_path):
    cache = tmp_path / ".desidera"
    cache.mkdir()
    (cache / "card-preview.png").write_bytes(b"not an actual image")
    (tmp_path / "image-r.fits").write_bytes(b"preserve")
    card = snapshot_from_record(
        record(tmp_path), read_project=lambda _: project(tmp_path)
    )
    assert card.thumbnail == cache / "card-preview.png"
    assert (tmp_path / "image-r.fits").read_bytes() == b"preserve"
