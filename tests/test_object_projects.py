"""Project Hub core contract: portable identity, migration, offline and safety."""
from pathlib import Path
import json

import pytest

from lsbpipe.surveys.legacy import BrickHit, SkyField, product_for
from lsbpipe.surveys.object_projects import (
    SCHEMA, PROJECT_FILENAME, REGISTRY_SCHEMA,
    create_object_project, import_object_project, library_index_path,
    read_object_project, recent_projects, register_object_project,
    relink_object_project, save_object_project, touch_object_project,
)

FIELD = SkyField(40.096208, 1.229250, 5.)


def _old_folder(root):
    data = root / "data"
    for region in ("north", "south"):
        hit = BrickHit("dr9", region, "0401p012", 0, 0, 0, 0, {"r": 1})
        for kind in ("image", "invvar", "nexp", "maskbits"):
            rel = product_for(hit, kind, None if kind == "maskbits" else "r")
            p = data / rel.relative_path
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_bytes((region + kind).encode())
    (data / "download-report.json").write_text(json.dumps({
        "schema": "lsbpipe.legacy-download-report/0.1",
        "field": {"ra_deg": FIELD.ra_deg, "dec_deg": FIELD.dec_deg,
                  "width_arcmin": FIELD.width_arcmin}, "results": [],
    }), encoding="utf-8")
    return data


def test_create_empty_galaxy_workspace_with_unresolved_position(tmp_path):
    p = create_object_project("LSB Candidate α", tmp_path / "LSB_Candidate")
    assert p.status == "NEEDS_TARGET_POSITION"
    assert (p.root / PROJECT_FILENAME).is_file()
    assert not (p.root / "data").exists()
    assert read_object_project(p.root).project_id == p.project_id
    assert read_object_project(p.root).field is None
    recents = recent_projects()
    assert len(recents) == 1
    assert recents[0]["available"]
    assert recents[0]["project_id"] == p.project_id
    assert library_index_path().is_file()


def test_create_known_target_roundtrip_and_identity_is_independent_of_name(tmp_path):
    p = create_object_project("UGC 2162", tmp_path / "UGC_2162", field=FIELD)
    assert p.field == FIELD
    assert p.status == "TARGET_DEFINED_NO_FITS"
    p.display_name = "UGC 2162 — independent photometry"
    p.aliases = ["UGC 2162", "test alias"]
    p.selected_dataset = ("dr9", "south", "0401p012", "r")
    p.catalog_match = {"provider": "SIMBAD", "confirmed": False,
                       "note": "no network lookup performed"}
    save_object_project(p)
    restored = read_object_project(p.root)
    assert restored.project_id == p.project_id
    assert restored.display_name == p.display_name
    assert restored.aliases == p.aliases
    assert restored.field == FIELD
    assert restored.selected_dataset == p.selected_dataset
    assert restored.catalog_match == p.catalog_match
    touch_object_project(restored)
    assert recent_projects()[0]["display_name"] == restored.display_name


def test_import_existing_old_ugc_2162_preserves_all_bytes_and_same_uuid(tmp_path):
    root = tmp_path / "UGC_2162"
    data = _old_folder(root)
    legacy = root / "lsbpipe-project.json"
    legacy.write_text(json.dumps({
        "schema": "lsbpipe.local-project/0.1",
        "field": {"ra_deg": FIELD.ra_deg, "dec_deg": FIELD.dec_deg,
                  "width_arcmin": FIELD.width_arcmin},
        "selected": ["dr9", "south", "0401p012", "r"],
    }), encoding="utf-8")
    originals = {p: p.read_bytes() for p in data.rglob("*.fits.fz")}
    old_metadata = legacy.read_bytes()
    p = import_object_project(root)
    assert p.display_name == "UGC_2162"
    assert p.field == FIELD
    assert p.data_root == data
    assert p.selected_dataset == ("dr9", "south", "0401p012", "r")
    again = import_object_project(root)
    via_data = import_object_project(data)
    assert again.project_id == via_data.project_id == p.project_id
    assert legacy.read_bytes() == old_metadata
    assert all(file.read_bytes() == payload for file, payload in originals.items())
    assert read_object_project(root).project_id == p.project_id


def test_portable_relative_root_and_relink_exact_uuid(tmp_path):
    old = tmp_path / "galaxy"
    project = create_object_project("Galaxy A", old, field=FIELD)
    project.data_root = old / "data"
    save_object_project(project)
    manifest = json.loads((old / PROJECT_FILENAME).read_text("utf-8"))
    assert manifest["data_root"] == {"relative": "data"}
    moved = tmp_path / "new-location" / "galaxy"
    moved.parent.mkdir()
    old.rename(moved)
    recent = recent_projects()
    assert recent[0]["available"] is False
    linked = relink_object_project(project.project_id, moved)
    assert linked.project_id == project.project_id
    assert linked.data_root == moved / "data"
    assert recent_projects()[0]["root"] == str(moved)
    other = create_object_project("Galaxy B", tmp_path / "other")
    with pytest.raises(ValueError, match="another project"):
        relink_object_project(project.project_id, other.root)


def test_explicit_separate_project_same_galaxy_has_different_identity(tmp_path):
    a = create_object_project("Galaxy X", tmp_path / "X1", field=FIELD)
    with pytest.raises(ValueError, match="already in the library"):
        create_object_project("GALAXY X", tmp_path / "X2", field=FIELD)
    b = create_object_project("Galaxy X independent", tmp_path / "X3", field=FIELD)
    assert a.project_id != b.project_id
    assert a.field == b.field
    assert len(recent_projects()) == 2


def test_refuse_nonempty_directory_or_unrecognized_metadata(tmp_path):
    root = tmp_path / "occupied"
    root.mkdir()
    (root / "notes.txt").write_text("keep")
    with pytest.raises(FileExistsError, match="not empty"):
        create_object_project("Test", root)
    assert (root / "notes.txt").read_text() == "keep"
    new = create_object_project("Empty", tmp_path / "empty")
    path = new.root / PROJECT_FILENAME
    path.write_text(json.dumps({"schema": "other-format"}))
    with pytest.raises(ValueError, match="Unsupported project schema"):
        read_object_project(new.root)
    with pytest.raises(ValueError, match="overwrite"):
        save_object_project(new)


def test_corrupted_registry_fails_closed_without_mutating_it(tmp_path):
    file = library_index_path()
    file.parent.mkdir(parents=True)
    file.write_text('{"schema": "foreign-index", "projects":[]}')
    old = file.read_bytes()
    with pytest.raises(ValueError, match="registry schema"):
        recent_projects()
    with pytest.raises(ValueError, match="registry schema"):
        create_object_project("New", tmp_path / "new")
    assert file.read_bytes() == old


def test_reject_path_traversal_in_saved_data_root(tmp_path):
    p = create_object_project("Test", tmp_path / "test")
    file = p.root / PROJECT_FILENAME
    raw = json.loads(file.read_text())
    raw["data_root"] = {"relative": "../outside"}
    file.write_text(json.dumps(raw))
    with pytest.raises(ValueError, match="Unsafe"):
        read_object_project(p.root)


def test_last_project_from_older_lsbpipe_is_visible_without_mutation(tmp_path):
    from lsbpipe.surveys.local_project import remember_project_directory
    root = tmp_path / "UGC_2162"
    data = _old_folder(root)
    remember_project_directory(root)
    rows = recent_projects()
    assert len(rows) == 1
    assert rows[0]["legacy"] is True
    assert rows[0]["display_name"] == "UGC_2162"
    assert rows[0]["available"] is True
    assert not (root / PROJECT_FILENAME).exists()
    project = import_object_project(root)
    updated = recent_projects()
    assert len(updated) == 1
    assert updated[0]["project_id"] == project.project_id
    assert not updated[0].get("legacy", False)
    assert all(p.is_file() for p in data.rglob("*.fits.fz"))


def test_object_passport_radii_survive_restart_and_old_schema_defaults(tmp_path):
    project = create_object_project("Radii", tmp_path / "Radii", field=FIELD)
    project.target_radius_arcsec = 35.
    project.sky_inner_arcsec = 65.
    project.sky_outer_arcsec = 110.
    save_object_project(project)
    reopened = read_object_project(project.root)
    assert (reopened.target_radius_arcsec, reopened.sky_inner_arcsec,
            reopened.sky_outer_arcsec) == (35., 65., 110.)
    meta = reopened.root / PROJECT_FILENAME
    record = json.loads(meta.read_text(encoding="utf-8"))
    for key in ("target_radius_arcsec", "sky_inner_arcsec", "sky_outer_arcsec"):
        del record[key]
    meta.write_text(json.dumps(record), encoding="utf-8")
    older = read_object_project(reopened.root)
    assert (older.target_radius_arcsec, older.sky_inner_arcsec,
            older.sky_outer_arcsec) == (60., 90., 180.)


def test_old_ugc_import_migrates_radii_while_preserving_legacy_bytes(tmp_path):
    root = tmp_path / "UGC_2162"
    _old_folder(root)
    legacy = root / "lsbpipe-project.json"
    old = {
        "schema": "lsbpipe.local-project/0.1",
        "field": {"ra_deg": FIELD.ra_deg, "dec_deg": FIELD.dec_deg,
                  "width_arcmin": FIELD.width_arcmin},
        "selected": ["dr9", "north", "0401p012", "r"],
        "target_radius_arcsec": 35.,
        "sky_inner_arcsec": 65.,
        "sky_outer_arcsec": 110.,
    }
    legacy.write_text(json.dumps(old), encoding="utf-8")
    before = legacy.read_bytes()
    project = import_object_project(root)
    assert (project.target_radius_arcsec, project.sky_inner_arcsec,
            project.sky_outer_arcsec) == (35., 65., 110.)
    assert legacy.read_bytes() == before
    again = read_object_project(root)
    assert (again.target_radius_arcsec, again.sky_inner_arcsec,
            again.sky_outer_arcsec) == (35., 65., 110.)


def test_authoritative_object_target_and_dataset_override_stale_legacy(tmp_path):
    from lsbpipe.surveys.local_project import open_local_project
    root = tmp_path / "UGC_2162"
    data = _old_folder(root)
    legacy = root / "lsbpipe-project.json"
    legacy.write_text(json.dumps({
        "schema": "lsbpipe.local-project/0.1",
        "field": {"ra_deg": 40.096208, "dec_deg": 1.22925,
                  "width_arcmin": 5.0},
        "selected": ["dr9", "north", "0401p012", "r"],
    }), encoding="utf-8")
    object_record = import_object_project(root)
    updated_field = SkyField(40.101234, 1.225678, 6.)
    object_record.field = updated_field
    object_record.selected_dataset = ("dr9", "south", "0401p012", "r")
    object_record.target_radius_arcsec = 35.
    object_record.sky_inner_arcsec = 65.
    object_record.sky_outer_arcsec = 110.
    save_object_project(object_record)
    before = {p: p.read_bytes() for p in root.rglob("*") if p.is_file()}
    inventory = open_local_project(
        data, authoritative_field=object_record.field,
        authoritative_selected=object_record.selected_dataset)
    assert inventory.field == updated_field
    assert inventory.field_source == "object-project"
    assert inventory.selected == object_record.selected_dataset
    assert inventory.hits[0].ra1 == updated_field.ra_deg
    assert legacy.read_bytes() == before[legacy]
    assert all(p.read_bytes() == content for p, content in before.items())
    # Older standalone imports still behave exactly as before.
    older = open_local_project(root)
    assert older.field == FIELD
    assert older.selected == ("dr9", "north", "0401p012", "r")


def test_stale_object_dataset_falls_back_to_local_inventory_not_old_selection(
        tmp_path):
    from lsbpipe.surveys.local_project import open_local_project
    root = tmp_path / "UGC_2162"
    _old_folder(root)
    legacy = root / "lsbpipe-project.json"
    legacy.write_text(json.dumps({
        "schema": "lsbpipe.local-project/0.1",
        "field": {"ra_deg": FIELD.ra_deg, "dec_deg": FIELD.dec_deg,
                  "width_arcmin": FIELD.width_arcmin},
        "selected": ["dr9", "south", "0401p012", "r"],
    }), encoding="utf-8")
    requested = ("dr9", "south", "9999p999", "r")
    inventory = open_local_project(
        root, authoritative_field=FIELD, authoritative_selected=requested)
    assert inventory.selected in inventory.complete
    assert inventory.selected != requested
    assert inventory.selected != ("dr9", "south", "0401p012", "r")



def test_forget_available_project_removes_only_registry_entry(tmp_path):
    from lsbpipe.surveys.object_projects import forget_recent_project
    first = create_object_project("UGC 2162", tmp_path / "real")
    other = create_object_project("Independent dwarf", tmp_path / "other")
    fits = first.root / "original-image-r.fits.fz"
    fits.write_bytes(b"unchanged original data")
    before = (first.root / PROJECT_FILENAME).read_bytes()
    other_before = (other.root / PROJECT_FILENAME).read_bytes()
    forget_recent_project(first.project_id, first.root)
    assert [r["project_id"] for r in recent_projects()] == [other.project_id]
    assert fits.read_bytes() == b"unchanged original data"
    assert (first.root / PROJECT_FILENAME).read_bytes() == before
    assert (other.root / PROJECT_FILENAME).read_bytes() == other_before
    assert library_index_path().is_file()
    # Undo by import: same UUID, no copy/download and no data deletion.
    imported = import_object_project(first.root)
    assert imported.project_id == first.project_id
    assert {r["project_id"] for r in recent_projects()} == {
        first.project_id, other.project_id}


def test_forget_missing_project_does_not_touch_moved_original(tmp_path):
    from lsbpipe.surveys.object_projects import forget_recent_project
    p = create_object_project("UGC_2162", tmp_path / "old_location")
    moved = tmp_path / "moved_elsewhere"
    p.root.rename(moved)
    assert recent_projects()[0]["available"] is False
    forget_recent_project(p.project_id, p.root)
    assert recent_projects() == []
    assert read_object_project(moved).project_id == p.project_id
    with pytest.raises(ValueError, match="no longer in the library"):
        forget_recent_project(p.project_id, p.root)


def test_forget_requires_exact_uuid_and_root_without_mutation(tmp_path):
    from lsbpipe.surveys.object_projects import forget_recent_project
    from uuid import uuid4
    a = create_object_project("One", tmp_path / "one")
    b = create_object_project("Two", tmp_path / "two")
    old_index = library_index_path().read_bytes()
    for pid, root in (
            (str(uuid4()), a.root),
            (a.project_id, b.root),
            (a.project_id, tmp_path / "imaginary"),
            ("not-a-uuid", a.root)):
        with pytest.raises(ValueError):
            forget_recent_project(pid, root)
        assert library_index_path().read_bytes() == old_index
    assert len(recent_projects()) == 2


def test_forget_legacy_suggestion_suppresses_last_session_pointer(tmp_path):
    from lsbpipe.surveys.local_project import (
        remember_project_directory, session_pointer_path)
    from lsbpipe.surveys.object_projects import forget_recent_project
    legacy = tmp_path / "historic_UGC_2162"
    _old_folder(legacy)
    remember_project_directory(legacy)
    original_pointer = session_pointer_path().read_bytes()
    suggestion = recent_projects()[0]
    assert suggestion["legacy"] is True
    assert suggestion["available"]
    forget_recent_project(suggestion["project_id"], suggestion["root"])
    assert recent_projects() == []
    assert session_pointer_path().read_bytes() == original_pointer
    assert not (legacy / PROJECT_FILENAME).exists()
    # A later normal registration must preserve the legacy dismissal.
    new = create_object_project("New candidate", tmp_path / "new")
    assert [p["project_id"] for p in recent_projects()] == [new.project_id]
    assert session_pointer_path().read_bytes() == original_pointer
    # Explicit import is allowed and makes the same old folder visible.
    imported = import_object_project(legacy)
    assert imported.project_id in {
        p["project_id"] for p in recent_projects()}


def test_forget_corrupt_registry_fails_closed(tmp_path):
    from lsbpipe.surveys.object_projects import forget_recent_project
    p = create_object_project("Candidate", tmp_path / "candidate")
    index = library_index_path()
    index.write_text('{"schema":"foreign","projects":[]}')
    original = index.read_bytes()
    with pytest.raises(ValueError, match="registry schema"):
        forget_recent_project(p.project_id, p.root)
    assert index.read_bytes() == original
    assert (p.root / PROJECT_FILENAME).is_file()
