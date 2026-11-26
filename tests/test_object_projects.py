"""Project Hub core contract: portable identity, migration, offline and safety."""
from pathlib import Path
import json
import os
import subprocess
import sys
import time

import pytest

from lsbpipe.surveys.legacy import BrickHit, SkyField, product_for
from lsbpipe.surveys.object_projects import (
    SCHEMA, PROJECT_FILENAME, REGISTRY_SCHEMA,
    create_object_project, import_object_project, library_index_path,
    read_object_project, read_registered_object_project, recent_projects,
    register_object_project,
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



def test_project_data_root_equal_project_root_roundtrips_and_can_journal(tmp_path):
    """Regression: a valid relative '.' must not block every worker task."""
    from lsbpipe.surveys.operation_journal import (
        start_operation, finish_operation)
    project = create_object_project(
        "UGC 2162", tmp_path / "UGC_2162", field=FIELD)
    original = project.root / "unmodified-image.fits.fz"
    original.write_bytes(b"original bytes are not changed")
    project.data_root = project.root
    save_object_project(project)
    manifest = project.root / PROJECT_FILENAME
    raw = json.loads(manifest.read_text(encoding="utf-8"))
    assert raw["data_root"] == {"relative": "."}
    again = read_object_project(project.root)
    assert again.data_root == project.root
    started = start_operation(project.root, project.project_id, "survey_search")
    finish_operation(started, "COMPLETED")
    assert original.read_bytes() == b"original bytes are not changed"


@pytest.mark.parametrize("unsafe", ["", "../elsewhere", "data/../../elsewhere"])
def test_relative_data_root_still_rejects_empty_and_traversal(tmp_path, unsafe):
    project = create_object_project("Target", tmp_path / "target")
    manifest = project.root / PROJECT_FILENAME
    state = json.loads(manifest.read_text(encoding="utf-8"))
    state["data_root"] = {"relative": unsafe}
    manifest.write_text(json.dumps(state), encoding="utf-8")
    with pytest.raises(ValueError, match="Unsafe relative data directory"):
        read_object_project(project.root)



def test_object_project_refuses_symlinked_manifest_on_read_and_save(
        tmp_path):
    """A symlinked project passport must not import foreign JSON evidence."""
    project = create_object_project("Protected", tmp_path / "Protected")
    manifest = project.root / PROJECT_FILENAME
    original = manifest.read_bytes()
    outside = tmp_path / "external-object-manifest.json"
    outside.write_bytes(original)
    manifest.unlink()
    try:
        manifest.symlink_to(outside)
    except (OSError, NotImplementedError):
        pytest.skip("Filesystem does not support file symlinks")
    assert manifest.is_symlink()
    with pytest.raises(ValueError, match="Unsafe symlinked object project"):
        read_object_project(project.root)
    with pytest.raises(ValueError, match="Unsafe symlinked object project"):
        save_object_project(project)
    assert outside.read_bytes() == original
    assert manifest.is_symlink()


def test_portable_relative_data_root_refuses_symlink_escape(
        tmp_path):
    """An unchanged lexical path can escape through a junction/symlink."""
    project = create_object_project("Target", tmp_path / "Target")
    outside = tmp_path / "external-original-fits"
    outside.mkdir()
    (outside / "ORIGINAL.fits.fz").write_bytes(b"original unchanged")
    alias = project.root / "data-alias"
    try:
        alias.symlink_to(outside, target_is_directory=True)
    except (OSError, NotImplementedError):
        pytest.skip("Filesystem does not support directory symlinks")
    passport = project.root / PROJECT_FILENAME
    record = json.loads(passport.read_text("utf-8"))
    record["data_root"] = {"relative": "data-alias"}
    passport.write_text(json.dumps(record), encoding="utf-8")
    with pytest.raises(ValueError, match="Unsafe relative data directory"):
        read_object_project(project.root)
    assert (outside / "ORIGINAL.fits.fz").read_bytes() == b"original unchanged"

    # Explicitly external ABSOLUTE paths remain valid for normal
    # researcher-managed original FITS directories.
    record["data_root"] = {"absolute": str(outside.resolve())}
    passport.write_text(json.dumps(record), encoding="utf-8")
    assert read_object_project(project.root).data_root == outside.resolve()


def test_library_registry_refuses_symlinked_external_index(
        tmp_path):
    project = create_object_project("One", tmp_path / "One")
    idx = library_index_path()
    original = idx.read_bytes()
    outside = tmp_path / "untrusted-projects.json"
    outside.write_bytes(original)
    idx.unlink()
    try:
        idx.symlink_to(outside)
    except (OSError, NotImplementedError):
        pytest.skip("Filesystem does not support file symlinks")
    with pytest.raises(ValueError, match="Unsafe symlinked project library"):
        recent_projects()
    with pytest.raises(ValueError, match="Unsafe symlinked project library"):
        register_object_project(project)
    assert outside.read_bytes() == original


def test_registry_registration_uses_single_consistent_index_snapshot(
        tmp_path, monkeypatch):
    """A second read could mix old projects with newer dismissed entries."""
    import lsbpipe.surveys.object_projects as module

    project = create_object_project("Single snapshot", tmp_path / "Single")
    original = module._registry_payload
    reads = []
    def observed_read():
        reads.append(True)
        return original()
    monkeypatch.setattr(module, "_registry_payload", observed_read)
    module.register_object_project(project)
    assert len(reads) == 1
    assert recent_projects()[0]["project_id"] == project.project_id



def test_stale_object_passport_cannot_overwrite_updated_target_or_qa(
        tmp_path):
    """Separate live GUI snapshots must not silently lose science pointers."""
    project = create_object_project(
        "Dwarf", tmp_path / "Dwarf", field=FIELD)
    observed_a = read_object_project(project.root)
    observed_b = read_object_project(project.root)
    assert observed_a.manifest_sha256 == observed_b.manifest_sha256
    assert observed_a.manifest_observed and observed_b.manifest_observed

    source_original = project.root / "coadd-original-image.fits.fz"
    source_original.write_bytes(b"keep original pixel stream")
    observed_a.selected_dataset = ("dr11", "south", "0401p012", "r")
    observed_a.last_science = {
        "action": "photometry", "path": "qa-runs/first.json",
        "scientific_status": "DIAGNOSTIC_NOT_ADOPTED"}
    observed_a.target_radius_arcsec = 45.
    observed_a.sky_inner_arcsec = 95.
    observed_a.sky_outer_arcsec = 180.
    manifest = save_object_project(observed_a)
    original_revision = manifest.read_bytes()
    observed_b.catalog_match = {"provider": "SIMBAD", "confirmed": False}
    with pytest.raises(RuntimeError, match="changed in another process"):
        save_object_project(observed_b)
    assert manifest.read_bytes() == original_revision
    assert source_original.read_bytes() == b"keep original pixel stream"

    observed_a.display_name = "Dwarf updated"
    save_object_project(observed_a)
    latest = read_object_project(project.root)
    assert latest.last_science == observed_a.last_science
    assert latest.selected_dataset == observed_a.selected_dataset
    assert latest.display_name == "Dwarf updated"
    assert latest.manifest_sha256 == observed_a.manifest_sha256
    assert latest.catalog_match is None
    assert json.loads(manifest.read_text())["science_status"] == (
        "DIAGNOSTIC_NOT_ADOPTED")


def test_two_process_object_passport_stale_save_conflict(
        tmp_path, monkeypatch):
    """Actual OS processes open the same passport before either can write."""
    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path / "session"))
    project = create_object_project("Galaxy", tmp_path / "Galaxy", field=FIELD)
    root = project.root
    first_ready = tmp_path / "first-ready"
    second_ready = tmp_path / "second-ready"
    first_go = tmp_path / "first-go"
    second_go = tmp_path / "second-go"
    code = (
        "import pathlib,sys,time\n"
        "from lsbpipe.surveys.object_projects import read_object_project,save_object_project\n"
        "root,ready,go=map(pathlib.Path,sys.argv[1:4])\n"
        "label=sys.argv[4]\n"
        "p=read_object_project(root)\n"
        "ready.write_text('OPENED')\n"
        "while not go.is_file(): time.sleep(.02)\n"
        "p.display_name=label\n"
        "try:\n"
        "    save_object_project(p)\n"
        "except RuntimeError as exc:\n"
        "    assert 'changed in another process' in str(exc),str(exc)\n"
        "    print('CONFLICT',flush=True)\n"
        "else:\n"
        "    print('SAVED',flush=True)\n"
    )
    def start(ready, go, label):
        return subprocess.Popen(
            [sys.executable, "-u", "-c", code,
             str(root), str(ready), str(go), label],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, env=dict(os.environ))
    a = start(first_ready, first_go, "First writer")
    b = start(second_ready, second_go, "Second writer")
    try:
        deadline = time.monotonic() + 18.
        while not (first_ready.is_file() and second_ready.is_file()):
            if a.poll() is not None or b.poll() is not None:
                pytest.fail("Object-project worker exited before opening")
            if time.monotonic() >= deadline:
                pytest.fail("Timed out opening two object-passport workers")
            time.sleep(.02)
        first_go.write_text("GO", encoding="utf-8")
        result_a, error_a = a.communicate(timeout=15)
        assert a.returncode == 0 and result_a.strip() == "SAVED", (
            result_a, error_a)
        second_go.write_text("GO", encoding="utf-8")
        result_b, error_b = b.communicate(timeout=15)
        assert b.returncode == 0 and result_b.strip() == "CONFLICT", (
            result_b, error_b)
        assert read_object_project(root).display_name == "First writer"
        assert (root / "lsbpipe-object.json.lock").is_file()
        assert read_object_project(root).project_id == project.project_id
    finally:
        for worker in (a, b):
            if worker.poll() is None:
                worker.kill()
                worker.communicate(timeout=10)


def test_never_observed_object_passport_cannot_replace_existing_identity(
        tmp_path):
    """Directly constructed project objects cannot guess a live revision."""
    from lsbpipe.surveys.object_projects import ObjectProject
    existing = create_object_project("Original", tmp_path / "original")
    record = existing.root / PROJECT_FILENAME
    prior = record.read_bytes()
    forged = ObjectProject(
        project_id=existing.project_id, display_name="Spoofed",
        root=existing.root, field=FIELD)
    assert not forged.manifest_observed
    with pytest.raises(RuntimeError, match="never opened"):
        save_object_project(forged)
    assert record.read_bytes() == prior


def test_legacy_import_concurrent_processes_reuse_one_project_uuid(
        tmp_path, monkeypatch):
    """Both importers scanned before either could create the object passport.

    A stable library transaction must recheck existence after the
    potentially slow FITS inventory. Independent imports cannot silently
    choose two different immutable project identifiers.
    """
    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path / "session"))
    root = tmp_path / "legacy"
    data = _old_folder(root)
    originals = {p: p.read_bytes() for p in data.rglob("*.fits.fz")}
    first_ready, second_ready = tmp_path / "first-ready", tmp_path / "second-ready"
    release = tmp_path / "allow-import"
    code = (
        "import pathlib,sys,time\n"
        "import lsbpipe.surveys.object_projects as m\n"
        "root,marker,release=map(pathlib.Path,sys.argv[1:])\n"
        "original=m.open_local_project\n"
        "def gated(*args,**kwargs):\n"
        "    inventory=original(*args,**kwargs)\n"
        "    marker.write_text('SCANNED',encoding='utf-8')\n"
        "    while not release.is_file(): time.sleep(.02)\n"
        "    return inventory\n"
        "m.open_local_project=gated\n"
        "project=m.import_object_project(root)\n"
        "print(project.project_id,flush=True)\n"
    )
    def launch(ready):
        return subprocess.Popen(
            [sys.executable, "-u", "-c", code,
             str(root), str(ready), str(release)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, env=dict(os.environ))
    a, b = launch(first_ready), launch(second_ready)
    try:
        deadline = time.monotonic() + 20.
        while not (first_ready.is_file() and second_ready.is_file()):
            if a.poll() is not None or b.poll() is not None:
                # If one importer exits prematurely, the OTHER importer
                # may still be deliberately blocked on the release marker.
                # Release it before communicate(): otherwise Windows CI
                # reports a misleading 2-second subprocess timeout and
                # hides the actual failure output from the early process.
                # The premature exit STILL fails this concurrency test.
                release.write_text("CONTINUE", encoding="utf-8")
                out_a, err_a = a.communicate(timeout=18)
                out_b, err_b = b.communicate(timeout=18)
                pytest.fail(
                    f"Importer stopped before rendezvous: "
                    f"{(a.returncode, out_a, err_a, b.returncode, out_b, err_b)}")
            if time.monotonic() > deadline:
                pytest.fail("Timed out waiting for two scanned FITS inventories")
            time.sleep(.02)
        release.write_text("CONTINUE", encoding="utf-8")
        out_a, err_a = a.communicate(timeout=18)
        out_b, err_b = b.communicate(timeout=18)
        assert a.returncode == b.returncode == 0, (out_a, err_a, out_b, err_b)
        assert out_a.strip() == out_b.strip()
        selected = read_object_project(root)
        assert selected.project_id == out_a.strip()
        rows = recent_projects()
        assert len(rows) == 1
        assert rows[0]["project_id"] == selected.project_id
        assert all(p.read_bytes() == original for p, original in originals.items())
    finally:
        for proc in (a, b):
            if proc.poll() is None:
                proc.kill()
                proc.communicate(timeout=10)


def test_interrupted_import_keeps_passport_uuid_and_recovers_by_reimport(
        tmp_path, monkeypatch):
    """Simulate process exit precisely between passport and library writes."""
    import lsbpipe.surveys.object_projects as module

    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path / "session"))
    root = tmp_path / "legacy"
    data = _old_folder(root)
    originals = {p: p.read_bytes() for p in data.rglob("*.fits.fz")}
    original_register = module._register_object_project_unlocked
    def abort_between_files(project):
        assert (root / PROJECT_FILENAME).is_file()
        raise SystemExit("simulate crash after passport publication")
    monkeypatch.setattr(module, "_register_object_project_unlocked",
                        abort_between_files)
    with pytest.raises(SystemExit, match="simulate crash"):
        import_object_project(root)
    orphan = read_object_project(root)
    assert orphan.project_id
    assert recent_projects() == []

    # The user explicitly imports the same known folder after restarting.
    # Recovery is an index repair, NOT a new project or FITS rewriting.
    monkeypatch.setattr(module, "_register_object_project_unlocked",
                        original_register)
    recovered = import_object_project(root)
    assert recovered.project_id == orphan.project_id
    assert recent_projects()[0]["project_id"] == orphan.project_id
    assert len(recent_projects()) == 1
    assert all(p.read_bytes() == original for p, original in originals.items())


def test_create_rejects_bad_name_before_making_empty_project_directory(
        tmp_path, monkeypatch):
    """Invalid metadata must not leave an unregistered new folder."""
    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path / "session"))
    root = tmp_path / "invalid-galaxy"
    with pytest.raises(ValueError, match="unsupported characters"):
        create_object_project("Bad|Name", root)
    assert not root.exists()
    assert not library_index_path().is_file()


def test_registry_path_replacement_is_not_accepted_as_same_galaxy(tmp_path):
    """The SAME pathname can now hold ANOTHER UUID: fail closed."""
    old = create_object_project("Original galaxy", tmp_path / "old", field=FIELD)
    foreign = create_object_project("Different galaxy", tmp_path / "foreign")
    original_data = old.root / "original-image-g.fits.fz"
    original_data.write_bytes(b"immutable original science image")
    original_registry = library_index_path().read_bytes()
    passport = old.root / PROJECT_FILENAME
    foreign_bytes = (foreign.root / PROJECT_FILENAME).read_bytes()
    passport.write_bytes(foreign_bytes)

    # The explicitly selected old library row cannot open another galaxy,
    # even though its path still exists and points to valid project JSON.
    with pytest.raises(ValueError, match="UUID mismatch"):
        read_registered_object_project(old.root, old.project_id)
    assert read_object_project(old.root).project_id == foreign.project_id

    by_id = {r["project_id"]: r for r in recent_projects()}
    assert not by_id[old.project_id]["available"]
    assert by_id[old.project_id]["availability_issue"] == "UUID_MISMATCH"
    assert by_id[old.project_id]["display_name"] == "Original galaxy"
    assert by_id[foreign.project_id]["available"]
    assert original_data.read_bytes() == b"immutable original science image"
    assert passport.read_bytes() == foreign_bytes
    assert library_index_path().read_bytes() == original_registry


def test_project_library_display_uses_current_passport_without_mutating_registry(
        tmp_path):
    """Project list is navigation; passport retains authoritative name."""
    project = create_object_project("Initial name", tmp_path / "initial")
    index_before = library_index_path().read_bytes()
    project.display_name = "Updated passport name"
    save_object_project(project)
    rows = recent_projects()
    assert len(rows) == 1
    assert rows[0]["available"] is True
    assert rows[0]["display_name"] == "Updated passport name"
    assert library_index_path().read_bytes() == index_before


def test_invalid_passport_is_not_offered_as_valid_registered_galaxy(
        tmp_path):
    project = create_object_project("Galaxy", tmp_path / "galaxy")
    index_before = library_index_path().read_bytes()
    passport = project.root / PROJECT_FILENAME
    passport.write_bytes(b'{"schema": "other.invalid/1"}')
    rows = recent_projects()
    assert len(rows) == 1
    assert rows[0]["available"] is False
    assert rows[0]["availability_issue"] == "INVALID_PASSPORT"
    assert library_index_path().read_bytes() == index_before
