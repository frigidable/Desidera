"""Non-measuring file-presence audit for the first photometry screen."""
from lsbpipe.photometry.readiness import (
    SCIENCE_STATE, inspect_photometry_inputs)
from lsbpipe.surveys.dataset_integrity import _expected_paths
from lsbpipe.surveys.object_projects import (
    create_object_project, save_object_project)


def test_unidentified_project_has_no_implicit_science_input(tmp_path, monkeypatch):
    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path / "session"))
    project = create_object_project("Draft", tmp_path / "draft")
    manifest = (project.root / "lsbpipe-object.json")
    original = manifest.read_bytes()
    info = inspect_photometry_inputs(project.root)
    assert info["project_id"] == project.project_id
    assert info["status"] == "NO_SELECTED_DATASET"
    assert info["selected_dataset"] is None
    assert info["products"] == ()
    assert info["science_state"] == SCIENCE_STATE
    assert manifest.read_bytes() == original
    assert not (project.root / "photometry-runs").exists()


def test_quartet_presence_never_certifies_fits_or_photometry(
        tmp_path, monkeypatch):
    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path / "session"))
    project = create_object_project("UGC 2162", tmp_path / "ugc")
    project.data_root = tmp_path / "original_data"
    project.data_root.mkdir()
    project.selected_dataset = ("dr9", "south", "0401p012", "r")
    save_object_project(project)
    paths = _expected_paths(project.data_root, project.selected_dataset)
    present = [paths["image"], paths["nexp"], paths["maskbits"]]
    for path in present:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"not a FITS file: existence check only")
    first = inspect_photometry_inputs(project.root)
    assert first["status"] == "INCOMPLETE_QUARTET"
    assert first["present_count"] == 3
    assert [p["kind"] for p in first["products"] if not p["present"]] == ["invvar"]
    paths["invvar"].write_bytes(b"still not a FITS")
    second = inspect_photometry_inputs(project.root)
    assert second["status"] == "FILES_PRESENT_UNVERIFIED"
    assert second["present_count"] == 4
    assert second["science_state"] == SCIENCE_STATE
    assert all(path.read_bytes().startswith(b"not a FITS") or
               path.read_bytes() == b"still not a FITS" for path in paths.values())
    assert not (project.root / "qa-runs").exists()


def test_missing_data_root_is_distinct_from_missing_products(
        tmp_path, monkeypatch):
    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path / "session"))
    project = create_object_project("Draft", tmp_path / "draft")
    project.selected_dataset = ("dr9", "north", "0401p012", "g")
    save_object_project(project)
    info = inspect_photometry_inputs(project.root)
    assert info["status"] == "DATA_ROOT_UNAVAILABLE"
    assert info["present_count"] == 0
    assert info["products"] == ()


def test_unsupported_optical_band_does_not_construct_fake_paths(
        tmp_path, monkeypatch):
    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path / "session"))
    project = create_object_project("Draft", tmp_path / "draft")
    project.selected_dataset = ("dr9", "north", "0401p012", "i")
    save_object_project(project)
    info = inspect_photometry_inputs(project.root)
    assert info["status"] == "UNSUPPORTED_SELECTION"
    assert not info["products"]
