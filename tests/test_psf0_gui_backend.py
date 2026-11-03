"""PSF-0 GUI contract: reuse saved Leo P galaxy project, no CLI directories."""
from __future__ import annotations

import pytest

from lsbpipe.photometry.psf_gui_backend import inspect_psf0_bands
from lsbpipe.surveys.object_projects import create_object_project, save_object_project
from lsbpipe.surveys.legacy import SkyField
from lsbpipe.surveys.dataset_integrity import _expected_paths


def _leo_project(tmp_path, monkeypatch, *, galaxy="Leo P"):
    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path / "session"))
    project = create_object_project(
        galaxy, tmp_path / "Leo_P",
        field=SkyField(155.438012, 18.088025, 5.))
    project.data_root = tmp_path / "Leo_P" / "originals"
    project.data_root.mkdir()
    project.selected_dataset = ("dr11", "south", "1554p180", "g")
    save_object_project(project)
    return project


def test_psf_gui_lists_r_and_g_independently_from_selected_g(tmp_path, monkeypatch):
    project = _leo_project(tmp_path, monkeypatch)
    rpaths = _expected_paths(project.data_root, ("dr11", "south", "1554p180", "r"))
    for path in rpaths.values():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"NOT A FITS - PATH PRESENCE ONLY")
    original = (project.root / "lsbpipe-object.json").read_bytes()
    info = inspect_psf0_bands(project.root)
    assert info["project_id"] == project.project_id
    assert info["dataset"] == ("dr11", "south", "1554p180", "g")
    bands = {row["band"]: row for row in info["bands"]}
    assert bands["r"]["paths_present_unverified"] is True
    assert bands["r"]["selected"] is False
    assert bands["g"]["selected"] is True
    assert bands["g"]["present_count"] == 1  # shared maskbits only
    assert bands["r"]["present_count"] == 4
    assert bands["r"]["files"] == {key: True for key in rpaths}
    assert info["status"] == "PATHS_ONLY_NOT_SCIENCE_VALIDATED"
    assert (project.root / "lsbpipe-object.json").read_bytes() == original
    assert not (project.root / "analyses").exists()


def test_psf_gui_refuses_other_project_without_touching_fits(tmp_path, monkeypatch):
    p = _leo_project(tmp_path, monkeypatch, galaxy="UGC 2162")
    with pytest.raises(ValueError, match="Leo P"):
        inspect_psf0_bands(p.root)
