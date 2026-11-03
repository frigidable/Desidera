"""P10g: selected-object survey readiness with bounded read-only file checks."""
from pathlib import Path
import pytest

from lsbpipe.surveys.galaxy_passport_readiness import inspect_passport_readiness
from lsbpipe.surveys.legacy import BrickHit, SkyField, product_for
from lsbpipe.surveys.object_projects import ObjectProject


def project(tmp_path, *, selection=("dr11", "south", "1554p180", "r"),
            linked=True, field=True):
    return ObjectProject(
        project_id="35a8b364-3819-427c-a835-83c1a542328f",
        display_name="Leo P", root=tmp_path,
        field=SkyField(155.438012, 18.088025, 4.) if field else None,
        selected_dataset=selection,
        data_root=tmp_path if linked else None)


def product_path(root, kind, band="r", *, data_child=False):
    hit = BrickHit("dr11", "south", "1554p180", 0, 0, 0, 0, {"r": None})
    p = product_for(hit, kind, band if kind != "maskbits" else None)
    return root / ("data" if data_child else "") / p.relative_path


@pytest.mark.parametrize("field,selection,linked,code", [
    (False, None, False, "TARGET_NOT_SET"),
    (True, None, False, "DATASET_NOT_SELECTED"),
    (True, ("dr11", "south", "1554p180", "r"), False,
     "DATA_ROOT_NOT_LINKED"),
])
def test_readiness_missing_stages_do_not_invent_original_fits(
        tmp_path, field, selection, linked, code):
    state = inspect_passport_readiness(
        project(tmp_path, selection=selection, linked=linked, field=field))
    assert state.code == code
    assert state.file_count == 0
    assert state.next_action_label.startswith("Open project")


def test_selected_but_empty_data_root_requires_original_files(tmp_path):
    status = inspect_passport_readiness(project(tmp_path))
    assert status.code == "NO_MATCHING_ORIGINAL_FITS"
    assert status.missing_kinds == ("image", "invvar", "nexp", "maskbits")
    assert "no original fits" in status.observation_label.lower()
    assert status.expected_count == 4


def test_incomplete_vs_full_quartet_are_presence_only(tmp_path):
    p = project(tmp_path)
    for kind in ("image", "invvar", "nexp"):
        path = product_path(p.root, kind)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"not-a-valid-FITS")
    partial = inspect_passport_readiness(p)
    assert partial.code == "INCOMPLETE_ORIGINAL_FITS"
    assert partial.file_count == 3
    assert partial.missing_kinds == ("maskbits",)
    assert "integrity" in partial.explanation.lower()

    mask = product_path(p.root, "maskbits")
    mask.write_bytes(b"not-a-valid-FITS")
    full = inspect_passport_readiness(p)
    assert full.code == "FOUR_ORIGINAL_FITS_PRESENT_UNVERIFIED"
    assert full.file_count == 4
    assert "have not been established" in full.explanation
    assert "inspect FITS" in full.next_action_label
    assert all(path.read_bytes() == b"not-a-valid-FITS"
               for path in [mask, product_path(p.root, "image")])


def test_nested_data_directory_source_and_zero_byte_does_not_count(tmp_path):
    p = project(tmp_path)
    image = product_path(tmp_path, "image", data_child=True)
    image.parent.mkdir(parents=True)
    image.write_bytes(b"original science image")
    invvar = product_path(tmp_path, "invvar", data_child=True)
    invvar.write_bytes(b"")
    status = inspect_passport_readiness(p)
    assert status.file_count == 1
    assert status.code == "INCOMPLETE_ORIGINAL_FITS"
    assert "invvar" in status.missing_kinds


def test_missing_root_and_symlinked_file_do_not_count(tmp_path):
    missing = tmp_path / "absent"
    p = project(tmp_path)
    p.data_root = missing
    assert inspect_passport_readiness(p).code == "DATA_ROOT_MISSING"
    existing = product_path(tmp_path, "image")
    existing.parent.mkdir(parents=True, exist_ok=True)
    other = tmp_path / "elsewhere"
    other.write_bytes(b"FITS")
    try:
        existing.symlink_to(other)
    except (OSError, NotImplementedError):
        pytest.skip("Symlinks unavailable under this Windows runner")
    p.data_root = tmp_path
    assert inspect_passport_readiness(p).file_count == 0


def test_inspection_never_calls_fits_open_or_hashes_files(tmp_path, monkeypatch):
    import astropy.io.fits
    p = project(tmp_path)
    path = product_path(tmp_path, "image")
    path.parent.mkdir(parents=True)
    path.write_bytes(b"not a FITS; only presence should be seen")
    monkeypatch.setattr(astropy.io.fits, "open", lambda *a, **k:
                        pytest.fail("Readiness must not decode science FITS"))
    assert inspect_passport_readiness(p).file_count == 1
