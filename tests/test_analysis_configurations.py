"""P0.8a immutable science setup regression tests (offline only)."""
import json
from uuid import uuid4

import pytest

from lsbpipe.surveys.legacy import SkyField
from lsbpipe.surveys.object_projects import (
    create_object_project, read_object_project, save_object_project,
)
from lsbpipe.surveys.dataset_registry import (
    read_dataset_index, registry_path, CHECKSUM_VERIFIED,
)
from lsbpipe.surveys.analysis_configurations import (
    capture_analysis_setup, freeze_analysis_configuration,
    read_analysis_configuration, list_analysis_configurations,
    compare_current_setup,
)


@pytest.fixture
def project(tmp_path):
    return create_object_project(
        "UGC 2162", tmp_path / "UGC_2162",
        field=SkyField(40.096208, 1.22925, 5))


def test_target_only_freeze_is_immutable_and_not_science(project):
    original = project.root / "original-r.fits.fz"
    original.write_bytes(b"REAL_ORIGINAL_BYTES_NOT_ALTERED")
    manifest = project.root / "lsbpipe-object.json"
    before = manifest.read_bytes()
    digest, item = freeze_analysis_configuration(project.root, "First exploration")
    assert item["status"] == "FROZEN_NOT_RUN"
    assert item["science_status"] == "NOT_SCIENTIFICALLY_ADOPTED"
    assert item["setup"]["input_evidence_state"] == "SETUP_ONLY_NO_DATASET"
    assert item["input_fits_rehashed_during_freeze"] is False
    assert read_analysis_configuration(project.root, digest) == item
    assert len(list_analysis_configurations(project.root)) == 1
    assert compare_current_setup(project.root, digest)["state"] == (
        "SETUP_METADATA_MATCHES_INPUTS_NOT_REHASHED")
    assert manifest.read_bytes() == before
    assert original.read_bytes() == b"REAL_ORIGINAL_BYTES_NOT_ALTERED"


def test_setup_changed_yields_stale_without_touching_old_snapshot(project):
    digest, item = freeze_analysis_configuration(project.root, "Baseline")
    project.target_radius_arcsec = 40
    save_object_project(project)
    assert compare_current_setup(project.root, digest)["state"] == (
        "PROJECT_OR_DATASET_METADATA_CHANGED")
    assert read_analysis_configuration(project.root, digest) == item
    digest2, newer = freeze_analysis_configuration(project.root, "Alternate radius")
    assert newer["setup_sha256"] != item["setup_sha256"]
    assert digest2 != digest
    assert len(list_analysis_configurations(project.root)) == 2


def test_recorded_sha_is_preserved_without_claiming_rehash(project):
    project.selected_dataset = ("dr9", "north", "0401p012", "r")
    save_object_project(project)
    index = read_dataset_index(project.root)
    source_id = str(uuid4())
    asset_id = str(uuid4())
    index["sources"].append({
        "source_id": source_id, "release": "dr9", "region": "north",
        "provider": "legacy_surveys"})
    index["assets"].append({
        "asset_id": asset_id, "source_id": source_id,
        "brick": "0401p012", "band": "r", "product_kind": "image",
        "product_family": "resampled_coadd", "size_bytes": 123,
        "sha256": "a" * 64, "checksum_status": CHECKSUM_VERIFIED,
        "fits_status": "NOT_VALIDATED"})
    file = registry_path(project.root)
    file.parent.mkdir(parents=True, exist_ok=True)
    file.write_text(json.dumps(index))
    digest, record = freeze_analysis_configuration(project.root, "DR9 baseline")
    assets = record["setup"]["selected_registered_assets"]
    assert len(assets) == 1
    assert assets[0]["asset_id"] == asset_id
    assert assets[0]["sha256_recorded"] == "a" * 64
    assert record["setup"]["input_evidence_state"] == (
        "RECORDED_SHA256_NOT_REVALIDATED")
    assert compare_current_setup(project.root, digest)["input_fits_rehashed"] is False
    # Different index status makes the old frozen setup stale.
    index["assets"][0]["checksum_status"] = "CHANGED_SINCE_CHECKSUM"
    file.write_text(json.dumps(index))
    assert compare_current_setup(project.root, digest)["state"] == (
        "PROJECT_OR_DATASET_METADATA_CHANGED")


def test_modified_configuration_fails_closed(project):
    digest, _ = freeze_analysis_configuration(project.root, "Baseline")
    path = project.root / "analyses" / "configurations" / f"config-{digest}.json"
    record = json.loads(path.read_text())
    record["setup"]["parameters"]["target_radius_arcsec"] = 9999
    path.write_text(json.dumps(record))
    with pytest.raises(ValueError, match="corrupted"):
        read_analysis_configuration(project.root, digest)
    with pytest.raises(ValueError):
        list_analysis_configurations(project.root)


def test_cross_project_configuration_does_not_adopt(project, tmp_path):
    other = create_object_project("Other", tmp_path / "Other")
    digest, _ = freeze_analysis_configuration(project.root, "Baseline")
    source = project.root / "analyses" / "configurations" / f"config-{digest}.json"
    target = other.root / "analyses" / "configurations"
    target.mkdir(parents=True)
    (target / source.name).write_bytes(source.read_bytes())
    with pytest.raises(ValueError, match="another project"):
        read_analysis_configuration(other.root, digest)


@pytest.mark.parametrize("bad_label", ["", "  ", "a" * 101, "newline\nname"])
def test_invalid_name_does_not_create_directory(project, bad_label):
    with pytest.raises(ValueError, match="label"):
        freeze_analysis_configuration(project.root, bad_label)
    assert not (project.root / "analyses" / "configurations").exists()
