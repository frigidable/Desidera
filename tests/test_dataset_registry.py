"""P0 Dataset Manager foundation: offline, identity-safe, non-science tests."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from lsbpipe.surveys.dataset_integrity import discover_dataset_inventory
from lsbpipe.surveys.dataset_registry import (
    SCHEMA, CHECKSUM_VERIFIED, CHECKSUM_CHANGED,
    read_dataset_index, registry_path, resolve_asset_path,
    import_legacy_fits, register_external_fits,
    verify_asset_checksum, science_cache_key, compatibility_policy,
)
from lsbpipe.surveys.object_projects import (
    create_object_project, read_object_project,
)
from lsbpipe.surveys.legacy import BrickHit, SkyField, product_for


@pytest.fixture
def project(tmp_path, monkeypatch):
    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path / "local"))
    return create_object_project(
        "UGC 2162", tmp_path / "UGC_2162",
        field=SkyField(40.096208, 1.22925, 5.))


def add_legacy_assets(root, *, release="dr9", region="north",
                      brick="0401p012", band="r"):
    hit = BrickHit(release, region, brick, 0, 0, 0, 0, {band: 1})
    files = {}
    for kind in ("image", "invvar", "nexp", "maskbits"):
        rel = product_for(hit, kind,
                          None if kind == "maskbits" else band).relative_path
        path = root / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(f"ORIGINAL_TEST_DATA:{release}/{brick}/{kind}/{band}".encode())
        files[kind] = path
    return files


def test_legacy_import_preserves_selection_uuid_and_files(project):
    project.selected_dataset = ("dr9", "north", "0401p012", "r")
    from lsbpipe.surveys.object_projects import save_object_project
    save_object_project(project)
    originals = add_legacy_assets(project.root / "data")
    before = {name: file.read_bytes() for name, file in originals.items()}
    result = import_legacy_fits(project.root, project.root / "data")
    assert result["sources_total"] == 1
    assert result["assets_total"] == 4
    assert result["science_status"] == "NOT_SCIENTIFICALLY_ADOPTED"
    idx = read_dataset_index(project.root)
    assert idx["schema"] == SCHEMA
    assert idx["project_id"] == project.project_id
    assert {a["product_kind"] for a in idx["assets"]} == {
        "image", "invvar", "nexp", "maskbits"}
    assert all(a["sha256"] is None and a["fits_status"] == "NOT_VALIDATED"
               for a in idx["assets"])
    assert read_object_project(project.root).project_id == project.project_id
    assert read_object_project(project.root).selected_dataset == project.selected_dataset
    assert {name: file.read_bytes() for name, file in originals.items()} == before
    again = import_legacy_fits(project.root, project.root / "data")
    assert again["assets_added"] == 0
    assert again["sources_added"] == 0
    assert [a["asset_id"] for a in idx["assets"]] == [
        a["asset_id"] for a in read_dataset_index(project.root)["assets"]]


def test_two_bricks_and_two_releases_are_separate_assets(project):
    root = project.root / "data"
    add_legacy_assets(root, brick="0401p012")
    add_legacy_assets(root, brick="0400p012")
    add_legacy_assets(root, release="dr11", region="south",
                      brick="0401p012", band="z")
    receipt = import_legacy_fits(project.root, root)
    assert receipt["sources_total"] == 2
    assert receipt["assets_total"] == 12
    idx = read_dataset_index(project.root)
    assert {s["release"] for s in idx["sources"]} == {"dr9", "dr11"}
    assert {a["source_id"] for a in idx["assets"]} == {
        s["source_id"] for s in idx["sources"]}
    assert all(s["observation_lineage"] == "UNKNOWN" for s in idx["sources"])
    assert compatibility_policy(idx["sources"], idx["assets"], mode="joint")["allowed"] is False
    assert compatibility_policy(idx["sources"], idx["assets"], mode="compare")["allowed"] is True


def test_maskbits_shared_between_bands_is_one_file_asset(project):
    root = project.root / "data"
    add_legacy_assets(root, band="g")
    add_legacy_assets(root, band="r")
    import_legacy_fits(project.root, root)
    idx = read_dataset_index(project.root)
    assert len(idx["assets"]) == 7
    assert len([a for a in idx["assets"]
                if a["product_kind"] == "maskbits"]) == 1
    assert len(idx["sources"]) == 1


def test_inventory_and_registry_support_more_than_64_datasets(project):
    root = project.root / "data"
    for n in range(65):
        add_legacy_assets(root, brick=f"{400+n:04d}p012")
    assert len(discover_dataset_inventory(root)) == 65
    summary = import_legacy_fits(project.root, root)
    assert summary["assets_total"] == 65 * 4
    # Import is metadata-only: contents are placeholders, not valid FITS.
    assert all(a["fits_status"] == "NOT_VALIDATED" for a in
               read_dataset_index(project.root)["assets"])


def test_checksum_hashes_are_explicit_and_changes_fail_closed(project):
    files = add_legacy_assets(project.root / "data")
    import_legacy_fits(project.root, project.root / "data")
    idx = read_dataset_index(project.root)
    asset = next(a for a in idx["assets"]
                 if a["product_kind"] == "image")
    receipt = verify_asset_checksum(project.root, asset["asset_id"])
    assert receipt["checksum_status"] == CHECKSUM_VERIFIED
    assert receipt["fits_status"] == "NOT_VALIDATED"
    assert receipt["sha256"] == hashlib.sha256(
        files["image"].read_bytes()).hexdigest()
    prior = science_cache_key(
        stage="sky", algorithm_version="1.0", parameters={"mode": "flat"},
        selected_assets=[next(a for a in read_dataset_index(project.root)["assets"]
                              if a["asset_id"] == asset["asset_id"])])
    files["image"].write_bytes(b"CHANGED_AFTER_CHECKSUM")
    changed = import_legacy_fits(project.root, project.root / "data")
    assert changed["assets_changed"] >= 1
    now = next(a for a in read_dataset_index(project.root)["assets"]
               if a["asset_id"] == asset["asset_id"])
    assert now["checksum_status"] == CHECKSUM_CHANGED
    with pytest.raises(ValueError, match="Unverified FITS fingerprint"):
        science_cache_key(
            stage="sky", algorithm_version="1.0",
            parameters={"mode": "flat"}, selected_assets=[now])
    with pytest.raises(RuntimeError, match="content changed"):
        verify_asset_checksum(project.root, asset["asset_id"])
    assert len(prior) == 64


def test_cache_dependencies_do_not_confuse_background_variants(project):
    add_legacy_assets(project.root / "data")
    import_legacy_fits(project.root, project.root / "data")
    idx = read_dataset_index(project.root)
    with pytest.raises(ValueError, match="Unverified FITS"):
        science_cache_key(stage="photometry", algorithm_version="0.1",
                          parameters={}, selected_assets=idx["assets"][:1])
    a = idx["assets"][0]
    a["sha256"] = hashlib.sha256(b"fixed").hexdigest()
    a["checksum_status"] = CHECKSUM_VERIFIED
    kwargs = {"stage": "sky", "algorithm_version": "v1",
              "selected_assets": [a]}
    key1 = science_cache_key(**kwargs, parameters={"order": 2})
    key2 = science_cache_key(**kwargs, parameters={"order": 3})
    key3 = science_cache_key(
        **kwargs, parameters={"order": 2}, dependency_keys=["mask-v2"])
    assert len({key1, key2, key3}) == 3
    assert science_cache_key(**kwargs, parameters={"order": 2}) == key1


def test_external_fits_not_treated_as_calibrated_survey(project, tmp_path):
    target = tmp_path / "hsc_frame.fits"
    target.write_bytes(b"FITS_PLACEHOLDER")
    result = register_external_fits(
        project.root, target, provider="HSC", band="i")
    idx = read_dataset_index(project.root)
    assert len(idx["sources"]) == len(idx["assets"]) == 1
    src, asset = idx["sources"][0], idx["assets"][0]
    assert result["source_id"] == src["source_id"]
    assert src["calibration_status"] == "UNASSESSED"
    assert asset["fits_status"] == "NOT_VALIDATED"
    assert resolve_asset_path(project.root, asset) == target
    assert compatibility_policy([src], [asset], mode="joint")["allowed"] is False


def test_registry_tampering_and_relative_path_traversal_fail_closed(project):
    add_legacy_assets(project.root / "data")
    import_legacy_fits(project.root, project.root / "data")
    path = registry_path(project.root)
    raw = json.loads(path.read_text())
    raw["project_id"] = "different-project"
    path.write_text(json.dumps(raw))
    with pytest.raises(ValueError, match="different galaxy"):
        read_dataset_index(project.root)
    with pytest.raises(ValueError, match="Unsafe"):
        resolve_asset_path(project.root,
                           {"path_ref": {"relative": "../other.fits"}})


def test_mosaic_geometry_plan_does_not_mean_coaddition_allowed(project):
    root = project.root / "data"
    add_legacy_assets(root, brick="0401p012")
    add_legacy_assets(root, brick="0400p012")
    import_legacy_fits(project.root, root)
    idx = read_dataset_index(project.root)
    source = idx["sources"][0]
    images = [a for a in idx["assets"] if a["product_kind"] == "image"]
    assert len(images) == 2
    decision = compatibility_policy([source], images, mode="mosaic")
    assert decision["allowed"] is True
    assert "Geometry plan only" in decision["reason"]
    assert compatibility_policy([source], images, mode="joint")["allowed"] is False


def test_science_cache_refuses_synthetic_unverifiable_asset_hash_labels():
    """A 64-character placeholder is not a canonical SHA-256 receipt."""
    import hashlib
    from lsbpipe.surveys.dataset_registry import science_cache_key
    actual = hashlib.sha256(b"immutable-source").hexdigest()
    asset = {"asset_id":"dataset-A","sha256":actual,
             "checksum_status":CHECKSUM_VERIFIED}
    kwargs = dict(stage="sky",algorithm_version="v1",
                  parameters={"plane":False},selected_assets=[asset])
    assert len(science_cache_key(**kwargs)) == 64
    for invalid in ("z"*64, actual.upper(), actual[:-1], "a"*65):
        with pytest.raises(ValueError, match="Unverified FITS fingerprint"):
            science_cache_key(
                **dict(kwargs,selected_assets=[dict(asset,sha256=invalid)]))
    with pytest.raises(ValueError, match="nonempty asset_id"):
        science_cache_key(
            **dict(kwargs,selected_assets=[dict(asset,asset_id="")]))
    with pytest.raises(ValueError, match="Stage, algorithm"):
        science_cache_key(**dict(kwargs,stage=42))
    with pytest.raises(ValueError, match="Stage, algorithm"):
        science_cache_key(**dict(kwargs,algorithm_version=False))
    with pytest.raises(ValueError, match="unique named"):
        science_cache_key(**kwargs, dependency_keys=["mask-v2","mask-v2"])


def test_registry_never_reads_external_symlinked_index_or_directory(
        project, tmp_path):
    """Project data lineage cannot be sourced from symlinked index metadata."""
    import_legacy_fits(project.root, project.root / "data")
    index = registry_path(project.root)
    external = tmp_path / "outside-dataset-index.json"
    index.rename(external)
    try:
        index.symlink_to(external)
    except (OSError, NotImplementedError):
        pytest.skip("Filesystem does not allow creating symlinks")
    assert index.is_symlink()
    with pytest.raises(ValueError, match="Unsafe symlinked dataset registry"):
        read_dataset_index(project.root)
    # Restoring the normal index must not rewrite the external record.
    index.unlink()
    external.replace(index)
    original = index.read_bytes()
    holder = tmp_path / "external-dataset-registry"
    holder.mkdir()
    moved = holder / "index.json"
    index.rename(moved)
    folder = index.parent
    folder.rmdir()
    try:
        folder.symlink_to(holder, target_is_directory=True)
    except (OSError, NotImplementedError):
        pytest.skip("Filesystem does not allow directory symlinks")
    assert folder.is_symlink()
    with pytest.raises(ValueError, match="Unsafe symlinked dataset registry"):
        read_dataset_index(project.root)
    from lsbpipe.surveys import dataset_registry
    with pytest.raises(ValueError, match="Unsafe symlinked dataset registry"):
        dataset_registry._save(
            project.root,
            {"schema":SCHEMA,"project_id":project.project_id})
    assert moved.read_bytes() == original
