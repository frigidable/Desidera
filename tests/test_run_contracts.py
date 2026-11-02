"""P0.3: immutable project-scoped QA execution and explicit staleness states."""
from __future__ import annotations

from dataclasses import replace
import json

import pytest

from lsbpipe.gui.qa_records import (
    ScienceContext, new_qa_run, original_quartet, snapshot_inputs)
from lsbpipe.gui.run_contracts import (
    canonical_digest, check_run_receipt, dependency_order,
    discover_run_history, input_content_digest, make_run_contract,
    science_parameters)
from lsbpipe.gui.science_workflows import run_science_qa
from lsbpipe.surveys.legacy import BrickHit, SkyField, product_for
from lsbpipe.surveys.object_projects import create_object_project

FIELD = SkyField(40.096208, 1.22925, 5.)
HIT = BrickHit("dr9", "south", "0401p012", 40, 40.2, 1.1, 1.4, {"r": 9})


def fixture_project(tmp_path):
    project = create_object_project(
        "UGC 2162 provenance", tmp_path / "UGC_2162", field=FIELD)
    root = project.root / "data"
    root.mkdir()
    for i, kind in enumerate(("image", "invvar", "nexp", "maskbits")):
        item = product_for(HIT, kind, None if kind == "maskbits" else "r")
        file = root / item.relative_path
        file.parent.mkdir(parents=True, exist_ok=True)
        file.write_bytes((kind + str(i)).encode() * 3)
    context = ScienceContext.capture(
        action="square", field=FIELD, hit=HIT, band="r",
        directory=project.root, project_id=project.project_id,
        target_radius_arcsec=35., sky_inner_arcsec=65.,
        sky_outer_arcsec=110.)
    inputs = snapshot_inputs({"source": original_quartet(root, HIT, "r")})
    return project, root, context, inputs


def test_semantic_parameters_ignore_irrelevant_qa_radii():
    field = FIELD
    first = ScienceContext.capture(
        action="square", field=field, hit=HIT, band="r", directory="C:/A")
    second = ScienceContext.capture(
        action="square", field=field, hit=HIT, band="r", directory="D:/B",
        target_radius_arcsec=40., sky_inner_arcsec=80.,
        sky_outer_arcsec=130.)
    assert science_parameters(first.to_dict()) == science_parameters(
        second.to_dict())
    third = ScienceContext.capture(
        action="census", field=field, hit=HIT, band="r", directory="C:/A")
    fourth = replace(third, target_radius_arcsec=40.)
    assert science_parameters(third.to_dict()) != science_parameters(
        fourth.to_dict())
    pair = replace(third, action="compare")
    assert science_parameters(pair.to_dict()) == science_parameters(
        replace(pair, target_radius_arcsec=999.).to_dict())


def test_content_fingerprint_is_portable_and_canonical(tmp_path):
    p, root, ctx, inputs = fixture_project(tmp_path)
    moved = json.loads(json.dumps(inputs))
    moved["source"]["image"]["path"] = "/foreign/new-folder/image.fits.fz"
    moved["source"]["image"]["mtime_ns"] = -1
    assert input_content_digest(moved) == input_content_digest(inputs)
    assert canonical_digest({"α": 2, "β": 1}) == canonical_digest(
        {"β": 1, "α": 2})
    a = make_run_contract(context=ctx.to_dict(), inputs=inputs,
                          project_id=p.project_id, code_revision="1"*40)
    b = make_run_contract(context=ctx.to_dict(), inputs=moved,
                          project_id=p.project_id, code_revision="1"*40)
    assert a["signature_sha256"] == b["signature_sha256"]
    assert a["dataset_ids"] == ["dr9/south/0401p012/r"]
    assert a["project_id"] == p.project_id
    assert a["status"] == "DIAGNOSTIC_NOT_ADOPTED"


def test_receipt_requires_explicit_current_input_snapshot(tmp_path):
    p, root, ctx, inputs = fixture_project(tmp_path)
    with new_qa_run(data_root=root, project_root=p.root,
                    context=ctx, inputs=inputs) as (stage, record):
        (stage / "square.json").write_text('{"valid": true}')
    manifest = p.root / "qa-runs" / record["run_id"] / "manifest.json"
    stored = json.loads(manifest.read_text(encoding="utf-8"))
    assert stored["context"]["project_id"] == p.project_id
    assert stored["run_contract"]["project_id"] == p.project_id
    assert stored["run_contract"]["parameter_sha256"]
    assert stored["run_contract"]["input_set_sha256"]
    assert (check_run_receipt(manifest, project_id=p.project_id)["state"] ==
            "RECORDED_NOT_REVALIDATED")
    check = check_run_receipt(
        manifest, project_id=p.project_id,
        current_context=ctx.to_dict(), current_inputs=inputs)
    assert check["state"] == "CURRENT_INPUTS_AND_PARAMETERS"
    assert check_run_receipt(
        manifest, project_id=str(__import__("uuid").uuid4())
    )["state"] == "PROJECT_MISMATCH"
    changed_field = replace(ctx, ra_deg=40.1012)
    assert check_run_receipt(
        manifest, current_context=changed_field.to_dict(),
        current_inputs=inputs)["state"] == "STALE_PARAMETERS"

    paths = original_quartet(root, HIT, "r")
    image = paths["image"]
    original = image.read_bytes()
    image.write_bytes(b"Z" * len(original))
    current = snapshot_inputs({"source": paths})
    assert check_run_receipt(
        manifest, current_context=ctx.to_dict(),
        current_inputs=current)["state"] == "STALE_INPUTS"
    assert len(discover_run_history(
        p.root, data_root=root, project_id=p.project_id)) == 1


def test_artifact_tampering_and_old_report_are_not_current(tmp_path):
    p, root, ctx, inputs = fixture_project(tmp_path)
    with new_qa_run(data_root=root, project_root=p.root,
                    context=ctx, inputs=inputs) as (stage, record):
        (stage / "square.json").write_text('{"valid": true}')
    manifest = p.root / "qa-runs" / record["run_id"] / "manifest.json"
    output = manifest.parent / "square.json"
    output.write_text('{"valid": null}')
    assert check_run_receipt(manifest)["state"] == "ARTIFACT_MODIFIED"
    # Historical v1 records without a science-run contract remain readable,
    # but may NOT be labeled validated against new source hashes.
    payload = json.loads(manifest.read_text())
    del payload["run_contract"]
    manifest.write_text(json.dumps(payload))
    assert check_run_receipt(manifest)["state"] == "ARTIFACT_MODIFIED"
    output.write_text('{"valid": true}')
    assert check_run_receipt(manifest)["state"] == "LEGACY_UNVERIFIED"
    assert check_run_receipt(
        manifest, current_context=ctx.to_dict(), current_inputs=inputs
    )["state"] == "LEGACY_UNVERIFIED"


def test_graph_rejects_missing_and_circular_dependencies():
    assert dependency_order({
        "square": [], "census": ["square"],
        "masks": ["square"], "photometry": ["census", "masks"],
    })[-1] == "photometry"
    with pytest.raises(ValueError, match="Missing scientific dependency"):
        dependency_order({"census": ["unknown"]})
    with pytest.raises(ValueError, match="cycle"):
        dependency_order({"a": ["b"], "b": ["a"]})


def test_dependency_must_be_existing_unchanged_same_project(tmp_path):
    p, root, ctx, inputs = fixture_project(tmp_path)
    with new_qa_run(data_root=root, project_root=p.root,
                    context=ctx, inputs=inputs) as (stage, prior):
        (stage / "square.json").write_text("signed")
    with new_qa_run(
            data_root=root, project_root=p.root,
            context=ctx, inputs=inputs,
            dependencies=(prior["run_id"],)) as (stage, derived):
        (stage / "square.json").write_text("derived")
    assert derived["run_contract"]["dependencies"] == [prior["run_id"]]
    (p.root / "qa-runs" / prior["run_id"] / "square.json").write_text("altered")
    with pytest.raises(ValueError, match="dependency"):
        with new_qa_run(
                data_root=root, project_root=p.root, context=ctx,
                inputs=inputs, dependencies=(prior["run_id"],)) as (stage, _):
            (stage / "square.json").write_text("unsafe")
    assert not list((p.root / "qa-runs").glob("*.partial"))
    with pytest.raises(ValueError, match="dependency"):
        with new_qa_run(data_root=root, project_root=p.root,
                        context=ctx, inputs=inputs,
                        dependencies=("20260101T000000000000Z-1234567890",)):
            pass


def test_run_science_qa_rejects_foreign_passport_before_io(tmp_path):
    p, root, ctx, inputs = fixture_project(tmp_path)
    another = create_object_project(
        "Another galaxy", tmp_path / "other", field=FIELD)
    with pytest.raises(ValueError, match="Project UUID mismatch"):
        run_science_qa(
            "square", field=FIELD, hit=HIT, band="r",
            directory=root, context=ctx, project_root=another.root)
    assert not (another.root / "qa-runs").exists()


def test_run_science_qa_publishes_in_project_not_external_data(tmp_path, monkeypatch):
    from lsbpipe.surveys import square_coverage
    p, data, ctx, inputs = fixture_project(tmp_path)
    data_root = tmp_path / "separate-datastore"
    data.rename(data_root)
    report = type("Report", (), {
        "pixel_centers_sampled": 10,
        "polygon_within_single_brick": True,
        "instrumentally_usable_fraction": 1.,
        "mean_nexp_where_positive": 2.,
        "maskbit_counts": {"SATUR_R": 0, "MEDIUM": 0, "GALAXY": 0},
    })()
    monkeypatch.setattr(square_coverage, "inspect_local_square",
                        lambda *a, **k: report)
    monkeypatch.setattr(square_coverage, "save_square_coverage",
                        lambda *a: a[-1].write_text('{"status":"diagnostic"}'))
    amended = replace(ctx, directory=str(data_root))
    output = run_science_qa(
        "square", field=FIELD, hit=HIT, band="r",
        directory=data_root, context=amended, project_root=p.root)
    receipt = p.root / "qa-runs" / output["run_id"] / "manifest.json"
    assert receipt.is_file()
    assert not (data_root / "qa-runs").exists()
    assert output["run_contract"]["project_id"] == p.project_id
    assert check_run_receipt(receipt, project_id=p.project_id)["state"] == (
        "RECORDED_NOT_REVALIDATED")
