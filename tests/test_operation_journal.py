"""Crash journal: never conflate attempted UI actions with scientific outputs."""
from __future__ import annotations

import json
import pytest

from lsbpipe.surveys.object_projects import create_object_project
from lsbpipe.surveys.operation_journal import (
    SCHEMA, finish_operation, interrupted_operations, request_stop,
    start_operation,
)


def test_uncertain_job_survives_reopening_without_mutating_original(tmp_path):
    project = create_object_project("UGC 2162", tmp_path / "UGC_2162")
    original = project.root / "notes.txt"
    original.write_bytes(b"existing original files untouched")
    marker = start_operation(project.root, project.project_id, "dataset_audit")
    assert marker.parent == project.root / "operations"
    before = marker.read_bytes()
    pending = interrupted_operations(project.root)
    assert len(pending) == 1
    assert pending[0]["schema"] == SCHEMA
    assert pending[0]["action"] == "dataset_audit"
    assert pending[0]["state"] == "RUNNING"
    assert pending[0]["scientific_status"] == "NOT_SCIENTIFIC_EVIDENCE"
    assert marker.read_bytes() == before
    assert original.read_bytes() == b"existing original files untouched"
    # New operations do not overwrite interrupted records; no replay.
    new = start_operation(project.root, project.project_id, "science_qa:square")
    assert new != marker
    assert len(interrupted_operations(project.root)) == 2
    finish_operation(new, "COMPLETED")
    assert len(interrupted_operations(project.root)) == 1
    assert marker.read_bytes() == before
    assert not (project.root / "qa-runs").exists()


def test_stop_request_and_terminal_status_are_atomic_and_explicit(tmp_path):
    project = create_object_project("Candidate", tmp_path / "candidate")
    marker = start_operation(project.root, project.project_id, "download")
    request_stop(marker)
    assert interrupted_operations(project.root)[0]["state"] == "STOP_REQUESTED"
    finish_operation(marker, "CANCELLED")
    assert interrupted_operations(project.root) == []
    state = json.loads(marker.read_text(encoding="utf-8"))
    assert state["state"] == "CANCELLED"
    assert state["scientific_status"] == "NOT_SCIENTIFIC_EVIDENCE"
    request_stop(marker)
    assert json.loads(marker.read_text()) == state
    with pytest.raises(ValueError, match="immutable"):
        finish_operation(marker, "COMPLETED")


@pytest.mark.parametrize("outcome", ["COMPLETED", "FAILED", "IGNORED"])
def test_terminal_state_never_claims_scientific_acceptance(tmp_path, outcome):
    project = create_object_project("Candidate", tmp_path / "candidate")
    marker = start_operation(project.root, project.project_id, "background")
    finish_operation(marker, outcome)
    assert interrupted_operations(project.root) == []
    record = json.loads(marker.read_text())
    assert record["state"] == outcome
    assert record["scientific_status"] == "NOT_SCIENTIFIC_EVIDENCE"


def test_identity_and_unsupported_actions_fail_closed(tmp_path):
    project = create_object_project("A", tmp_path / "A")
    other = create_object_project("B", tmp_path / "B")
    with pytest.raises(ValueError, match="UUID"):
        start_operation(project.root, other.project_id, "download")
    for action in ("", "../escape", "an operation", ".hidden"):
        with pytest.raises(ValueError, match="Invalid operation"):
            start_operation(project.root, project.project_id, action)
    assert interrupted_operations(project.root) == []


def test_corrupt_journal_is_ignored_without_auto_cleanup(tmp_path):
    project = create_object_project("A", tmp_path / "A")
    marker = start_operation(project.root, project.project_id, "background")
    marker.write_text("bad json", encoding="utf-8")
    assert interrupted_operations(project.root) == []
    assert marker.read_text() == "bad json"
    with pytest.raises((ValueError, json.JSONDecodeError)):
        finish_operation(marker, "FAILED")
