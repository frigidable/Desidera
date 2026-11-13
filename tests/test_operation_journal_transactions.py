"""Cross-process regression for journal state-machine terminal safety."""
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import time

import pytest

from lsbpipe.surveys.object_projects import create_object_project
from lsbpipe.surveys.operation_journal import (
    finish_operation, interrupted_operations, request_stop, start_operation,
)


def _start_child(script, *args):
    return subprocess.Popen(
        [sys.executable, "-u", "-c", script, *map(str, args)],
        env=dict(os.environ), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True)


def _await_file(path: Path, process: subprocess.Popen):
    deadline = time.monotonic() + 18.
    while not path.is_file():
        if process.poll() is not None:
            out, err = process.communicate(timeout=3)
            pytest.fail("Child exited before rendezvous: " +
                        repr((process.returncode, out, err)))
        if time.monotonic() > deadline:
            process.kill()
            out, err = process.communicate(timeout=8)
            pytest.fail("Child stalled before rendezvous: " + repr((out, err)))
        time.sleep(.02)


def test_stop_finish_race_cannot_resurrect_terminal_operation(
        tmp_path, monkeypatch):
    """A stale STOP_REQUESTED writer must never overwrite COMPLETED."""
    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path / "user"))
    project = create_object_project("Journal race", tmp_path / "galaxy")
    journal = start_operation(
        project.root, project.project_id, "dataset_audit")
    marker = tmp_path / "stop-has-read-running"
    ready = tmp_path / "finish-ready"
    go = tmp_path / "finish-go"
    stopped_script = (
        "import pathlib,sys,time\n"
        "import lsbpipe.surveys.operation_journal as j\n"
        "path=pathlib.Path(sys.argv[1]);marker=pathlib.Path(sys.argv[2])\n"
        "original=j._read\n"
        "def delayed(path):\n"
        "    state=original(path)\n"
        "    marker.write_text('READ')\n"
        "    time.sleep(1.0)\n"
        "    return state\n"
        "j._read=delayed\n"
        "j.request_stop(path)\n"
    )
    finish_script = (
        "import pathlib,sys,time\n"
        "import lsbpipe.surveys.operation_journal as j\n"
        "path,ready,go=map(pathlib.Path,sys.argv[1:])\n"
        "ready.write_text('READY')\n"
        "while not go.is_file(): time.sleep(.01)\n"
        "j.finish_operation(path,'COMPLETED')\n"
    )
    finishing = _start_child(finish_script, journal, ready, go)
    stopping = None
    try:
        _await_file(ready, finishing)
        stopping = _start_child(stopped_script, journal, marker)
        _await_file(marker, stopping)
        go.write_text("GO", encoding="utf-8")
        for p in (stopping, finishing):
            out, err = p.communicate(timeout=18)
            assert p.returncode == 0, (out, err)
        state = json.loads(journal.read_text(encoding="utf-8"))
        assert state["state"] == "COMPLETED"
        assert state["scientific_status"] == "NOT_SCIENTIFIC_EVIDENCE"
        assert interrupted_operations(project.root) == []
        before = journal.read_bytes()
        request_stop(journal)
        assert journal.read_bytes() == before
        with pytest.raises(ValueError, match="immutable"):
            finish_operation(journal, "FAILED")
        assert (journal.parent / ".operation-journal.lock").is_file()
    finally:
        for p in (stopping, finishing):
            if p is not None and p.poll() is None:
                p.kill()
                p.communicate(timeout=8)


def test_journal_symlink_lock_fails_closed_without_rewriting_state(
        tmp_path, monkeypatch):
    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path / "user"))
    project = create_object_project("Journal lock", tmp_path / "galaxy")
    journal = start_operation(project.root, project.project_id, "download")
    original = journal.read_bytes()
    lock = journal.parent / ".operation-journal.lock"
    old_lock = tmp_path / "old-lock-inode"
    lock.rename(old_lock)
    outside = tmp_path / "external-original"
    outside.write_bytes(b"do not touch")
    try:
        lock.symlink_to(outside)
    except (OSError, NotImplementedError):
        pytest.skip("Filesystem does not permit symlinks")
    with pytest.raises(ValueError, match="Unsafe metadata lock"):
        request_stop(journal)
    with pytest.raises(ValueError, match="Unsafe metadata lock"):
        finish_operation(journal, "COMPLETED")
    assert journal.read_bytes() == original
    assert outside.read_bytes() == b"do not touch"
