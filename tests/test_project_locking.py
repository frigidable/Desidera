"""Real OS-process regression for crash-safe local Desidera library locks.

These tests operate on temporary metadata only. No original FITS, GUI,
network, user profile or scientific products are touched.
"""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys

import pytest

from lsbpipe.surveys.object_projects import (
    REGISTRY_SCHEMA, create_object_project, library_index_path,
    read_object_project, recent_projects,
)
from lsbpipe.surveys.project_locking import exclusive_metadata_lock


def _child(script: str, *args: str, env: dict | None = None):
    return subprocess.Popen(
        [sys.executable, "-u", "-c", script, *map(str, args)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True, env=env)


def test_lock_is_stable_and_usable_after_nested_failure(tmp_path):
    lock_path = tmp_path / "projects.json.lock"
    with pytest.raises(RuntimeError, match="simulated"):
        with exclusive_metadata_lock(lock_path):
            assert lock_path.is_file()
            raise RuntimeError("simulated crash in transaction")
    before = lock_path.stat()
    with exclusive_metadata_lock(lock_path):
        assert lock_path.is_file()
    after = lock_path.stat()
    assert (before.st_dev, before.st_ino) == (after.st_dev, after.st_ino)
    assert lock_path.read_bytes() in (b"", b"0")


def test_os_process_holds_lock_and_other_writer_times_out(tmp_path):
    """A separate interpreter must not enter another process's critical section."""
    lock_path = tmp_path / "projects.json.lock"
    script = (
        "import sys, time\n"
        "from lsbpipe.surveys.project_locking import exclusive_metadata_lock\n"
        "with exclusive_metadata_lock(sys.argv[1]):\n"
        "    print('HELD', flush=True)\n"
        "    time.sleep(0.7)\n"
    )
    child = _child(script, lock_path)
    try:
        assert child.stdout.readline().strip() == "HELD"
        with pytest.raises(TimeoutError, match="Metadata lock busy"):
            with exclusive_metadata_lock(lock_path, timeout_seconds=0.08):
                pytest.fail("second process incorrectly entered lock")
        stdout, stderr = child.communicate(timeout=12)
        assert child.returncode == 0, (stdout, stderr)
        with exclusive_metadata_lock(lock_path, timeout_seconds=0.5):
            pass
    finally:
        if child.poll() is None:
            child.kill()
            child.communicate(timeout=12)


def test_process_exit_without_unlock_releases_os_lock(tmp_path):
    """A killed/crashed instance cannot strand a permanent lock claim."""
    lock_path = tmp_path / "projects.json.lock"
    script = (
        "import os, sys\n"
        "from lsbpipe.surveys.project_locking import exclusive_metadata_lock\n"
        "with exclusive_metadata_lock(sys.argv[1]):\n"
        "    print('HELD', flush=True)\n"
        "    os._exit(0)\n"
    )
    child = _child(script, lock_path)
    stdout, stderr = child.communicate(timeout=12)
    assert child.returncode == 0 and stdout.strip() == "HELD", stderr
    with exclusive_metadata_lock(lock_path, timeout_seconds=1):
        assert lock_path.is_file()


def test_symlinked_lock_is_rejected_without_changing_external_bytes(tmp_path):
    lock_path = tmp_path / "projects.json.lock"
    outside = tmp_path / "original-other-file"
    outside.write_bytes(b"external data unchanged")
    try:
        lock_path.symlink_to(outside)
    except (OSError, NotImplementedError):
        pytest.skip("Operating system does not permit symlinks")
    with pytest.raises(ValueError, match="Unsafe metadata lock"):
        with exclusive_metadata_lock(lock_path):
            pytest.fail("symlinked metadata lock accepted")
    assert outside.read_bytes() == b"external data unchanged"


def test_simultaneous_registry_updates_keep_both_project_identities(
        tmp_path, monkeypatch):
    """Force one reader to pause while another process attempts update.

    Without an OS lock, the second process publishes a competing snapshot;
    the paused writer then loses the second registry entry. The final
    registry must retain both independently created UUIDs.
    """
    state = tmp_path / "session"
    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(state))
    first = create_object_project("Galaxy A", tmp_path / "galaxyA")
    second = create_object_project("Galaxy B", tmp_path / "galaxyB")
    registry_path = library_index_path()
    registry_path.write_text(
        '{"schema":"' + REGISTRY_SCHEMA +
        '","projects":[],"dismissed_legacy_roots":[]}',
        encoding="utf-8")
    env = dict(os.environ)
    script = (
        "import sys,time\n"
        "import lsbpipe.surveys.object_projects as module\n"
        "p=module.read_object_project(sys.argv[1])\n"
        "if sys.argv[2]=='pause':\n"
        "    old=module._registry_payload\n"
        "    def paused():\n"
        "        result=old()\n"
        "        print('READ', flush=True)\n"
        "        time.sleep(.55)\n"
        "        return result\n"
        "    module._registry_payload=paused\n"
        "module.register_object_project(p)\n"
        "print('SAVED', flush=True)\n"
    )
    a = _child(script, first.root, "pause", env=env)
    b = None
    try:
        assert a.stdout.readline().strip() == "READ"
        b = _child(script, second.root, "normal", env=env)
        out_b, err_b = b.communicate(timeout=15)
        out_a, err_a = a.communicate(timeout=15)
        assert a.returncode == 0, (out_a, err_a)
        assert b.returncode == 0, (out_b, err_b)
        assert "SAVED" in out_b
        # Both processes used the SAME metadata index and stable lock.
        recents = recent_projects()
        assert {row["project_id"] for row in recents} == {
            first.project_id, second.project_id}
        assert all(row["available"] for row in recents)
        assert read_object_project(first.root).project_id == first.project_id
        assert read_object_project(second.root).project_id == second.project_id
        assert registry_path.with_name("projects.json.lock").is_file()
    finally:
        for process in (a, b):
            if process is not None and process.poll() is None:
                process.kill()
                process.communicate(timeout=12)
