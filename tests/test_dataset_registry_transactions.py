"""Independent process and conflict tests for dataset-index transactions.

These use temporary text placeholders with FITS filenames, never source FITS,
download services, adopted science or Qt. SHA is of exact placeholder bytes.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import pytest

from lsbpipe.surveys.dataset_registry import (
    CHECKSUM_VERIFIED, read_dataset_index, register_external_fits,
    registry_path, verify_asset_checksum,
)
from lsbpipe.surveys.object_projects import create_object_project


@pytest.fixture
def temporary_project(tmp_path, monkeypatch):
    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path / "session"))
    return create_object_project("Dataset race", tmp_path / "galaxy")


def _start(script, *args):
    return subprocess.Popen(
        [sys.executable, "-u", "-c", script, *map(str, args)],
        env=dict(os.environ), stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, text=True)


def _await_marker(marker: Path, process: subprocess.Popen):
    deadline = time.monotonic() + 18.
    while not marker.is_file():
        if process.poll() is not None:
            stdout, stderr = process.communicate(timeout=2)
            pytest.fail(
                "Child exited before checkpoint: " + repr(
                    (process.returncode, stdout, stderr)))
        if time.monotonic() > deadline:
            process.kill()
            out, err = process.communicate(timeout=8)
            pytest.fail(
                "Child timed out before checkpoint: " + repr((out, err)))
        time.sleep(.02)


def _finish(process):
    output, error = process.communicate(timeout=18)
    assert process.returncode == 0, (output, error)


def test_two_process_external_fits_imports_do_not_lose_asset_records(
        temporary_project, tmp_path):
    """One writer pauses after reading an old index; another races to save."""
    root = temporary_project.root
    first = tmp_path / "hsc-original.fits"
    second = tmp_path / "decals-original.fits"
    first.write_bytes(b"first original preserved")
    second.write_bytes(b"second original preserved")
    signal = tmp_path / "first-read"
    release = tmp_path / "go"
    ready = tmp_path / "second-ready"
    script = (
        "import pathlib,sys,time\n"
        "import lsbpipe.surveys.dataset_registry as d\n"
        "mode=sys.argv[1]\n"
        "project,asset,flag=map(pathlib.Path,sys.argv[2:])\n"
        "if mode=='pause':\n"
        "    original=d.read_dataset_index\n"
        "    def delayed(root):\n"
        "        snapshot=original(root)\n"
        "        flag.write_text('READ')\n"
        "        time.sleep(1.0)\n"
        "        return snapshot\n"
        "    d.read_dataset_index=delayed\n"
        "else:\n"
        "    flag.write_text('READY')\n"
        "    while not pathlib.Path(sys.argv[5]).is_file(): time.sleep(.01)\n"
        "d.register_external_fits(project,asset,provider='HSC',band='i')\n"
    )
    # The 'wait' worker is fully imported and ready BEFORE the
    # paused worker reaches its stale registry read.
    b = _start(script, "wait", root, second, ready, release)
    a = None
    try:
        _await_marker(ready, b)
        a = _start(script, "pause", root, first, signal)
        _await_marker(signal, a)
        release.write_text("GO", encoding="utf-8")
        _finish(a)
        _finish(b)
        index = read_dataset_index(root)
        assert len(index["sources"]) == len(index["assets"]) == 2
        assert {s["provider"] for s in index["sources"]} == {"HSC"}
        assert {a["path_ref"]["absolute"] for a in index["assets"]} == {
            str(first), str(second)}
        assert all(a["sha256"] is None for a in index["assets"])
        assert first.read_bytes() == b"first original preserved"
        assert second.read_bytes() == b"second original preserved"
    finally:
        for child in (a, b):
            if child is not None and child.poll() is None:
                child.kill()
                child.communicate(timeout=8)


def test_sha_verification_reloads_index_after_concurrent_import(
        temporary_project, tmp_path):
    """A long file hash must not overwrite a new import in the registry."""
    root = temporary_project.root
    image = tmp_path / "original-r.fits"
    added = tmp_path / "other-i.fits"
    image.write_bytes(b"exact optical bytes" * 4096)
    added.write_bytes(b"second original")
    item = register_external_fits(root, image, provider="HSC", band="r")
    sha_marker = tmp_path / "sha-started"
    ready = tmp_path / "writer-ready"
    release = tmp_path / "writer-go"
    verify_script = (
        "import pathlib,sys,time\n"
        "import lsbpipe.surveys.dataset_registry as d\n"
        "root,asset,marker=map(pathlib.Path,sys.argv[1:])\n"
        "original=d.hashlib.sha256\n"
        "class Tracked:\n"
        "    def __init__(self): self.h=original();self.first=True\n"
        "    def update(self,chunk):\n"
        "        if self.first:\n"
        "            self.first=False\n"
        "            marker.write_text('HASH')\n"
        "            time.sleep(1.0)\n"
        "        self.h.update(chunk)\n"
        "    def hexdigest(self): return self.h.hexdigest()\n"
        "d.hashlib.sha256=lambda *a,**k: Tracked()\n"
        "d.verify_asset_checksum(root,str(asset))\n"
    )
    writer_script = (
        "import pathlib,sys,time\n"
        "import lsbpipe.surveys.dataset_registry as d\n"
        "root,asset,ready,go=map(pathlib.Path,sys.argv[1:])\n"
        "ready.write_text('READY')\n"
        "while not go.is_file(): time.sleep(.01)\n"
        "d.register_external_fits(root,asset,provider='HSC',band='i')\n"
    )
    writer = _start(writer_script, root, added, ready, release)
    worker = None
    try:
        _await_marker(ready, writer)
        worker = _start(verify_script, root, item["asset_id"], sha_marker)
        _await_marker(sha_marker, worker)
        release.write_text("GO", encoding="utf-8")
        _finish(writer)
        _finish(worker)
        index = read_dataset_index(root)
        assert len(index["sources"]) == len(index["assets"]) == 2
        initial = next(a for a in index["assets"]
                       if a["asset_id"] == item["asset_id"])
        assert initial["checksum_status"] == CHECKSUM_VERIFIED
        assert initial["sha256"] == hashlib.sha256(image.read_bytes()).hexdigest()
        newcomer = next(a for a in index["assets"]
                        if a["asset_id"] != item["asset_id"])
        assert newcomer["sha256"] is None
        assert added.read_bytes() == b"second original"
    finally:
        for child in (worker, writer):
            if child is not None and child.poll() is None:
                child.kill()
                child.communicate(timeout=8)


def test_sha_commit_rejects_concurrent_asset_identity_change(
        temporary_project, tmp_path, monkeypatch):
    """Never publish an old hash over a later import metadata generation."""
    import lsbpipe.surveys.dataset_registry as module
    from lsbpipe.surveys.local_project import _atomic_json

    file = tmp_path / "original.fits"
    file.write_bytes(b"unchanged FITS payload")
    record = register_external_fits(
        temporary_project.root, file, provider="HSC", band="g")
    index_path = registry_path(temporary_project.root)
    original_digest = module.hashlib.sha256
    changed = []
    class TamperDigest:
        def __init__(self):
            self.hash = original_digest()
        def update(self, data):
            self.hash.update(data)
            if not changed:
                payload = json.loads(index_path.read_text(encoding="utf-8"))
                payload["assets"][0]["mtime_ns"] += 1
                _atomic_json(index_path, payload)
                changed.append(True)
        def hexdigest(self):
            return self.hash.hexdigest()
    monkeypatch.setattr(
        module.hashlib, "sha256", lambda *a, **kw: TamperDigest())
    with pytest.raises(RuntimeError, match="metadata changed"):
        verify_asset_checksum(temporary_project.root, record["asset_id"])
    after = read_dataset_index(temporary_project.root)["assets"][0]
    assert after["sha256"] is None
    assert after["checksum_status"] != CHECKSUM_VERIFIED
    assert file.read_bytes() == b"unchanged FITS payload"
    assert changed


def test_dataset_index_lock_rejects_symlink_without_touching_original(
        temporary_project, tmp_path):
    file = tmp_path / "original.fits"
    file.write_bytes(b"unmodified source")
    index = registry_path(temporary_project.root)
    index.parent.mkdir(parents=True, exist_ok=True)
    lock = index.with_name("index.json.lock")
    outsider = tmp_path / "external-metadata"
    outsider.write_bytes(b"untouched")
    try:
        lock.symlink_to(outsider)
    except (OSError, NotImplementedError):
        pytest.skip("Filesystem does not allow symlinks")
    with pytest.raises(ValueError, match="Unsafe metadata lock"):
        register_external_fits(
            temporary_project.root, file, provider="HSC", band="i")
    assert outsider.read_bytes() == b"untouched"
    assert file.read_bytes() == b"unmodified source"
    assert not index.exists()
