"""Immutable human mask drafts retain per-object origin and history."""
import json

import pytest

from lsbpipe.photometry.mask_drafts import (
    read_active_mask_draft, save_mask_draft)
from lsbpipe.photometry.sky_masks import (
    SkyMaskEdit, SkyMaskShape)
from lsbpipe.surveys.object_projects import create_object_project


def test_explicit_save_reload_undo_and_concurrent_head_guard(tmp_path, monkeypatch):
    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path / "session"))
    galaxy = create_object_project("UGC 2162", tmp_path / "ugc")
    draft, digest = read_active_mask_draft(galaxy.root)
    assert digest is None
    assert not (galaxy.root / "analyses").exists()
    circle = SkyMaskShape("circle", (40.096208, 1.229250), radius_arcsec=8.)
    edited = draft.with_edit(SkyMaskEdit("manual-star-1", "exclude", circle))
    result_hash = save_mask_draft(galaxy.root, edited, expected_digest=digest)
    loaded, current = read_active_mask_draft(galaxy.root)
    assert loaded == edited
    assert current == result_hash == edited.digest
    old = galaxy.root / "analyses" / "mask-drafts" / ("draft-" + current + ".json")
    bytes_before = old.read_bytes()
    clean = edited.without_edit("manual-star-1")
    with pytest.raises(ValueError, match="modified elsewhere"):
        save_mask_draft(galaxy.root, clean, expected_digest=None)
    new_hash = save_mask_draft(
        galaxy.root, clean, expected_digest=current)
    assert new_hash != current
    assert old.read_bytes() == bytes_before
    assert not (galaxy.root / "qa-runs").exists()


def test_corrupt_snapshot_is_refused_instead_of_invisibly_reset(
        tmp_path, monkeypatch):
    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path / "session"))
    galaxy = create_object_project("Draft", tmp_path / "draft")
    initial, _ = read_active_mask_draft(galaxy.root)
    digest = save_mask_draft(galaxy.root, initial, expected_digest=None)
    file = galaxy.root / "analyses" / "mask-drafts" / f"draft-{digest}.json"
    content = json.loads(file.read_text("utf8"))
    content["project_id"] = "another-project"
    file.write_text(json.dumps(content), encoding="utf8")
    with pytest.raises(ValueError, match="Modified/cross-project"):
        read_active_mask_draft(galaxy.root)


def test_cross_project_save_refused(tmp_path, monkeypatch):
    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path / "session"))
    a = create_object_project("GalaxyA", tmp_path / "a")
    b = create_object_project("GalaxyB", tmp_path / "b")
    draft, _ = read_active_mask_draft(a.root)
    with pytest.raises(ValueError, match="different galaxy"):
        save_mask_draft(b.root, draft, expected_digest=None)
    assert not (b.root / "analyses").exists()
