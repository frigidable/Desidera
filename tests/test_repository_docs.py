"""Protect public documentation links and current maintenance branch policy."""
from __future__ import annotations

from pathlib import Path
import re
from urllib.parse import unquote


ROOT = Path(__file__).resolve().parents[1]
LANDING_DOCS = (
    ROOT / "README.md",
    ROOT / "docs" / "README.md",
)


def test_public_documentation_links_point_to_existing_tracked_files():
    broken = []
    pattern = re.compile(r"\]\(([^)]+)\)")
    for doc in LANDING_DOCS:
        for raw in pattern.findall(doc.read_text(encoding="utf-8")):
            href = raw.split("#", 1)[0].split("?", 1)[0]
            if not href or re.match(r"^[a-zA-Z][a-zA-Z0-9+.-]*:", href):
                continue
            target = (doc.parent / unquote(href)).resolve()
            if not target.is_relative_to(ROOT) or not target.exists():
                broken.append(f"{doc.relative_to(ROOT)}: {raw}")
    assert not broken, "Broken landing-page links: " + "; ".join(broken)


def test_active_branch_instructions_do_not_point_to_historical_ui_branch():
    guide = (ROOT / "docs" / "ux" / "AUTONOMY_PLAYBOOK.md").read_text(
        encoding="utf-8")
    assert "feature/survey-explorer-v04" not in guide
    assert "one short-lived branch" in guide
