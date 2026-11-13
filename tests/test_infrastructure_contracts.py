"""Repository infrastructure and installation contracts (no network).

These tests deliberately inspect the actual repository configuration that
would run in GitHub Actions. They cannot establish branch protections or
external service availability and do not replace Actions validation.
"""
from __future__ import annotations

from importlib import metadata
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
WORKFLOWS = ROOT / ".github" / "workflows"


def _workflow(name: str) -> str:
    return (WORKFLOWS / name).read_text(encoding="utf-8")


def test_desidera_installed_package_exposes_both_apis_and_gui_test_extra():
    """The declared test tool belongs to the installable distribution."""
    import desidera
    import lsbpipe

    dist = metadata.distribution("desidera")
    assert dist.version == desidera.__version__ == lsbpipe.__version__
    provided = set(dist.metadata.get_all("Provides-Extra") or ())
    assert {"test", "test-gui", "exact", "gui"} <= provided
    requirements = dist.requires or ()
    assert any(
        "pytest-qt" in item.lower() and "test-gui" in item.lower()
        for item in requirements)


def test_science_ci_observes_source_tests_scripts_and_builds_a_real_wheel():
    """An editable-only test matrix is not sufficient to vet packaging."""
    text = _workflow("survey-tests.yml")
    assert re.search(r"(?m)^  pull_request:\s*$", text)
    assert all(token in text for token in (
        '"src/**"', '"tests/**"', '"scripts/**"',
        '"pyproject.toml"', 'python -m pip check',
        'package-wheel:', 'pip wheel --no-deps',
        'desidera-wheels/desidera-*.whl',
        'INSTALLED_WHEEL_IMPORT_AND_ENTRYPOINTS_OK',
        'timeout-minutes: 20'))
    assert 'cancel-in-progress: true' in text
    assert 'permissions:\n  contents: read' in text


def test_gui_ci_covers_core_changes_and_cancels_stale_pr_runs():
    """Background/scene edits must not silently bypass the full GUI suite."""
    text = _workflow("ui-baseline.yml")
    assert all(token in text for token in (
        '"src/**"', '"tests/**"', '"scripts/**"',
        '"pyproject.toml"', '"docs/ux/**"',
        '"[exact,gui,test,test-gui]"',
        'python -m pip check',
        'pytest -q tests/gui',
        'cancel-in-progress: true',
        'QT_QPA_PLATFORM: offscreen',
        'actions/upload-artifact@v4'))
    assert 'permissions:\n  contents: read' in text


def test_actual_archive_network_is_manual_only_but_offline_pr_tests_run():
    """Network-dependent original-FITS acceptance cannot be an implicit PR step."""
    text = _workflow("legacy-ccd-real-archive.yml")
    assert 'workflow_dispatch:' in text
    assert 'pull_request:' in text
    assert 'python -m pytest -q tests/test_real_legacy_ccd_acceptance.py' in text
    assert re.search(
        r"- name: Fetch official DR9 DR11 original FITS and audit\n"
        r"\s*if: github.event_name == 'workflow_dispatch'\n", text)
    assert re.search(
        r"- name: Print verified source hashes and cross-release overlap\n"
        r"\s*if: github.event_name == 'workflow_dispatch'\n", text)
    assert "if: always() && github.event_name == 'workflow_dispatch'" in text
