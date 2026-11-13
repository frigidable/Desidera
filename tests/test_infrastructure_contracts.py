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
        '"pyproject.toml"', '"README.md"', '"AGENTS.md"',
        '"docs/architecture/**"', '".github/workflows/**"',
        'python -m pip check',
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
        '"pyproject.toml"', '"README.md"', '"AGENTS.md"',
        '"docs/architecture/**"', '"docs/ux/**"',
        '".github/workflows/**"',
        '".[exact,gui,test,test-gui]"',
        'python -m pip check',
        'pytest -q tests/gui',
        'cancel-in-progress: true',
        'QT_QPA_PLATFORM: offscreen',
        'actions/upload-artifact@ea165f8d65b6e75b540449e92b4886f43607fa02'))
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



def test_all_external_github_actions_are_sha_pinned():
    """Floating tags must not silently change the Actions runtime."""
    for name in (
            "survey-tests.yml", "ui-baseline.yml",
            "legacy-ccd-real-archive.yml"):
        source = _workflow(name)
        uses = re.findall(
            r"(?m)^\s*- uses:\s+([^\s#]+)", source)
        assert uses, name
        for action in uses:
            assert re.fullmatch(
                r"[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+@[0-9a-f]{40}",
                action), (name, action)
        assert "actions/checkout@" in source
        assert "actions/setup-python@" in source



def test_scientific_and_survey_core_do_not_import_gui_or_qt():
    """Enforce the source-level GUI/science dependency inversion.

    Parse the real Python files rather than import Qt in the test runner;
    this remains valid in minimal headless installations and tests direct
    import regressions even on otherwise untested CLI paths.
    """
    import ast

    core = ROOT / "src" / "lsbpipe"
    paths = [
        core / "background.py",
        core / "scene.py",
        core / "survey.py",
        core / "mask_engine.py",
        *sorted((core / "photometry").rglob("*.py")),
        *sorted((core / "surveys").rglob("*.py")),
    ]
    forbidden = (
        "PySide6", "PyQt5", "PyQt6", "qtpy", "lsbpipe.gui",
        "desidera.gui",
    )
    assert len(paths) >= 55
    for path in paths:
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                names = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom):
                names = [node.module or ""]
                # Relative imports inside the science package remain legal;
                # Python resolves these against the local science module.
                if node.level:
                    continue
            else:
                continue
            for name in names:
                assert not any(
                    name == prefix or name.startswith(prefix + ".")
                    for prefix in forbidden
                ), (path.relative_to(ROOT), node.lineno, name)
