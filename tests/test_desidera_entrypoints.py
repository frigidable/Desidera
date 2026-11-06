"""Headless regression tests for DESIDERA's public compatibility layer."""
from __future__ import annotations

import sys
from types import ModuleType, SimpleNamespace


def test_public_facade_preserves_legacy_objects_and_version():
    import desidera
    import lsbpipe

    assert desidera.__version__ == lsbpipe.__version__
    assert set(lsbpipe.__all__).issubset(desidera.__all__)
    for name in lsbpipe.__all__:
        assert getattr(desidera, name) is getattr(lsbpipe, name)


def test_cli_is_the_existing_legacy_implementation():
    from desidera.cli import main as desidera_cli
    from lsbpipe.cli import main as legacy_cli

    assert desidera_cli is legacy_cli


def test_gui_package_can_be_imported_without_optional_qt():
    import desidera.gui as gui

    assert callable(gui._configure_webengine_for_launch)
    assert callable(gui.__getattr__)


def test_gui_resolves_exact_legacy_main_without_importing_qt(monkeypatch):
    import desidera.gui as gui

    stub = ModuleType("lsbpipe.gui.project_hub")
    stub.main = lambda: 37
    monkeypatch.setitem(sys.modules, "lsbpipe.gui.project_hub", stub)
    monkeypatch.delenv("DESIDERA_WEBENGINE_SOFTWARE", raising=False)
    assert gui.main is stub.main
    assert gui.main() == 37


def test_gui_unknown_attribute_is_not_a_launcher():
    import pytest
    import desidera.gui as gui

    with pytest.raises(AttributeError):
        gui.__getattr__("unknown")


def test_windows_webengine_diagnostic_opt_in_is_idempotent(monkeypatch):
    import desidera.gui as gui

    env = {
        "DESIDERA_WEBENGINE_SOFTWARE": "1",
        "QTWEBENGINE_CHROMIUM_FLAGS": "--single-process",
    }
    monkeypatch.setattr(gui, "os", SimpleNamespace(name="nt", environ=env))
    gui._configure_webengine_for_launch()
    gui._configure_webengine_for_launch()
    tokens = env["QTWEBENGINE_CHROMIUM_FLAGS"].split()
    assert tokens.count("--single-process") == 1
    assert tokens.count("--disable-gpu") == 1
    assert tokens.count("--disable-gpu-compositing") == 1


def test_non_windows_launcher_does_not_modify_webengine(monkeypatch):
    import desidera.gui as gui

    env = {"DESIDERA_WEBENGINE_SOFTWARE": "1"}
    monkeypatch.setattr(gui, "os", SimpleNamespace(name="posix", environ=env))
    gui._configure_webengine_for_launch()
    assert "QTWEBENGINE_CHROMIUM_FLAGS" not in env
