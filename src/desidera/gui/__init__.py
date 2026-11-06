"""DESIDERA desktop launcher with optional pre-Qt Windows diagnostics."""
from __future__ import annotations

import os


def _configure_webengine_for_launch() -> None:
    """Set opt-in WebEngine flags before importing any PySide6 widgets."""
    if os.name != "nt" or os.environ.get("DESIDERA_WEBENGINE_SOFTWARE") != "1":
        return
    flags = os.environ.get("QTWEBENGINE_CHROMIUM_FLAGS", "")
    for option in ("--disable-gpu", "--disable-gpu-compositing"):
        if option not in flags.split():
            flags = f"{flags} {option}".strip()
    os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] = flags


def __getattr__(name: str):
    """Expose the exact legacy entry point, but configure Qt before loading it.

    Both `python -m desidera.gui` and the `desidera-gui` console script
    resolve `main` through this hook. Ordinary imports of the GUI package
    remain safe when optional PySide6 dependencies are not installed.
    """
    if name != "main":
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    _configure_webengine_for_launch()
    from lsbpipe.gui.project_hub import main
    return main
