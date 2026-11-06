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


def main() -> int:
    """Launch the existing Project Hub without importing Qt prematurely."""
    _configure_webengine_for_launch()
    from lsbpipe.gui.project_hub import main as launch_project_hub
    return launch_project_hub()
