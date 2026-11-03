"""Launch DESIDERA Project Hub.

Set DESIDERA_WEBENGINE_SOFTWARE=1 before launch to test whether an
overlaid Windows artifact is caused by the Qt WebEngine GPU compositor.
This is an opt-in diagnostic, not a default performance regression.
"""
import os

if os.name == "nt" and os.environ.get("DESIDERA_WEBENGINE_SOFTWARE") == "1":
    flags = os.environ.get("QTWEBENGINE_CHROMIUM_FLAGS", "")
    os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] = (
        flags + " --disable-gpu --disable-gpu-compositing").strip()

from lsbpipe.gui.project_hub import main

if __name__ == "__main__":
    raise SystemExit(main())
