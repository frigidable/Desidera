"""Never let automated tests read/write a user's real recent-project pointer."""
import pytest


@pytest.fixture(autouse=True)
def isolate_local_project_settings(tmp_path, monkeypatch):
    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path / "session-config"))
    # GUI tests opt into auto-fetch with mocked transports explicitly.
    # Keep older offline fixtures isolated from real remote JPEG requests.
    monkeypatch.setenv("DESIDERA_DISABLE_AUTOMATIC_RGB", "1")
