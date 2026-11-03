"""P0.4 real Qt lifecycle: cooperative stop and stale GUI callback isolation."""
from __future__ import annotations

from threading import Event
import pytest

pytest.importorskip("pytestqt")
pytest.importorskip("PySide6.QtWebEngineWidgets")

from PySide6.QtCore import QThread
from PySide6.QtWidgets import QApplication
from scripts.offline_gui_fixture import offline_legacy_viewer
from lsbpipe.gui.survey_explorer import SurveyExplorer
from lsbpipe.gui.task_execution import TaskToken, TaskCancelled


@pytest.fixture
def workspace(qtbot):
    with offline_legacy_viewer():
        window = SurveyExplorer(restore_session=False, launch_remote=False)
    qtbot.addWidget(window)
    window.show()
    yield window
    # A real Qt lifecycle must never destroy an executing QThread.
    if window.thread is not None:
        window._request_stop()
        qtbot.waitUntil(lambda: window.thread is None, timeout=15000)


def test_cancel_token_checkpoint_is_explicit_not_thread_termination():
    ticket = TaskToken()
    ticket.checkpoint()
    assert not ticket.cancel_requested
    ticket.request_cancel()
    assert ticket.cancel_requested
    with pytest.raises(TaskCancelled, match="Stop requested"):
        ticket.checkpoint()


def test_stop_waits_for_safe_checkpoint_without_result_or_error_dialog(
        workspace, qtbot, monkeypatch):
    w = workspace
    started = Event()
    release = Event()
    results = []
    alerts = []
    monkeypatch.setattr(
        "lsbpipe.gui.task_execution.QMessageBox.warning",
        lambda *args: alerts.append(args))
    def work(progress):
        started.set()
        release.wait(timeout=3.)
        progress((1, 2, "after-stop"))
        return "must-not-deliver"
    try:
        w._start(work, lambda result: results.append(result))
        qtbot.waitUntil(started.is_set, timeout=8000)
        assert w.stop_task_btn.isEnabled()
        w.stop_task_btn.click()
        assert not w.stop_task_btn.isEnabled()
        assert w.thread is not None
        assert "safe checkpoint" in w.message.text()
    finally:
        release.set()
    qtbot.waitUntil(lambda: w.thread is None, timeout=15000)
    assert results == []
    assert alerts == []
    assert w._active_task_token is None
    assert "not automatically" in w.science_qa.result.toPlainText()
    assert w.search_btn.isEnabled()
    assert not w.stop_task_btn.isEnabled()


def test_stop_noncooperative_work_drops_late_success_and_recovers(
        workspace, qtbot):
    w = workspace
    started = Event()
    release = Event()
    side_effect = []
    completions = []
    def work(progress):
        started.set()
        release.wait(timeout=3.)
        side_effect.append("work-finished")
        return "success-after-stop"
    try:
        w._start(work, completions.append)
        qtbot.waitUntil(started.is_set, timeout=8000)
        w._request_stop()
        assert w.thread is not None
    finally:
        release.set()
    qtbot.waitUntil(lambda: w.thread is None, timeout=15000)
    assert side_effect == ["work-finished"]
    assert completions == []
    assert "not automatically" in w.science_qa.result.toPlainText()
    # A new independent job can start, with its own token and GUI callback.
    callbacks = []
    w._start(lambda progress: "fresh", callbacks.append)
    qtbot.waitUntil(
        lambda: w.thread is None and callbacks == ["fresh"], timeout=12000)
    assert w._active_task_token is None


def test_stop_signal_delivered_on_gui_thread(workspace, qtbot):
    w = workspace
    started, release = Event(), Event()
    threads = []
    original = w.science_qa.clear_result

    def observe():
        threads.append(QThread.currentThread() == QApplication.instance().thread())
        return original()

    w.science_qa.clear_result = observe

    def work(progress):
        started.set()
        release.wait(timeout=3.)
        progress("checkpoint")
    try:
        w._start(work, lambda x: None)
        qtbot.waitUntil(started.is_set, timeout=8000)
        w._request_stop()
    finally:
        release.set()
    qtbot.waitUntil(lambda: w.thread is None and bool(threads), timeout=15000)
    assert all(threads)


def test_close_while_running_is_blocked_without_destroying_qthread(
        workspace, qtbot, monkeypatch):
    w = workspace
    started, release = Event(), Event()
    notices = []
    monkeypatch.setattr(
        "lsbpipe.gui.survey_explorer.QMessageBox.information",
        lambda *args: notices.append(args))
    def work(progress):
        started.set()
        release.wait(timeout=3.)
        return "complete"
    try:
        w._start(work, lambda result: None)
        qtbot.waitUntil(started.is_set, timeout=8000)
        w.close()
        assert w.thread is not None
        assert w.isVisible()
        assert notices
    finally:
        release.set()
    qtbot.waitUntil(lambda: w.thread is None, timeout=15000)



def test_project_switch_is_refused_until_worker_is_safely_finished(
        qtbot, monkeypatch, tmp_path):
    from lsbpipe.gui.project_hub import ProjectHub
    from lsbpipe.surveys.object_projects import create_object_project
    project = create_object_project(
        "Active-worker project", tmp_path / "active-worker")
    with offline_legacy_viewer():
        hub = ProjectHub()
        hub._enter_project(project)
    qtbot.addWidget(hub)
    qtbot.addWidget(hub.workspace)
    window = hub.workspace
    started, release = Event(), Event()
    notices = []
    monkeypatch.setattr(
        "lsbpipe.gui.survey_explorer.QMessageBox.information",
        lambda *args: notices.append(args))
    monkeypatch.setattr(
        "lsbpipe.gui.project_hub.QMessageBox.information",
        lambda *args: notices.append(args))
    def work(progress):
        started.set()
        release.wait(timeout=3.)
        progress((1, 1, "finished"))
        return "done"
    try:
        window._start(work, lambda _: None)
        qtbot.waitUntil(started.is_set, timeout=8000)
        window._go_back_to_projects()
        assert hub.workspace is window
        assert notices
        assert not hub.isVisible()
        window._request_stop()
    finally:
        release.set()
    qtbot.waitUntil(lambda: window.thread is None, timeout=15000)
    window._go_back_to_projects()
    assert hub.workspace is None
    assert hub.isVisible()


@pytest.mark.parametrize("thread_first", [True, False])
def test_terminal_and_thread_end_both_acknowledged_before_idle(
        qtbot, thread_first):
    from PySide6.QtCore import QObject
    from lsbpipe.gui.task_execution import TaskDispatch
    class Owner(QObject):
        def __init__(self):
            super().__init__()
            self._closing_task_ui = False
            self._active_task_token = None
            self.idle_calls = 0

        def _idle(self):
            self.idle_calls += 1

    owner = Owner()
    token = TaskToken()
    owner._active_task_token = token
    results = []
    adapter = TaskDispatch(owner, results.append, None, token)
    if thread_first:
        adapter.on_thread_finished()
        assert owner.idle_calls == 0
        adapter.deliver_completed("delivered")
    else:
        adapter.deliver_completed("delivered")
        assert owner.idle_calls == 0
        adapter.on_thread_finished()
    assert results == ["delivered"]
    assert adapter.terminal_outcome == "COMPLETED"
    assert owner.idle_calls == 1
    adapter.on_thread_finished()
    assert owner.idle_calls == 1


def test_object_workspace_worker_journal_reopens_with_crash_warning(
        qtbot, tmp_path):
    from lsbpipe.gui.project_hub import ProjectHub
    from lsbpipe.surveys.object_projects import create_object_project, read_object_project
    from lsbpipe.surveys.operation_journal import (
        interrupted_operations, start_operation)
    p = create_object_project("NGC 2146sat", tmp_path / "NGC2146sat")
    with offline_legacy_viewer():
        hub = ProjectHub()
        hub._enter_project(p)
    qtbot.addWidget(hub)
    window = hub.workspace
    done = []
    window._start(lambda progress: "ok", done.append,
                  operation_name="dataset_audit")
    qtbot.waitUntil(
        lambda: window.thread is None and done == ["ok"], timeout=15000)
    assert interrupted_operations(p.root) == []
    journals = list((p.root / "operations").glob("*.json"))
    assert len(journals) == 1
    import json
    assert json.loads(journals[0].read_text())["state"] == "COMPLETED"
    # Emulate a previous process dying with an active job; no manual
    # flag flips or changes to original datasets are required.
    stale = start_operation(p.root, p.project_id, "download")
    window._go_back_to_projects()
    assert hub.workspace is None
    with offline_legacy_viewer():
        hub._enter_project(read_object_project(p.root))
    reopened = hub.workspace
    assert reopened is not None
    assert "uncertain prior operation" in reopened.project_status.text()
    assert len(interrupted_operations(p.root)) == 1
    assert stale.is_file()
    reopened._go_back_to_projects()
