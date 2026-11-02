"""Qt-thread-safe job execution primitives with cooperative stop requests.

Workers NEVER access Qt widgets. All completion, failure and progress callbacks
are delivered through GUI-owned QObject slots on the GUI thread. Cancellation
does not terminate threads or pretend that already committed reports disappeared.
"""
from __future__ import annotations

import logging
from threading import Event
from typing import Callable

from PySide6.QtCore import QObject, QThread, Qt, Signal, Slot
from PySide6.QtWidgets import QMessageBox

# Preserve the public GUI exception while allowing pure scientific modules
# to propagate cooperative cancellation without importing any Qt bindings.
from lsbpipe.surveys.task_control import OperationCancelled as TaskCancelled


class TaskToken:
    def __init__(self):
        self._cancel = Event()

    def request_cancel(self):
        self._cancel.set()

    @property
    def cancel_requested(self) -> bool:
        return self._cancel.is_set()

    def checkpoint(self):
        if self.cancel_requested:
            raise TaskCancelled("Stop requested; worker reached a safe checkpoint")


class TaskWorker(QObject):
    completed = Signal(object)
    progress = Signal(object)
    failed = Signal(str)
    cancelled = Signal(str)
    finished = Signal()

    def __init__(self, task: Callable, token: TaskToken | None = None):
        super().__init__()
        self.task = task
        self.token = token if token is not None else TaskToken()

    @Slot()
    def run(self):
        try:
            self.token.checkpoint()

            def checked_progress(value):
                self.token.checkpoint()
                self.progress.emit(value)
                self.token.checkpoint()

            # The work function can pass this hook into long-running,
            # Qt-independent FITS/HTTP APIs without emitting GUI progress
            # for every 1 MiB read or every raster slab.
            checked_progress.checkpoint = self.token.checkpoint
            result = self.task(checked_progress)
            self.token.checkpoint()
            self.completed.emit(result)
        except TaskCancelled as exc:
            self.cancelled.emit(str(exc))
        except Exception as exc:
            if self.token.cancel_requested:
                self.cancelled.emit(
                    "Stop requested; original work ended with "
                    f"{type(exc).__name__}. Check saved receipts before reuse.")
            else:
                self.failed.emit(f"{type(exc).__name__}: {exc}")
        finally:
            self.finished.emit()


class TaskDispatch(QObject):
    """GUI-owned adapter; ignores events from superseded jobs."""

    def __init__(self, owner, completed: Callable,
                 on_error: Callable | None, token: TaskToken):
        super().__init__(owner)
        self.owner_window = owner
        self.completed_callback = completed
        self.error_callback = on_error
        self.token = token
        self._terminal_seen = False
        self._thread_done = False
        self._finalized = False
        self.terminal_outcome = "IGNORED"

    def _acknowledge_terminal(self):
        self._terminal_seen = True
        self._maybe_finalize()

    @Slot()
    def on_thread_finished(self):
        self._thread_done = True
        self._maybe_finalize()

    def _maybe_finalize(self):
        # The worker's terminal result is queued independently of
        # QThread.finished. Never reset the active token or delete this
        # receiver until BOTH events have reached this GUI thread.
        if not (self._terminal_seen and self._thread_done):
            return
        if self._finalized:
            return
        self._finalized = True
        self.owner_window._idle()
        self.deleteLater()

    def _current(self) -> bool:
        return (getattr(self.owner_window, "_active_task_token", None)
                is self.token and
                not getattr(self.owner_window, "_closing_task_ui", False))

    @Slot(object)
    def deliver_completed(self, value):
        try:
            if not self._current():
                self.terminal_outcome = "IGNORED"
                return
            if self.token.cancel_requested:
                self.terminal_outcome = "CANCELLED"
                self._show_cancelled(
                    "Stop requested. Work may have written immutable receipts; "
                    "inspect saved run history before reuse.")
                return
            try:
                self.completed_callback(value)
                self.terminal_outcome = "COMPLETED"
            except Exception as exc:
                self.terminal_outcome = "FAILED"
                logging.exception("LSBpipe GUI completion callback failed")
                self._show_failure(f"{type(exc).__name__}: {exc}")
        finally:
            self._acknowledge_terminal()

    def _show_failure(self, message: str):
        self.owner_window.message.setText(message)
        if self.error_callback is not None:
            try:
                self.error_callback(message)
            except Exception:
                logging.exception("LSBpipe GUI error callback failed")
        QMessageBox.warning(self.owner_window, "Survey operation failed", message)

    @Slot(str)
    def deliver_failed(self, message: str):
        try:
            if self._current():
                self.terminal_outcome = "FAILED"
                self._show_failure(message)
        finally:
            self._acknowledge_terminal()

    def _show_cancelled(self, message: str):
        self.owner_window.message.setText(
            f"Stop request acknowledged. {message} "
            "No new result has been displayed as accepted science.")
        if hasattr(self.owner_window, "science_qa"):
            self.owner_window.science_qa.clear_result()
            self.owner_window.science_qa.result.setPlainText(
                "Task stopped or its late result discarded. Any existing "
                "immutable scientific receipts remain on disk, and are "
                "not automatically scientifically adopted.")

    @Slot(str)
    def deliver_cancelled(self, message: str):
        try:
            if self._current():
                self.terminal_outcome = "CANCELLED"
                self._show_cancelled(message)
        finally:
            self._acknowledge_terminal()

    @Slot(object)
    def deliver_progress(self, value):
        if self._current() and not self.token.cancel_requested:
            self.owner_window._progress(value)

class WorkspaceTaskController(QObject):
    """Owns one QThread, GUI delivery, cancellation and crash journal.

    Public SurveyExplorer._start/_request_stop/_idle remain thin compatibility
    delegates. GUI controls and Qt object lifetimes are managed in one place;
    science callables are always passed immutable work inputs.
    """

    def __init__(self, owner):
        super().__init__(owner)
        self.owner = owner

    def start(self, action, finished, on_error=None, *,
              operation_name="background"):
        w = self.owner
        if w.thread is not None:
            QMessageBox.information(
                w, "Busy", "Another survey operation is running.")
            return
        if w._object_project is not None:
            from lsbpipe.surveys.operation_journal import start_operation
            try:
                w._operation_journal_path = start_operation(
                    w._object_project.root, w._object_project.project_id,
                    operation_name)
            except (OSError, ValueError) as exc:
                QMessageBox.warning(
                    w, "Cannot journal operation",
                    "No background work was started: " + str(exc))
                return
        else:
            w._operation_journal_path = None

        w._task_serial += 1
        token = TaskToken()
        w._active_task_token = token
        w.search_btn.setEnabled(False)
        w.download_btn.setEnabled(False)
        w.verify_btn.setEnabled(False)
        w.qa_btn.setEnabled(False)
        w.stop_task_btn.setEnabled(True)
        w.science_qa.set_busy(True)

        thread = QThread(w)
        worker = TaskWorker(action, token=token)
        dispatch = TaskDispatch(w, finished, on_error, token)
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.completed.connect(dispatch.deliver_completed, Qt.QueuedConnection)
        worker.failed.connect(dispatch.deliver_failed, Qt.QueuedConnection)
        worker.cancelled.connect(dispatch.deliver_cancelled, Qt.QueuedConnection)
        worker.progress.connect(dispatch.deliver_progress, Qt.QueuedConnection)
        worker.finished.connect(thread.quit)
        worker.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)
        # The delivery barrier lives on the main Qt thread and outlives
        # QThread.finished until the final result/cancellation is handled.
        thread.finished.connect(
            dispatch.on_thread_finished, Qt.QueuedConnection)
        w._task_dispatch = dispatch
        w.worker = worker
        w.thread = thread
        thread.start()

    def request_stop(self):
        w = self.owner
        token = w._active_task_token
        if w.thread is None or token is None:
            return
        token.request_cancel()
        w.thread.requestInterruption()
        if w._operation_journal_path is not None:
            from lsbpipe.surveys.operation_journal import request_stop
            try:
                request_stop(w._operation_journal_path)
            except (OSError, ValueError) as exc:
                w.project_status.setText(
                    f"Stop requested but journal write failed: {exc}")
        w.stop_task_btn.setEnabled(False)
        w.message.setText(
            "Stop requested. Waiting for a safe checkpoint; the running "
            "operation may finish writing an immutable report before stopping.")

    def idle(self):
        w = self.owner
        if w._operation_journal_path is not None:
            from lsbpipe.surveys.operation_journal import finish_operation
            try:
                finish_operation(
                    w._operation_journal_path,
                    w._task_dispatch.terminal_outcome
                    if w._task_dispatch is not None else "IGNORED")
            except (OSError, ValueError) as exc:
                w.project_status.setText(
                    f"Operation finished; journal is uncertain: {exc}")
        w._operation_journal_path = None
        w.thread = None
        w.worker = None
        w._task_dispatch = None
        w._active_task_token = None
        w.stop_task_btn.setEnabled(False)
        w.search_btn.setEnabled(True)
        w.download_btn.setEnabled(True)
        w.verify_btn.setEnabled(True)
        w.qa_btn.setEnabled(True)
        w.science_qa.set_busy(False)
