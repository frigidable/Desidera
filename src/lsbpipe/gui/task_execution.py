"""Qt-thread-safe job execution primitives with cooperative stop requests.

Workers NEVER access Qt widgets. All completion, failure and progress callbacks
are delivered through GUI-owned QObject slots on the GUI thread. Cancellation
does not terminate threads or pretend that already committed reports disappeared.
"""
from __future__ import annotations

import logging
from threading import Event
from typing import Callable

from PySide6.QtCore import QObject, Signal, Slot
from PySide6.QtWidgets import QMessageBox


class TaskCancelled(Exception):
    """Raised at a cooperative checkpoint on the worker thread."""


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
                return
            if self.token.cancel_requested:
                self._show_cancelled(
                    "Stop requested. Work may have written immutable receipts; "
                    "inspect saved run history before reuse.")
                return
            try:
                self.completed_callback(value)
            except Exception as exc:
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
                self._show_cancelled(message)
        finally:
            self._acknowledge_terminal()

    @Slot(object)
    def deliver_progress(self, value):
        if self._current() and not self.token.cancel_requested:
            self.owner_window._progress(value)
