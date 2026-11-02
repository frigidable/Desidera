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

    def _current(self) -> bool:
        return (getattr(self.owner_window, "_active_task_token", None)
                is self.token and
                not getattr(self.owner_window, "_closing_task_ui", False))

    @Slot(object)
    def deliver_completed(self, value):
        if not self._current():
            return
        if self.token.cancel_requested:
            self.deliver_cancelled(
                "Stop requested. Work may have written immutable receipts; "
                "inspect saved run history before reuse.")
            return
        try:
            self.completed_callback(value)
        except Exception as exc:
            logging.exception("LSBpipe GUI completion callback failed")
            self.deliver_failed(f"{type(exc).__name__}: {exc}")

    @Slot(str)
    def deliver_failed(self, message: str):
        if not self._current():
            return
        self.owner_window.message.setText(message)
        if self.error_callback is not None:
            try:
                self.error_callback(message)
            except Exception:
                logging.exception("LSBpipe GUI error callback failed")
        QMessageBox.warning(self.owner_window, "Survey operation failed", message)

    @Slot(str)
    def deliver_cancelled(self, message: str):
        if not self._current():
            return
        self.owner_window.message.setText(
            f"Stop request acknowledged. {message} "
            "No new result has been displayed as accepted science.")
        if hasattr(self.owner_window, "science_qa"):
            self.owner_window.science_qa.clear_result()
            self.owner_window.science_qa.result.setPlainText(
                "Task stopped or its late result discarded. Any existing "
                "immutable scientific receipts remain on disk, and are "
                "not automatically scientifically adopted.")

    @Slot(object)
    def deliver_progress(self, value):
        if self._current() and not self.token.cancel_requested:
            self.owner_window._progress(value)
