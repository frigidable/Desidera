"""Pure-Python cooperative cancellation for science and survey I/O.

The scientific core never imports Qt. Stop requests are checked only at
documented safe points and never forcibly interrupt atomic publication.
"""
from __future__ import annotations
from typing import Callable


class OperationCancelled(Exception):
    """The caller requested cancellation at a safe processing checkpoint."""


def checkpoint(cancel_check: Callable[[], None] | None = None) -> None:
    if cancel_check is not None:
        cancel_check()
