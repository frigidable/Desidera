"""Cross-platform, crash-released advisory lock for local metadata transactions.

A stable sidecar file is intentionally NOT deleted after unlock. Removing a
lock file while another process is waiting on its inode can split the lock
domain and cause lost updates. These locks protect Desidera processes that
follow this protocol; they do not constrain arbitrary third-party file writes.

No Qt, network, FITS or additional directories beyond the existing metadata
parent. A lock guards the entire read-modify-atomic-replace interval.
"""
from __future__ import annotations

from contextlib import contextmanager
import errno
import math
import os
from pathlib import Path
import stat
import time


@contextmanager
def exclusive_metadata_lock(path: str | Path, *, timeout_seconds: float = 10.0):
    """Hold one same-filesystem lock until the guarded transaction completes.

    Unix: flock(LOCK_EX|LOCK_NB); Windows: msvcrt.locking(LK_NBLCK, 1).
    The OS releases the byte-range/flock lock when the process dies.
    The path persists to preserve lock identity between processes.
    """
    if (isinstance(timeout_seconds, bool) or
            not isinstance(timeout_seconds, (int, float)) or
            not math.isfinite(timeout_seconds) or timeout_seconds < 0):
        raise ValueError("Invalid metadata lock timeout")
    file = Path(path).expanduser().absolute()
    if (file.parent.is_symlink() or file.is_symlink() or
            (file.exists() and not file.is_file())):
        raise ValueError("Unsafe metadata lock path")
    file.parent.mkdir(parents=True, exist_ok=True)
    flags = os.O_RDWR | os.O_CREAT
    if hasattr(os, "O_BINARY"):
        flags |= os.O_BINARY
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    fd = os.open(file, flags, 0o600)
    acquired = False
    try:
        # Detect a swapped pathname, a malicious symlink or a nonregular
        # device. Linux O_NOFOLLOW also rejects the symlink atomically.
        opened = os.fstat(fd)
        if file.is_symlink() or not stat.S_ISREG(opened.st_mode):
            raise ValueError("Unsafe metadata lock file")
        visible = file.stat()
        if (opened.st_dev, opened.st_ino) != (visible.st_dev, visible.st_ino):
            raise ValueError("Metadata lock path changed during open")
        # Windows CRT locks a byte range starting at the file offset.
        # The stable sentinel is not metadata and need not be interpreted.
        if os.name == "nt":
            if opened.st_size == 0:
                os.lseek(fd, 0, os.SEEK_SET)
                os.write(fd, b"0")
            import msvcrt
        else:
            import fcntl
        deadline = time.monotonic() + float(timeout_seconds)
        while True:
            try:
                if os.name == "nt":
                    os.lseek(fd, 0, os.SEEK_SET)
                    msvcrt.locking(fd, msvcrt.LK_NBLCK, 1)
                else:
                    fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
                acquired = True
                break
            except OSError as exc:
                # Locked by another process, not an I/O error or a bad fd.
                if exc.errno not in (errno.EACCES, errno.EAGAIN, errno.EDEADLK):
                    raise
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    raise TimeoutError(
                        f"Metadata lock busy: {file}") from exc
                time.sleep(min(0.05, remaining))
        # The pathname must still refer to our opened inode. No rename,
        # overwrite or unlink of the permanent lock file is permitted.
        visible = file.stat()
        if (opened.st_dev, opened.st_ino) != (visible.st_dev, visible.st_ino):
            raise ValueError("Metadata lock file was replaced during wait")
        yield
    finally:
        try:
            if acquired:
                if os.name == "nt":
                    os.lseek(fd, 0, os.SEEK_SET)
                    msvcrt.locking(fd, msvcrt.LK_UNLCK, 1)
                else:
                    fcntl.flock(fd, fcntl.LOCK_UN)
        finally:
            os.close(fd)
