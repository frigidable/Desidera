"""Bounded, *display-only* previews of native optical IMAGE arrays.

No science data, WCS, noise maps or pixel readouts are downsampled.
An area-mean display overview preserves mean surface brightness of finite
samples within each nonoverlapping source-pixel block; it is NOT a
PSF-matched image, replacement FITS, or photometric measurement.
"""
from __future__ import annotations

import math
import numpy as np

MAX_PREVIEW_SIDE = 1024
MAX_NATIVE_BYTES = 512 * 1024 * 1024


def native_dataset_bytes(dataset: dict) -> int:
    """Nominal uncompressed native pixel memory; conservative for memmaps."""
    return sum(int(np.asarray(arr).nbytes)
               for arr in dataset["data"].values())


def display_preview(data: np.ndarray, *, max_side: int = MAX_PREVIEW_SIDE
                    ) -> tuple[np.ndarray, int]:
    """Finite-pixel block average, bounded raster only; source untouched.

    Returns a finite/NaN float32 overview and integer bin factor.
    Each output sample averages a nonoverlapping square of source pixels,
    with edge blocks averaging only existing pixels. The caller must
    scale the display item back to the native scene size so WCS hover and
    native pixel readout do NOT change coordinates.
    """
    arr = np.asarray(data)
    if arr.ndim != 2 or 0 in arr.shape:
        raise ValueError("Overview requires a nonempty 2D native IMAGE")
    if not isinstance(max_side, int) or max_side < 32:
        raise ValueError("Preview side limit must be an integer >= 32")
    ny, nx = arr.shape
    factor = max(1, math.ceil(max(ny, nx) / max_side))
    if factor == 1:
        # No conversion of the native input is required; the normal
        # rendering transfer function still creates a separate RGB raster.
        return arr, factor
    oy, ox = math.ceil(ny / factor), math.ceil(nx / factor)
    expanded = np.full((oy * factor, ox * factor),
                       np.nan, dtype=np.float32)
    expanded[:ny, :nx] = arr
    blocks = expanded.reshape(oy, factor, ox, factor)
    finite = np.isfinite(blocks)
    totals = np.where(finite, blocks, 0.).sum(
        axis=(1, 3), dtype=np.float64)
    counts = finite.sum(axis=(1, 3), dtype=np.int64)
    out = np.full((oy, ox), np.nan, dtype=np.float32)
    np.divide(totals, counts, out=out, where=counts > 0)
    return out, factor
