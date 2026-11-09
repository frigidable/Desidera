"""Read-only, explicitly selected Siena Galaxy Atlas initial-geometry priors.

SGA is *not* a galaxy extent mask or a proof of source membership. The
catalogues may be derived from Tractor/Legacy mosaics and must never be
counted as statistically independent evidence from native Tractor data.
This reader does not download, copy, index or modify catalogue FITS.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import math
from pathlib import Path

import numpy as np


MAX_SGA_BYTES = 2_600_000_000
RELEASE_HDU = {"SGA2020": "ELLIPSE", "SGA2025": "SGA2025"}


@dataclass(frozen=True)
class SGAGeometryPrior:
    release: str
    sga_id: int
    name: str
    ra_deg: float
    dec_deg: float
    offset_arcsec: float
    d26_arcmin: float | None
    d26_ref: str
    ba: float | None
    pa_deg: float | None
    ellipsebit: int
    ellipsemode: int
    flags: tuple[str, ...]
    source_sha256: str
    source_path: str
    status: str = "SGA_CORE_PRIOR_NOT_SCIENCE_ADOPTED"

    def core_shape(self):
        """D26 is a *central reference scale*, not an outer LSB boundary."""
        from lsbpipe.photometry.sky_masks import SkyMaskShape
        if (self.d26_arcmin is None or self.ba is None or
                self.pa_deg is None):
            return None
        # Diameter in arcmin -> semi-major axis in arcsec.
        return SkyMaskShape(
            "ellipse", center=(self.ra_deg, self.dec_deg),
            a_arcsec=self.d26_arcmin*30.,
            q=self.ba, pa_deg=self.pa_deg)


@dataclass(frozen=True)
class SGAPriorLookup:
    release: str
    source_path: str
    source_sha256: str
    matches: tuple[SGAGeometryPrior, ...]
    status: str = "SGA_LOCAL_MATCHES_REQUIRE_RESEARCHER_SELECTION"


def _scalar(row, key, names, default=None):
    actual = names.get(key)
    if actual is None:
        return default
    val = row[actual]
    if np.ma.is_masked(val):
        return default
    if isinstance(val, (bytes, np.bytes_)):
        return bytes(val).decode("utf-8", errors="replace").strip()
    if isinstance(val, str):
        return val.strip()
    if isinstance(val, (np.integer, int)):
        return int(val)
    try:
        n = float(val)
    except (TypeError, ValueError, OverflowError):
        return default
    return n if math.isfinite(n) else default


def _physical(value, *, upper=None):
    if value is None:
        return None
    try:
        n = float(value)
    except (ValueError, TypeError, OverflowError):
        return None
    if not math.isfinite(n) or n <= 0 or (upper is not None and n > upper):
        return None
    return n


def _sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8*1024*1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _distance_arcsec(ra, dec, target_ra, target_dec):
    dra = np.deg2rad(((ra-target_ra+180.) % 360.)-180.)
    lat = np.deg2rad(dec)
    target = math.radians(target_dec)
    chord = (np.sin((lat-target)*.5)**2
             + np.cos(lat)*math.cos(target)*np.sin(dra*.5)**2)
    return np.rad2deg(2.*np.arcsin(np.sqrt(np.clip(chord, 0, 1))))*3600.


def read_local_sga_priors(path, target_icrs, *, radius_arcmin=3.):
    """Inspect one user-chosen local SGA2020/SGA2025 FITS, in a worker.

    Rows are cone-matched by position, *not* automatically adopted. The
    researcher must disambiguate all candidates before the central SGA
    ellipse influences even the provisional Tractor-screening masks.
    """
    from astropy.io import fits
    file = Path(path).expanduser()
    if not file.is_file() or file.is_symlink():
        raise FileNotFoundError("Choose an existing non-symlink SGA FITS table")
    before = file.stat()
    if not 100 < before.st_size <= MAX_SGA_BYTES:
        raise ValueError("SGA file outside bounded catalogue size")
    try:
        target_ra, target_dec = (float(x) for x in target_icrs)
    except (ValueError, TypeError) as exc:
        raise ValueError("Invalid target sky coordinates") from exc
    if (not 0 <= target_ra < 360 or not -90 <= target_dec <= 90
            or not .1 <= float(radius_arcmin) <= 10.
            or not all(math.isfinite(x) for x in (target_ra, target_dec))):
        raise ValueError("Invalid target cone search")
    with fits.open(file, memmap=True, lazy_load_hdus=True) as hdus:
        available = {h.name.upper(): h for h in hdus
                     if isinstance(h, fits.BinTableHDU)}
        if "SGA2025" in available:
            release, hdu = "SGA2025", available["SGA2025"]
        elif "ELLIPSE" in available:
            release, hdu = "SGA2020", available["ELLIPSE"]
        else:
            raise ValueError("Neither SGA2025 nor SGA-2020 ELLIPSE extension")
        names = {name.upper(): name for name in hdu.columns.names}
        if not {"RA","DEC","D26","PA","BA"} <= names.keys():
            raise ValueError("SGA table lacks required sky/geometry fields")
        id_col = "SGAID" if release == "SGA2025" else "SGA_ID"
        if id_col not in names:
            raise ValueError("SGA table has no expected release-specific ID")
        rows = hdu.data
        ra = np.asarray(rows[names["RA"]], dtype=float)
        dec = np.asarray(rows[names["DEC"]], dtype=float)
        good = (np.isfinite(ra) & np.isfinite(dec)
                & (ra >= 0) & (ra < 360)
                & (dec >= -90) & (dec <= 90))
        separation = _distance_arcsec(ra, dec, target_ra, target_dec)
        indices = np.flatnonzero(good & (separation <= radius_arcmin*60.))
        if len(indices) > 32:
            raise ValueError(
                "More than 32 plausible SGA matches; narrow the cone "
                "before selecting a specific galaxy")
        indices = sorted(indices, key=lambda i: float(separation[i]))
        raw = []
        for index in indices:
            row = rows[int(index)]
            d26 = _physical(_scalar(row, "D26", names), upper=60.)
            ba = _physical(_scalar(row, "BA", names), upper=1.)
            pa = _scalar(row, "PA", names)
            if pa is None or not 0 <= pa <= 360:
                pa = None
            ident = _scalar(row, id_col, names)
            if ident is None or ident < 0:
                continue
            flags = []
            if d26 is None:
                flags.append("NO_VALID_D26")
            if pa is None or ba is None:
                flags.append("NO_VALID_SKY_ELLIPSE")
            ref = str(_scalar(row, "D26_REF", names, "") or "")
            if release == "SGA2020" and ref != "SB26":
                flags.append("D26_NOT_DIRECT_SB26")
            bit = int(_scalar(row, "ELLIPSEBIT", names, 0) or 0)
            mode = int(_scalar(row, "ELLIPSEMODE", names, 0) or 0)
            if release == "SGA2020":
                # SGA-2020 ELLIPSEBIT is a bit-index mask; not the
                # similarly named SGA2025 bit-value/status vocabulary.
                if bit & ((1 << 1) | (1 << 2)):
                    flags.append("SGA2020_TOO_SMALL_FOR_ELLIPSE_FIT")
                if bit & (1 << 3):
                    flags.append("SGA2020_ELLIPSE_FIT_FAILED")
                if bit & (1 << 4):
                    flags.append("SGA2020_ELLIPSE_NOT_FITTED")
                if bit & (1 << 5):
                    flags.append("SGA2020_ELLIPSE_REJECTED")
            if release == "SGA2025":
                flags.append("SGA2025_DRAFT_NOT_FINAL_RELEASE")
                if bit & 16384:
                    flags.append("SGA_FIT_GEOMETRY_FAILED")
                if bit & 16:
                    flags.append("SGA_LARGE_CENTRE_SHIFT")
                if bit & 8:
                    flags.append("SGA_BLENDED")
                if mode & 1 or bit & 4:
                    flags.append("SGA_FIXED_GEOMETRY")
                if ref in ("mom","fix"):
                    flags.append("D26_FALLBACK_OR_FIXED")
            raw.append(dict(
                release=release, sga_id=int(ident),
                name=str(_scalar(row, "GALAXY", names, "") or
                         _scalar(row, "SGANAME", names, "") or ident),
                ra_deg=float(ra[index]), dec_deg=float(dec[index]),
                offset_arcsec=float(separation[index]),
                d26_arcmin=d26, d26_ref=ref, ba=ba,
                pa_deg=(float(pa) % 180. if pa is not None else None),
                ellipsebit=bit, ellipsemode=mode, flags=tuple(flags)))
    # Confirm file identity after FITS is closed. No persistent cache is
    # silently trusted if the user later modifies/replaces the file.
    source_hash = _sha256(file)
    after = file.stat()
    if (before.st_size, before.st_mtime_ns) != (
            after.st_size, after.st_mtime_ns):
        raise RuntimeError("SGA FITS changed while being read")
    canonical = str(file.resolve())
    matches = tuple(SGAGeometryPrior(
        **row, source_path=canonical, source_sha256=source_hash)
        for row in raw)
    return SGAPriorLookup(
        release=release, source_path=canonical,
        source_sha256=source_hash, matches=matches)
