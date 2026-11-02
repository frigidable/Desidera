"""Strict read-only adapter: original Legacy Surveys quartet to P0.5 frame.

Optical MASKBITS is HDU1, never the WISE extensions. Existing P0.2
validation provides HDU/WCS/checksum gates, not science acceptance.
"""
from __future__ import annotations

from contextlib import ExitStack
from pathlib import Path
import json

from lsbpipe.surveys.dataset_integrity import (
    KINDS, _expected_paths, _sha256_file, _validate_open_hdus)
from lsbpipe.surveys.legacy import BRICK_RE, regions_for_release
from lsbpipe.surveys.task_control import checkpoint
from .contracts import MaskPolicy, NativeBandFrame

_ALLOWED_IMAGE_UNITS = {
    "nanomaggy", "nanomaggies", "nmgy",
    "nanomaggy/pixel", "nanomaggies/pixel", "nmgy/pixel",
    "nanomaggy/pix", "nanomaggies/pix", "nmgy/pix"}
_ALLOWED_IVAR_UNITS = {
    "1/nanomaggies^2", "1/nanomaggy^2", "1/nmgy^2",
    "nanomaggies^-2", "nanomaggy^-2", "nmgy^-2",
    "1/(nanomaggies)^2", "1/(nanomaggy)^2",
}


def _norm(value):
    return str(value).lower().replace(" ", "").replace("**", "^")


def load_native_legacy_quartet(
    data_root: str | Path, *, release: str, region: str,
    brick: str, band: str, mask_policy: MaskPolicy | None = None,
    require_explicit_bunit: bool = True, hash_inputs: bool = True,
    cancel_check=None,
) -> NativeBandFrame:
    """Load one unmodified DR9/10/11 optical quartet.

    For missing IMAGE BUNIT an explicit caller opt-in
    require_explicit_bunit=False is required, documented in evidence;
    no hidden arbitrary conversion or zeropoint adjustment is performed.
    No FITS file is written, padded, interpolated or reprojected.

    File integrity SHA-256 is local only, not an archive authenticity claim.
    """
    from astropy.io import fits
    if region not in regions_for_release(release):
        raise ValueError("Invalid survey release/region")
    if band not in ("g", "r", "i", "z") or (region == "north" and band == "i"):
        raise ValueError("Unsupported optical band in this region")
    if not isinstance(brick, str) or not BRICK_RE.fullmatch(brick):
        raise ValueError("Invalid native brick")
    # The root must contain <release>/<region>/coadd, not the output folder.
    root = Path(data_root).expanduser().resolve()
    paths = _expected_paths(root, (release, region, brick, band))
    for name, path in paths.items():
        if not path.is_file() or path.is_symlink():
            raise FileNotFoundError(f"Missing safe original {name}: {path}")
    checkpoint(cancel_check)
    stats = {name: path.stat() for name, path in paths.items()}
    with ExitStack() as stack:
        hdus = {
            name: stack.enter_context(fits.open(path, memmap=True))
            for name, path in paths.items()}
        arrays, wcs, indices, units = _validate_open_hdus(hdus, band)
        image_unit = units["image"]
        if image_unit is None or not str(image_unit).strip():
            if require_explicit_bunit:
                raise ValueError(
                    "IMAGE BUNIT missing: explicit documented survey-unit "
                    "assumption required, not silent calibration")
            unit_evidence = "LEGACY_SURVEY_DOCUMENTATION_ASSUMPTION"
        elif _norm(image_unit) not in _ALLOWED_IMAGE_UNITS:
            raise ValueError(f"Unexpected IMAGE BUNIT: {image_unit!r}")
        else:
            unit_evidence = "IMAGE_FITS_BUNIT_CONFIRMED"
        iv_unit = units["invvar"]
        if iv_unit is not None and str(iv_unit).strip():
            if _norm(iv_unit) not in _ALLOWED_IVAR_UNITS:
                raise ValueError(f"Unexpected inverse-variance BUNIT: {iv_unit!r}")
            ivar_evidence = "INVVAR_FITS_BUNIT_CONFIRMED"
        else:
            ivar_evidence = "LEGACY_SURVEY_DOCUMENTATION_ASSUMPTION"
        checkpoint(cancel_check)
        frame = NativeBandFrame(
            image_nMgy=arrays["image"],
            invvar_nMgy2=arrays["invvar"],
            nexp=arrays["nexp"],
            maskbits=arrays["maskbits"],
            native_wcs=wcs,
            release=release, region=region, brick=brick, band=band,
            policy=mask_policy if mask_policy is not None else MaskPolicy(),
            source_evidence={
                "data_root": str(root),
                "files": {k: str(v) for k, v in paths.items()},
                "hdu_indices": indices,
                "observed_header_units": units,
                "image_unit_evidence": unit_evidence,
                "invvar_unit_evidence": ivar_evidence,
                "data_status": "NATIVE_FITS_WCS_VALIDATED_NOT_SCIENCE_ADOPTED",
                "local_hashes": {},
            })
        # All pixel arrays are now immutable copies detached from the HDU.
    checkpoint(cancel_check)
    hashes = {}
    for name, path in paths.items():
        checkpoint(cancel_check)
        before = stats[name]
        if hash_inputs:
            info = _sha256_file(path, cancel_check=cancel_check)
            hashes[name] = info["sha256_local"]
            size, modified = info["size_bytes"], info["mtime_ns"]
        else:
            final = path.stat()
            size, modified = final.st_size, final.st_mtime_ns
        if (before.st_size, before.st_mtime_ns) != (size, modified):
            raise RuntimeError(
                "Original FITS changed during photometric input preparation: "
                + str(path))
    # Rebuild only metadata; pixel data/WCS remain their detached immutable
    # snapshots. Avoid hashing data while holding FITS file handles open.
    evidence = {
        "data_root": str(root),
        "files": {k: str(v) for k, v in paths.items()},
        "hdu_indices": indices,
        "observed_header_units": units,
        "image_unit_evidence": unit_evidence,
        "invvar_unit_evidence": ivar_evidence,
        "data_status": "NATIVE_FITS_WCS_VALIDATED_NOT_SCIENCE_ADOPTED",
        "local_hashes": hashes,
        "hash_provenance": ("LOCAL_SHA256_NOT_REMOTE_AUTHENTICITY"
                            if hash_inputs else "INPUT_HASHES_NOT_COMPUTED"),
    }
    checkpoint(cancel_check)
    return NativeBandFrame(
        image_nMgy=frame.image_nMgy,
        invvar_nMgy2=frame.invvar_nMgy2,
        nexp=frame.nexp,
        maskbits=frame.maskbits,
        native_wcs=frame.wcs,
        release=release, region=region, brick=brick, band=band,
        policy=frame.policy,
        source_evidence=evidence)
