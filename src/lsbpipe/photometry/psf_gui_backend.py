"""GUI orchestration for the Leo P-first PSF-0 experiment.

No Qt imports; work is executed by existing SurveyExplorer TaskWorker.
Uses the saved object-project manifest and the original legacy FITS tree.
Never creates a second project, moves a source selection, downloads, or
changes main project state. The selected analysis band is independent of
the saved Galaxy Workspace *default band*.
"""
from __future__ import annotations

from hashlib import sha256
from dataclasses import replace
import json
import math
import os
import tempfile
from pathlib import Path
from typing import Callable

from lsbpipe.photometry.native_fits import load_native_legacy_quartet
from lsbpipe.photometry.psf_empirical import (
    PSF0Config, PSF0Result, PSF0QualityError, build_empirical_coadd_psf,
)
from lsbpipe.photometry.psf_export import export_psf0_diagnostic
from lsbpipe.surveys.dataset_integrity import _expected_paths
from lsbpipe.surveys.object_projects import read_object_project

LEO_P_CENTER = (155.438012, 18.088025)
BANDS = ("g", "r", "i", "z")
KINDS = ("image", "invvar", "nexp", "maskbits")

# Reproducible diagnostic alternatives. Wing-smoothing and iteration
# modes each change one PSFBuilder parameter; the larger-sky mode must
# jointly increase cutout, sky-annulus and isolation footprints. Its
# stellar sample can therefore differ: do not claim a controlled
# single-variable causal inference or select a "best PSF" from these QA
# summaries. Every experiment remains non-adopted.
PSF0_DIAGNOSTIC_EXPERIMENTS = {
    "baseline": {
        "label": "Baseline · 31 px, sky 10–15, 8 iterations",
        "kwargs": {},
    },
    "outer_sky": {
        "label": "Extended cutouts · 51 px, sky 19–24",
        "kwargs": {
            "cutout_size_pix": 51,
            "sky_annulus_inner_pix": 19.,
            "sky_annulus_outer_pix": 24.,
            "min_neighbor_separation_pix": 55.,
        },
    },
    "wing_smoothing_off": {
        "label": "Disable final wing smoothing · 31 px",
        "kwargs": {"wing_smoothing": False},
    },
    "longer_convergence": {
        "label": "Extended convergence · 16 iterations",
        "kwargs": {"max_iterations": 16},
    },
}


def psf0_experiment_recipe(guessed_fwhm_pix: float, experiment: str) -> PSF0Config:
    if experiment not in PSF0_DIAGNOSTIC_EXPERIMENTS:
        raise ValueError("Unknown controlled PSF-0 diagnostic experiment")
    return PSF0Config(
        guessed_fwhm_pix=float(guessed_fwhm_pix),
        **PSF0_DIAGNOSTIC_EXPERIMENTS[experiment]["kwargs"])




def _leo_p_project(project_root):
    """Use existing UUID and saved source identity, never a guessed folder."""
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    project = read_object_project(project_root)
    if project.display_name.strip().casefold().replace(" ", "") != "leop":
        raise ValueError("PSF-0 real-object pilot currently supports Leo P only")
    if project.field is None or project.data_root is None or project.selected_dataset is None:
        raise ValueError("Leo P needs its saved coordinates, selected dataset and data root")
    center = SkyCoord(project.field.ra_deg*u.deg, project.field.dec_deg*u.deg)
    expected = SkyCoord(LEO_P_CENTER[0]*u.deg, LEO_P_CENTER[1]*u.deg)
    if center.separation(expected).arcsec > 45:
        raise ValueError("Leo P project center differs by over 45 arcsec; review galaxy identity")
    return project


def inspect_psf0_bands(project_root: str | Path) -> dict:
    """Metadata/path-only availability for the current Galaxy Workspace."""
    project = _leo_p_project(project_root)
    release, region, brick, selected = project.selected_dataset
    entries = []
    for band in BANDS:
        if region == "north" and band == "i":
            continue
        locations = _expected_paths(project.data_root, (release, region, brick, band))
        present = {kind: path.is_file() and not path.is_symlink()
                   for kind, path in locations.items()}
        entries.append({
            "band": band,
            "selected": band == selected,
            "present_count": sum(present.values()),
            "paths_present_unverified": all(present.values()),
            "files": present,
        })
    return {
        "project_id": project.project_id,
        "project_root": str(project.root),
        "dataset": (release, region, brick, selected),
        "bands": entries,
        "target_protection_default_arcsec": max(project.target_radius_arcsec, 90.),
        "status": "PATHS_ONLY_NOT_SCIENCE_VALIDATED",
    }


def _progress(progress: Callable | None, text: str):
    if progress is not None:
        progress((0, 1, text))
        checkpoint = getattr(progress, "checkpoint", None)
        if checkpoint is not None:
            checkpoint()


def run_gui_psf0(
    project_root: str | Path, *,
    band: str, config: PSF0Config,
    protection_arcsec: float, progress: Callable | None = None,
) -> dict:
    """Read, validate and build in background, return result without saving."""
    snapshot = inspect_psf0_bands(project_root)
    chosen = next((x for x in snapshot["bands"] if x["band"] == band), None)
    if chosen is None or not chosen["paths_present_unverified"]:
        raise ValueError(f"Selected {band} native FITS quartet is incomplete or unsupported")
    if not isinstance(config, PSF0Config):
        raise TypeError("PSF-0 configuration is required")
    if not math.isfinite(protection_arcsec) or not 0 < protection_arcsec <= 1800:
        raise ValueError("Invalid galaxy-protection radius")
    release, region, brick, _ = snapshot["dataset"]
    project = _leo_p_project(project_root)
    _progress(progress, f"Checking {release.upper()} {band}: original FITS and source SHA-256")
    frame = load_native_legacy_quartet(
        project.data_root, release=release, region=region, brick=brick,
        band=band, hash_inputs=True, require_explicit_bunit=True)
    _progress(progress, f"Building provisional {band} PSF · star selection, sky, holdout")
    try:
        result = build_empirical_coadd_psf(
            frame, target_icrs_deg=LEO_P_CENTER,
            target_protection_arcsec=float(protection_arcsec), config=config)
    except PSF0QualityError as exc:
        # This is a *normal scientific rejection* with useful candidate-star
        # evidence, not an exception that should lose context in a GUI popup.
        return {
            "status": "PSF0_QUALITY_REJECTED_NOT_ADOPTED",
            "reason": str(exc),
            "diagnostics": exc.diagnostics,
            "star_review": exc.review_snapshot,
            "source_hashes": dict(frame.source_evidence.get("local_hashes", {})),
            "project_id": snapshot["project_id"],
            "project_root": snapshot["project_root"],
            "dataset": snapshot["dataset"],
        }
    _progress(progress, f"Completed provisional {band} PSF · not science-adopted")
    return {
        "status": "EXPERIMENTAL_PSF_NOT_ADOPTED",
        "result": result,
        "star_review": result.review_snapshot,
        "source_hashes": dict(result.source_hashes),
        "project_id": snapshot["project_id"],
        "project_root": snapshot["project_root"],
        "dataset": snapshot["dataset"],
    }


def _hash_original_file(path: Path) -> str:
    sha = sha256()
    with path.open("rb") as reader:
        for chunk in iter(lambda: reader.read(4 * 1024 * 1024), b""):
            sha.update(chunk)
    return sha.hexdigest()


def save_gui_psf0(preview: dict, progress: Callable | None = None) -> str:
    """Explicit user-initiated diagnostic export with source-change detection."""
    if not isinstance(preview, dict) or not isinstance(preview.get("result"), PSF0Result):
        raise TypeError("Run a successful PSF-0 preview before saving")
    result = preview["result"]
    if result.lineage_id is None or not result.source_hashes:
        raise ValueError("Verified source hashes required before diagnostic save")
    project = _leo_p_project(preview["project_root"])
    if (project.project_id != preview["project_id"]
            or tuple(project.selected_dataset) != tuple(preview["dataset"])):
        raise ValueError("Galaxy project or selected dataset changed; rerun the diagnostic")
    release, region, brick, _ = project.selected_dataset
    if result.dataset_id != f"{release}/{region}/{brick}/{result.dataset_id.split('/')[-1]}":
        raise ValueError("PSF result does not match the Galaxy Workspace source")
    band = result.dataset_id.split("/")[-1]
    expected_paths = _expected_paths(project.data_root, (release, region, brick, band))
    _progress(progress, "Rechecking original FITS bytes against PSF source fingerprints")
    for kind in KINDS:
        p = expected_paths[kind]
        if not p.is_file() or p.is_symlink():
            raise ValueError(f"Original {kind.upper()} FITS missing/unsafe; cannot save")
        if _hash_original_file(p) != result.source_hashes[kind]:
            raise ValueError(f"Original {kind.upper()} FITS changed since PSF build; rerun")
    _progress(progress, "Publishing immutable experimental PSF FITS + QA report")
    output = export_psf0_diagnostic(
        result, project.root / "analyses" / "psf-pilots", explicitly_requested=True)
    return str(output)
