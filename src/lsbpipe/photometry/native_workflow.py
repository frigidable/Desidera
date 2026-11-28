"""Staged, opt-in native annular photometry workflow.

A candidate scene or geometry is never adopted implicitly. No network calls,
FITS writes, automatic mask restoration or calibrated uncertainties.
"""
from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Mapping
import copy
import numpy as np

from .contracts import NativeBandFrame, PSFDescription
from lsbpipe.scene import SceneState
from .stellar_contamination import fit_native_independent_sky_halo
from .native_annular_diagnostic import (
    NativeCorrection, common_native_validity, prepare_native_annuli,
    fit_native_covariance_sersic,
)
from .geometry_inference import fit_native_geometry_candidates, reviewed_geometry_candidate
from .radial_models import NestedEllipse
from .block_covariance import native_scene_block_covariance
from .annular_export import export_annular_diagnostic
from .asymptotic_diagnostic import asymptotic_flux_diagnostic
from .uncertainty_budget import STATUS


class WorkflowBlocked(ValueError):
    """An explicit scientific/input gate stopped the diagnostic workflow."""

    def __init__(self, stage, reason):
        self.stage = stage
        self.reason = str(reason)
        super().__init__(f"{stage}: {self.reason}")


def _stage(stage, fn, *args, **kwargs):
    try:
        return fn(*args, **kwargs)
    except (ValueError, np.linalg.LinAlgError) as exc:
        raise WorkflowBlocked(stage, exc) from exc


@dataclass(frozen=True)
class SceneCandidates:
    frames: tuple
    original_scenes: tuple
    proposed_scenes: tuple
    halo_fits: tuple
    corrections: tuple
    psfs: tuple
    psf_identities: tuple
    status: str = STATUS


def prepare_scene_candidates(frames, scenes, *, halo_settings, compact_sources,
                             psfs, psf_identities):
    """Fit independent sky/halo candidates and propose additional layers.

    halo_settings are explicit ordered per-band keyword mappings for the
    existing fit_native_independent_sky_halo API. Compact-source corrections
    may be None. PSF identities are ordered (input SHA256, grid fingerprint)
    pairs. Neither scene nor FITS is mutated.
    """
    frames, scenes, settings, compact, psfs, identities = map(
        tuple, (frames, scenes, halo_settings, compact_sources, psfs, psf_identities))
    if not frames or any(len(v)!=len(frames) for v in
                         (scenes, settings, compact, psfs, identities)):
        raise WorkflowBlocked("INPUT", "One explicit input set per band required")
    if any(not isinstance(f,NativeBandFrame) for f in frames):
        raise WorkflowBlocked("INPUT", "Validated native frames required")
    if len({f.band for f in frames})!=len(frames):
        raise WorkflowBlocked("INPUT", "Duplicate native band")
    _stage("INPUT",common_native_validity,frames,scenes)
    proposed=[]; fits=[]; corrections=[]
    for f,s,k,c,p,identity in zip(frames,scenes,settings,compact,psfs,identities):
        if (not isinstance(k,Mapping) or "compact_sources" in k
                or (c is not None and (not isinstance(c,NativeCorrection)
                                        or c.kind!="compact_sources"))):
            raise WorkflowBlocked("INPUT", "Explicit halo settings and compact-source correction required")
        if (not isinstance(p,PSFDescription) or p.state!="EMPIRICAL_KERNEL"
                or tuple(identity)!=(f.native_input_sha256,f.native_grid_fingerprint)):
            raise WorkflowBlocked("INPUT", "PSF identity mismatch")
        h=_stage(f"SKY_HALO_{f.band}",fit_native_independent_sky_halo,
                 f,s,compact_sources=c,**dict(k))
        candidate=s
        for layer in h["proposed_layers"]:
            candidate=_stage(f"MASK_PROPOSAL_{f.band}",candidate.with_layer,layer)
        proposed.append(candidate);fits.append(h)
        corrections.append((h["sky"],h["halo"])+(() if c is None else (c,)))
    return SceneCandidates(frames,scenes,tuple(proposed),tuple(fits),
                           tuple(corrections),psfs,tuple(map(tuple,identities)))


def geometry_candidates(prepared, *, masks_confirmed, anchor_icrs, region,
                        region_grid_id, geometry_settings):
    """Fit geometry hypotheses only after explicit review of proposed masks."""
    if not isinstance(prepared,SceneCandidates) or masks_confirmed is not True:
        raise WorkflowBlocked("MASK_REVIEW", "Review and confirm proposed stellar masks first")
    return _stage("GEOMETRY",fit_native_geometry_candidates,
        prepared.frames,prepared.proposed_scenes,anchor_icrs=anchor_icrs,
        corrections_by_band=prepared.corrections,psfs=prepared.psfs,
        psf_identities=prepared.psf_identities,region=region,
        region_grid_id=region_grid_id,**dict(geometry_settings))


def run_reviewed_annular_workflow(prepared, *, masks_confirmed, geometry_confirmed,
                                   geometry, anchor_icrs, annular_settings,
                                   covariance_settings, fit_settings,
                                   output=None, render=False, asymptotic_settings=None,
                                   publication=False):
    """Measure, estimate conditional covariance, fit and optionally export.

    geometry is an explicitly reviewed NestedEllipse (which may come from
    geometry_candidates and reviewed_geometry_candidate). There is no
    automatic fallback to another mask, radial limit, noise or geometry
    hypothesis. A failed or boundary/rank-deficient fit is not exported
    as a successful result. Training/holdout statistics are conditional.
    """
    if not isinstance(prepared,SceneCandidates) or masks_confirmed is not True:
        raise WorkflowBlocked("MASK_REVIEW", "Explicit review of the proposed scene required")
    if geometry_confirmed is not True or not isinstance(geometry,NestedEllipse):
        raise WorkflowBlocked("GEOMETRY_REVIEW", "Explicit reviewed nested geometry required")
    if not isinstance(publication,bool) or (publication and output is None):
        raise WorkflowBlocked("EXPORT", "Six-figure publication export requires an explicit output path")
    frames, scenes=prepared.frames,prepared.proposed_scenes
    common=_stage("COMMON_SUPPORT",common_native_validity,frames,scenes)
    if not common.any():
        raise WorkflowBlocked("COMMON_SUPPORT", "No common valid pixels")
    ann=dict(annular_settings)
    if any(k in ann for k in ("common_valid","retain_pixel_support","anchor_icrs",
                               "corrections","psf","psf_input_sha256","psf_grid_id",
                               "retain_full_aperture_operator")):
        raise WorkflowBlocked("INPUT", "Annular settings cannot override native identities or support")
    records=[]
    for f,s,c,p,identity in zip(frames,scenes,prepared.corrections,
                                  prepared.psfs,prepared.psf_identities):
        rec=_stage(f"ANNULI_{f.band}",prepare_native_annuli,
            f,s,geometry,anchor_icrs=anchor_icrs,corrections=c,psf=p,
            psf_input_sha256=identity[0],psf_grid_id=identity[1],
            common_valid=common,retain_pixel_support=True,
            retain_full_aperture_operator=asymptotic_settings is not None,**ann)
        # Require finite PSF support around measured pixels, separately
        # from the in-frame ring coverage (not a full contour certificate).
        yy,xx=np.nonzero(rec["ring_index"]>=0)
        hy,hx=np.asarray(p.kernel.shape)//2
        if (not len(yy) or yy.min()<hy or xx.min()<hx
                or yy.max()+hy>=f.shape[0] or xx.max()+hx>=f.shape[1]):
            raise WorkflowBlocked(f"ANNULI_{f.band}", "Measured support lacks full finite PSF footprint")
        records.append(rec)
    records=tuple(records)
    tr=np.asarray(fit_settings["train_rings"])
    te=np.asarray(fit_settings["test_rings"])
    nr=len(records[0]["ring_edges_arcsec"])-1
    if (tr.shape!=(nr,) or te.shape!=(nr,) or tr.dtype.kind!="b"
            or te.dtype.kind!="b" or not tr.any() or np.any(tr&te)):
        raise WorkflowBlocked("COVERAGE", "Explicit disjoint training and test rings required")
    if any(np.any((tr|te)&~r["coverage"]["eligible"]) for r in records):
        raise WorkflowBlocked("COVERAGE", "Requested rings fail contiguous coverage; no automatic truncation")
    covariance=_stage("COVARIANCE",native_scene_block_covariance,
                      frames,records,prepared.halo_fits,**dict(covariance_settings))
    fit=_stage("GLS",fit_native_covariance_sersic,records,covariance,**dict(fit_settings))
    if (not fit["success"] or fit["boundary"]
            or fit["conditional_parameter_covariance"] is None):
        raise WorkflowBlocked("GLS", "Fit failed, reached a boundary, or lacks full-rank conditional covariance")
    if te.any() and fit["holdout"] is None:
        raise WorkflowBlocked("HOLDOUT", "Requested predictive check unavailable")
    fit["workflow"]={"status":STATUS,"masks_explicitly_reviewed":True,
                     "geometry_explicitly_reviewed":True,
                     "observed_finite_psf_support":True,
                     "full_geometry_containment_verified":False,
                     "contamination":[{"band":f.band,
                        "sky_coefficients":h["sky_coefficients"].tolist(),
                        "amplitude":h["amplitude"],"sky_condition":h["sky_condition"],
                        "holdout_rmse":h["holdout_rmse"],"block_counts":list(h["counts"])}
                        for f,h in zip(frames,prepared.halo_fits)],
                     "science_adopted":False}
    if asymptotic_settings is not None:
        fit["asymptotic_diagnostic"]=_stage("ASYMPTOTIC",asymptotic_flux_diagnostic,
            records,covariance["covariance"],fit,**dict(asymptotic_settings))
    manifest=None
    if output is not None:
        context=None
        if publication:
            context={"frames":frames,"scenes":scenes,"psfs":prepared.psfs,
                     "corrections_by_band":prepared.corrections}
        manifest=_stage("EXPORT",export_annular_diagnostic,output,records,
                        covariance["covariance"],fit,render=render,publication_context=context)
    return {"status":STATUS,"stage":"DIAGNOSTIC_COMPLETE","records":records,
            "covariance":covariance,"fit":fit,"manifest":manifest,
            "science_adopted":False,"covariance_calibrated":False}
