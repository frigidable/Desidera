"""Differential native-pixel injections for non-adopted recovery tests.

Add an explicitly known PSF-convolved signal to a separate in-memory
frame. Never edit FITS, change native masks or manufacture new noise.
A differential recovery is not an independent observation or a full
injection/recovery calibration.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
import numpy as np
from scipy.signal import fftconvolve

from .contracts import NativeBandFrame, PSFDescription
from .radial_models import NestedEllipse, sersic_intensity, sersic_model_total
from .native_annular_diagnostic import NativeCorrection, numeric_sha256
from .native_workflow import prepare_scene_candidates, run_reviewed_annular_workflow
from .uncertainty_budget import STATUS


@dataclass(frozen=True)
class NativeInjection:
    frame: NativeBandFrame
    parent_sha256: str
    signal_nMgy_per_pixel: np.ndarray
    metadata: dict
    status: str = STATUS


def inject_native_sersic(frame, psf, *, psf_identity, geometry, anchor_icrs,
                         profile, max_work_pixels=500000):
    """Inject positive intrinsic Sérsic intensity sampled at native pixels.

    Profile = (central intensity nMgy/pixel, scale arcsec, index).
    The entire finite model raster is convolved before selecting pixels.
    Nonfinite original pixels remain nonfinite; all native masks and
    IVAR are unchanged. No extra photon noise or exposure is simulated.
    """
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    from hashlib import sha256
    if not isinstance(frame,NativeBandFrame) or not isinstance(geometry,NestedEllipse):
        raise ValueError("Validated native frame and explicit geometry required")
    grid=frame.native_grid_fingerprint
    parent=frame.native_input_sha256
    if (not isinstance(psf,PSFDescription) or psf.state!="EMPIRICAL_KERNEL"
            or tuple(psf_identity)!=(parent,grid)):
        raise ValueError("Injection PSF must refer to the parent native observation")
    if (not isinstance(max_work_pixels,int) or max_work_pixels<1
            or np.prod(frame.shape)>max_work_pixels):
        raise ValueError("Explicit injection work-pixel budget exceeded")
    ra,dec=map(float,anchor_icrs)
    if not np.isfinite([ra,dec]).all() or not 0<=ra<360 or not -90<=dec<=90:
        raise ValueError("Invalid injection ICRS anchor")
    amp,re,n=map(float,profile)
    center=SkyCoord(ra*u.deg,dec*u.deg,frame="icrs")
    yy,xx=np.indices(frame.shape)
    coords=frame.wcs.pixel_to_world(xx,yy).icrs
    sep=center.separation(coords)
    if not np.isfinite(sep.deg).all() or np.max(sep.deg)>.5:
        raise ValueError("Small-field injection restricted to 0.5 degree")
    angle=center.position_angle(coords).rad
    radius=geometry.radius(sep.arcsec*np.sin(angle),sep.arcsec*np.cos(angle))
    intrinsic=sersic_intensity(radius,amp,re,n)
    signal=fftconvolve(intrinsic,psf.kernel,mode="same")
    # Numerical FFT noise only. Do not clip negative observational pixels.
    tolerance=1e-12*amp
    if not np.isfinite(signal).all() or np.min(signal)<-tolerance:
        raise ValueError("Invalid injected PSF model")
    signal=np.maximum(signal,0.)
    image=np.array(frame.image_nMgy,copy=True)
    finite=np.isfinite(image)
    image[finite]+=signal[finite]
    if not np.isfinite(image[finite]).all():
        raise ValueError("Injection overflow")
    model=sersic_model_total(geometry,amp/frame.local_pixel_area_arcsec2,re,n)
    metadata={"kind":"DIFFERENTIAL_SYNTHETIC_SIGNAL_ON_NATIVE_FRAME",
              "parent_input_sha256":parent,"parent_dataset_id":frame.dataset_id,
              "grid_id":grid,"psf_sha256":sha256(psf.kernel.tobytes()).hexdigest(),
              "anchor_icrs":[ra,dec],"geometry":asdict(geometry),
              "profile":[amp,re,n],"profile_units":["nMgy/pixel","arcsec","dimensionless"],
              "signal_sha256":numeric_sha256(signal),
              "intrinsic_model_total_nMgy":model["model_total"],
              "finite_intrinsic_raster_sum_nMgy":float(intrinsic.sum()),
              "finite_convolved_raster_sum_nMgy":float(signal.sum()),
              "additional_noise_simulated":False,"science_adopted":False}
    child=NativeBandFrame(image,frame.invvar_nMgy2,frame.nexp,frame.maskbits,frame.wcs,
                          release=frame.release,region=frame.region,brick=frame.brick,
                          band=frame.band,policy=frame.policy,
                          extra_exclusion=frame.extra_exclusion,
                          source_evidence={"synthetic_injection":metadata})
    if (child.native_grid_fingerprint!=grid or not np.array_equal(child.valid,frame.valid)
            or not np.array_equal(child.maskbits,frame.maskbits)):
        raise ValueError("Injection changed native observation support")
    return NativeInjection(child,parent,
        np.frombuffer(signal.tobytes(),dtype=signal.dtype).reshape(signal.shape),metadata)


def run_differential_sersic_injection(frames, scenes, *, psfs, psf_identities,
                                      halo_settings, compact_sources,
                                      injection_geometry, injection_anchor_icrs,
                                      injection_profiles, workflow_settings,
                                      max_work_pixels=500000):
    """Run matched original/injected diagnostics with independent sky/halo refits.

    Both runs use the same reviewed mask hypotheses, geometry, ring split
    and covariance settings, while covariance and halo/sky coefficients
    are re-estimated. The geometry and compact-source amplitudes stay
    frozen. This is a local differential response, not recovery of a
    separately identifiable physical component. No archives or FITS are
    written by this experiment.
    """
    frames,scenes,psfs,ids,hs,cs,profiles=map(tuple,
        (frames,scenes,psfs,psf_identities,halo_settings,compact_sources,injection_profiles))
    if not frames or any(len(x)!=len(frames) for x in (scenes,psfs,ids,hs,cs,profiles)):
        raise ValueError("One explicit matched injection specification per band required")
    options=dict(workflow_settings)
    if any(k in options for k in ("output","render")):
        raise ValueError("Differential experiment does not write science archives")
    if options.get("masks_confirmed") is not True or options.get("geometry_confirmed") is not True:
        raise ValueError("Explicit review of the shared mask and geometry hypotheses required")
    original=prepare_scene_candidates(frames,scenes,halo_settings=hs,
        compact_sources=cs,psfs=psfs,psf_identities=ids)
    baseline=run_reviewed_annular_workflow(original,**options)
    injections=tuple(inject_native_sersic(f,p,psf_identity=identity,
        geometry=injection_geometry,anchor_icrs=injection_anchor_icrs,profile=profile,
        max_work_pixels=max_work_pixels)
        for f,p,identity,profile in zip(frames,psfs,ids,profiles))
    children=tuple(i.frame for i in injections)
    child_ids=tuple((f.native_input_sha256,f.native_grid_fingerprint) for f in children)
    child_cs=tuple(None if c is None else NativeCorrection(
        c.kind,f.native_grid_fingerprint,f.native_input_sha256,
        f"Frozen parent correction {c.sha256}: {c.source}",c.values)
        for c,f in zip(cs,children))
    injected=prepare_scene_candidates(children,scenes,halo_settings=hs,
        compact_sources=child_cs,psfs=psfs,psf_identities=child_ids)
    after=run_reviewed_annular_workflow(injected,**options)
    rows={}
    for i,f in enumerate(frames):
        a,b=baseline["records"][i],after["records"][i]
        if (a["ring_index_sha256"]!=b["ring_index_sha256"]
                or a["valid_sha256"]!=b["valid_sha256"]):
            raise ValueError("Injection altered measurement support; differential comparison blocked")
        if (not np.array_equal(a["ring_index"],b["ring_index"])
                or a["geometry"]!=b["geometry"]):
            raise ValueError("Injection changed ring geometry")
        ring=a["ring_index"]
        use=ring>=0
        nr=len(a["observed_flux_nMgy"])
        truth=np.bincount(ring[use],weights=injections[i].signal_nMgy_per_pixel[use],minlength=nr)
        recovered=b["observed_flux_nMgy"]-a["observed_flux_nMgy"]
        delta_model=(after["fit"]["model_extrapolations"][f.band]["model_total"]
                     -baseline["fit"]["model_extrapolations"][f.band]["model_total"])
        total=injections[i].metadata["intrinsic_model_total_nMgy"]
        rows[f.band]={"expected_valid_ring_flux_nMgy":truth,
                      "recovered_ring_flux_change_nMgy":recovered,
                      "ring_flux_bias_nMgy":recovered-truth,
                      "intrinsic_injected_model_total_nMgy":total,
                      "fitted_model_total_change_nMgy":delta_model,
                      "fitted_model_total_differential_fractional_bias":float(delta_model/total-1),
                      "sky_coefficients_change":injected.halo_fits[i]["sky_coefficients"]-original.halo_fits[i]["sky_coefficients"],
                      "halo_amplitude_change":injected.halo_fits[i]["amplitude"]-original.halo_fits[i]["amplitude"],
                      "injection_metadata":injections[i].metadata}
    return {"status":STATUS,"kind":"DIFFERENTIAL_PIXEL_INJECTION_NOT_INDEPENDENT_RECOVERY",
            "band_order":tuple(f.band for f in frames),"bands":rows,
            "baseline_record_sha256":baseline["fit"]["record_sha256"],
            "injected_record_sha256":after["fit"]["record_sha256"],
            "geometry_refitted":False,"compact_sources_refitted":False,
            "sky_halo_covariance_refitted":True,"noise_realizations":0,
            "covariance_calibrated":False,"science_adopted":False}
