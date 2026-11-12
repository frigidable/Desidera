"""Stage 04.2: project-bound, fixed-ellipse exact curves of growth.

All measurements refer to *observed valid pixels* on one unchanged native
coadd. Annuli are independently integrated from exact ellipse annulus
weights, not inferred by subtracting errors or replacing masked data.
This is an experimental diagnostic, never an adopted galaxy total.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from hashlib import sha256
import json
import math
from collections.abc import Mapping

import numpy as np

from .aperture_flux import (
    ApertureGeometryICRS, ApertureProvenance, ApertureSkyScenario,
    measure_native_aperture, project_icrs_ellipse_to_native, _finite_float,
)
from .contracts import NativeBandFrame


@dataclass(frozen=True)
class GrowthHypothesis:
    """Frozen ICRS angle/shape hypothesis; center never re-fitted."""
    name: str
    q: float
    pa_deg: float
    source: str

    def __post_init__(self):
        if not isinstance(self.name,str) or not self.name:
            raise ValueError("Named growth geometry hypothesis required")
        if not isinstance(self.source,str) or not self.source:
            raise ValueError("Hypothesis provenance required")
        # Frozen hypotheses must preserve canonical numeric geometry,
        # never a Boolean masquerading as a valid q=1 or PA=0.
        for name in ("q", "pa_deg"):
            object.__setattr__(
                self, name, _finite_float(name, getattr(self, name)))
        if not 0 < self.q <= 1:
            raise ValueError("Invalid ellipse axis ratio")
        if not 0 <= self.pa_deg < 180:
            raise ValueError("Position angle must be in [0,180)")


@dataclass(frozen=True)
class AnnulusFluxRecord:
    a_inner_arcsec: float
    a_outer_arcsec: float
    exact_geometric_area_pixel2: float
    inframe_area_pixel2: float
    valid_area_pixel2: float
    valid_fraction_of_geometric: float
    signed_flux_nMgy: float | None
    formal_diagonal_sigma_nMgy: float | None
    signed_mean_intensity_nMgy_per_arcsec2: float | None
    ab_surface_brightness_mag_arcsec2: float | None
    no_valid_pixels: bool

    def to_dict(self):
        record=asdict(self)
        json.dumps(record,allow_nan=False)
        return record


def geometry_hypotheses(*, sga_q=None, sga_pa_deg=None):
    """Proposals only; never recenter or declare an SGA member as truth."""
    probes=[GrowthHypothesis("anchor_circle",1.,0.,
                             "PROJECT_TARGET_ANCHOR_NEUTRAL_CIRCLE")]
    if sga_q is not None and sga_pa_deg is not None:
        q=float(sga_q); pa=float(sga_pa_deg)
        if math.isfinite(q) and math.isfinite(pa) and .05 < q <= 1:
            probes.append(GrowthHypothesis(
                "sga_q_pa_candidate",q,pa%180,
                "SELECTED_SGA_SHAPE_ONLY_ANGLE_AND_AXIS_RATIO_NOT_VERIFIED"))
    if len(probes)==1:
        probes.extend([
            GrowthHypothesis("exploratory_q0p7_north",.7,0.,
                             "PREDECLARED_SHAPE_SENSITIVITY_NOT_A_FIT"),
            GrowthHypothesis("exploratory_q0p7_east",.7,90.,
                             "PREDECLARED_SHAPE_SENSITIVITY_NOT_A_FIT"),
        ])
    return tuple(probes)


def predeclared_growth_radii(
        frame:NativeBandFrame, anchor_icrs, *,
        target_protection_native_pixels=None):
    """Anchored same-band radial grid, not a measured galaxy outer radius."""
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    wcs=frame.wcs
    ny,nx=frame.shape
    ra,dec=map(float,anchor_icrs)
    px,py=map(float,wcs.all_world2pix(ra,dec,0))
    if not np.isfinite([px,py]).all() or not (0 <= px < nx and 0 <= py < ny):
        raise ValueError("The verified galaxy project anchor is outside native image")
    # Use minimum pixel-axis extent only as a *hard cap*. The genuine
    # circular/elliptic sky footprint is later checked by coverage.
    edge_margin=min(px+0.5,nx-0.5-px,py+0.5,ny-0.5-py)
    pixel_scale=math.sqrt(frame.local_pixel_area_arcsec2)
    if edge_margin < 12:
        raise ValueError("Insufficient native image margin for radial diagnostics")
    if target_protection_native_pixels is None:
        # A deliberately labelled image-scale fallback, not a light radius.
        scale_arcsec=.35*edge_margin*pixel_scale
        source="NATIVE_FIELD_SCALE_ONLY_NO_GALAXY_EXTENT_INFERRED"
    else:
        if (isinstance(target_protection_native_pixels,bool)
                or not isinstance(target_protection_native_pixels,(int,np.integer))
                or target_protection_native_pixels < 25):
            raise ValueError("Nonempty bound target-protection mask count required")
        scale_arcsec=math.sqrt(
            target_protection_native_pixels/math.pi)*pixel_scale
        source="PROVISIONAL_TARGET_PROTECTION_AREA_EQUIVALENT_NOT_ISOPHOTE"
    # Keep 10 fixed nested radii, with enough center-to-edge space for
    # most on-sky orientations, but DO NOT pretend all radii fully covered.
    fractions=(.12,.20,.30,.42,.56,.72,.9,1.12,1.40,1.75)
    cap=.85*edge_margin*pixel_scale
    if scale_arcsec*fractions[-1] > cap:
        scale_arcsec=cap/fractions[-1]
        source += "_CLIPPED_BY_IMAGE_MARGIN"
    radii=tuple(float(scale_arcsec*v) for v in fractions)
    if not all(math.isfinite(a) and a>0 for a in radii):
        raise ValueError("Invalid frozen radial grid")
    return radii,{
        "strategy":"FROZEN_ANCHOR_RADIAL_GRID_V1",
        "basis":source,
        "radius_scale_arcsec":float(scale_arcsec),
        "edge_cap_arcsec":float(cap),
        "target_protection_is_not_photometry_mask":True,
        "radii_are_not_measured_isophotes":True,
    }


def measure_native_annulus(frame:NativeBandFrame,
                           outer:ApertureGeometryICRS,
                           a_inner_arcsec:float,
                           provenance:ApertureProvenance, *,
                           effective_source_mask:np.ndarray,
                           manual_exclude:np.ndarray,
                           sky:ApertureSkyScenario):
    """Exact independent fractional annular sum, same frozen WCS/Scene.

    An annulus exists only for concentric, same-q, same-PA nested ellipses.
    Photutils computes fractional geometric intersections, not pixel centers.
    """
    from photutils.aperture import EllipticalAnnulus

    if not (isinstance(outer,ApertureGeometryICRS) and
            isinstance(frame,NativeBandFrame) and
            isinstance(provenance,ApertureProvenance) and
            isinstance(sky,ApertureSkyScenario)):
        raise TypeError("Validated native annular science inputs required")
    ain=_finite_float("inner semimajor axis", a_inner_arcsec)
    if not 0 < ain < outer.semimajor_arcsec:
        raise ValueError("Inner semimajor axis must be strictly inside outer")
    # Validate the scientific mask grid BEFORE an early off-frame
    # return. Otherwise a stale/mismatched Scene mask is silently
    # accepted whenever the particular annulus misses the native FITS.
    source=np.asarray(effective_source_mask)
    manual=np.asarray(manual_exclude)
    if (source.shape!=frame.shape or manual.shape!=frame.shape or
        source.dtype.kind!="b" or manual.dtype.kind!="b"):
        raise ValueError("Annulus Scene masks must be same-grid Boolean arrays")
    # The enclosing growth-run already freezes both full-size masks once.
    # Avoid duplicating two entire native coadd grids for every radial ring.
    e=project_icrs_ellipse_to_native(outer,frame.wcs)
    factor=ain/outer.semimajor_arcsec
    aperture=EllipticalAnnulus(
        (e.x_center,e.y_center),a_in=e.a_pix*factor,
        a_out=e.a_pix,b_out=e.b_pix,b_in=e.b_pix*factor,
        theta=e.theta_rad_pixel_ccw)
    mask=aperture.to_mask(method="exact")
    large,small=mask.get_overlap_slices(frame.shape)
    ideal=math.pi*e.a_pix*e.b_pix*(1-factor**2)
    if large is None:
        return AnnulusFluxRecord(
            ain,outer.semimajor_arcsec,ideal,0.,0.,0.,
            None,None,None,None,True)
    weights=np.asarray(mask.data[small],dtype=float)
    # Exact photutils overlap weights must be finite fractions. Silently
    # clipping an invalid -0.2/1.2 weight would change signed flux,
    # surface brightness and the diagonal-variance estimate without
    # revealing a numerical failure to Stage 04.3.
    if (not np.isfinite(weights).all() or
            np.any((weights < -1e-10) | (weights > 1.+1e-10))):
        raise ValueError("Invalid exact annulus fractional overlap weights")
    weights=np.clip(weights,0.,1.)
    valid=(frame.valid[large]&~source[large]&~manual[large])
    used=valid&(weights>0)
    weighted_area=float(np.sum(weights[used]))
    inframe=float(np.sum(weights))
    if weighted_area<=0.:
        return AnnulusFluxRecord(
            ain,outer.semimajor_arcsec,ideal,inframe,0.,0.,
            None,None,None,None,True)
    dy,dx=np.nonzero(used)
    x=dx+large[1].start
    y=dy+large[0].start
    w=weights[dy,dx]
    pixels=frame.image_nMgy[large][dy,dx]
    invvar=frame.invvar_nMgy2[large][dy,dx]
    xref,yref=sky.reference_xy_pix
    b=(sky.b0_nMgy_per_pixel+
       sky.bx_nMgy_per_pixel2*(x-xref)+
       sky.by_nMgy_per_pixel2*(y-yref))
    flux=float(np.sum(w*(pixels-b),dtype=np.float64))
    sigma=float(np.sqrt(np.sum(w*w/invvar,dtype=np.float64)))
    intensity=float(flux/(weighted_area*e.local_area_arcsec2_per_native_pixel))
    mu=(float(22.5-2.5*math.log10(intensity))
        if intensity>0 else None)
    return AnnulusFluxRecord(
        ain,outer.semimajor_arcsec,float(ideal),float(inframe),
        weighted_area,float(min(1.,weighted_area/ideal)),
        flux,sigma,intensity,mu,False)


def measure_growth_diagnostic(
        frame:NativeBandFrame, *,
        anchor_icrs:tuple[float,float],
        provenance:ApertureProvenance,
        effective_source_mask:np.ndarray,
        manual_exclude:np.ndarray,
        radii_arcsec:tuple[float,...],
        hypothesis:GrowthHypothesis,
        sky:ApertureSkyScenario|None=None):
    """Fixed-geometry finite cumulative flux and independent rings.

    A report never estimates total light or extrapolates missing pixels.
    """
    if not isinstance(hypothesis,GrowthHypothesis):
        raise TypeError("Explicit frozen geometry hypothesis required")
    if not 2 <= len(radii_arcsec) <= 40:
        raise ValueError("Growth curve needs 2–40 radii")
    radii=tuple(_finite_float("growth radius",v) for v in radii_arcsec)
    if (any(not math.isfinite(v) or v<=0 for v in radii)
            or any(b<=a for a,b in zip(radii[:-1],radii[1:]))):
        raise ValueError("Growth radii must be strictly increasing")
    sky=sky if sky is not None else ApertureSkyScenario()
    # Explicitly copy user-editable draft masks once for the full run.
    source=np.array(effective_source_mask,copy=True)
    manual=np.array(manual_exclude,copy=True)
    geometry_base={
        "ra_deg":float(anchor_icrs[0]),"dec_deg":float(anchor_icrs[1]),
        "axis_ratio":hypothesis.q,
        "pa_deg_east_of_north":hypothesis.pa_deg,
    }
    cumulative=[]
    annuli=[]
    for i,outer in enumerate(radii):
        geom=ApertureGeometryICRS(
            semimajor_arcsec=outer,**geometry_base)
        cumulative.append(measure_native_aperture(
            frame,geom,provenance,effective_source_mask=source,
            manual_exclude=manual,sky=sky))
        if i:
            annuli.append(measure_native_annulus(
                frame,geom,radii[i-1],provenance,
                effective_source_mask=source,manual_exclude=manual,
                sky=sky))
    residuals=[]
    n_verified_closure=0
    for i,ann in enumerate(annuli,1):
        outer=cumulative[i].signed_valid_flux_nMgy
        inner=cumulative[i-1].signed_valid_flux_nMgy
        if outer is None or inner is None:
            residuals.append(None)
            continue
        # A fully masked ring is not a *measured* zero-flux ring.
        # Nevertheless two measured nested cumulative apertures must
        # then integrate the same original usable pixels.
        ring_flux = ann.signed_flux_nMgy
        if ring_flux is None and not ann.no_valid_pixels:
            raise ValueError("Stage 04.2 annulus missing measured flux")
        deviation=float(outer-inner-(ring_flux if ring_flux is not None else 0.))
        scale=max(1.,abs(outer),abs(inner),
                  abs(ring_flux) if ring_flux is not None else 0.)
        # An exact geometric annulus must equal a difference of
        # concentric exact apertures. This is a numerical identity,
        # not a model-dependent science tolerance or chi^2 gate.
        if (not math.isfinite(deviation) or
                abs(deviation) > 1.e-6 + 1.e-9*scale):
            raise ValueError(
                "Stage 04.2 independent annulus and cumulative flux "
                "closure mismatch at radial index " + str(i))
        n_verified_closure+=1
        residuals.append(deviation if ring_flux is not None else None)
    # CoG is not forcibly monotonic. A sequence of entirely masked
    # apertures is UNDERCONSTRAINED, not a valid partial measurement.
    n_measured=sum(
        rec.signed_valid_flux_nMgy is not None for rec in cumulative)
    complete=n_measured >= 2
    qa={
        "schema":"desidera.photometry.native-growth-diagnostic/1",
        "project_id":provenance.project_id,
        "dataset_id":frame.dataset_id,
        "grid_fingerprint":cumulative[0].grid_fingerprint,
        "native_input_sha256":frame.native_input_sha256,
        "native_input_sha_is_local_array_identity_not_fits_authenticity":True,
        "draft_sha256":provenance.scene_draft_sha256,
        "catalogue_sha256":provenance.catalogue_sha256,
        "source_policy_id":provenance.source_policy_id,
        "geometry_hypothesis":asdict(hypothesis),
        "fixed_project_anchor_icrs":list(anchor_icrs),
        "radii_arcsec":list(radii),
        "sky_method":sky.method,
        "sky_evidence_id":sky.evidence_id,
        "cumulative":[r.to_dict() for r in cumulative],
        "independently_integrated_annuli":[r.to_dict() for r in annuli],
        "annulus_minus_cumulative_difference_nMgy":residuals,
        "max_abs_annulus_closure_nMgy":(
            max(abs(v) for v in residuals if v is not None)
            if any(v is not None for v in residuals) else None),
        "n_measured_direct_apertures":int(n_measured),
        "n_measured_independent_annuli":int(sum(
            rec.signed_flux_nMgy is not None for rec in annuli)),
        "n_comparable_nested_radial_closures_verified":int(
            n_verified_closure),
        "all_comparable_annular_flux_closures_within_numerical_tolerance":
            True,
        "closure_absolute_tolerance_nMgy":1.e-6,
        "closure_relative_tolerance":1.e-9,
        "closure_is_not_a_noise_or_science_quality_test":True,
        "independent_coadd_covariance_calibrated":False,
        "masked_galaxy_light_restored":False,
        "aperture_total_magnitude_validated":False,
        "source_membership_scientifically_validated":False,
        "science_status":(
            "DIAGNOSTIC_MEASURED_PARTIAL_NOT_ADOPTED" if complete
            else "UNDERCONSTRAINED_INSUFFICIENT_VALID_RADIAL_SAMPLES"),
        "adopted_for_science":False,
        "warning_flags":([
            "NO_USABLE_RADIAL_CURVE_MEASUREMENT"
        ] if not complete else []) + [
            "FIXED_ELLIPSE_HYPOTHESIS_NOT_MEASURED_ISOPHOTE",
            "NEIGHBOUR_MEMBERSHIP_PROVISIONAL",
            "COADD_NOISE_COVARIANCE_UNCALIBRATED",
            "DR11_UPSTREAM_SKY_SUBTRACTION_BIAS_UNKNOWN",
            "NO_MASKED_FLUX_RECONSTRUCTION",
        ],
    }
    qa["receipt_sha256"]=sha256(json.dumps(
        {k:v for k,v in qa.items() if k!="receipt_sha256"},
        sort_keys=True,separators=(",",":"),allow_nan=False).encode()).hexdigest()
    json.dumps(qa,allow_nan=False)
    return qa
