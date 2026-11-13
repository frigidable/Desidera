"""Stage 04.0–04.1: finite, signed native-pixel aperture flux (diagnostic).

Exact fractional *pixel* overlap comes from photutils. The local on-sky
ellipse -> pixel ellipse projection is validated against WCS boundary
samples, rather than equating astronomical PA with photutils theta.
No data writing, background fitting, mask interpolation, sky adoption,
PSF deconvolution, or total-galaxy flux claim.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, fields
from hashlib import sha256
import json
import math
import re
from types import MappingProxyType

import numpy as np

from lsbpipe.photometry.contracts import NativeBandFrame
from lsbpipe.scene import grid_fingerprint


_SCHEMA = "desidera.photometry.native-aperture/1"
_SHA256 = re.compile(r"^[0-9a-f]{64}$")


def _finite_float(name, value):
    if isinstance(value, (bool, str)) or not isinstance(
            value, (int, float, np.integer, np.floating)):
        raise ValueError(name + " must be a finite number")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(name + " must be finite")
    return result


@dataclass(frozen=True)
class ApertureGeometryICRS:
    ra_deg: float
    dec_deg: float
    semimajor_arcsec: float
    axis_ratio: float = 1.0
    pa_deg_east_of_north: float = 0.0

    def __post_init__(self):
        for key in self.__dataclass_fields__:
            value = _finite_float(key, getattr(self, key))
            object.__setattr__(self, key, value)
        if not 0 <= self.ra_deg < 360 or not -90 <= self.dec_deg <= 90:
            raise ValueError("Invalid ICRS target centre")
        if not 0 < self.semimajor_arcsec <= 3600:
            raise ValueError("Semimajor axis must be 0 < a <= 3600 arcsec")
        if not 0 < self.axis_ratio <= 1:
            raise ValueError("Axis ratio must be 0 < q <= 1")
        if not 0 <= self.pa_deg_east_of_north < 360:
            raise ValueError("Astronomical PA must be 0 <= PA < 360 degrees")


@dataclass(frozen=True)
class ApertureSkyScenario:
    method: str = "zero_additional_correction"
    b0_nMgy_per_pixel: float = 0.0
    bx_nMgy_per_pixel2: float = 0.0
    by_nMgy_per_pixel2: float = 0.0
    reference_xy_pix: tuple[float, float] = (0., 0.)
    evidence_id: str = "DR11_UPSTREAM_SKY_SUBTRACTED"

    def __post_init__(self):
        if self.method not in (
                "zero_additional_correction", "constant", "tilted_plane"):
            raise ValueError("Unsupported residual sky scenario")
        for name in ("b0_nMgy_per_pixel", "bx_nMgy_per_pixel2",
                     "by_nMgy_per_pixel2"):
            object.__setattr__(
                self, name, _finite_float(name, getattr(self, name)))
        if not isinstance(self.reference_xy_pix, (tuple, list)) or (
                len(self.reference_xy_pix) != 2):
            raise ValueError("Sky reference requires two native pixel coordinates")
        xy = tuple(_finite_float("sky reference pixel", v)
                   for v in self.reference_xy_pix)
        object.__setattr__(self, "reference_xy_pix", xy)
        if not isinstance(self.evidence_id, str) or not self.evidence_id:
            raise ValueError("Sky evidence identifier required")
        if self.method == "zero_additional_correction" and any(
                (self.b0_nMgy_per_pixel, self.bx_nMgy_per_pixel2,
                 self.by_nMgy_per_pixel2)):
            raise ValueError("Zero additional sky correction must be identically zero")
        if self.method == "constant" and any(
                (self.bx_nMgy_per_pixel2, self.by_nMgy_per_pixel2)):
            raise ValueError("Constant sky cannot contain a slope")


@dataclass(frozen=True)
class ApertureProvenance:
    """Immutable project/Scene/catalogue lineage; never self-certifies QA."""
    project_id: str
    scene_draft_sha256: str
    catalogue_sha256: str
    source_policy_id: str

    def __post_init__(self):
        if not isinstance(self.project_id, str) or not self.project_id.strip():
            raise ValueError("Aperture must be tied to a galaxy project")
        for key in ("scene_draft_sha256", "catalogue_sha256"):
            v = getattr(self, key)
            if not isinstance(v, str) or _SHA256.fullmatch(v) is None:
                raise ValueError(key + " needs full SHA-256")
        if not isinstance(self.source_policy_id, str) or not self.source_policy_id:
            raise ValueError("Explicit source/membership mask policy required")


@dataclass(frozen=True)
class PixelEllipse:
    x_center: float
    y_center: float
    a_pix: float
    b_pix: float
    theta_rad_pixel_ccw: float
    max_projection_error_pix: float
    local_area_arcsec2_per_native_pixel: float


@dataclass(frozen=True)
class ApertureFluxRecord:
    schema: str
    project_id: str
    dataset_id: str
    grid_fingerprint: str
    source_draft_sha256: str
    catalogue_sha256: str
    source_policy_id: str
    recipe_sha256: str
    sky_method: str
    sky_evidence_id: str
    geometry_icrs: dict
    geometry_native: dict
    area_native_pixel2: dict
    area_by_reason_native_pixel2: dict
    frame_invalid_by_reason_native_pixel2: dict
    sector_qa: tuple[dict, ...]
    signed_valid_flux_nMgy: float | None
    formal_diagonal_sigma_nMgy: float | None
    ab_magnitude_valid_area_only: float | None
    measured_mean_intensity_nMgy_per_arcsec2: float | None
    science_state: str
    adopted_for_science: bool
    covariance_calibrated: bool
    psf_validated: bool
    sky_science_adopted: bool
    warnings: tuple[str, ...]

    def __post_init__(self):
        # frozen=True is only shallow. Freeze nested QA tables as well so
        # callers cannot rewrite a science result after its provenance was
        # computed (to_dict gives a separate, mutable export copy).
        for name in (
                "geometry_icrs", "geometry_native", "area_native_pixel2",
                "area_by_reason_native_pixel2",
                "frame_invalid_by_reason_native_pixel2"):
            object.__setattr__(
                self, name, MappingProxyType(dict(getattr(self, name))))
        object.__setattr__(
            self, "sector_qa",
            tuple(MappingProxyType(dict(item)) for item in self.sector_qa))
        object.__setattr__(self,"warnings",tuple(self.warnings))

    def to_dict(self):
        result = {entry.name:getattr(self,entry.name)
                  for entry in fields(self)}
        for key in ("geometry_icrs","geometry_native",
                    "area_native_pixel2","area_by_reason_native_pixel2",
                    "frame_invalid_by_reason_native_pixel2"):
            result[key] = dict(result[key])
        result["sector_qa"] = [dict(item) for item in self.sector_qa]
        result["warnings"] = list(self.warnings)
        json.dumps(result, allow_nan=False)
        return result


def project_icrs_ellipse_to_native(geometry: ApertureGeometryICRS, wcs):
    """Locally project a true ICRS PA/q ellipse to a photutils pixel ellipse.

    The WCS may be rotated/reflected or cross RA=0. Fail closed where
    its local affine approximation deviates >0.15 px around the aperture.
    """
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    if not isinstance(geometry, ApertureGeometryICRS):
        raise TypeError("ApertureGeometryICRS required")
    if not getattr(wcs, "has_celestial", False):
        raise ValueError("Celestial native WCS required")
    wcs = wcs.celestial
    if getattr(wcs, "has_distortion", False):
        raise ValueError("Distorted WCS requires nonaffine aperture rasterizer")
    target = SkyCoord(ra=geometry.ra_deg*u.deg,
                      dec=geometry.dec_deg*u.deg, frame="icrs")
    x0, y0 = wcs.all_world2pix(geometry.ra_deg, geometry.dec_deg, 0)
    east = target.spherical_offsets_by(1*u.arcsec, 0*u.arcsec)
    north = target.spherical_offsets_by(0*u.arcsec, 1*u.arcsec)
    xe, ye = wcs.all_world2pix(east.ra.deg, east.dec.deg, 0)
    xn, yn = wcs.all_world2pix(north.ra.deg, north.dec.deg, 0)
    xy = np.array([x0, y0], dtype=float)
    xy_e = np.array([xe, ye], dtype=float)-xy
    xy_n = np.array([xn, yn], dtype=float)-xy
    basis = np.column_stack([xy_e, xy_n])
    if not (np.all(np.isfinite(basis)) and np.all(np.isfinite(xy))):
        raise ValueError("Nonfinite sky-to-native WCS Jacobian")
    if abs(np.linalg.det(basis)) < 1.e-12:
        raise ValueError("Degenerate native WCS Jacobian")
    pa = np.deg2rad(geometry.pa_deg_east_of_north)
    a = geometry.semimajor_arcsec
    q = geometry.axis_ratio
    major = a*(np.sin(pa)*xy_e+np.cos(pa)*xy_n)
    minor = a*q*(np.cos(pa)*xy_e-np.sin(pa)*xy_n)
    transform = np.column_stack([major, minor])
    u_matrix, singular, _ = np.linalg.svd(transform)
    if np.any(~np.isfinite(singular)) or singular[1] <= 0:
        raise ValueError("Unresolved ellipse in native WCS")
    if singular[0] > 8192:
        raise ValueError("Aperture too large for bounded exact-pixel integration")
    theta = float(math.atan2(u_matrix[1, 0],u_matrix[0, 0]))
    # 16 explicit on-sky ellipse boundary checkpoints. An approximate
    # local ellipse is not silently treated as exact sky geometry.
    phases = np.linspace(0., 2*np.pi, 16, endpoint=False)
    east_offset = a*(np.sin(pa)*np.cos(phases) +
                     q*np.cos(pa)*np.sin(phases))
    north_offset = a*(np.cos(pa)*np.cos(phases) -
                      q*np.sin(pa)*np.sin(phases))
    perimeter = target.spherical_offsets_by(
        east_offset*u.arcsec, north_offset*u.arcsec)
    real_x, real_y = wcs.all_world2pix(
        perimeter.ra.deg, perimeter.dec.deg, 0)
    model_x = float(x0)+major[0]*np.cos(phases)+minor[0]*np.sin(phases)
    model_y = float(y0)+major[1]*np.cos(phases)+minor[1]*np.sin(phases)
    max_error = float(np.max(np.hypot(
        np.asarray(real_x)-model_x,np.asarray(real_y)-model_y)))
    if not math.isfinite(max_error) or max_error > .15:
        raise ValueError(
            "Projected ICRS ellipse exceeds WCS affine error budget (0.15 px)")
    area = float(1./abs(np.linalg.det(basis)))
    return PixelEllipse(
        x_center=float(x0), y_center=float(y0),
        a_pix=float(singular[0]), b_pix=float(singular[1]),
        theta_rad_pixel_ccw=theta,
        max_projection_error_pix=max_error,
        local_area_arcsec2_per_native_pixel=area)


def _snapshot_native_boolean_mask(name, mask, shape):
    """Own one full-grid GUI mask snapshot in an immutable bytes buffer."""
    values = np.asarray(mask)
    if values.dtype.kind != "b" or values.shape != shape:
        raise ValueError(name + " requires an explicit same-grid boolean array")
    # Copy only once, at the boundary of a complete growth operation.
    # An ordinary read-only NumPy *view* is not an immutable snapshot:
    # its owning GUI array can be mutated or flags reset.
    return np.frombuffer(values.tobytes(order="C"), dtype=np.bool_).reshape(shape)


def _mask_for(name, mask, shape):
    values = np.asarray(mask)
    if values.dtype.kind != "b" or values.shape != shape:
        raise ValueError(name+" requires an explicit same-grid boolean array")
    # Stage04.2 takes an immutable BYTES-backed snapshot once. A byte
    # root cannot be modified through NumPy (including setflags), so
    # each nested aperture may safely reuse it without 2× full-image
    # copies per radius. Ordinary writable masks, including read-only
    # views backed by a writable ndarray, STILL get defensive snapshots.
    parent = values
    safe = not values.flags.writeable
    while safe and isinstance(parent, np.ndarray):
        if parent.flags.writeable:
            safe = False
            break
        parent = parent.base
    if safe and (isinstance(parent, bytes) or
                 (isinstance(parent, memoryview) and
                  parent.readonly and isinstance(parent.obj, bytes))):
        return values
    return np.array(values, dtype=bool, copy=True)


def measure_native_aperture(
        frame: NativeBandFrame, geometry: ApertureGeometryICRS,
        provenance: ApertureProvenance, *,
        effective_source_mask: np.ndarray,
        manual_exclude: np.ndarray,
        sky: ApertureSkyScenario | None = None) -> ApertureFluxRecord:
    """Read-only valid-area aperture sum, no interpolation or total claim.

    Manual restorations MUST be applied upstream to effective_source_mask,
    not to the immutable NativeBandFrame.valid. Target protection is not
    a photometric exclusion. Hard defects always take precedence.
    """
    if not isinstance(frame, NativeBandFrame):
        raise TypeError("Validated NativeBandFrame required")
    if not isinstance(provenance, ApertureProvenance):
        raise TypeError("Verified project Scene provenance required")
    sky = ApertureSkyScenario() if sky is None else sky
    if not isinstance(sky, ApertureSkyScenario):
        raise TypeError("Versioned sky scenario required")
    shape = frame.shape
    source = _mask_for("effective source mask",effective_source_mask,shape)
    manual = _mask_for("manual exclusion",manual_exclude,shape)
    wcs = frame.wcs
    ellipse = project_icrs_ellipse_to_native(geometry,wcs)

    try:
        from photutils.aperture import EllipticalAperture
    except ImportError as exc:
        raise RuntimeError(
            "Exact aperture photometry requires the installed "
            "photutils extra: pip install -e '.[exact]'") from exc
    aperture = EllipticalAperture(
        (ellipse.x_center,ellipse.y_center),
        a=ellipse.a_pix,b=ellipse.b_pix,
        theta=ellipse.theta_rad_pixel_ccw)
    aperture_mask = aperture.to_mask(method="exact")
    large,small = aperture_mask.get_overlap_slices(shape)
    ideal = float(math.pi*ellipse.a_pix*ellipse.b_pix)
    flags = [
        "DIAGNOSTIC_NOT_SCIENCE_ADOPTED",
        "SOURCE_MEMBERSHIP_POLICY_NOT_SCIENTIFICALLY_VALIDATED",
        "RESIDUAL_SKY_SYSTEMATIC_NOT_CALIBRATED",
        "EMPIRICAL_PSF_SCIENCE_QA_NOT_PASSED",
        "COADD_PIXEL_COVARIANCE_UNKNOWN_FORMAL_IVAR_ONLY",
        "SECTOR_BOUNDARIES_BY_PIXEL_CENTRE_DIAGNOSTIC",
        "LOCAL_AFFINE_ICRS_ELLIPSE_APPROXIMATION",
    ]
    if sky.method != "zero_additional_correction":
        flags.append("CONDITIONAL_NONZERO_SKY_SCENARIO_NOT_ADOPTED")
    hashes = frame.source_evidence.get("local_hashes", {})
    if not isinstance(hashes, dict) and not hasattr(hashes,"get"):
        hashes = {}
    if not hashes or any(
            not isinstance(hashes.get(k),str) or
            _SHA256.fullmatch(hashes.get(k,"")) is None
            for k in ("image","invvar","nexp","maskbits")):
        flags.append("FULL_SOURCE_FITS_HASHES_NOT_VERIFIED")
    reasons = {
        "off_frame": 0., "frame_hard_or_invalid": 0.,
        "source_excluded_after_valid": 0.,
        "manual_excluded_after_source": 0.,
    }
    native_reasons = {
        name: 0. for name in (
            "NONFINITE_IMAGE", "BAD_INVERSE_VARIANCE", "BAD_NEXP",
            "SATUR_" + frame.band.upper(),
            "ALLMASK_" + frame.band.upper(),
            "NPRIMARY", "BRIGHT", "GALAXY",
            "USER_EXCLUSION", "OTHER_INVALID")
    }
    sectors = []
    inframe = 0.
    valid_area = 0.
    flux = None
    sigma_diag = None
    mean_intensity = None
    magnitude = None
    if large is not None:
        weights = np.asarray(
            aperture_mask.data[small],dtype=float)
        if np.any(~np.isfinite(weights)) or np.any(
                (weights < -1.e-10) | (weights > 1.+1.e-10)):
            raise ValueError("Invalid exact fractional pixel overlap")
        weights = np.clip(weights,0.,1.)
        valid_native = frame.valid[large]
        source_local = source[large]
        manual_local = manual[large]
        frame_invalid = ~valid_native
        source_ex = valid_native & source_local
        manual_ex = valid_native & ~source_local & manual_local
        usable = valid_native & ~source_local & ~manual_local
        inframe = float(np.sum(weights))
        reasons["frame_hard_or_invalid"] = float(
            np.sum(weights[frame_invalid]))
        remaining = frame_invalid.copy()
        for key in tuple(native_reasons)[:-1]:
            layer = frame.mask_reasons.get(key)
            if layer is None:
                continue
            chosen = remaining & layer[large]
            native_reasons[key] = float(np.sum(weights[chosen]))
            remaining &= ~chosen
        native_reasons["OTHER_INVALID"] = float(
            np.sum(weights[remaining]))
        if not math.isclose(
                sum(native_reasons.values()),
                reasons["frame_hard_or_invalid"],
                abs_tol=1e-7,rel_tol=1e-10):
            raise ValueError("Native invalidity accounting is not exclusive")
        reasons["source_excluded_after_valid"] = float(
            np.sum(weights[source_ex]))
        reasons["manual_excluded_after_source"] = float(
            np.sum(weights[manual_ex]))
        valid_area = float(np.sum(weights[usable]))
        if not math.isclose(
                inframe, valid_area+sum(reasons.values()),
                rel_tol=1.e-10, abs_tol=1.e-7):
            raise ValueError("Overlapping mask priority violates area accounting")
        if valid_area > 0.:
            ys, xs = np.nonzero(usable & (weights > 0.))
            yy = ys+large[0].start
            xx = xs+large[1].start
            w = weights[ys,xs]
            pix = frame.image_nMgy[large][ys,xs]
            ivar = frame.invvar_nMgy2[large][ys,xs]
            xref,yref = sky.reference_xy_pix
            b = (sky.b0_nMgy_per_pixel +
                 sky.bx_nMgy_per_pixel2*(xx-xref) +
                 sky.by_nMgy_per_pixel2*(yy-yref))
            flux = float(np.sum(w*(pix-b),dtype=np.float64))
            variance = float(np.sum(w*w/ivar,dtype=np.float64))
            sigma_diag = float(math.sqrt(variance))
            mean_intensity = float(
                flux/(valid_area*ellipse.local_area_arcsec2_per_native_pixel))
            if flux > 0:
                magnitude = float(22.5-2.5*math.log10(flux))
        # Eight *diagnostic* sectors use pixel-centre membership; flux
        # itself uses exact fractional elliptical pixel overlap.
        yy,xx = np.indices(weights.shape,dtype=float)
        xx += large[1].start-ellipse.x_center
        yy += large[0].start-ellipse.y_center
        # Rebuild local east/north from the same SVD-free WCS sampling.
        # Pixel position angle can have arbitrary handedness.
        from astropy import units as u
        from astropy.coordinates import SkyCoord
        anchor = SkyCoord(ra=geometry.ra_deg*u.deg,
                          dec=geometry.dec_deg*u.deg,frame="icrs")
        east = anchor.spherical_offsets_by(1*u.arcsec,0*u.arcsec)
        north = anchor.spherical_offsets_by(0*u.arcsec,1*u.arcsec)
        ex,ey = wcs.all_world2pix(east.ra.deg,east.dec.deg,0)
        nx,ny = wcs.all_world2pix(north.ra.deg,north.dec.deg,0)
        matrix = np.column_stack([
            [ex-ellipse.x_center, ey-ellipse.y_center],
            [nx-ellipse.x_center, ny-ellipse.y_center]])
        tangent = np.linalg.inv(matrix) @ np.stack(
            [xx.ravel(),yy.ravel()])
        angles = np.mod(np.arctan2(tangent[0],tangent[1]),2*np.pi)
        sec = (np.floor(angles/(np.pi/4)).astype(int)%8).reshape(
            weights.shape)
        for k in range(8):
            selection = sec == k
            sector_in = float(np.sum(weights[selection]))
            sector_valid = float(np.sum(weights[selection & usable]))
            sectors.append({
                "sector_id_pa_45deg": k,
                "inframe_area_pixel2_by_centre_sector": sector_in,
                "valid_area_pixel2_by_centre_sector": sector_valid,
                "valid_fraction_of_inframe_sector": (
                    sector_valid/sector_in if sector_in > 0 else None),
            })
    reasons["off_frame"] = max(0.,float(ideal-inframe))
    frac_frame = min(1.,inframe/ideal)
    frac_usable = min(1.,valid_area/ideal)
    if frac_frame < .999:
        flags.append("APERTURE_PARTLY_OUTSIDE_NATIVE_FITS")
    if frac_usable < .999:
        flags.append("PARTIAL_VALID_APERTURE_NO_MASK_FLUX_RECONSTRUCTION")
    if valid_area <= 0.:
        flags.append("NO_VALID_OBSERVED_APERTURE_PIXELS")
    geometry_dict = asdict(geometry)
    native_dict = asdict(ellipse)
    area = {
        "ideal_native_pixel2": ideal,
        "inside_image_native_pixel2": inframe,
        "valid_native_pixel2": valid_area,
        "outside_image_native_pixel2": reasons["off_frame"],
        "inframe_fraction_of_ideal": frac_frame,
        "valid_fraction_of_ideal": frac_usable,
        "valid_arcsec2_local_projection_approx": (
            valid_area*ellipse.local_area_arcsec2_per_native_pixel),
    }
    meta = {
        "schema":_SCHEMA,"dataset_id":frame.dataset_id,
        "project_id":provenance.project_id,
        "source":asdict(provenance),
        "geometry":geometry_dict,"native":native_dict,
        "sky":asdict(sky),
        "grid_fingerprint":grid_fingerprint(wcs,shape),
        "mask_sha256": {
            "effective_source":sha256(source.tobytes()).hexdigest(),
            "manual_exclude":sha256(manual.tobytes()).hexdigest(),
        },
        "original_fits_hashes":dict(hashes),
        "aperture_method":"photutils_exact_native_pixel_overlap_v1",
    }
    recipe_hash = sha256(json.dumps(
        meta,sort_keys=True,separators=(",",":"),
        allow_nan=False).encode("utf8")).hexdigest()
    status = ("MEASURED_PARTIAL" if valid_area > 0.
              else "UNDERCONSTRAINED_NO_VALID_PIXELS")
    record = ApertureFluxRecord(
        schema=_SCHEMA,
        project_id=provenance.project_id,
        dataset_id=frame.dataset_id,
        grid_fingerprint=meta["grid_fingerprint"],
        source_draft_sha256=provenance.scene_draft_sha256,
        catalogue_sha256=provenance.catalogue_sha256,
        source_policy_id=provenance.source_policy_id,
        recipe_sha256=recipe_hash,
        sky_method=sky.method,
        sky_evidence_id=sky.evidence_id,
        geometry_icrs=geometry_dict,
        geometry_native=native_dict,
        area_native_pixel2=area,
        area_by_reason_native_pixel2=reasons,
        frame_invalid_by_reason_native_pixel2=native_reasons,
        sector_qa=tuple(sectors),
        signed_valid_flux_nMgy=flux,
        formal_diagonal_sigma_nMgy=sigma_diag,
        ab_magnitude_valid_area_only=magnitude,
        measured_mean_intensity_nMgy_per_arcsec2=mean_intensity,
        science_state=status,
        adopted_for_science=False,
        covariance_calibrated=False,
        psf_validated=False,
        sky_science_adopted=False,
        warnings=tuple(flags),
    )
    record.to_dict()
    return record
