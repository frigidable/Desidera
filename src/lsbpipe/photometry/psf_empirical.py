"""PSF-0: experimental, independently held-out empirical *coadd* PSF.

Operates exclusively on an existing NativeBandFrame (strict original FITS
adapter), not a remote preview. Native pixels are never written or changed.
Local stellar cutout sky is distinct from the LSB galaxy sky solution.

No measurement here is an adopted, spatially validated, wide-angle or
SBF-grade PSF. In particular, small-cutout EPSFs do not constrain wings.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from hashlib import sha256
import json
import math

import numpy as np
from scipy.spatial import cKDTree

from lsbpipe.photometry.contracts import NativeBandFrame
from lsbpipe.photometry.workflow_dependencies import PSFLineageKey
from lsbpipe.sky_geometry import ellipse_from_wcs

STATUS = "EXPERIMENTAL_PSF_NOT_ADOPTED"
RECIPE = "desidera.psf0-coadd-epsf/1"


@dataclass(frozen=True)
class PSF0Config:
    """Reproducible, conservative core-PSF recipe (not a science calibration)."""
    guessed_fwhm_pix: float = 4.8
    detection_sigma: float = 7.
    cutout_size_pix: int = 31
    sky_annulus_inner_pix: float = 10.
    sky_annulus_outer_pix: float = 15.
    max_invalid_fraction: float = .01
    min_neighbor_separation_pix: float = 32.
    min_train_stars: int = 6
    min_holdout_stars: int = 2
    holdout_every: int = 5
    max_training_stars: int = 75
    max_iterations: int = 8
    psf_oversampling: int = 1
    max_star_gradient_sigma: float = 6.
    seed: int = 20261103

    def __post_init__(self):
        floats = (self.guessed_fwhm_pix, self.detection_sigma,
                  self.sky_annulus_inner_pix, self.sky_annulus_outer_pix,
                  self.max_invalid_fraction, self.min_neighbor_separation_pix,
                  self.max_star_gradient_sigma)
        if not all(math.isfinite(x) for x in floats):
            raise ValueError("Nonfinite PSF recipe value")
        if not 2.5 <= self.guessed_fwhm_pix <= 20:
            raise ValueError("PSF-0 requires sampled stars with guessed FWHM >=2.5 px")
        if not 3 <= self.detection_sigma <= 30:
            raise ValueError("Detection threshold out of conservative range")
        if type(self.cutout_size_pix) is not int or not 15 <= self.cutout_size_pix <= 91 or self.cutout_size_pix % 2 != 1:
            raise ValueError("Star cutout must be odd and 15..91 px")
        h = self.cutout_size_pix // 2
        if not 2 < self.sky_annulus_inner_pix < self.sky_annulus_outer_pix <= h:
            raise ValueError("Local stellar sky annulus must lie inside cutout")
        if not 0 <= self.max_invalid_fraction <= .1 or self.min_neighbor_separation_pix < self.cutout_size_pix:
            raise ValueError("Invalid pixel fraction or unsafe neighbor isolation")
        if (type(self.min_train_stars) is not int or self.min_train_stars < 4 or
                type(self.min_holdout_stars) is not int or self.min_holdout_stars < 2 or
                type(self.holdout_every) is not int or self.holdout_every < 3 or
                type(self.max_training_stars) is not int or
                self.max_training_stars < self.min_train_stars or
                type(self.max_iterations) is not int or not 2 <= self.max_iterations <= 25 or
                type(self.seed) is not int):
            raise ValueError("Invalid deterministic star split/build settings")
        if self.psf_oversampling != 1:
            raise ValueError("PSF-0 currently supports only a native-grid EPSF; oversampling requires separate calibration")
        if self.max_star_gradient_sigma <= 0:
            raise ValueError("Gradient rejection threshold must be positive")

    @property
    def digest(self):
        return sha256(json.dumps(asdict(self), sort_keys=True,
                                 allow_nan=False).encode("utf8")).hexdigest()


@dataclass(frozen=True)
class PSF0Star:
    """Audit record: selection, background, and sky coordinates of a star."""
    star_id: str
    x_pix: float
    y_pix: float
    ra_deg: float
    dec_deg: float
    peak_nMgy: float
    cutout_sky_nMgy: float
    cutout_rms_nMgy: float
    gradient_sigma_span: float
    role: str

    def __post_init__(self):
        if self.role not in ("train", "holdout"):
            raise ValueError("Star role must be declared")


def _immutable_kernel(value):
    arr = np.asarray(value, dtype=np.float64)
    if (arr.ndim != 2 or min(arr.shape) < 15 or
            any(n % 2 != 1 for n in arr.shape) or
            not np.all(np.isfinite(arr)) or np.any(arr < 0) or
            not np.isclose(float(arr.sum()), 1., atol=1e-9, rtol=1e-8)):
        raise ValueError("Experimental PSF kernel must be finite, nonnegative, odd and normalized")
    return np.frombuffer(arr.tobytes(order="C"), dtype=np.float64).reshape(arr.shape)


@dataclass(frozen=True)
class PSF0Result:
    dataset_id: str
    instrument: str
    kernel: np.ndarray = field(repr=False, compare=False)
    training: tuple[PSF0Star, ...]
    holdout: tuple[PSF0Star, ...]
    rejected_reasons: dict
    detection_count: int
    metrics: dict
    flags: tuple[str, ...]
    recipe: PSF0Config
    source_hashes: dict
    lineage_id: str | None
    center_target_icrs: tuple[float, float]
    target_protection_arcsec: float
    algorithm_versions: dict
    status: str = STATUS

    def __post_init__(self):
        if self.status != STATUS:
            raise ValueError("PSF-0 cannot produce adopted PSF results")
        object.__setattr__(self, "kernel", _immutable_kernel(self.kernel))
        object.__setattr__(self, "training", tuple(self.training))
        object.__setattr__(self, "holdout", tuple(self.holdout))
        object.__setattr__(self, "flags", tuple(self.flags))
        if (set(s.star_id for s in self.training) &
                set(s.star_id for s in self.holdout)):
            raise ValueError("PSF train/holdout contamination")
        if len(self.training) < self.recipe.min_train_stars or len(self.holdout) < self.recipe.min_holdout_stars:
            raise ValueError("PSF lacks disjoint holdout/training stars")
        if self.lineage_id is not None and (
                len(self.lineage_id) != 64 or "SYNTHETIC_TEST_INPUT" in self.flags):
            raise ValueError("Invalid science provenance for PSF lineage")

    def report(self):
        return {
            "schema": "desidera.psf0-experimental/1",
            "status": self.status,
            "science_adopted": False,
            "dataset_id": self.dataset_id, "instrument": self.instrument,
            "kernel_shape": list(self.kernel.shape),
            "kernel_native_pixel_sampling": True,
            "kernel_sum": float(self.kernel.sum()),
            "kernel_sha256": sha256(self.kernel.tobytes(order="C")).hexdigest(),
            "source_hashes": dict(self.source_hashes),
            "lineage_id": self.lineage_id,
            "recipe": asdict(self.recipe),
            "recipe_digest": self.recipe.digest,
            "target_icrs_deg": list(self.center_target_icrs),
            "target_protection_arcsec": self.target_protection_arcsec,
            "n_detected": self.detection_count,
            "n_train": len(self.training), "n_holdout": len(self.holdout),
            "training": [asdict(s) for s in self.training],
            "holdout": [asdict(s) for s in self.holdout],
            "rejections": dict(self.rejected_reasons),
            "diagnostic_metrics": dict(self.metrics),
            "flags": list(self.flags),
            "algorithm_versions": dict(self.algorithm_versions),
            "psf_wings_validated": False,
            "psf_spatial_representativity_validated": False,
            "noise_covariance_validated": False,
            "stellar_contamination_model_validated": False,
        }


def _verified_source_hashes(frame):
    info = frame.source_evidence
    hashes = info.get("local_hashes", {})
    if (info.get("hash_provenance") != "LOCAL_SHA256_NOT_REMOTE_AUTHENTICITY"
            or info.get("synthetic", False) or not isinstance(hashes, dict)
            or set(hashes) != {"image", "invvar", "nexp", "maskbits"}
            or any(not isinstance(v, str) or len(v) != 64 or
                   any(char not in "0123456789abcdef" for char in v)
                   for v in hashes.values())):
        return None
    return dict(hashes)


def _target_exclusion(frame, center_icrs, radius_arcsec):
    if (len(center_icrs) != 2 or not all(math.isfinite(float(c)) for c in center_icrs)
            or not math.isfinite(radius_arcsec) or not 0 < radius_arcsec <= 1800):
        raise ValueError("Target ICRS and conservative protection radius required")
    ra, dec = map(float, center_icrs)
    if not 0 <= ra < 360 or not -90 <= dec <= 90:
        raise ValueError("Invalid ICRS coordinates")
    ellipse = ellipse_from_wcs(frame.wcs, ra, dec, sky_q=1., sky_pa_deg=0.)
    ny, nx = frame.shape
    yy, xx = np.indices((ny, nx))
    # Local-affine approximation for one survey brick/cutout only, not
    # a spherical sky polygon. Guard this limitation via flags.
    return ellipse.pixel.elliptical_radius(xx, yy) < radius_arcsec * ellipse.pixel_a_per_sky_arcsec


def _star_cutout(frame, x, y, cfg, protected):
    """Fit a local planar background on the annulus of this one star."""
    h = cfg.cutout_size_pix // 2
    ix, iy = int(np.rint(x)), int(np.rint(y))
    ny, nx = frame.shape
    if ix-h < 0 or iy-h < 0 or ix+h >= nx or iy+h >= ny:
        return None, "CUTOUT_EDGE"
    rows = slice(iy-h, iy+h+1)
    cols = slice(ix-h, ix+h+1)
    valid = np.array(frame.valid[rows, cols], dtype=bool, copy=True)
    for key in (f"SATUR_{frame.band.upper()}", f"ALLMASK_{frame.band.upper()}"):
        if key in frame.mask_reasons:
            valid &= ~frame.mask_reasons[key][rows, cols]
    valid &= ~protected[rows, cols]
    if (not valid[h, h] or
            np.count_nonzero(~valid) / valid.size > cfg.max_invalid_fraction):
        return None, "PIXEL_QUALITY_OR_TARGET_OVERLAP"
    src = np.asarray(frame.image_nMgy[rows, cols], dtype=float)
    ivar = np.asarray(frame.invvar_nMgy2[rows, cols], dtype=float)
    gy, gx = np.indices(src.shape)
    dx, dy = gx-(x-(ix-h)), gy-(y-(iy-h))
    rr = np.hypot(dx, dy)
    annulus = (valid & (rr >= cfg.sky_annulus_inner_pix)
               & (rr < cfg.sky_annulus_outer_pix))
    if np.count_nonzero(annulus) < 35:
        return None, "INSUFFICIENT_LOCAL_SKY"
    # Three robust least-squares iterations, independent of the LSB
    # galaxy-wide background. Even accepted fit is PSF-wing sensitive.
    design = np.column_stack((np.ones(np.count_nonzero(annulus)),
                              dx[annulus], dy[annulus]))
    observed = src[annulus]
    include = np.ones(len(observed), dtype=bool)
    for _ in range(4):
        if include.sum() < 25:
            return None, "INSUFFICIENT_LOCAL_SKY"
        coeff = np.linalg.lstsq(design[include], observed[include], rcond=None)[0]
        resid = observed - design @ coeff
        scale = 1.4826 * np.median(np.abs(resid[include] - np.median(resid[include])))
        if not np.isfinite(scale) or scale <= 0:
            scale = np.std(resid[include])
        if not np.isfinite(scale) or scale <= 0:
            return None, "NO_LOCAL_NOISE_ESTIMATE"
        include = np.abs(resid-np.median(resid[include])) < 3.5*scale
    full_sky = coeff[0] + coeff[1]*dx + coeff[2]*dy
    residual = src - full_sky
    core = (rr <= max(2.5, cfg.guessed_fwhm_pix))
    star_peak = float(np.max(residual[valid & core])) if np.any(valid & core) else float("nan")
    if not np.isfinite(star_peak) or star_peak <= cfg.detection_sigma*scale:
        return None, "LOW_LOCAL_STAR_CONTRAST"
    gradient_span = float(np.hypot(coeff[1], coeff[2]) * 2*h / scale)
    if gradient_span > cfg.max_star_gradient_sigma:
        return None, "STRONG_LOCAL_SKY_GRADIENT"
    # A nonzero masked cutout cannot silently become a reconstructed
    # original pixel: those samples receive zero statistical weight.
    cutout = np.where(valid, residual, 0.)
    weights = np.where(valid, np.sqrt(np.maximum(ivar, 0.)), 0.)
    if not np.isfinite(weights).all():
        return None, "BAD_CUTOUT_WEIGHTS"
    return {
        "image": cutout, "weights": weights,
        "center_in_cutout": (x-(ix-h), y-(iy-h)),
        "sky_nMgy": float(coeff[0]), "rms_nMgy": float(scale),
        "gradient_sigma_span": gradient_span,
        "peak_nMgy": star_peak,
        "valid": valid,
    }, None


def _holdout_diagnostic(star, cutout, kernel):
    """Native-pixel weighted flux fit: diagnostic only, not calibrated chi^2."""
    from scipy.ndimage import shift
    cy, cx = (np.array(kernel.shape)-1)/2.
    local_x, local_y = cutout["center_in_cutout"]
    model = shift(kernel, (local_y-cy, local_x-cx),
                  order=3, mode="constant", cval=0., prefilter=True)
    usable = cutout["weights"] > 0
    # Exclude any interpolation ringing: tiny signed numerical artifacts
    # outside the model support are *not* interpreted as stellar flux.
    model = np.maximum(model, 0.)
    numerator = float(np.sum(cutout["weights"][usable]**2 *
                             model[usable] * cutout["image"][usable]))
    denominator = float(np.sum(cutout["weights"][usable]**2 * model[usable]**2))
    if denominator <= 0 or not math.isfinite(denominator):
        raise ValueError("Holdout PSF normalization is undefined")
    flux = numerator/denominator
    if flux <= 0 or not math.isfinite(flux):
        raise ValueError("Holdout star fitted flux is invalid")
    err = cutout["image"]-flux*model
    inner = np.indices(kernel.shape)
    r = np.hypot(inner[1]-cx, inner[0]-cy)
    local = usable & (r <= min(cx, cy, 9))
    if not np.any(local):
        raise ValueError("Empty PSF holdout aperture")
    return {
        "source_id": star.star_id,
        "relative_core_abs_residual": float(np.median(np.abs(err[local])) /
                                            max(cutout["peak_nMgy"], 1e-12)),
        "flux_nMgy_diagnostic": float(flux),
        "background_rms_nMgy": cutout["rms_nMgy"],
    }


def build_empirical_coadd_psf(
    frame: NativeBandFrame, *,
    target_icrs_deg: tuple[float, float],
    target_protection_arcsec: float,
    config: PSF0Config | None = None,
) -> PSF0Result:
    """Build an exploratory native-grid ePSF from independent local star cutouts.

    No automatic archive access, source substitution, external Gaia query,
    wide-wing extrapolation, or scientific adoption. Fails when the stellar
    sample is too small or insufficiently isolated.
    """
    if not isinstance(frame, NativeBandFrame):
        raise TypeError("Pass NativeBandFrame, never a remote/preview image")
    cfg = config if config is not None else PSF0Config()
    if not isinstance(cfg, PSF0Config):
        raise TypeError("Expected PSF0Config")
    try:
        import photutils
        import astropy
        from astropy.stats import sigma_clipped_stats
        from photutils.detection import DAOStarFinder
        from photutils.psf import EPSFBuilder, EPSFStar, EPSFStars
    except ImportError as exc:
        raise RuntimeError("Empirical PSF requires photutils and astropy; install -e '.[exact]'") from exc
    data = np.asarray(frame.image_nMgy)
    protected = _target_exclusion(frame, target_icrs_deg, target_protection_arcsec)
    valid = np.array(frame.valid, dtype=bool, copy=True)
    for key in (f"SATUR_{frame.band.upper()}", f"ALLMASK_{frame.band.upper()}"):
        if key in frame.mask_reasons:
            valid &= ~frame.mask_reasons[key]
    sample = valid & ~protected
    if sample.sum() < 400:
        raise ValueError("Insufficient background pixels outside protected Leo P/target region")
    _, sky0, std0 = sigma_clipped_stats(data[sample], sigma=3., maxiters=6)
    if not np.isfinite(std0) or std0 <= 0:
        raise ValueError("Cannot estimate an initial detection noise scale")
    detection_image = np.where(sample, data-sky0, 0.)
    finder = DAOStarFinder(
        threshold=cfg.detection_sigma*std0,
        fwhm=cfg.guessed_fwhm_pix,
        exclude_border=True)
    detections = finder(detection_image)
    if detections is None or len(detections) == 0:
        raise ValueError("No stellar candidates detected; PSF-0 cannot be built")
    # Split source identification and star-quality selection. All positive
    # detections are used for isolation even if a candidate is later rejected.
    # Photutils 2.x and 3.x use different finder centroid column names.
    # Resolve the documented schema explicitly instead of silently losing
    # candidates or assuming the installed release.
    names = set(detections.colnames)
    xkey = "x_centroid" if "x_centroid" in names else "xcentroid"
    ykey = "y_centroid" if "y_centroid" in names else "ycentroid"
    if not {xkey, ykey, "peak"}.issubset(names):
        raise ValueError("Unsupported photutils DAOStarFinder column schema")
    candidates = sorted(
        [(float(row[xkey]), float(row[ykey]), float(row["peak"]))
         for row in detections
         if np.isfinite(row[xkey]) and np.isfinite(row[ykey])
         and np.isfinite(row["peak"])],
        key=lambda item: (-item[2], item[1], item[0]))
    if len(candidates) < cfg.min_train_stars+cfg.min_holdout_stars:
        raise ValueError("Too few detected star candidates for an independent PSF holdout")
    positions = np.array([item[:2] for item in candidates])
    neighbors = cKDTree(positions)
    rejected = {}
    def reject(reason):
        rejected[reason] = rejected.get(reason, 0)+1
    accepted = []
    cutouts = {}
    max_candidates = max(cfg.max_training_stars*5, 100)
    # Sort by brightness for capping, then spatially sort/split below.
    for ix, (x, y, peak) in enumerate(candidates[:max_candidates]):
        if len(neighbors.query_ball_point([x,y], r=cfg.min_neighbor_separation_pix)) != 1:
            reject("NEIGHBOR_WITHIN_ISOLATION_RADIUS")
            continue
        data_cut, cause = _star_cutout(frame, x, y, cfg, protected)
        if cause:
            reject(cause)
            continue
        star_id = f"star-{ix+1:04d}"
        ra, dec = frame.wcs.all_pix2world(x, y, 0)
        if not np.isfinite(ra) or not np.isfinite(dec):
            reject("INVALID_STAR_WCS")
            continue
        accepted.append((star_id, x, y, float(ra)%360, float(dec),
                         peak, data_cut))
        cutouts[star_id] = data_cut
        if len(accepted) >= cfg.max_training_stars + max(
                cfg.min_holdout_stars, cfg.max_training_stars//(cfg.holdout_every-1)):
            break
    # Prevent the brightest section of one brick monopolizing the holdout.
    accepted.sort(key=lambda t: (int(t[2]//(cfg.cutout_size_pix*2)),
                                 int(t[1]//(cfg.cutout_size_pix*2)), t[1], t[0]))
    if len(accepted) < cfg.min_train_stars + cfg.min_holdout_stars:
        raise ValueError(
            f"Insufficient isolated, unsaturated, sky-screened stars: "
            f"{len(accepted)} (minimum {cfg.min_train_stars+cfg.min_holdout_stars})")
    # Deterministic holdout every fifth source after spatial ordering.
    holdout_idx = {i for i in range(len(accepted)) if i % cfg.holdout_every == 0}
    if len(holdout_idx) < cfg.min_holdout_stars:
        holdout_idx.update(range(cfg.min_holdout_stars))
    training_entries = [v for i,v in enumerate(accepted) if i not in holdout_idx]
    holdout_entries = [v for i,v in enumerate(accepted) if i in holdout_idx]
    if len(training_entries) < cfg.min_train_stars:
        raise ValueError("Independent holdout would leave too few training stars")
    training_entries = training_entries[:cfg.max_training_stars]
    def metadata(entry, role):
        ident, x, y, ra, dec, peak, c = entry
        return PSF0Star(ident, x, y, ra, dec, c["peak_nMgy"],
                        c["sky_nMgy"], c["rms_nMgy"],
                        c["gradient_sigma_span"], role)
    train = tuple(metadata(v, "train") for v in training_entries)
    holdout = tuple(metadata(v, "holdout") for v in holdout_entries)
    epsf_stars = []
    for entry in training_entries:
        c = entry[-1]
        epsf_stars.append(EPSFStar(
            c["image"], weights=c["weights"],
            cutout_center=c["center_in_cutout"],
            id_label=entry[0]))
    fitted = EPSFBuilder(
        oversampling=1, shape=(cfg.cutout_size_pix, cfg.cutout_size_pix),
        maxiters=cfg.max_iterations, progress_bar=False)(EPSFStars(epsf_stars))
    psf_model = fitted.epsf if hasattr(fitted, "epsf") else fitted[0]
    raw_kernel = np.asarray(psf_model.data, dtype=float)
    if raw_kernel.shape != (cfg.cutout_size_pix, cfg.cutout_size_pix):
        raise ValueError("Unexpected ePSF output grid/shape; refuse silent resampling")
    if not np.isfinite(raw_kernel).all():
        raise ValueError("Nonfinite ePSF solution")
    negative_fraction = float(np.sum(np.abs(raw_kernel[raw_kernel<0])) /
                              max(np.sum(np.abs(raw_kernel)), 1e-20))
    if negative_fraction > .025:
        raise ValueError("Strong negative EPSF ringing; inspect star selection/sky")
    kernel = np.maximum(raw_kernel, 0.)
    if kernel.sum() <= 0:
        raise ValueError("Nonpositive empirical PSF normalization")
    kernel /= kernel.sum()
    checks = []
    for entry in holdout_entries:
        checks.append(_holdout_diagnostic(metadata(entry, "holdout"), entry[-1], kernel))
    med_core_resid = float(np.median([x["relative_core_abs_residual"] for x in checks]))
    flags = [
        "CORE_PSF_ONLY_EXTENDED_WINGS_UNVALIDATED",
        "SPATIAL_REPRESENTATIVITY_NOT_VALIDATED",
        "CORRELATED_COADD_NOISE_UNMODELED",
        "PSF_STAR_IDENTITIES_NOT_GAIA_VERIFIED",
        "CUTOUT_SKY_MAY_ABSORB_WINGS",
        "FITTED_STAR_FLUX_NOT_TOTAL_FLUX",
        "HOLDOUT_RESIDUAL_NOT_CALIBRATED_CHI_SQUARE",
        "TARGET_PROTECTION_LOCAL_WCS_APPROXIMATION",
    ]
    if negative_fraction > 1e-4:
        flags.append("EPSF_NEGATIVE_RINGING_CLIPPED")
    if hasattr(fitted, "converged") and not fitted.converged:
        flags.append("EPSF_NOT_CONVERGED")
    source_hashes = _verified_source_hashes(frame)
    lineage = None
    if source_hashes is None:
        flags.append("SYNTHETIC_TEST_INPUT" if frame.source_evidence.get("synthetic") else
                     "SOURCE_HASHES_NOT_VERIFIED")
        source_hashes = {}
    else:
        digest = sha256(json.dumps(source_hashes, sort_keys=True,
                                   separators=(",", ":")).encode("utf8")).hexdigest()
        selection_id = sha256(json.dumps(
            [(s.star_id, s.x_pix, s.y_pix, s.role) for s in train+holdout],
            sort_keys=True).encode("utf8")).hexdigest()
        star_sky_id = sha256(json.dumps(
            [(s.star_id, s.cutout_sky_nMgy, s.cutout_rms_nMgy,
              s.gradient_sigma_span) for s in train+holdout],
            sort_keys=True).encode("utf8")).hexdigest()
        lineage = PSFLineageKey(
            observation_id=frame.dataset_id, instrument=(
                "BASS_90PRIME" if frame.region == "north" and frame.band in ("g","r")
                else "MZLS_MOSAIC3" if frame.region == "north" and frame.band == "z"
                else "DECAM"),
            band=frame.band, source_sha256=digest,
            star_selection_sha256=selection_id,
            cutout_sky_sha256=star_sky_id,
            recipe_sha256=cfg.digest, model_scope="global").digest
    instrument = ("BASS_90PRIME" if frame.region == "north" and frame.band in ("g","r")
                  else "MZLS_MOSAIC3" if frame.region == "north" and frame.band == "z"
                  else "DECAM")
    # FWHM here is only the guessed finder width, NOT measured EPSF FWHM.
    metrics = {
        "holdout_median_abs_core_residual_over_peak": med_core_resid,
        "holdout_diagnostics": checks,
        "epsf_negative_mass_fraction_before_clipping": negative_fraction,
        "initial_detection_background_nMgy_per_pixel": float(sky0),
        "initial_detection_rms_nMgy_per_pixel": float(std0),
        "assumed_detection_fwhm_pix_not_measured": cfg.guessed_fwhm_pix,
        "epsf_converged_reported": (bool(fitted.converged)
                                    if hasattr(fitted, "converged") else None),
    }
    return PSF0Result(
        dataset_id=frame.dataset_id, instrument=instrument, kernel=kernel,
        training=train, holdout=holdout,
        rejected_reasons=rejected, detection_count=len(candidates),
        metrics=metrics, flags=tuple(flags), recipe=cfg,
        source_hashes=source_hashes, lineage_id=lineage,
        center_target_icrs=(float(target_icrs_deg[0]),float(target_icrs_deg[1])),
        target_protection_arcsec=float(target_protection_arcsec),
        algorithm_versions={
            "recipe": RECIPE, "numpy":np.__version__,
            "scipy":__import__("scipy").__version__,
            "photutils":photutils.__version__,
            "astropy":astropy.__version__,
        })
