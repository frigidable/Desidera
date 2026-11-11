"""Robust *residual* background estimators for already sky-subtracted coadds.

Fit constant or tilted-plane candidates to explicitly labelled sky-training
pixels; evaluate separate holdouts. No automatic adoption, masked-pixel
reconstruction, uncertainty calibration, or mutation of source images.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Literal
import numpy as np

from .scene import Operation, SceneState, evaluate_mask


class FitStatus(str, Enum):
    FITTED = 'FITTED'
    INSUFFICIENT_SKY = 'INSUFFICIENT_SKY'
    ILL_CONDITIONED = 'ILL_CONDITIONED'
    NUMERICAL_FAILURE = 'NUMERICAL_FAILURE'


class HoldoutStatus(str, Enum):
    NOT_EVALUATED = 'NOT_EVALUATED'
    INSUFFICIENT_HOLDOUT = 'INSUFFICIENT_HOLDOUT'
    INSUFFICIENT_SPATIAL_COVERAGE = 'INSUFFICIENT_SPATIAL_COVERAGE'
    MEASURED_UNCALIBRATED = 'MEASURED_UNCALIBRATED'
    REJECTED_BY_USER_THRESHOLD = 'REJECTED_BY_USER_THRESHOLD'
    WITHIN_USER_THRESHOLD = 'WITHIN_USER_THRESHOLD'


@dataclass(frozen=True)
class BackgroundSettings:
    min_training_pixels: int = 12
    min_holdout_pixels: int = 12
    min_holdout_quadrants: int = 2
    max_design_condition: float = 1e8
    huber_t: float = 1.345
    max_iterations: int = 50
    convergence_relative: float = 1e-10
    max_abs_holdout_tile_median_nMgy: float | None = None

    def __post_init__(self):
        if self.min_training_pixels < 3 or self.min_holdout_pixels < 1:
            raise ValueError('Pixel counts must be positive (training at least 3)')
        if not 1 <= self.min_holdout_quadrants <= 4:
            raise ValueError('min_holdout_quadrants must be between 1 and 4')
        for value in (self.max_design_condition, self.huber_t, self.convergence_relative):
            if not np.isfinite(value) or value <= 0:
                raise ValueError('Numerical settings must be finite and positive')
        if self.max_iterations < 1:
            raise ValueError('max_iterations must be >= 1')
        limit = self.max_abs_holdout_tile_median_nMgy
        if limit is not None and (not np.isfinite(limit) or limit < 0):
            raise ValueError('User-supplied holdout threshold must be non-negative finite')


@dataclass(frozen=True)
class HoldoutTile:
    quadrant: int
    pixels: int
    median_residual_nMgy: float
    mean_residual_nMgy: float
    rms_residual_nMgy: float


@dataclass(frozen=True)
class BackgroundCandidate:
    method: Literal['constant', 'tilted_plane']
    fit_status: FitStatus
    holdout_status: HoldoutStatus
    reference_xy_pix: tuple[float, float]
    coefficients_nMgy: tuple[float, float, float] | None
    training_pixels: int
    holdout_pixels: int
    design_condition: float | None
    iterations: int
    training_median_residual_nMgy: float | None
    training_rms_residual_nMgy: float | None
    holdout_median_residual_nMgy: float | None
    holdout_rms_residual_nMgy: float | None
    tiles: tuple[HoldoutTile, ...]
    warnings: tuple[str, ...] = ()
    grid_id: str | None = None
    band: str | None = None
    image_shape: tuple[int, int] | None = None

    def evaluate_at(self, x_pix, y_pix):
        """Evaluate a residual-sky candidate in nanomaggies per pixel."""
        if self.coefficients_nMgy is None or self.fit_status != FitStatus.FITTED:
            raise ValueError('Cannot evaluate an unsuccessful background candidate')
        b0, bx, by = self.coefficients_nMgy
        x, y = np.asarray(x_pix, dtype=float), np.asarray(y_pix, dtype=float)
        return b0 + bx * (x - self.reference_xy_pix[0]) + by * (y - self.reference_xy_pix[1])

    def image(self, shape: tuple[int, int]) -> np.ndarray:
        """Materialize a fresh model map without modifying the input."""
        if self.image_shape is not None and tuple(shape) != tuple(self.image_shape):
            raise ValueError('Candidate source image shape mismatch')
        yy, xx = np.indices(shape, dtype=float)
        return np.asarray(self.evaluate_at(xx, yy), dtype=float)


def _weighted_median(values, weights):
    order = np.argsort(values)
    v, w = values[order], weights[order]
    return float(v[np.searchsorted(np.cumsum(w), 0.5 * np.sum(w), side='left')])


def _robust_scale(residuals):
    mad = 1.4826 * np.median(np.abs(residuals - np.median(residuals)))
    # Scale must be data-derived. A formal pixel error can exceed scatter by
    # orders of magnitude and thereby defeat Huber downweighting.
    floor = 64. * np.finfo(float).eps * max(1., float(np.max(np.abs(residuals))))
    return max(float(mad), floor)


def _weighted_lstsq(design, values, weights, max_condition):
    sqrtw = np.sqrt(weights / np.median(weights))
    matrix = design * sqrtw[:, None]
    rhs = values * sqrtw
    singular = np.linalg.svd(matrix, compute_uv=False)
    condition = float(singular[0] / singular[-1]) if singular[-1] > 0 else np.inf
    if not np.isfinite(condition) or condition > max_condition:
        return None, condition
    coefficients, _, rank, _ = np.linalg.lstsq(matrix, rhs, rcond=None)
    if rank != design.shape[1] or not np.all(np.isfinite(coefficients)):
        return None, condition
    return coefficients, condition


def _fit_robust(design, values, ivar, method, settings):
    if method == 'constant':
        coefficients = np.array([_weighted_median(values, ivar)])
        condition = 1.0
    else:
        coefficients, condition = _weighted_lstsq(
            design, values, ivar, settings.max_design_condition)
        if coefficients is None:
            return None, condition, 0
    for iteration in range(1, settings.max_iterations + 1):
        residuals = values - design @ coefficients
        centered = residuals / _robust_scale(residuals)
        huber_weights = np.ones_like(centered)
        outliers = np.abs(centered) > settings.huber_t
        huber_weights[outliers] = settings.huber_t / np.abs(centered[outliers])
        updated, condition = _weighted_lstsq(
            design, values, ivar * huber_weights, settings.max_design_condition)
        if updated is None:
            return None, condition, iteration
        converged = (np.max(np.abs(updated - coefficients)) <=
                     settings.convergence_relative *
                     max(1., float(np.max(np.abs(coefficients)))))
        coefficients = updated
        if converged:
            break
    return coefficients, condition, iteration


def _rms(values):
    return float(np.sqrt(np.mean(values**2)))


def fit_residual_background(
    image_nMgy: np.ndarray,
    inverse_variance_nMgy2: np.ndarray,
    scene: SceneState,
    *,
    method: Literal['constant', 'tilted_plane'],
    grid_id: str,
    observation_valid: np.ndarray | None = None,
    settings: BackgroundSettings | None = None,
) -> BackgroundCandidate:
    """Estimate one sky candidate, retaining fitting and QA status separately.

    Fitting samples are selected through the SKY_FIT operation; holdout samples
    through SKY_HOLDOUT. Their spatial separation beyond disjoint membership
    is the researcher's responsibility, and covariance is not modelled.
    """
    if method not in ('constant', 'tilted_plane'):
        raise ValueError('Only constant and tilted_plane candidates are implemented')
    settings = settings or BackgroundSettings()
    image = np.asarray(image_nMgy, dtype=float)
    ivar = np.asarray(inverse_variance_nMgy2, dtype=float)
    train = evaluate_mask(scene, image, ivar, grid_id=grid_id,
                          operation=Operation.SKY_FIT,
                          observation_valid=observation_valid)
    held = evaluate_mask(scene, image, ivar, grid_id=grid_id,
                         operation=Operation.SKY_HOLDOUT,
                         observation_valid=observation_valid)
    ntrain, nheld = train.n_valid, held.n_valid
    ny, nx = scene.shape
    xref, yref = (nx-1)/2., (ny-1)/2.
    ref = (xref, yref)

    def failure(status, cond=None, iterations=0):
        return BackgroundCandidate(method, status, HoldoutStatus.NOT_EVALUATED,
                                   ref, None, ntrain, nheld, cond, iterations,
                                   None, None, None, None, (), (status.value,),
                                   grid_id=grid_id, band=scene.band, image_shape=scene.shape)

    if ntrain < settings.min_training_pixels:
        return failure(FitStatus.INSUFFICIENT_SKY)
    yy, xx = np.nonzero(train.valid)
    values, weights = image[yy, xx], ivar[yy, xx]
    xs, ys = max(1., xref), max(1., yref)
    if method == 'constant':
        design = np.ones((ntrain, 1))
    else:
        design = np.column_stack(
            (np.ones(ntrain), (xx-xref)/xs, (yy-yref)/ys))
    try:
        fit, cond, iterations = _fit_robust(
            design, values, weights, method, settings)
    except (ValueError, np.linalg.LinAlgError, FloatingPointError, OverflowError):
        return failure(FitStatus.NUMERICAL_FAILURE)
    if fit is None:
        return failure(FitStatus.ILL_CONDITIONED, cond, iterations)

    coeff = ((float(fit[0]), 0., 0.) if method == 'constant'
             else (float(fit[0]), float(fit[1])/xs, float(fit[2])/ys))
    provisional = BackgroundCandidate(
        method, FitStatus.FITTED, HoldoutStatus.NOT_EVALUATED, ref, coeff,
        ntrain, nheld, cond, iterations, None, None, None, None, (),
        grid_id=grid_id, band=scene.band, image_shape=scene.shape)
    residual_train = values - provisional.evaluate_at(xx, yy)
    hy, hx = np.nonzero(held.valid)
    residual_held = image[hy, hx] - provisional.evaluate_at(hx, hy)
    tiles = []
    for quadrant in range(4):
        selection = ((hx > xref).astype(int) +
                     2 * (hy > yref).astype(int)) == quadrant
        r = residual_held[selection]
        if r.size:
            tiles.append(HoldoutTile(
                quadrant, int(r.size), float(np.median(r)),
                float(np.mean(r)), _rms(r)))
    warnings = []
    if nheld < settings.min_holdout_pixels:
        qa = HoldoutStatus.INSUFFICIENT_HOLDOUT
        warnings.append('HOLDOUT_PIXELS_BELOW_MINIMUM')
    elif len(tiles) < settings.min_holdout_quadrants:
        qa = HoldoutStatus.INSUFFICIENT_SPATIAL_COVERAGE
        warnings.append('HOLDOUT_SPATIAL_COVERAGE_LIMITED')
    elif settings.max_abs_holdout_tile_median_nMgy is None:
        qa = HoldoutStatus.MEASURED_UNCALIBRATED
        warnings.append('NO_CALIBRATED_HOLDOUT_THRESHOLD')
    elif max(abs(tile.median_residual_nMgy) for tile in tiles) > settings.max_abs_holdout_tile_median_nMgy:
        qa = HoldoutStatus.REJECTED_BY_USER_THRESHOLD
        warnings.append('HOLDOUT_RESIDUAL_EXCEEDS_USER_LIMIT')
    else:
        qa = HoldoutStatus.WITHIN_USER_THRESHOLD

    return BackgroundCandidate(
        method, FitStatus.FITTED, qa, ref, coeff, ntrain, nheld,
        cond, iterations, float(np.median(residual_train)),
        _rms(residual_train),
        float(np.median(residual_held)) if nheld else None,
        _rms(residual_held) if nheld else None,
        tuple(tiles), tuple(warnings),
        grid_id=grid_id, band=scene.band, image_shape=scene.shape)


def fit_observation_background(observation, scene: SceneState, *, method,
                               settings: BackgroundSettings | None = None) -> BackgroundCandidate:
    """Fit from v0.2 ObservationSet, preserving its existing validity map."""
    from .mask_engine import MaskEngine
    MaskEngine(scene).evaluate(observation, Operation.SKY_FIT)
    return fit_residual_background(
        observation.image_nMgy, observation.inverse_variance_nMgy2,
        scene, method=method, grid_id=scene.grid_id,
        observation_valid=observation.valid, settings=settings)



def _automatic_scene_sky_pixels(
        dataset, preview, contamination, *,
        tile_edge_native_pix=None, grid_offset_xy_native_pix=(0, 0),
        target_guard_factor=1.0, source_guard_factor=1.0):
    """Deterministic spatially disjoint source-vetoed tiles, no image-based picking.

    Selection is solely based on geometry, validity, existing Scene
    provenance and fixed design constants. We do NOT choose tiles after
    comparing their residual intensities, or assume they are galaxy-free.
    """
    from scipy.ndimage import distance_transform_edt
    from astropy.wcs.utils import proj_plane_pixel_scales

    shape = tuple(dataset["pixel_shape"])
    ny, nx = shape
    if (not np.isfinite(target_guard_factor) or target_guard_factor < 1.
            or not np.isfinite(source_guard_factor) or source_guard_factor < 1.):
        raise ValueError("Sensitivity guards may only extend original exclusions")
    if (len(grid_offset_xy_native_pix) != 2 or
            any(not isinstance(k, (int, np.integer)) or k < 0
                for k in grid_offset_xy_native_pix)):
        raise ValueError("Fixed sky tile-grid offset must be nonnegative integers")
    if (tile_edge_native_pix is not None and
            (not isinstance(tile_edge_native_pix, (int, np.integer))
             or not 16 <= tile_edge_native_pix <= 256)):
        raise ValueError("Invalid predeclared sky tile edge")
    hard = np.asarray(dataset["native_hard"], dtype=bool)
    image = np.asarray(dataset["data"]["image"])
    invvar = np.asarray(dataset["data"]["invvar"])
    auto = np.asarray(contamination["auto_mask"], dtype=bool)
    manual = np.asarray(preview["manual_exclude"], dtype=bool)
    target = (np.asarray(contamination["protected"], dtype=bool)
              | preview["target_protect"] | preview["isophote_candidate"])
    source = auto | manual
    scales = np.asarray(proj_plane_pixel_scales(
        dataset["wcs"].celestial), dtype=float)
    if scales.shape != (2,) or not np.all(np.isfinite(scales)) or np.any(scales <= 0):
        raise ValueError("Cannot determine Scene native pixel scale")
    arcsec_per_pixel = float(np.sqrt(scales[0]*scales[1])*3600.)
    # This is an exclusion-only heuristic, not validated galaxy extent.
    equivalent_target_radius = np.sqrt(np.count_nonzero(target)/np.pi)
    target_guard_px = int(np.ceil(target_guard_factor*max(
        15./arcsec_per_pixel, .30*equivalent_target_radius)))
    source_guard_px = int(np.ceil(source_guard_factor*max(
        2., 6./arcsec_per_pixel)))

    def guarded(mask, radius):
        if not np.any(mask):
            return np.zeros(shape, dtype=bool)
        return distance_transform_edt(~mask) <= radius

    target_veto = guarded(target, target_guard_px)
    catalog_veto = guarded(source, source_guard_px)

    # Existing Tractor/Gaia sources only. These radii are intentionally
    # conservative pilot vetoes and NOT a measured bright-star halo model.
    inventory = contamination.get("source_inventory")
    entries = inventory.get("entries", ()) if isinstance(inventory, dict) else ()
    bright_veto = np.zeros(shape, dtype=bool)
    bright_n = 0
    for entry in entries:
        mag = entry.get("gaia_g_mag")
        if mag is None:
            continue
        try:
            mag = float(mag)
            sx, sy = float(entry["x"]), float(entry["y"])
        except (TypeError, ValueError, KeyError, OverflowError):
            continue
        if not (np.isfinite(mag) and np.isfinite(sx) and np.isfinite(sy)
                and 0. < mag <= 18.):
            continue
        guard_arcsec = (60. if mag <= 13. else
                        40. if mag <= 15. else
                        30. if mag <= 16.5 else 18.)
        radius = guard_arcsec/arcsec_per_pixel
        x0, x1 = max(0, int(np.floor(sx-radius))), min(
            nx, int(np.ceil(sx+radius))+1)
        y0, y1 = max(0, int(np.floor(sy-radius))), min(
            ny, int(np.ceil(sy+radius))+1)
        if x1 <= x0 or y1 <= y0:
            continue
        xx, yy = np.meshgrid(np.arange(x0,x1), np.arange(y0,y1))
        bright_veto[y0:y1,x0:x1] |= (xx-sx)**2+(yy-sy)**2 <= radius**2
        bright_n += 1

    valid = ((~hard) & np.isfinite(image) &
             np.isfinite(invvar) & (invvar > 0))
    safe = valid & ~target_veto & ~catalog_veto & ~bright_veto
    # Full spatial tiles, with a 2-native-pixel inner edge guard.
    # The same physical sky pixels are never assigned to train and holdout.
    tile = (max(16, min(96, int(min(shape)/12)))
            if tile_edge_native_pix is None else int(tile_edge_native_pix))
    offset_x, offset_y = map(int, grid_offset_xy_native_pix)
    if offset_x >= tile or offset_y >= tile:
        raise ValueError("Grid offset must be smaller than tile edge")
    border = 2
    train = np.zeros(shape, dtype=bool)
    holdout = np.zeros(shape, dtype=bool)
    tiles_train, tiles_holdout = [], []
    n_candidate = 0
    n_eligible = 0
    # The fixed hash is wholly geometry-based, decided before pixel fluxes
    # are considered. No resampling or cherry-picking after a fit.
    for grid_y0 in range(-offset_y, ny, tile):
        for grid_x0 in range(-offset_x, nx, tile):
            y0, y1 = max(0, grid_y0), min(ny, grid_y0+tile)
            x0, x1 = max(0, grid_x0), min(nx, grid_x0+tile)
            if min(y1-y0, x1-x0) < .75*tile:
                continue
            n_candidate += 1
            yy0, yy1 = y0+border, y1-border
            xx0, xx1 = x0+border, x1-border
            footprint = safe[yy0:yy1, xx0:xx1]
            count = int(np.count_nonzero(footprint))
            total = int(footprint.size)
            if total < 25 or count < 40 or count / total < .65:
                continue
            n_eligible += 1
            yy, xx = np.nonzero(footprint)
            # Cap the nominal pixel weighting per tile. This is a
            # deterministic spatial subsample, not sigma clipping.
            if len(xx) > 384:
                index = np.linspace(0, len(xx)-1, 384, dtype=np.int64)
                yy, xx = yy[index], xx[index]
            ty, tx = ((grid_y0+offset_y)//tile,
                      (grid_x0+offset_x)//tile)
            is_holdout = (tx+2*ty) % 5 == 0
            selected = holdout if is_holdout else train
            selected[yy0+yy,xx0+xx] = True
            record = {
                "tile_index": [int(tx),int(ty)],
                "tile_bounds_native_xyxy": [int(x0),int(y0),int(x1),int(y1)],
                "n_available_valid_pixels": count,
                "n_sampled_pixels": int(len(xx)),
                "eligible_fraction_of_interior": float(count/total),
                "quadrant": int((x0+x1)/2 > (nx-1)/2) +
                            2*int((y0+y1)/2 > (ny-1)/2),
            }
            (tiles_holdout if is_holdout else tiles_train).append(record)

    def coverage(rows):
        return sorted({r["quadrant"] for r in rows})

    enough = (len(tiles_train) >= 8 and len(tiles_holdout) >= 4
              and len(coverage(tiles_train)) >= 3
              and len(coverage(tiles_holdout)) >= 3
              and int(train.sum()) >= 100 and int(holdout.sum()) >= 100)
    qa = {
        "strategy": "SCENE_MASKS_FIXED_SPATIAL_TILES_V1",
        "selection_predeclared_geometry_only_not_flux_optimized": True,
        "sky_regions_drawn_by_user_required": False,
        "manual_sky_annotations_ignored_not_deleted": True,
        "native_pixel_scale_arcsec": arcsec_per_pixel,
        "tile_edge_native_pix": tile,
        "tile_grid_offset_xy_native_pix": [offset_x, offset_y],
        "target_guard_factor": float(target_guard_factor),
        "source_guard_factor": float(source_guard_factor),
        "tile_inner_border_native_pix": border,
        "min_eligible_fraction_per_tile": .65,
        "max_sampled_pixels_per_tile": 384,
        "target_outskirts_guard_native_pix": target_guard_px,
        "source_neighbour_guard_native_pix": source_guard_px,
        "bright_gaia_veto_is_unvalidated_radius_heuristic": True,
        "n_existing_bright_gaia_sources_vetoed": bright_n,
        "bright_gaia_inventory_present": isinstance(inventory, dict),
        "candidate_tiles": n_candidate,
        "eligible_tiles": n_eligible,
        "training_tiles": len(tiles_train),
        "holdout_tiles": len(tiles_holdout),
        "training_quadrants": coverage(tiles_train),
        "holdout_quadrants": coverage(tiles_holdout),
        "training_tile_evidence": tiles_train,
        "holdout_tile_evidence": tiles_holdout,
        "native_valid_pixels": int(np.count_nonzero(valid)),
        "geometry_eligible_pixels_before_tile_sampling":
            int(np.count_nonzero(safe)),
        "target_guard_veto_pixels": int(np.count_nonzero(target_veto)),
        "catalogue_guard_veto_pixels": int(np.count_nonzero(catalog_veto)),
        "bright_star_veto_pixels": int(np.count_nonzero(bright_veto)),
        "training_pixels": int(np.count_nonzero(train)),
        "holdout_pixels": int(np.count_nonzero(holdout)),
        "geometry_sufficient_for_diagnostic": bool(enough),
        "galaxy_free_certified": False,
        "psf_wings_calibrated": False,
        "independent_exposure_holdout": False,
    }
    assert not np.any(train & holdout)
    assert not np.any((train | holdout) & ~safe)
    return train, holdout, qa



def _audit_automatic_sky_tile_residuals(image, sky_train, sky_holdout,
                                        selection, candidates):
    """Equal-area-tile descriptive checks without treating pixels as iid.

    Post-fit statistics are never used to reselect training/holdout tiles,
    tune veto radii or replace the predeclared model choice.
    """
    from collections import Counter

    edge = selection["tile_edge_native_pix"]
    ny, nx = image.shape
    checks = {}
    warnings = []
    for role, records in (
            ("training", selection["training_tile_evidence"]),
            ("holdout", selection["holdout_tile_evidence"])):
        row_counts = Counter(r["tile_index"][1] for r in records)
        dominant = max(row_counts.values(), default=0)
        fraction = (float(dominant / len(records))
                    if records else None)
        missing_quadrants = sorted(set(range(4)) -
                                   set(selection[f"{role}_quadrants"]))
        clipped_edge_tiles = sum(
            (r["tile_index"][1]+1)*edge > ny or
            (r["tile_index"][0]+1)*edge > nx
            for r in records)
        checks[role] = {
            "n_tiles": len(records),
            "n_distinct_tile_rows": len(row_counts),
            "max_tiles_in_any_single_row": dominant,
            "dominant_single_row_fraction": fraction,
            "n_partial_edge_tiles": clipped_edge_tiles,
            "partial_edge_tile_fraction": (
                float(clipped_edge_tiles/len(records))
                if records else None),
            "missing_quadrants": missing_quadrants,
            "spatial_mask_selection_unchanged": True,
        }
        if fraction is not None and fraction > .40:
            warnings.append(f"{role.upper()}_TILES_CLUSTERED_IN_SINGLE_ROW")
        if missing_quadrants:
            warnings.append(f"{role.upper()}_QUADRANTS_NOT_ALL_REPRESENTED")
        if clipped_edge_tiles >= max(2, int(np.ceil(.40*len(records)))):
            warnings.append(f"{role.upper()}_TILES_DOMINATED_BY_PARTIAL_IMAGE_EDGE")

    method_details = {}
    for method, candidate in candidates.items():
        if candidate.fit_status != FitStatus.FITTED:
            continue
        data = {}
        for role, records, valid in (
                ("training", selection["training_tile_evidence"], sky_train),
                ("holdout", selection["holdout_tile_evidence"], sky_holdout)):
            tiles = []
            medians = []
            for record in records:
                ix, iy = record["tile_index"]
                x0, x1 = ix*edge, min(nx, (ix+1)*edge)
                y0, y1 = iy*edge, min(ny, (iy+1)*edge)
                yy_local, xx_local = np.nonzero(valid[y0:y1, x0:x1])
                xx = xx_local+x0
                yy = yy_local+y0
                if len(xx) != record["n_sampled_pixels"]:
                    raise ValueError("Sky sample tile geometry/provenance changed")
                residual = np.asarray(
                    image[yy,xx] - candidate.evaluate_at(xx,yy),
                    dtype=float)
                median = float(np.median(residual))
                medians.append(median)
                tiles.append({
                    "tile_index": list(record["tile_index"]),
                    "quadrant": int(record["quadrant"]),
                    "n_samples": int(len(xx)),
                    "median_residual_nMgy_per_pixel": median,
                    "rms_residual_nMgy_per_pixel":
                        float(np.sqrt(np.mean(residual**2))),
                })
            data[role] = {
                "tiles": tiles,
                "equal_tile_median_of_medians_nMgy_per_pixel": (
                    float(np.median(medians)) if medians else None),
                "tile_median_peak_to_peak_nMgy_per_pixel": (
                    float(np.ptp(medians)) if medians else None),
                "tile_median_mad_nMgy_per_pixel": (
                    float(np.median(np.abs(np.asarray(medians) -
                                            np.median(medians))))
                    if medians else None),
                "equal_tile_median_abs_residual_nMgy_per_pixel": (
                    float(np.median(np.abs(medians))) if medians else None),
            }
        train_rms = candidate.training_rms_residual_nMgy
        hold_rms = candidate.holdout_rms_residual_nMgy
        ratio = (float(train_rms/hold_rms)
                 if (train_rms is not None and hold_rms is not None
                     and hold_rms > 0.) else None)
        data["training_to_holdout_pixel_rms_ratio"] = ratio
        data["tile_residuals_are_correlated_same_coadd_not_errors"] = True
        data["pixel_weighted_rms_does_not_validate_sky_level"] = True
        if ratio is not None and ratio > 2.:
            warnings.append(f"{method.upper()}_TRAINING_RMS_OVER_TWICE_HOLDOUT")
        method_details[method] = data
    return {
        "schema": "desidera.sky.automatic-tile-spatial-qa/1",
        "geometry": checks,
        "methods": method_details,
        "warning_flags": sorted(set(warnings)),
        "spatial_support_complete_for_science": False,
        "science_holds_even_without_warning": True,
        "no_tile_reselection_after_inspecting_residuals": True,
    }



def _audit_zero_sky_baseline(image, training, holdout, selection):
    """No additional sky subtraction: a mandatory DR11 reference candidate.

    Uses exact same fixed sky tiles, before any fitted model is adopted.
    Pixel-RMS comparisons are descriptive because coadd pixels correlate.
    """
    shape = image.shape
    ny, nx = shape
    tile_edge = int(selection["tile_edge_native_pix"])
    records = []
    for record in selection["holdout_tile_evidence"]:
        ix, iy = record["tile_index"]
        x0, x1 = ix*tile_edge, min(nx,(ix+1)*tile_edge)
        y0, y1 = iy*tile_edge, min(ny,(iy+1)*tile_edge)
        values = np.asarray(
            image[y0:y1,x0:x1][holdout[y0:y1,x0:x1]],
            dtype=float)
        if values.size != record["n_sampled_pixels"]:
            raise ValueError("Zero-background holdout sample lineage differs")
        records.append({
            "tile_index": list(record["tile_index"]),
            "quadrant": int(record["quadrant"]),
            "n_samples": int(values.size),
            "median_nMgy_per_pixel": float(np.median(values)),
            "rms_nMgy_per_pixel": float(np.sqrt(np.mean(values**2))),
        })
    vals = np.asarray(image[holdout],dtype=float)
    train_vals = np.asarray(image[training],dtype=float)
    if not vals.size or not train_vals.size:
        raise ValueError("Zero baseline requires same nonempty fixed sky samples")
    medians = np.array([r["median_nMgy_per_pixel"] for r in records])
    q = {}
    for k in range(4):
        ys, xs = np.nonzero(holdout)
        chosen = (xs > (nx-1)/2).astype(int) + (
            2*(ys > (ny-1)/2).astype(int)) == k
        if np.any(chosen):
            q[str(k)] = float(np.median(image[ys[chosen],xs[chosen]]))
    return {
        "method": "zero_additional_correction_DR11_coadd",
        "same_original_holdout_as_fitted_candidates": True,
        "background_subtracted_from_image": False,
        "additional_correction_nMgy_per_pixel": 0.0,
        "training_rms_nMgy_per_pixel":
            float(np.sqrt(np.mean(train_vals**2))),
        "holdout_rms_nMgy_per_pixel":
            float(np.sqrt(np.mean(vals**2))),
        "holdout_median_nMgy_per_pixel": float(np.median(vals)),
        "n_training_samples": int(train_vals.size),
        "n_holdout_samples": int(vals.size),
        "holdout_tile_median_peak_to_peak_nMgy_per_pixel":
            float(np.ptp(medians)),
        "holdout_tile_median_of_medians_nMgy_per_pixel":
            float(np.median(medians)),
        "holdout_tile_medians_nMgy_per_pixel": records,
        "holdout_quadrant_medians_nMgy_per_pixel": q,
        "validates_original_sky_level": False,
        "original_DR11_may_already_be_oversubtracted": True,
        "not_an_independent_noise_realization": True,
    }




def _make_automatic_sky_scene(
        *, project_id, grid, band, shape, revision,
        hard, source_mask, manual_mask, target_mask,
        sky_train, sky_holdout):
    """Construct diagnostic masks; never change supplied native arrays."""
    from lsbpipe.scene import LayerKind, MaskLayer, SceneState

    def layer(name, kind, values):
        return MaskLayer(
            name, kind, values, grid, band, "auto-sky-fixed-evidence")

    return SceneState(
        target_id=project_id, grid_id=grid, band=band, shape=shape,
        layers=(
            layer("native_hard", LayerKind.INVALID_DATA, hard),
            layer("catalogue_provisional", LayerKind.NEIGHBOR, source_mask),
            layer("manual_exclusion", LayerKind.MANUAL_EXCLUDE, manual_mask),
            layer("target_provisional", LayerKind.TARGET_PROTECT, target_mask),
            layer("sky_train_automatic_tiles", LayerKind.SKY_TRAIN, sky_train),
            layer("sky_holdout_automatic_tiles", LayerKind.SKY_HOLDOUT, sky_holdout),
        ), revision=revision,
    )


# These sensitivity trials are predeclared, not chosen from the best result.
# The baseline is the default 96px grid on Leo P; none can reduce its
# original source/galaxy protection. No new bricks, downloads or files.
_SKY_FIXED_SENSITIVITY_TRIALS = (
    ("tile_72_origin", {"tile_edge_native_pix": 72}),
    ("tile_96_shift48", {
        "tile_edge_native_pix": 96,
        "grid_offset_xy_native_pix": (48,48)}),
    ("tile_128_origin", {"tile_edge_native_pix": 128}),
    ("expanded_target_1p5", {"target_guard_factor": 1.5}),
    ("expanded_sources_1p5", {"source_guard_factor": 1.5}),
)


def _audit_current_project_sky_grid_sensitivity(
        dataset, preview, contamination, *, grid, band, shape,
        revision, original_hard, source_mask, manual_mask, target_mask,
        original_image, inverse_variance, anchor_native_xy,
        original_candidates, settings):
    """Predeclared same-project sky-grid perturbations, all non-adopted."""
    result = []
    ax, ay = map(float, anchor_native_xy)
    for name, overrides in _SKY_FIXED_SENSITIVITY_TRIALS:
        train, hold, selection = _automatic_scene_sky_pixels(
            dataset, preview, contamination, **overrides)
        item = {
            "name": name,
            "parameters": {
                "tile_edge_native_pix": selection["tile_edge_native_pix"],
                "tile_grid_offset_xy_native_pix": selection[
                    "tile_grid_offset_xy_native_pix"],
                "target_guard_factor": selection["target_guard_factor"],
                "source_guard_factor": selection["source_guard_factor"],
            },
            "training_tiles": selection["training_tiles"],
            "holdout_tiles": selection["holdout_tiles"],
            "training_quadrants": selection["training_quadrants"],
            "holdout_quadrants": selection["holdout_quadrants"],
            "training_pixels": selection["training_pixels"],
            "holdout_pixels": selection["holdout_pixels"],
            "geometry_sufficient": selection[
                "geometry_sufficient_for_diagnostic"],
            "status": "NOT_FITTED",
            "models": {},
            "no_reselection_based_on_pixel_brightness": True,
            "sky_areas_are_same_project_not_independent_data": True,
        }
        centers = []
        edge = selection["tile_edge_native_pix"]
        for record in selection["training_tile_evidence"]:
            x0,y0,x1,y1 = record["tile_bounds_native_xyxy"]
            centers.append(((x0+x1)/2.,(y0+y1)/2.))
        if centers:
            center_array = np.asarray(centers,dtype=float)
            item["target_inside_training_tile_center_bbox"] = bool(
                center_array[:,0].min() <= ax <= center_array[:,0].max()
                and center_array[:,1].min() <= ay <= center_array[:,1].max())
            item["nearest_training_tile_center_distance_native_pix"] = float(
                np.min(np.hypot(
                    center_array[:,0]-ax,center_array[:,1]-ay)))
        else:
            item["target_inside_training_tile_center_bbox"] = False
            item["nearest_training_tile_center_distance_native_pix"] = None
        if not selection["geometry_sufficient_for_diagnostic"]:
            item["status"] = "SKIPPED_INSUFFICIENT_FIXED_GEOMETRY"
            result.append(item)
            continue
        if not selection["bright_gaia_inventory_present"]:
            item["status"] = "SKIPPED_MISSING_SOURCE_INVENTORY"
            result.append(item)
            continue
        trial_scene = _make_automatic_sky_scene(
            project_id=dataset["project_id"], grid=grid, band=band,
            shape=shape, revision=revision, hard=original_hard,
            source_mask=source_mask, manual_mask=manual_mask,
            target_mask=target_mask, sky_train=train, sky_holdout=hold)
        for method in ("constant","tilted_plane"):
            fit = fit_residual_background(
                original_image, inverse_variance, trial_scene,
                method=method, grid_id=grid,
                observation_valid=~original_hard, settings=settings)
            if (fit.fit_status == FitStatus.FITTED
                    and fit.holdout_status == HoldoutStatus.MEASURED_UNCALIBRATED):
                baseline = original_candidates[method]
                item["models"][method] = {
                    "status": "FITTED_HOLDOUT_UNCALIBRATED",
                    "background_at_target_nMgy_per_pixel":
                        float(fit.evaluate_at(ax,ay)),
                    "change_vs_original_grid_at_target_nMgy_per_pixel":
                        float(fit.evaluate_at(ax,ay) -
                              baseline.evaluate_at(ax,ay)),
                    "holdout_pixel_rms_nMgy_per_pixel":
                        float(fit.holdout_rms_residual_nMgy),
                    "within_trial_zero_correction_rms_nMgy_per_pixel":
                        float(np.sqrt(np.mean(
                            np.asarray(original_image)[hold]**2))),
                }
            else:
                item["models"][method] = {
                    "status": "FIT_OR_HOLDOUT_REJECTED",
                    "fit_status": fit.fit_status.value,
                    "holdout_status": fit.holdout_status.value,
                }
        item["status"] = (
            "BOTH_MODELS_MEASURED_UNCALIBRATED"
            if all(v["status"] == "FITTED_HOLDOUT_UNCALIBRATED"
                   for v in item["models"].values())
            else "INCOMPLETE_TRIAL_NOT_ADOPTED")
        result.append(item)

    summary = {}
    for method in ("constant","tilted_plane"):
        numbers = [
            row["models"][method]["background_at_target_nMgy_per_pixel"]
            for row in result
            if row["models"].get(method,{}).get("status") ==
            "FITTED_HOLDOUT_UNCALIBRATED"
        ]
        summary[method] = {
            "n_measurable_alternatives": len(numbers),
            "min_at_target_nMgy_per_pixel": (
                float(min(numbers)) if numbers else None),
            "max_at_target_nMgy_per_pixel": (
                float(max(numbers)) if numbers else None),
            "peak_to_peak_nMgy_per_pixel": (
                float(np.ptp(numbers)) if len(numbers) > 1 else None),
            "not_an_uncertainty_interval": True,
        }
    return {
        "schema": "desidera.sky.internal-grid-sensitivity/1",
        "status": "FIXED_TRIALS_DIAGNOSTIC_NOT_ADOPTED",
        "n_predeclared_variants": len(_SKY_FIXED_SENSITIVITY_TRIALS),
        "trials": result,
        "summary": summary,
        "original_tile_selection_unchanged": True,
        "galaxy_extinction_or_missing_upstream_flux_not_identified": True,
        "no_absolute_sky_calibration": True,
        "no_science_adoption": True,
    }


def _audit_current_project_sky_injection_recovery(
        image, inverse_variance, scene, candidates, *,
        anchor_native_xy, protected, source_mask, manual_mask, hard,
        settings):
    """Artificial exponential-light sensitivity on copied image arrays.

    Freeze the exact Scene, train/holdout and aperture for each experiment.
    This tests sky-fit-induced flux bias in a fixed *unmasked* aperture,
    NOT the actual galaxy photometry pipeline or historic DR11 sky removal.
    """
    from math import hypot, pi

    original = np.asarray(image)
    ivar = np.asarray(inverse_variance)
    ny,nx = original.shape
    ax,ay = map(float, anchor_native_xy)
    yy,xx = np.indices(original.shape,dtype=float)
    rr = np.hypot(xx-ax,yy-ay)
    effective_target_radius = float(
        np.sqrt(np.count_nonzero(protected)/pi))
    if effective_target_radius < 3.:
        return {
            "schema": "desidera.sky.frozen-mask-exponential-injection/1",
            "status": "NOT_MEASURED_NO_TARGET_PROTECTION",
            "scenarios": [],
            "science_adopted": False,
        }
    aperture_radius = 2.*effective_target_radius
    aperture = ((rr <= aperture_radius) & (~hard) & (~source_mask)
                & (~manual_mask) & np.isfinite(original)
                & np.isfinite(ivar) & (ivar > 0))
    n_geometric = int(np.count_nonzero(rr <= aperture_radius))
    n_usable = int(np.count_nonzero(aperture))
    if n_usable < 100:
        return {
            "schema": "desidera.sky.frozen-mask-exponential-injection/1",
            "status": "NOT_MEASURED_APERTURE_TOO_MASKED",
            "scenarios": [],
            "geometric_aperture_pixels": n_geometric,
            "unmasked_aperture_pixels": n_usable,
            "science_adopted": False,
        }
    ax_ap, ay_ap = xx[aperture],yy[aperture]
    scenarios = []
    # Fixed by protocol and labelled as uncalibrated, not a Leo P
    # morphological/stellar population model. q=0.65 test ellipse.
    for factor in (0.5, 1.0, 2.0):
        h = factor*effective_target_radius
        elliptical_radius = np.hypot(xx-ax,(yy-ay)/0.65)
        signal = 0.005*np.exp(-elliptical_radius/h)
        simulated = np.array(original,copy=True)
        simulated += signal
        aperture_input = float(np.sum(signal[aperture]))
        scenario = {
            "scale_length_as_protected_area_radius": factor,
            "scale_length_native_pix": h,
            "axis_ratio": .65,
            "central_amplitude_nMgy_per_pixel": .005,
            "psf_convolution_applied": False,
            "fixed_aperture_radius_native_pix": aperture_radius,
            "geometric_aperture_pixels": n_geometric,
            "unmasked_aperture_pixels": n_usable,
            "unmasked_aperture_fraction": float(n_usable/n_geometric),
            "injected_total_image_flux_nMgy":
                float(np.sum(signal)),
            "injected_flux_in_unmasked_aperture_nMgy": aperture_input,
            "zero_additional_correction_recovery_fraction": 1.,
            "models": {},
        }
        for method, base in candidates.items():
            fitted = fit_residual_background(
                simulated, ivar, scene,
                method=method, grid_id=scene.grid_id,
                observation_valid=~hard, settings=settings)
            if (fitted.fit_status != FitStatus.FITTED or
                    fitted.holdout_status != HoldoutStatus.MEASURED_UNCALIBRATED):
                scenario["models"][method] = {
                    "status": "INJECTED_FIT_NOT_VALIDATED",
                    "fit_status": fitted.fit_status.value,
                    "holdout_status": fitted.holdout_status.value,
                }
                continue
            fit_change = fitted.evaluate_at(ax_ap,ay_ap) - (
                base.evaluate_at(ax_ap,ay_ap))
            modeled_flux_loss = float(np.sum(fit_change))
            fraction = (1.-modeled_flux_loss/aperture_input
                        if aperture_input > 0 else None)
            scenario["models"][method] = {
                "status": "INJECTED_FIT_MEASURED_UNCALIBRATED",
                "sky_fit_shift_at_target_nMgy_per_pixel": float(
                    fitted.evaluate_at(ax,ay)-base.evaluate_at(ax,ay)),
                "induced_sky_flux_bias_in_unmasked_aperture_nMgy":
                    modeled_flux_loss,
                "injection_flux_recovery_fraction_fixed_aperture":
                    float(fraction) if fraction is not None else None,
                "zero_sky_correction_control_fraction": 1.,
                "absolute_galaxy_flux_bias_identified": False,
            }
        scenarios.append(scenario)
    return {
        "schema": "desidera.sky.frozen-mask-exponential-injection/1",
        "status": "SAME_COADD_FROZEN_MASK_INJECTION_DIAGNOSTIC_ONLY",
        "same_native_original_scene_and_heldout_tiles": True,
        "original_arrays_mutated": False,
        "model_selection_and_scene_masks_frozen": True,
        "no_noise_or_psf_realization_synthesized": True,
        "not_full_galaxy_detections_or_photometric_pipeline": True,
        "does_not_test_historic_DR11_sky_subtraction": True,
        "unmasked_aperture_only_not_total_galaxy_flux": True,
        "scenarios": scenarios,
        "science_adopted": False,
    }


def diagnose_current_project_residual_sky(dataset, draft, contamination):
    """Automatically select Scene-safe sky tiles and compare residual models.

    No manual sky drawing, FITS reads, network access or new survey bricks.
    All fitted results remain explicitly non-adopted research diagnostics.
    """
    from dataclasses import asdict, replace
    from lsbpipe.photometry.sky_masks import evaluate_mask_draft
    from lsbpipe.scene import grid_fingerprint

    if (not isinstance(dataset, dict) or
            not isinstance(contamination, dict) or draft is None):
        raise ValueError("Loaded Scene, draft and prepared masks are required")
    pid = dataset["project_id"]
    band = dataset["band"]
    source = dataset["source"]
    if (draft.project_id != pid or
            contamination.get("project_id") != pid or
            contamination.get("band") != band or
            contamination.get("dataset_source") != source or
            contamination.get("draft_digest") != draft.digest or
            contamination.get("status") !=
            "CATALOGUE_CONTAMINATION_MASK_NOT_ADOPTED"):
        raise ValueError("Stale project/band/draft/Tractor masking provenance")
    image = np.asarray(dataset["data"]["image"])
    ivar = np.asarray(dataset["data"]["invvar"])
    hard = np.asarray(dataset["native_hard"])
    shape = tuple(dataset["pixel_shape"])
    if (len(shape) != 2 or image.shape != shape or ivar.shape != shape
            or hard.shape != shape or hard.dtype.kind != "b"):
        raise ValueError("Original native quartet has inconsistent pixel grids")
    auto = np.asarray(contamination["auto_mask"])
    protected = np.asarray(contamination["protected"])
    prepared_hard = np.asarray(contamination["hard"])
    for array in (auto, protected, prepared_hard):
        if array.shape != shape or array.dtype.kind != "b":
            raise ValueError("Prepared masking does not match original native grid")
    if not np.array_equal(hard, prepared_hard):
        raise ValueError("Prepared native-hard mask differs from original source")
    if not isinstance(contamination.get("catalog_sha256"), str):
        raise ValueError("Prepared catalogue fingerprint missing")

    # Sky annotations from previous experiments are optional and are
    # ignored by the automatic selector; no draft edits are deleted.
    working_draft = replace(
        draft, edits=tuple(edit for edit in draft.edits
                           if edit.purpose not in ("sky_train","sky_holdout")))
    preview = evaluate_mask_draft(
        working_draft, dataset["wcs"], shape, band,
        native_hard=hard, automatic_contaminants=auto)
    sky_train, sky_holdout, selection = _automatic_scene_sky_pixels(
        dataset, preview, contamination)
    sky_protected = protected | preview["target_protect"] | preview["isophote_candidate"]
    grid = grid_fingerprint(dataset["wcs"], shape)
    warnings = [
        "AUTOMATIC_SKY_SAMPLE_SELECTION_INTERNAL_SAME_COADD_ONLY",
        "PROVISIONAL_CATALOGUE_MASK_NOT_PSF_WING_CALIBRATED",
        "GALAXY_OUTSKIRTS_MAY_EXTEND_BEYOND_PROTECTED_REGION",
        "BRIGHT_GAIA_EXCLUSION_HEURISTIC_NOT_A_PSF_WING_MODEL",
        "DESI_COADD_PIXEL_COVARIANCE_NOT_CALIBRATED",
        "UPSTREAM_SURVEY_SKY_SUBTRACTION_IRREVERSIBLE",
    ]
    if not selection["bright_gaia_inventory_present"]:
        warnings.append("BRIGHT_GAIA_SOURCE_INVENTORY_ABSENT")
    base = {
        "schema": "desidera.sky.project-native-residual-diagnostic/2",
        "status": "NOT_MEASURED",
        "science_status": "UNADOPTED_DIAGNOSTIC_ONLY",
        "project_id": pid,
        "band": band,
        "dataset_source": source,
        "grid_fingerprint": grid,
        "draft_sha256": draft.digest,
        "catalog_sha256": contamination["catalog_sha256"],
        "native_pixel_shape": list(shape),
        "native_hard_pixels": int(np.count_nonzero(hard)),
        "provisional_catalogue_mask_pixels": int(np.count_nonzero(auto)),
        "protected_target_pixels": int(np.count_nonzero(sky_protected)),
        "sky_train_effective_pixels": selection["training_pixels"],
        "sky_holdout_effective_pixels": selection["holdout_pixels"],
        "automatic_sky_selection_qa": selection,
        "preview_gaia_inner_mask_used":
            bool(contamination.get("gaia_inner_preview_enabled", False)),
        "preview_red_extended_mask_used":
            bool(contamination.get("red_extended_preview_enabled", False)),
        "uses_only_existing_project_loaded_native_band": True,
        "additional_bricks_or_catalogues_used": False,
        "original_fits_or_project_draft_modified": False,
        "background_subtracted_image_created": False,
        "sky_model_adopted": False,
        "psf_wings_validated": False,
        "noise_covariance_validated": False,
        "sky_galaxy_degeneracy_resolved": False,
        "source_masks_scientifically_validated": False,
        "methods": {},
        "method_comparison": None,
        "provisional_automatic_model": None,
        "automatic_sky_spatial_qa": None,
        "zero_additional_correction_baseline": None,
        "model_improvement_vs_zero_baseline": None,
        "recommended_science_action": "NO_ADDITIONAL_SKY_SUBTRACTION_UNTIL_VALIDATED",
        "automatic_sky_background_correction_authorized": False,
        "sky_grid_sensitivity_qa": None,
        "sky_injection_recovery_qa": None,
        "warnings": warnings,
    }
    if not selection["geometry_sufficient_for_diagnostic"]:
        return dict(
            base, status="INSUFFICIENT_AUTOMATIC_SCENE_SKY_GEOMETRY",
            warnings=warnings+["NO_SPATIALLY_DISTRIBUTED_CLEAN_SKY_TILES"])
    if not selection["bright_gaia_inventory_present"]:
        return dict(
            base, status="SOURCE_INVENTORY_REQUIRED_FOR_AUTOMATIC_SKY",
            warnings=warnings+["BRIGHT_SOURCE_HALO_VETO_UNAVAILABLE"])

    scene = _make_automatic_sky_scene(
        project_id=pid, grid=grid, band=band, shape=shape,
        revision=draft.revision, hard=hard, source_mask=auto,
        manual_mask=preview["manual_exclude"],
        target_mask=sky_protected,
        sky_train=sky_train, sky_holdout=sky_holdout)
    settings = BackgroundSettings(
        min_training_pixels=100, min_holdout_pixels=100,
        min_holdout_quadrants=3)
    candidates = {}
    methods = {}
    for method in ("constant","tilted_plane"):
        candidate = fit_residual_background(
            image, ivar, scene, method=method, grid_id=grid,
            observation_valid=~hard, settings=settings)
        candidates[method] = candidate
        serialized = asdict(candidate)
        serialized["fit_status"] = candidate.fit_status.value
        serialized["holdout_status"] = candidate.holdout_status.value
        serialized["warnings"] = list(candidate.warnings)
        cond = serialized.get("design_condition")
        if cond is not None and not np.isfinite(cond):
            serialized["design_condition"] = None
            serialized["warnings"].append("NONFINITE_CONDITION_NUMBER")
        methods[method] = serialized
    base["methods"] = methods
    baseline = _audit_zero_sky_baseline(
        image, sky_train, sky_holdout, selection)
    base["zero_additional_correction_baseline"] = baseline
    zero_rms = baseline["holdout_rms_nMgy_per_pixel"]
    # Numeric differences, not significance, are the only permissible
    # interpretation of a single-correlated-coadd pixel holdout.
    base["model_improvement_vs_zero_baseline"] = {
        method: {
            "rms_improvement_fraction": (
                float((zero_rms-c.holdout_rms_residual_nMgy)/zero_rms)
                if zero_rms > 0. and c.holdout_rms_residual_nMgy is not None
                else None),
            "improvement_is_not_statistical_significance": True,
        }
        for method,c in candidates.items()
    }
    tile_audit = _audit_automatic_sky_tile_residuals(
        image, sky_train, sky_holdout, selection, candidates)
    base["automatic_sky_spatial_qa"] = tile_audit
    base["warnings"] = warnings + tile_audit["warning_flags"]
    if not all(c.fit_status == FitStatus.FITTED for c in candidates.values()):
        base["status"] = "SKY_CANDIDATE_FIT_INCOMPLETE_NOT_ADOPTED"
        return base
    if not all(c.holdout_status == HoldoutStatus.MEASURED_UNCALIBRATED
               for c in candidates.values()):
        base["status"] = "SKY_HOLDOUT_GEOMETRY_FAILED_NOT_ADOPTED"
        return base

    anchor = contamination.get("anchor_icrs")
    if anchor is None:
        raise ValueError("Protected target ICRS anchor is unavailable")
    x, y = dataset["wcs"].all_world2pix(*map(float, anchor), 0)
    x, y = float(x), float(y)
    if not (np.isfinite(x) and np.isfinite(y)):
        raise ValueError("Project target WCS anchor is not finite")
    const, plane = candidates["constant"], candidates["tilted_plane"]
    b_const = float(const.evaluate_at(x,y))
    b_plane = float(plane.evaluate_at(x,y))
    ny, nx = shape
    corners = ((0,0),(nx-1,0),(0,ny-1),(nx-1,ny-1))
    differences = [
        float(plane.evaluate_at(cx,cy)-const.evaluate_at(cx,cy))
        for cx,cy in corners]
    # Select ONLY the provisional numerical candidate from the spatial
    # same-coadd holdout; never assert a scientific / noise-calibrated PASS.
    score_const = float(const.holdout_rms_residual_nMgy)
    score_plane = float(plane.holdout_rms_residual_nMgy)
    preferred = ("tilted_plane" if score_plane < .75*score_const
                 else "constant")
    base["status"] = "AUTOMATIC_SKY_CANDIDATES_MEASURED_NOT_ADOPTED"
    base["provisional_automatic_model"] = {
        "method": preferred,
        "selection_rule": (
            "PLANE_IF_HELDOUT_RMS_LESS_THAN_75_PERCENT_OF_CONSTANT_ELSE_CONSTANT"),
        "same_coadd_holdout_not_independent_validation": True,
        "not_a_science_adoption": True,
        "choice_only_between_two_fitted_nonzero_candidates": True,
        "zero_correction_reference_was_not_in_the_original_selection_rule": True,
        "model_selection_25percent_rule_is_heuristic_not_calibrated": True,
        "scientifically_recommended_for_subtraction": False,
        "spatial_support_warning_flags":
            tile_audit["warning_flags"],
        "spatial_quality_not_certified": True,
        "coefficients_nMgy": list(candidates[preferred].coefficients_nMgy),
        "reference_xy_pix": list(candidates[preferred].reference_xy_pix),
        "predicted_residual_sky_at_target_nMgy_per_pixel":
            float(candidates[preferred].evaluate_at(x,y)),
    }
    base["method_comparison"] = {
        "target_anchor_native_xy_pix":[x,y],
        "constant_at_target_nMgy_per_pixel":b_const,
        "plane_at_target_nMgy_per_pixel":b_plane,
        "plane_minus_constant_at_target_nMgy_per_pixel":b_plane-b_const,
        "max_abs_plane_minus_constant_corners_nMgy_per_pixel":
            max(abs(v) for v in differences),
        "both_methods_share_exact_same_sky_train_and_holdout":True,
        "not_a_galaxy_flux_bias_estimate":True,
        "not_a_model_selection_p_value":True,
        "holdout_rms_constant_nMgy_per_pixel":score_const,
        "holdout_rms_plane_nMgy_per_pixel":score_plane,
    }
    base["sky_grid_sensitivity_qa"] = (
        _audit_current_project_sky_grid_sensitivity(
            dataset, preview, contamination, grid=grid, band=band,
            shape=shape, revision=draft.revision, original_hard=hard,
            source_mask=auto, manual_mask=preview["manual_exclude"],
            target_mask=sky_protected, original_image=image,
            inverse_variance=ivar, anchor_native_xy=(x,y),
            original_candidates=candidates, settings=settings))
    base["sky_injection_recovery_qa"] = (
        _audit_current_project_sky_injection_recovery(
            image, ivar, scene, candidates,
            anchor_native_xy=(x,y), protected=sky_protected,
            source_mask=auto, manual_mask=preview["manual_exclude"],
            hard=hard, settings=settings))
    if sum(row["status"] == "BOTH_MODELS_MEASURED_UNCALIBRATED"
           for row in base["sky_grid_sensitivity_qa"]["trials"]) < 2:
        base["warnings"].append(
            "FEWER_THAN_TWO_PREDECLARED_GRID_VARIANTS_MEASURABLE")
    if base["sky_injection_recovery_qa"]["status"] != (
            "SAME_COADD_FROZEN_MASK_INJECTION_DIAGNOSTIC_ONLY"):
        base["warnings"].append(
            "SKY_INJECTION_RECOVERY_NOT_MEASURABLE")
    return base
