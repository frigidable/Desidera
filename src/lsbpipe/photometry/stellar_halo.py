"""SH1-A: native-pixel sector measurements for a saturated foreground star.

No PSF extrapolation, sky subtraction, halo fit or pixel modification.
Target-facing and protected pixels cannot train the future halo model.
The angular sectors are defined relative to the fixed star-target vector.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.ndimage import binary_dilation

EDGES_ARCSEC = (12., 16., 20., 25., 30., 40., 50., 60., 70., 85., 100., 115.)
MIN_PIXELS = 100
MIN_COVERAGE = .25


def _stats(image, select):
    v = image[select]
    if not v.size:
        return {"n": 0, "median_nMgy_pix": None, "mad_nMgy_pix": None,
                "q16_nMgy_pix": None, "q84_nMgy_pix": None}
    m = float(np.median(v))
    q = np.quantile(v, [.16, .84])
    return {"n": int(v.size), "median_nMgy_pix": m,
            "mad_nMgy_pix": float(np.median(np.abs(v-m))),
            "q16_nMgy_pix": float(q[0]), "q84_nMgy_pix": float(q[1])}


def measure_saturated_halo_sectors(image, invvar, native_hard, *,
        source_xy, target_xy, pixel_scale_arcsec, source_mask, protected,
        unresolved_sky_veto, neighbor_centres=()):
    a = np.asarray(image, dtype=float)
    iv = np.asarray(invvar, dtype=float)
    masks = [np.asarray(m, dtype=bool) for m in
             (native_hard, source_mask, protected, unresolved_sky_veto)]
    if a.ndim != 2 or iv.shape != a.shape or any(m.shape != a.shape for m in masks):
        raise ValueError("SH1 requires aligned native image, INVVAR and scene masks")
    sx, sy = map(float, source_xy)
    tx, ty = map(float, target_xy)
    scale = float(pixel_scale_arcsec)
    if (not all(math.isfinite(v) for v in (sx,sy,tx,ty,scale))
            or not .1 <= scale <= 1.
            or not (0 <= sx < a.shape[1] and 0 <= sy < a.shape[0])
            or math.hypot(tx-sx,ty-sy)*scale < 1.):
        raise ValueError("SH1 invalid source/target/scale")
    yy, xx = np.ogrid[:a.shape[0], :a.shape[1]]
    dx, dy = xx-sx, yy-sy
    distance = np.hypot(dx,dy)*scale
    direction = math.atan2(ty-sy, tx-sx)
    angle = (np.arctan2(dy,dx)-direction+np.pi/8)%(2*np.pi)
    sector = np.floor(angle/(np.pi/4)).astype(np.int8)
    hard, source, protection, veto = masks
    # Explicit engineering guard around invalid pixels. Do not allow
    # the protection mask to restore saturation/bleed pixels.
    guarded_hard = binary_dilation(hard, iterations=3)
    native_valid = (~guarded_hard & np.isfinite(a) & np.isfinite(iv) & (iv > 0))
    # Mask catalogue cores even when automatic photometry screening
    # left them protected or unfitted. This is only a halo measurement mask.
    cores = np.zeros(a.shape, dtype=bool)
    n_neighbours = 0
    for cx,cy in neighbor_centres:
        cx,cy = float(cx),float(cy)
        if not math.isfinite(cx) or not math.isfinite(cy):
            continue
        radius = 4./scale
        x0=max(0, math.floor(cx-radius)); x1=min(a.shape[1],math.ceil(cx+radius)+1)
        y0=max(0, math.floor(cy-radius)); y1=min(a.shape[0],math.ceil(cy+radius)+1)
        if x1 <= x0 or y1 <= y0:
            continue
        cyy,cxx=np.ogrid[y0:y1,x0:x1]
        cores[y0:y1,x0:x1] |= (cxx-cx)**2+(cyy-cy)**2 <= radius**2
        n_neighbours += 1
    descriptive = native_valid & ~source & ~cores
    train_valid = descriptive & ~protection & ~veto
    rows=[]
    for lo,hi in zip(EDGES_ARCSEC[:-1], EDGES_ARCSEC[1:]):
        annulus = (distance >= lo) & (distance < hi)
        nominal = math.pi*(hi**2-lo**2)/(8*scale**2)
        for k in range(8):
            geom = annulus & (sector == k)
            # only opposite three sectors train; side sectors 2,6
            # remain independent checks; 0,1,7 face the galaxy.
            role = "AWAY_TRAIN_CANDIDATE" if k in (3,4,5) else (
                "SIDE_CONTROL" if k in (2,6) else "TARGET_FACING_DESCRIPTIVE")
            desc = _stats(a, geom & descriptive)
            safe = _stats(a, geom & train_valid)
            fraction = min(1.,safe["n"]/nominal)
            rows.append({"r_inner_arcsec":lo, "r_outer_arcsec":hi,
                         "sector_index":k,"relative_angle_center_deg":45*k,
                         "role":role,"nominal_full_sector_pixels_approx":nominal,
                         "in_frame_pixels":int(geom.sum()),
                         "native_valid_pixels":int((geom & native_valid).sum()),
                         "descriptive":desc,"outside_protection_and_veto":safe,
                         "safe_full_sector_fraction_approx":fraction,
                         "adequate_for_descriptive_comparison":
                             safe["n"] >= MIN_PIXELS and fraction >= MIN_COVERAGE})
    return {"schema":"desidera.sh1-saturated-halo-sector-measurement/1",
            "status":"SH1_SECTORS_MEASURED_NOT_ADOPTED",
            "source_xy_native_pix":[sx,sy],"target_xy_native_pix":[tx,ty],
            "pixel_scale_arcsec":scale,"star_target_distance_arcsec":
                math.hypot(tx-sx,ty-sy)*scale,
            "radial_edges_arcsec":list(EDGES_ARCSEC),
            "native_invalid_guard_pix":3,"catalogue_core_radius_arcsec":4.,
            "catalogue_neighbour_count":n_neighbours,
            "min_pixels":MIN_PIXELS,"min_full_sector_fraction_approx":MIN_COVERAGE,
            "rows":rows,
            "catalog_flux_used":False,"compact_psf_used":False,
            "sky_subtracted":False,"halo_fitted":False,"science_adopted":False,
            "original_pixels_modified":False,"noise_covariance_calibrated":False,
            "warning":"Signed coadd medians include unknown sky; coverage and MAD are descriptive, not calibrated uncertainties. Target-facing values may include galaxy light. No halo amplitude or galaxy correction inferred."}
