"""Explicit stellar mask and compact-halo diagnostics on a native WCS.

No parameters are inferred from galaxy photometry; no automatic adoption.
The halo is a finite effective model, not a physical PSF wing measurement.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np

from lsbpipe.scene import LayerKind, MaskLayer, grid_fingerprint

STATUS = "DIAGNOSTIC_NOT_ADOPTED"


def _number(x):
    if isinstance(x, (bool, np.bool_, str)):
        raise ValueError("Expected a finite numerical parameter")
    x = float(x)
    if not math.isfinite(x):
        raise ValueError("Expected a finite numerical parameter")
    return x


@dataclass(frozen=True)
class RayMask:
    """Directed ray; PA east of north, knots (distance, full width) in arcsec."""
    pa_deg: float
    knots: tuple[tuple[float, float], ...]

    def __post_init__(self):
        pa = _number(self.pa_deg)
        knots = tuple(tuple(map(_number, k)) for k in self.knots)
        if len(knots) < 2 or any(len(k) != 2 for k in knots):
            raise ValueError("At least two distance/width knots required")
        d, w = np.asarray(knots).T
        if d[0] < 0 or np.any(np.diff(d) <= 0) or np.any(w <= 0):
            raise ValueError("Increasing nonnegative distances and positive widths required")
        object.__setattr__(self, "pa_deg", pa % 360)
        object.__setattr__(self, "knots", knots)


def stellar_exclusion_layers(frame, *, source_id, ra_deg, dec_deg,
                             core_arcsec, rays=(), row_chunk=128):
    """Create proposed STAR_CORE/INSTRUMENTAL_DEFECT layers, without unmasking.

    Sky offsets are computed on this frame's WCS (not copied across bands).
    Ray widths interpolate linearly between explicit knots; no extrapolation.
    Core/ray radii are user-supplied hypotheses, not automatic recommendations.
    """
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    ra, dec, radius = map(_number, (ra_deg, dec_deg, core_arcsec))
    if not 0 <= ra < 360 or not -90 <= dec <= 90 or not 0 < radius <= 1800:
        raise ValueError("Invalid stellar center or core radius")
    if not isinstance(row_chunk, int) or isinstance(row_chunk, bool) or row_chunk < 1:
        raise ValueError("Invalid chunk size")
    rays = tuple(rays)
    if any(not isinstance(ray, RayMask) for ray in rays):
        raise ValueError("Expected explicit RayMask objects")
    wcs = frame.wcs
    shape = frame.shape
    grid = grid_fingerprint(wcs, shape)
    center = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame="icrs")
    core = np.zeros(shape, dtype=bool)
    spikes = np.zeros(shape, dtype=bool)
    for start in range(0, shape[0], row_chunk):
        stop = min(shape[0], start + row_chunk)
        y, x = np.mgrid[start:stop, :shape[1]]
        sky = wcs.pixel_to_world(x, y).icrs
        separation = center.separation(sky).arcsec
        pa = center.position_angle(sky).rad
        core[start:stop] = separation <= radius
        for ray in rays:
            angle = pa - np.deg2rad(ray.pa_deg)
            along = separation*np.cos(angle)
            cross = separation*np.sin(angle)
            d, width = np.asarray(ray.knots).T
            half = np.interp(along, d, width)/2
            spikes[start:stop] |= (along >= d[0]) & (along <= d[-1]) & (np.abs(cross) <= half)
    source = f"Explicit stellar hypothesis {source_id}; {STATUS}"
    return (
        MaskLayer(f"{source_id}_core", LayerKind.STAR_CORE, core, grid, frame.band, source),
        MaskLayer(f"{source_id}_rays", LayerKind.INSTRUMENTAL_DEFECT, spikes, grid, frame.band, source),
    )


def compact_halo_template(radius_arcsec, *, cutoff_arcsec, power,
                          pivot_arcsec, core_arcsec):
    """Unit-pivot template outside the excluded core, C1 at outer cutoff.

    Returns NaN inside the core: those pixels must stay excluded, never filled.
    Normalization is an analytic convention, not a catalogued stellar flux.
    """
    cutoff, p, pivot, core = map(_number, (cutoff_arcsec, power, pivot_arcsec, core_arcsec))
    if not 0 < core <= pivot < cutoff or p < 0:
        raise ValueError("Require 0 < core <= pivot < cutoff and nonnegative power")
    r = np.asarray(radius_arcsec, dtype=float)
    if not np.isfinite(r).all() or np.any(r < 0):
        raise ValueError("Invalid radial coordinates")
    out = np.zeros_like(r)
    out[r < core] = np.nan
    use = (r >= core) & (r < cutoff)
    out[use] = (r[use]/pivot)**(-p)*((1-r[use]/cutoff)/(1-pivot/cutoff))**2
    if not np.isfinite(out[r >= core]).all():
        raise ValueError("Nonfinite halo outside excluded core")
    return out


def _robust_fit(x, y, *, max_condition, max_iter=100):
    """Huber IRLS with an explicit rank/condition gate and frozen influence."""
    n, k = x.shape
    if n <= k or np.linalg.matrix_rank(x) != k:
        raise ValueError("Insufficient independent design support")
    w = np.ones(n)
    beta = np.linalg.lstsq(x, y, rcond=None)[0]
    for _ in range(max_iter):
        resid = y-x@beta
        scale = max(1.4826*np.median(np.abs(resid-np.median(resid))),
                    np.finfo(float).eps*max(1., np.max(np.abs(y))))
        w = np.minimum(1., 1.5*scale/np.maximum(np.abs(resid), np.finfo(float).tiny))
        xx = np.sqrt(w)[:, None]*x
        if np.linalg.cond(xx) > max_condition:
            raise ValueError("Unreliable spatial extrapolation: ill-conditioned design")
        update = np.linalg.lstsq(xx, np.sqrt(w)*y, rcond=None)[0]
        if np.linalg.norm(update-beta) <= 1e-10*max(1., np.linalg.norm(beta)):
            beta = update
            break
        beta = update
    else:
        raise ValueError("Robust fit did not converge")
    operator = np.linalg.pinv(np.sqrt(w)[:, None]*x)*np.sqrt(w)
    return beta, operator, float(np.linalg.cond(np.sqrt(w)[:, None]*x))


def fit_independent_sky_halo(block_values, sky_design, halo_template, *,
                             far, near, holdout, max_condition=100.):
    """Fit far sky then freeze it while fitting a single nonnegative halo amplitude.

    Inputs are means of DISJOINT same-frame blocks; masks must be predeclared.
    Sky design and template must be averaged over identical valid pixels.
    Caller must exclude galaxy light, rays, invalid pixels and saturated core.
    No implicit subtraction of a galaxy, source, or previously removed sky.
    Returned influence is a local fixed-Huber-weight derivative, not calibrated.
    """
    y, x, h = map(lambda v: np.asarray(v, dtype=float),
                  (block_values, sky_design, halo_template))
    n = y.size
    if y.ndim != 1 or x.ndim != 2 or x.shape[0] != n or x.shape[1] < 1 or h.shape != (n,):
        raise ValueError("Incompatible block arrays")
    if not all(np.isfinite(v).all() for v in (y, x, h)) or np.any(h < 0):
        raise ValueError("Nonfinite block input or negative halo template")
    masks = tuple(np.asarray(v) for v in (far, near, holdout))
    if any(m.shape != (n,) or m.dtype.kind != "b" or not m.any() for m in masks):
        raise ValueError("Explicit nonempty Boolean block selections required")
    far, near, holdout = masks
    if np.any(far & near) or np.any((far | near) & holdout):
        raise ValueError("Sky, halo and holdout blocks must be disjoint")
    if np.any(h[far] != 0):
        raise ValueError("Far sky must be outside effective halo support")
    if np.count_nonzero(h[near]) < 3:
        raise ValueError("Insufficient measured halo support")
    condition = _number(max_condition)
    if condition <= 1:
        raise ValueError("Invalid condition threshold")
    beta, sky_op, cond = _robust_fit(x[far], y[far], max_condition=condition)
    residual = y-x@beta
    alpha, halo_op, _ = _robust_fit(h[near, None], residual[near], max_condition=condition)
    if alpha[0] <= 0:
        raise ValueError("Halo amplitude unconstrained at nonnegative boundary")
    # beta derivative in all input blocks; amplitude includes sky coupling.
    b = np.zeros((x.shape[1], n))
    b[:, far] = sky_op
    a = np.zeros(n)
    a[near] = halo_op[0]
    a -= (halo_op[0]@x[near])@b
    prediction = x@beta+h*alpha[0]
    e = y[holdout]-prediction[holdout]
    return {"status": STATUS, "sky_coefficients": beta, "amplitude": float(alpha[0]),
            "sky_influence": b, "halo_influence": a,
            "sky_condition": cond, "holdout_residual": e,
            "holdout_rmse": float(np.sqrt(np.mean(e**2))),
            "counts": tuple(int(m.sum()) for m in masks),
            "covariance_calibrated": False}


def fit_native_independent_sky_halo(frame, scene, *, star_icrs, galaxy_icrs,
                                     core_arcsec, cutoff_arcsec, power,
                                     pivot_arcsec, galaxy_guard_arcsec,
                                     far_star_guard_arcsec, rays, block_size,
                                     minimum_block_fraction, holdout_pa_deg,
                                     holdout_width_deg, buffer_deg,
                                     sky_degree=1, compact_sources=None,
                                     max_condition=100., max_work_pixels=2_000_000):
    """Strict independent block fit from one native observation.

    Entire block footprints obey both geometric guards; only usable pixels
    contribute to means. Galaxy light is never fitted/subtracted implicitly.
    A spatial sector and buffer are excluded from BOTH training stages.
    The caller declares the geometry/guards and repeats folds explicitly.
    No automatic family selection, correction adoption, or calibrated errors.
    """
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    from lsbpipe.scene import SceneState, evaluate_mask, Operation
    from .contracts import NativeBandFrame
    from .native_annular_diagnostic import NativeCorrection

    if not isinstance(frame,NativeBandFrame) or not isinstance(scene,SceneState):
        raise ValueError("Validated frame and explicit scene required")
    grid=frame.native_grid_fingerprint
    digest=frame.native_input_sha256
    if scene.grid_id!=grid or scene.shape!=frame.shape or scene.band!=frame.band:
        raise ValueError("Frame/scene identity mismatch")
    if not isinstance(max_work_pixels,int) or max_work_pixels<1 or np.prod(frame.shape)>max_work_pixels:
        raise ValueError("Explicit work-pixel budget exceeded")
    if not isinstance(block_size,int) or isinstance(block_size,bool) or block_size<2:
        raise ValueError("Invalid block size")
    core,cutoff,pivot,p,gg,fg,fraction,hpa,hw,buf=map(_number,
        (core_arcsec,cutoff_arcsec,pivot_arcsec,power,galaxy_guard_arcsec,
         far_star_guard_arcsec,minimum_block_fraction,holdout_pa_deg,holdout_width_deg,buffer_deg))
    if (not 0<core<=pivot<cutoff<=fg or p<0 or gg<=0 or not 0<fraction<=1
            or not 0<hw<180 or buf<0 or hw/2+buf>=180 or sky_degree not in (0,1)):
        raise ValueError("Invalid guards, template, holdout or sky degree")
    ra,dec=map(_number,star_icrs)
    gra,gdec=map(_number,galaxy_icrs)
    if not 0<=ra<360 or not -90<=dec<=90 or not 0<=gra<360 or not -90<=gdec<=90:
        raise ValueError("Invalid ICRS center")
    # Existing scene and new core/ray hypotheses can only reduce validity.
    core_layer,ray_layer=stellar_exclusion_layers(frame,source_id="diagnostic_star",
        ra_deg=ra,dec_deg=dec,core_arcsec=core,rays=rays)
    base=frame.valid & ~frame.mask_reasons[f"SATUR_{frame.band.upper()}"]
    base &= ~frame.mask_reasons[f"ALLMASK_{frame.band.upper()}"]
    valid=evaluate_mask(scene,frame.image_nMgy,frame.invvar_nMgy2,grid_id=grid,
                        operation=Operation.PHOTOMETRY,observation_valid=base).valid
    # Respect any explicitly protected target region as well as the
    # independent radial guard used below.
    for layer in scene.layers:
        if layer.kind==LayerKind.TARGET_PROTECT:
            valid=valid & ~layer.data
    valid=valid & ~core_layer.data & ~ray_layer.data
    values=np.array(frame.image_nMgy,copy=True)
    if compact_sources is not None:
        c=compact_sources
        if (not isinstance(c,NativeCorrection) or c.kind!="compact_sources"
                or c.grid_id!=grid or c.input_sha256!=digest or c.values.shape!=frame.shape
                or not np.isfinite(c.values[valid]).all()):
            raise ValueError("Invalid explicitly supplied compact-source correction")
        values[valid]-=c.values[valid]
    yy,xx=np.indices(frame.shape)
    world=frame.wcs.pixel_to_world(xx,yy).icrs
    star=SkyCoord(ra*u.deg,dec*u.deg,frame="icrs")
    galaxy=SkyCoord(gra*u.deg,gdec*u.deg,frame="icrs")
    sr=star.separation(world)
    if np.max(sr.deg)>.5:
        raise ValueError("Small-field sky plane exceeds 0.5 degree")
    r=sr.arcsec
    angle=star.position_angle(world).rad
    gr=galaxy.separation(world).arcsec
    delta=np.abs((np.rad2deg(angle)-hpa+180)%360-180)
    x=r*np.sin(angle)/fg
    y=r*np.cos(angle)/fg
    template=compact_halo_template(r,cutoff_arcsec=cutoff,power=p,
                                   pivot_arcsec=pivot,core_arcsec=core)
    means=[]; designs=[]; halos=[]; masks=[]; tiles=[]; counts=[]
    ny,nx=frame.shape
    for j in range(0,ny-block_size+1,block_size):
        for i in range(0,nx-block_size+1,block_size):
            sl=np.s_[j:j+block_size,i:i+block_size]
            ok=valid[sl]
            if ok.sum()<fraction*block_size**2 or gr[sl].min()<gg or r[sl].min()<core:
                continue
            far=bool(r[sl].min()>=fg)
            near=bool(r[sl].max()<=cutoff)
            training=bool(delta[sl].min()>=hw/2+buf)
            held=bool(delta[sl].max()<=hw/2)
            mask=(far and training,near and training,near and held)
            if not any(mask):
                continue
            means.append(float(values[sl][ok].mean()))
            designs.append([1.] if sky_degree==0 else [1.,float(x[sl][ok].mean()),float(y[sl][ok].mean())])
            halos.append(float(template[sl][ok].mean()))
            masks.append(mask);tiles.append((i,j));counts.append(int(ok.sum()))
    if not masks:
        raise ValueError("No eligible independent blocks")
    far,near,holdout=np.asarray(masks,dtype=bool).T
    result=fit_independent_sky_halo(means,designs,halos,far=far,near=near,
                                   holdout=holdout,max_condition=max_condition)
    beta=result["sky_coefficients"]
    sky=np.full(frame.shape,beta[0]) if sky_degree==0 else beta[0]+beta[1]*x+beta[2]*y
    halo=result["amplitude"]*template
    provenance=(f"strict independent block fit; core={core}; cutoff={cutoff}; "
                f"power={p}; pivot={pivot}; galaxy_guard={gg}; far_guard={fg}; "
                f"holdout_PA={hpa}; width={hw}; buffer={buf}; block={block_size}; "
                f"valid_fraction={fraction}; {STATUS}")
    # NaN inside the candidate halo core is intentional. Applying these
    # corrections requires adding the same core/ray layers to the scene.
    result.update({"dataset_id":frame.dataset_id,"grid_id":grid,"native_input_sha256":digest,
        "sky":NativeCorrection("residual_sky",grid,digest,provenance,sky),
        "halo":NativeCorrection("stellar_halo",grid,digest,provenance,halo),
        "proposed_layers":(core_layer,ray_layer),"block_xy":np.asarray(tiles),
        "block_counts":np.asarray(counts),"far":far,"near":near,"holdout":holdout,
        "block_means":np.asarray(means),"block_sky_design":np.asarray(designs),
        "block_halo_template":np.asarray(halos),"compact_sources_sha256":
            None if compact_sources is None else compact_sources.sha256,
        "conditional_on_strict_galaxy_guard":True})
    return result
