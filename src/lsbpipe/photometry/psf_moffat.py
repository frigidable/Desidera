"""PSF-M1: finite, native pixel-integrated Moffat diagnostic.

Fit shared elliptical shape and individual centroids, fluxes and residual
planes to TRAINING pixels only. Diagonal weights are not a calibrated
coadd likelihood. Never use the analytic extrapolation as a bright-star halo.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.optimize import least_squares
from scipy.sparse import lil_matrix

# Gauss-Legendre quadrature on one native pixel, not image resampling.
_nodes, _weights = np.polynomial.legendre.leggauss(3)
_offsets = [(u/2, v/2, wu*wv/4)
            for u, wu in zip(_nodes, _weights)
            for v, wv in zip(_nodes, _weights)]


def pixel_moffat(x, y, shape):
    """shape=(log axial scale FWHM x,y, cross term rho, beta)."""
    wx, wy = np.exp(shape[:2])
    rho, beta = shape[2:]
    factor = 2*np.sqrt(np.expm1(np.log(2)/beta))
    ax, ay = wx/factor, wy/factor
    norm = (beta-1)*np.sqrt(1-rho*rho)/(np.pi*ax*ay)
    z = np.zeros_like(x, dtype=float)
    for dx, dy, w in _offsets:
        u, v = (x+dx)/ax, (y+dy)/ay
        z += w*norm*(1+u*u+v*v+2*rho*u*v)**(-beta)
    return z


def fit_joint_moffat(training, guessed_fwhm, stamp_size, fit_radius):
    """Profile linear nuisance parameters and jointly fit nonlinear shape/centres.

    Each input is (star_id, cutout). Fit only fixed valid pixels inside
    fit_radius from the DAO centroid; do not change membership on residuals.
    Flux is constrained nonnegative; sky correction is an unconstrained plane.
    """
    if len(training) < 4:
        raise ValueError("Moffat diagnostic requires at least four training stars")
    prepared = []
    for sid, c in training:
        image, weight = np.asarray(c["image"], float), np.asarray(c["weights"], float)
        if image.shape != (stamp_size, stamp_size) or weight.shape != image.shape:
            raise ValueError("Moffat training stamp shape mismatch")
        yy, xx = np.indices(image.shape, dtype=float)
        cx, cy = c["center_in_cutout"]
        x, y = xx-cx, yy-cy
        good = (weight > 0) & np.isfinite(weight) & np.isfinite(image)
        good &= x*x+y*y <= fit_radius**2
        if good.sum() < 35:
            raise ValueError("Insufficient weighted Moffat training pixels")
        x, y, w, data = x[good], y[good], weight[good], image[good]
        plane = np.column_stack((np.ones(len(x)), x/fit_radius, y/fit_radius))*w[:,None]
        prepared.append((sid, x, y, w, data*w, plane))
    starts = np.cumsum([0]+[len(p[1]) for p in prepared])
    def solve_one(p, shape, offset):
        sid, x, y, w, data, plane = p
        star = pixel_moffat(x-offset[0], y-offset[1], shape)*w
        design = np.column_stack((star, plane))
        # Scale each column for a well-conditioned four-parameter solve.
        scales = np.sqrt(np.sum(design**2, axis=0))
        mat = design/scales
        coef = np.linalg.lstsq(mat, data, rcond=None)[0]/scales
        if coef[0] <= 0:
            sky = np.linalg.lstsq(plane, data, rcond=None)[0]
            coef = np.r_[0., sky]
        residual = np.einsum("ij,j->i", design, coef)-data
        return residual, coef
    def residual(par):
        return np.concatenate([solve_one(p, par[:4], par[4+2*i:6+2*i])[0]
                               for i,p in enumerate(prepared)])
    n = len(prepared)
    jac = lil_matrix((int(starts[-1]), 4+2*n), dtype=np.int8)
    for i in range(n):
        jac[starts[i]:starts[i+1], :4] = 1
        jac[starts[i]:starts[i+1], 4+2*i:6+2*i] = 1
    init = np.r_[np.log([guessed_fwhm, guessed_fwhm]), 0., 3., np.zeros(2*n)]
    lo = np.r_[np.log([.5*guessed_fwhm]*2), -.65, 1.2, np.full(2*n,-1.5)]
    hi = np.r_[np.log([1.6*guessed_fwhm]*2), .65, 12., np.full(2*n,1.5)]
    fit = least_squares(residual, init, bounds=(lo,hi), jac_sparsity=jac.tocsr(),
                        x_scale="jac", ftol=1e-5, xtol=1e-5, gtol=1e-5,
                        max_nfev=100)
    if not fit.success or not np.isfinite(fit.x).all():
        raise ValueError("Joint Moffat fit did not converge: "+str(fit.message))
    wx, wy = np.exp(fit.x[:2])
    rho, beta = fit.x[2:4]
    q = np.array([[1/wx**2, rho/(wx*wy)], [rho/(wx*wy),1/wy**2]])
    eig, vec = np.linalg.eigh(q)
    widths = 1/np.sqrt(eig)
    angle = math.degrees(math.atan2(vec[1,0],vec[0,0]))%180
    per_star = []
    for i,p in enumerate(prepared):
        delta = fit.x[4+2*i:6+2*i]
        err, coef = solve_one(p, fit.x[:4], delta)
        if not np.isfinite(coef).all() or coef[0] <= 0:
            raise ValueError("Joint Moffat training flux is nonpositive: "+p[0])
        per_star.append({"star_id":p[0], "centroid_offset_xy_pix":delta.tolist(),
                         "analytic_flux_nMgy_model_dependent":float(coef[0]),
                         "residual_sky_plane_nMgy":[float(coef[1]),float(coef[2]/fit_radius),float(coef[3]/fit_radius)],
                         "n_pixels":len(err), "diagonal_weighted_rss":float(err@err)})
    yy, xx = np.indices((stamp_size,stamp_size), dtype=float)
    h = (stamp_size-1)/2
    raw = pixel_moffat(xx-h, yy-h, fit.x[:4])
    finite_sum = float(raw.sum())
    if not np.isfinite(raw).all() or finite_sum <= 0 or np.any(raw < 0):
        raise ValueError("Invalid finite Moffat kernel")
    return raw/finite_sum, {
        "method":"PSF_M1_JOINT_ELLIPTICAL_MOFFAT/1",
        "shape_parameters":fit.x[:4].tolist(), "beta":float(beta),
        "continuous_fwhm_major_pix":float(widths[0]),
        "continuous_fwhm_minor_pix":float(widths[1]),
        "continuous_axis_ratio_minor_over_major":float(widths[1]/widths[0]),
        "major_axis_angle_native_x_toward_y_deg":float(angle),
        "finite_stamp_fraction_of_analytic_integral":finite_sum,
        "fit_radius_pix":float(fit_radius), "pixel_quadrature_order":3,
        "n_function_evaluations":int(fit.nfev), "solver_status":int(fit.status),
        "solver_message":str(fit.message), "optimality":float(fit.optimality),
        "active_bounds":fit.active_mask.tolist(),
        "training_only":True, "training_stars":per_star,
        "diagonal_weighted_rss":float(2*fit.cost),
        "noise_covariance_calibrated":False, "wide_wings_validated":False,
        "analytic_flux_is_not_adopted_total_flux":True,
        "science_adopted":False,
    }
