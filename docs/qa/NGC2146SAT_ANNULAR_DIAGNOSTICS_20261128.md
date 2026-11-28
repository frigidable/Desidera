# NGC2146sat-derived experimental annular diagnostics (2026-11-28)

Status: **DIAGNOSTIC_NOT_ADOPTED / SCIENCE HOLD**. The new numerical
components are opt-in research APIs, not a replacement for existing
photometry or a new default GUI workflow.

## Implemented components

| Module under src/lsbpipe/photometry | Role |
| --- | --- |
| uncertainty_budget.py | Conditional magnitude covariance with shared extinction and explicit zeropoint covariance; finite tested scenario extrema |
| stellar_contamination.py | Native-WCS stellar core and directed tapered-ray proposals; compact radial halo template; rank/condition-gated independent far-sky and near-halo block fits with buffered held-out sectors |
| radial_models.py | Explicit smooth nested ellipse hypotheses with conservative non-crossing gate; radial Sérsic intensity and untruncated model integration; signed observed growth and conditional plateau check |
| annular_inference.py | Contiguous total/ring/sector support; finite PSF adjoint; linear covariance propagation; bounded GLS and correlated holdout predictions; ring-domain control comparison |
| native_annular_diagnostic.py | Explicit native-frame/scene/correction/PSF identity checks; signed ring measurements and joint conditional Sérsic fit on exact common observed support |

All use the existing NativeBandFrame and SceneState contracts where applicable.
Original FITS remain immutable. No new external dependency is introduced.
No default star radii, sky guards, Sérsic bounds, calibration errors or
covariance kernel parameters are taken from NGC2146sat.

## Opt-in workflow

1. Load native data using the existing validated FITS loader. Preserve the
   frame input SHA-256, WCS fingerprint, release, region, brick and band.
2. Supply explicit star positions and core/ray hypotheses. Review the
   MaskLayer proposals before adding them to a SceneState. Native hard
   exclusions are always preserved by the new diagnostic.
3. Supply a strict galaxy guard and a far-star guard outside the chosen
   compact template. fit_native_independent_sky_halo returns provisional
   NativeCorrection rasters, the proposed core/ray layers, disjoint block
   samples and frozen linear influence matrices. Repeat with predeclared
   held-out sectors; never select a model by galaxy curve-of-growth shape.
   The halo contains NaN inside its excluded core by design. The strict
   implementation does not implicitly subtract a galaxy model.
4. Supply a reviewed NestedEllipse hypothesis, exact native scene,
   residual sky and any explicitly chosen source/halo corrections, and an
   explicitly same-frame empirical PSF. prepare_native_annuli returns
   signed observed fluxes, formal diagonal variance, per-sector coverage,
   and a finite PSF-response matrix. Missing rings are NaN, not zero.
   Coverage is in-frame only: geometric containment needs a separate check.
   Row chunks and an explicit default 2 million work-pixel gate limit
   memory; use validated native cutouts without changing original FITS.
5. Supply the full conditional covariance in explicit band/ring order.
   fit_common_support_sersic requires identical WCS, release/brick,
   geometry and exact valid-pixel support across bands. It uses caller
   supplied log-parameter bounds. Model extrapolations and measured
   fluxes have different meanings. No covariance rescaling by reduced
   chi-square is performed.
6. For external matched control vectors use control_model_comparison.
   Controls must have unique identities, identical declared support and
   nonoverlapping full footprints including sky-training regions. Counts
   of exceedances are not a calibrated false-positive probability.
7. Report conditional covariance and finite scenario extrema separately.
   Do not infer a probability distribution from sensitivity cases.

## Regression anchor, not scientific acceptance

NGC2146sat, DR11 North, brick 0965p780, reference
core26_compact_refit_50_L15_A50. Working hypotheses: 26 arcsec core,
tapered north ray, compact g/r halo cutoff 70/100 arcsec, independent
plane, conditional L15_A50 covariance and eligible radius 50 arcsec.
These are object-specific and must never become general defaults.

The numerical calibration regression uses the supplied
NGC2146sat_uncertainty_budget.zip manifest:

- core-mask-gls.zip SHA256:
  b677b20fa4d53beea80bf315fbeceeea1db0e00c36c2c7cc7cb125f684e38cb7
- core-publication-qa.zip SHA256:
  2c695197c2f679fb569383cc4211b6dcf31f81538e44a23c34eb2331798b775f
- outer-source-gls.zip SHA256:
  f27db2ad782425a9344b28c67bb4e76dfe9d6ab3346f2eaf23cfb8bcaee3e0f0
- full-brick-controls.zip SHA256:
  2153ed9c0b5b3a2255e8e05d356d5080ee518f43d92a8e6472b765d6aa879e84

Working color 0.313993 mag; conditional plus assumed calibration
sigma 0.095326 mag. Alternate halo changes the color substantially.
No independent CCD zeropoint calibration; no full confidence interval.
The observed outer excess is diagnostic, not an adopted second component.
Only two strict nonoverlapping controls at 80% coverage were available,
and none at 95%. Do not claim a confirmed physical envelope.

## Verification and outstanding gates

Tests are in tests/test_uncertainty_budget.py,
tests/test_stellar_annular_diagnostics.py, tests/test_radial_models.py
and tests/test_native_annular_diagnostic.py. They include synthetic
geometry/masks, direct vs adjoint PSF, signed flux, correlated holdout,
a real numerical covariance excerpt and a synthetic end-to-end native
scene. Execution status must be read from exact-commit CI, not inferred
from test presence.

Prior exact commit 724b1db9c581f271e5fdeff10a88464a3ad28df4 passed Survey
Explorer tests but failed the Windows offline GUI smoke suite in
test_sga_informed_radial_initial_ellipse_requires_explicit_review.
See Actions runs 46084871418 and 46084871407. No GUI changes in this series.

Still required before any default workflow/science adoption:
full exact-commit CI; full real-FITS replay against archived object
diagnostics; review of memory/edge support and sky/halo transfer;
calibration of joint block covariance and independent noise;
explicit fit/selection and regularization of varying geometry (the new
geometry class evaluates supplied hypotheses only); wider empirical
pixel injection/control selection; independent CCD calibration; full
publication/export and GUI workflow review. Existing diagnostic figures
remain external research products, not generated by this new API.
Never modify or merge main without approval.
