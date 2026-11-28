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


## Follow-up: spatial block covariance, archive replay and export

- `block_covariance.py` builds the conditional spatial covariance from
  supplied ring-by-tile counts, design sums and fixed fitted sky/halo
  influence matrices. It rejects overlaps between training tile footprints
  and observed aperture tiles, supports matched far-control residuals for
  robust marginal scales and cross-band correlations, and evaluates
  explicit exponential kernel length/mixing hypotheses. It does not select
  a covariance scenario automatically or certify independent controls.
  Tile arrays must still be assembled from the exact native pixel support;
  the automatic adapter from native scenes is not yet implemented.
- `annular_export.py` writes a new ZIP with signed ring measurements,
  fit JSON, numeric NPZ and a SHA-256 manifest. Optional matplotlib
  figures show observed-support curves of growth and signed mean radial
  intensity (not reconstructed surface brightness). Exact fit input and
  covariance digests must match. It is not the complete six-figure
  publication exporter and is not wired into the GUI.
- `scripts/replay_ngc2146sat_annular.py` reads the historical template NPZ
  member in full_brick_controls_script.zip, core-mask-gls.zip and
  core-publication-qa.zip with the original SHA-256 checks. It tests the
  new GLS solver and correlated holdout against the archived 20-ring
  observed vector and PSF operators. It never downloads FITS or reruns
  scene preparation.

Example (PowerShell, from repository checkout with the package installed;
archives stay in the existing project thesis_figures directory):

```powershell
$d = "C:\Users\frigi\DESIDERAProjects\NGC2146sat\thesis_figures"
python scripts/replay_ngc2146sat_annular.py --template-archive "$d\full_brick_controls_script.zip" --reference-archive "$d\core-mask-gls.zip" --holdout-archive "$d\core-publication-qa.zip" --output "$d\desidera-annular-replay.json"
```

An independent SciPy calculation on the archived arrays reproduced the
reference chi-square 65.35565473 (12 degrees of freedom), maximum ring
prediction difference 2.27e-6 nMgy, and holdout chi-square 18.76479051
with residuals [3.63703817, 1.46130443] nMgy. This is numerical evidence
for the ring-domain reference, not a claim that the exact-commit package
suite or the full native FITS replay has passed.

New tests: tests/test_block_covariance.py and archive identity/no-overwrite
checks in tests/test_native_annular_diagnostic.py. CI must be checked for
the final head; intermediate runs can be canceled by successive commits.
The prior GUI failure is a Windows native access violation during
test_sga_informed_radial_initial_ellipse_requires_explicit_review, with
the worker stack in ellipse_seed._outer_boundary. Root cause is not
established; do not disable the test or change scientific algorithms
merely to make CI pass.


## Native covariance adapter (follow-up)

The previously missing native adapter is now implemented as
`native_scene_block_covariance` in block_covariance.py. Opt-in sequence:

1. Fit each native band with fit_native_independent_sky_halo; use a
   common block lattice and explicit holdout angular sector. Far held
   blocks are retained separately, not used to fit sky or halo.
2. Review and add the proposed core/ray layers, then intersect native
   validity with common_native_validity. Use that mask and
   retain_pixel_support=True when preparing each annular record.
3. Pass the exact frames, records and independent halo fit results to
   native_scene_block_covariance. It verifies grid, data and correction
   identity, exact observed ring support, and disjoint training/control
   tiles. Far controls must match both tile positions and pixel-mask
   fingerprints across bands. Block means, counts and formal IVAR
   variances are retained in the halo-fit results.
4. Use fit_native_covariance_sersic to bind the conditional covariance
   to the exact annular records and export as before.

This is a NEW conditional covariance hypothesis, not an exact
reimplementation of the archived L15_A50 matrix. It models coherent
tile fluctuations as additive to independent pixel-IVAR noise. It
subtracts the median formal block-mean variance from the robust
paired far-control covariance, rejects unresolved or indefinite
coherent covariance, and propagates independent IVAR through ring sums
and fixed sky/halo influence matrices. No covariance eigenvalue
clipping, arbitrary jitter or automatic fit acceptance is used.
The control residuals still include uncertainty of predicted sky;
robust scales and clipped correlation are approximate, and spatial
stationarity, independence of pixel IVAR and transfer from held-out
control tiles to aperture tiles remain hypotheses. The model does not
propagate compact-source or geometry refits. An insufficient control
sample is a failure, not permission to relax the support criteria.

New tests exercise the noise decomposition, exact formal IVAR
propagation, and a synthetic native two-band pipeline. Their presence
does not establish execution. Real FITS replay, comparison against
archived covariance, empirical calibration, geometry optimization,
and GUI integration remain open. Do not update adopted photometry.


## Regularized geometry inference (follow-up)

`geometry_inference.py` adds three PSF-forward hypotheses on the same
predeclared native block support: fixed reference, free constant ellipse,
and two smoothly varying nested endpoints. The geometry midpoint and
endpoint differences are bounded; the entire variable parameter box
must satisfy the conservative non-crossing certificate of NestedEllipse.
The user supplies explicit regularization scales, profile bounds, block
size, minimum usable fraction, holdout buffer and block intensity scales.
There is no automatic geometry selection or calibration of the diagonal
block score as a chi-square. The penalty is recorded separately.
A finite PSF convolution footprint is required inside the native frame.

The native adapter verifies frame, WCS, mask and correction/PSF identities,
uses common_native_validity, and preserves the split and geometry
configuration. The selected geometry can be returned only after explicit
review via reviewed_geometry_candidate; this does not adopt a measurement.
Geometry candidates with failed or boundary-constrained optimization
are rejected by that review helper. The entire result is diagnostic.

Historical NGC2146sat_reggeom_fits.csv, variable, all, g, gives endpoint
q=0.5604268448317696 and 0.67487168 (rounded), PA=44.867813 and
50.866228 degrees, and transition 10–70 arcsec. The conservative nesting
margin for these archived parameters is 0.14809569. An independent
calculation on 100000 uniform tangent-plane test positions reproduced
the original 52-step bisection using the new endpoint/transition
shortcut to max absolute 2.73e-11 arcsec. This checks the radius algorithm,
not an end-to-end real FITS geometry fit or a package CI result.

Regression tests in tests/test_geometry_inference.py cover synthetic
PSF-convolved geometry recovery, radius equality, holdout leakage,
bounds, native identity and explicit review. Full CI is still required.
At exact prior commit cac0414fcc7b7f9a3f81c20c4c3dbdd1328c1c2e,
Ubuntu Python 3.12 reported 1058 passed, 49 skipped and one failure:
test_native_joint_workflow_on_synthetic_frames had insufficient inner
ring sector support. The test fixture now samples the same 40 arcsec
field at 0.5 arcsec/pixel instead of 1 arcsec/pixel, preserving the
coverage gate and a hard exclusion. This repair requires exact-commit CI.
The prior Windows GUI access violation remains unresolved.

Outstanding: evaluate geometry candidates on real native frames with
predeclared bounds/splits, joint empirical covariance calibration,
full real-data replay, profile family and pixel injection studies, six
publication figures, end-to-end workspace orchestration and GUI review.
No new method is science-adopted.
