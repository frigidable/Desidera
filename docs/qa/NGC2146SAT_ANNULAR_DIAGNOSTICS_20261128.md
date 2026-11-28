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


## Staged native workflow (follow-up)

The exact prior commit b1edfc217850dd065202be62e0e74a6834c5772e
passed all six Survey Explorer jobs (Ubuntu and Windows Python 3.10
and 3.12, package wheel and GUI import) and the Windows UI baseline.
Actions runs 46139519136 and 46139519113. The earlier native access
violation did not recur; its root cause has not been established.

`native_workflow.py` composes the previously separate research APIs:

- prepare_scene_candidates fits each band on its existing native scene,
  constructs candidate additional mask layers, and retains correction
  rasters, PSF identities and block fits. Original scenes stay unchanged.
- geometry_candidates optionally evaluates PSF-forward fixed, constant
  and varying geometry after explicit mask review.
- run_reviewed_annular_workflow requires mask confirmation and an
  explicitly reviewed NestedEllipse, computes common valid support,
  measures signed annuli, checks finite PSF footprint around the
  measured pixels and contiguous coverage of requested rings, builds
  conditional covariance, fits the joint Sérsic model and optionally
  writes a new diagnostic archive.

There is no fallback to a different radial limit, mask, geometry or
noise hypothesis. WorkflowBlocked carries the failed stage and reason
for expected scientific/input gate failures. Failed, boundary or
rank-deficient fits are not exported by this workflow as complete.
The PSF footprint check does not certify geometric containment of an
entire contour; in-frame coverage is still explicitly labelled as such.
The workflow uses an independently fitted single-star compact halo
per band and fixed compact-source corrections. It does not estimate
geometry uncertainty, select a covariance model, calibrate the noise,
or validate independent full-brick controls. No scientific acceptance
is implied by an archive or the DIAGNOSTIC_COMPLETE stage.

The in-memory candidate bundle is a research API, not a persistent
review token or an immutable serialized provenance artifact. Callers
must preserve the prepared results unchanged. Geometry provenance
from geometry_candidates should be preserved separately; passing an
explicit NestedEllipse is also supported. A dedicated GUI integration
and persistent review state remain outstanding.

tests/test_native_workflow.py adds a synthetic two-band end-to-end
fixture with different PSFs, bright star and correlated block noise,
plus review, coverage and archive tests. The new test execution
requires exact-commit CI. No real native FITS replay or additional
scientific adoption was performed in this follow-up.


## Asymptotic flux decomposition and differential injections (follow-up)

The exact prior commit 5e05ef5b77f7c47fd5fbba04e8c8aac336f3ca9a
had successful Ubuntu Python 3.10/3.12 scientific jobs and package wheel
in run 46145869317. The remaining jobs were canceled on subsequent
pushes. Do not attribute these results to the new head.

`asymptotic_diagnostic.py` requires the optional full in-frame
aperture PSF operator from prepare_native_annuli. It reports separately:
the signed valid-pixel aperture sum, model missing-pixel increment,
hybrid in-frame aperture flux, intrinsic Sérsic enclosed flux and
outer tail, PSF/sampling transfer (intrinsic enclosed minus full
in-frame convolved model), pure model total and hybrid corrected total.
The latter equals observed valid flux plus model total minus predicted
valid flux. No masked observational pixels are filled. Full finite PSF
footprint around in-frame aperture pixels is required; full geometric
containment of the actual contour remains unverified. Model magnitudes
are before extinction. Negative measured flux is retained.

The outer tail is integrated directly, with the regularized upper
incomplete gamma function for the constant-q outer domain. Conditional
component covariance uses the Gauss-Newton response of the parameters
to the training measurements and propagates the full ring covariance
through the total data-plus-model derivative. This retains the cross
term and does not add a separate independent parameter variance.
Neither the model total nor a plateau check establishes an empirical
asymptote or constrains arbitrary additional light outside the fit.
Optional workflow asymptotic_settings (including predeclared
aperture_count and optionally all three plateau settings) includes
asymptotic.json and full operators in arrays.npz. No automatic adoption.

`pixel_injection.py` creates a separate NativeBandFrame with an explicit
positive native-sampled Sérsic signal convolved with the empirical PSF.
It retains original invalid pixels, native masks, IVAR and NEXP, records
parent and model identities, and does not simulate additional noise.
Original FITS and in-memory observations stay unchanged. Annular records
carry the synthetic parent identity and the observational exporter rejects
them. The differential wrapper re-runs both the original and injected
workflow with matched hypotheses, checks identical measurement support,
refits sky/halo and covariance, and reports expected valid ring signal,
measured ring-flux change, biases and fitted model-total response.
Geometry and compact-source amplitudes remain fixed. This is a
differential stress test, not independent recovery of a component or
a calibrated coverage/false-positive experiment. The wrapper returns
a diagnostic result only; it does not write an observational archive.

Added tests in test_asymptotic_diagnostic.py, test_pixel_injection.py
and extended end-to-end native workflow tests. New tests still require
exact-commit CI. No additional real NGC2146sat measurements have been
made in this follow-up. Independent control positioning, replacement
injections and full noise/geometry/source refits remain outstanding,
along with real-FITS replay, full publication figures and GUI integration.


## Six-figure diagnostic export (follow-up)

`publication_products.py` provides an opt-in native-grid figure set.
Call run_reviewed_annular_workflow with publication=True and an explicit
new output ZIP path. Matplotlib is the existing optional demo extra.
The workflow passes the exact frames, reviewed scenes, fixed corrections,
PSFs, annular records, fitted parameters and covariance to the exporter.

Before rendering, publication_arrays independently reconstructs the
common validity, native annulus labels, corrected data and the finite
PSF-convolved model using the fitted radial interpolation. It verifies
native input/grid/scene identities, correction and PSF checksums, and
both observed and model ring-sum closure against the fit. A mismatch
blocks the export. No masked observational pixel is reconstructed.
The displayed model is masked to the same observed support as the data.

Six PNG/SVG pairs:
- 01_stellar_halo: before halo removal, halo model, corrected image;
- 02_prepared_images: corrected native images and geometry contours;
- 03_curve_of_growth: signed measured cumulative flux and same-support model;
- 04_surface_brightness: positive eligible measurements, before extinction;
- 05_model_residual: observed, PSF-convolved model, signed residual;
- 06_radial_model: signed mean radial intensity and raw model residual.

Native image axes are pixels, not claimed north-up sky coordinates.
Nonpositive flux stays in numerical tables but cannot be expressed as
a finite magnitude. Positive bins whose conditional lower flux interval
crosses zero are marked separately, without fabricated symmetric
magnitude errors. All error bars are conditional and all figures are
marked DIAGNOSTIC_NOT_ADOPTED. Raw residuals are not the correlated
conditional held-out residuals.

publication_arrays.npz contains the actual plotted rasters, common mask,
radius/annulus index, signed ring fluxes, conditional covariance and
support metadata. publication_qa.json records closure errors and
provenance. All members are covered by the existing ZIP manifest.
The full contour containment and empirical uncertainty calibration
remain unverified. These are publication-format diagnostics, not
scientifically approved final article figures.

Extended test_native_workflow.py checks all twelve image members when
Matplotlib is installed, manifest hashes, signed residual identity,
masked pixels and rejection of a substituted correction. CI installs
the existing exact,demo,test extras. New code has not yet passed
exact-commit CI and no real FITS publication figures were rendered
or visually inspected in this follow-up. Prior run 46146235229 for
742e2829 was canceled before completing. GUI integration and real
NGC2146sat replay remain separate open gates.


## Photometry Workspace: archive review acceptance plan

Source review: existing PhotometryWorkspace has four functional stages,
and stage/band navigation is persisted elsewhere. Introduce an explicit
Joint QA dialog from the Apertures stage, not another automatic analysis
or an incompatible navigation key. Use the existing project UUID,
selected release/region/brick and target position as the review context.

Acceptance cases before UI approval:
- Unbound/no target/no dataset: the action is disabled with a clear reason.
- Explicit local ZIP selection: a bounded read-only worker validates the
  complete diagnostic manifest, member hashes, dataset/anchor identity
  and signed ring tables before enabling review. Opening the dialog or
  canceling a file chooser performs no science operation.
- Missing, malformed, incompatible, changed and injected evidence must
  not appear as a successful review. Retain explicit recovery actions.
- Show band/radius/observed/model/conditional sigma and coverage without
  converting signed flux to magnitudes. Existing verified PNG figures
  can be inspected; absence of figures must be explicit.
- A project-bound review receipt is an explicit new JSON file selected
  by the researcher. It binds the archive digest, UUID, target and source
  selection; reopening verifies them again. It is not saved science,
  automatic adoption, or a checkpoint of a running numerical calculation.
- Empty, loading, populated, error, and resized states need offline
  regression coverage. Actual screenshot capture and inspection, Windows
  high-font/DPI and keyboard review remain separate GUI gates.

The archive schema does not currently encode a project UUID, and a
matching brick or nearby anchor does not prove galaxy identity or that
original FITS are still unchanged. Display this limitation; an explicit
researcher association is necessary before creating a receipt. This
increment intentionally does not launch the new multiband optimizer,
edit saved masks or automatically convert earlier PSF/Scene drafts.


### Implemented archive review increment

The existing fourth Photometry Workspace tab (Apertures) now provides
**Joint photometry QA…**. No fifth navigation key or new project
directory is introduced. The action is disabled without saved target
and dataset; on click the project passport is reread and UUID/source
identity rechecked. No FITS, network or science worker is opened.

`joint_archive_review.py` validates the supported ZIP schema, member
set and sizes, all manifest SHA-256 entries, exact selected source,
band order, target-anchor proximity within the existing project target
radius, signed ring table and fit/optional asymptotic/publication QA
provenance. It does not extract ZIP members or load NPZ arrays. PNGs
have size limits. These are integrity and association checks, not
authentication, numerical reanalysis, original FITS verification or
scientific validation.

`joint_archive_dialog.py` displays a loading/error status, provenance
summary, signed ring table and archived PNG selector. Existing figures
are viewed at fit-to-viewport or native display size. A read-only worker
is used for archive validation; new selection invalidates prior review.
A separate explicit confirmation enables a new JSON association receipt.
The receipt stores the exact archive SHA-256, absolute path and current
project context. Reopen is explicit, checks the current context and
revalidates the archive; moved/changed archives require a new review.
No project passport or adopted result is changed.

Tests: tests/test_joint_archive_review.py checks archive-format fixtures,
negative signed flux, association/receipt safeguards, corruption and Qt
empty/loading/populated/error/narrow states. test_native_workflow.py
also reads the actual exporter result from the synthetic native workflow.
scripts/capture_gui_baseline.py now captures four additional screenshots
with an explicitly synthetic archive-format fixture and action bounds.
No real research-object images are included in these fixtures.

Evidence at implementation time is E0 (source review). Exact new-head
CI and screenshot capture/inspection are pending; no visual approval
or on-device claim. The prior 230a07a CI was not observed completing
successfully before these edits. This increment is an archive review
interface, not yet a GUI launcher for the new numerical pipeline.
Next: integrate reviewed native Scene/PSF input snapshots, explicit
parameter setup and staged cancellable worker runs. Full resumable
computation state, original FITS verification, and real NGC2146sat
replay remain open. No changes to the numerical photometry algorithms
or science-adoption status in this GUI increment.
