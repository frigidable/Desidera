# Desidera PSF Core & Wings Validation — fixed completion contract

Status: **active validation track, NOT_ADOPTED**. This document defines a
finite, auditable exit from PSF-0 research; it supersedes the idea of
adding an indefinite sequence of P15 neighbor-deblending diagnostics.
It does **not** grant science use to the empirical PSF.

Existing repository: no additional branch or folder. Initial work is
on the existing draft PR #43. Local originals stay under the user's
existing per-object Desidera data project.

## Scientific target and definition of “complete”

Build a **band-, field-, observation- and use-case-specific** normalized,
positive, reproducible coadd PSF with an uncertainty/bias budget small
enough for the planned LSB galaxy surface-photometry measurements.
This is not a promise of a mathematically “perfect” PSF, a universal
PSF across all DECam exposures, or an automatic SBF-quality PSF.

There are only **four milestones**. Later code/diagnostics must directly
satisfy one of these unmet exit conditions. Diagnostic improvements
alone are not a reason to rename or extend P15 indefinitely.

| Milestone | Required scientific evidence | Exit decision |
| --- | --- | --- |
| 1. Native core | Accurate, reproducible FWHM and encircled-energy core on the final, physically normalized kernel; brightness-/position-matched independent stars with a predeclared core-shape/flux error budget | HOLD until independent core residuals pass that budget |
| 2. Wings and sky | Original-pixel broad-radius bright-star wings measured beyond the core model support (radius chosen from the **galaxy measurement aperture**, not a fixed convenient threshold), separated from scattered light and background trends; finite-support loss and stitching uncertainty quantified | HOLD until wing+sky model agrees with independent stars and has bounded additive flux error |
| 3. Spatial/noise representativity | Independent stars distributed across the relevant coadd footprint and magnitude range; comparison to survey PSF maps and spatial PSF variations; resampling-induced coadd covariance measured or validated against blank-sky statistics | HOLD until position/brightness and covariance tests satisfy the predeclared error budget |
| 4. Galaxy injection-recovery | Known-input LSB galaxy image models injected into representative **original** survey-image backgrounds; full photometry pipeline rerun blind to injected truth; recovered surface-brightness profiles, sizes, ellipticity and integrated flux assessed for systematic bias and coverage | HOLD until **externally chosen, science-driven** bias/uncertainty tolerances pass on prospective tests |

**No numbers will be retroactively tuned to Leo P's three holdout stars.**
Core/wing/bias tolerances must be written **before** selecting recipes
or rerunning reference experiments, justified by a specified scientific
observable (surface-brightness profile to a declared radius, faint
diffuse flux, integrated magnitude, or later SBF). They must then
be evaluated on a separate, prospective sample. If the science
requires stricter accuracy than these data can constrain, status
remains HOLD and the pipeline must disclose the limitation.

## First deliverable: P16 Core & Wings validation preflight

The currently implemented PSF kernel remains the **original**
81×81 native-pixel ePSF. No source selection, PSF building, clipping,
Tractor fitting, pixel data, or downstream photometry is modified.

Added pure diagnostic module:

\`src/lsbpipe/photometry/psf_core_wings_validation.py\`.

The new audit is serialized in the existing:

\`diagnostic_metrics.core_wings_validation_qa\`

with method schema
\`desidera.psf0-core-wings-validation/1\`.
The complete QA and existing GUI summary report the following:

* FWHM along the **two native pixel axes**, measured by linear
  interpolation of the first half-maximum crossings on each side
  of the observed **final-kernel peak**. Their geometric mean is
  \`measured_core_fwhm_axis_geometric_mean_native_pix\`. This is
  not a two-dimensional fitted FWHM, not a variance estimate, not an
  optical PSF spatial model, and **does not replace**
  \`assumed_detection_fwhm_pix_not_measured\` used for earlier
  predeclared star selection and P14/P15 tests.
* The native-peak offset from the grid center and the x/y width
  ratio. Extreme asymmetry needs separate optical, subpixel and
  centroid validation; measured axis width alone is not sufficient.
* Encircled **fraction of finite kernel** at 1, 2, 3 and 4
  **measured** FWHM, counting only complete circular apertures.
  Any aperture extending past the rectangular kernel's minimum
  edge distance returns \`null\` rather than summing a clipped
  circle. Normalizing by finite-kernel sum is **not true PSF
  total-flux calibration**, even when the number is close to 1.
* Existing, separately fit, same-three-star original-native
  holdout residual and core-normalized curve-of-growth evidence,
  explicitly labelled descriptive. Coadd covariance and the
  stellar sky-wing degeneracy have not been calibrated.
* The original pre-clipping negative absolute-mass fraction,
  **not** a negative physical probability and not the post-clipped
  PSF error. It requires a scientifically justified model
  selection/treatment; clipping alone is not a fidelity guarantee.

### Gate semantics

A \`PASS\` means the **required basic evidence exists** for a
numerical/provenance audit. It never means the PSF is scientifically
validated. \`HOLD\` means more independent scientific evidence is
required. There is no success path from this initial module to
scientific adoption.

| QA stage | Initial condition | Science decision |
| --- | --- | --- |
| Original native input provenance | Local image/invvar/nexp/maskbit SHA256 quartet, with provenance | PASS/HOLD on availability only |
| Finite kernel numerical geometry | Nonnegative unit-normalized kernel and bracketed half maxima on both native axes | PASS/HOLD on numerical measurability only |
| Predeclared holdout coverage | All **original** holdouts with all original zones scored | PASS/HOLD on coverage only |
| Core accuracy | External, predeclared, unbiased error budget and prospective stars | HOLD |
| Wide-angle stellar wings | Real large-radius empirical measurements with quantified truncation | HOLD |
| Sky/wing identifiability | Constrained sky–wing degeneracy and outer background | HOLD |
| Spatial representativity | Independent field- and brightness-matched sources, survey PSF comparison | HOLD |
| Resampled coadd covariance | Measured or validated correlated-noise model | HOLD |
| LSB galaxy injection-recovery | Predeclared bounds on profile, structure and integrated flux bias | HOLD |

Global status must remain \`SCIENCE_GATE_HOLD\` and
\`science_adopted=false\` until a **separate reviewed science
adoption implementation** has demonstrated every required milestone.
The existing \`PSF0Result\` cannot produce an adopted PSF result.
P15-type star/galaxy/sky nuisance fits remain QA-only and are
*not* promoted to wide-angle PSF measurements.

## Baseline and immediate test

User-provided original Leo P dataset
\`dr11/south/1554p180/g\`:

* Frozen ePSF SHA256
  \`8d976b91b4760fc8c51f607ccaf1cc59d4cac622c4c9992e5df11044bed2852a\`.
* 10 preselected training stars and 3 fixed Gaia holdouts.
* Finder FWHM **6.57 px assumed**, still **not** a
  scientific native-kernel width measurement.
* Three final-kernel holdout annuli assessed; at 4×**guessed**
  FWHM the median absolute core-normalized growth difference
  is **0.095869648**, which is NOT automatically 9.6% bias in
  galaxy surface brightness.
* Pre-clipping negative absolute mass fraction
  **0.016537987**, concentrated outside the bright core.
* Previous P15G observed native astrophysical deblending
  is incomplete; the real outer PSF wings, noise covariance
  and spatial representativity are not validated.

First prospective run should produce
\`core_wings_validation_qa\` with measured core geometry.
The **existing kernel SHA256, train/holdout split, source
FITS hashes and science_not_adopted status must match** the
previous same-recipe Leo P run. Difference means an input/lineage
change, not a free "improvement" to adopt.

## Stop rule and future work

After the first real core-width measurement, priority is to
study *independent bright stars with genuinely larger native
cutouts* (not simply extrapolate the 81×81 kernel). Distinguish
PSF structure from coadd sky; then bound spatial/coadd covariance
and perform galaxy injection-recovery. Work on source deblending
only if it removes a **demonstrated blocker** of the above
tests; there is no requirement to make P15 perfectly reproduce
Tractor before measuring native PSF wings.

If all milestones pass to a reviewed band/use-case-specific
error budget, a separately reviewed production/adoption path
may be implemented. If not, retain a scientifically defensible
core-only experimental PSF and report the limits honestly.


## First original Leo P g P16 result (2026-11-09)

From the first real original-FITS output,
\`dr11/south/1554p180/g\`, P16 method
\`627222e7f3f5193ffb10db274112ef4201d58ca50b7e493fcc4ba6b3c2ad361a\`:

- The exact original 81×81 kernel hash remains
  \`8d976b91b4760fc8c51f607ccaf1cc59d4cac622c4c9992e5df11044bed2852a\`.
  The original sample remains 10 training/3 fixed Gaia holdouts.
- Its **peak-axis half-maximum widths** are **4.963218 px**
  (native X) and **5.316666 px** (native Y); geometric mean
  **5.136903 native pix**, compared with a **6.57 native pix
  detection FWHM guess**. The measured width is **21.81% lower**
  than the guess. The axis max/min ratio is **1.071213**.
  The native kernel peak is at exact 81×81 grid center (40,40).
- Fractions of the 81×81 **finite-normalized** kernel within
  1, 2, 3, 4 **measured** FWHM circular radii are
  **0.8035703, 0.9403957, 0.9720075, 0.9831932**.
  These are *not* true wide-angle PSF encircled-energy
  fractions; there may be substantial flux outside the 81×81
  support and/or absorbed into sky fitting.
- On three fixed independent holdouts the previous
  4×**guessed** FWHM median absolute core-normalized
  growth mismatch is **0.095869648**; this uses a different
  aperture and normalization and is NOT a calibrated
  galaxy-photometry flux bias.
- Exactly three numerical/provenance **availability** stages
  PASS; six scientific accuracy/wing/noise/spatial/injection
  stages HOLD. The global gate is
  \`SCIENCE_GATE_HOLD\`, not adopted. The native negative
  absolute ePSF precursor fraction remains **0.016537987**.

The wide-angle PSF is therefore **not validated** by the
finite-kernel encircled fraction or native FWHM.

## P16B — original-image outer-profile feasibility, before any PSF wing fit

A separate **read-only** follow-up
\`src/lsbpipe/photometry/psf_wing_profile.py\` uses **original
native coadd pixels** and the same already frozen 10+3 stars;
it does NOT retrain, change or accept any PSF. For each selected
source, use an original large odd-pixel native stamp that fully
contains the disjoint sky-hypothesis rings through **30 measured
native-kernel FWHM**. A star outside this entire footprint is
marked unavailable, never silently truncated.

Require the verified same-brick native Tractor catalog and
exactly one counterpart within 2.5 native pixels. Mask
the **cores** of every other original-catalog position inside
the patch at a fixed **3 measured-FWHM radius**, and
also preserve the already predeclared **Leo P target-exclusion
mask** on the whole larger stamp and both sky annuli. This
does NOT remove wings of catalogued neighboring stars,
unmodelled broad halo/deblending errors, diffraction,
background galaxies or scattered light. Original valid,
inverse-variance and SATUR/ALLMASK masks are also retained.

Record original-pixel robust annular medians and eight-sector
coverage for the immutable annular edges **1,2,3,4,6,8,
12,16,20,22** measured-FWHM units. Two **non-overlapping**
outer sky hypotheses use measured-FWHM rings **22–26**
and **26–30**, respectively. Compare the same annular
original-pixel median after subtracting either *conditional*
constant sky and dividing by the star's original measured
peak. These sky estimates are deliberately **not** claimed
to be free of the target star's outer PSF wings or of
spatial sky gradients; differences are a sky/wing
degeneracy diagnostic only.

The original catalog mask can leave incomplete azimuthal
coverage and cannot remove other source wings. Per-annulus
usable pixel counts, geometric coverage and octant counts
must be retained, with an explicit cannot-measure status
for missing data/counterpart/sky. The program must **never**
use the outer-profile evidence to alter original training,
holdout roles, current PSF kernel, assumed finder FWHM,
legacy fixed-ring diagnostics or science-adoption status.

QA serialization:
\`diagnostic_metrics.original_native_outer_stellar_profile_qa\`
with schema
\`desidera.psf0-original-native-wide-wing-profile/1\`.
Expected success status
\`MASKED_OUTER_PROFILE_DESCRIPTIVE_SKY_DEGENERATE\`
means *raw data can be inspected*, not stellar wings have
been measured without bias. This is not an independent
bright-star sample, a certified sky model, a calibrated
PSF wing, or a science pass. The initial scientific
stage \`wide_angle_wing_measurement\` remains **HOLD**.

**Next** after checking the Leo P P16B real output: if
selected stars are confused or limited by spatial
coverage, design a *separately fixed* prospective
bright-star wing sample. A PSF core and true wide-angle
profile cannot be inferred by choosing only the
best-looking posthoc wing residuals from the same
ten training and three holdout stars.
