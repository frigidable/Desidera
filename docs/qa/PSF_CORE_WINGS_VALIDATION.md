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


## Real Leo P g P16B result and P16C prospective test (2026-11-09)

The user's first real native P16B QA has algorithm fingerprint
\`46d99bfee698f53fdea02c13526117a683ee19e155828c69a1eb15c5f9544501\`.
The unchanged physical-finite PSF-kernel SHA256 is
\`8d976b91b4760fc8c51f607ccaf1cc59d4cac622c4c9992e5df11044bed2852a\`.

- **10/13** preselected Gaia-star profiles reached both conditional sky
  rings, including **3/3** fixed holdouts. Three training stars,
  \`star-0011\`, \`star-0036\`, and \`star-0097\`, were correctly withheld
  for absent complete native 30-FWHM cutouts; no rectangular
  clipping/imaginary ring coverage was accepted.
- All ten samples are explicitly sky-degenerate; the native-pixel
  catalog-neighbor core-mask counts range **28–72** per broad stamp.
  This masks only selected-size object cores, not their scattered
  wings or uncatalogued emission. Some original profile rings have
  severely incomplete azimuthal coverage.
- The brightest available training \`star-0002\` has an original
  native-pixel peak approximately **23.317 nMgy/native pixel** and
  conditional-outer-sky-subtracted median brightness near
  **0.0037 nMgy/native pixel** at 4–6 measured FWHM. This is
  descriptive original-pixel light in the 20.55–30.82 px ring,
  **not** a demonstrated unblended, sky-calibrated stellar PSF wing.
- The faint holdout \`star-0110\` (peak approximately
  **0.1454 nMgy/native pixel**) has a large inner-versus-outer
  sky-hypothesis difference **+0.0013014 of peak**; its conditional
  radial profiles change sign repeatedly. Do not include its
  profile in a claimed calibrated outer PSF.
- Frozen final ePSF has not changed. All six physical accuracy/
  wing/background/noise/spatial/injection stages remain HOLD.

### P16C: predeclared distinct bright Gaia-star wing-probe cohort

\`audit_disjoint_bright_gaia_wing_probes\` reads the same
**original** DR11 native FITS and **same-brick** versioned Tractor
Gaia metadata and generates a *new catalog-only cohort* that does
not contain any of the original ePSF training/holdout star
identities or native centroid matches (≤2.5 pix). The
selection is **before reading any new profile or native core-peak
pixels**; it NEVER uses P16B faint-star residual quality, an
outer slope, the desired science outcome or holdout performance.

Predeclared catalog/geometry requirements:

- \`REF_CAT=G3\` with positive unique \`REF_ID\`, Tractor \`TYPE=PSF\`,
  non-forced/non-frozen fitted model, valid maskbits (including
  band-specific saturation and ALLMASK), no Gaia duplication,
  strictly positive native-band Tractor \`FLUX\` and \`FLUX_IVAR\`,
  and diagonal Tractor S/N ≥20 (NOT covariance-calibrated).
- A *complete* native **30-measured-FWHM** stamp inside the brick;
  the source center must lie outside the fixed protected galaxy
  mask, and not have another catalog centroid within the fixed
  **3-measured-FWHM neighbor-core-mask radius**. Reject
  ambiguous centroids within 2.5 pix. These restrictions are
  initial engineering sample choices, not astrophysically
  calibrated cleanliness guarantees.
- Sort the qualifying cohort **only by original-band Tractor
  fitted total flux descending**, with deterministic \`REF_ID\`/
  \`OBJID\` tie breaks; choose the **first 16**, with no later
  backfill/replacement, regardless of native peak, original
  wing/sky profile, saturation missed by catalog, or
  unmodelled extended contamination.
- After fixing all candidate identities and their catalog-flux
  ranking, measure **original** native-pixel 3×3 local peak for
  descriptive normalization. If no positive valid peak exists
  in the original coadd, fail that candidate; **do not** pick
  a different one. Profile the accepted native pixels with
  precisely the same radial geometry, original galaxy mask,
  catalog-neighbor core exclusions and sky hypotheses as P16B.

The resulting science-HOLD QA
\`diagnostic_metrics.disjoint_bright_gaia_outer_profile_qa\`
has schema
\`desidera.psf0-disjoint-bright-tractor-gaia-native-profile/1\`.
It reports **all** catalog rejection reason counts and the
prospective inventory. This cohort is **identity-disjoint**
from the existing ePSF fit; it is *not* an independent
coadd, seeing realization, star/sky covariance calibration
or a proven clean wide-angle PSF. No per-star measurement
changes the sample size, training, holdouts, finite kernel,
existing P16B measures, science adoption, original WCS,
image pixels or inverse variance.

The real **P16C** Leo P `g` experiment was executed on
2026-11-09: 5363 original Tractor rows, 100 catalog/geometry-eligible
Gaia candidates before the fixed cap, 16 fixed identities, 16 measured
native outer profiles. All 16 share the DR11 coadd's calibration and sky.
The new identities do not validate wide wings. The faintest stars and
neighbor-rich fields remain background/contamination limited.
Compare its observed native outer radial profiles and
two **conditional** sky solutions before any PSF core+wing
stitching. If these brighter sources are contaminated or
insufficient, STOP: obtaining a separate adequately
isolated optical field or independent bright-star
calibration will be required. No new empirical PSF can
be scientifically adopted without subsequent independent
validation and quantitative LSB injection/recovery.


## P16D: one bounded, preregistered neighbor and sky stress test

**Purpose:** decide whether the same fixed 16 P16C Gaia sources support
a stable outer profile in this one DR11 coadd, or whether the Leo P
field must be set aside in favor of an independently chosen field.
This is the intended one bounded P16D experiment, NOT a sequence of
retuned deblending or selection attempts.

P16D is read-only and attaches to the **existing** PSF-0 QA as
\`diagnostic_metrics.fixed_bright_gaia_neighbor_sky_stress_qa\`, schema
\`desidera.psf0-fixed-bright-gaia-neighbor-sky-stress/1\`.

It keeps P16C's native catalog hash, all fixed Gaia IDs and order,
native coadd images/invvar, original protected target mask, unchanged
81×81 ePSF, FWHM measured on that kernel, train/holdout split, and
scientific HOLD. A missing/mismatched catalog, core geometry, stamp,
or original P16C result fails closed, rather than constructing a
replacement star list.

Predeclared nuisance-model comparison (fixed before viewing P16D results):

1. Original catalog-neighbor core exclusion at **3 measured FWHM**;
   conservative stress exclusion at **6 measured FWHM**. Both use the
   same original Tractor centroids and identical native-pixel masks,
   apart from the fixed radius. These masks deliberately do *not* model
   the actual wings of unmasked neighbors or uncataloged sources.
2. The original non-overlapping constant background hypotheses
   \`22–26 FWHM\` and \`26–30 FWHM\`. Additionally fit a robust
   first-order spatial plane using the **pooled 22–30 FWHM pixels**
   (three iterations with 3-MAD clipping); this is a conditional
   alternative, not a calibrated sky truth. It can absorb real wings.
3. Measure original-pixel medians in the same six outer bins,
   \`4–6\`, \`6–8\`, \`8–12\`, \`12–16\`, \`16–20\`, \`20–22\`
   measured FWHM. Record usable-area fractions, both constant-sky
   residuals, plane-subtracted residuals and normalized differences
   between the two masking radii.
4. Split each ring into eight fixed azimuthal sectors; only report
   a sector median with at least **25 unmasked original pixels**.
   Missing sectors stay null. Neither sector-to-sector spread nor a
   positive ring median is a covariance-calibrated detection.
5. **Internal invariance check:** at exactly 3 FWHM, P16D must
   reproduce the original P16C sky-ring medians and each matching
   radial-bin residual within a fixed absolute tolerance \`1e-10\`
   in its native units. If not, label the source
   \`BASELINE_P16C_REPRODUCTION_MISMATCH\`; do not interpret the
   stress difference. No reranking or backfilling is permitted.

**What P16D does NOT establish:** a wing detection at 8–30 FWHM, a
wide-PSF total flux, noise covariance, the optimal neighbor mask,
photometric uncertainty bounds, a science-ready PSF, or a distance.
The fixed mask radii are sensitivity hypotheses, not choices to be
optimized on the same stars.

**Decision after the real Leo P run:**

- If both masking radii leave usable, reasonably well-populated annuli
  and the conditional sky alternatives yield comparable *outer*
  profiles without large octant anisotropy, advance to a separately
  defined core/wing modeling and independent-star validation experiment.
- If outer ring residuals are comparable to mask/sky perturbations,
  or most annuli become poorly covered, **do not coadd the 16 profiles
  or fit a wide PSF**. Prefer another, prospectively chosen DESI field
  with better isolated bright stars, or independently calibrated
  exposure-level/bright-star information.
- No numerical scientific acceptance thresholds may be tuned after
  inspecting this pilot; a subsequent adoption must use a declared
  galaxy-photometry error budget and prospective injection/recovery.

For the local Leo P rerun, preserve the original options
\`gaia_strict\`, detection FWHM \`6.57 px\`,
\`outer_sky_isolation_footprint\`. Confirm the new schema exists
alongside unchanged \`kernel_sha256\`, \`n_train=10\`, \`n_holdout=3\`,
and \`science_adopted=false\`. Keep the full JSON as the experiment
record. The stage remains \`SCIENCE_GATE_HOLD\`.
