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
| Core accuracy | Predeclared core error budget with within-brick heldout stars and injections | HOLD |
| Wide-angle stellar wings | Real large-radius empirical measurements with quantified truncation | HOLD |
| Sky/wing identifiability | Constrained sky–wing degeneracy and outer background | HOLD |
| Spatial representativity | Position-separated, brightness-matched same-brick source controls; global PSF scope limited to this brick | HOLD |
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
insufficient, STOP the physical wing inference and preserve
the measurable conditional profile bounds for THIS brick.
An external field is not an input to this brick-specific PSF. No new empirical PSF can
be scientifically adopted without subsequent independent
validation and quantitative LSB injection/recovery.


## P16D: one bounded, preregistered neighbor and sky stress test

**Purpose:** decide whether the same fixed 16 P16C Gaia sources support
a stable outer profile in this one DR11 coadd, or whether the
outer-wing model must remain unresolved within this brick.
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


## P16D original Leo P g result — 2026-11-09

Input: the user-provided unmodified scientific QA export for
\`dr11/south/1554p180/g\`; method SHA256
\`e984d2560da95c0dd8409d4a59d9723e790e441d6f86a531a3f5bf410307a5af\`.
Its original ePSF SHA256 remains
\`8d976b91b4760fc8c51f607ccaf1cc59d4cac622c4c9992e5df11044bed2852a\`,
81×81 pixels and measured peak-axis geometric-mean FWHM 5.13690296 px;
10 training and 3 original independent-identity holdout stars, unchanged.

- 16/16 P16C catalog-selected stars and sky/radial profiles exactly
  reproduced in the 3-FWHM P16D baseline; 16/16 have both conditional
  sky-ring values and conditional planar backgrounds for both radii.
  Those are **computational completeness measures**, not scientific
  wing detections or adequate unmasked ring coverage.
- Median available annular area across these fixed 16 stars (3-FWHM
  masks -> 6-FWHM masks):

  | Ring (measured FWHM) | 3-FWHM masked | 6-FWHM masked |
  | --- | ---: | ---: |
  | 4–6 | 85.7% | 40.1% |
  | 6–8 | 72.9% | 27.7% |
  | 8–12 | 67.9% | 25.8% |
  | 12–16 | 72.7% | 20.3% |
  | 16–20 | 70.2% | 23.1% |
  | 20–22 | 66.1% | 18.8% |

- Illustrative *post-result diagnostic*, not a preregistered adoption
  gate: requiring >=50% unmasked area AND >=6 populated octants under
  the 6-FWHM neighbor mask gives 4/16 usable samples at 4–6 FWHM,
  1/16 at 6–8 FWHM, and 0/16 at every measured ring from 8–22 FWHM.
- At 4–6 FWHM 14/16 of the 6-FWHM-mask ring medians are positive
  against the outer 26–30 FWHM constant-sky hypothesis; 13/16
  positive after subtracting the alternative outer-sky plane.
  This does not establish wing detection. In the 8–12 FWHM ring,
  just 10/16 stars preserve their ring-residual sign between 3-FWHM
  and 6-FWHM masks, and the median absolute mask-induced residual
  change is 1.72e-5 of the native peak, compared with a median
  6-FWHM-mask ring residual of -8.42e-6 of peak. The median
  absolute change induced by plane-versus-constant-sky subtraction
  in that ring is 4.51e-6 of peak.
- \`wing-cat-0014\` and \`wing-cat-0015\` retain just 18 and 7 native
  pixels, respectively, in their 6-FWHM-mask 4–6 FWHM rings.
  Both have 0 populated octants. In particular, the large positive
  4–6 FWHM residual of \`wing-cat-0015\` must not be interpreted
  as wing evidence.

**Scientific decision at P16D:** STOP wide-angle PSF estimation from the
16-star Leo P cohort under the 6-FWHM *stress mask* with the available coadd alone. The stress test
demonstrates that the nominal profiles cannot reliably constrain
the **outer** wings with conservative neighbor masking; the two
background hypotheses and the plane remain degenerate with possible
real stellar flux, and covariance has not been calibrated. This does
not prove the true stellar wings are absent, nor does it invalidate
the native coadd or justify adopting the unvalidated 81×81 core.
Keep \`science_adopted=false\`, \`SCIENCE_GATE_HOLD\`; preserve the
unchanged original P16C cohort and the P16D report as a documented
negative feasibility pilot.

**Superseded planning note — fixed Leo P brick only:** use the
existing catalog-qualified stars, native coadd pixels and same-brick
blank-sky controls. If the stellar population cannot support calibrated
wings, report an unresolved wing component; do not silently source
stars or background from a different brick. Define galaxy-aperture-dependent
outer radial extent and quantitative predeclared stability/uncertainty
requirements; compare field-wise wing profile measurements to the
existing Leo P data only after the new selection is fixed. Do not
relax Leo P's isolation cuts post hoc, average its 16 incomplete
profiles, or manufacture a wide-angle core+wing PSF.

**Additional caution:** a 6-FWHM circular neighbor exclusion can
cut into a target's own intermediate radial annulus when another
Tractor source lies within 6 FWHM of the target; its severity is
a feature of this deliberately conservative stress test, not a
physically calibrated optimal masking policy.


## P17: one preregistered 48-star Leo P mask sensitivity experiment

P16D's 6-FWHM mask is a deliberately pessimistic geometric test.
Its severe annular area loss **does not** establish that a
3-FWHM exclusion fails or that Leo P contains too few stars.
The original P16C catalog screening found 100 Gaia/Tractor candidates
before the 16-star *quota*. Therefore, repeat this one **bounded**
feasibility experiment in Leo P before requesting a new DR11 field.
This does not undo P16D's refusal to adopt or extrapolate PSF wings.

The immutable catalog-only ranking remains the exact P16C sequence:
Gaia G3, Tractor PSF, valid native-band flux and diagonal-only
nominal SNR >= 20, clean mask/fit/duplicated-source flags, complete
30-FWHM native cutout, unambiguous catalog core, at least 3 measured
FWHM from the nearest OTHER Tractor centroid, outside the protected
target, and excluding the original 10 training + 3 holdout identities.
No profile signal, image peak, annular completeness, or sky-model fit
can change eligibility or the original flux ordering.

The fixed cap is **48 catalog identities**. Identity roles are
assigned strictly by this rank, before reading new native profiles:

- Ranks 1–16: **legacy P16C** (previously inspected, never validation).
- Ranks 17–48 not divisible by four: **24 new exploratory stars**.
- Ranks 20,24,28,32,36,40,44,48: **8 new identity-heldout stars**.

A missing/invalid profile stays in the rank list and is *never*
replaced. All 16 legacy Gaia IDs, Tractor IDs, native coordinates
and original ranking must match the P16C report exactly. P17 must
also reproduce the legacy 3-FWHM sky and radial medians at absolute
1e-10 tolerance; mismatch is marked, never fit away. A new
identity-heldout star is **not** independent of the coadd, its sky,
exposure-level PSF, or neighboring sources. All four masks are
evaluated on all stars for descriptive comparison; there is no
outcome-driven best-mask selection.

Four fixed hypotheses based on **measured** core FWHM:

1. Catalog-neighbor circular mask of 3 FWHM.
2. Same centroids with mask radius 4.5 FWHM.
3. **Catalog-flux-tiered** mask only: relative Tractor native-band
   neighbor/target flux below 0.01 -> 3 FWHM; in [0.01, 0.1)
   -> 4.5 FWHM; at least 0.1 -> 6 FWHM. A non-PSF neighbor,
   unknown, negative or nonfinite catalog flux receives 6 FWHM.
   This is a *hypothesis* rather than a validated PSF-based
   contamination calculation; extended morphology may exceed
   any selected mask radius.
4. Uniform 6-FWHM mask, retained as an intentionally conservative
   stress condition, not an automatic quality gate.

Output:
\`diagnostic_metrics.extended_bright_gaia_mask_experiment_qa\`,
schema \`desidera.psf0-fixed-48-gaia-catalog-neighbor-masks/1\`.
For every strategy record both conditional sky-ring medians
(22–26 and 26–30 measured FWHM), outer 22–30 FWHM planar
background, radial medians in the six fixed 4–22 FWHM bins,
unmasked area, and eight fixed-sector medians requiring at least
25 original pixels in an octant. The original coadd/variance/masks,
the 81×81 native ePSF, training/holdout and all science gates
are unchanged.

**Predeclared interpretation for this feasibility pilot:**
Compare the 24 new exploratory stars to the 8 new identity-heldout
stars *without tuning masks on either split*. Report sky-model and
neighbor-mask induced scatter against the actual residual amplitude,
as well as available annular coverage and populated octants. Use
the 16 old P16C stars for continuity, not for confirming a new
scientific detection. Do not infer a significance without estimating
the image/noise covariance and background controls. Neither residual
sign consistency nor a large population size substitutes for
morphological/brightness-matched independent PSF validation.

**Decision discipline:** P17 can demonstrate *feasibility* or
*infeasibility* for selected mask hypotheses, not validate or
authorize a wide-angle PSF model. If no hypothesis yields
stable profiles with adequate spatial coverage across both
new exploratory and holdout stars, retain the residual
limits as an explicit unresolved component for this brick.
If promising, the NEXT
separately declared phase must independently constrain
background/noise covariance, neighbors and spatial PSF
variation, and perform blind galaxy-injection recovery
before adopting any core+wing PSF.

This adds no new branch, data directory, alternative PSF
product, empirical wing fit, or automatic science gate.


## P17 real Leo P g result — 2026-11-09

Source: user-provided full native DR11 QA export,
\`dr11/south/1554p180/g\` (JSON text file SHA256
\`8b9821162935cf1d9b92799c72eaafae0d5fa6eff79c98a65519edf65cbe5875\`).
Implementation-method SHA256:
\`86dec8f63d050561a1aaed17612b918bccdd24184f0f10d97229ec50d63c9ca2\`.
The original 81×81 ePSF SHA256
\`8d976b91b4760fc8c51f607ccaf1cc59d4cac622c4c9992e5df11044bed2852a\`,
10 training + 3 holdout stars, and \`science_adopted=false\`
are unchanged.

**Lineage and execution: PASS**. The 100 eligible catalog-only
candidates remain available. All fixed 48 have measurable native
profiles under the P17 diagnostic, with 16 legacy, 24 new exploratory,
8 new identity-heldout stars. The 16 legacy 3-FWHM sky and radial
P16C baselines reproduce exactly within the declared tolerance.
All 48 have four recorded mask variants. Seven of the 1,152
annular mask×star×radius entries have *no* native pixels;
all seven are under the 6-FWHM stress mask. "48 measurable stars"
must **not** be reported as 48 independently adequate wing profiles.

**Unmasked area** (median across all 48; 4 mask strategies in
3, 4.5, catalog-flux-tiered, 6 measured FWHM order):

| Ring (FWHM) | 3 FWHM | 4.5 FWHM | Tiered | 6 FWHM |
| --- | ---: | ---: | ---: | ---: |
| 4–6 | 72.3% | 50.3% | 37.8% | 28.0% |
| 6–8 | 72.9% | 46.1% | 34.7% | 25.3% |
| 8–12 | 71.0% | 49.2% | 37.5% | 26.0% |
| 12–16 | 73.9% | 48.0% | 36.0% | 24.9% |
| 16–20 | 72.5% | 47.9% | 38.0% | 25.6% |
| 20–22 | 72.3% | 47.8% | 38.2% | 27.1% |

For comparison, a **post-result descriptive diagnostic** requiring
at least 50% unmasked area AND at least six octants each with at
least 25 pixels has these qualifying star counts (out of 48):

| Ring | 3 FWHM | 4.5 FWHM | Tiered | 6 FWHM |
| --- | ---: | ---: | ---: | ---: |
| 4–6 | 43 | 20 | 14 | 9 |
| 8–12 | 44 | 23 | 8 | 3 |
| 12–16 | 47 | 21 | 6 | 0 |

The corresponding **8 new heldout identities** at 8–12 FWHM
qualify 7, 5, 2, **0** times, respectively. These are
*illustrative spatial/area diagnostics* and not preregistered
statistical acceptance or proof that a 3-FWHM mask cleans
neighboring PSF wings. The catalog-flux-tiered variant is
more spatially destructive than uniform 4.5 FWHM in this
field; by design it assigns 6-FWHM masks to non-PSF
neighbors regardless of their flux. The rule is not
validated by these results.

**Source-brightness selection effect:** median catalog Tractor
\`g\` flux of the 16 legacy stars is 176.52 nMgy,
compared with 15.08 nMgy for 24 new exploratory stars and
13.39 nMgy for eight new controls. Increasing the count
does not increase the flux-limited information in faint
outer wings: the newly added probes are roughly an order
of magnitude fainter.

**Normalized residual and sky sensitivity** (3-FWHM mask,
all values relative to each star's original native peak):

- At 4–6 FWHM, 15/16 legacy stars have a positive residual
  under each of the 3 and 4.5 FWHM masks; legacy
  median residual under 3-FWHM mask: +9.24e-5, median
  absolute deviation 1.54e-5. But only 5/8 new
  heldout stars have a positive residual with the
  3-FWHM mask (median +4.16e-5; MAD 1.73e-4).
  This is **not** independent validation of a detected wing.
- At 8–12 FWHM under the 3-FWHM mask, the all-star median
  signed residual is +5.22e-6 of peak but the median
  **absolute** residual is 5.42e-5. The median absolute
  change from 3 to 4.5 FWHM is 3.70e-5,
  and 37/48 source-residual signs agree across those masks;
  from 3 to 6 FWHM, the median absolute change reaches
  9.25e-5 (over 47 sources with non-null residuals).
  The spatial coverage of the 6-FWHM results is poor.
- For the 32 NEW stars at 8–12 FWHM with 3-FWHM
  masking, the median absolute stellar residual is
  1.41e-4 of peak, while the median absolute difference
  between the two disjoint sky-ring hypotheses is
  1.18e-4 of peak, and median absolute constant-versus-plane
  offset is 4.79e-5. Sky subtraction is therefore a
  limiting nuisance at the flux level being investigated.
- At 4–6 FWHM under the 3-FWHM mask, the new 24
  exploratory stars have median normalized azimuthal
  max-minus-min octant spread of 2.86e-3 of peak,
  compared with a median signed radial residual
  +1.71e-4. The eight new controls have median octant
  spread 2.50e-3 and median signed residual +4.16e-5.
  These are angular-scatter diagnostics, **not** calibrated
  errors or formal detection significances.

**Interpretation and decision:** P17 successfully demonstrates
the *number* of available Gaia stars is not the binding catalog
constraint. The 3-FWHM mask retains usable geometry for most
of the 48 objects; the uniform 6-FWHM and the current
catalog-flux-tiered masks erase too much annular area.
However, the additional stars are much fainter,
show large sector-to-sector variation, and the 8-star
identity holdout does not provide brightness-matched,
sky/covariance-independent evidence of broad wings.
There is no validated 4–22-FWHM stellar PSF; the
original ePSF and \`SCIENCE_GATE_HOLD\` persist.

Do **not** choose a "winning" P17 mask post hoc from this
report, or extend the cohort to 100 simply to increase N.
Retain Leo P as a **diagnostic** field for developing
source-specific neighbor footprints and explicit sky
uncertainty; require a separately selected, sufficiently
bright and spatially isolated **independent** DESI
field (and ultimately individual-exposure/covariance
checks) for any scientific confirmation. Predeclare
morphology-aware contamination and recovery tests
before implementing them; do not assume a model-subtracted
neighbor is correctly decontaminated from a coadd
without injection/recovery validation.


## P18: native Tractor-shape neighbor footprint and azimuthal sky sensitivity (2026-11-09)

**Purpose:** P17 confirmed that uniform 6 measured-FWHM masks and
the first flux-tiered masks erase too much Leo P annular area,
while a 3-FWHM circular mask has not been shown to remove
extended galaxies or bright stellar halos. P18 introduces ONE
additional *predeclared* catalog-shape exclusion hypothesis,
without reselecting stars or adopting a preferred P17 method
based on the observed residuals.

Input lineage: exactly the same fixed 48 native Gaia/Tractor star
identities, brightness rank and 16/24/8 legacy/exploratory/new
identity-heldout roles as P17; exact native FITS catalog hash.
The same native original image, inverse variance, masked-pixel
geometry, measured ePSF core FWHM, target mask and sky rings.
P18 fails closed when WCS is missing/degenerate, a star's Gaia
REF_ID/Tractor OBJID/native coordinate is inconsistent with its
P17 identity or the P17 cohort cannot be confirmed. There is
NO adaptive cherry-picking, replacement or refitting of stars.

**Catalog-only source footprint:**

- For Tractor TYPE=PSF the circular exclusion radius is 3,
  3.5 or 4.5 measured FWHM, according as native Tractor
  neighbor/target flux ratio is <0.1, [0.1,1), or >=1.
  These fixed *hypotheses* are neither encircled-energy
  measurements nor physical-scattering guarantees.
- For TYPE=REX, EXP, DEV, SER, use the original release/brick
  Tractor \`SHAPE_R\` in arcsec, \`SHAPE_E1\`, \`SHAPE_E2\`
  with \(\epsilon=(a-b)/(a+b)\). The sky tangent major
  axis is \(\theta = \tfrac12\operatorname{atan2}(e_2,e_1)\)
  measured from local east toward north. The ellipse's
  semi-axes in arcsec are \(a^2=r_{\rm core}^2+(4r_e)^2\)
  and \(b^2=r_{\rm core}^2+(4qr_e)^2\) for
  \(q=(1-|\epsilon|)/(1+|\epsilon|)\).
  For REX, \(\epsilon=0\). A local WCS Jacobian
  converts these oriented tangent-plane ellipses into
  native pixel masks; no hardcoded DECam plate scale.
- Missing/invalid non-PSF radius or ellipticity, unknown
  morphological type or missing/negative band flux triggers
  an explicit **6-FWHM circular fallback** for that
  source, rather than guessing an ellipticity or
  claiming an unmasked clean region. A very large
  galaxy model may obscure an entire radial ring.
  Full fallback reason counts are recorded per star.
- The factor four, core radii and fallback are intentionally
  provisional, catalogue-derived stress geometry.
  These analytic elliptical footprints are NOT
  exact Tractor source rendering or convolved profiles;
  real galaxy tails and stellar scattering may remain.
  A footprint does not establish absence of
  light beyond its boundary.

**Sky diagnostics:** Preserve constant-sky medians in *two disjoint*
22–26 and 26–30 measured FWHM rings plus a pooled conditional
planar fit. Within the outer 26–30 ring, separately delete each
of eight fixed angular octants and recompute the masked native
median, recording all eight replacements and their maximum
absolute change divided by the *unchanged original 3×3
native stellar peak*. These eight estimates share most pixels;
they ARE NOT independent jackknife replicates, calibrated
standard deviations, error bars or a coadd covariance model.
Keep the same six 4–22 FWHM radial bins with coverage,
at least-25-pixel octant count and constant/planar residuals.
The stronger background requirement cannot create false
"sky-free" PSF wings.

Output:
\`diagnostic_metrics.catalog_shape_neighbor_sky_qa\`, schema
\`desidera.psf0-catalog-shape-neighbor-sector-sky/1\`.
Compare geometry, normalized residuals and sky sensitivity
with P17 3-FWHM baseline by the SAME fixed star identities.
The readout must separately report 16 previously inspected
legacy stars, 24 new exploratory stars and eight predeclared
identity-heldout stars. Preserve all missing-profile outcomes;
do not use a post-hoc residual criterion to promote a star,
mask or PSF model. In particular do not turn the eight
intra-coadd identity holdouts into fictitious independent
exposure/sky validation.

**Fixed-brick scientific scope:** Do not select another field.
Within `1554p180`, the measured bright-star census,
original coadd, Tractor morphology, same-brick blank sky,
spatially separated stars and within-brick injections are
all admissible sources of evidence. A high number of faint
stars alone is insufficient for wings. A calibrated wing
correction requires explicit background/noise-covariance
bounds and morphology-matched injection/recovery in this
same brick, and remains on HOLD if these are not obtained.

No new branch or directory, no changes to legacy P16B–P17
QA, no PSF kernel change, and \`science_adopted=false\`.


## P18 real Leo P g result — 2026-11-09

Source: two identical user-provided full native DR11 QA exports
(\`Вставленный текст(20261109-101744).txt\`,
\`Вставленный текст(20261109-101827).txt\`),
each 5,074,178 bytes, bytewise SHA256
\`4fa0482c10200e20d09dfdf08fda482cf4346b9be5d2224f0fe2b84b57bcbbf5\`.
Exact dataset \`dr11/south/1554p180/g\`; native
PSF kernel SHA256
\`8d976b91b4760fc8c51f607ccaf1cc59d4cac622c4c9992e5df11044bed2852a\`.
P18 implementation-method SHA256
\`bfc0343deb5b316ab9f34f6d7354a60ecdaa6deb67fb2b51aec6db2dd3078085\`.
Source catalog hash matches P17:
\`8e8e3f75af868f3ef7cf4c4c9d7ca28210cbca87ef27c98351a6ab237010562c\`.
Original core training=10, holdout=3, with
\`science_adopted=false\`, \`psf_wings_validated=false\`
and uncalibrated coadd-noise covariance.

**Execution and frozen identities PASS:** all fixed 48 P17 identities
(16 legacy, 24 new exploratory, 8 new identity-heldout) are preserved.
P18 reports 48/48 with BOTH conditional sky rings and conditional
plane; this count **does not imply** all annular measurements are
geometrically complete. Heldout rank 20 has exactly zero usable
native pixels in the 4–6 and 6–8 FWHM bins under P18 despite
both 22–30-FWHM sky rings being measurable. No star was replaced.

**Coverage comparison:** medians across the same 48 stars. The
"coverage criterion" is POST-RESULT descriptive: >=50% unmasked
annular area and >=6 azimuthal octants each containing >=25 pixels.

| Native annulus | P17 3FWHM area / qualifying N | P18 shape area / qualifying N | P17 4.5FWHM area / qualifying N | P17 tiered area / qualifying N | P17 6FWHM area / qualifying N |
| --- | ---: | ---: | ---: | ---: | ---: |
| 4–6 FWHM | 72.3% / 43 | 67.0% / 36 | 50.3% / 20 | 37.8% / 14 | 28.0% / 9 |
| 6–8 FWHM | 72.9% / 40 | 64.0% / 36 | 46.1% / 20 | 34.7% / 10 | 25.3% / 6 |
| 8–12 FWHM | 71.0% / 44 | 63.4% / 37 | 49.2% / 23 | 37.5% / 8 | 26.0% / 3 |
| 12–16 FWHM | 73.9% / 47 | 68.2% / 42 | 48.0% / 21 | 36.0% / 6 | 24.9% / 0 |
| 16–20 FWHM | 72.5% / 46 | 65.5% / 38 | 47.9% / 21 | 38.0% / 9 | 25.6% / 2 |
| 20–22 FWHM | 72.3% / 47 | 64.5% / 40 | 47.8% / 23 | 38.2% / 8 | 27.1% / 2 |

At 8–12 FWHM by the identical star identities, P18 vs P17
3-FWHM geometrically qualifies **13 vs 15 /16 legacy**,
**18 vs 22 /24 newly exploratory**, and
**6 vs 7 /8 newly identity-heldout**. P18 conserves
more useful area than the old blanket 4.5- or 6-FWHM
assumptions but less than the P17 3-FWHM circles.
Mask footprints alone cannot establish the purity of
the remaining pixels or that the narrower strategy is superior.

**8–12 FWHM normalized stellar residuals:** across 48 stars
the P18 median signed Sky-B-subtracted residual is
\`+2.12e-6\` of native stellar peak, while the median
*absolute* residual is \`4.45e-5\`. Under P17 fixed
3-FWHM the corresponding values are \`+5.22e-6\`
and \`5.42e-5\`. For the same stars, the median
absolute P18-versus-P17 residual shift is \`2.43e-5\`;
39/48 preserve the residual sign. These values show
a *mask-model-dependent weak signal*, NOT a wing
detection or measured PSF profile. Comparing
absolute-residual medians is not a quality metric:
a smaller residual may reflect extra sky absorption.

**Correlated background diagnostics:** median absolute
P18 Sky-A versus Sky-B background difference, and
median maximum absolute leave-one-octant Sky-B
change, each normalized by native stellar peak:

| Predefined fixed identity cohort | Sky-A vs Sky-B | Sky-B sector deletion maximum | Median absolute 8–12FWHM P18 residual |
| --- | ---: | ---: | ---: |
| 16 legacy | 8.39e-6 | 4.46e-6 | 6.78e-6 |
| 24 new exploratory | 8.81e-5 | 6.54e-5 | 7.80e-5 |
| 8 new identity-heldout | 1.21e-4 | 7.00e-5 | 2.33e-4 |
| All 48 | 2.49e-5 | 4.49e-5 | 4.45e-5 |

The eight leave-one-octant medians overlap strongly:
their maximum deviation is **not** a calibrated standard
error, variance, jackknife covariance, or independent
statistical test. The Sky-A/Sky-B discrepancy is
also a *systematic-sensitivity* comparison rather
than a unique sky truth. For the NEW 24 exploratory
stars, the median sky-sensitivity is of the
same order as the median absolute measured residual.
For the eight new intra-coadd identity holdouts,
median angular max-minus-min sector variation at
8–12 FWHM is \`1.77e-3\` of peak, versus the median
signed radial signal \`+1.43e-4\`, so there is no
brightness-matched spatially robust outer wing detection.

**Bright legacy inner annulus:** P18's 4–6-FWHM
Sky-B residual is positive for 16/16 legacy stars
(median \`+8.81e-5\` of native peak), versus 15/16
under P17 fixed 3-FWHM (median \`+9.24e-5\`).
However the median P18 legacy angular sector spread
there is \`1.60e-4\`, exceeding the median radial
signal, and the old legacy cohort has already been
inspected. A positive median or 16/16 signs is
**not a blind validation**, given unresolved
foreground/background and neighbor contamination.

**Source-footprint and reproducibility check:** P18 evaluated
the Tractor neighbor \`PSF\` and extended
\`REX/EXP/DEV/SER\` shapes from the exact native DR11 brick,
without altering the core PSF. The DR11 catalog documentation
defines \`shape_e1/e2\` as
\(\epsilon=(a-b)/(a+b)\exp(2i\phi)\) and
\(q=(1-|\epsilon|)/(1+|\epsilon|)\), agreeing
with the scalar conversion implemented here.
This does not verify the complete WCS orientation,
per-source 4-half-light-radius masking hypothesis,
or the contamination outside each mask.

**Decision:** P18 was executed and its catalog-only
morphological mask is a geometrically viable diagnostic
alternative. It is NOT proven to outperform P17
3-FWHM on neighbor-cleanliness. Neither strategy
resolves the dominant sky/outer-profile degeneracy,
and the extra faint 32 stars have not supplied
independent-of-star-identity but same-coadd spatial/sky calibration.
Maintain \`SCIENCE_GATE_HOLD\` and all PSF wings,
noise covariance, neighbor-model and spatial
representativity science gates as unvalidated.
Do not tune the footprint multiplier/flux tiers
using this same Leo P outcome and then claim
confirmation from the eight same-coadd stars.

**Brick-only closure:** no further field selection is allowed.
Validate all reachable native-coadd PSF evidence and report
unknown wide-angle wings and remaining sky/contamination
uncertainty instead of constructing an ungrounded model.
Any subsequent source-injection or covariance work must
be strictly within this original native brick.


## Final single-brick empirical PSF closure — scoped decision, 2026-11-09

**Scope frozen:** \`dr11/south/1554p180/g\` only, using the original
locally verified \`image\`, \`invvar\`, \`nexp\`, \`maskbits\` FITS and
same native-brick Tractor/Gaia catalog. No additional brick, synthetic
sky from an adjacent brick, fabricated independent control field,
changes in source identities or refits of the original kernel.
The P16–P18 reports remain an archived exhaustive feasibility history;
their earlier suggestions to select another field are superseded.

### What the brick actually constrains

1. **Empirical coadd PSF core, numerical product:** 81×81 native pixels,
   frozen kernel SHA256
   \`8d976b91b4760fc8c51f607ccaf1cc59d4cac622c4c9992e5df11044bed2852a\`.
   From 10 train and 3 untouched core holdouts. Axial mean
   FWHM = 5.136902957677135 native pixels = about 1.346 arcsec
   at the original 0.262 arcsec/native-pixel scale. Full
   circular support of this finite kernel is 40 native
   pixels (about 10.48 arcsec); its central finite-kernel
   encircled fractions at 1/2/3/4 FWHM are approx
   0.80357, 0.94040, 0.97201 and 0.98319.
   These fractions **refer to a kernel normalized over its
   81×81 grid**, and may NOT be cited as fractions of
   all physical stellar flux. The native coadd core is
   reconstructed, but its science photometric bias and
   spatial representativity have not passed quantitative gates.
2. **Outer original-pixel constraints, descriptive product:** 48
   catalog-only selected Gaia stars (16 previously inspected
   bright + 24 new exploratory + 8 new same-coadd heldout
   identities), with six 4–22 measured-FWHM native
   profile bins, sky comparisons at 22–26 and 26–30
   FWHM, four P17 neighbor masks and one P18
   morphology+flux ellipse mask. P18 8-sector deletion
   is a correlated sky sensitivity test, NOT a covariance
   estimate or statistical significance. The 16 original
   P16C profiles reproduce and all star identities match.
3. **Limits on what the brick can tell us:** P18 4–6 FWHM
   medians are positive for all 16 bright stars with median
   +8.81e-5 of star peak, yet a typical azimuthal
   spread exceeds this level and the profile is
   environment/sky-dependent. At 8–12 FWHM the
   48-star median signed residual is +2.12e-6;
   median absolute residual 4.45e-5, median absolute
   P18-vs-P17(3 FWHM) change 2.43e-5. Sky-B
   maximum sector-deletion sensitivity across the 48
   has median 4.49e-5: comparable to the
   absolute residual. These figures cannot be
   fitted into a unique nonnegative, total-flux-normalized
   physical PSF. In particular, 48/48 successful P18
   sky estimates mask the fact that heldout rank 20
   has **zero** valid pixels in both 4–6 and
   6–8 FWHM stellar annuli. This is explicitly
   counted separately by the terminal QA.
4. **Single brick defines effective coadd PSF, not
   raw-telescope scattering.** Official DR11 processing
   documentation states that per-CCD PSFEx functions
   vary spatially, and since DR9 a modified extended
   PSF model is used for subtraction of bright-star
   wings from DECam images. An outer residual in
   DR11 coadd is therefore not automatically the
   untouched instrumental PSF. Sources:
   https://www.legacysurvey.org/dr11/description/
   and https://www.legacysurvey.org/dr11/catalogs/.

### Machine-readable terminal QA

\`diagnostic_metrics.single_brick_empirical_psf_closure_qa\`
with schema \`desidera.psf0-single-brick-empirical-closure/1\`
is derived **only from existing matched P17/P18 reports**
and the frozen finite core. It includes six ring-by-ring
medians (signed/absolute), P17-to-P18 mask shifts,
area/octant counts, 16/24/8 brightness cohorts,
sky A/B and azimuthal sensitivity, and a catalog-only
four-quadrant (same brick) spatial summary.
A strict provenance gate requires verified original
input hashes, matching native Tractor hash, all
48 immutable ranks and Gaia identities, exact radial
bins and the preexisting science-gate result. Missing
profiles are counted rather than set to zero.
No weighting or best-mask choice is tuned from Leo P,
and no raw pixels or kernel values are changed.

**Separate completion states:**
\`engineering_stage_finished=true\` means the single
brick's empirical core and all available native
outer-profile evidence have been reconstructed,
quantified and archived. This is a legitimate stopping
point for the *empirical PSF development stage*.
It does **not** imply \`science_adopted=true\`.
The existing \`SCIENCE_GATE_HOLD\` remains because
total physical PSF wing flux, sky/wing separation,
coadd covariance and quantified LSB photometry bias
remain unknown. Publish the finite diagnostic ePSF
with explicit radius/norm limits and retain the
non-adopted wing-constraint table. Never fabricate
PSF flux outside the measured region, force positive
wing tails from background-subtracted medians,
or apply a wide-wing scattering correction in
scientific photometry.

**Deliverable:** the existing "Save diagnostic FITS + report"
action preserves the frozen 81×81 kernel and embeds
this terminal QA in its complete report; status in
the PSF0 GUI reads "single-brick characterization
CLOSED" while uncalibrated wing science remains
explicitly withheld. No other brick is required.


## P20 same-brick joint sky/wing and retrospective spatial validation

Scope remains the original native \`dr11/south/1554p180/g\` brick.
This stage is NOT a galaxy-injection/recovery experiment; galaxy
injections, synthetic-galaxy structural parameter tests and
photometric correction bias thresholds belong to the later
photometry stage. P20 cannot silently satisfy them or switch
\`science_adopted\` to \`true\`.

**Reason to add P20:** the 48-star P18 pilot demonstrated
that the median absolute residual at 8–12 measured FWHM
is comparable to sky-ring and sector-deletion differences.
A median positive profile or a good nominal mask does not
suffice to identify the absolute wide-angle PSF. P20
separates three predeclared questions: (a) what radial
profile is compatible with a shared star component plus
each star's unknown local constant sky offset, (b) do
catalog-selected blank positions in the SAME brick
produce similar spurious "wings", and (c) is the
conditional profile predictive across independent
positions/identities within that brick?

**Input and controls:**

- Verify P17 and P18 exact Tractor native catalog SHA256,
  all 48 fixed star ranks and Gaia IDs, frozen measured
  ePSF-core FWHM, exact same-coadd native original
  \`image\`/\`invvar\`/\`valid\`/\`mask_reasons\` and target
  protection. Retain any unmeasurable ring and all
  historical selection decisions. No refilling.
- Predeclared brighter 16-star cohort:
  **12 model-building identities** at ranks
  \`[1,2,3,5,6,7,9,10,11,13,14,15]\`,
  and **4 retained identity holdouts** at catalog
  ranks \`[4,8,12,16]\`. They were included in past
  Leo P diagnostic exploration: despite withholding
  them in the P20 fit, these are RETROSPECTIVE
  controls, not previously uninspected blind
  independent observations. The eight fainter new
  identity holdouts from P17 remain recorded and
  cannot be promoted into a brightness-matched
  bright-star validation set.
- Fixed regular nonoverlap grid across this same
  native brick, spacing
  \`2*ceil(30*measured FWHM)+5\` native pixels
  (implemented with a bounded grid-stamp margin).
  Candidate blanks require >=3 FWHM distance
  from every native Tractor centroid and >=40
  FWHM from each fixed original bright star;
  target mask and full 30-FWHM original cutout
  must be available. At most four catalog-eligible
  sites are selected in EACH of four brick
  pixel-coordinate quadrants; catalog-only
  and geometric decisions are made BEFORE
  measuring blank-patch image residuals.
  A blank candidate that later fails a native sky
  or mask check is recorded as failed; not replaced.
  Entire blank original patch is subject to the
  same existing P18 catalog shape/flux mask,
  two outer sky rings and six profile bins,
  except it has no central PSF source.
  Coadd-adjacent controls are NOT independent
  exposures, and the Tractor catalog can miss
  diffuse sources, artifacts or large PSF halos.
- For each blank radial bin record raw physical
  \`median(original_pixels)-sky_B\` in nMgy/native
  pixel, coverage and missingness; report
  median, median absolute and maximum absolute
  spurious ring residual across all fixed blank
  sites. No normalization to a nonexistent
  blank-star peak.

**Joint nuisance-profile model:**

- Use the unchanged P18 **shape-mask** star
  profiles and three fixed alternative sky
  treatments: outer constant sky B (26–30 FWHM),
  inner constant sky A (22–26 FWHM), and
  pooled sky plane. The model in each of six
  broad 4–22-FWHM bins is
  \[
  y_{ik}=w_k+b_i+\varepsilon_{ik},
  \]
  where \`y_ik\` is the observed native stellar
  median residual divided by the original native
  stellar peak, \`w_k\` the shared conditional
  wing candidate, and \`b_i\` a separate constant
  background-nuisance term for each bright star.
  A bin enters only if it retains at least
  100 native valid pixels and four populated
  25-pixel azimuthal sectors (both criteria are
  fixed pre-observation).
- Fit with robust \`soft_l1\` least squares
  using the previously observed star-specific
  sky A/B difference, maximum leave-one-sector-out
  sky sensitivity, and independently chosen
  same-brick blank median absolute residual
  normalized by that star's peak as *nuisance
  scaling*. The \`b_i\` sky offsets use a
  zero-centered soft penalty of the same
  *diagnostic* magnitude. **This penalty
  anchors the otherwise nonidentifiable common
  additive wing/sky zero point.** The fitted
  absolute wing offset depends on that prior.
  No diagonal inverse-variance assumption, no
  Gaussian error bar, no covariance-based
  significance and no physically normalized
  PSF halo is claimed.
- Compare the six \`w_k\` estimates across
  the THREE background hypotheses; a smaller
  residual alone is not a positive quality gate.
  A physically nonnegative wide-wing PSF is
  NOT enforced, stitched to the finite
  normalized core, saved as a science kernel,
  or used to correct galaxies.

**Spatial and holdout QA:** predict each of
four fixed bright rank holdouts using the joint
profile trained on the 12 bright training IDs.
Report signed/absolute prediction differences
and per-star sky sensitivity, without refitting
the withheld star. Additionally run four
same-brick *leave-one-pixel-quadrant-out*
fits; omit every bright star in that
quadrant from the corresponding training
and evaluate the fixed profile at those
positions. A quadrant with too few
bright training stars or missing coverage
retains a blocked outcome; there is no
substitution with another brick or cohort.

**Output:** \`diagnostic_metrics.same_brick_joint_wing_sky_spatial_qa\`,
schema \`desidera.psf0-same-brick-joint-wing-sky-spatial/1\`,
including every blank control, ring
and withheld-star prediction. The status
\`SAME_BRICK_WING_SKY_AND_SPATIAL_FEASIBILITY_ONLY\`
signifies that the analysis was executed,
not that all data support a unique wide
PSF model. \`psf_wings_validated=false\`,
\`science_adopted=false\`, and
\`wide_angle_kernel_constructed=false\` remain
immutable for this diagnostic.

**Scientific decision:** Within this brick
one may ultimately report a restricted,
radius-bounded coadd PSF model if actual
stars, blank nulls, sky-model variations
and separated spatial positions support it,
with quantitative uncertainty. P20 supplies
those controls and measures identifiability;
its result must be reviewed before a
separate explicit science-validation decision.
If the null controls are empty/inadequate
or the sky-conditional profiles disagree,
the correct result is an unresolved outer
wing—not a false calibrated detection.
No galaxy injections are required *for P20*
or performed here; they remain necessary
for scientific validation of later
galaxy photometry.
