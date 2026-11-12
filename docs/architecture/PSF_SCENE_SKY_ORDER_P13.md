# P1.3 — PSF-first scene, sky and photometry architecture

Date: 2026-11-03. Status: **approved architectural requirements expressed as
a design specification, not implementation or scientific validation**.
This supersedes the *logical stage ordering* in P1.1 and P1.2a, but
does not change the numerical claims, saved science statuses or existing
source/FITS handling. Base: P1.1 -> P1.2a -> P1.2c stacked draft PRs.

This document incorporates two detailed researcher-supplied PSF/methodology
notes (NGC2146sat empirical-PSF approach and the corrected PSF/scene/sky
sequence). The source notes are design input. The official documentation
links at the end independently support the specific survey/library behavior.

## 0. Non-negotiable scientific premise

The PSF, source contamination and residual galaxy sky are **coupled**
inverse problems. An apparently successful optimizer or a visually smooth
image cannot establish that a faint outer galaxy profile is preserved.

- A final contaminant mask needs information from an appropriately tested PSF.
- A reliable PSF needs selected stars, valid original pixels and a local sky
  estimate on their cutouts, not a scientifically final galaxy sky model.
- A final galaxy sky needs a PSF-aware contaminant model and protected
  faint galactic light. During structural photometry galaxy and residual
  sky parameters can still require **joint fitting**.
- A preliminary isophotal or interactive user mask remains useful before
  the PSF has been built, but must be identified as provisional.
- Avoid a literal recursive execution graph. Scientific refinement
  creates versioned artifacts and selectively invalidates descendants.

## 1. Adopted order: eight logical stages, not eight separate GUI screens

| # | Stage | Artifacts and meaning | Prohibited automatic claim |
|---|---|---|---|
| 01 | Data Preparation | Native FITS quartet/CCD, quality masks, WCS, band, instrument, source fingerprints | File presence means accepted science |
| 02 | Scene Initialization | Scene_v0, tentative SourceCatalog_v0, rough target protection, preliminary segmentation, candidate PSF stars | Definitive target/source classification |
| 03 | Preliminary Sky | Initial detection noise/sky; distinct PSF-star cutout backgrounds after star selection | Final faint-galaxy background |
| 04 | PSF Reconstruction & Validation | Per-band empirical PSF_v0; optional CCD PSFEx/independent coadd reconstruction; coverage/held-out QA | Locally valid PSF from global stars without test |
| 05 | Scene & Smart Masks | Scene_v1, PSF-aware contaminants, manual shared ICRS anchors with band-local masks, neighbor models | “Best mask” based only on source threshold |
| 06 | Sky Reconstruction | Residual Sky_v1 and sky-model sensitivity/holdout evidence | Upstream removed galaxy halo recovered |
| 07 | Refinement Loop | Conditional assessment and revisions PSF_vN/Scene_vN/Sky_vN; selective reexecution | Optimization converged implies scientific precision |
| 08 | Galaxy Photometry & QA | Fixed/free isophotes, parametric forward models, colors, uncertainties and adoption gate | Final result with unchecked PSF/sky/masks |

The UI maintains a single Photometry Workspace and supports an
interactive preliminary/manual mask editor from stage 02. PSF-informed
mask proposals are enabled only when the relevant per-band PSF evidence
is available at stage 05. Stage 04 belongs **before** stage 05.
SBF-grade PSF remains a downstream specialized science-gated extension.

## 2. Define boundaries precisely

### 2.1. Scene Initialization (Scene_v0)

Inputs: immutable observational source references, WCS, policy-specific
native validity masks and user/project ICRS target location. Local robust
noise can be iterated internally with tentative source detection.

Outputs:
- source candidates with explicit **uncertain classification**
  (star/galaxy/unknown) and reasons/evidence;
- rough protected target footprint, preferably larger than the immediate
  visible bright component so genuine diffuse outskirts are not fit away;
- basic contaminant/bright-source and stellar-training masks by purpose;
- candidate positions for PSF construction, with cuts not predicated on
  having already removed or masked *all* useful stars.

Do not try to complete a PSF-aware scene here. A preliminary user seed
in Smart Mask can mark a suspect source, but a connected threshold
from the coadd does not certify its extended wings.

### 2.2. Preliminary Sky (Sky_detection and Sky_star_cutouts)

Keep these as **different data products** with different footprints
and dependency identifiers. Local star-cutout background is usually
constant or low-order over a small uncontaminated field, but must exclude
neighbours and must not absorb significant stellar PSF wings.

A prelim detection sky can inform star selection. The actual
`EPSFBuilder` input cutouts require background subtraction *before*
PSF reconstruction; photutils does not estimate this inside the builder.

The coadd has already undergone survey processing. In DR9+ the
Community Pipeline and Legacypipe perform background treatments,
so Desidera analyzes **residual** backgrounds and potential oversubtraction.
A locally fitted sky cannot restore unresolved signal removed upstream.

### 2.3. PSF Reconstruction & Validation (PSF_v0)

Use the empirical star-stack/ePSF approach previously employed for
NGC2146sat in DR9 North as the baseline coadd method. It must be
possible to process g, r, z separately (i only where actually available).

**StarSelector** should use Gaia (optional evidence, not a sole gate)
and local image morphology, S/N and saturation, neighbor isolation,
sky quality and spatial distribution. Selection criteria are scale-
dependent: a star suitable for a compact core cannot be assumed to
measure PSF wings out to tens of arcsec.

**EmpiricalPSFBuilder**:
- use individual background-subtracted stellar cutouts;
- refine subpixel centers and flux normalization iteratively;
- use robust stack/fit with documented outlier reasons;
- choose oversampling adaptively, not always 4x; retain native pixel
  integration/normalization conventions;
- preserve absolute/encircled-flux normalization and report unmeasured
  power beyond the empirically constrained PSF extent.

**Spatial representativity** must be measured, not assumed:
- Global empirical: stars across an observational brick/coadd, a sound
  baseline if holdout-star residuals show sufficient spatial stability;
- Local empirical: a model near the target only if enough clean local
  stars support lower expected predictive error;
- Spatial model: e.g. smooth basis/low-dimensional coefficients varying
  with detector/sky position, complexity limited by quality and
  geographical spread of training stars;
- Comparison on independent holdout stars/leave-one-out, including
  ellipticity, residual structure, radial profile, wings, systematic
  flux errors and evidence for position dependence.
- If spatial interpolation has weak support, do not silently force it.

Three meaningful PSF scales must be represented, with actual
boundaries established empirically for instrument/seeing:
1. Core: FWHM, ellipticity, asymmetry, accurate subpixel response.
2. Intermediate: a few-to-tens of arcsec, star contamination near
   the galaxy, multi-band matching, compact structures.
3. Extended/wide-angle: far wings, scattering and star halos, using
   unsaturated outer pixels and independent normalization constraints.

A complete model `P(x,y; position, band, observation)` should maintain
flux normalization on a clearly specified finite region and flag
unknown flux beyond it. The extended PSF can retain asymmetric
diffraction/reflection artifacts; a radial average alone is insufficient
when those residuals affect the target.

#### Survey-specific PSF evidence

- `SurveyPSFAdapter` for separate CCD PSFEx (or compatible source model);
  the actual per-exposure/per-position model should be preferred for
  CCD-forward modelling in Hybrid Analysis, **with an independent
  stellar validation**.
- Coadd PSF is affected by underlying CCD PSFs, weights, resampling,
  coverage and source masks. An effective coadd PSF reconstructed from
  CCDs is an *approximation* until checked on coadd stars.
- Empirical coadd PSF and PSFEx-derived coadd PSF are methodologically
  distinct but **not statistically independent observations**; they
  derive from related imaging. Do not combine them as independent likelihoods.
- DR9 North g/r observations are BASS/90Prime and z MzLS/Mosaic-3;
  DECam-specific bright-star halo models must **never** be copied
  directly to North data.
- Starting DR9, Legacy Survey applies an extended-wing bright-star
  subtraction model for DECam. Desidera must model the **remaining
  contamination in the delivered processed image**, not subtract the
  *already corrected* total star halo again.

Products: `StarSelection`, `StarCutoutSky`, `PSFCore`,
`PSFIntermediate`, optional `PSFWings`, `PSFQA`, and an immutable
`PSFModelReceipt`. These are **target interfaces only**. Do not
claim the Python classes or algorithms exist in P1.3.

PSF retrieval should be associated with the **observation dataset**,
instrument, release, band and source hashes, not merely a galaxy
project. Global PSF for the exact same coadd may be reused for several
galaxies only when spatial applicability has been checked at each
target. Local evaluation also depends on position.

### 2.4. Scene & Smart Masks (Scene_v1)

This is where the first PSF-informed mask/refinement engine runs.
- The researcher points at a source once in ICRS; preserve the shared
  source identity across g/r/z while creating **independent** local
  pixels/contours and uncertainty in each band.
- For ordinary isolated stars, an anchored connected-isophote mask is
  a quick preview. It must not masquerade as a validated
  full-wing mask or post-PSF `optimal` solution.
- For bright stars, use fitted star amplitude times evaluated PSF and
  choose a masking boundary from predicted impact on recovered galaxy
  sky/flux/profile, with a wing/halo model uncertainty envelope.
- If a bright/saturated star is ambiguous or PSF wings are unresolved,
  use a conservative mask and explicit `UNCONSTRAINED_WINGS` flag.
- For blends overlapping legitimate galaxy light, compare masking,
  joint star+galaxy+sky modelling and explicitly validated subtraction
  alternatives. A smooth residual is not a guarantee of unbiased galaxy
  morphology.
- Make each proposal visible, reversible and reviewable; keep
  native hard masks immutable, retain manual edits and per-band
  exceptions, and invalidate downstream products after accepted edits.

### 2.5. Sky Reconstruction (Sky_v1)

Fit the remaining sky around the protected, PSF-aware source model.
Compare plausible flat/gradient/smooth templates using disjoint
training and holdout sky regions and assess degeneracy with large
effective radius/extrapolated galaxy flux. Do not make low-order
spline flexibility a substitute for identifiability. Output includes
per-band residual maps, background-parameter uncertainties, correlated
noise and systematics due to star wings and survey oversubtraction.

### 2.6. Refinement and Photometry

Run additional loops **only if evidence demands it**:
- previously used PSF star found blended by Scene_v1;
- significant control-star PSF residuals/spatial mismatch;
- alternative PSF wings change the adopted star masks or galaxy halo;
- sky-model uncertainty changes target structural inference.

A newly estimated PSF_v1 can trigger new Scene_v2 and Sky_v2, while
keeping PSF_v0/Sky_v1 available as historical artifacts. Stopping
criteria include holdout residuals, encircled-energy/radial profiles,
mask stability, sky quality and downstream galaxy parameters—not FWHM
alone. Non-identifiability across plausible solutions is an explicit
systematic uncertainty, not silent successful convergence.

The final galaxy fit may jointly refine nuisance sky with target
Sérsic/other models. Coadd inverse variance is not necessarily a
complete independent-pixel likelihood due to resampling correlation.

## 3. PSF product reuse, ownership and reproducibility

**Dataset Registry/Observation Registry owns PSF identity**; each galaxy
project holds a reference to a reviewed dataset-bound PSF realization.

The initial `PSFLineageKey` contract is deliberately a **cache-key
design primitive**, not a validated PSF product, a new parallel
registry, or a proof that a hash was checked against downloaded bytes.
Store:
- exact source product (coadd vs CCD), release, region/brick/exposure,
  instrument/camera and band;
- the original observational source checksum and source evidence;
- PSF star-selection version, cutout-background version, recipe digest,
  applicable spatial scope and position if local;
- stellar candidates used/held out, rejection reasons and residuals;
- center convention, pixel response, normalization extent, flux outside
  region and per-scale model uncertainty;
- independent QA and full kernel/model as a derivative FITS artifact;
- links to original unmodified files and the native observation lineage.

Do not reuse a DECam model on BASS/MzLS, a coadd PSF for a CCD
likelihood, or a global PSF as a local PSF without spatial QA.
PSFEx availability must be verified in each particular release and
exposure; catalogue mentions alone do not prove data were downloaded.

## 4. Versioned dependency graph (minimal reprocessing)

New pure-Python contract `workflow_dependencies.py` enumerates
artifact-stage dependencies for planning, with no jobs executed.

| Changed artifact | Invalidate downstream | Do **not** invalidate solely for this change |
|---|---|---|
| Original FITS/WCS/source identity | Scene_v0, preliminary sky, PSF input/fit, final scene, galaxy sky, photometry | Stored/manual sky-vector shape itself may remain, but must be reprojected/checked |
| Star selection or PSF cutout background | PSF, final scene, galaxy sky, refinement, photometry | Unrelated external catalogue snapshots |
| PSF model or wing QA | Final scene, galaxy sky, refinement, photometry | Original FITS |
| Manual contamination-mask edit | Final scene, galaxy sky, refinement, photometry | PSF star sample/cutout sky if unchanged |
| Final galaxy residual sky | Refinement/photometry | Empirical PSF cutout sky, star sample or PSF by default |

The graph is intentionally conservative; a later real scheduler should
use **exact hashes of the consumed artifacts** so a version change
that does not alter underlying star selection/cutout background need not
rebuild PSF. The source identity graph must stay consistent with the
existing P0.3 run infrastructure, not invent a second results register.

## 5. Explicit acceptance and PSF QA gates

The GUI must expose PSF status per band and task. **No fake completion**.
Critical gates before publishable results:

- sufficiently constrained core and intermediate response, with holdout
  stellar residuals and spatial representativity at the galaxy position;
- if faint outskirts depend on star halo leakage, wing model measured or
  a conservative, propagated uncertainty/masking limitation reported;
- full distinction among unsaturated core, outer-wing fit and residual
  post-survey-halo-subtraction model;
- bootstrap/star-sample sensitivity, optional alternative PSFEx model,
  model-vs-observation residual and encircled energy;
- effect of PSF alternatives on `m_tot`, `R_e`, `n`, `mu_0` and
  outer surface brightness/colour gradients, not just FWHM agreement;
- if no suitable stars, mark PSF unavailable/underconstrained and
  provide an explicit user option to evaluate a documented survey PSF,
  rather than inventing an empirical model.

Per-use policy:
- isophotal photometry: remain on original pixels with documented
  PSF contamination and seeing/ellipticity effects;
- forward Sérsic and joint CCD likelihood: convolve each observation
  with its own appropriately localized PSF;
- matched-aperture colors: PSF matching/homogenization or joint
  multi-band forward models; noise covariance changes after matching;
- faint diffuse features: wide-angle PSF and scattered-light QA;
- SBF: area-matched PSF power-spectrum template combined with masks,
  non-white noise, alternative PSFs and injection/recovery;
- WISE: **different instrument/passband PSFs**, never reuse Legacy PSF.

## 6. GUI contract and acceptance tests

In `Photometry Workspace`, display the ordered high-level stages
without making a placeholder look like a runnable validated algorithm.
At PSF stage show source type, band, model family and completeness.
Future PSF inspector: candidate-star map, train/holdout classes, 2D
kernel, radial log profile, core/intermediate/wings, FWHM/ellipticity,
residuals, spatial stability, dependence on selected stars, local
representativity and export/QA provenance.

Manual mask editor is available for Scene_v0, but a default `Smart
Mask` sourced only from thresholded coadds must clearly say
`PROVISIONAL_PRE_PSF`. After PSF construction a **different**
`PSF_AWARE_REVIEW_REQUIRED` candidate can be offered; both require
explicit researcher acceptance. It is incorrect to hard-code a
fixed angular radius as “optimal” for all bands.

Validate both scientific and UX behavior:
1. Parent offline source/GUI regression suite unchanged;
2. dependency tests show manual-mask/galaxy-sky edits do not
   automatically rebuild a stellar PSF;
3. PSF lineage hash changes for instrument, source, band, star sample,
   cutout sky, recipe, coadd/CCD and local position;
4. no new automatic downloads, source-FITS writes or claimed science
   products from opening the Photometry Workspace;
5. first run the Leo P control on its validated local DR11 South source;
   inspect held-out stars and core/wing QA, then compare methodologically
   with the older NGC2146sat empirical PSFs as a secondary control;
   subsequently UGC 2162 and crowded/bright-star/PSF-sparse fields;
6. use injection-recovery for recovered galaxy magnitude, R_e, n,
   mu_0, color gradients and SBF when appropriate.

## 7. Recommended implementation sequence, revised

- P1.3a: this architecture, conditional dependency contract and GUI
  stage order; no fitted sky/PSF or photometry claimed.
- PSF-0: Leo P is the **first** real control; construct an independent
  empirical coadd PSF per available band, review star selection,
  cutout-sky fit, provenance/QA and two-dimensional diagnostic FITS.
  NGC2146sat remains a **second** comparison baseline.
- P1.2b: research-quality FITS mask editor with provisional masks,
  sky-coordinate linked edits, undo/redo and review.
- PSF-1: local/global/spatial evaluation + CCD PSFEx adapter with
  independent holdouts and camera-specific policy.
- P1.3b: final scene and PSF-aware Smart Mask; weak wing handling,
  model-versus-mask decision diagnostics.
- P1.3c: galaxy residual sky, heldout tests, conditional iterations
  and joint sky–galaxy fits.
- Later P1.4+: isophotal/structural photometry, band matching,
  scientific QA, then specialized SBF implementation.

Ordering in this plan distinguishes a **development sequence** from
the logical **execution order**. In particular, GUI drawing tools
can be developed while PSF-0 is implemented, but they cannot falsely
assert PSF-informed mask optimality before its validation.

## 8. Official methodological references

- Photutils EPSFBuilder requires background-subtracted stellar cutouts:
  https://photutils.readthedocs.io/en/3.1.0/api/photutils.psf.EPSFBuilder.html
- Legacy DR11 PSFEx description (per-CCD spatially varying PSF)
  and previous sky processing:
  https://www.legacysurvey.org/dr11/description/
- Legacy DR9 BASS 90Prime g/r vs MzLS Mosaic-3 z:
  https://www.legacysurvey.org/dr9/description/
- Legacy DR9 DECam bright-star extended-halo subtraction:
  https://www.legacysurvey.org/dr9/psf/

Future review should also examine Anderson & King (2000) ePSF,
AutoProf, elderflower wide-angle PSFs and injection-recovery work
before selecting exact solvers. Bibliographic discussion does **not**
by itself validate a PSF implementation.


## 2026-11-06 amendment: mandatory empirical PSF in every optical science band

**New non-negotiable researcher requirement:** Desidera must obtain and
scientifically validate an **empirical, observed-star-derived effective PSF
for every band being used in subsequent galaxy photometry**, with no
filter-specific exceptions. In the Leo P DR11 south native coadd analysis
this means **g, r, i and z**; a scientifically incomplete `g`
means **all-band galaxy photometry cannot start**. The current r/i/z
PSF-0 reconstructions remain `EXPERIMENTAL_PSF_NOT_ADOPTED` and
also **cannot** satisfy the gate. A diagnostic coverage report, a
numerically converged PSF, or a PSFSIZE seeing estimate do not satisfy
the scientific requirement. The prior pilot's honest `g` refusal
closes **only a diagnostic experiment matrix**, not the empirical PSF
development objective.

Define the required input to a future numerical photometry executor
as one *trusted internal science receipt per exact band/observational
dataset*, recording:
- empirical origin from **observed stellar pixels in that band**;
- relevant coadd/CCD observation lineage and verified original hashes;
- properly normalized final model with explicit immutable kernel hash;
- independent stellar evidence (at least 3 Gaia-identified holdouts in
  the current primary policy), extended-wing validation, spatial
  representativity, stellar/sky separability, correlated coadd noise QA;
- explicit **scientific adoption** distinct from experimental preview.

`inspect_empirical_psf_coverage` and
`require_empirical_psf_coverage` in
`src/lsbpipe/photometry/workflow_dependencies.py` express the
**fail-closed planning/entry contract** for later photometric execution.
They accept only future *trusted* internally verified science receipts;
the contract cannot authenticate arbitrary user-entered dictionaries
or replace the actual measurement QA. There is **no numerical
photometry executor yet** to invoke this guard, so these additions
must not be misrepresented as scientific validation or as a
retroactive upgrade of PSF-0. Subsequent executor implementation
**must invoke the guard before any scientific measurement is written**.

An external single-CCD PSFEx kernel, a Gaussian from PSFSIZE, and
an uncalibrated analytically imposed Moffat profile are **not
acceptable substitutes** for the mandated empirical coadd PSF.
PSFEx and analytical shapes remain useful for benchmark comparisons,
initialization and contamination models; the *final* adopted PSF
must be constrained and independently verified by real observed
stars in the required filter and observational product.

### How to recover the missing Leo P g-band stellar PSF without fake passes

**Priority 1 — improve the sample selection, not the acceptance cut.**
Read the complete `g` Gaia/Tractor first-failure inventory (972
screened). Specifically inspect *Gaia-confirmed* sources rejected by
the 85-pixel isolation condition or source-neighbour tests, rather than
reclassifying rejected Tractor galaxy models as stars. For each
candidate inspect native cutouts, weighted segmentation and companion
flux; assess whether *PSF-model-based neighbour subtraction or
neighbour-only masking* leaves a clean, symmetric stellar core and
sufficient azimuthal coverage for the 30–38-pixel local-sky
annulus. Any new neighbour-aware regime must be independently
versioned, predeclared, and verified using **g-only stellar pixels**
and fixed heldouts; modified acceptance criteria are never
retrospectively compared as if the same recipe.

**Priority 2 — leverage more observed g-band star pixels while
preserving coadd compatibility.** Explore supplementary DR11 nearby
brick data **only** when the contributing exposure mix, seeing,
position and coadd resampling can be shown to support transferring
or spatially interpolating the empirical PSF to Leo P. Do not claim
adjacent-brick PSFs are automatically interchangeable. For severe
star scarcity consider observing-frame/CCD-based **empirical**
per-exposure ePSFs, transformed and combined using actual
coadd weights and astrometry to obtain the effective g-band
PSF at the target. Validate that synthesized coadd-equivalent
kernel with stars from the **actual target coadd**, or report
unvalidated transfer and continue blocking photometry.

**Priority 3 — transparent nonstandard validation if samples remain
sparse.** Leave-one-star-out or K-fold refitting can provide
cross-fitted diagnostics of star-derived g PSF stability, but these
folds **do not magically supply three independent Gaia holdouts
of the final all-stars model** and cannot be silently counted as
the original 10-train+3-heldout gate. Any alternative policy must
define a separate statistically justified validation contract,
include injection/sky/covariance tests and be explicitly reviewed
before becoming eligible for scientific adoption.

If real g observations cannot support a scientifically validated
empirical coadd-effective PSF by any supported path, the pipeline
must continue to **block scientific galaxy photometry** and clearly
record why. The requirement does not imply that a trustworthy
PSF can be mathematically guaranteed in arbitrary poor-quality
data; it implies that missing PSF can never be silently ignored.

Primary sources: Legacy DR11 PSF and CCD-model description at
https://www.legacysurvey.org/dr11/description/ and Photutils
`EPSFBuilder` at
https://docs.astropy.org/projects/photutils/en/stable/user_guide/epsf_building.html.


---

## 2026-11-11 implementation: 03 Background built automatically from 02 Scene

**Corrected operational contract** (supersedes the earlier same-date
manual-sky-region requirement). Users must **not** draw `sky_train` or
`sky_holdout` regions to run background modelling. That was an error in
translating an optional low-level `SceneState` API constraint into the GUI.
The third UI tab is a *post-Scene residual-sky diagnostic*, not a final
scientific sky acceptance gate.

Photometry Workspace: `1 · PSF-0`, `2 · Scene & Masks`,
`3 · Background`. The third tab shares the **currently loaded** native
per-band Scene image/INVVAR/native validity and the prepared project-local
catalogue source masks, WCS, draft and target protection. It does not
fetch another brick, download catalogues, open other FITS, alter masks or
create a new directory. A project-scoped resume pointer may select
`stage="sky"`, but only an explicit `Build automatic sky candidate`
action runs the numerical worker. Changes to project/band/draft/source
invalidate previous in-memory evidence.

**Automatic selection, current implementation `SCENE_MASKS_FIXED_SPATIAL_TILES_V1`:**

1. Begin with original pixels with positive finite inverse variance and
   finite image values, excluding the immutable native-hard mask and
   prepared Stage 02 Tractor/manual source exclusions.
2. Exclude an expanded target-protection region; also use saved
   `isophote_candidate` annotations conservatively **as sky vetoes only**,
   not as adopted science isophotes. The target expansion is the larger
   of 15 arcsec and 30% of the area-equivalent target radius.
3. Expand catalogue/manual exclusions by 6 arcsec and, where the
   *already prepared* Source Registry provides Gaia G photometry,
   exclude bright-star neighbourhoods using provisional fixed radii
   from 18 to 60 arcsec for `G ≤ 18`. These heuristics do **not**
   establish actual PSF halo limits. Without the current source registry,
   fail closed instead of claiming bright-star cleaning.
4. Partition the current native Scene crop into a fixed rectangular
   grid with tile edges from 16 to 96 native pixels. Tiles with
   insufficient uncontaminated interior (<65% or <40 valid pixels)
   are rejected *using mask geometry alone*. Choose disjoint spatial
   training and holdout tiles via the fixed, reproducible rule
   `(tile_x + 2 * tile_y) % 5 == 0` for holdout; all other
   eligible tiles train. Use a 2-pixel interior border, and at most
   384 deterministically spread native pixels per tile to prevent
   dense coadd sampling from dominating the fit. No selection based
   on the eventual residual brightness.
5. Require >=8 training tiles, >=4 holdout tiles, >=100 usable
   native samples on each side, and at least three full-image
   quadrants represented on each side. If these are unavailable,
   **no sky model is fitted**; report
   `INSUFFICIENT_AUTOMATIC_SCENE_SKY_GEOMETRY`.
6. Fit both the existing robust constant and tilted-plane
   residual sky estimators to the same training samples, then
   evaluate both on the same heldout spatial tiles. A simple
   predeclared descriptive preference selects the tilted plane
   only when its heldout native-pixel RMS is below 75% of the
   constant model's; otherwise the constant level is preferred.
   This is **numerical provisional model selection**, not
   `SCIENCE PASS` or calibrated evidence of preserved galaxy outskirts.

Researcher-drawn `sky_train` and `sky_holdout` draft annotations
are ignored by this automatic selector without editing or removing
them; the core `fit_residual_background` API still supports explicit
training/holdout layers for isolated algorithm tests. Source protection
and manual exclusions from 02 remain binding. The numerical output is
`desidera.sky.project-native-residual-diagnostic/2`, containing
the tile identities, eligibility fractions, chosen pixel counts,
spatial quadrant coverage, exclusion guard radii, Gaia veto counts,
both sky models and holdout QA, their disagreement at the galaxy
anchor, and the provisional automated choice when supported.
`Copy background QA` exports this complete provenance.
No model-subtracted FITS or derived image is written.

**Science limitations:** these are patches of **the same DESI
coadd**, not independent exposures or noise-only calibration.
The detection of faint diffuse galaxy light outside the
provisional exclusion is unproven; the wide PSF wings and
mask completeness remain uncalibrated. A bright-star veto
radius is an engineering precaution, not a physical wing
model. No recovered galaxy flux, background systematic
error or final sky adoption can be inferred. Coadd spatial
covariance and upstream sky oversubtraction remain unresolved.
Scientific status remains `HOLD` until injection/recovery,
spatial and mask/sky-sensitivity QA justify adoption.

Automated selection must be tested on Leo P's **already downloaded
native brick** before reporting observed values. Synthetic unit/Qt
tests separately exercise no-manual-region success, protected
galaxy exclusion, withheld tile selection independent of pixel
brightness, fixed provenance, no usable tiles, absence of source
registry, GUI invalidation, and unchanged original arrays.


---

## Actual Leo P DR11 r Stage 03 first automatic QA — 2026-11-11

**Observed source:** user-exported complete
`desidera.sky.project-native-residual-diagnostic/2` at
`2026-11-11T14:19:25.949366+00:00`,
project `0876bbc4-7405-481b-88a8-0e3486c7be98`,
`DR11 south · 1554p180 · r`,
native image 1324×1466 pixels at 0.262 arcsec/pixel.
WCS fingerprint
`5eee1d9581324371fd042bf848b2cffbbcbef081e3edb8302262324906b000e2`
(prefixed by `wcs-sha256:`); draft
`2727f8959a9bcd15338f5a3c6abbf12eafe569e34450a14c4c59643cfd502619`;
Tractor catalogue
`8e8e3f75af868f3ef7cf4c4c9d7ca28210cbca87ef27c98351a6ab237010562c`.
No new bricks, catalogs, branches or image changes.

The automatic sky sampler returned **31 eligible tiles of 210
geometrically examined**, of which 24 train and 7 are
spatial holdouts; **9216 training native pixels** and
**2688 holdout native pixels**. Masks: 265 native-hard pixels,
393,573 preliminary auto-catalogue-mask pixels and
164,759 protected pixels. There are 540,925 geometrically
available sky pixels prior to sparse tile subsampling;
8 existing bright Gaia sources were vetoed.
The approximate fixed target guard was 69 native pixels;
catalogue exclusion extension 23 pixels. The radii
remain heuristics rather than calibrated extended PSF wings.

**Actual fitted comparisons (nMgy per native pixel):**

| Metric | Constant | Tilted plane |
|---|---:|---:|
| Residual sky at Leo P anchor | +0.00014247858 | -0.00005439155 |
| Training pixel RMS residual | 0.03963203 | 0.03962365 |
| Holdout pixel RMS residual | 0.007283561 | 0.007272366 |
| Holdout pixel median residual | +0.000005215 | +0.000175846 |
| Holdout status | MEASURED_UNCALIBRATED | MEASURED_UNCALIBRATED |

The planar model only reduces the holdout RMS by **0.154%**
relative to constant, so the predeclared 25%-improvement
gate provisionally chooses the **constant** model.
The local difference in their inferred sky at the target
is -0.00019687012 nMgy/pixel and the maximum model
disagreement at original image corners is 0.00088015471.
These differences cannot be interpreted as galaxy flux
loss or calibration uncertainty.

**Unresolved spatial quality:**
- The RMS of training pixels is **5.44×** the holdout
  RMS, compatible with substantial mixture/outlier or
  spatial-coverage differences but not diagnostic of
  their cause from summary data.
- Holdout quadrant median residuals for the constant
  model are Q0=-0.00051869, Q1=-0.00016206,
  Q2=+0.00111337, Q3=+0.00019200 nMgy/pixel.
  Their peak-to-peak span is **0.00163206**,
  approximately 11.5× the provisional fitted constant.
  The absolute fitted level is therefore much smaller
  than the spatial variation of the observed heldout
  medians. This ratio does NOT establish uncertainty
  without a calibrated sampling/noise model.
- Training spatial coverage omits quadrant Q1 while
  holdout covers all four quadrants. On the clipped
  final tile row y=13, **12/24 (50%)** of all
  training tiles and **3/7 (43%)** of holdout
  tiles are concentrated at the image edge.
  Most training information is therefore
  spatially unbalanced despite the 3-quadrant minimum.
- Only **2.20%** of the 540,925 geometry-eligible
  pixels are selected as sampled pixels; the
  dominant concern is representativity, not number
  of pixels or a conventional standard error.

**Immediate P3 follow-up (implemented, observational results pending):**
retain `SCENE_MASKS_FIXED_SPATIAL_TILES_V1` unchanged
and append read-only, method-specific **tile-level residual QA**
to the same current-band `Copy background QA`.
For each source-selected training and heldout tile,
measure the residual median and pixel RMS on exactly
the previously selected original pixels, and summarize
the unweighted median and peak-to-peak dispersion of
tile medians. Add exact row-concentration fractions,
partial-edge-tile counts, missing quadrants and
training/holdout pixel RMS ratio. Flag >40% of tiles
in one native image row or >2× training/heldout RMS,
without treating these heuristic warnings as
statistical significance tests. Never reselect or
trim tiles in response to their residuals: preserve
the predeclared holdout. A targeted future automatic
sensitivity study may compare *predeclared alternative*
tile scales/offsets in the **same project's** existing
native image, rather than searching additional
DESI bricks or asking researchers to draw sky samples.

The fitted constant is a **preliminary descriptive
candidate**, not a validated sky model. The earliest
safe next step is to run this augmented QA on Leo P
and then quantify tile-scale/sky-area sensitivity
and injection–recovery within the same native project
before allowing any sky subtraction of LSB emission.
`SCIENCE HOLD` remains binding.


---

## 2026-11-11 external literature and logic audit of Stage 03

**Basis: verified primary survey documentation and methodological research**
(accessed 2026-11-11). Reference-level provenance, with limits:

1. Legacy Surveys DR11 Data Release Description:
   https://www.legacysurvey.org/dr11/description/
   DR11 coadds already contain Community Pipeline sky removal and
   per-CCD 512-pixel-box spline sky model subtraction. DECam
   bright-star **extended-wing subtraction** is present since DR9.
   The image stacks were resampled through Lanczos-3 and are
   expressly not designed for precision work; the Tractor instead
   fits individual exposures. Do NOT frame coadd residual sky
   modelling as reconstruction of the raw night sky, recovery of
   signal already subtracted from source pixels, or a proof of
   the raw PSF's far wings.
2. Legacy Surveys DR11 File specification:
   https://www.legacysurvey.org/dr11/files/
   Coadd IMAGE is in *nanomaggies/pixel*, INVVAR is the
   corresponding *inverse variance in nanomaggies squared*.
   Pixel sampling is 0.262 arcsec/pixel. Any per-pixel surface
   brightness conversion requires division by 0.262²
   arcsec²/pixel, and signed/negative residuals cannot be
   assigned ordinary magnitudes without an explicit convention.
3. Kelvin, Hasan & Tyson (2023), MNRAS 520, 2484–2516,
   https://doi.org/10.1093/mnras/stad180 .
   Faint undetected galaxies and extended stellar/galaxy light
   leak through masks, creating positively biased sky models
   and over-subtraction of LSB flux; dilated and modelled
   source masks can help, but model mismatch is another
   source of systematic error.
4. Watkins et al. (2024), MNRAS 528, 4289–4306,
   https://doi.org/10.1093/mnras/stae236 .
   PSF-scattered light and undetected sources limit
   sky models; deeper masks and as-simple-as-possible
   parametric sky modelling are preferred. Their
   synthetic experiments do NOT validate fixed 18–60
   arcsec Gaia exclusion radii or 65% tile occupancy
   for DR11/Leo P; these are **our heuristics**.
   Errors of sky-subtracted coadds propagate and
   cannot be repaired merely by repeated adjustment
   of the already-subtracted coadd.
5. Photutils documentation, background methodology and
   Background2D:
   https://docs.astropy.org/projects/photutils/en/stable/user_guide/background.html
   https://docs.astropy.org/projects/photutils/en/stable/api/photutils.background.Background2D.html
   Even robust clipping/medians remain biased under
   undetected sources. 2D mesh size must exceed
   typical sources yet resolve *actual background*
   structure. Their `mask` for source pixels differs
   from `coverage_mask` for areas with no observation.
   Desidera's no-coverage/native-hard and source
   exclusion masks must remain distinct.
6. Gnuastro NoiseChisel manual, sky and correlated noise:
   https://www.gnu.org/software/gnuastro/manual/html_node/Sky-value.html
   https://www.gnu.org/software/gnuastro/manual/gnuastro.html
   Tile-based sky maps around diffuse extended galaxies
   can directly include galaxy flux; resampling introduces
   spatial pixel correlations. A sky image must be
   scrutinized for residual emission tracking the target.
7. Astropy CCD Photometry Guide:
   https://www.astropy.org/ccd-reduction-and-photometry-guide/v/3.0.0a1/notebooks/photometry/01.01.01-Background-removal-simulated-data.html
   The best background for *source detection* need
   not be the best background for *precision photometry*.

**Audit verdicts on our current Stage 03 implementation:**

* CORRECT with stated limits: source-bound same-project
  pixels and WCS, non-destructive masks, disjoint
  fixed spatial training/holdout tiles, candidate-only
  constant/plane models, native coadd status
  `UNADOPTED`, no external bricks/catalogues.
* P0 LOGIC DEFECT **remediated in code**: the original
  `provisional_automatic_model` ranked constant
  vs tilted plane but omitted the physical
  `zero_additional_correction_DR11_coadd`
  baseline. This confused ranking fitted residual
  candidates with a scientific decision to subtract.
  The code now evaluates **exactly the same heldout
  original native pixels without any correction**,
  publishes per-tile and full-image original RMS
  and median, and reports descriptive holdout RMS
  improvement for each fitted candidate vs zero.
  It explicitly sets
  `recommended_science_action =
  NO_ADDITIONAL_SKY_SUBTRACTION_UNTIL_VALIDATED`
  and
  `automatic_sky_background_correction_authorized =
  false`, regardless of the numerical ranking.
  The `provisional_automatic_model` field remains
  for backwards-compatible engineering provenance:
  it is **only the choice between the two nonzero
  fitted models**, and its 25%-improvement heuristic
  is **not** calibrated by the cited literature.
* P0 UNKNOWN: preexisting DR11 sky over-subtraction,
  scattered-light subtraction, missing/hidden LSB flux,
  unmapped background astrophysical emission, and
  true noise covariance are **not identifiable** from
  one already-subtracted, resampled coadd and
  its original Tractor/Gaia mask, even if a
  fitted background is spatially flat.
* P1 OPEN: native 96px sky-tile geometry, 65% source-
  free fraction, max 384 sampled pixels/tile,
  15+ arcsec target guard, fixed bright-star
  radii, and 25% model-choice rule are **engineering
  heuristics**, not validated scientific thresholds
  or literature standards for Leo P. They require
  predeclared multi-grid/mask-radius sensitivity
  experiments and injection–recovery on copies of
  the already downloaded project's same native image.
* P1 OPEN: the same-coadd 7-tile holdout is
  geographically disjoint but **not an independent
  exposure/noise realization**; adjacent tiles can
  share coadd correlations. Pixel RMS cannot estimate
  uncertainty of a global sky pedestal by
  `RMS/sqrt(number_of_pixels)`. Require
  spatial block/tile sampling or empirically calibrated
  apertures, not a false pixel-count significance.
* P1 OPEN: the sky estimate in the center of the
  target is a conditional interpolation/extrapolation
  of sparsely placed **outer** tiles. A low
  coefficient covariance/condition number does
  not certify true galaxy-free training pixels
  or absence of spatial gradients under Leo P.
  A low-confidence Scene isophote can be used as
  a conservative veto, **never** as established
  outer radius of the galaxy.
* P1 OPEN: the current pilot uses the coadd
  `IMAGE` and `INVVAR`; varying exposure depth,
  source masking, interpolated correlated noise
  and single-image residual structures can create
  different train/holdout distributions. Compare
  spatial tile medians and large aperture
  fluctuations, rather than pixel RMS alone.
* P1 OPEN: multi-filter consistency must compare
  each band using its own native image/PSF/source
  and science mask, without importing a fitted sky
  surface from r into g/i/z.

**Actual Leo P DR11 r QA reconsidered with ZERO baseline:**

The user-provided `2026-11-11T14:19:25.949366+00:00`
QA contains constant holdout RMS `0.007283560776944524`,
plane RMS `0.0072723658775994925`, constant
`b0=0.00014247857884609668` nMgy/pixel, and
its *holdout pixel-weighted mean residual*
`0.0003411811903179568` nMgy/pixel as reconstructed
from published quadrant means and counts.
Using the identity
`E[I^2]=E[(I-b0)^2]+2*b0*E[I-b0]+b0^2`,
the **exact no-additional-subtraction holdout RMS**
can be reconstructed from summary precision as
`0.007291623945262521` nMgy/pixel.
The constant reduces holdout RMS by just **0.11058%**,
plane by **0.26411%**, vs no further correction.
The difference between constant and plane is only
**0.15370%** of the constant RMS. None of these
is a p-value, a true-sky recovery, or enough
to authorize extra subtraction.

Other live QA warnings: 24 train tiles with 12/24
(50%) in the final clipped y=13 row; 7 holdout
tiles with 3/7 (42.9%) in the same partial row.
Training covers Q0/Q2/Q3 but not Q1. Pixel-level
training RMS 0.039632 vs holdout 0.007284
nMgy/pixel (~5.44x). Constant heldout quadrant
medians vary by ~0.001632 nMgy/pixel,
11.45 times larger than the fitted constant
level. This is descriptive evidence of highly
heterogeneous spatial sampling, not a
calibrated 11.45-sigma discrepancy or
a measurement of the true sky error.

**Mandatory pre-adoption gates, all internal to the
already downloaded Leo P project:**

1. Immutable versioned provenance and visible
   original mask/source/target protection.
2. No-correction baseline alongside fitted
   candidates and full per-spatial-tile residual
   maps/summary, with no post-hoc holdout replacement.
3. Predeclared fixed alternative tile grids and
   protection radii; flag unstable fit at fixed
   target position, inadequately supported/clustered
   sky samples and sensitivity to the low-confidence
   galaxy candidate mask.
4. Empirical covariance on existing image via
   multiple **large**, spatially separated
   source-vetoed apertures; treat these as
   same-coadd controls, not independent exposures.
5. Injection/recovery of artificial faint extended
   sources into **copies** of current images, with
   the exact same auto-sky selection and masking
   recipe, to estimate processing-induced flux
   loss. This measures method bias on that image,
   not the unknown historic DR11 sky-subtraction
   bias.
6. No automatic subtraction for Leo P if
   thresholds are unmet. For now **SCIENCE HOLD**,
   keep original FITS untouched. Do not interpret
   any numeric fit or diagnostic 'preferred method'
   as a final sky solution.


---

## 2026-11-11 same-brick Stage 03 grid and exponential-injection experiments

**Status: implemented as automatically executed research diagnostics; no
scientific sky correction is adopted.** This extension operates solely
on the current loaded galaxy project's native DR11 brick, WCS,
`IMAGE`, `INVVAR`, `native_hard`, provisional Tractor/Gaia
source registry, original Scene draft and target protection.
No other brick, independent star field, downloaded source data,
new directory, changed original FITS or changed adopted mask is
created or required. The worker remains background-threaded
in the existing `3 · Background` tab.

**1. Prespecified spatial-support and mask-sensitivity trials.**
The original algorithm, its 0.262-arcsec/pixel coadd and its
existing training/holdout tiles are **frozen**, then the same
automatic mask-driven sampler is rerun from scratch under
exactly five predeclared perturbations:

- `tile_72_origin`: 72-native-pixel tiles, original origin
- `tile_96_shift48`: 96px tiles shifted by (48,48) original
  native pixels (edge-partial tiles excluded at the original
  75%-coverage gate)
- `tile_128_origin`: 128px tiles, original origin
- `expanded_target_1p5`: 1.5× original target-veto dilation
- `expanded_sources_1p5`: 1.5× original source-veto dilation

All catalogue/target/bright-Gaia exclusions still apply
and remain **at least** as conservative as the original.
Each experiment constructs its own **fixed** spatially
separate training and holdout pixels from predeclared
geometry without inspecting their brightness or picking
a better-performing set after fit. If the stricter
geometry fails the >=8 training/4 holdout tile or
three-quadrant conditions, record
`SKIPPED_INSUFFICIENT_FIXED_GEOMETRY` and do **not**
replace the trial with another mask or grid.

Both constant and tilted-plane candidates are refitted
to the same current project's ORIGINAL unmodified
sky-subtracted coadd, with the zero-correction control
reported on that variant's own original heldout pixels.
Report each trial's target-centre prediction, its
difference from the original fixed-grid model,
number of accepted regions, approximate nearest
training-tile distance, and whether the target falls
within the bounding box of the training-tile centres.
The min/max/peak-to-peak sensitivity envelope
**includes the original grid**, even if all alternative
grids are unavailable. It is neither an uncertainty
interval nor significance because the trial predictions
are correlated and the model-family bias is unknown.
Holdout RMS on different variants uses DIFFERENT
spatial samples, so cross-variant RMS is not a
controlled paired score. The sensitivity experiments
cannot alone authorize model subtraction.

**2. Controlled exponential artificial-galaxy injection.**
On *temporary in-memory copies* of the current native
`IMAGE` (never written to disk), add three deterministic
two-dimensional exponential light distributions with
axis ratio 0.65, central amplitude
0.005 nMgy/native pixel and scale lengths
0.5, 1 and 2 times the area-equivalent radius of
the current *provisional* target-protected region.
The injection is centred at the same verified native
project target anchor. Each artificial profile is
**intentionally NOT a PSF-convolved physical galaxy
model**; its morphology and luminosity are test
parameters, not estimates for Leo P.

Hold the EXACT original Source/Scene masks, original
training/heldout sky pixel identities and an unmasked
aperture of radius twice the original target-area-
equivalent radius fixed. Fit the same constant and
tilted-plane estimators to each injected image. Compare
the fitted sky field to the pre-injection fitted sky at
the *same* aperture pixel coordinates and compute

\[
F_{\rm in,usable}=\sum_{p\in A_{\rm valid}}S(p),\qquad
\widehat{F}_{\rm rec,usable} =
F_{\rm in,usable}-\sum_{p\in A_{\rm valid}}
\left[B_{\rm injected}(p)-B_{\rm original}(p)\right].
\]

The ratio \(\widehat{F}_{\rm rec,usable}/F_{\rm in,usable}\)
is explicitly a **fixed-unmasked-aperture estimate of
the sky-estimator's response**, not a full image
reduction, detection, segmentation, PSF/model subtraction,
uncertainty or unbiased real-galaxy flux recovery.
Pixels excluded by original neighbour/hard masks are
reported as missing from the synthetic aperture rather
than silently assigned interpolated galaxy flux.
The **zero-additional-sky-correction** control always
recovers exactly 1 by linearity, which is a definition
of this controlled experiment, not evidence that the
already sky-subtracted DR11 input is unbiased.
No upstream lost galaxy light can be reconstructed.

**3. Scientific acceptance barrier.** New
`sky_grid_sensitivity_qa` and
`sky_injection_recovery_qa` fields are nested in
the existing current-band `Copy background QA`,
with full fixed recipes, per-trial skipped statuses,
source/mask provenance and non-adoption flags.
Neither exercise supplies a reliable **absolute sky
accuracy**, true galaxy morphology or native-coadd
noise covariance. The automatic numerical preferred
constant/plane estimator remains an **engineering
candidate only**. The science recommendation stays
`NO_ADDITIONAL_SKY_SUBTRACTION_UNTIL_VALIDATED`,
with `automatic_sky_background_correction_authorized=false`
and `SCIENCE HOLD`.

**4. The next possible methodological step** requires
a comprehensive, predeclared error budget over spatially
separated sky blocks (effective large-aperture
fluctuations), signal-to-noise and flux biases under
more realistic PSF-convolved dwarf galaxy injections,
and g/r/i/z-specific independent PSF and mask
dependencies. Do not expand project bricks to perform
these checks; where existing data cannot identify true
far wings or galaxy-free sky, continue HOLD.


---

## 2026-11-11 Leo P r observed fixed-grid and exponential recovery

**Researcher-provided GUI readout, not a new independent FITS run:**
`DR11 south · 1554p180 · r`, project
`0876bbc4-7405-481b-88a8-0e3486c7be98`;
draft fingerprint
`2727f8959a9bcd15338f5a3c6abbf12eafe569e34450a14c4c59643cfd502619`.
Source has 24 training/7 holdout tiles (9216/2688
original native pixels). No manual sky annotation was required.
The complete machine-readable `Copy background QA` (as opposed to
this GUI readout) is still needed to verify all source-mask,
trial-sample and injection provenance in a single object.

**Actual conditional sensitivity:**
`tile_72_origin`: 41/13 train/holdout tiles, both models
measured; `tile_96_shift48`: 18/4, both measured;
`tile_128_origin`: only 7/2, appropriately
`SKIPPED_INSUFFICIENT_FIXED_GEOMETRY`;
`expanded_target_1p5`: 24/7, both measured;
`expanded_sources_1p5`: 13/6, both measured.
Over the original grid plus four measurable perturbations,
constant target-sky predictions span
`+0.00002335535..+0.00041074240` nMgy/pixel,
a full range of `0.00038738704` nMgy/pixel;
tilted plane spans
`-0.00007173299..+0.00054140047`,
range `0.00061313346`. These are **not** statistically
calibrated confidence limits. The first constant estimate,
`+0.00014247858`, lies inside a perturbation range
almost 2.72 times its amplitude. The original 96px
sky training footprint has 12/24 clipped bottom-row
tiles, and the holdout has 3/7 clipped bottom-row
tiles; Q1 has no training region.
The model-free DR11 coadd remains the safer operational
default while these structural weaknesses remain.

**Actual frozen-mask artificial-source stress:**
The user reports unmasked-aperture recovered fractions
for injected deterministic exponential profiles:

| Injected scale length \(h/R_{\rm protect}\) | Constant fitted sky | Tilted fitted sky |
|---|---:|---:|
| 0.5 | 0.96615471 | 0.94854787 |
| 1.0 | 0.83932383 | 0.78059261 |
| 2.0 | 0.61203313 | 0.53651548 |

Equivalently, sky-fit-induced *conditional*
unmasked-aperture losses are 3.38/5.15%,
16.07/21.94%, and 38.80/46.35%. The source
protected-mask area in the earlier current-draft
Leo P QA was 164,759 native pixels, corresponding
to area-equivalent `R = 229.0` native pixels
(`≈60 arcsec` at 0.262 arcsec/pixel). Provided
this same protected mask is unchanged, the three
injection h are approximately **30, 60 and 120
arcsec**, and the fixed test-aperture radius is
approximately **120 arcsec**. These are **not**
Leo P structural parameters. The 120-arcsec
profile is intentionally broad enough that a
nonzero synthetic flux can extend far into
the nominal sky controls.

**Critical limits:** the experiment injects light
*after* Legacy's original sky subtraction; it tests
only additional local residual-sky-fitting bias.
The fixed scene/training masks, unconvolved
profile (q=0.65) and unmasked aperture omit
PSF effects, full source detection, galaxy
segmentation, neighbour flux filling, photometric
aperture uncertainty, native coadd covariance
and historic DR11 source over-subtraction.
The reported recovered fractions must **never**
be described as measured Leo P total-galaxy
flux loss or inverted as a photometric correction.
The 4/4 apparently measurable variants share one
original correlated coadd and have unequal sky
sample geometry. Large differences at the
target centre are a robust *sensitivity
warning*, not four independent measurements.

**Follow-up implementation after this QA:**
the frozen-mask injection now reports
`injected_signal_on_frozen_sky_samples` with
original training and holdout sample counts
and explicit injected-signal medians, maxima,
and training 90th percentile, computed from
the *same fixed native Scene masks used by the
robust fitted sky*. This directly tests the
source-to-background-sample leakage driving
the synthetic biases, rather than only the
end recovery fractions. For float32 input FITS,
the additive test image is held in float64,
not quantized back to native float32; original
image, native mask and inventory remain
unchanged. No external brick is used. These
leakage values still need a new Leo P
`Copy background QA` export; they cannot
be reconstructed from the current text-only
readout without the original native pixels.

**Scientific decision:** maintain
`recommended_science_action =
NO_ADDITIONAL_SKY_SUBTRACTION_UNTIL_VALIDATED`,
`automatic_sky_background_correction_authorized=false`,
and `SCIENCE HOLD`. Further development should
use predeclared large-aperture spatial covariance
and several PSF-aware artificial dwarf profiles
on copies of *the same already attached project*.


### PSF-0: explicit internal-holdout selection leakage and spatial coverage QA (2026-11-13)

The empirical ePSF fitting cohort and its fixed held-out
stars are disjoint, but in catalogue-vetted modes the image
core-width reference used to screen stars is learned from a
population of Gaia-confirmed candidates **before** the final
training/holdout split. Thus the holdout is independent of
the *kernel fit*, not a completely blinded population with
quality thresholds learned exclusively from training stars.
Spatially distant stars and a narrow magnitude range may also
prevent reliable transfer of one finite kernel to a galaxy.

PSF-0 now records explicit `stellar_spatial_and_brightness_coverage_qa`
with per-cohort native-quadrant counts, spatial bounding boxes,
Gaia-G magnitude ranges and unambiguous flags distinguishing
fit independence from preselection independence. The existing
star selection, ePSF kernel, heldout numerical residuals and
all non-adoption gates are unchanged. Synthetic fixture
regressions verify exact cohort/quadrant accounting and
that no fictitious external-validation or PSF-wing claim appears.
This supports designing a future genuinely independent stellar
validation cohort and does not itself cure PSF wing/sky degeneracy.
