# Photometry Quality & Validation — scientific product and adoption contract

Date: 2026-11-10. Status: **PROPOSED SCIENTIFIC SPECIFICATION / IMPLEMENTATION HOLD**.
Scope: DESI Legacy optical (native DECam/BASS/MzLS bands, DR9–DR11), with
explicit extensibility to independent CCD and external datasets.
This file defines intended products and scientific adoption criteria;
it does **not** certify existing measurements, enable automatic total
photometry or claim any real Leo P acceptance result.

Related specifications:
- [Photometry input contract](PHOTOMETRY_INPUT_CONTRACT_P05.md)
- [PSF / scene / sky dependency order](PSF_SCENE_SKY_ORDER_P13.md)
- [SC1 / SA1 / SA2 source handling](SCENE_GEOMETRY_P12.md)
- [PSF validation](../qa/PSF_CORE_WINGS_VALIDATION.md)

## 0. Central scientific question

For a low-surface-brightness galaxy, a measured finite-aperture flux,
a profile integral, a model-extrapolated flux and a "total" flux
are **different estimands**. Agreement among them does not validate any
if they inherit the same sky or source-mask error.

Desidera must record what was actually constrained by the observations
and refuse an unqualified total magnitude when the missing-light
fraction, masking correction, large-scale sky covariance, calibration
or source membership are not sufficiently constrained.

The target galaxy is not identified with a single smooth Sérsic
component. Its **observed broad-band total** may include diffuse
stellar light, genuine compact stellar associations, resolved members
and nebular line/continuum emission. A smooth-component measurement
is a different product. A foreground/background contaminant is not
defined merely by being compact or by overlapping an isophote.

## 1. Non-negotiable input/lineage contracts (P0)

Every measurement refers to a versioned **Galaxy Project** and a
specific analysis run, not an implicit global "current best" FITS.

Record the following for each band/source/fit:
- Object/project UUID and deliberately confirmed target coordinates.
- Survey, release, instrument, native filter-throughput identifier,
  field/brick/CCD, input file paths **and content hashes**.
- Native WCS, pixel solid angle or area map, photometric scale,
  units, AB/Vega convention, observational vs dereddened flux convention.
- IMAGE, INVVAR, NEXP, MASKBITS and additional instrument validity
  with explicit bit-policy version and spatial coverage; missing
  IMAGE/INVVAR or NEXP cannot be silently substituted with zeros.
- Scene/segmentation revision, per-source membership hypotheses
  and chosen policy; immutable hard invalidity distinct from
  reversible provisional contaminant masking.
- PSF product, finite normalization extent, band/spatial applicability,
  wings QA and uncertainty; whether used for forward modelling or
  PSF-matched image colours.
- Sky estimator and training/holdout footprints, residual large-scale
  statistics, background scenario IDs, and record of upstream sky
  subtraction that cannot be undone by local fitting.
- Aperture geometry in ICRS, centre, PA (north through east), q,
  radial grid, reference-image definition, and rasterization version.
- Calculation method, code/recipe/config revision, random seeds,
  dependency DAG hashes and scientific release state.

No cross-release stacking or joint-likelihood use may treat reused CCD
exposures as independent data. The presence of g/r/i/z FITS, a saved
draft, or a successful GUI test does not imply a science-grade run.

## 2. Definitions: keep the estimands separate

### 2.1 Native signed finite-aperture flux

Let A be a specified sky aperture and V_b its **valid observed**
pixel subset in band b. Define, in a consistently calibrated linear
flux system, the explicitly observed valid-area sum

    F_valid,b(A) = sum_{p in A intersect V_b} w_geom,p [I_b(p) - S_b(p)]

where w_geom,p is the fractional pixel-aperture overlap for the native
grid and S_b is the selected residual sky model in the SAME flux/pixel
units. Pixel invalidity and missing footprint are not replaced by
a fictional zero signal. Overlapping aperture fractions must be
applied consistently to flux and uncertainty.

F_valid is a *partial-area measured flux* whenever valid coverage is
less than complete. It must NEVER be called whole-aperture or total
flux. Do not automatically multiply by A_geom/A_valid: that assumes
azimuthal homogeneity and can bias irregular dwarfs.

Keep signed fluxes and meaningful upper limits. Negative or near-zero
flux is not an error and must not be clipped before modelling.
Magnitude m_AB = 22.5 - 2.5 log10(F_nMgy) is defined only for
positive F_nMgy, and its uncertainty is not generally symmetric in
low-S/N cases. All downstream inference should use linear flux and
covariance rather than silently dropping nondetections.

### 2.2 Spatial coverage and correction

Report for each band AND each annulus/aperture:
- Geometric area, exposure-covered area, instrument-valid area,
  foreground-removed area, target-protected area, and unrecoverable
  saturated/invalid area as distinct quantities.
- Coverage fraction(s) with clear denominators and WCS/pixel scale.
- Whether missing pixels are clustered azimuthally or concentrated
  in a region that could host galaxy light; sector-level diagnostics.
- A separate estimated flux for masked/missing galaxy areas
  F_reconstructed,b(A), with model and uncertainty recorded.

An optional corrected full-aperture estimator is

    F_aperture_corrected,b(A) = F_valid,b(A) + F_reconstructed,b(A).

It is **model-conditional** and cannot inherit F_valid's purely
observational interpretation. If the underlying stellar light in
missing regions is unidentifiable, report an interval or a rejection,
not an unconditional correction.

### 2.3 Radial and model-based fluxes

Maintain independent estimators:
- F_isophotal,b(R): area-integrated fitted isophotal surface-brightness
  profile out to the **declared finite** R, including uncertainty
  in centre, PA, q(R), sky and masked-area reconstruction.
- F_growth,b(R): directly accumulated, valid-area and/or explicitly
  reconstructed curve of growth, with a covariance across radii.
- F_model,b(<R): forward PSF-convolved structural model integrated to
  a stated radius, with projection and neighbour treatment recorded.
- F_model,b(infty): formally asymptotic analytic/model flux.
  This is an extrapolation **even if the mathematical integral
  converges**, not a direct observation.
- F_smooth,b and F_compact,b: optional separately modelled
  morphological components, with their overlap/cross-covariance.
  Neither may be presented as the old stellar population alone.

Never replace an observed aperture measurement with a fitted model
without changing both its estimator label and its uncertainty budget.
A "plateau" or apparent zero annular gradient is not proof that light
outside the aperture is negligible; noise and over-subtracted sky
can create such a plateau.

### 2.4 What qualifies as a total

A science total F_total,b is a **qualified estimate** of galaxy light
over a declared target-component definition, with an explicitly
bounded remainder outside the measured area and inside excluded
source/invalid regions.

The term TOTAL_VALIDATED is reserved for a run passing Section 5.
Its record must identify:
- how much flux is observed, reconstructed within the aperture,
  and extrapolated outside its boundary;
- the effective outer radial domain and its selection method;
- the adopted definition of target galaxy light and any alternative
  source-membership scenarios;
- a quantitative interval/bound for missing light, including model
  and sky uncertainty;
- statistical and systematic uncertainties and the validation domain.

If an outer profile is consistent with zero but an extended halo
remains allowed, output a bound or UNDERCONSTRAINED, not a fake total.

### 2.5 Three distinct colour products

1. MATCHED_APERTURE_COLOR: same sky aperture in both bands; use
   per-band WCS and a validated PSF-matching or joint forward model.
   Only computed as a logarithmic colour for acceptable positive
   flux distributions; otherwise carry flux likelihood / limits.
2. TOTAL_COLOR: difference of two separately **qualified** total
   magnitudes for the same galaxy-component definition. It is
   not automatically equal to a central matched-aperture colour.
3. MODEL_COLOR: explicitly model-dependent colour of one
   spatial component, with assumptions and posterior uncertainty.

Do not use a galaxy central aperture colour as an unlabelled global
SED colour while using total luminosity from another aperture.
For SSP/CSP/SED fits preserve the same spatial domain and pass
joint flux errors and model assumptions.

A common reference ellipse defines comparable aperture coordinates;
it does not prove morphology is wavelength-invariant. Compare
single-band geometry, multiband joint detection and union-of-band
detection, avoiding arbitrary truncation in the weaker bands.
Significant blue outskirts or asymmetric H II knots must be
represented in total-light alternatives, not removed to make
isophotes smoother.

## 3. Uncertainty model: what must be carried downstream

Maintain a signed multiband flux vector f and a joint covariance or
validated equivalent posterior samples. The uncertainty model must
allow the following *correlated* contributors:
- pixel/shot/read noise and coadd resampling covariance;
- residual local and large-scale sky, including sky/R_e/n/flux
  non-identifiability and upstream reduction sensitivity;
- photometric zero-points, native filter responses, WCS/pixel areas,
  Galactic extinction correction when used;
- per-band PSF core/intermediate/wing errors, neighbour model error,
  source membership and variable masks;
- aperture geometry selection, radial plateau selection and
  projection/coverage reconstruction;
- model family and regularisation priors for extrapolated flux.

Errors cannot be quadrature-summed unless the components are
actually independent. In particular, fixed total geometry, PSF
matching, common sky footprints and common source decisions can
create nonzero inter-band covariance.

For a colour c=m_g-m_r the uncertainty contains
Var(c)=Var(m_g)+Var(m_r)-2 Cov(m_g,m_r), in a regime where magnitude
linearization is justified. For low S/N use joint flux samples.

**Selection-aware propagation:** in each Monte Carlo / bootstrap
replicate rerun detection, sky, source policy, geometry, radius
selection, masks and measurement as appropriate. Recomputing
flux inside one radius selected on the original noisy CoG does
not capture the look-elsewhere/selection bias.

Return statistical intervals, major systematic scenarios, posterior
covariance (or reproducible draws), and the contribution of each
named scenario. A narrow formal posterior alone must not establish
identifiability.

## 4. Scientific quality record per product

Proposed versioned record schema:
    desidera.photometry-quality/1

Required attributes (conceptual schema, NOT yet implemented):
- input_fingerprints, code_version, recipe_hash, object_id,
  dataset_ids, band, instrument_filter, extinction_state;
- estimand, target_component_definition, radius_domain,
  sky_aperture_definition, geometry_version, psf_version,
  scene_version and source_policy_revision;
- signed_flux_nMgy, flux_uncertainty_or_samples,
  covariance_reference, calibration_system, valid_areas,
  masked_and_unobserved_areas, model_correction;
- observed_flux_fraction_or_bound and finite-aperture
  missing-light estimate with interval (if identifiable);
- statuses[] and hard_blockers[], diagnostics,
  independent_validation_runs, benchmark_population/domain;
- adopted_for_science (boolean, default FALSE), adopted_by,
  adoption_timestamp, review rationale and provenance.

The same run may produce several estimator records in different
states. Never erase or overwrite rejected estimates; preserve them
as diagnostic evidence. A scientific adoption is a separate,
explicit, auditable action from completing an algorithm.

## 5. Stage gates and adoption semantics (P0)

Use an auditable, configurable **science policy** defined and frozen
BEFORE checking the results for a target. Do not claim universal
thresholds such as "70% coverage is sufficient", "0.01 mag CoG
change establishes a total", or "all sources within the target
ellipse are members". The numerical budget depends on galaxy
morphology, scale, instrument and intended science.

### Gate Q0 — input validity
Verified source hashes, compatible native photometric scales and
filter throughput, valid WCS, explicit band/pixel mask policy, native
coverage and no double application of extinction/zero-points.
Failure: INPUT_INVALID / CALIBRATION_UNVERIFIED.

### Gate Q1 — scene and valid pixels
Original hard-invalid pixels immutable; explicit treatment of
foreground, background, uncertain, target compact, target diffuse
and instrument artefacts. Tractor, Gaia, SGA and MASKBITS are
**evidence**, not automatically physical memberships. In particular
DR11 GALAXY and RESOLVED context bits are not blanket invalidity.
Compare reasonable alternative source assignments where needed.
Failure: CONTAMINATION_UNRESOLVED / MASK_POLICY_UNVERIFIED.

### Gate Q2 — PSF and residual sky
Band-local PSF applicability and independent stellar holdout checks,
including wing leakage at the measured surface-brightness limit.
Sky trained away from diffuse target and audited on sky holdouts,
with large-scale correlated/structured residuals.
If upstream subtraction could have removed relevant outer light,
flag it; later local fitting cannot guarantee recovery.
Failure: PSF_UNVERIFIED / SKY_UNIDENTIFIABLE /
UPSTREAM_LIGHT_LOSS_UNBOUNDED.

### Gate Q3 — coverage and morphology
Band/annulus coverage reported, including sector structure and
asymmetric or off-centre light. Any extrapolation over invalid
or unobserved regions must have independently tested recovery
uncertainty. One deep reference band cannot automatically cap
the target extent in another band.
Failure: INCOMPLETE_COVERAGE / GEOMETRY_UNDERCONSTRAINED.

### Gate Q4 — flux definition and missing light
Keep direct signed flux, isophotal integral, matched-aperture and
model extrapolation distinct. Evaluate posterior/alternative
sky-and-profile bounds on masked plus exterior missing light.
Require a science-policy-defined tolerated bias/uncertainty for
the claimed estimand. A plateau alone never passes.
Failure: NO_PLATEAU / MISSING_LIGHT_UNBOUNDED /
MODEL_EXTRAPOLATION_DOMINATED.

### Gate Q5 — calibrated end-to-end validation
Test *the exact pipeline and decision policy*, not merely the final
integrator, on independent synthetic and observed controls spanning
the intended morphology, angular size, surface brightness, colours,
crowding, sky gradients, position on bricks, PSF wings, and footprint
gaps. Evaluate bias, scatter, coverage of nominal 68%/95% intervals,
false discovery/false rejection, and accepted-sample selection bias.
Require held-out objects and/or genuinely independent deeper data,
not only agreement of correlated methods on the same coadd.
Failure: RECOVERY_NOT_VALIDATED / OUT_OF_VALIDATION_DOMAIN.

### Gate Q6 — reproducibility and release review
All upstream artifacts versioned; invalidate downstream products
when any input, mask, PSF, sky or source policy changes. A human
or separately specified science-release policy must explicitly
adopt a validated run. Failures do not silently adopt "best
available" results.
Failure: STALE_DEPENDENCY / ADOPTION_NOT_APPROVED.

**States (per estimand, separate from worker success):**
- DIAGNOSTIC_ONLY: quick measurements/preview, no adoption claims;
- MEASURED_PARTIAL: signed valid-area flux is available;
- MODEL_CONDITIONAL: a model correction/extrapolation is reported;
- UNDERCONSTRAINED: uncertainty or alternative hypotheses prevent
  a unique inference;
- SCIENCE_REVIEW_READY: Q0–Q5 passed against a documented policy;
- TOTAL_VALIDATED: only an explicitly adopted qualifying TOTAL
  estimand after Q0–Q6, with its science definition and domain;
- REJECTED: failed mandatory scientific checks; keep audit record.

An estimator can have several flags at once, but TOTAL_VALIDATED and
REJECTED/UNDERCONSTRAINED cannot apply to the same adopted estimand.

For each downstream task (colours, SED fitting, M_star, SBF, UDG
classification), the eligibility engine decides whether **that
specific product and its dependencies** suffice. SBF requires
additional dedicated PSF/noise/calibration gates beyond this document.

## 6. Numerical acceptance budgets — configuration, not invented truth

A science policy must include:
- max tolerated absolute/radial-flux bias and coverage-dependent
  systematic effects for each estimand;
- max tolerated residual missing-light fraction and upper-tail
  confidence level, or an alternate justified observable-domain
  definition; exact domain of any total-light claim;
- permitted model extrapolation sensitivity and dependence on priors;
- required recovery-test ranges and finite-sample confidence intervals
  on measured bias, coverage and false acceptance rate;
- required spatial PSF/sky scale, target membership assumptions,
  completeness/purity requirements for source recovery.

Set science thresholds **before** seeing object outputs. No single
ad hoc number is hard-coded as universally valid at this stage.
For illustration ONLY, an analysis seeking better than 0.05 mag total
accuracy would need independent evidence that its combined bias,
missing-light corrections and uncertainty calibration satisfy that
budget; a small formal optimizer error is not sufficient.

If coverage or corrected flux cannot be constrained to the declared
budget, return MEASURED_PARTIAL, MODEL_CONDITIONAL or UNDERCONSTRAINED
rather than pretending a total was validated.

## 7. Verification matrix before release

Design synthetic and real-data checks with explicit expected outputs:

| ID | Experiment | Must demonstrate |
|---|---|---|
| PQ01 | Signed near-zero/negative annulus | Flux retained, magnitude omitted, likelihood usable |
| PQ02 | Annular coverage + wedge-shaped gaps | Geometry/covered/valid/masked area separated; no hidden 1/f coverage inflation |
| PQ03 | Target RESOLVED/GALAXY contextual bits | Genuine target photons not globally invalidated |
| PQ04 | Central foreground star, overlapping H II knot | Alternative membership scenarios; no automatic erasure of target light |
| PQ05 | Blue outskirts absent in z-band | Pairwise matched colours do not silently shrink entire galaxy to z detection |
| PQ06 | Different PSFs + WCS handedness/RA wrap | Same ICRS aperture, PSF-aware colour and complete provenance |
| PQ07 | Spatial sky gradient + oversubtracted coadd | Correct identifiability warning and sky-dependent missing-light limits |
| PQ08 | Compact vs extended galaxy at same per-pixel noise | Scale-dependent limiting surface-brightness/coverage behaviour |
| PQ09 | Noisy apparent CoG plateau | Selection-aware MC prevents falsely validated total |
| PQ10 | Sérsic vs broken/asymmetric outer profile | Extrapolation model uncertainty appears, no unjustified asymptotic truth |
| PQ11 | Multi-band correlated geometry/sky error | Joint covariance changes colour uncertainties appropriately |
| PQ12 | Joint multi-release with reused CCD exposure | No double-counted likelihood/inflated precision |
| PQ13 | Injection before sky subtraction vs into final coadd | Separates upstream lost-light bias from downstream recovery |
| PQ14 | Independent holdout test of nominal intervals | Honest uncertainty calibration / failure decision with confidence bounds |
| PQ15 | Modify mask/PSF/sky after saved run | Downstream result marked stale, not silently adopted |
| PQ16 | Unresolved extra-galaxy component | Total rejected/conditional; finite-aperture product still retained |

Tests must cover resolved, semi-resolved, diffuse and irregular dwarfs,
not only a single smooth Sérsic profile. They should exercise actual
background and mask operations, not just analytical flux recovery.

## 8. Leo P acceptance sequence

Leo P is a **stress-test**, not a calibrator whose apparent colour,
H II morphology or SGA ellipse can serve as independent truth.

1. Retain exact native DR11 IMAGE/INVVAR/NEXP/MASKBITS and Tractor SHA,
   project/scene revision and per-band source inventory. Review all
   protected-overlap entries; report any inventory truncation.
2. Complete SA2 independent native-band detections, association
   ambiguity audits and independent Gaia quality checks; reserve
   HST/MUSE as external morphology/population checks.
3. Validate band-local PSF core/wings on held-out stars and explicitly
   reject/retain its applicability at the galaxy position.
4. Measure background alternatives with independent sky samples,
   spatial covariance and sensitivity to an extended galaxy profile.
5. Report signed valid-area curves of growth, isophotal integrals,
   matched-aperture colours, and separately model-conditioned total
   estimates under multiple source/sky hypotheses.
6. Run injection/recovery and benchmark tests across realistic
   morphology and extended low-SB outskirts, including source
   association errors and coadd reduction loss where testable.
7. Only then consider science adoption under a predeclared policy.
   Do not promote PR #84 or a passing SA2 test to scientific
   total-photometry validation.

## 9. Implementation order without expanding the GUI prematurely

A. Specify and test a small pure-data product schema and fail-closed
   adoption/status validator; existing fixed-aperture integration remains
   DIAGNOSTIC_ONLY / MEASURED_PARTIAL until coverage and correlated
   errors are validated.
B. Extend coverage/annulus ledger and signed-flux covariance
   without changing existing FITS or previous project manifests.
C. Add alternative source/PSF/sky scenario runner with science
   fingerprints and downstream invalidation; keep the GUI as a
   review surface for real measured evidence.
D. Implement independent multiband geometry/CoG/model estimators
   and matched-aperture outputs with explicit PSF controls.
E. Run blinded injection/recovery and quantitative acceptance,
   then enable manual adoption for qualifying products.

No new project directory structure or alternate feature branch is
required for this specification. Treat this file as a scientific
requirement, not a claim that any of these gates already pass.

## References and source-status caveat

This contract synthesizes the project methodological review and
existing architecture. The following are **methodological orientation**
rather than claims that a particular threshold is externally validated:

- DESI Legacy Surveys DR11 documentation:
  https://www.legacysurvey.org/dr11/description/
  https://www.legacysurvey.org/dr11/catalogs/
  https://www.legacysurvey.org/dr11/bitmasks/
- Stone et al. (2021), AutoProf, MNRAS 508, 1870:
  https://doi.org/10.1093/mnras/stab2709
- Moustakas et al. (2023), Siena Galaxy Atlas 2020,
  ApJS 269, 3: https://doi.org/10.3847/1538-4365/acfaa2
- PHOTOMETRY_INPUT_CONTRACT_P05.md and
  PSF_SCENE_SKY_ORDER_P13.md contain project-specific engineering
  conventions and may be stricter than generic literature methods.

Where documentation of SGA2025 products remains preliminary, it is
comparative context, not the scientific truth catalogue. The current
code base has no externally validated TOTAL_VALIDATED implementation.
