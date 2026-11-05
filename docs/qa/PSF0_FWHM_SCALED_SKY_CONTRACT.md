# PSF-0: FWHM-scaled local-star sky geometry (experimental)

**Scope:** Leo P-first DR11 south / g,r,i,z native coadd PSF-0. The original FITS are read-only. This change is a diagnostic **recipe change**, not a statement that PSF wings/sky are solved, and not a scientifically adopted PSF. PSFSIZE values are survey-weighted seeing estimates, **not** empirical ePSF measurements.

## Diagnosed problem

For Leo P DR11 `g`, the original native coadd `psfsize-g` map gives `1.721 arcsec` at `0.2620 arcsec/pix`, a detection FWHM prior of `6.57 pix`. The former fixed Extended sky annulus `19–24 pix` samples only **2.89–3.65 × FWHM**, where stellar wings may bias the fitted local sky plane. Sigma clipping does not eliminate smooth wing light. This invalidates any automatic claim that a low stellar residual confirms correct background removal.

Official photutils guidance: https://docs.astropy.org/projects/photutils/en/stable/user_guide/epsf_building.html . DR11 seeing/FWHM: https://www.legacysurvey.org/dr11/files/ .

## Active, exactly four, deliberately non-adopted recipes

| Mode | Cutout / annulus / isolation | Interpretation |
|---|---|---|
| `baseline` | fixed 31 pix / 10–15 px / 32 px | Retained *as an explicitly wing-sensitive reference only*, NOT correct sky |
| `outer_sky` | adaptive `r_inner=ceil(4.5 FWHM)`, `r_outer=ceil(5.75 FWHM)`, `size=max(51,2*(r_outer+2)+1)` odd, `isolation=size+4` | Conservative expanded local sky sampling, **not proof of sky free from wings** |
| `outer_sky_robust` | same geometry, plus ≤2 passes of *guarded* train-only residual screening | Screen is still disabled automatically when it tracks star amplitude, and holdouts are not used to choose rejections |
| `outer_sky_convergence` | **exactly the same** adaptive sky annulus, stamp, isolation and untrimmed star sample as `outer_sky`, but `max_iterations=16` instead of 8 | Controlled Photutils convergence A/B (especially for Leo P `i`), never automatic PSF promotion |

At `FWHM=6.57 pix`, both extended modes have `81 × 81` cutouts, `30–38 pix` sky annulus, `85 px` minimum neighbour isolation; ratios **4.57–5.78 × FWHM**. At `FWHM=4.8 pix`, extended geometry is `61 × 61` and `22–28 pix`. Both examples are based on a *guessed seeing prior*, not verified physical ePSF.

If an adequate annulus would require a cutout larger than the PSF-0 safety limit of `91 × 91`, **raise an actionable error** rather than silently moving the sky closer to the star. A baseline reference remains possible but cannot be presented as a validated replacement. The enlarged cutout and isolation constraints may reject many stars; the number of accepted train/holdout sources must be evaluated on real Leo P FITS, not simulated as a pass.

## Scientific restrictions and next validation

1. Per-star 2D signed residual, core/outer profiles, radial trend of sky annulus and octant spread remain **descriptive**, not covariance-calibrated chi-square.
2. `EPSFBuilder` requires background-subtracted input. The Desidera robust planar sky fit is upstream of the ePSF, but it is not a joint star-wing + residual sky model.
3. A `>=4.5 FWHM` annulus may still contain the wings of bright sources, especially if nearby sources or wide-angle scatter are present.
4. Changing seeing or cutout size changes the detection/selection footprint. Cross-recipe comparisons require matched Gaia source IDs, consistent input hashes, and a scientifically justified common population.
5. The batch dispatcher uses distinct `PSFSIZE_{G,R,I,Z}` priors, except explicit per-band manual overrides. A manual `star-XXXX` exclusion from another FWHM is **rejected**, not silently transferred.
6. The detailed `griz` debug JSON includes exact effective recipe + SHA256 and original source hashes for rejected and completed outputs. Strict project UUID, release/brick, and optical-band checks prevent cross-band displays.
7. Future science QA: injection–recovery under a simulated smooth wing, blend, bright/saturated star and sky-gradient grid; multi-annulus background sensitivity; stratified heldout star comparisons (magnitude, colour, position), photutils subpixel-phase diagnostics, and covariance-aware goodness of fit. None are completed here.

**Scientific status:** `EXPERIMENTAL_NOT_ADOPTED` / `DIAGNOSTIC_NOT_ADOPTED`. Input source files have not been modified.


## Leo P griz evidence (2026-11-04) and two scientifically distinct catalog policies

The real DR11 south brick `1554p180` four-band QA was inspected on the local
native FITS, with independently verified Tractor source SHA-256 and per-band
PSFSIZE priors. Under the FWHM-scaled `outer_sky` recipe:

| Band | Prior FWHM (native px) | Strict Gaia/Tractor train+holdout | Signed negative mass | ePSFBuilder converged |
|---|---:|---:|---:|---|
| g | 6.57 | selection underflow: 10 accepted, minimum 13 | not fitted | not fitted |
| r | 4.88 | 16+4 | 1.235% | yes |
| i | 5.12 | 10+3 | 2.133% | **no** (8 iterations) |
| z | 3.99 | 39+10 | 1.430% | yes |

This is observational diagnostic evidence, not a calibrated model adoption.
The g run examined 969 DAO candidates, of which 692 were first rejected by
`NOT_GAIA_DR3`; this category includes **unverified sources of all types**,
not automatically usable stars. Other leading g reasons were 109 strict
Tractor neighbours and 109 missing centroid crossmatches. The i model is
not converged even though its final fitted-star centroid convergence
fraction is 100%: the Photutils build stopped at the configured maximum,
so the complete convergence status must be used. The negative-mass figures
are measured on **signed pre-clipping kernels**, not photometric flux losses.

Two noninterchangeable candidate-training policies are now provided:

- `gaia_strict` (unchanged): original native DR11 Tractor `TYPE=PSF`
  and `REF_CAT=G3`, valid Gaia ID, masks, fitted flux, isolation, shape
  and local sky. Both training and holdout stars have Gaia-confirmed IDs.
- `tractor_psf_validated` (new **EXPERIMENTAL TRAINING POLICY**): same
  original DR11 Tractor brick and the same conservative source/cutout
  gates, with the relative core-width reference measured exclusively
  from Gaia DR3-confirmed stars, but `TYPE=PSF` stars lacking Gaia IDs
  are eligible **only for training**. The expanded mode rejects the
  sample if fewer than five Gaia reference cores establish that prior. Gaia-confirmed stars are independently reserved for the
  holdout **before** building or iteratively screening the training PSF.
  It fails closed if there are fewer than the required Gaia holdouts,
  even if many non-Gaia Tractor stars are available.

Neither `TYPE=PSF` nor a fitted-flux S/N from a coadd is an independent
spectroscopic or high-resolution proof of stellar morphology. Gaia
holdouts are independent of the *ePSF fitting and train-only trimming*, but
the initial candidate preselection (e.g. masks, core moment width) still
uses their cutout pixels. Consequently, Gaia-only holdout residuals are
not an unbiased validation of the preselection classifier. Different
training policies may yield different PSF star samples. Compare
source-ID overlap and report brightness, position, and noise dependence
before attributing residual changes to any individual gate.

A per-candidate `candidate_screening_audit` in the full griz QA JSON now
records each inspected DAO detection's position, peak, nearest DAO
neighbour, catalog-match evidence (including ref-cat, optional Gaia ID,
Tractor type, fitted flux/S/N, mask/fit bits and nearest catalog neighbour),
first refusal reason, passed gates, local-sky QA and shape QA where they
were actually reached, and its pre-model train/holdout role. **Missing
metrics mean that the earlier failure prevented that test**; they are
not imputed. The audit reports `n_screened_candidates` and
`n_unexamined_candidates`, so DAO-only diagnostic caps cannot be
mistaken for complete census. No full native pixel arrays are exported.

The training-only policy is a research control, not a recommendation to
replace Gaia-strict as the default, and it cannot be science-adopted.
The native FITS/Tractor sources remain read-only. Further acceptance
requires explicit per-band injection/recovery, sky-wing and spatial
PSF validation and covariance-aware heldout residual assessment.


## Independent real-image follow-up: Tractor-vetted pilot (2026-11-04 15:09 UTC filename)

The uncommitted local researcher export
`Desidera_LeoP_PSF0_griz_debug_20261104_150909.json` uses the exact
DR11 south `1554p180` Tractor SHA-256
`8e8e3f75af868f3ef7cf4c4c9d7ca28210cbca87ef27c98351a6ab237010562c`.
These figures are **observed diagnostics only**, not clean PSF measurements:

| Band | Total DAO examined | Train (Gaia/Tractor-only) | Gaia-only holdout | Signed negative absolute mass | Solver/QA outcome |
|---|---:|---:|---:|---:|---|
| g | 972 | 11 candidates before split (10 Gaia, 1 Tractor-only) | no model | n/a | Insufficient total 11 vs 13 required |
| r | 882 | 16/10 | 4 | 1.707% | Converged, still NOT_ADOPTED |
| i | 1427 | 10/12 | 3 | 4.825% | Not converged after 8 iterations, blocked |
| z | 675 | 36/39 | 10 | 2.719% | Converged but 2.5% experimental screen failed |

In r/i/z the **training** median native positive-moment core sigmas for
Gaia and Tractor-only are respectively
2.30/3.15, 2.34/3.18, and 2.19/2.68 pixels. The corresponding
median native **peak** fluxes in nMgy/pix are
2.879/0.094, 2.654/0.106, and 2.970/0.214. Comparing unmatched
medians is **not a causal test for galaxy contamination**:
positive-clipped moments are S/N dependent. For example, in the
z 0.2–0.5 nMgy/pix peak bin, already observed median widths
were 2.46 pixels for six Gaia-associated objects and 2.58 for
19 Tractor-only objects. The other bands have poor faint-end Gaia
coverage and hence insufficient brightness-matched controls.

The negative power remains concentrated in the model's **outer**
radii: i contributes about 4.70% of total unsigned mass beyond
2 assumed FWHM and z about 2.72% beyond 2 assumed FWHM. Negative
absolute mass is an uncalibrated shape heuristic, NOT an error rate
or negative physical light. Photutils `converged` is the actual
builder-level Boolean, not just the fraction of recentered stars;
i had a final converged fraction of 1.0 but `converged=False`
because an additional center-shift criterion was unmet.

**Instrumentation repair in this PR:** preserve exact candidate and
training-cap roles in the full QA JSON (including underflow and
cap-unused rows), record median peak flux, native peak/RMS and
positive-core moment width separately for Gaia and Tractor-only,
and provide fixed, **descriptive** peak-flux strata with sample-size
gates on the interpretation. No new automated classifier, PSF
acceptance criterion or loosened selection threshold is introduced.

**Next physical validation:** match Gaia holdout `REF_ID` and source
hashes between strict and pilot runs, perform matched-brightness
train-set comparisons where actual overlap exists, test cutout
background/ring sensitivity, and then investigate PSF spatial
variation and broad wings on real data. Do not use more Tractor-only
stars as evidence of better ePSF by itself.


## Approved PSF-page visualization contract (2026-11-04)

The approved GUI prototype at the project level requires a **persistent**
per-selected-band scientific view: a native PSF model image, a radial
profile/residual diagram and a quality/status summary. These three science
surfaces must never be hidden behind the detailed-QA disclosure. The left
Galaxy Workspace navigation and independent g/r/i/z states remain intact;
per-band configuration and griz batch/export follow the science results.

Data binding is deliberately stricter than a cosmetic illustration:

- **Successful experimental fit:** the 2D model is the *actual in-memory*
  normalized, nonnegative \`PSF0Result.kernel\`. It remains labelled as
  provisional/non-adopted; its signed pre-clipping diagnostic is preserved
  separately in QA. The radial diagram is generated from the real,
  in-memory, local-sky-subtracted native cutouts and the *signed aligned*
  model from the existing star-review worker snapshot. Each training star
  uses its measured native peak as a display normalization and identical
  1-pixel annular masks in observed and model arrays. Each displayed radius
  pools valid stars using the median, for observed intensity, matched-model
  intensity and **observed−model signed residual**.
- **Signed-kernel quality rejection:** an actual computed but
  scientifically rejected signed model may be shown in diagnostic colors
  with explicit REJECTED/non-adopted status. Radial plots are permitted
  only if the worker also supplied the actual paired cutouts/model.
- **Selection underflow / no computed fit / missing inputs:** model and
  radial cards show explicit unavailable explanations. The GUI does **not**
  invent a PSF model, normalize an unavailable model, fabricate residual
  curves or reuse a prior-band model. The complete selection failure
  remains accessible in per-band full QA.
- The radial diagram is **TRAINING-IN-SAMPLE**, not an independent Gaia
  holdout, calibrated chi-square, full PSF-wing validation, or spatial
  PSF constraint. Only stars still marked \`train\` are aggregated;
  \`excluded\` and \`auto_rejected\` cutouts cannot silently contribute.
  The original independent heldout metrics and immutable source
  provenance remain accessible through the same QA surface.

Implementation reuses \`psf_star_review.py\` rather than introducing an
additional image pipeline or directories. All images are rendered in the Qt
GUI from previously computed arrays, without network access, FITS writes,
or changes to any scientific fitting/selection thresholds.

### 2026-11-04: real dark-theme screenshot readability repair

A local Leo P quality-rejected (`G`, training 75, negative-mass 3.23%)
screenshot exposed a **GUI-only** regression: the 540×315 radial chart
was placed into a narrower/shorter QLabel with native-size pixmap
semantics. Qt clipped the header, parts of the axes and residual labels.
The dedicated `_AspectFitProfileLabel` preserves the original plot as an
immutable presentation source and scales it on every resize with
`Qt.KeepAspectRatio`; the fixed-height viewport always displays the
entire scientific diagram. The chart itself now labels the signed
residual as `observed − aligned model / measured stellar peak`, reports
its actual plotting scale and the native-pixel radius ticks, and clips
curves to their respective upper/log and lower/signed plot rectangles.
A Qt regression test covers multiple widths and stale-plot invalidation.
No photometric, training/holdout or science-adoption logic changes.


## PSF-0 methodological hardening: signed holdouts, solver convergence and lineage

The controlled solver experiment \`outer_sky_convergence\` changes **only**
the Photutils EPSFBuilder iteration limit from 8 to 16. This avoids
confounding iteration budget with seeing/FWHM, candidate selection,
local-sky annulus, star masks, and optional star-level trimming. Both
models remain diagnostics, including if the 16-iteration fit converges.
Photutils documentation cautions that a non-converged ePSF can have
percent-level shape/centroid error:
https://docs.astropy.org/projects/photutils/en/stable/user_guide/epsf_building.html

**Nonconvergence is now a scientifically blocking PSF-0 outcome:** when
Photutils returns \`converged=False\`, Desidera preserves the signed kernel,
the real star-review cutouts, the full iteration history and fixed holdout
residuals in a structured \`PSF0QualityError\`. It never labels the model a
successful provisional PSF merely because the
\`final_converged_fraction\` was 100%. More iterations do **not**
guarantee convergence or improved scientific quality. With older
Photutils versions returning only a two-tuple (no convergence Boolean),
the preview remains explicitly `EPSF_CONVERGENCE_NOT_REPORTED_BY_LIBRARY`
and must not be interpreted as convergence-verified. Real Leo P runs
should use a pinned Photutils version exposing this result field.

**Holdout completeness is now a blocking integrity requirement:** all
originally reserved holdout candidates must have a valid *signed*
core-amplitude fit and valid median signed/absolute residuals in **each**
of the three radial ranges (at least four usable native pixels per zone).
Failures and missing annuli cannot be silently omitted from the
denominator. An incomplete test results in
\`INCOMPLETE_FIXED_HOLDOUT_RESIDUAL_QA\`, with the original holdout IDs
intact and the signed diagnostic model available only for review.
This integrity gate does not make the residual metrics covariance
calibrated, independently validate Gaia preselection, or prove PSF
spatial representativity.

In `gaia_strict` and `tractor_psf_validated` the fixed holdouts have
Gaia DR3 identifications (not absolute proof that their cutouts are clean).
In `dao_diagnostic`, the holdouts are **unverified DAO detections**;
their independence is limited to the ePSF training pixels. They must not
be counted as Gaia-vetted or Tractor-only sources. The QA now distinguishes
`n_train_unverified_dao`, `n_holdout_unverified_dao`, and
`holdout_identity_validation`, and uses
`DAO_UNVERIFIED_HOLDOUT` for their final selection roles.

A separate existing diagnostic of the nonnegative clipped-kernel core
now uses the actual **subpixel star centroid** as the centre of its
residual-evaluation aperture, consistent with shifting the model to
that same location. It reports the aperture center, radius and number
of native pixels. Previous code mistakenly centered that aperture on
the unshifted model grid.

**Reproducibility:** the method contract is versioned
\`desidera.psf0-coadd-epsf/2\`. The observation-bound PSF lineage key
now includes a digest of the configuration, effective Photutils solver
defaults, recentering function, relevant package versions and the new
holdout measurement policies. A version change cannot silently reuse
the same lineage key, although a fully reproducible scientific release
still requires a pinned software environment and archived original
FITS. The original \`recipe_digest\` remains separately available.
This change intentionally makes new lineage keys incompatible with
the old diagnostic keys; original FITS and prior experiment JSON remain
unchanged.

## Criteria to close the PSF-0 pilot branch, not to adopt a science PSF

1. Linux and Windows Python test matrices, Qt snapshots, holdout
   integrity and lineage regressions pass on the same HEAD.
2. Re-run real **Leo P g/r/i/z** with Gaia-strict and experimental
   Tractor-vetted selection as appropriate; separately test
   \`outer_sky\` (8 iterations) versus \`outer_sky_convergence\`
   (16 iterations) on the same native FITS, sample IDs, sky and masks.
   Do not assert a result from synthetic fixtures as real-data success.
3. Export and inspect complete original-FITS source SHA-256,
   candidate-selection QA, immutable heldout IDs, convergence history,
   signed negative mass **by radius**, holdout completeness,
   brightness-dependent star morphologies and radial residuals.
4. Do not bypass an underflowing g sample (11 vetted, minimum 13)
   through invented stars; do not accept i solely for convergence; do
   not lower the arbitrary 2.5% signed-negative-mass diagnostic
   screen to force r/z to pass.
5. Merge/close the draft only after the above reproducibility and
   review evidence is accepted. **Even a fully passing PSF-0 milestone
   is a non-adopted, global, core-only diagnostic.** The eventual
   science PSF requires spatial representativity, external PSFEx
   comparison, extended-wing/sky modelling, injection recovery and
   correlated-coadd-noise validation. Those are distinct subsequent
   PSF-science milestones, not validations supplied by PSF-0.


## Additional real Leo P evidence, 2026-11-05 16:35 UTC

Two uploaded griz JSON snapshots have matching `i`-band hashes for the
original DR11 `image`, `invvar`, `maskbits`, and `nexp` inputs.
They do **not** constitute an eight-versus-sixteen iteration experiment:

- At `16:35:00Z`, the selected `i` recipe was
  `outer_sky_robust` (8 iterations, `tractor_psf_validated`,
  `training_outlier_passes=2`), but `science_report_or_rejection`
  was null. Its QUALITY REJECTED status alone provides no usable model
  or convergence metrics.
- At `16:35:18Z`, the completed `i` run used
  `outer_sky_convergence` (16 iterations, **`dao_diagnostic`**,
  `training_outlier_passes=0`). It reported **56 train + 15
  computational holdout DAO candidates**, with **zero Gaia matches and
  no Tractor catalogue used in source screening**. Of 1427 detected
  DAO candidates, only 375 were screened (298 rejected for proximity,
  5 stamp edges, 1 masked stellar core; 71 accepted). The other
  1052 detections were not screened under the diagnostic cap.
- The 16-iteration signed kernel has negative absolute mass
  **3.7075%** (2.50% diagnostic threshold), **entirely beyond
  2 guessed FWHMs**. Builder `converged=False`; final stable-center
  fraction **14.29%** versus required 95%. The iterative fractions
  fluctuated, reaching 92.86% in iteration 10 but declining to 14.29%
  by iteration 16; the last model change was 0.87% of peak. This is
  **not** evidence that adding iterations improves convergence.
- The fixed 15/15 DAO holdout candidates had successfully reported
  residuals in all radial zones, but their stellar identities were
  not independently validated and coadd noise correlations are
  unmodeled. Median absolute signed residuals divided by measured
  peak were **0.0388 / 0.0175 / 0.0120** in core/intermediate/outer
  zones, respectively. These are not calibrated confidence intervals
  or scientific model-acceptance metrics.
- The **g, r and z bands were not computed in these two exports**:
  their original FITS paths were present, but their per-band
  source-specific scientific reports were null. The fact that the
  top-level selector is `dr11/south/1554p180/g` does not make an
  `i`-band result into a `g`-band measurement.

The source-level audit had previously labelled unverified DAO holdouts
`GAIA_HOLDOUT` and counted them as `Tractor-only`, which falsely
suggested independent catalogue-based selection. This metadata bug is
corrected here without relaxing any cuts or changing the fitted PSF.

**Outstanding closure requirement:** run `outer_sky` (8 iterations)
and `outer_sky_convergence` (16 iterations) with exactly the same
explicit `gaia_strict` or `tractor_psf_validated` policy and
unmodified input files, and export full non-null reports for both.
Compare exact training and holdout star IDs and fixed selection
fingerprints, not only the headline residual percentages. Keep this
branch DRAFT and the scientific adoption gate closed.


## New real Leo P i-band diagnostic: 2026-11-05 16:52–16:53 UTC

Two full-griz metadata exports were inspected directly:

- `16:52:40.531Z`: `i` records `outer_sky_robust`
  (`max_iterations=8`, `tractor_psf_validated`, 2 guarded training
  outlier passes), but `science_report_or_rejection=null`. The current
  UI policy was `dao_diagnostic`, contradicting the old
  `exact_result_recipe.star_catalog_policy`. This is a stale
  *provenance* bug in a reset/unfinished report, not an 8-iteration
  scientific result. Original `i`-band image/invvar/maskbits/nexp SHA256
  match the next export.
- `16:53:45.141Z`: the completed `i` result uses
  `outer_sky_convergence` (`max_iterations=16`),
  `tractor_psf_validated`, no train residual auto-trim,
  `PSFSIZE=5.12` guessed FWHM and a `65×65` pixel stamp. All 1427
  DAO detections were screened, 1402 rejected before the sample split:
  **22 training** (10 Gaia-matched + 12 Tractor-PSF without Gaia) plus
  **3 Gaia-matched independent holdouts**. The signed holdout integrity
  report was `COMPLETE` (3/3 stars with core, intermediate and outer
  annuli); none of their IDs was replaced based on fit residuals.
  The selected instrumental source quartet is unchanged from 16:52.

The second run is a useful numerical convergence milestone. Photutils
3.1.1 reached `converged=True` in **12/16** iterations,
`final_converged_fraction=1.0`, `max_center_shift=0.000301 px`;
the last `max_epsf_change` was **0.002732** of the ePSF peak.
These conditions establish the library's *centroid-based*
convergence, not wing fidelity or PSF applicability.

The PSF is still **`STRONG_NEGATIVE_EPSF_RINGING`**:
negative absolute mass **0.053711** versus blocking pilot 0.025.
Of the *total signed-kernel absolute mass*, 0% is negative in the
core, **0.1739%** is negative at 1–2 guessed FWHM and **5.1972%**
is negative outside 2 guessed FWHM; within that last outer region
**40.95% of the absolute ring mass itself is negative**.
The 3 Gaia holdouts have median per-pixel absolute signed-model
residual / stellar peak **0.009757 / 0.003151 / 0.002241** across
the three annuli. These medians do not constrain integrated outer
ringing: pixels with opposite residual signs can offset, the ring
contains thousands of near-sky pixels, coadd noise is correlated, and
the sample contains only 3 independently selected stars.

**Strong population/brightness confound:** all 12 Tractor-only
training sources are significantly fainter than the Gaia training
sources. The two groups have no overlapping measured-peak range:
Tractor-only peaks 0.0798–0.1804 nMgy/pixel versus Gaia 1.371–17.323.
The groups' median peak/RMS is **17.83 versus 484.95** and their
median positive native core moments are **3.182 versus 2.339 px**
(≈36% larger for the faint population). These are *not PSF FWHM*;
at this brightness disparity they cannot establish that the faint
population is nonstellar. Bright Gaia holdouts cannot directly
validate potential population-dependent biases introduced by the
faint Tractor-only training sample. Spatial sampling may also be
different. Do not automatically remove the faint sources based on
an unmatched width or on training-set residuals.

**Required controlled tests before deciding whether the sample or the
method is at fault:**

1. Run `outer_sky` (8) and `outer_sky_convergence` (16), both with
   `tractor_psf_validated`, the same original `i` native FITS,
   `5.12 px` FWHM guess, and no automatic training outlier passes.
   Both runs must have non-null full result QA; compare exact train/
   holdout Gaia IDs and exact effective recipe fields.
2. Repeat the same 8/16 pair with `gaia_strict`, retaining 10+3
   minima. This is a population ablation, **not** a direct fixed-star
   comparison between selection policies. Record selection underflow
   if it occurs; never lower the minimum.
3. Only if necessary, evaluate `wing_smoothing=False` as a **separate
   controlled ablation** at oversampling=1 with a fixed training
   population, because Photutils 3.1 documentation cautions wing
   smoothing may reduce wing fidelity in undersampled ePSFs.
   This is not evidence that smoothing caused the current rejection:
   the existing outer annulus (>2 FWHM) combines the inside/outside
   of the library's wing-smoothing radius.
4. Independently inspect rejected signed ePSF, per-star cutouts,
   2D residuals, annular flux integrals, sky and field-position
   variation. Keep PSF0 non-adopted; do not weaken the 2.5% screening
   heuristic to manufacture a pass.

**Export provenance repair:** an old `_last_reported_payload` could
survive a filter/catalog policy/experiment reset, creating a non-null
`exact_result_recipe` beside `science_report_or_rejection=null`.
The new exporter records `result_provenance_state`, clears stale
recipe state upon reset, and refuses to attach old source hashes or
effective recipe when no scientific result exists. This repairs
report integrity, without modifying the PSF-fitting algorithm, saved
FITS or historical QA files.

Reference: official Photutils 3.1 ePSF builder documentation:
https://docs.astropy.org/projects/photutils/en/stable/user_guide/epsf_building.html
