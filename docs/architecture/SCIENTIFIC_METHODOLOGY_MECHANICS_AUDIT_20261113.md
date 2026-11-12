# Desidera: end-to-end scientific methodology and mechanism reconciliation

**Audited revision:** `75334183cdcd440933e0dd241af1322121c59246`, draft PR #84 (`feature/scene-wcs-editor-sc1`).  
**Date:** 2026-11-13.  
**Assessment type:** code, scientific-contract, documented Leo P QA, and external methodological comparison. **No new original-FITS reduction or independently reproduced Leo P pixel measurement was performed in this audit.**
**Scope:** current optical DESI Legacy Surveys DR9–DR11 coadd pipeline, conditional PSF, Scene, sky, direct aperture photometry, covariance diagnostic, staged GUI and future science gates. Existing specifications are *intended requirements*, not evidence that every named algorithm has been implemented.

## Executive verdict

1. **Engineering/data-lineage validation: substantially implemented and synthetically tested.** Original quartets, IMAGE/INVVAR native units, optical MASKBITS, celestial WCS, native-source hashes, isolated project/draft identity, exact aperture overlaps, hard pixel validity and cross-stage QA gates exist. SHA-256 is a local byte-identity check, not proof of archive provenance or photometric quality.
2. **Current astronomical measurement: signed valid-area *finite-aperture* photometry, diagnostic only.** It is not a total, a decontaminated galaxy flux or an isophotal profile. Its formal variance includes only pixel-diagonal inverse variances, and it presently uses zero additional residual-sky correction. Negative pixel and aperture fluxes are retained.
3. **Main scientific failure modes are upstream/inferential, not FITS parsing:** imperfect initial Legacy sky removal, poorly identified foreground-versus-target light, finite/space-varying PSF wings, correlated/coherent sky noise, missing areas, choice of photometric geometry and extrapolated outskirts. Passing the implemented synthetic/GUI test matrix does not bound these errors.
4. **No production science acceptance:** the `SCIENCE HOLD` status is appropriate. Do not publish total magnitudes, $R_e$, colours, SBF distances or physical stellar-population parameters from the current Stage 04 product.

## Implemented mechanics versus scientific estimand

| Layer / implementation | Mechanic verified in code | Scientific interpretation supported | Missing prerequisite / verdict |
|---|---|---|---|
| Object project / data manager; `surveys/{object_projects,dataset_registry,dataset_integrity}.py` | Project-bound native source, local index, complete coadd quartet, header/WCS validation, source SHA, explicit band/release | A reproducible *candidate* optical input and sampled sky footprint | Neither all CCD exposures nor an independent photometric validation; **engineering usable** |
| P0.5 `photometry/{native_fits,contracts}.py` | Original IMAGE, INVVAR, NEXP and MASKBITS are read-only, native-grid, in nMgy/pixel and nMgy^-2, with hard invalid-pixel gate | Signed flux-bearing native samples with documented metadata | INVVAR lacks off-diagonal resampling covariance; no absolute/foreground-extinction or PSF correction; **input usable** |
| PSF-0 `psf_{gui_backend,empirical,holdout_qa,core_wings_validation,...}.py` | DAO/Gaia/Tractor-assisted candidate selection, star-cutout local sky, native ePSF construction, held-out stellar core/ring diagnostics | A **finite empirical PSF candidate** conditional on stellar sample and local sky | Specific GUI backend is **Leo P-only**; far wings, brightness/position dependence, total encircled energy and galaxy-position transfer uncalibrated; **HOLD** |
| Scene `scene_contaminants.py`, `sky_masks.py`, `scene_source_inventory.py` | Exact native raster masks, semantic target-protection and manual exclusion, catalogue screening, Gaia/SGA labels and fixed-source inventory | Provisional contaminant screening; target protection is **not** bad-pixel masking | Tractor classification does not establish foreground membership; bounding radii are heuristics (not wing flux predictions); H II knots/resolved members can be lost or spurious sources retained; **HOLD** |
| SA2a `scene_source_association.py` | Native-band detections provisionally linked by angular distance; nontransitive/multi-source ambiguity flagged | Positional candidate group, not an astrophysical identity or membership assignment | Relative astrometry/crowding/independent segmentation completeness/purity unvalidated; **HOLD** |
| Preliminary isophote `ellipse_seed.py` | Mask-dependent Gaussian-smoothed contour, SGA-seeded morphology, diagnostic perturbation/sensitivity | Geometry *hypothesis* and on-screen contour only | Not a fitted variable isophote, validated outer extent, asymptotic radius, measured $R_e$ or target/foreground classifier; **HOLD** |
| Stage 03 `background.py` | Explicitly selected, disjoint spatial training/holdout tiles; robust constant and tilted-plane residual-sky candidate fits; fixed-mask exponential injection | A *conditional additional sky level* and within-one-coadd diagnostic | Prior Legacy subtraction may already have removed galaxy light; tiles can contain unrecognized outskirts/scattered star wings, have unequal spatial sampling and correlated noise; no adopted sky or calibrated posterior; **HOLD** |
| Stage 04.1 `aperture_flux.py` | `photutils` exact fractional overlap, source/manual/hard-area ledger, signed $F_{\rm valid}$, positive-only AB magnitude, $w^2/\mathrm{INVVAR}$ variance | Observed valid-area flux in a **declared finite ICRS ellipse** | Not galaxy-only, missing-light-corrected, PSF-matched, residual-sky-corrected or total; reported $\sigma_{\rm diag}$ is not complete error; **diagnostic MEASURED_PARTIAL** |
| Stage 04.2 `aperture_growth.py` | Fixed set of nested reference ellipses; direct cumulative and independently integrated annular sums; closure check | Conditional valid-area curve of growth over a **predeclared grid** | Radii depend on provisional protected area + image margin; one fixed geometry, no fitted outer isophotes, plateau significance or reliable extrapolation; **diagnostic** |
| Stage 04.3 `aperture_systematics.py` | Analytic derivatives to constant/plane sky, source-mask counterfactual, geometry span, spatial blank-ring covariance (when available) | Conditional sensitivity and local same-coadd empirical scatter **only** | Blank controls are not independent exposures; full scale may lack enough placements, sky structures remain, inverse covariance not calibrated; **diagnostic** |
| Old `profiles.py`, `sensitivity.py`, `simulation.py` | Fixed pixel-scale prototype profiles and simple PSF/noise-free Sérsic mock | Useful numerical unit-test tools, **not** an adopted native-WCS isophotal/total pipeline | Avoid conflating them with the Stage 04 operator or end-to-end injection/recovery |
| Optical multi-band, SED, SBF, CCD-hybrid | Scientific methods and contracts specified in docs/notes | **Planned** quantities only | No validated matched-PSF colours, end-to-end CCD-level coadd redo, SBF likelihood/calibration, or science-adopted galaxy total in the current GUI |

### Methodology is not yet equivalent to the intended eight-stage P1.3 sequence

`PSF_SCENE_SKY_ORDER_P13.md` prescribes Data Preparation → provisional Scene → preliminary star-cutout sky → PSF → final PSF-aware Scene → final sky → joint refinement → photometry. The **four visible GUI tabs** implement overlapping exploratory portions (PSF-0, Scene, background, direct aperture) but are not evidence that all eight *scientific operators* execute.

In particular Stage 04 deliberately **does not apply** the Stage 01 ePSF or Stage 03 fitted residual sky; zero additional correction is hard-coded in the integrated GUI Stage 04 worker (`gui/aperture_workspace.py`). This is a defensible non-adoption policy while the models remain unvalidated, and a decisive limit on physical interpretation.

## Mathematical invariants that must govern implementation

**Forward model (coadd schematic only):**

$$
I_b(p) =
 \left[ G_b \otimes P_{b,\mathrm{eff}} \right](p) +
 \sum_j \left[S_{j,b}\otimes P_{j,b,\mathrm{eff}}\right](p)
 + B_{b,\mathrm{res}}(p)+\epsilon_b(p).
$$

The already-subtracted survey sky component is not directly observable from the sky-subtracted coadd. No unique solution for true extended galaxy light follows merely by fitting a new constant or plane. Individual CCD modeling can provide better leverage, but still needs independently tested sky/PSF/galaxy identifiability.

**Finite masked-aperture flux:**

$$
F_{i,b}^{\mathrm{valid}}
=\sum_p w_{i,b,p}\,v_{b,p}[I_b(p)-B_b(p)],
\qquad 0\le w_{i,b,p}\le 1.
$$

`v` is hard validity AND NOT provisional-source/manual-excluded, not a probabilistic galaxy membership. A valid-area mean $F_{\rm valid}/A_{\rm valid}$ is not automatically an unbiased estimate of the **whole** elliptical annulus if masked regions preferentially contain knots or foreground stars. Negative signed flux remains legitimate.

**Covariance:**

$$
\operatorname{Cov}(F_i,F_j)
=\mathbf w_i^\mathsf T C_{\rm pix}\mathbf w_j
+\mathbf d_i^\mathsf T\Sigma_{\rm sky}\mathbf d_j
+\text{PSF/membership/model/selection terms and cross-terms}.
$$

Here $\mathbf d_i$ is the signed flux derivative with respect to sky-plane parameters. The implemented formal aperture variance $\sum_p w_p^2/\mathrm{INVVAR}_p$ accounts only for the **diagonal** pixel term; a resampled Legacy coadd has nonzero pixel correlations. The current blank-ring `np.cov` is a useful empirical same-image diagnostic, not $\Sigma_{\rm sky}$ or a universally invertible science covariance. Nested cumulative apertures share pixels, so independent-point curve-of-growth least squares is not justified.

**Galaxy total:**

$$
F_{b,\mathrm{total}}\equiv F_{b,\mathrm{observed}}
+\Delta F_{b,\mathrm{masked}}
+\Delta F_{b,\mathrm{outside\ finite\ domain}},
$$

with an explicitly defined target-component convention. The last two terms cannot be set to zero because an observed curve is flat. Statistical and model uncertainty must be propagated; the current Stage 04 only supplies the first term for surviving pixels.

**Matched-aperture colour:** same ICRS sky domain and consistent galaxy component, with explicitly matched effective PSF or joint forward model and joint $g/r/i/z$ covariance. Merely using the same angular radius on dissimilar native coadd PSFs does not certify colour. Separate per-band asymptotic radii, common-aperture colours and total colours are different estimands.

## Independent Leo P findings already documented, not independently reproduced here

Existing researcher-supplied Leo P DR11 South r-band QA, written to `PSF_SCENE_SKY_ORDER_P13.md`, reported 24/7 training/holdout spatial tiles and a **frozen-mask**, post-Legacy-sky-subtraction exponential injection. The conditional recovered injected fractions are:

| Exponential scale length relative to provisional target-protection area-equivalent radius | Constant additional sky | Tilted additional sky |
|---:|---:|---:|
| 0.5 | 0.966 | 0.949 |
| 1.0 | 0.839 | 0.781 |
| 2.0 | 0.612 | 0.537 |

These numbers diagnose **additional local sky fitting absorbing artificially injected galaxy light**. They are not empirical loss fractions for Leo P, not uncertainty estimates, and do not test historical DR11 sky loss, a re-run of segmentation, a convolved PSF or foreground subtraction. This nevertheless directly falsifies any unconditional rule that a clean holdout RMS guarantees preservation of LSB outskirts.

In `SCENE_LOGIC_AUDIT_20261110.md`, the researcher-supplied P21/2 r-band source/sector experiment for bright Gaia star / Tractor OBJID 2035 found **zero of five** radius bins satisfying the predeclared sector-coverage quality criteria, and a low-confidence galaxy contour overlapped a large fraction of source-facing distant rings. Apparent anisotropic excess cannot therefore be attributed uniquely to the star's far PSF wing or to Leo P's outer galaxy light. The finite ePSF stamp did not reach the region of concern. P21/3 conditional segmentation comparison still requires a new original-project GUI run and cannot itself establish the physical decomposition.

## External scientific-methodology cross-check

- **Legacy Surveys DR11 data description:** coadds are inverse-variance-weighted and resampled; a survey background has already been subtracted (including spatial spline correction); Tractor operates on individual exposures using spatially dependent PSFEx models. Thus a single uniform coadd kernel and a post-coadd residual plane are approximations, not their joint CCD forward model. https://www.legacysurvey.org/dr11/description/ and https://www.legacysurvey.org/dr11/files/
- **Legacy DR11 bitmasks:** `GALAXY` and `RESOLVED` are source-context flags, not default exclusions of intrinsic dwarf-galaxy light. `TYPE=PSF` in masked regions is not stellar-membership proof; Gaia-retained entries can have zero fitted band inverse variance. https://www.legacysurvey.org/dr11/bitmasks/ and https://www.legacysurvey.org/dr11/description/
- **Watkins et al. (2024):** diffuse unrecognized source flux and extended stellar scattered light bias modelled sky levels. Local sky holdout and PSF-star annuli alone cannot certify the true broad-band background. https://arxiv.org/abs/2401.12297
- **Stone et al. (2021), AutoProf:** published machinery includes sky, centre, azimuthal ellipse fitting, radial profiles, direct cumulative measurements and fit checks. Desidera currently has geometry *proposals* plus fixed exact-aperture diagnostics; do not label them already equivalent to a full AutoProf-like profile. https://doi.org/10.1093/mnras/stab2709
- **Moustakas et al. (2023), SGA-2020:** custom coadds, other-source Tractor models, fixed isophote profiles and explicit CoG extrapolation demonstrate a mature comparison method. **Do not use the affected SGA-2020 `[G,R,Z]_COG_MAG`/`COG_MAGERR` columns as uncritical external truth**: the release documentation explicitly warns of a coding bug, with some other derived surface-brightness quantities remeasured from unaffected profiles. https://www.legacysurvey.org/sga/sga2020/
- **SBF:** normalized galaxy-residual Fourier power depends on the PSF *and* image/mask transfer; resampled correlated coadd noise need not be white. No fitted SBF power spectrum, contamination variance, suitable colour-calibration applicability or validated SBF distance currently exists. Example methodological comparison: https://academic.oup.com/mnras/article/527/2/1656/7326767

## Prioritized scientific defects / uncertainties

**P0: block scientific adoption until empirically bounded**
1. **Upstream sky/light loss and residual-sky degeneracy.** Predefine multiple sky models, nearby outer-light scales, and observational control; use injection at both final coadd and as far upstream in CCD/reduction as feasible. Separately quantify "cannot recover from current input" regimes.
2. **Foreground-membership / target-light preservation.** Audit bright Gaia sources (Leo P OBJID 2035 first), resolved stars, H II knots, red extended background objects and overlapping masks; compare mutually exclusive plausible identity treatments at fixed geometry. A single colour/TYPE threshold is not an astrophysical classifier.
3. **PSF wing/spatial transfer.** Validate core, intermediate and wide-angle wings with *disjoint* star cohorts, stellar brightness ranges, subpixel phase, position and local sky nuisance; test impact on recovered galaxy flux. Coadd finite-ePSF normalization is not infinite-wing normalization.
4. **Full multi-scale covariance.** Derive aperture/ring-dependent empirical noise on truly appropriate controls, model coherent residual sky covariance, include cross-band/geometry dependencies, preserve rank checks, do not extrapolate a small-scale blank covariance to a large galaxy aperture.
5. **Galaxy observables/uncertainty calibration.** Build actual isophote and structural model alternatives, missing-light and outer-flux bounds, selection-aware end-to-end injection/recovery. Do not adopt a model total because a fitted profile has a plateau.

**P1: capability needed for general surveys, not evidence of P0 acceptance**
6. Generalize the **Leo P-specific** PSF GUI/backend to arbitrary project-target identity without making relaxed defaults.
7. Add explicit two-level products: *coadd diagnostic* and *CCD-aware science inference*; treat DR9/DR10/DR11 derived from reused CCD frames as correlated, never independent measurements by release alone.
8. Implement multi-band common native-sky geometry, band-specific WCS/PSF, matched-pixel/PSF colour likelihood, and separate per-band total estimates. Do not force an identical adopted outer radius or blue outskirts to disappear.
9. Introduce an immutable, typed scientific output schema and predeclared adoption policy for a **specific estimand**, not a global "PSF passed" toggle. Keep dependencies, masks, observational source and selection decisions versioned.
10. Leave WISE/GALEX population inference and SBF **downstream** of the optical total/colour and SBF-specific validation gates.

## Minimal preregistered end-to-end acceptance experiment

Choose a small but stratified set of *independent* nearby dwarfs and blank/sparse controls, retaining Leo P as a high-risk case, not the sole calibrator. Predeclare the science claim (e.g. finite matched-aperture flux vs fully model-corrected total), target surface brightness, angular size, morphologies, survey footprints, crowding and tolerance on bias/coverage.

Run the **entire configured chain**, rather than reusing masks/centre/sky/radius chosen on the unmodified image:

1. Independent synthetic images with known truth, nonlinear WCS/RA wrap and area-integration controls; numerical agreement of signed pixels, WCS, annular closure and IVAR diagonal errors.
2. Synthetic galaxies inserted at controlled positions into **copies of real native coadds**, covering several exponential/Sérsic indices, irregular/clumpy/resolved systems, colour gradients, orientation and effective radii, multiple brightness levels, scattered PSF light and several sky backgrounds. Re-run source detection, foreground classification, PSF, sky, geometry, aperture selection and estimation. Report missed and wrongly assigned target light, not just surviving-pixel flux.
3. Repeat, as data allow, with injected simulated galaxy signal into **individual CCD exposures before survey sky subtraction**; propagate real per-CCD PSFs and sampling through the coadd. Quantify differential loss relative to post-coadd injection. If upstream data inaccessible, record the resulting **unbounded upstream loss** instead of inferring a correction.
4. Hold out galaxies and stellar cohorts from choice of model/threshold. Recover bias and scatter as a function of true structural parameters; audit nominal 68%/95% interval coverage, sky/PSF/membership contamination and false adoption rate.
5. Compare with genuinely independent deeper imaging or validated survey measurements, taking care not to ingest the known faulty SGA-2020 CoG columns. Sensitivity agreement of two fits to the same coadd is not an independent check.

**Acceptance rule:** remain `DIAGNOSTIC_ONLY`/`MEASURED_PARTIAL`/`MODEL_CONDITIONAL`/`UNDERCONSTRAINED` until Q0–Q5 in `PHOTOMETRY_QUALITY_VALIDATION.md` have *measured* results under a frozen object-class-specific policy. Only an explicitly reviewed Q6 run may be `TOTAL_VALIDATED`. No thresholds are chosen retroactively from Leo P's own outputs.

## Audit provenance and limitations

Repository review inspected the live GitHub SHA above, including `surveys/dataset_integrity.py`, `photometry/native_fits.py`, `photometry/contracts.py`, `photometry/scene_{contaminants,source_inventory,source_association}.py`, `photometry/psf_{empirical,joint_scene,core_wings_validation}.py`, `background.py`, `photometry/aperture_{flux,growth,systematics}.py`, `gui/{scene_editor,aperture_workspace,photometry_workspace,project_qa_bundle}.py`, `profiles.py`, `simulation.py` and the existing architecture/QA documents.

The most recent passing synthetic CI suite validates engineering behavior *as encoded*, not the astrophysical accuracy of those choices. No remote original Leo P FITS were fetched or newly reduced, no real-galaxy remeasurement was repeated and no external science benchmark was executed for this audit. The present report intentionally makes **no** declaration of a science-ready PSF, total magnitude, galaxy colour, $R_e$ or SBF distance.


---

## Follow-on forensic audit — current Stage 01–04 mechanics and concrete remediation (2026-11-13)

**Purpose:** identify *specific changes to the running stages*, not restate the future architecture. Examined at code revision `c342ac01121c00f621a426e8c82c8acfda4902ba`. All proposals below are explicitly separated into: **CONFIRMED DEFECT** (a demonstrable discrepancy between stated and calculated quantity), **OBSERVED/REPORTED LIMIT** (source code plus recorded QA, without independent FITS rerun), and **PROPOSED TEST** (not yet shown to improve real data). The term "P0" here is an implementation/scientific-validation priority, **not** a scientific status or claim of a universally valid threshold.

### F01. Stage 04.3 sky derivative used the wrong parameter pivot — CONFIRMED DEFECT, FIXED (P0)

**Location:** `photometry/aperture_systematics.py::diagnose_aperture_systematics`, versus `aperture_flux.py::measure_native_aperture`.

For the stored model `b0 + bx*(x-xref) + by*(y-yref)`, the published gradient derivative must be `-sum(w_usable*(x-xref))`. Stage 04.3 formerly reported `-sum(w_usable*(x-x_aperture))` regardless of `sky.reference_xy_pix`, even though its field name implied an exact derivative of the *frozen* scenario. The difference is `-A_valid*(x_aperture-xref)` for `bx`, similarly for `by`. It can be large for a large native elliptical aperture and a sky plane anchored at the cutout rather than the galaxy.

**Remedy committed:** `c342ac0` makes `reference_native_xy` the actual frozen sky-scenario reference, uses it for `dflux_dbx_native_pixel3` and `dflux_dby_native_pixel3`, and preserves differently named target-centred derivatives as alternative *parameterizations*. An independent regression perturbs all three sky-plane parameters with an off-centre reference and checks against new direct exact native aperture integration.

**Gate:** finite-difference or algebraic sky response must match the *same* sky-parameterized flux operator for asymmetric masks, nonzero gradients, non-central `(xref,yref)`, several WCS handednesses and radii. An internally self-consistent derivative is still **not** an estimate of sky-parameter uncertainty.

### F02. Stage 03 target-protection footprint is insufficient to establish galaxy-free sky — OBSERVED/REPORTED LIMIT (P0)

**Location:** `background.py::_automatic_scene_sky_pixels` lines around target veto and fixed tile assignment; `_audit_current_project_sky_injection_recovery`.

The target veto expands a provisional mask by an area-equivalent radius and pixel-distance guard. Neither procedure constrains the unknown, potentially asymmetric low-surface-brightness outskirts. It can protect all currently *detected* pixels yet admit real faint halo/stream light into the sky training sample. The current frozen-mask injection freezes source detection, geometry, masks and sky tiles and uses an **unconvolved** ideal exponential injected after the original Legacy reduction; by design it cannot calibrate the full procedure.

**Existing Leo P r researcher-supplied evidence:** for increasing injected $h/R_\mathrm{protect}=0.5,1,2$, the fitted-sky recovered fractions in the frozen unmasked aperture were about 0.966/0.839/0.612 (constant) and 0.949/0.781/0.537 (plane). These numbers are **additional-sky-fitting response to artificial signal**, not the true lost fraction of Leo P light. The automatic no-subtraction recommendation is therefore defensible.

**Improve:** make a *protected-radius and asymmetric-footprint ensemble* a first-class, separately versioned Stage 03 sensitivity product. Compare several target-protection shapes, bright-star wing veto alternatives, sky-grid offsets and one genuinely held-out outer field; report $b_\mathrm{target}$, full-aperture $A_\mathrm{valid}\Delta b$, tile-by-tile residual medians and injected-signal leakage into the exact selected sky pixels. Do not select a model using only a lower same-coadd RMS. If independent sky cannot be found, surface `SKY_UNIDENTIFIABLE` rather than silently subtracting.

**Acceptance test:** on synthetic extended targets and independent observed controls, evaluate distribution of recovered flux conditioned on *retained/rejected* galaxy extensions after **re-running** the mask/sky selector, not merely frozen-mask sky fitting. Compare constant, plane, zero correction with a policy frozen before seeing target results. No fixed percentage is science-approved yet.

### F03. Stage 03 train/holdout tile design has a spatial-selection and correlated-sampling weakness — CONFIRMED CODE POLICY / PROPOSED IMPROVEMENT (P1)

**Location:** `background.py::_automatic_scene_sky_pixels`.

The source-aware tile layout is geometry-deterministic, with `(tx+2*ty)%5==0` for holdout and at most 384 pixels selected via `np.linspace` on the row-major list of valid indices. Equal tile *counts* do not ensure equal radial/azimuthal coverage about the galaxy, matched contamination distributions or statistical independence at the coadd's correlation length. The 384-index linear selection can exhibit a regular spatial lattice and sensitivity to chip/coadd patterns. This is a *possible bias mechanism*, **not** a proven Leo P numerical error. Reported Leo P tile distributions were spatially unbalanced and training RMS differed strongly from holdout.

**Improve:** predeclare several *spatially stratified folds*, balanced by target radius, quadrant and coverage, not by measured image brightness. Record sampled positions and their 2D pair-separation histogram, tile completeness, clipped edges, per-fold sky shifts and target-centre response. Compare row-major, seeded 2D stratified and spatial-block sampling **only on the same fixed eligible geometry**, then test out of sample. Version the sky sampling recipe: changing it changes the sample and invalidates downstream evidence.

**Acceptance test:** synthetic stationary white/correlated fields, gradients and galaxy-outskirts contamination; require demonstrably unbiased constant/plane recovery and calibrated *block-level* uncertainty before assigning any scientific preference.

### F04. Stage 01 empirical PSF star-cutout sky and PSF wings are non-identifiable from the same ring — CONFIRMED ALGORITHM / PHYSICAL HOLD (P0)

**Location:** `psf_empirical.py::_star_cutout` and `PSF0Config`; `psf_survey_defaults.py`.

The per-star local sky uses a clipped planar fit in an annulus around each selected star. A real extended wing, PSF spatial variation, close companion or broad halo can contribute coherent signal that is not removed by sigma clipping and is then subtracted from the ePSF stamp. A guessed FWHM or DR11 PSFSIZE is a *seeing-scale prior*, not a calibrated encircled-energy function; a normalized finite kernel cannot predict light beyond its cutout. The current kernel is conditionally evaluated on held-out stars and its reported `psf_wings_validated` is **false**; that is the correct semantics.

**Improve:** for each native band, retain at least two stellar sky-radius hypotheses with equal fixed star cohort; compare radial flux fractions and local background planes on genuinely held-out stars in several spatial/brightness strata. Measure intermediate and wide-angle wings in suitably sized original-star fields, with neighbour modelling, and track the **finite-PSF normalization radius** explicitly. Reject a "global" ePSF as representative if spatial/brightness holdouts do not support transfer to the galaxy position. The official Photutils EPSFBuilder guide warns both about star cleanliness and about wings biasing a local annular background: https://docs.astropy.org/projects/photutils/en/stable/user_guide/epsf_building.html .

**Acceptance test:** known injected stars with varied seeing, phase, brightness and wing profiles; recover a galaxy+star composite at the galaxy's relevant limiting surface brightness. Do **not** make a global far-wing claim merely from low heldout core residual.

### F05. Stage 01 PSF holdout can share sample-selection information with the training cohort — CODE MECHANISM / VALIDATION RISK (P1)

**Location:** `psf_empirical.py` near `accepted`, `field_relative_core_width` and subsequent deterministic split.

The reference population of Gaia-matched stellar **core widths** and other preliminary quality screens is computed over candidates before the final training/holdout split; the held-out stars are therefore independent of the final ePSF **fit**, but not necessarily independent of *all empirical decisions* determining the final eligible population. This does not invalidate every residual check, but it limits interpretation as a fully blinded accuracy estimate.

**Improve:** freeze quality-cut thresholds from a separate training/calibration star cohort (or from survey-independent criteria), then evaluate a disjoint untouched *outer* validation cohort. Keep one internal training split for EPSF iteration and a separate external holdout for stellar population / transfer calibration; record actual Gaia IDs and spatial locations. Report holdout completeness by angular radius and magnitude rather than just count.

**Acceptance test:** compare internal holdout vs genuinely disjoint external validation on known star simulations, measuring centroid, encircled energy, core-to-wing residual and spatial transfer bias. No ad hoc reduction of the minimum star count to manufacture a pass.

### F06. Stage 02 catalogue masks do not implement a physical source decomposition — CONFIRMED CODE POLICY / PHYSICAL HOLD (P0)

**Location:** `scene_contaminants.py::conservative_catalogue_mask`, `scene_source_inventory.py`, `scene_source_association.py`.

The automatic source envelope uses magnitude/Tractor-shape-dependent *heuristic* radii. Even after the geometrical WCS fixes, it neither models PSF wings nor decides which unresolved or compact sources belong to the dwarf galaxy. An interior Gaia object can be foreground, an unresolved target member, a blend or an erroneous crossmatch. The optional red-extended/inner-Gaia exclusion rules are diagnostic and **not** proven classifiers. A positional SA2a match likewise does not establish astrophysical identity.

**Improve:** create an explicit *per-source competing-treatment ledger* (unmasked; excluded but *no* flux reconstruction; forward PSF/Tractor model-subtracted; target compact member retained) and a WCS-identical image-difference comparison. Use group-matched Gaia astrometry with quality controls, individual CCD/residual shape, multi-band colour and external resolved-star data where available; no automatic membership probability until calibrated. Keep hard invalidity separate from source policy and always quantify the light lost from each treatment using artificial-galaxy injection.

**Acceptance test:** inject partially overlapping foreground PSFs, blue HII knots, resolved upper RGB/AGB stars and diffuse red background galaxies into native coadd copies; report false removals, false retained contaminants, integrated galaxy flux bias and identified ambiguous regions. Repeat the actual source finder, not just fixed masks.

### F07. Leo P foreground source OBJID 2035 is a distinct identifiability test, not yet a PSF measurement — OBSERVED/REPORTED LIMIT (P0)

**Location:** `psf_joint_scene.py`, P21 directional/sector QA and `SCENE_LOGIC_AUDIT_20261110.md`.

Existing user-supplied P21/2 QA reports that **all five** directional rings fail coverage criteria and part of the source-facing signal overlaps a provisional low-confidence galaxy contour. A positive facing-versus-away contrast cannot be uniquely interpreted as the star's PSF wing or Leo P light; same-image "blank" controls are not independent stars or sky realizations.

**Improve:** evaluate the existing P21/3 provisional-*component*-aware split on a fresh exact user-project QA export (not the padded display contour), then compare four *conditional* scenarios at frozen sky/geometry: star kept; star masked; empirical PSF with an independently constrained broad tail; separate star+galaxy+sky forward fit. No star model or mask may be adopted on this comparison alone.

**Acceptance test:** held-out bright stars of matched magnitude/environment and injected artificial dwarf light around them; calibrate the false stellar-wing assignment rate before using subtraction in Leo P.

### F08. Stage 02 preliminary isophote is an anchored smoothed-component detector, not a fitted light profile — CONFIRMED MECHANIC / IMPROVEMENT (P1)

**Location:** `ellipse_seed.py::propose_local_isophote`, `sga_informed_radial_seed`.

The algorithm derives a provisional scalar sky from a distant cutout annulus, convolves the masked image, thresholds by the *empirical smoothed-annulus scatter*, finds a connected region near the target and uses intensity-truncated image moments for one $q$/PA suggestion. It does **not** fit $q(a)$, PA$(a)$, centre drift, ellipticity or azimuthal Fourier harmonics; a contaminated sky annulus and bridged neighbouring structure can shift the contour. A small trial-to-trial perturbation does not supply a confidence interval.

**Improve:** separate detection mask from *isophote fitting*; introduce a controlled radial solver in original native sky coordinates with explicit exposure/invalid fraction, coherent PSF-aware contamination alternatives, constrained centre and smooth $q(a)$ and PA$(a)$, recorded covariance. Keep a freely varying centre as a *separate hypothesis* to detect drift and blends. Validate on irregular blue knots, inclined diffuse disks and low-S/N outer isophotes; measure angular sector scatter and low Fourier modes, not just global moments. Compare against AutoProf on the **same** sky/mask/PSF inputs; agreement of profiles alone is not truth: https://arxiv.org/abs/2106.13809 .

### F09. Stage 04 direct flux is a partial valid-area estimator; missing-light correction is not available — CORRECT CODE SEMANTICS / P0 SCIENCE LIMIT

**Location:** `aperture_flux.py::measure_native_aperture`, `aperture_growth.py`.

The implementation correctly integrates exact fractional overlap with signed pixels, respects hard invalidity and manual/source exclusions, and does not multiply by the geometric-to-valid area ratio. A masked area can contain a diffuse central galaxy gradient or HII knots, so **neither zero replacement nor simple area scaling** is an unbiased galaxy flux reconstruction. A positive-only AB magnitude of surviving pixels cannot be renamed galaxy total. Current outputs are appropriately labelled MEASURED_PARTIAL/DIAGNOSTIC.

**Improve:** publish distinct `F_valid`, `F_model_in_mask`, `F_extrapolated_outside` and `F_total_conditional` in an eventually adopted schema with uncertainty for each. Until then report a *quantitative lower/upper model-conditional envelope* of missed flux under alternative source/sky/geometry models, without silently modifying `F_valid`.

**Acceptance test:** synthetic nonaxisymmetric galaxies and bright superposed stars; by mask policy, compare exact injected true flux inside known geometric/missing domains and distribution of reconstructed errors. Do not simply test ring-vs-curve-of-growth numerical closure.

### F10. Stage 04 sky response must be radial and tied to covariance — IMPROVEMENT TO CURRENT DIAGNOSTIC (P1)

**Location:** `aperture_systematics.py::diagnose_aperture_systematics`.

After F01, the largest-aperture derivative is exact for the frozen sky reference, but the routine still reports sensitivities **only for the largest radius and baseline geometry**. A tiny coherent gradient may create a false outer curve-of-growth plateau, change ring signs and bias any $R_e$ selected from the same noisy curve. A single derivative and an uncalibrated same-coadd scatter cannot quantify a radial uncertainty envelope.

**Improve:** for *each* frozen radius $a_i$, calculate the exact area and first two *usable weighted pixel moments* about the **same sky pivot**, then publish the Jacobian matrix $D_{i\alpha}=-\sum_p w_{ip}v_p(1,x-x_\mathrm{ref},y-y_\mathrm{ref})$. Show the covariance-induced CoG envelope $D\Sigma_\mathrm{sky}D^T$ **only if** a separately calibrated $\Sigma_\mathrm{sky}$ exists; otherwise show signed response to predeclared sky perturbations and mark no confidence interpretation. Derive ring responses as differences of cumulative Jacobians, subject to the same masks and geometry.

**Acceptance test:** finite differences of independently integrated flux at every radius for asymmetric and wedge-masked apertures; response under sky-pivot changes must transform algebraically with the coefficients and leave the *physical sky image* and measured flux unchanged. Failure at any radius blocks use of a sky-error bar.

### F11. Stage 04.3 blank-ring covariance cannot be promoted to galaxy-scale correlated-noise error — CONFIRMED PROTOCOL LIMIT (P0)

**Location:** `aperture_systematics.py::estimate_blank_annulus_covariance`.

The routine measures flux scatter from geometrically non-overlapping, source-masked, fixed-size annular placements in one resampled coadd; it can run on smaller scales if the full target aperture does not fit. Large-scale sky patterns, unmasked faint sources and angularly anisotropic correlations may remain. Spatially separated positions in *one exposure/coadd* are not independent noise realizations, and the number of available large blank placements can be zero. Rank/PSD of a sample covariance do not prove astrophysical completeness or accurate uncertainties.

**Improve:** add a scale-by-scale block-control ledger with number of accepted spatial locations, all-excluded reasons, sky-area, exposure/NEXP structure, resampling correlation length, mean offsets and independent/overlapping CCD evidence where available. Do **not** extrapolate a 2–5 arcsec blank sigma to a 100 arcsec galaxy aperture by $\sqrt{N_\mathrm{pix}}$. Until independently calibrated, keep covariance only as a diagnostic and do not invert it for a formal chi-squared isophotal fit.

**Acceptance test:** inject known correlated Gaussian fields *and* coherent large-scale sky modes in synthetic original coadds, compare predicted vs measured sample covariance at each relevant aperture radius and verify holdout coverage intervals. Report failure when the full radial domain lacks uncontaminated blank positions.

### F12. Upstream sky subtraction is an input limitation, not recoverable by a post-coadd constant fit — EXTERNALLY CONFIRMED SCIENCE LIMIT (P0)

The DESI Legacy data description documents spline-based sky estimation after source masking, with medians computed in sliding **512-pixel boxes**, and confirms that delivered coadds already have that level removed. The DR11 file specification states that the native image is a Lanczos-3 resampled inverse-variance-weighted coadd in nMgy/pixel, whereas Tractor models individual exposures. This warns that a post-coadd test cannot quantify flux irreversibly lost to *upstream* sky subtraction or reconstruct the complete original covariance.

**Improve:** maintain explicitly separate coadd-only and CCD-aware processing modes. For difficult LSB objects, compare custom per-exposure background/PSF reduction from photometrically calibrated CCD files, not a local re-fit interpreted as raw-sky reconstruction. This is consistent with the independently published Stellar Stream Legacy Survey's CCD-level processing strategy. Do not treat a nominal 512 native-*coadd*-pixel angular scale as identical to an upstream CCD-pixel scale across instruments.

**Acceptance test:** same synthetic galaxy injected (A) after standard coadd sky subtraction and (B) into original calibrated CCDs before sky modeling/coaddition, with matched truth/PSF/filter response, then the same downstream Scene→sky→photometry chain. A difference isolates part of the upstream loss. If source CCDs unavailable, identify the loss as *unbounded from the current coadd*, not zero.

Official sources: https://www.legacysurvey.org/dr11/description/ ; https://www.legacysurvey.org/dr11/files/ . Independent treatment: https://doi.org/10.1051/0004-6361/202245011 .

### F13. Multiband matched colour and PSF transfer remain absent science operators — DOCUMENTED CAPABILITY GAP (P1)

The current same-brick band selector does not jointly fit per-band backgrounds/PSFs and does not make the per-band finite-aperture profiles share a validated target-light definition, compatible angular footprint and noise covariance. An identical ellipse in sky coordinates is **necessary but insufficient** for a physically comparable colour if PSFs and masks differ. A deep $r$-band geometry may also truncate blue outskirts visible only in $g$.

**Improve:** define a union-of-band candidate target segmentation and common ICRS geometry; compare source memberships and native valid coverage **separately** per band. Evaluate two independent estimators: band-by-band forward PSF-convolved galaxy models, and a PSF-homogenized *matched-aperture* measurement with traceable convolution kernel/induced covariance. Selectively report per-band integrated fluxes and corresponding alternative colours, not silently equate with total colour. Galactic extinction correction and cross-band calibration uncertainties must be recorded once and propagated.

**Acceptance test:** injection with a known $g-r$ colour gradient, variable seeing and mask holes, checking recovered matched-aperture colour on a predeclared common sky footprint. With insufficient overlapping valid coverage or wide-wing kernel, keep colour UNDERCONSTRAINED.

### F14. Injection-recovery currently validates sky sensitivity under frozen masks, not the pipeline — CONFIRMED TEST-SCOPE LIMIT (P0)

**Location:** `background.py::_audit_current_project_sky_injection_recovery` and simple Sérsic `simulation.py`.

Existing tests are valuable controlled mechanisms but hold the previous image-detected source inventory, sky training/holdout pixels, mask, aperture, PSF and reduction history fixed. There is no demonstrated full-image recovery of input galaxy total and size across the algorithm's **own** source detection, mask assignment, sky selection, isophotal geometry and aperture-radius choice. This can substantially understate selection and missing-light bias.

**Improve:** stratified and blinded full-chain native-coadd injections into untouched science data **copies**. Rerun source detection, ambiguous membership policy, PSF identification, sky tiles, centre/ellipse, aperture selection and flux integration. Preserve precise truth components: diffuse galactic light, genuine compact members, injected foreground, extrapolated outside footprint and pre/post-survey-sky differences. Separate tests with PSF-convolved sources and irregular/resolved dwarfs; simple Sérsic injections alone do not establish representativeness.

**Acceptance test:** for every predeclared morphological/size/crowding/brightness stratum, report median bias, scatter, interval coverage, fraction of missed/false total claims and accepted-sample selection bias; use genuinely held-out galaxies and external deeper data. Keep *all* rejected and unsuccessful trials in denominators, not just successful runs.

### F15. Current generalization and science release must remain gated — IMPLEMENTATION ORDER (P1)

- PSF GUI/backend is Leo P-specific despite the more general lower-level `build_empirical_coadd_psf` function; do not infer galaxy-independent experimental PSF readiness from a successful Leo P run.
- Separate the **measurement operator** from its adoption validator. P0.5/P1.1 already has useful array/data contracts; scientific statuses must be *per estimand* and include full source/mask/PSF/sky versions, joint noise, declared spatial domain, model family and acceptance evidence.
- No current GUI output should silently become an asymptotic total, stellar-population constraint or SBF distance merely because a numerical solver converges. SBF additionally needs galaxy-residual power-spectrum modeling and PSF/noise calibration specific to its Fourier estimator.

### Execution sequence by expected scientific gain and implementation risk

| Batch | Stage(s) | Do first | Evidence required before next batch |
|---|---|---|---|
| A — **implemented now** | 04.1→04.3 | Repair actual frozen-sky pivot derivative; preserve target-centred reparameterization and direct-aperture finite-difference regression | All current CI suites pass; no regression in native fixed-aperture flux |
| B — safest next diagnostic | 03→04 | Freeze source/sky scenarios, calculate **per-radius** sky Jacobian and direct exact annular response; display *conditional flux shifts* not guessed sigma | Synthetic affine-plane perturbations, masked wedges, native WCS, versioned receipts and heldout blank QA |
| C — physical degeneracy | 01→02→03 | Joint source/bright-star wing/galaxy protection ablation; no automatic source membership or sky subtraction | Leo P OBJID2035 source-specific experiment + injections with separate star/galaxy truth |
| D — sky and measurement | 02→03→04 | Full-chain post-coadd PSF-aware artificial galaxy recovery with detection/mask/geometry re-run, spatially balanced sky-holdout folds | Finite-sample bias/coverage curves across size, SB, crowding and masks; failure statuses retained |
| E — total and multiband | 04+ | Fit variable-isophote and/or structural models; reconstruct missing light; PSF-matched native multi-band colours with covariance | Independent deeper images/CCD-aware sky test and predeclared accuracy/false-adoption limits |
| F — downstream | SBF/SED | Only after an adopted flux/colour/PSF/sky domain exists, implement SBF-specific and SED-specific calibrations | SBF Fourier transfer/noise and colour-calibration eligibility or SED filter/metallicity/systematics validated |

**Stop/go condition:** A technical CI pass demonstrates software contracts, not any physical performance in Leo P. Stage 01–04 should continue returning DIAGNOSTIC_ONLY, MEASURED_PARTIAL or UNDERCONSTRAINED until source membership, large-angle PSF, residual/upstream sky bias, missing flux and full-scale covariance have *measured*, independently held-out accuracy evidence. PR #84 remains Draft and `SCIENCE HOLD` is unchanged.


## Implemented methodological changes after F01–F15 review (2026-11-13)

**Audit trail.** The original F01–F15 findings above are a frozen
methodological assessment at their stated revision. The following
is the subsequent **implemented** research-QA code, not a claim
that physical PSF wings, foreground membership or sky level have
been independently calibrated. All work remains in the existing
draft PR #84 and its existing source/test/docs directories.

| Audit issue | Implementation commit / existing operator | What changed | Scientific limit retained |
|---|---|---|---|
| F10 / F01 (radial sky response) | `ad468537`, `aperture_systematics.py` | Every frozen Stage04.2 radius now has an exact original-native-pixel Jacobian `(-sum w, -sum w*(x-xref), -sum w*(y-yref))`; adjacent-ring differences are recorded and all radial valid areas verified. Independent synthetic direct-aperture plane perturbations cover arbitrary off-centre pivots, asymmetric masks and closure. | No `Sigma_sky` posterior, correlated-coadd error, adopted subtraction or total flux. |
| F03 (sampling structure) | `e6b0f3e`, `background.py` | Added a fixed-seed 8x8 two-dimensional stratified *alternative* for the 384-native-pixel/tile sky cap; comparable to unchanged row-major V1 under a sixth predeclared sensitivity experiment. Proportional quotas honor irregular mask occupancy. | Neither policy is galaxy-free, statistically independent or a preferred science sky. |
| F06/F09 (masked light) | `df96dfd`, `aperture_systematics.py` | Every growth radius now exposes exact newly available native area and *signed* flux when only provisional source exclusion is disabled. Manual and hard exclusions stay frozen; outer result independently verified against direct integration; ring differences retained. | Recovered flux is NOT inferred galaxy membership, star subtraction or missing-light correction. |
| F02/F04 (bright Gaia veto) | `55920ad`, `background.py` | A complete count alone no longer certifies an inventory: nonmapping rows, invalid positions, invalid G magnitudes and Gaia IDs missing G force incomplete source-veto evidence. Stage03 refuses fitting rather than silently skipping potentially bright contaminating stars. | Catalog completeness/halo radii are not independently validated and do not capture uncatalogued or truly broad wings. |
| F11 (conditional covariance) | `b182faa`, `aperture_systematics.py` | Both full and smaller-size blank-ring covariance diagnostics now subtract **exactly the same frozen** constant/tilted `ApertureSkyScenario` as Stage04.2, on original unmodified pixels. Schema version 2 records the sky scenario and rejects unsuitable numerical inputs. Controlled plane+noise experiment verifies residual covariance. | Same-coadd spatial sample covariance is not calibrated science covariance, even when a known synthetic gradient is removed. |
| F02/F03 (zero baseline) | `6cc88c2`, `background.py` | Stage03 now shows a three-way numerical heldout RMS ranking with `zero_additional_correction` on the exact same samples, in addition to the preserved historical two-fit heuristic. Synthetic training/holdout discrepancy puts zero first. | Neither low RMS nor model rank authorizes subtracting real LSB outskirts or adopting a sky model. |
| F05 / F04 (PSF validation cohort) | `9b411e0`, `psf_empirical.py` | PSF0 receipts now distinguish disjoint *ePSF-fit* holdouts from population selection informed by pre-split Gaia core widths; train/holdout positions by native quadrant and Gaia magnitude ranges are explicit. | No new independent external stellar cohort or calibrated field/brightness transfer; wing/sky degeneracy unresolved. |
| Regression quality | `462ef18`, `9e0f453` | Ensured smaller-scale blank-sky scenario test uses an actually measurable noisy synthetic field rather than vacuous `all()` assertions; clarified zero-model selection flags. | Still synthetic only. |

**Immediate follow-up validation:** after verifying all scientific/GUI CI,
the first **new real-data run** should be a single read-only Leo P
project-bound QA export with both Stage03 sampling policies,
the Stage04 radial sky Jacobian and radial source-mask ablation,
the Stage04 blank covariance for the frozen conditional sky
and Stage01's PSF cohort coverage. Differences must be shown
by radius, angular sector, actual native source and mask
identity, not only one integrated flux. Preserve the full list
of underconstrained/rejected trials.

**Not implemented in this patch series:** real DR11 Leo P
end-to-end remeasurement; full chain re-detection injection/recovery;
CCD-level upstream-sky recovery; physically calibrated far-PSF wing
and brightness/position transfer; Gaia astrometric member/foreground
classification; galaxy/PSF/sky joint inference; variable native
isophotes and missing-light/total estimators; covariance calibrated
over actual galaxy angular scales; matched-PSF multiband total
colours; SBF. None may be inferred from a successful CI result.
`SCIENCE HOLD` remains mandatory.


## Follow-on method consistency review: PSF fit domain and sky-only ambiguous-source veto

This follow-on forensic pass at the current draft branch focuses on
the *same estimand under different nuisance-parameter and source-mask
operators* rather than new photometric products. Two actionable
methodological mismatches were verified in live code and changed,
while major physical calibration claims remain on hold.

### M16 — PSF held-out core test used a whole-stamp amplitude fit — FIXED

**Source:** `photometry/psf_empirical.py::_holdout_diagnostic`.
Before this change, the held-out ePSF star's amplitude was derived
by inverse-variance least squares over **all** valid pixels of the
31-pixel stamp. The published residual/normalization score was
then evaluated only over the core (radius ≤9 pixels).
An unmodelled outer stellar wing or locally subtracted background
could alter the amplitude used to judge the core, even when its
pixels were outside the reported core-quality region.

**Numerical fix:** the independent held-out core amplitude and
residual now use precisely the same valid central pixel domain.
The old whole-stamp amplitude remains *descriptive sensitivity*
evidence, explicitly kept separate from the core result and
allowed to be undefined/non-positive without invalidating an
otherwise measurable core. Synthetic fixed-PSF tests add very
large positive and negative outer-halo perturbations while
preserving the original core, showing that the core residual and
its fitted amplitude remain unchanged. No PSF kernel, star
training selection, science PSF acceptance, galaxy image or
reference FITS was changed. This operation still fits only
core *amplitude*, not the PSF's unknown far-wing flux.

**Next experiment:** blinded independent stars across positions
and Gaia magnitudes, with several **local sky annuli**, and
explicit wide-angle PSF encircled-energy curves. A good core
amplitude must not be extrapolated into a calibrated total
stellar flux or a galaxy decontamination correction.

Photutils and DAOPHOT document restricted spatial fitting
regions for stellar PSF normalization:
https://docs.astropy.org/projects/photutils/en/stable/user_guide/psf.html
and https://iraf.readthedocs.io/en/latest/tasks/noao/digiphot/daophot/daopars.html .

### M17 — protected-overlap Tractor sources could enter sky pixels — FIXED

**Source:** `photometry/scene_contaminants.py::conservative_catalogue_mask`,
`background.py::_automatic_scene_sky_pixels`.
The Scene policy rightly did *not* automatically remove a
Tractor footprint intersecting provisional dwarf protection,
because parts of it might be intrinsic H II regions, stellar
associations or diffuse target light. Yet Stage 03 previously
vetoed only actual source-masked regions and a separate
size-based target buffer. The **exterior fragment** of a
contested source footprint could remain eligible for sky
training or holdout despite being already catalogued as an
ambiguous overlap.

**Numerical/mechanistic fix:** Stage 02 now exports a separate
native Boolean `unresolved_target_overlap_sky_veto_mask`,
the count of contributing sources, and per-conflict
on-target/off-target footprint-pixel accounting, with a SHA-256
receipt bound to native pixels. Stage 03 checks the grid,
Boolean type and receipt before additionally vetoing the
*entire* unresolved source footprint solely in the selection
of sky candidates. The original source photometry mask,
protected target light, manual edits, original FITS, object
association and PSF treatment remain unchanged.

This is a one-sided conservative policy: accepting fewer sky
pixels is preferable to falsely calling a known ambiguous
source empty sky, but it can make a sky plane **unidentifiable**
in a small or crowded field. A successful sky fit still does
not prove the exterior diffuse target light is absent, or that
Tractor's truncated PSF envelopes cover all scattered light.
Tests verify source-mask invariance, source-overlap accounting,
fewer geometrically eligible sky pixels and fail-closed
detection of a tampered sky-only-mask receipt.

### M18 — unresolved risks found but not "fixed" by synthetic checks

1. **Local PSF star-cutout sky** currently estimates a robust
   plane in one inner annular region. It does not yet jointly
   fit a stellar wide-wing component, neighbour light and sky;
   varying the sky annulus radius may produce a different
   PSF. A fixed-reference star-core selection is not an
   independent stellar-population validation.
2. **Stage 02 membership and star-vs-dwarf decomposition**
   remains an observational inference. The new sky-only mask
   is not evidence that ambiguous Tractor sources are
   foreground; it explicitly preserves them in photometric
   apertures.
3. **Stage 03 original Legacy coadd sky loss** remains
   unrecoverable from an additional fitted constant/plane;
   the protected-radius buffer cannot be assumed galaxy-free.
   Restricting sky pixels does not infer the missing light.
4. **Stage 04.3 blank same-coadd covariance** is conditional on
   the sky/reference and mask assumptions; no empirical
   large-scale covariance calibration or per-CCD independent
   noise exists. Mathematical exactness of fractionally
   overlapped apertures does not license an asymptotic total.

These remaining tasks require original coadd/CCD science
inputs, star-cohort validation and end-to-end blinded
galaxy+stellar contamination injections. None of them is
satisfied by this synthetic/tested numerical revision.

**Release gate:** `SCIENCE HOLD` remains. The affected
modules produce finite-aperture/PSF/sky diagnostic evidence
only; no automatic sky or PSF model is accepted and no
galaxy total/SBF distance is created. All updates remain
on the existing draft branch PR #84 without creating
new directories or modifying `main`.
