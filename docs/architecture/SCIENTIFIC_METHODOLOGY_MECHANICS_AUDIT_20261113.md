# Desidera: end-to-end scientific methodology and mechanism reconciliation

**Audited revision:** \`75334183cdcd440933e0dd241af1322121c59246\`, draft PR #84 (\`feature/scene-wcs-editor-sc1\`).  
**Date:** 2026-11-13.  
**Assessment type:** code, scientific-contract, documented Leo P QA, and external methodological comparison. **No new original-FITS reduction or independently reproduced Leo P pixel measurement was performed in this audit.**
**Scope:** current optical DESI Legacy Surveys DR9–DR11 coadd pipeline, conditional PSF, Scene, sky, direct aperture photometry, covariance diagnostic, staged GUI and future science gates. Existing specifications are *intended requirements*, not evidence that every named algorithm has been implemented.

## Executive verdict

1. **Engineering/data-lineage validation: substantially implemented and synthetically tested.** Original quartets, IMAGE/INVVAR native units, optical MASKBITS, celestial WCS, native-source hashes, isolated project/draft identity, exact aperture overlaps, hard pixel validity and cross-stage QA gates exist. SHA-256 is a local byte-identity check, not proof of archive provenance or photometric quality.
2. **Current astronomical measurement: signed valid-area *finite-aperture* photometry, diagnostic only.** It is not a total, a decontaminated galaxy flux or an isophotal profile. Its formal variance includes only pixel-diagonal inverse variances, and it presently uses zero additional residual-sky correction. Negative pixel and aperture fluxes are retained.
3. **Main scientific failure modes are upstream/inferential, not FITS parsing:** imperfect initial Legacy sky removal, poorly identified foreground-versus-target light, finite/space-varying PSF wings, correlated/coherent sky noise, missing areas, choice of photometric geometry and extrapolated outskirts. Passing the implemented synthetic/GUI test matrix does not bound these errors.
4. **No production science acceptance:** the \`SCIENCE HOLD\` status is appropriate. Do not publish total magnitudes, $R_e$, colours, SBF distances or physical stellar-population parameters from the current Stage 04 product.

## Implemented mechanics versus scientific estimand

| Layer / implementation | Mechanic verified in code | Scientific interpretation supported | Missing prerequisite / verdict |
|---|---|---|---|
| Object project / data manager; \`surveys/{object_projects,dataset_registry,dataset_integrity}.py\` | Project-bound native source, local index, complete coadd quartet, header/WCS validation, source SHA, explicit band/release | A reproducible *candidate* optical input and sampled sky footprint | Neither all CCD exposures nor an independent photometric validation; **engineering usable** |
| P0.5 \`photometry/{native_fits,contracts}.py\` | Original IMAGE, INVVAR, NEXP and MASKBITS are read-only, native-grid, in nMgy/pixel and nMgy^-2, with hard invalid-pixel gate | Signed flux-bearing native samples with documented metadata | INVVAR lacks off-diagonal resampling covariance; no absolute/foreground-extinction or PSF correction; **input usable** |
| PSF-0 \`psf_{gui_backend,empirical,holdout_qa,core_wings_validation,...}.py\` | DAO/Gaia/Tractor-assisted candidate selection, star-cutout local sky, native ePSF construction, held-out stellar core/ring diagnostics | A **finite empirical PSF candidate** conditional on stellar sample and local sky | Specific GUI backend is **Leo P-only**; far wings, brightness/position dependence, total encircled energy and galaxy-position transfer uncalibrated; **HOLD** |
| Scene \`scene_contaminants.py\`, \`sky_masks.py\`, \`scene_source_inventory.py\` | Exact native raster masks, semantic target-protection and manual exclusion, catalogue screening, Gaia/SGA labels and fixed-source inventory | Provisional contaminant screening; target protection is **not** bad-pixel masking | Tractor classification does not establish foreground membership; bounding radii are heuristics (not wing flux predictions); H II knots/resolved members can be lost or spurious sources retained; **HOLD** |
| SA2a \`scene_source_association.py\` | Native-band detections provisionally linked by angular distance; nontransitive/multi-source ambiguity flagged | Positional candidate group, not an astrophysical identity or membership assignment | Relative astrometry/crowding/independent segmentation completeness/purity unvalidated; **HOLD** |
| Preliminary isophote \`ellipse_seed.py\` | Mask-dependent Gaussian-smoothed contour, SGA-seeded morphology, diagnostic perturbation/sensitivity | Geometry *hypothesis* and on-screen contour only | Not a fitted variable isophote, validated outer extent, asymptotic radius, measured $R_e$ or target/foreground classifier; **HOLD** |
| Stage 03 \`background.py\` | Explicitly selected, disjoint spatial training/holdout tiles; robust constant and tilted-plane residual-sky candidate fits; fixed-mask exponential injection | A *conditional additional sky level* and within-one-coadd diagnostic | Prior Legacy subtraction may already have removed galaxy light; tiles can contain unrecognized outskirts/scattered star wings, have unequal spatial sampling and correlated noise; no adopted sky or calibrated posterior; **HOLD** |
| Stage 04.1 \`aperture_flux.py\` | \`photutils\` exact fractional overlap, source/manual/hard-area ledger, signed $F_{\rm valid}$, positive-only AB magnitude, $w^2/\mathrm{INVVAR}$ variance | Observed valid-area flux in a **declared finite ICRS ellipse** | Not galaxy-only, missing-light-corrected, PSF-matched, residual-sky-corrected or total; reported $\sigma_{\rm diag}$ is not complete error; **diagnostic MEASURED_PARTIAL** |
| Stage 04.2 \`aperture_growth.py\` | Fixed set of nested reference ellipses; direct cumulative and independently integrated annular sums; closure check | Conditional valid-area curve of growth over a **predeclared grid** | Radii depend on provisional protected area + image margin; one fixed geometry, no fitted outer isophotes, plateau significance or reliable extrapolation; **diagnostic** |
| Stage 04.3 \`aperture_systematics.py\` | Analytic derivatives to constant/plane sky, source-mask counterfactual, geometry span, spatial blank-ring covariance (when available) | Conditional sensitivity and local same-coadd empirical scatter **only** | Blank controls are not independent exposures; full scale may lack enough placements, sky structures remain, inverse covariance not calibrated; **diagnostic** |
| Old \`profiles.py\`, \`sensitivity.py\`, \`simulation.py\` | Fixed pixel-scale prototype profiles and simple PSF/noise-free Sérsic mock | Useful numerical unit-test tools, **not** an adopted native-WCS isophotal/total pipeline | Avoid conflating them with the Stage 04 operator or end-to-end injection/recovery |
| Optical multi-band, SED, SBF, CCD-hybrid | Scientific methods and contracts specified in docs/notes | **Planned** quantities only | No validated matched-PSF colours, end-to-end CCD-level coadd redo, SBF likelihood/calibration, or science-adopted galaxy total in the current GUI |

### Methodology is not yet equivalent to the intended eight-stage P1.3 sequence

\`PSF_SCENE_SKY_ORDER_P13.md\` prescribes Data Preparation → provisional Scene → preliminary star-cutout sky → PSF → final PSF-aware Scene → final sky → joint refinement → photometry. The **four visible GUI tabs** implement overlapping exploratory portions (PSF-0, Scene, background, direct aperture) but are not evidence that all eight *scientific operators* execute.

In particular Stage 04 deliberately **does not apply** the Stage 01 ePSF or Stage 03 fitted residual sky; zero additional correction is hard-coded in the integrated GUI Stage 04 worker (\`gui/aperture_workspace.py\`). This is a defensible non-adoption policy while the models remain unvalidated, and a decisive limit on physical interpretation.

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

\`v\` is hard validity AND NOT provisional-source/manual-excluded, not a probabilistic galaxy membership. A valid-area mean $F_{\rm valid}/A_{\rm valid}$ is not automatically an unbiased estimate of the **whole** elliptical annulus if masked regions preferentially contain knots or foreground stars. Negative signed flux remains legitimate.

**Covariance:**

$$
\operatorname{Cov}(F_i,F_j)
=\mathbf w_i^\mathsf T C_{\rm pix}\mathbf w_j
+\mathbf d_i^\mathsf T\Sigma_{\rm sky}\mathbf d_j
+\text{PSF/membership/model/selection terms and cross-terms}.
$$

Here $\mathbf d_i$ is the signed flux derivative with respect to sky-plane parameters. The implemented formal aperture variance $\sum_p w_p^2/\mathrm{INVVAR}_p$ accounts only for the **diagonal** pixel term; a resampled Legacy coadd has nonzero pixel correlations. The current blank-ring \`np.cov\` is a useful empirical same-image diagnostic, not $\Sigma_{\rm sky}$ or a universally invertible science covariance. Nested cumulative apertures share pixels, so independent-point curve-of-growth least squares is not justified.

**Galaxy total:**

$$
F_{b,\mathrm{total}}\equiv F_{b,\mathrm{observed}}
+\Delta F_{b,\mathrm{masked}}
+\Delta F_{b,\mathrm{outside\ finite\ domain}},
$$

with an explicitly defined target-component convention. The last two terms cannot be set to zero because an observed curve is flat. Statistical and model uncertainty must be propagated; the current Stage 04 only supplies the first term for surviving pixels.

**Matched-aperture colour:** same ICRS sky domain and consistent galaxy component, with explicitly matched effective PSF or joint forward model and joint $g/r/i/z$ covariance. Merely using the same angular radius on dissimilar native coadd PSFs does not certify colour. Separate per-band asymptotic radii, common-aperture colours and total colours are different estimands.

## Independent Leo P findings already documented, not independently reproduced here

Existing researcher-supplied Leo P DR11 South r-band QA, written to \`PSF_SCENE_SKY_ORDER_P13.md\`, reported 24/7 training/holdout spatial tiles and a **frozen-mask**, post-Legacy-sky-subtraction exponential injection. The conditional recovered injected fractions are:

| Exponential scale length relative to provisional target-protection area-equivalent radius | Constant additional sky | Tilted additional sky |
|---:|---:|---:|
| 0.5 | 0.966 | 0.949 |
| 1.0 | 0.839 | 0.781 |
| 2.0 | 0.612 | 0.537 |

These numbers diagnose **additional local sky fitting absorbing artificially injected galaxy light**. They are not empirical loss fractions for Leo P, not uncertainty estimates, and do not test historical DR11 sky loss, a re-run of segmentation, a convolved PSF or foreground subtraction. This nevertheless directly falsifies any unconditional rule that a clean holdout RMS guarantees preservation of LSB outskirts.

In \`SCENE_LOGIC_AUDIT_20261110.md\`, the researcher-supplied P21/2 r-band source/sector experiment for bright Gaia star / Tractor OBJID 2035 found **zero of five** radius bins satisfying the predeclared sector-coverage quality criteria, and a low-confidence galaxy contour overlapped a large fraction of source-facing distant rings. Apparent anisotropic excess cannot therefore be attributed uniquely to the star's far PSF wing or to Leo P's outer galaxy light. The finite ePSF stamp did not reach the region of concern. P21/3 conditional segmentation comparison still requires a new original-project GUI run and cannot itself establish the physical decomposition.

## External scientific-methodology cross-check

- **Legacy Surveys DR11 data description:** coadds are inverse-variance-weighted and resampled; a survey background has already been subtracted (including spatial spline correction); Tractor operates on individual exposures using spatially dependent PSFEx models. Thus a single uniform coadd kernel and a post-coadd residual plane are approximations, not their joint CCD forward model. https://www.legacysurvey.org/dr11/description/ and https://www.legacysurvey.org/dr11/files/
- **Legacy DR11 bitmasks:** \`GALAXY\` and \`RESOLVED\` are source-context flags, not default exclusions of intrinsic dwarf-galaxy light. \`TYPE=PSF\` in masked regions is not stellar-membership proof; Gaia-retained entries can have zero fitted band inverse variance. https://www.legacysurvey.org/dr11/bitmasks/ and https://www.legacysurvey.org/dr11/description/
- **Watkins et al. (2024):** diffuse unrecognized source flux and extended stellar scattered light bias modelled sky levels. Local sky holdout and PSF-star annuli alone cannot certify the true broad-band background. https://arxiv.org/abs/2401.12297
- **Stone et al. (2021), AutoProf:** published machinery includes sky, centre, azimuthal ellipse fitting, radial profiles, direct cumulative measurements and fit checks. Desidera currently has geometry *proposals* plus fixed exact-aperture diagnostics; do not label them already equivalent to a full AutoProf-like profile. https://doi.org/10.1093/mnras/stab2709
- **Moustakas et al. (2023), SGA-2020:** custom coadds, other-source Tractor models, fixed isophote profiles and explicit CoG extrapolation demonstrate a mature comparison method. **Do not use the affected SGA-2020 \`[G,R,Z]_COG_MAG\`/\`COG_MAGERR\` columns as uncritical external truth**: the release documentation explicitly warns of a coding bug, with some other derived surface-brightness quantities remeasured from unaffected profiles. https://www.legacysurvey.org/sga/sga2020/
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

**Acceptance rule:** remain \`DIAGNOSTIC_ONLY\`/\`MEASURED_PARTIAL\`/\`MODEL_CONDITIONAL\`/\`UNDERCONSTRAINED\` until Q0–Q5 in \`PHOTOMETRY_QUALITY_VALIDATION.md\` have *measured* results under a frozen object-class-specific policy. Only an explicitly reviewed Q6 run may be \`TOTAL_VALIDATED\`. No thresholds are chosen retroactively from Leo P's own outputs.

## Audit provenance and limitations

Repository review inspected the live GitHub SHA above, including \`surveys/dataset_integrity.py\`, \`photometry/native_fits.py\`, \`photometry/contracts.py\`, \`photometry/scene_{contaminants,source_inventory,source_association}.py\`, \`photometry/psf_{empirical,joint_scene,core_wings_validation}.py\`, \`background.py\`, \`photometry/aperture_{flux,growth,systematics}.py\`, \`gui/{scene_editor,aperture_workspace,photometry_workspace,project_qa_bundle}.py\`, \`profiles.py\`, \`simulation.py\` and the existing architecture/QA documents.

The most recent passing synthetic CI suite validates engineering behavior *as encoded*, not the astrophysical accuracy of those choices. No remote original Leo P FITS were fetched or newly reduced, no real-galaxy remeasurement was repeated and no external science benchmark was executed for this audit. The present report intentionally makes **no** declaration of a science-ready PSF, total magnitude, galaxy colour, $R_e$ or SBF distance.
