# Desidera: evidence-checked scientific photometry core roadmap
**Status:** execution plan, not an adopted scientific method. **Date:** 2026-11-09.
**Source baseline:** \`main\` merge \`30d2ac64ffba43ef86f5fc2f49ceb00c23bca9b7\` (PR #83); no new branches or directories for this planning record.

## 0. Scope and source priority

This is the **next implementation plan after the Leo P PSF-0 engineering pilot**, not another PSF-0 pilot. Preserve its frozen finite 81×81 coadd ePSF evidence, all historical P15–P20 reports, original source SHA, and \`science_adopted=false\`/unvalidated wing flags. Additional galaxy injections belong to photometry validation **after** a defined measured-flux method, not to PSF-0.

Canonical requirement index: [Issue #44](https://github.com/frigidable/Desidera/issues/44) (117 atomic requirements, 35 WPs); authoritative relevant issue details: [#46 inputs](https://github.com/frigidable/Desidera/issues/46), [#47 calibration](https://github.com/frigidable/Desidera/issues/47), [#50 resumable analysis](https://github.com/frigidable/Desidera/issues/50), [#53 eligibility](https://github.com/frigidable/Desidera/issues/53), [#56 scene–sky](https://github.com/frigidable/Desidera/issues/56), [#58 noise](https://github.com/frigidable/Desidera/issues/58), [#59 empirical PSF](https://github.com/frigidable/Desidera/issues/59), [#62 observed photometry](https://github.com/frigidable/Desidera/issues/62), [#66 uncertainty](https://github.com/frigidable/Desidera/issues/66), [#67 injections](https://github.com/frigidable/Desidera/issues/67), [#65 benchmarks](https://github.com/frigidable/Desidera/issues/65), [#68 science release](https://github.com/frigidable/Desidera/issues/68).

Do **not** count all 35 entire WP scopes as fully blocking one another. Implement atomic bootstrap increments. The original issue #44 states 253 atomic dependency edges, but the complete machine-readable edge list was *not* independently verified from its issue body. Keep this claim provisional until #45 traceability work makes it reproducible.

Source-of-truth hierarchy for this plan: current \`main\` source and exact CI; November 6 [PSF–scene–sky amendment](../architecture/PSF_SCENE_SKY_ORDER_P13.md); the November 4 authoritative issue execution contracts; [photometry architecture](../architecture/PHOTOMETRY_CORE_P11.md); official survey/library documentation; cited research. Old issue text referring to #39–#43 as open drafts is HISTORICAL (they merged; #83 also merged), not present status.

**Unresolved decision that must not be silently overridden:** The November 6 PSF–scene–sky amendment explicitly requires scientifically validated observed-star empirical PSF in *each optical science band* before adopted galaxy photometry. #68 instead proposes the smallest per-quantity release (M1 observed signed fixed-aperture flux) and conditional gates. Until a recorded science decision reconciles them, build F02/F03 as **ENGINEERING/METHOD_CANDIDATE diagnostics** and block \`SCIENCE_VALIDATED\`/adoption without the required per-band PSF receipts. Do not weaken this by labelling observed non-parametric flux a science-certified total.

## 1. Current implementation, not merely roadmap

| Source | Implemented, reusable | Missing for accepted photometry |
| --- | --- | --- |
| \`src/lsbpipe/photometry/native_fits.py\`, \`contracts.py\` | Immutable native IMAGE/INVVAR/NEXP/MASKBITS, strict FITS/WCS, local source SHA, nMgy native units and mask reasons | Actual calibration/area edge cases, scientifically accepted real-field input gate |
| \`src/lsbpipe/scene.py\`, \`photometry/sky_masks.py\`, \`smart_masks.py\` | Typed grid-linked semantic masks and candidate/manual draft geometry | One reviewed final-scene selection, foreground vs target-light QA and durable decision |
| \`src/lsbpipe/background.py\` | Robust residual constant/tilted-plane candidates, separate sky training/holdout | Independent spatial covariance and tested galaxy–sky identifiability; no claim of recovered upstream light |
| \`src/lsbpipe/profiles.py\`, \`geometry.py\`, \`sensitivity.py\` | Signed valid-pixel annular flux, geometric/imaged/valid areas, non-forced COG, exact/subpixel mask options, fixed-geometry sky perturbation identities | NativeBandFrame/WCS bridge; complete uncertainty; science-grade results and QA receipts |
| \`src/lsbpipe/photometry/workflow_dependencies.py\` | Staged immutable dependency design, invalidation planning, PSF receipt requirements | Persisted runnable hash-addressed scheduler and selective resume |
| \`src/lsbpipe/photometry/psf_empirical.py\` | Finite coadd core, unchanged holdout checks, separate opt-in P15–P20 audit | E02 per-band adoption, wide wings, spatial domain, stellar/sky disentanglement still unresolved |

Do not duplicate these implementations under new names or delete historical numerical pilot modules. No science release is established by existing successful unit/CI tests alone.

## 2. Measurement estimands: freeze before algorithms

Let \`v_i\` be immutable valid-pixel selection, \`w_{Ai}\` the *fraction of native pixel i* in intended fixed aperture A, \`I_i\` original survey-coadd flux in nMgy/native pixel, \`B_i\` a separately named residual-sky candidate and \`N_i\` an explicitly validated neighbour model. The first measured quantity is:

\[
F_{\rm valid}(A)=\sum_i v_i\,w_{Ai}(I_i-B_i-N_i),\qquad
A_{\rm valid}=\sum_i v_i w_{Ai}\,\Omega_i .
\]

Also output \`A_intended\`, \`A_imaged\`, unmasked fractions, original signed \`F_valid\`, sky hypothesis/receipt and whether \`N_i\` was **absent**, **masked**, or **model-subtracted**. An actually missing \`N_i\` model is not zero physical contamination by assumption. Preserve signed annular measurements even where \`F_valid<=0\`, and give a flux-space uncertainty/limit; never compute a finite logarithmic magnitude for a nonpositive flux.

Keep four distinct outputs and units:
- \`F_valid(A)\`: sum on **actually observed valid pixels**. Not total flux, not corrected missing-mask flux.
- \`I_iso(a)\`: sampled/robust isophote surface brightness along a fitted ellipse. Not an integral.
- \`F_infill(A)\`: separately labelled *model-derived* estimate for unseen/contaminated area, with assumptions/uncertainty.
- \`F_total_model\`: separately labelled intrinsic or PSF-convolved extrapolated structural-model integral, with explicit model radial domain.

If \`F_valid>0\`, \`m_AB=22.5-2.5log10(F_valid/nMgy)\` is an **observed valid-area aperture magnitude** and is not the masked-area-corrected or total galaxy magnitude. Surface brightness \`mu_AB=22.5-2.5log10[(F_valid/nMgy)/(A_valid/arcsec²)]\` requires positive flux and a validated local WCS pixel area. For irregular masks, this describes sampled pixels, not necessarily unbiased galaxy mean intensity; report validity fractions. Original/dereddened photometry and native passbands remain separate.

The core covariance definition is \(\operatorname{Cov}(\mathbf F)=W C_{\rm pix}W^\top\) plus *joint* sky, zero-point, calibration, PSF, geometry and neighbour-nuisance propagation, with cross-terms where source/sky pixels overlap. For independent constant-sky estimate b, sky flux variance is \(A_w^2\operatorname{Var}(b)\), where \(A_w=\sum_i v_iw_{Ai}\) is the weighted **pixel** count; never replace with linear-area scaling. The coadd \`INVVAR\` does not supply its off-diagonal covariance. An early diagonal error field must be named \`diagnostic_diagonal_only\`, not \`sigma_total\`.

## 3. Implementation sequence: one focused, testable increment at a time

| Increment | Atomic requirements / Issues | Deliverable and acceptance oracle | Science state |
| --- | --- | --- | --- |
| **P21a — contract reconciliation** | M01/M02/M03 #45; A04/G01 #53; #59/#68 | One quantity-by-quantity scientific-state table and PSF-eligibility decision; clarify P13 vs #68 before releasing measurements; provenance/decision mapping to tests | Specification only |
| **P21b — native measurement adapter** | B04 #46; C01 #47; F02.bootstrap #62 | Map \`NativeBandFrame\` + evaluated scene mask + one residual sky + WCS-registered ellipse to existing signed flux engine; exact fraction geometry; immutable source hash/recipe; no new global pipeline | ENGINEERING_READY after tests |
| **P21c — scene/sky candidate measurement** | D05 #56; B06.bootstrap #50; D10.bootstrap #58 | Predeclare target/contaminant and sky train/holdout blocks; compare constant vs plane on identical frozen masks/apertures; explicit insufficient/holdout failure; rerun *only* sky descendants on mask updates | METHOD_CANDIDATE |
| **P21d — empirical noise/sky covariance** | C06 #58; G02.bootstrap #66 | Blank-aperture scaling, coadd pixel-lag/PSD, full annulus/COG covariance with shared-sky nuisance; real sky holdouts frozen by geometry before reading residual values | METHOD_CANDIDATE |
| **P21e — observed fixed-aperture photometry v1** | F02/F03/F10 #62; #68 M1 | g/r/i/z *where native data actually exist*: each band independent, signed flux, intended/imaged/valid area, raw COG, angular radii, AB flux-based quantities if positive, explicit masks/sky/quality and original FITS unchanged | First defensible measurement candidate, not release |
| **P22 — identifiable isophote extension** | F02 advanced #62; D06 #56; optional A07 #55 | Reference ellipse then conditional centre/q/PA variation; diagnostic radial-gradient/fit stop code, crossing-isophote rejection; inner result survives outer fit failure; no artificial total | M2 candidate only |
| **P23 — structural model / multi-band** | F04 #63; E08 #60; F05 #64; C02 #47 | Single PSF-convolved exponential/Sérsic with explicit sky nuisance; separately calibrated colours via matched common resolution or validated native-PSF forward modelling and measured band covariance | M3/M4 separately gated |
| **P24 — independent bias/coverage validation and publication QA** | G04/G09/D08 #67; G03/G05 #66; F09/G07 #65; G06/G11 #68 | Fixed synthetic and observational challenge sets, negative blanks, matched AutoProf/Photutils/AstroPhot input tests, true scientific bias/interval-coverage/failure ledger | \`SCIENCE_VALIDATED\` only per evidenced quantity/domain |

The sequence expresses **incremental dependencies**. Begin basic analytic tests, simple source/sky eligibility and preregistration of held-out samples early; do not insist that every feature in entire #58/#66/#67/#68 issue be complete before implementing F02.bootstrap. The first usable scientific article table needs M1 and perhaps a limited M2; Sérsic/decompositions, full colour gradients, Hybrid exposure fitting, SBF, WISE/GALEX and stellar population analysis are later independent scientific capabilities.

## 4. Next reviewable implementation: P21b (start here)

Work in existing \`src/lsbpipe/photometry/\` and tests, reusing existing \`profiles.py\` and \`geometry.py\` where numerically valid. Prefer a small, typed native-frame measurement adapter rather than a parallel second photometry stack. The first PR should:

1. Load/validate one native brick/band; require trusted byte/HDU/WCS/area/units/mask evidence. Verify IMAGE nMgy/native pixel and IVAR nMgy⁻², bitfield masks, effective valid area, local pixel solid angle, extinction state and camera/passband provenance; do not infer a per-CCD zeropoint for an already calibrated coadd.
2. Define a pinned target centre (ICRS), reference major-axis angular size, axis ratio \`0<q<=1\`, position angle north-through-east; transform through local WCS to native pixels. Validate +/− RA wrap, rotated/mirrored orientation, and reference-area distortion; when a scalar-area/affine ellipse is inadequate, fail closed or apply a validated per-pixel WCS footprint.
3. Use *the same frozen fractional weights* for annular and cumulative flux/area. Prefer Photutils exact elliptical overlap where available, otherwise test the existing subpixel method against exact-area oracle. Preserve no-overlap/fully-masked/partial-edge flags and HII knots. Do not treat a segmentation-masked source as flux that has been recovered.
4. Apply one explicitly named residual-sky candidate and one explicit neighbour policy. Output signed per-ring flux and raw COG, source/recipe hashes, missing area, sky sensitivity, diagonal-only variance flag and honest quantity state. No PSF-deconvolved total or calibrated outer radius.
5. Validate analytic uniform intensity, exponential source, negative rings, edge clipping, zero area, single masked core, different pixel scales, WCS flips and RA wrap, masked bright clump, pure sky shift identity \(\Delta F=-\sum_i vw_i\Delta B_i\), and comparison against Photutils reference implementation.
6. Test original Leo P DR11 South native FITS as observational **method diagnostic**. Unit tests on synthetic FITS are distinct evidence. Record exact input file SHA, dependency version, selected band, native physical grid and full failure cases; keep original data unchanged.

**Potential dependency mismatch found in \`main\`:** \`pyproject.toml\` advertises \`photutils>=2.0\`, while the newer \`AperturePhotometry\` class with mask-aware area and bitwise flags was introduced in **Photutils 3.1**. Choose **one** explicit tested strategy: (a) require \`photutils>=3.1\` for the new class, with minimum-version/Windows CI; or (b) retain \`EllipticalAperture.to_mask(method="exact")\` on supported earlier versions and independently reproduce mask-area and quality flags. Do not import an API unavailable under the declared minimum package version.

## 5. Subsequent precision and falsification rules

### Sky/scene and P21d
The native Legacy optical coadd was resampled (Lanczos-3) and sky-subtracted upstream. The P20 Leo P pixel-lag correlations are local *diagnostics*, not a reusable, whole-brick noise covariance or an upstream sky-loss correction. Use a fixed catalogue/geometric stratification for sky/null blocks with no outcome-driven replacement; report complete block attempts and spatial inhomogeneity. Fit constant/tilted-plane first; a more flexible spline is NOT automatically preferable. Evaluate heldout sky regions, flux shifts versus sky family, star-wing masks, contamination and positive/negative nulls. A sky fit trained on the outer galaxy must not be claimed independent of that galaxy. Unidentifiable galaxy–sky profiles return \`UNDERCONSTRAINED\`.

### P22 free isophotes
Start fixed geometry; release q(a), PA(a), and centre only where signal/gradient and valid sample support them. Photutils stop codes are not interchangeable with convergence. Adjacent fitted ellipses can cross: forbid naive annulus subtraction of crossing sky footprints. Isophote *sample intensity* and masked flux integral are different estimands. Outer \`A4/B4\` and twists are not certified tidal evidence when local gradient is ill-conditioned.

### P23 colours/forward modelling
Band filters are distinct physical passbands (DECam, BASS, MzLS, later WISE/UV). A same-angle aperture on images with different PSFs is not automatically an intrinsic colour profile; either validate an angular-resolution matching kernel (including propagated coadd covariance) or fit per-band native observations with declared PSF and forward-model assumptions. Never treat missing DR11 North \`i\` as measured. Preserve negative/weak linear fluxes and upper limits without artificial colour logarithms. Never use a Gaussian seeing scalar in place of the empirical science PSF when the PSF gate requires validation.

### P24 truth and external controls
Photometry injection–recovery is **separate from PSF-0**. Freeze injection families and science thresholds before examining heldout outcomes. Distinguish (a) forced known-target parameter recovery, (b) blind-detection completeness, (c) fit-success fraction, and (d) upstream-reduction sky-loss recovery: post-processed coadd injections **cannot** determine flux removed upstream. Hold out independent fields/galaxies, not overlapping brick coadds or many annuli from one object. Log **every attempted** injection, including failures, and evaluate grouped 68/95% coverage, bias, false positives and domain limitations. Leo P, NGC2146sat and UGC 2162 are development fields, not future untouched blind certification samples once used to tune decisions.

AutoProf (non-parametric), Photutils (exact apertures/isophotes) and AstroPhot (forward galaxy+sky+PSF models) are existing implementations; benchmark on **identical images, source masks, PSF, sky and estimands**. A visually deeper profile is not evidence of lower flux bias. SGA-2020 has documented unusable original ellipse FITS COG fields; SGA-2025 is a possible newer control only after its field definitions, flags and source-release provenance are verified.

## 6. Quantity states, evidence, stop rules

Separate **capability evidence** (\`ENGINEERING_READY\`, \`METHOD_CANDIDATE\`, \`SCIENCE_VALIDATED\`, \`RESEARCH/INELIGIBLE\`) from **one measurement's result state** (\`VALID\`, \`LIMIT\`, \`UNDERCONSTRAINED\`, \`FAILED_QA\`, \`UNAVAILABLE\`). A numerical \`VALID\` measurement under a local acceptance convention is NOT automatically a \`SCIENCE_VALIDATED\` capability. Default result/export remains scientifically non-adopted until the approved per-band PSF and independent-science requirements pass.

Suggested immutable result record (not yet implemented):
\`project_uuid, dataset_id, original_image/ivar/nexp/maskbits_sha256, native_bandpass_id, measured_quantity, observation_unit, sky_version, mask_version, psf_receipt_id_or_missing_reason, geometry_wcs_id, pixel_area_model_id, estimate_signed_nMgy, area_valid/area_intended, covariance_method, uncertainty_components, source/selection/iteration_history, algorithm_version, random_seed, result_state, capability_evidence_state, reviewer_decision\`.
When PSF science is unavailable, mark PSF-gated results \`METHOD_CANDIDATE\` and keep adoption blocked; do not erase otherwise well-measured observed flux diagnostics.

No unreviewed numeric "universal" \(\mu_{\lim}\), S/N threshold, 0.1 mag bias tolerance, or heldout 68%/95% coverage claim should be embedded by assumption. Design and freeze domain-specific test thresholds **before** blind evaluation with independent field and sample-size uncertainty. A partially valid *inner* aperture may be useful while outer profiles remain underconstrained. Do not overwrite original results on changed masks/sky/distance; use versioned selective invalidation.

## 7. Verified primary references (accessed 2026-11-09)

- DESI Legacy DR11 optical coadds: nMgy/pixel, INVVAR and Lanczos-3; see https://www.legacysurvey.org/dr11/files/ and https://www.legacysurvey.org/dr11/description/.
- Photutils aperture exact overlap and diagonal error assumption: https://docs.astropy.org/projects/photutils/en/stable/user_guide/aperture.html .
- Photutils 3.1 \`AperturePhotometry\` area/flag API: https://docs.astropy.org/projects/photutils/en/stable/whats_new/3.1.html .
- Photutils ellipse gradient and early termination limitations: https://docs.astropy.org/projects/photutils/en/stable/api/photutils.isophote.Ellipse.html .
- Astropy reproject surface brightness vs integrated-flux distinction: https://docs.astropy.org/projects/reproject/en/latest/explanation/algorithms.html .
- Watkins et al. (2024), *Strategies for optimal sky subtraction in the low surface brightness regime*, MNRAS 528, 4289: https://academic.oup.com/mnras/article/528/3/4289/7603727 .
- Stone et al. (2021), *AutoProf I*, MNRAS 508, 1870: https://academic.oup.com/mnras/article/508/2/1870/6373481 .
- Stone et al. (2023), *AstroPhot*, MNRAS 525, 6377: https://doi.org/10.1093/mnras/stad2477 .
- SGA-2020 published COG caveat: https://escholarship.org/uc/item/5cj0x5sm ; SGA-2025 data model: https://sga.readthedocs.io/en/latest/sga2025.html .
- Issue and code discrepancies intentionally **not hidden**: old #44/#59 draft-PR status superseded by merged main; PSF scientific gate vs minimal M1 release requires an explicit decision; the scientific library already contains signed fixed-geometry photometry but lacks native-WCS measurement acceptance; newer Photutils class and currently specified version floor differ.

**Next actual PR recommendation:** a single, small **P21b fixed native aperture/annulus measurement adapter** with an immediately preceding PSF-eligibility/result-state contract note. No new roadmap Issues; update existing #46/#47/#53/#62 with test/commit evidence only when code and observational results exist. Avoid broad scene/sky GUI refactoring until the numerical acceptance oracles are passing.
