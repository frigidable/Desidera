# Desidera: evidence-checked scientific photometry core roadmap
**Status:** corrected execution plan, not an adopted scientific method. **Date:** 2026-11-09. **Sequence correction:** SC1–SC4 Scene & Smart Masks precede residual galaxy sky and PH1 measured apertures.
**Source baseline:** `main` merge `30d2ac64ffba43ef86f5fc2f49ceb00c23bca9b7` (PR #83); no new branches or directories for this planning record.

## 0. Scope and source priority

This is the **next implementation plan after the Leo P PSF-0 engineering pilot**, not another PSF-0 pilot. Preserve its frozen finite 81×81 coadd ePSF evidence, all historical P15–P20 reports, original source SHA, and `science_adopted=false`/unvalidated wing flags. Additional galaxy injections belong to photometry validation **after** a defined measured-flux method, not to PSF-0.

Canonical requirement index: [Issue #44](https://github.com/frigidable/Desidera/issues/44) (117 atomic requirements, 35 WPs); authoritative relevant issue details: [#46 inputs](https://github.com/frigidable/Desidera/issues/46), [#47 calibration](https://github.com/frigidable/Desidera/issues/47), [#50 resumable analysis](https://github.com/frigidable/Desidera/issues/50), [#53 eligibility](https://github.com/frigidable/Desidera/issues/53), [#56 scene–sky](https://github.com/frigidable/Desidera/issues/56), [#58 noise](https://github.com/frigidable/Desidera/issues/58), [#59 empirical PSF](https://github.com/frigidable/Desidera/issues/59), [#62 observed photometry](https://github.com/frigidable/Desidera/issues/62), [#66 uncertainty](https://github.com/frigidable/Desidera/issues/66), [#67 injections](https://github.com/frigidable/Desidera/issues/67), [#65 benchmarks](https://github.com/frigidable/Desidera/issues/65), [#68 science release](https://github.com/frigidable/Desidera/issues/68).

Do **not** count all 35 entire WP scopes as fully blocking one another. Implement atomic bootstrap increments. The original issue #44 states 253 atomic dependency edges, but the complete machine-readable edge list was *not* independently verified from its issue body. Keep this claim provisional until #45 traceability work makes it reproducible.

Source-of-truth hierarchy for this plan: current `main` source and exact CI; November 6 [PSF–scene–sky amendment](../architecture/PSF_SCENE_SKY_ORDER_P13.md); the November 4 authoritative issue execution contracts; [photometry architecture](../architecture/PHOTOMETRY_CORE_P11.md); official survey/library documentation; cited research. Old issue text referring to #39–#43 as open drafts is HISTORICAL (they merged; #83 also merged), not present status.

**Unresolved decision that must not be silently overridden:** The November 6 PSF–scene–sky amendment explicitly requires scientifically validated observed-star empirical PSF in *each optical science band* before adopted galaxy photometry. #68 instead proposes the smallest per-quantity release (M1 observed signed fixed-aperture flux) and conditional gates. Until a recorded science decision reconciles them, build F02/F03 as **ENGINEERING/METHOD_CANDIDATE diagnostics** and block `SCIENCE_VALIDATED`/adoption without the required per-band PSF receipts. Do not weaken this by labelling observed non-parametric flux a science-certified total.

## 1. Current implementation, not merely roadmap

| Source | Implemented, reusable | Missing for accepted photometry |
| --- | --- | --- |
| `src/lsbpipe/photometry/native_fits.py`, `contracts.py` | Immutable native IMAGE/INVVAR/NEXP/MASKBITS, strict FITS/WCS, local source SHA, nMgy native units and mask reasons | Actual calibration/area edge cases, scientifically accepted real-field input gate |
| `src/lsbpipe/scene.py`, `photometry/sky_masks.py`, `smart_masks.py` | Typed grid-linked semantic masks and candidate/manual draft geometry | One reviewed final-scene selection, foreground vs target-light QA and durable decision |
| `src/lsbpipe/background.py` | Robust residual constant/tilted-plane candidates, separate sky training/holdout | Independent spatial covariance and tested galaxy–sky identifiability; no claim of recovered upstream light |
| `src/lsbpipe/profiles.py`, `geometry.py`, `sensitivity.py` | Signed valid-pixel annular flux, geometric/imaged/valid areas, non-forced COG, exact/subpixel mask options, fixed-geometry sky perturbation identities | NativeBandFrame/WCS bridge; complete uncertainty; science-grade results and QA receipts |
| `src/lsbpipe/photometry/workflow_dependencies.py` | Staged immutable dependency design, invalidation planning, PSF receipt requirements | Persisted runnable hash-addressed scheduler and selective resume |
| `src/lsbpipe/photometry/psf_empirical.py` | Finite coadd core, unchanged holdout checks, separate opt-in P15–P20 audit | E02 per-band adoption, wide wings, spatial domain, stellar/sky disentanglement still unresolved |

Do not duplicate these implementations under new names or delete historical numerical pilot modules. No science release is established by existing successful unit/CI tests alone.

## 2. Measurement estimands: freeze before algorithms

Let `v_i` be immutable valid-pixel selection, `w_{Ai}` the *fraction of native pixel i* in intended fixed aperture A, `I_i` original survey-coadd flux in nMgy/native pixel, `B_i` a separately named residual-sky candidate and `N_i` an explicitly validated neighbour model. The first measured quantity is:

\[
F_{\rm valid}(A)=\sum_i v_i\,w_{Ai}(I_i-B_i-N_i),\qquad
A_{\rm valid}=\sum_i v_i w_{Ai}\,\Omega_i .
\]

Also output `A_intended`, `A_imaged`, unmasked fractions, original signed `F_valid`, sky hypothesis/receipt and whether `N_i` was **absent**, **masked**, or **model-subtracted**. An actually missing `N_i` model is not zero physical contamination by assumption. Preserve signed annular measurements even where `F_valid<=0`, and give a flux-space uncertainty/limit; never compute a finite logarithmic magnitude for a nonpositive flux.

Keep four distinct outputs and units:
- `F_valid(A)`: sum on **actually observed valid pixels**. Not total flux, not corrected missing-mask flux.
- `I_iso(a)`: sampled/robust isophote surface brightness along a fitted ellipse. Not an integral.
- `F_infill(A)`: separately labelled *model-derived* estimate for unseen/contaminated area, with assumptions/uncertainty.
- `F_total_model`: separately labelled intrinsic or PSF-convolved extrapolated structural-model integral, with explicit model radial domain.

If `F_valid>0`, `m_AB=22.5-2.5log10(F_valid/nMgy)` is an **observed valid-area aperture magnitude** and is not the masked-area-corrected or total galaxy magnitude. Surface brightness `mu_AB=22.5-2.5log10[(F_valid/nMgy)/(A_valid/arcsec²)]` requires positive flux and a validated local WCS pixel area. For irregular masks, this describes sampled pixels, not necessarily unbiased galaxy mean intensity; report validity fractions. Original/dereddened photometry and native passbands remain separate.

The core covariance definition is \(\operatorname{Cov}(\mathbf F)=W C_{\rm pix}W^\top\) plus *joint* sky, zero-point, calibration, PSF, geometry and neighbour-nuisance propagation, with cross-terms where source/sky pixels overlap. For independent constant-sky estimate b, sky flux variance is \(A_w^2\operatorname{Var}(b)\), where \(A_w=\sum_i v_iw_{Ai}\) is the weighted **pixel** count; never replace with linear-area scaling. The coadd `INVVAR` does not supply its off-diagonal covariance. An early diagonal error field must be named `diagnostic_diagonal_only`, not `sigma_total`.

## 3. Corrected next executable GUI stage: Scene & Smart Masks (after PSF-0)

**2026-11-09 researcher correction:** The immediately preceding revision of this document incorrectly prioritized an F02 fixed-aperture adapter *before* post-PSF galaxy scene preparation. This conflicts with the already approved P1.3 sequence in [PSF_SCENE_SKY_ORDER_P13.md](../architecture/PSF_SCENE_SKY_ORDER_P13.md), #55 and #56. The following ordering supersedes that premature priority. This change is about **implementation and user-workflow order**; it does not erase the original P15–P20 measurements or change the science readiness of PSF-0.

**One Photometry Workspace, sequential real actions, no decorative empty stages:**

`Native observation & preliminary scene → PSF-0 → Scene & Smart Masks → residual galaxy sky/noise → measured aperture/isophotal photometry → optional structural/multiband inference → quantity-specific science release`.

The early `Scene_v0` and temporary sky **before PSF-0** exist to identify candidate PSF stars, mask gross contaminants and estimate *stellar-cutout* sky. They are **not** the final galactic scene or background. After PSF-0, `Scene_v1` is the next actionable GUI workbench. Although a PSF kernel is physically measured in this pilot, it remains `science_adopted=false`; therefore every PSF-dependent masking proposal must be labelled diagnostic/limited, not adopted contamination correction.

### 3.1 Scene-v1 user objective and minimal GUI contract

The astronomer should be able to:
1. Inspect one immutable native image in the exact selected band, zoom/pan/inspect WCS coordinates, and switch g/r/i/z while preserving one **sky-coordinate (ICRS) target identity** and band-local pixel geometry.
2. Seed the approximate galaxy region with a catalogue/user centre, an ellipse or a polygon. Preview a proposed isophotal/connected light region under a **named provisional residual sky and noise estimate**; distinguish bright knots/H II regions and resolved galaxy members from foreground stars. The initial ellipse is a **geometry seed**, not a measured isophotal surface-brightness profile or true outer boundary.
3. Show the protected-target footprint, alternative provisional isophotal contours, native invalid pixels (MASKBITS/IVAR/NEXP), candidate contaminant footprints and proposed sky train/holdout separately, using clear colours/legend and reason tooltips. Target protection prevents fitting out suspected galaxy outskirts but **must not restore bad detector pixels**.
4. Suggest masks for stars and neighbouring objects from WCS/Tractor/Gaia candidates, band-specific segmentation and the actual observed image, with source classification `STAR/GALAXY/UNKNOWN` and uncertainty. The existing `propose_seeded_star_mask` is a useful heuristic, **not** evidence of full star-wing coverage. A bounded PSF-core preview may help propose compact masks, but unvalidated broad wings prohibit claims of a physically optimal halo boundary or automatic scattered-light correction. Bright stars/bleeds/ghosts require conservative manual/uncertain treatment.
5. Distinguish **automatic candidate**, **accepted researcher mask draft**, **instrument-hard exclusion**, **target-protected region**, **uncertain blend**, **sky-training** and **sky-holdout**. A visual candidate must never be published as reviewed automatically.
6. Offer manual draw/edit circles, ellipses and polygons, band-specific exceptions, source identity and annotation, undo/redo, per-edit review/revert, mask overlay on/off and a meaningful side-by-side "original / proposed masks / approved scene" comparison. Shared ICRS geometry must rasterize separately on each band's actual native WCS, not copy a Boolean mask between bands.
7. Freeze an immutable `Scene_v1` receipt after explicit review, with source FITS hashes, dataset/PSF reference and limitation, WCS+band, target centre and anchor ellipse, segmentation recipe, mask/pixel counts, source IDs, each decision, reviewer, timestamp/revision, and downstream invalidation. No native FITS rewritten; other project analyses remain intact.

**GUI entry:** after PSF-0 the Photometry Workspace should expose a single functional `Scene & Masks` step (not one new screen for every internal algorithm). The scene canvas is actionable even when PSF science HOLD; PSF-dependent automatic claims remain disabled/flagged and a non-adopted scene can still be prepared. Do not require an unvalidated PSF wing model to let the astronomer draw/protect a target and conservatively mask obvious foreground sources. The scene may be saved as *prepared but science-unaccepted*, with explicit reasons.

**Minimal normal controls:** band/source, `Propose galaxy region`, `Propose contaminant masks`, `Edit masks`, `Review & save scene`. Advanced controls: sky/segmentation detection levels, smoothing scales, candidate-source policies, PSF-sensitive expansions and method/QA details. Never expose arbitrary scientific thresholds as silently calibrated S/N.

### 3.2 Precise scene data products

Separate and version:
- `TargetFootprintDraft`: tentative protected diffuse light and internal substructures, potentially multi-component; **not a hard mask excluding the galaxy from its future photometry**.
- `GalaxyGeometrySeed`: ICRS centre, angular semimajor axis, q, PA, local WCS/native ellipse and evidence/contrast regime; derived from `ellipse_seed.py` or manually supplied, subject to user review.
- `SourceMaskProposal`: source ID/category, band-specific pixel footprint and PSF-wing uncertainty; generated from `smart_masks.py` or independent catalogue/segmentation evidence.
- `SceneMaskDraft`: immutable series of ICRS shape edits and manual-review receipts, reusing `sky_masks.py`, `mask_drafts.py`, `scene.py` and `MaskEngine`; no second incompatible mask registry.
- `Scene_v1`: reviewed scene-layer composition **per band** with target protection, hard validity, contaminants, sky exclusions, preliminary sky-block candidates, QA flags and downstream dependencies. "Reviewed" refers to user selection, not scientific proof of completeness.
- `ContaminationAssessment`: unresolved star wings/blends, central knot classification uncertainty, nearby-source masked area and possible impact on future galaxy/sky flux, including `UNCONSTRAINED_WINGS` when appropriate.

The **final residual galaxy sky** is a subsequent stage: only after a protected region and contaminant scene exist should independent `SKY_TRAIN`/`SKY_HOLDOUT` blocks be frozen and physically assessed. Preliminary isophotal contours are thus useful both to protect diffuse light during scene preparation and to initialize later full radial photometry, but are not an interchangeable substitute for later free-isophote sampling/annular-flux measurement.

### 3.3 Reuse existing code rather than implementing a parallel pipeline

Existing `src/lsbpipe/photometry/ellipse_seed.py` provides connected-isophote *moment seeds* with provisional background input; `smart_masks.py` provides researcher-click anchored band-specific stellar-mask proposals; `sky_masks.py` and `mask_drafts.py` provide ICRS circle/ellipse/polygon drafts and persistence; `scene.py`/`mask_engine.py` provide immutable semantic layers and operation-specific precedence; `workflow_dependencies.py` defines the initial post-PSF stage invalidation graph. These are research primitives with targeted tests, **not an existing integrated final-scene GUI**. `photometry_workspace.py` currently shows per-band PSF0 panels and specifically states that final background and masks are not implemented.

**Known missing integration:** WCS-linked live editable canvas; project-centred target-footprint proposal and safeguarded geometry editing; source-list integration and candidate classification; mask-review/undo/redo receipts and accurate WCS preview across optical bands; band-local scene composition; robust scene QA/uncertainty handoff to galaxy sky; a real Scene stage inside Photometry Workspace. Reuse existing FITS Explorer display/WCS machinery where it meets source-safety contracts; avoid independent pixel raster/UI frameworks.

### 3.4 Scene-phase increments, in strict order

| Increment | Owning existing Issues | Actual deliverable and exit gate |
| --- | --- | --- |
| **SC1 — interactive WCS scene canvas** | A07 #55; B04 #46; A04 #53 | Single native-band image with native invalid map, project target ICRS anchor, zoom/pan, editable sky circle/ellipse/polygon, multi-band position mapping and screenshot/keyboard controls. Offline Windows GUI and WCS flip/RA-wrap tests. No auto classification yet. |
| **SC2 — galaxy footprint & isophotal-seed assistant** | D02/D04 #40/#41 integrated; A07 #55; D05 #56 | Existing isophotal connected-component/moment seed applied with provisional sky and explicit contrast recipe; multi-threshold and scale sensitivity shown; protect faint outskirts and target bright knots, flag disconnected/edge-clipped/no-signal/ambiguous cases, allow ellipse+polygon hand correction. Seed is not a final galaxy measurement. |
| **SC3 — smart star/neighbor mask proposal** | A07 #55; D07/D09/E05 #57; D05 #56 | Candidate source identities and band-specific mask footprints, Gaia/Tractor cross-check, local image segmentation and partial observed-PSF evidence. Preview/compare mask choices; saturated bright stars and unknown extended wings produce conservative/excluded or underconstrained status, no inferred subtraction. |
| **SC4 — reviewed scene and immutable handoff** | A07 #55; B06 #50; A04 #53; D05 #56 | Manual accept/reject/edit, undo/redo and atomic scene draft revisions; operation-specific precedence; project-local scene receipt/source fingerprint/mask source ID; only downstream sky/photometry invalidated. GUI "Continue to residual sky" becomes actionable with honest eligibility checks; never auto-adopts science. |
| **SK1 — residual sky after reviewed scene** | D10/C06 #58; D05/D06 #56; G02 #66 | Frozen spatially disjoint train/holdout sky on evaluated scene, robust constant/plane comparisons and null/correlation diagnostics; insufficient margin or contamination causes explicit refusal. |
| **PH1 — WCS native fixed-aperture/COG** | B04 #46; C01 #47; F02/F03 #62 | Use reviewed `Scene_v1` plus qualified `Sky_v1` to bridge `NativeBandFrame` into existing signed valid-area `profiles.py`; separate observed flux and missing area, no total-flux inference. |
| **PH2 — free isophotes / uncertainty / deeper validation** | F02 #62; G02/G03 #66; G04 #67; #65/#68 | Only after PH1, test free q/PA/centre and nonnested ellipses, correlated covariance, model bias and blind controls; later structural and multiband releases have independent acceptance. |

SC1–SC4 are the **actual next GUI development**. PH1 was mistakenly named the next immediate task in the earlier version of this plan; it is now explicitly downstream of reviewed Scene_v1 and qualified residual-sky diagnostics. Analytical test utilities for PH1 can still be developed independently without presenting user-facing galaxy flux as ready.

### 3.5 Scene scientific counterexamples and fail-closed acceptance

- **Galaxy with multiple disconnected star-forming knots:** preserve genuine target H II substructures rather than automatically masking every compact high-S/N detection; allow non-elliptical diffuse body and uncertain member stars.
- **LSB halo below an initial segmentation threshold:** target protection/radius sensitivity cannot be defined solely by one bright isophote; vary smoothing, local sky and detection threshold, show unstable outskirts and flag uncertainty.
- **Foreground star on diffuse galaxy:** compare masked-valid-area treatment, explicit star/galaxy model, and purely predicted unobserved flux as *different* estimands; never claim the star-overlap pixels were measured. Initial phase chooses mask/flag, not unsupported subtraction.
- **Bright/saturated star with unknown wide wings:** present finite-core PSF diagnostic only in its validated numeric support; unvalidated extended halo cannot produce a "clean" scene status.
- **Rotated/mirrored WCS, RA wrap and source near brick edge:** sky ellipse and manual polygon must map correctly, show lost coverage, preserve original world-coordinate edit. No transferring a raster mask from one native band into another.
- **Hard bad pixels inside protected target:** native invalidity always wins; manual `restore_auto` never restores instrumental saturation/zero-IVAR/no-coverage.
- **Sky holdout accidentally including target or neighbouring wings:** reject/flag contamination; do not automatically optimize training masks against heldout residuals.
- **Change of manual mask after PSF run:** invalidate only affected Scene/Sky/Photometry descendants, **not** original PSF star training unless that star/cutout-sky was explicitly changed. On cancel/reopen, preserve signed drafts and preceding runs without writing original FITS.
- **All-band policy:** one source ICRS identity, independent per-band observed masks and statuses; no inferred missing `i` in North. Source/PSF SCIENCE HOLD may restrict downstream adoption, not basic review and scene drafting.

Scene release is **ENGINEERING_READY / REVIEWED_DRAFT** after source-and-WCS regressions and real Leo P UI inspection. It is **not** `SCIENCE_VALIDATED` galaxy segmentation, PSF-wing removal or unbiased photometry. Any numerical recovery/completeness claims require independent field-level injection/recovery later (#67), separate from PSF-0.


## 4. Next reviewable implementation: SC1 — native WCS scene editor

**Do this next instead of the old P21b aperture measurement adapter.** Use existing `src/lsbpipe/gui/photometry_workspace.py`, current FITS Explorer raster/WCS interaction and current `src/lsbpipe/photometry/{sky_masks,mask_drafts,ellipse_seed,smart_masks,workflow_dependencies}.py` with `scene.py`; do not create another mask representation, stage-card mockup, persistent parallel project or unnecessary top-level module structure.

1. After the per-band PSF section, expose a real `Scene & Masks` tab/step in the same object-bound Photometry Workspace. The page must show the **original immutable native image**, clear band/source/version metadata and PSF science status. It must remain accessible if the PSF is diagnostic/HOLD, while PSF-dependent mask corrections remain limited/flagged.
2. Reuse current native FITS source handling and its WCS view; enforce source project UUID, band/grid identity, native invalidity, mask bits and independent band-native pixel raster. No screenshot/mock-RGB may masquerade as scientific image.
3. Enable ICRS-positioned target ellipse/circle and freehand/polygon drafts with visible drag handles, an unobscured original image, live overlay and undo/redo. Each action must round-trip via existing `SkyMaskDraft` and scene semantic layers, with explicit purpose (`protect_target`, `exclude`, `restore_auto`, `sky_train`, `sky_holdout`); do not implement actual final sky fitting or photometry in this increment.
4. Present draft-only candidate galaxy footprint initialized from `propose_isophotal_seed` and native user anchor. A connected bright isophote must be marked **preliminary**; if field light is disconnected/clumpy, edges invalid or sky contrast unstable, show uncertainty and preserve user adjustment instead of forcing a falsely precise ellipse.
5. Regression-test angular coordinates, rotated and reflected WCS, RA wrap, exact band-local masks, hard-invalid precedence, empty/no-coverage/edge region, project restart without lost edits, and no source FITS writes. Inspect real Qt screenshots at laptop resolution, use accessible labels and demonstrate actual scene interaction; synthetic GUI captures do not establish astronomical segmentation accuracy.

**First visible success:** an astronomer opens Leo P in Photometry Workspace, goes from PSF-0 to **Scene & Masks**, sees its original g-band native image, marks/protects the true galaxy, identifies contaminants and saves a versioned scene **draft**, without performing or claiming measured galaxy photometry. Next SC2/SC3 enhance the proposals, then SC4 reviews/freezes a scene and hands it to residual-sky measurement.

**Separate later work:** under #62, the existing `profiles.py` gets the native-WCS/valid-area measured-flux adapter only *after* this Scene & Masks sequence and sky candidate acceptance. Preserve the future PH1 analytic oracles (constant/exponential galaxy, irregular H II clumps, negative rings, edge truncation, non-monotonic COG, exact Photutils overlap, sky-flux identity), but they are **not** SC1 UI acceptance tests. The Photutils 3.1 versus `pyproject.toml photutils>=2.0` version mismatch remains a PH1 API selection decision, not a reason to defer SC1.


## 5. Subsequent precision and falsification rules (after SC1–SC4)

### Sky/scene and P21d
The native Legacy optical coadd was resampled (Lanczos-3) and sky-subtracted upstream. The P20 Leo P pixel-lag correlations are local *diagnostics*, not a reusable, whole-brick noise covariance or an upstream sky-loss correction. Use a fixed catalogue/geometric stratification for sky/null blocks with no outcome-driven replacement; report complete block attempts and spatial inhomogeneity. Fit constant/tilted-plane first; a more flexible spline is NOT automatically preferable. Evaluate heldout sky regions, flux shifts versus sky family, star-wing masks, contamination and positive/negative nulls. A sky fit trained on the outer galaxy must not be claimed independent of that galaxy. Unidentifiable galaxy–sky profiles return `UNDERCONSTRAINED`.

### P22 free isophotes
Start fixed geometry; release q(a), PA(a), and centre only where signal/gradient and valid sample support them. Photutils stop codes are not interchangeable with convergence. Adjacent fitted ellipses can cross: forbid naive annulus subtraction of crossing sky footprints. Isophote *sample intensity* and masked flux integral are different estimands. Outer `A4/B4` and twists are not certified tidal evidence when local gradient is ill-conditioned.

### P23 colours/forward modelling
Band filters are distinct physical passbands (DECam, BASS, MzLS, later WISE/UV). A same-angle aperture on images with different PSFs is not automatically an intrinsic colour profile; either validate an angular-resolution matching kernel (including propagated coadd covariance) or fit per-band native observations with declared PSF and forward-model assumptions. Never treat missing DR11 North `i` as measured. Preserve negative/weak linear fluxes and upper limits without artificial colour logarithms. Never use a Gaussian seeing scalar in place of the empirical science PSF when the PSF gate requires validation.

### P24 truth and external controls
Photometry injection–recovery is **separate from PSF-0**. Freeze injection families and science thresholds before examining heldout outcomes. Distinguish (a) forced known-target parameter recovery, (b) blind-detection completeness, (c) fit-success fraction, and (d) upstream-reduction sky-loss recovery: post-processed coadd injections **cannot** determine flux removed upstream. Hold out independent fields/galaxies, not overlapping brick coadds or many annuli from one object. Log **every attempted** injection, including failures, and evaluate grouped 68/95% coverage, bias, false positives and domain limitations. Leo P, NGC2146sat and UGC 2162 are development fields, not future untouched blind certification samples once used to tune decisions.

AutoProf (non-parametric), Photutils (exact apertures/isophotes) and AstroPhot (forward galaxy+sky+PSF models) are existing implementations; benchmark on **identical images, source masks, PSF, sky and estimands**. A visually deeper profile is not evidence of lower flux bias. SGA-2020 has documented unusable original ellipse FITS COG fields; SGA-2025 is a possible newer control only after its field definitions, flags and source-release provenance are verified.

## 6. Quantity states, evidence, stop rules

Separate **capability evidence** (`ENGINEERING_READY`, `METHOD_CANDIDATE`, `SCIENCE_VALIDATED`, `RESEARCH/INELIGIBLE`) from **one measurement's result state** (`VALID`, `LIMIT`, `UNDERCONSTRAINED`, `FAILED_QA`, `UNAVAILABLE`). A numerical `VALID` measurement under a local acceptance convention is NOT automatically a `SCIENCE_VALIDATED` capability. Default result/export remains scientifically non-adopted until the approved per-band PSF and independent-science requirements pass.

Suggested immutable result record (not yet implemented):
`project_uuid, dataset_id, original_image/ivar/nexp/maskbits_sha256, native_bandpass_id, measured_quantity, observation_unit, sky_version, mask_version, psf_receipt_id_or_missing_reason, geometry_wcs_id, pixel_area_model_id, estimate_signed_nMgy, area_valid/area_intended, covariance_method, uncertainty_components, source/selection/iteration_history, algorithm_version, random_seed, result_state, capability_evidence_state, reviewer_decision`.
When PSF science is unavailable, mark PSF-gated results `METHOD_CANDIDATE` and keep adoption blocked; do not erase otherwise well-measured observed flux diagnostics.

No unreviewed numeric "universal" \(\mu_{\lim}\), S/N threshold, 0.1 mag bias tolerance, or heldout 68%/95% coverage claim should be embedded by assumption. Design and freeze domain-specific test thresholds **before** blind evaluation with independent field and sample-size uncertainty. A partially valid *inner* aperture may be useful while outer profiles remain underconstrained. Do not overwrite original results on changed masks/sky/distance; use versioned selective invalidation.

## 7. Verified primary references (accessed 2026-11-09)

- DESI Legacy DR11 optical coadds: nMgy/pixel, INVVAR and Lanczos-3; see https://www.legacysurvey.org/dr11/files/ and https://www.legacysurvey.org/dr11/description/.
- Photutils aperture exact overlap and diagonal error assumption: https://docs.astropy.org/projects/photutils/en/stable/user_guide/aperture.html .
- Photutils 3.1 `AperturePhotometry` area/flag API: https://docs.astropy.org/projects/photutils/en/stable/whats_new/3.1.html .
- Photutils ellipse gradient and early termination limitations: https://docs.astropy.org/projects/photutils/en/stable/api/photutils.isophote.Ellipse.html .
- Astropy reproject surface brightness vs integrated-flux distinction: https://docs.astropy.org/projects/reproject/en/latest/explanation/algorithms.html .
- Watkins et al. (2024), *Strategies for optimal sky subtraction in the low surface brightness regime*, MNRAS 528, 4289: https://academic.oup.com/mnras/article/528/3/4289/7603727 .
- Stone et al. (2021), *AutoProf I*, MNRAS 508, 1870: https://academic.oup.com/mnras/article/508/2/1870/6373481 .
- Stone et al. (2023), *AstroPhot*, MNRAS 525, 6377: https://doi.org/10.1093/mnras/stad2477 .
- SGA-2020 published COG caveat: https://escholarship.org/uc/item/5cj0x5sm ; SGA-2025 data model: https://sga.readthedocs.io/en/latest/sga2025.html .
- Issue and code discrepancies intentionally **not hidden**: old #44/#59 draft-PR status superseded by merged main; PSF scientific gate vs minimal M1 release requires an explicit decision; the scientific library already contains signed fixed-geometry photometry but lacks native-WCS measurement acceptance; newer Photutils class and currently specified version floor differ.

**Next actual implementation recommendation (corrected by researcher):** begin with **SC1, the interactive post-PSF Scene & Masks workbench** in the current Photometry Workspace, then SC2 galaxy footprint/isophotal seed, SC3 smart source masks and SC4 human-reviewed immutable scene receipt. Only afterward implement the SK1 residual-sky and PH1 measured-aperture execution paths. Use existing #55/#56/#57/#50/#53 and current source modules; do not open duplicate roadmap Issues or add unnecessary branches/directories. Do not reinterpret a provisional isophotal seed as fitted galaxy photometry or a PSF-core pilot as validated bright-star wings.
