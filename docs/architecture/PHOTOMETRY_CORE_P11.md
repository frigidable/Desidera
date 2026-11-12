# Photometric core architecture — design baseline P1.1

> **Stage-order update (P1.3, 2026-11-03):** the logical workflow is now
> Data Preparation → Scene Initialization → Preliminary Sky →
> PSF Reconstruction → Final Scene & Smart Masks → Sky Reconstruction →
> Conditional Refinement → Galaxy Photometry & QA. See
> `PSF_SCENE_SKY_ORDER_P13.md` for authoritative ordering, PSF reuse
> and validation requirements. The milestone labels below are historical
> and must not be interpreted as executable order. P1.1 remains a
> non-measuring preparation screen.

Date: 2026-11-03. Status: scientific architecture plus a **non-measuring**
project-bound GUI preparation screen. This document is a design contract,
not evidence that the following numerical methods exist or are validated.

## 1. Navigation and authority boundaries

Galaxy Projects -> Galaxy Workspace -> Photometry Workspace.

- **Projects** owns the object identity and library entry, not a photometric
  result or a FITS copy.
- **Workspace** owns the project state, selected observational sources,
  catalogues, data manager, source QA, frozen analysis configurations, and
  history. Existing survey/FITS exploration remains available.
- **Photometry** is a dedicated task surface. Its entry point requires an
  object project. Navigation must not issue network requests, download files,
  adopt previous QA, or write results.
- External catalogues (SIMBAD/NED/Tractor) provide literature/comparison
  values, never an adopted result or hidden prior on the fit.
- Scientific results belong to a unique immutable **analysis run**. A galaxy
  can have many runs with distinct releases, sky models, masks and algorithms.

P1.1 is a preparation surface: reads the *saved* project selection and checks
whether four native optical paths exist. This is neither source integrity nor
FITS/WCS validation. It never reports `SCIENCE_READY`, displays invented
fluxes, or launches analysis. A complete native quartet is necessary for
the first intended coadd path, but not sufficient for photometry.

## 2. Scientific invariants

1. Preserve original IMAGE/INVVAR/NEXP/MASKBITS bytes and pixel coordinates.
   Use read-only references with file hash / survey identities in run receipts.
2. All photometric source measurements are signed **linear fluxes** first.
   Never log-transform non-positive fluxes; maintain flux likelihoods and
   limits distinct from unavailable measurements.
3. Make explicit all units, ICRS/WCS, calibration source, filter transmission,
   AB/Vega convention, extinction, pixel solid angle and coordinate origin.
4. The normal Legacy optical coadd IMAGE is nMgy/native pixel;
   `m_AB = 22.5 - 2.5 log10(F_nMgy)` applies to positive integrated flux
   *after* the justified background/contamination treatment. Surface
   brightness needs an explicit area correction and pixel geometry.
5. The coadd INVVAR is **diagonal evidence**, not an aperture flux
   covariance matrix. Resampling, background variation and masking can
   induce correlated errors. NEXP is a coverage/count map, not inverse noise.
6. MASKBITS is an unsigned semantic bit field, not a reject-all-nonzero mask.
   Preserve named invalid/contaminant/target-protected/user-exclusion layers.
7. Keep observed/processed/model/residual and neighbouring-source models
   separate; no silent subtraction, interpolation, PSF matching or mosaicking.
8. Overlapping bricks/releases and CCDs may reuse original exposures.
   Record shared origin and do not treat such products as independent.
9. A locally fitted residual sky cannot generally recover light removed
   upstream by the survey reduction; report identifiability limits.
10. Any automatic selection of method/parameters must be reviewable,
    versioned and scientifically justified; user overrides are recorded.

## 3. Core data flow

```text
ObjectProject + Dataset Registry
          |
Input Contract: NativeBandFrame / future CCDFrame
          |
Scene Definition: region, target center, neighbors, mask policy
          |
Sky / PSF / Noise preparation (independent versioned artifacts)
          |
   +------+------+
   |             |
Isophotal       Parametric forward fitting
non-parametric  galaxy + sky + contaminant sources
   |             |
   +------+------+
          |
Joint error model + identifiability + injection/recovery tests
          |
PhotometryRun -> ScientificMeasurements -> QA -> Export
          |
Downstream: colors, SBF eligibility, distances, populations, UDG
```

**Coadd-only (first deliverable):** independent image per
release/region/brick/band. Multi-brick analysis must wait for a validated
mosaic or sky-domain joint-fitting policy. Do not combine duplicate pixels
or average independent-looking releases.

**Hybrid (target architecture):** coadd provides detection and initialization.
Independent CCD frames then constrain a common sky-coordinate galaxy model
while retaining exposure-specific PSF, zero point, residual background,
valid-pixel map, calibration and noise. The CCD likelihood should operate
on native pixels when possible; image reprojection is an explicit derivative.
Coadd and CCD fits are not statistically independent evidence if coadd
contains those CCDs.

## 4. Reconstruction and photometric estimators

### 4.1 Observation forward model

Conceptually for calibrated exposure e and band b:

```text
y[e,b] = a[e,b] * W[e,b]((G[b] + N[b]) convolved P[e,b])
       + B[e,b] + noise[e,b]
```

G is the target's astrophysical light; N its neighbours; P the PSF
(including wings as appropriate); W the WCS/pixel integration operator;
B the residual instrumental/sky contribution. A cirrus component may
be explicitly modelled only with defensible constraints. Exact operator
order, pixel response, units and likelihood depend on the camera/product.

For Gaussian errors with known covariance `C`:
`-2 log L = r.T @ C^-1 @ r + log|C| + constant`.
The practical implementation may use empirically calibrated local
correlation kernels, blank apertures or structured operators rather
than an infeasible dense global covariance. These approximations
need recovery/coverage testing.

### 4.2 Sky & contamination

- Candidate constant/plane/low-order smooth sky models, with independent
  blank regions, target protection and masked sources.
- Explicitly fit or marginalize sky along with extended galaxy light.
  Require a held-out sky diagnostic, robust spread across plausible models
  and tests of degeneracy with effective radius and extrapolated flux.
- Extended PSF wings and halos of bright stars may require a larger field
  than the galaxy measurement window.
- A flexible spline that reproduces observed pixels is **not** proof that
  it preserved the galaxy outskirts.

### 4.3 Non-parametric photometry (baseline)

- WCS-registered elliptical geometry: center, q=b/a, PA (with documented
  on-sky convention), semimajor axis and annular boundaries in arcsec.
- First fixed-geometry, band-specific/forced-geometry aperture estimates;
  subsequently iterative free isophotes with explicit stopping criteria.
- Use exact or sufficiently converged subpixel fractional integration
  for masks and elliptical annuli, preserving native flux.
- Per-annulus *signed* flux, effective unmasked area, local surface
  brightness when positive, intensity uncertainty, flags, valid fraction.
- Curves of growth are accumulated from **measured** aperture fluxes;
  do not enforce monotonicity by modifying noisy points.
- Define a measured truncation radius from signal/background evidence;
  record extrapolated model totals separately.

Useful references:
- photutils aperture methods:
  https://docs.astropy.org/projects/photutils/en/stable/user_guide/aperture.html
- photutils isophotes:
  https://docs.astropy.org/projects/photutils/en/stable/user_guide/isophote.html

### 4.4 Parametric photometry

Initial target model: single Sérsic, convolved with an empirical PSF;
fit galaxy plus sky and selected contaminants. Store parameter covariances
and model-extrapolation sensitivity. Fit in native flux and pixels.
Alternative profiles/nuclei/disk breaks only after a model-selection and
identifiability check. A success flag from an optimizer is not enough.

Reuse and benchmark third-party solvers where scientifically suitable.
AstroPhot already supports joint multi-image modelling, PSFs, neighbours
and joint parameter covariance; Desidera should not build a competing
general optimizer without demonstrating a concrete unmet requirement.
https://academic.oup.com/mnras/article/525/4/6377/7243402
https://astrophot.readthedocs.io/en/v0.17.2/tutorials/JointModels.html

### 4.5 Multi-band analysis

Preserve individual filter-system response (BASS/MzLS vs DECam).
Fit or measure each band under its PSF and correct spatial footprint;
colors need matched or forward-model-consistent apertures.
Flux covariance across bands includes shared geometry/sky priors where
applicable. A missing or failed z/i band must not be silently imputed.

### 4.6 SBF and physical inference

Downstream only. SBF requires independently tested PSF, spatial noise,
background/masks, validated passband/population calibration and recovery
tests. Detection of SBF power does not imply a valid distance.
Padova SSP/CSP, WISE, HI and distance-based physical parameters consume
photometric products *with correlations and provenance*, never replace
the core measurement or its QA.

## 5. Typed scientific objects (planned)

- `ObservationSelection`: object UUID, source IDs, release, region, brick,
  band, exposure identities, instrument, archive path, source hashes.
- `PreparedScene`: source WCS, ROI, protected target, mask layers,
  neighbouring sources and explicit versioned footprint.
- `SkySolution`: model coefficients, units, fit/holdout regions, priors,
  posterior/covariance, residual maps, validation state.
- `PSFSolution`: empirical model/position, normalization, wings,
  residuals, source-star evidence, validation state.
- `NoiseSolution`: diagonal weights, effective correlated-noise model,
  scale-dependent blank-sky sensitivity with validation evidence.
- `PhotometryConfiguration`: frozen data, method and all settings
  with a digest. No hidden mutable defaults in published science.
- `PhotometryRun`: immutable input dependency IDs, code version, full
  environment, seeds, completion manifest, artifacts and QA.
- `PhotometricMeasurement`: method, band, aperture/model definition,
  signed flux, uncertainty and covariance link, units, effective area,
  optional magnitude, masks, provenance and scientific status.

These are **new conceptual interfaces**, not classes claimed to exist in P1.1.
Reuse P0.5 `NativeBandFrame` and P0.3 run contracts rather than building
a second registry or overwriting their schema.

## 6. Eligibility, scientific statuses and adoption

Input readiness and scientific acceptance must remain distinct:
`NO_SELECTED_DATASET`, `DATA_ROOT_UNAVAILABLE`,
`INCOMPLETE_QUARTET`, `FILES_PRESENT_UNVERIFIED` describe *file paths*,
not photometry. Subsequent validation should add explicit states for WCS,
units, masks, sky, PSF, noise and geometric coverage.

Measurement results should support `VALID`, `LIMIT`,
`UNDERCONSTRAINED`, `FAILED_QA` and `UNAVAILABLE` with a versioned
definition of their criteria. Exploratory products must retain existing
`DIAGNOSTIC_NOT_ADOPTED` / `EXPERIMENTAL_NOT_ADOPTED` labels until a
separately reviewed scientific adoption protocol is implemented.

Require flags for background–galaxy degeneracy, unmodelled wings/neighbor
contamination, incomplete footprint, inconsistent band calibration and
posterior/prior domination. No single photometric quality score may conceal
a failed blocking gate.

## 7. Screen contract and intended UX

P1.1 shows a distinct Photometry workspace only inside an object project,
with return navigation and Data Manager access. On entry it reads the saved
selection *offline*, showing the four expected input file paths and an
unambiguous non-adopted status. Any input existence check is read-only.

The eventual screen should use:
- Left: pipeline stage navigation and progress, clear blocked conditions.
- Center: large local FITS canvas (image/model/residual/masks with WCS),
  plots and cross-band controls; no remote Viewer pixels in science.
- Right: contextual inputs: region, source/mask editor, sky models, PSF,
  isophote geometry, model selection, nuisance parameters.
- Bottom: immutable run selector, stage-level error budget, QA evidence,
  audit logs, reproducible export.

A real-time preview can display exploratory derivatives but must not
write over source FITS. A new accepted run must be explicit, attributable
and reproducible, with previous runs remaining accessible.

## 8. Release sequence and acceptance gates

| Stage | Deliverable | Blocking scientific evidence |
| --- | --- | --- |
| P1.1 | Workspace navigation, read-only data snapshot, design contract | Offline GUI path, missing-file cases, no network or FITS writes |
| P1.2 | Native ROI, region editor, evidence-preserving masks | WCS/area accuracy and mask semantics synthetic + real |
| P1.3 | Residual sky candidate engine | Injected diffuse light, holdout sky, sky/galaxy degeneracy |
| P1.4 | Coadd fixed-ellipse flux/profile engine | Analytic images, signed flux, fractional masks, blank-aperture tests |
| P1.5 | Free isophotes + growth curves | Ellipse recovery and failure cases, measured vs extrapolated flux |
| P1.6 | PSF and coadd noise model | Empirical stellar PSF, wings, scale-dependent correlation |
| P1.7 | Sérsic forward fit | Model recovery, residuals, structural parameter coverage |
| P1.8 | Multi-band colors and QA | Filter compatibility, PSF effects, shared uncertainty |
| P1.9 | End-to-end runs, provenance, reporting | Immutable dependency DAG, 68/95% coverage, science gate |
| H1 | CCD exposure adapter and joint fitting | Shared-CCD lineage, native-pixel model, compare with AstroPhot |
| S1 | Experimental SBF | Masked power spectra, correlated noise, injection/TRGB tests |

Do not claim completion from synthetic-only tests; require real-FITS
acceptance and independent benchmarking before publication-grade status.
For initial *photometry* pilots use Leo P first, followed by available UGC 2162 and a diverse control set
(extended faint dwarf, compact dwarf, sky-crowded/star-contaminated field,
brick boundary, known bright comparison galaxy). Leo P can stress-test
resolved stars and irregular morphology but is not a universal Sérsic
ground truth.

**Validation questions:** bias in total magnitude, R_e, n, q, sky and
outer-profile slopes; repeatability under sky/mask changes; empirical
interval coverage; recovery across surface brightness, size, PSF and
crowding; and whether hybrid genuinely improves over existing solvers.
A claim of methodological superiority is conditional on those tests.


---

## 2026-11-11 — Stage 04 aperture-photometry research audit and implementation plan

**Document status: AUDITED DESIGN / IMPLEMENTATION NOT STARTED.**
This section is the development contract for the GUI task
`04 · Object Aperture Photometry`; it is **not** a report that any
real galaxy aperture flux or integrated magnitude has been measured.
The existing `PHOTOMETRY_QUALITY_VALIDATION.md` controls all scientific
adoption decisions and is not superseded by this plan.

### A. Primary-source methodology and scope of evidence

1. Legacy Surveys DR11 data products and reduction:
   https://www.legacysurvey.org/dr11/files/
   https://www.legacysurvey.org/dr11/description/
   Native coadd `IMAGE` pixels are **nanomaggies per pixel**,
   `INVVAR` per-pixel inverse variance in `nMgy^-2`,
   0.262 arcsec nominal scale for optical bricks.
   The stacks are already sky-subtracted and Lanczos-3-resampled,
   not the single-epoch images used in Tractor precision fitting.
   No residual correction can recover signal removed upstream.
2. SGA-2020 methods (Moustakas et al. 2023,
   DOI 10.3847/1538-4365/acfaa2), with public DR9 precedent:
   https://a.legacysurvey.org/sga/sga2020/
   https://www.desi.lbl.gov/2021/06/29/siena-galaxy-atlas-2020-publicly-released-by-desi-legacy-imaging-surveys-team/
   Source models/masks, fixed-ellipse geometry, surface-brightness
   profiles, *separately measured* curves of growth and
   asymptotic fits are useful comparator strategies; these
   products do not supply independent photometric truth for
   DR11 and should not override source membership in Leo P.
3. SGA-2025 documentation (v1.0/draft publication caveat):
   https://sga.readthedocs.io/en/latest/sga2025.html
   Describes integrated `FLUX_AP`, `FLUX_ERR_AP` and
   `FMASKED_AP` for [0.5,1.0,1.25,1.5,2.0]×moment SMA,
   `ELLIPSEPHOT`, curve-of-growth parameters, fit/mask flags
   and a crucial distinction between fixed and fitted geometry.
   The documentation currently displays a **draft/not-yet-public**
   warning despite its release-history note; do not assume a
   particular catalog/geometry is immutable or validated.
4. AutoProf (Stone et al. 2021, MNRAS 508, 1870;
   DOI 10.1093/mnras/stab2709):
   https://academic.oup.com/mnras/article/508/2/1870/6373481
   Distinguishes direct unmasked-pixel flux summation from
   radial-profile integration, records masked pixels and
   permits variable isophotes; a fitted ellipse or model is
   not equivalent to the direct observed finite-aperture sum.
5. Photutils aperture API (maintained documentation):
   https://docs.astropy.org/projects/photutils/en/stable/user_guide/aperture.html
   https://docs.astropy.org/projects/photutils/en/stable/api/photutils.aperture.AperturePhotometry.html
   `AperturePhotometry` is preferred over its legacy function.
   `method="exact"` provides fractional pixel/ellipse overlaps.
   `flux_err` is a *diagonal-pixel* propagation and does not
   include unknown coadd off-diagonal covariance or sky/geometry.
   `area` is the effective unmasked overlap, whereas the
   aperture's ideal area is a different quantity.
6. Photutils ellipse/isophote algorithm (Jedrzejewski 1987):
   https://docs.astropy.org/projects/photutils/en/stable/api/photutils.isophote.Ellipse.html
   Iterative free-isophote fits can fail for masked/low-S/N
   radii. Fixing center/PA/ellipticity or switching to
   noniterative sampling is an explicitly reported method
   change, **not** a silent recovery from an invalid solution.
7. Correlated-noise aperture controls, *lambdar*:
   https://doi.org/10.1093/mnras/stw832
   and DustPedia/CAAPR:
   https://doi.org/10.1051/0004-6361/201731419
   Empirical blank/source-masked apertures of matched size
   and shape are useful checks on resampling covariance and
   structured sky, but may be unavailable at Leo P's large
   scales. Never substitute tiny apertures or treat overlapping
   same-coadd draws as independent.
8. Schombert, "The Structure of Galaxies I: Surface Photometry
   Techniques":
   https://connectsci.au/as/article/29/2/174/117537/The-Structure-of-Galaxies-I-Surface-Photometry
   Sky-level uncertainty dominates large, faint apertures;
   an apparent CoG plateau or unregularized asymptotic fit
   is not evidence of all true outer-galaxy light being present.

### B. Repository architecture audit — what exists, what does not

**Existing and reusable, do not duplicate:**

- Native strict array/WCS/NEXP/MASKBITS contracts in
  `src/lsbpipe/photometry/contracts.py` and
  `src/lsbpipe/photometry/native_fits.py`;
  `NativeBandFrame` keeps immutable source pixels and
  contextual MASKBITS separated from hard invalidity.
- `src/lsbpipe/scene.py` has typed `LayerKind`,
  `SceneState`, `Operation.PHOTOMETRY`, and mask
  evaluation. Crucially, `TARGET_PROTECT` is **excluded
  for sky fitting but is not automatically rejected in
  photometry**. Preserve this.
- `scene_editor.py`, `scene_contaminants.py` and
  `sky_masks.py` have the current project-native dataset,
  revisioned ICRS draft, an original Tractor/Gaia inventory,
  `effective_auto` after explicit restoration and
  `manual_exclude`. Mask and membership evidence remain
  **provisional**, not a final source-subtracted image.
- `ellipse_seed.py` contains initial/local diagnostic
  geometry, sensitivity trials and an SGA-informed seed;
  **not** a scientifically fitted free-isophote series.
  The Leo P prior failed a prior centeredness test: a
  provisional `trial 0` had ≈66.6 native-pixel centroid
  drift. It cannot define the science centre or outer radius.
- `background.py` and `gui/sky_workspace.py` deliver
  project-bound non-adopted zero/constant/plane sky
  scenarios, fixed sky-sample provenance, grid sensitivity
  and frozen-mask artificial-source tests.
- `workflow_dependencies.py` already declares PSF,
  final scene, sky and refinement as *science*
  dependencies and has fail-closed per-band star-derived
  PSF gate; `PHOTOMETRY_QUALITY_VALIDATION.md` defines
  `MEASURED_PARTIAL`, `MODEL_CONDITIONAL`,
  `UNDERCONSTRAINED` and `TOTAL_VALIDATED`.
- `PhotometryWorkspace` has three live tabs and
  project navigation/resume; do not build a second
  photometry application or scatter independent
  project-directory owners across tabs.

**Missing/at-risk before Stage 04:**

1. No current real native-aperture numerical engine,
   aperture quality object, radial coverage ledger,
   covariance estimate, model-conditional missing-light
   correction or adopted `PhotometryRun`.
2. Current Sky v2 output has no scientifically
   accepted local background/covariance solution;
   `zero_additional_correction_DR11_coadd` is the
   **default analysis of the already sky-subtracted
   pixels**, not proof upstream subtraction is unbiased.
3. Current Scene's `protected` region and
   `isophote_candidate` are *sky exclusion candidates*,
   **NOT the galaxy's exact geometry or missing-pixel
   mask**. Source `auto_mask` is not a definitive
   classification of foreground vs genuine H II knots,
   resolved target stars or associations.
4. In Scene v0.3, `Operation.PHOTOMETRY` can exclude
   hard defects and contaminant layers; a Stage 04
   adapter must construct these layers from the
   **effective** Scene draft (`effective_auto`,
   `manual_exclude`) with source membership QA,
   rather than blindly using a conservative background
   mask, original `auto_mask` or treating the sky
   target-protection disk as invalid data.
5. `PhotometryWorkspace` navigation currently only
   knows `psf0/scene/sky`; a fourth stage needs
   schema-valid persisted resume, guarded source-band
   handoff, stale input invalidation and a background
   worker with no FITS writes on merely opening.
6. `photutils` lives in the optional `exact` extra
   (`photutils>=2.0`); GUI-only installation does
   not currently guarantee it. Stage 04 must
   either require an explicit version-compatible
   optional dependency with a clear disabled state
   or use a separately validated exact overlap engine.
   No silent substitution to center-only summation.
7. DR11 target masks/Tractor/SGA flags do not
   independently certify source membership; the
   close Gaia foreground OBJID 2035 is a required
   overlap/wing stress test, not a source for
   automatically interpolated flux into Leo P.

### C. Stage 04 science contract and minimal estimator

**Primary Stage 04 output:** a signed, measured **finite
aperture observed-valid-area flux**, by filter, radius
and fully specified geometry, with separate masks
and errors. It is NOT a total magnitude.

Let `G_p(A)` be the exact geometric fraction (0..1)
of native pixel `p` covered by a specified on-sky
aperture. Let `V_p` be 1 only for original data-valid
and currently retained pixels under the declared
source/membership policy. `I_p` is the original DR11
coadd in nMgy/pixel and `B_p^{(s)}` is one of
`zero_additional`, `constant`, or `tilted_plane`
residual-sky scenarios **with declared scientific state**.

\[
F_{\rm valid}^{(s)}(A)=
\sum_p G_p(A)V_p\,[I_p-B_p^{(s)}]\quad{\rm nMgy}.
\]

Default `s=zero_additional` sets `B_p=0` and
leaves original FITS unchanged. Constant and
plane give *conditional sensitivity alternatives*,
not scientific adoption. The sky model at each
pixel must be evaluated with its original
reference `(x_ref,y_ref)` in the same native grid;
never use its intercept as a constant substitute
for a tilted field.

**Four separate area denominators** must be tracked:

\[
A_{\rm ideal}(a,q)=\pi q a^2,\qquad
A_{\rm inframe}=\sum_{p\in\mathrm{FITS}}G_p,\qquad
A_{\rm valid}=\sum G_pV_p,\qquad
f_{\rm valid}=A_{\rm valid}/A_{\rm ideal}.
\]

The first three are in native pixel-area units
when `a` is in native pixels; for angular areas
multiply by an accurate local per-pixel area
or integrate the WCS solid-angle map as needed.
Report inframe fraction, hard-invalid fraction,
provisional source/membership mask fraction and
`f_valid` **separately**. Mask reason accounting
must not double-count overlapping reasons:
store raw overlapping reason footprints *and*
an explicitly prioritized mutually exclusive
area ledger summing to `A_inframe`.
A high total masked fraction is not automatically
repaired or normalized away.

**Geometric model.** Freeze a single project-bound
ICRS center and candidate semimajor axis `a`
in arcsec, `0<q=b/a<=1`, and PA measured
north-through-east (east of north).
At the current target position compute a
local celestial-to-native-pixel Jacobian:
RA periodicity, cos(dec), handedness, off-axis
distortion and pixel-coordinate convention
(`origin=0`) must be handled explicitly.
`photutils.EllipticalAperture.theta` is an angle
in the pixel coordinate frame and must NOT be
set equal to the astronomer's PA without
the local WCS rotation/reflection transform.
A fixed ICRS aperture should be reprojected to
each filter's own native WCS; never reuse another
filter's mask array by shape alone.

**Geometry selection for the first run:**
automatic, target-anchor-centered, fixed ellipses,
with a preliminary axis-ratio/PA candidate from
SGA or validated local moment seed **marked as
a hypothesis**; also supply an anchor-centered
circular comparison and predeclared alternative
q/PA values when the seed is low-confidence.
No manual sky drawing is required and the
user is not required to manually draw apertures.
An optional geometry override is a versioned
new hypothesis, not an in-place edit of an
adopted run. Never recenter the galaxy by
chasing the brightest resolved star or
the 66.6px-offset Leo P trial.

**Integral versus surface brightness.** The
primary sum uses actual native flux per pixel.
For an annulus `a_in<a_out` with frozen
geometry, independently integrate its exact
fractional overlap `W_p=G_p(a_out)-G_p(a_in)`
(with numerical tolerance for 0..1). The
observed-valid mean intensity per arcsec² is

\[
\langle I\rangle_{\rm ann} =
\frac{\sum_p W_pV_p(I_p-B_p)}
     {\sum_p W_pV_p\Omega_p},
\]

where \(\Omega_p\) is the physical area in
arcsec²/pixel. Reporting its AB surface
brightness \(\mu=22.5-2.5\log_{10}
\langle I\rangle\) is allowed only for positive
intensity; always preserve the signed native
intensity and its uncertainty/limit. For integrated
positive `F_nMgy`, `m_AB=22.5-2.5log10(F)`;
do not claim a logarithmic magnitude for a
negative/zero flux or use symmetric magnitude
errors at extremely low S/N.

**Photutils sum/statistics asymmetry:** the exact overlap
weighting of `AperturePhotometry` and the sum properties
of `ApertureStats` must not be confused with the
latter class's center-based median/standard-deviation
statistics. Annular mean intensity is explicitly
the fractional-area-weighted native sum over
fractional-area-weighted valid area; any clipped
median or robust sector statistic is a separate
labelled estimator with its own geometry and
error model. Never substitute a median value
into an exact integrated flux without an
explicit reconstruction model.

**Direct curve of growth (CoG):** a nested
sequence of **observed** `F_valid(a)` and
`A_valid(a)`, with the *same fixed ellipse*
and band-local Scene/sky scenarios. Do not
force cumulative flux monotonicity or infer
an asymptote from an apparent plateau.
Annular values and direct aperture values
should agree under frozen geometry to
documented numerical precision. If
`q(a)`, PA(a) or center(a) vary later,
annulus widths/overlaps and missing coverage
must be redefined explicitly: subtracting
fluxes from geometrically **nonnested** variable
ellipses does not produce an annulus.

**Mask loss and source membership:** no default
azimuthal scaling by `A_ideal/A_valid`, no
interpolation of masked data treated as measured,
no automatic deletion of all Tractor or Gaia
sources within a proposed galaxy contour.
Save `F_valid`, signed `F_reconstructed` as
separate **model-conditioned** estimates only
in a later phase with tested recovery. Keep
resolved target members, star-forming knots,
ambiguous objects and foreground stars as
independent membership hypotheses where relevant.
A full "galaxy light" product must specify
whether target compact associations and genuine
emission line regions are included in its definition.

### D. Quantitative uncertainty, covariance and scientific limits

The coadd `INVVAR` allows a *labelled formal*
diagonal-only estimate

\[
\sigma^2_{\mathrm{diag}}(A)=
\sum_p\frac{(G_pV_p)^2}{\mathrm{INVVAR}_p},
\]

which is **not** the calibrated total variance
after Lanczos-3 resampling, nor necessarily a
mathematical bound on it. The general flux
covariance is \(\boldsymbol{w}^{\mathsf T}
C_{\rm image}\boldsymbol{w}\);
cross-radius covariance is nonzero for nested
apertures. Do not derive \(\sigma_{\rm sky}\)
as pixel RMS divided by \(\sqrt{N_{\rm pix}}\).

The effect of a coherent additive sky offset
\(\delta b\) in nMgy/native pixel is

\[
\delta F(A)=-A_{\rm valid,pix}\,\delta b,
\]

where `A_valid,pix=sum(GV)` is the
weighted usable area. For a tilted plane
and uncertain coefficients \(\boldsymbol{\beta}\),
use the aperture's design-integral vector
\(\boldsymbol{H}_A=\sum_pG_pV_p[1,
x_p-x_{\rm ref},y_p-y_{\rm ref}]\).
Then the sky contribution to two aperture
covariances is
\(\boldsymbol{H}_A C_\beta
\boldsymbol{H}_B^\mathsf T\), **only if**
the sky coefficient covariance has actually
been calibrated. Cross-correlations with
galaxy/neighbor components must also be
accounted for when estimating a total budget.

**Required diagnostic noise programme:**
(1) retain the diagonal formal term with its
explicit limitations; (2) measure flux sums
for source-vetoed, separated **matched-shape
blank apertures** on the same project/brick
at multiple radii; (3) report their number,
mask fraction, positions, overlap, depth,
spatial heterogeneity and empirical flux
scatter, along with inability to estimate
a scale when there are too few nonoverlapping
controls; (4) use radial *joint* blank
vectors or validated covariance simulations
to estimate cross-radius covariance; (5)
propagate predeclared zero/constant/plane
sky and mask/geometry alternatives separately.
Do not present `1/sqrt(Nnative_pixels)`
or a handful of overlapping controls as
calibrated precision.

**Separate numerical uncertainty from major
systematics:** formal coadd-diagonal, resampling/
structured sky, sky/galaxy degeneracy, target
geometry, missing footprint and neighbour/source
ownership, PSF core and wings, zero point,
Galactic extinction, and outer extrapolation.
For a masked or very low-S/N aperture, retain
linear flux intervals and limits; magnitude
intervals are optional derived products.
Never quadrature-sum scenario ranges as if
they were independent Gaussian errors.

**Leo P observed sky stress evidence** (researcher
2026-11-11 Stage 03 readout): no fitted additional
sky is adopted; constant vs zero-coadd holdout
RMS improvement 0.111%, tilted plane vs zero
0.264%; sensitivity of the constant target
prediction over the original grid plus
four measurable alternatives spans
`+0.00002335535..+0.00041074240`
nMgy/pixel, which is **not** a confidence interval.
A former frozen-mask artificial aperture of
radius ≈120 arcsec had only 78.3% usable area
and fitted-sky recoveries ranging from 0.966
to 0.612 (constant), and 0.949 to 0.537 (plane)
as synthetic profile scale grew from 30 to
120 arcsec. These are **not measured Leo P
light-loss fractions**. The extra diagnostic
leakage within original sky training pixels
grew from ≈2.76e-6 to ≈7.67e-4
nMgy/native pixel. Stage 04 must therefore
never substitute those synthetic ratios
as aperture corrections, and must report
actual radius-dependent valid-area fractions.

### E. Mandatory scientific status separation

Stage 04 is permitted to compute and show
`DIAGNOSTIC_ONLY` signed finite-aperture
fluxes and `MEASURED_PARTIAL` valid-area
fluxes under full source lineage; this **does
not mean** all-band science photometry passed.
`MODEL_CONDITIONAL` is a separate later
product for reconstructed missing pixels
or finite-radius profile integrals. Neither
diagnostic pixels nor a SGA estimate can
become `TOTAL_VALIDATED` without the
existing Q0–Q6 gates in
`PHOTOMETRY_QUALITY_VALIDATION.md`.

For Leo P the per-band empirical PSF remains
experimental, and g-band remains insufficient
under the predeclared stellar-holdout policy.
Existing `require_empirical_psf_coverage`
**must be invoked at the trusted science
adoption/execution boundary**, never passed
with a manually filled `True`. This does
not prevent a clearly labelled, standalone
single-band *diagnostic sum* in native pixels:
such a sum does not physically require PSF
deconvolution, but its aperture/outer-light
interpretation and final science corrections
do require band-local PSF QA and the other gates.
No implicit g->r PSF borrowing, made-up
sky posterior, zero-point correction or
`SCIENCE PASS` from finished UI.

The GUI's Stage 03 background result is
currently an in-memory *diagnostic* result,
not an immutable adopted `SkySolution`.
Stage 04 must freeze a precise source-identity
snapshot and reject stale project/draft/sky
objects; default to the original DR11
coadd (zero additional sky) even if Stage
03's numerical preference is `constant`.
Model-subtracted full-frame FITS should not
be generated to feed Stage 04.

### F. Concrete implementation phases (no new directories/branches)

**04.0 — Contracts and science receipts (first code milestone).**
Reuse existing `contracts.py`, Scene and
`PHOTOMETRY_QUALITY_VALIDATION.md`. Implement
small frozen schemas for
`ApertureGeometryICRS`, `ApertureScenario`,
`ApertureFluxRecord`, `AnnulusFluxRecord`
and `ApertureQA`, with API version,
input SHA/WCS and Scene/source/sky hashes,
band, mask policy, signed nMgy and all
areas. For missing advanced QA, explicitly
store `None/NOT_CALIBRATED`, never zero
as if calibrated. Reviewable statuses:
`DIAGNOSTIC_ONLY`, `MEASURED_PARTIAL`,
`UNDERCONSTRAINED`, `MODEL_CONDITIONAL`;
`adopted_for_science=false`.
Acceptance: JSON round-trip, invalid
units/WCS/mask/sky/dataset fail closed,
no source array writes, deterministic hash.

**04.1 — WCS and exact-pixel aperture integrator
(single-band `r` first).** Use the
existing optional photutils `exact`
extra, with explicit dependency/ABI QA;
implement original-pixel ICRS ellipse to
native transform and exact fractional
elliptical, circular and annular overlaps.
Construct hard/neighbor/manual masks from
the current *effective* Scene policy.
Return signed `F_valid`, ideal/inframe/
covered/valid/masked ledger, per-sector
coverage, `sigma_diag` as uncalibrated,
and optional positive-only AB magnitude.
Acceptance: analytic constant/plane images,
subpixel shifts, rotated/negative-CD WCS,
RA wrap, exact masked-area integration,
NaN/zero INVVAR and edge fail-closed.
No interpolation and no full aperture
rescaling.

**04.2 — Automatic geometry proposals and
multi-radius direct curves of growth.**
Use project anchor as frozen center,
fixed candidate q/PA and a neutral circle.
Use SGA/local moments as *proposals*
with provenance, not physical truth.
Record predeclared radial spacing and
outer valid-coverage constraints.
Produce independently computed direct
elliptical-annulus surface brightness,
signed aperture CoG and area/sector
plots. Initial geometry defaults can be
reviewed/overridden but never require
a hand-drawn aperture.
Acceptance: nested-flux/annulus identity
and no forced monotonicity, stable positive
F for analytic exponential test images,
no hidden extrapolation. For low-confidence
isophote or star overlap, tag
`GEOMETRY_UNDERCONSTRAINED` / 
`SOURCE_MEMBERSHIP_UNRESOLVED`.

**04.3 — Sky/geometry/source scenarios and
spatial-scale noise.** Always retain
zero additional sky. Show constant/plane
and predeclared sky-grid envelope as
*conditional alternative flux curves*;
never present envelope min/max as 68%
confidence intervals. Integrate
aperture-matched blank-source controls
at multiple scales with minimum
geometric-spread requirements;
report unmeasurable covariance rather
than `sigma_diag` as full `sigma`.
Introduce source-policy sensitivity,
especially the Leo P Gaia near-centre
star and uncertain resolved compact
members. Run controlled injection/restore
without changing source FITS. Acceptance:
known correlated-noise simulated fields,
blank-aperture scale dependence,
flux response \(\partial F/\partial b
=-A_{\rm valid,pix}\), scenario leakage,
and fail-closed controls when all blanks
overlap the galaxy or have high mask loss.

**04.4 — Existing GUI integration.**
Append `4 · Object Aperture Photometry`
inside current `PhotometryWorkspace`
with the same project and band owner,
existing scene editor and Stage 03
diagnostic snapshot (no competing project
state). Project-bound worker and
no-run-on-entry rules remain. Main
view: original native image with
clearly distinguished candidate
ellipse(s), hard/neighbor masks and
valid coverage. Adjacent views:
signed curve of growth, intensity/
\(\mu\) vs semimajor axis, coverage
vs radius/sector, sky-scenario
comparison and uncertainty flags.
Controls: `Build diagnostic apertures`,
optional reviewed geometry choice,
`Copy complete aperture QA`,
existing-project saved run/export with
explicit `UNADOPTED` names.
No compulsory user region drawing.
Changing scene draft/band/sky/source
invalidates stale photometry and
restoring Stage 03 never clobbers
Stage 02's band/masks.
Acceptance: Qt offline screenshot,
responsive narrow window, project
switch, worker cancellation/staleness,
no FITS modification or network access.

**04.5 — Band-local `griz` and comparable
aperture colours (later, separately gated).**
Project the same *ICRS* ellipses independently
into the four native WCS grids, audit
coverage per band and PSF matching or
joint forward models before an adopted
colour. Maintain genuine band-local
outer structure: do not cap g-band
light at r or z isophotal radius.
Source-component ownership consistent
across bands but may be uncertain.
Acceptance: WCS handedness, different
PSFs/depths, shared geometry covariance,
signed low-S/N flux-vector handling,
and refusal to adopt a mixed-PSF
matched colour without checks.
The g empirical PSF remains a blocking
science issue, not a reason to fabricate
a g-band kernel.

**04.6 — Separate corrected/aperture-total
research (not the first Stage 04 MVP).**
Only after the direct measurements and
their independent validation: optional
profile integration over masked areas
with explicit model posterior, free
isophotes `photutils.isophote`, partial
pixel reconstruction, model extrapolation
beyond a finite radius and asymptotic
CoG estimates. Each changes the estimand
to `MODEL_CONDITIONAL`; a zero derivative
of a noisy CoG is never proof that
remaining light is negligible. Source
membership, sky and PSF uncertainties
must propagate jointly and include
model-family alternatives.

**04.7 — Science-adoption gate and external
benchmarking.** Execute Q0–Q6 against a
policy frozen prior to target inspection;
cross-check diagnostic values against
the published **methodology** of
SGA/AutoProf and against synthetic
truths independently generated from
the exact Stage 04 instrument/PSF
models. Same DR11 galaxy products
are correlated comparisons, not
independent evidence of zero bias.
Test other already available native
galaxy projects only within authorized
data scope; Leo P alone is insufficient
for a general published accuracy claim.
Publishing `TOTAL_VALIDATED` requires
independent proof of missing-light
bounds, adequate sky/PSF/noise
calibration and explicit adoption.

### G. Regression/acceptance matrix for the first implementation

| ID | Required experiment | Mandatory observable / gate |
|---|---|---|
| A04-01 | Constant-image exact ellipse, fractional edges | `F=c*A_valid` in original nMgy; absolute area tested |
| A04-02 | Original zero sky vs constant and tilted alternative | Same pixel weights; exact `deltaF=-sum(w*deltaB)` |
| A04-03 | Annuli and nested CoG on analytic Gaussian/exponential | `F_outer-F_inner=F_annulus` only with fixed nested geometry |
| A04-04 | Centre shifted by fractions of a pixel | No step-wise binary aperture bias; compare photutils exact |
| A04-05 | WCS flip, rotation, RA=359.99/0.01 wrap, varying scale | Same ICRS ellipse, correct orientation, sky area |
| A04-06 | Foreground hole, large wedge and off-frame segment | Geometric/inframe/valid/excluded reasons distinguished; no invented flux |
| A04-07 | Negative annulus, near-zero cumulative flux | Preserve signed nMgy; no logarithmic magnitude |
| A04-08 | Zero-IVAR, nonfinite image, saturation/NEXP defect | Immutable hard mask wins over any manual restore |
| A04-09 | Tractor masked H II knot and central Gaia star | Target membership alternative, no unconditional source deletion |
| A04-10 | Missing candidate isophote or ≈66px shifted Leo P trial | Anchor circle/explicit low-confidence geometry, not false fitted ellipse |
| A04-11 | Nested covariance and correlated synthetic noise field | Diagonal only flagged; no `1/sqrt(Npix)` scientific claim |
| A04-12 | Masked/highly crowded blank controls | Empirical large-aperture noise refuses precision when sample unavailable |
| A04-13 | Injection of circular/elongated/irregular faint galaxy | Measure partial flux bias separately from model correction and PSF wing bias |
| A04-14 | Noisy false CoG plateau and negative outer rings | No forced monotonicity or automatic total flux |
| A04-15 | Readout/QA serialization, project/sky/PSF/mask mutation | Full immutable lineage, stale fail-closed, no source FITS write |
| A04-16 | Four different band PSFs and native WCSs | No cross-band mask reuse or unvalidated science colour |
| A04-17 | Narrow-window GUI and interrupted worker | Visible warning, safe cancellation, zero network/file writes on entry |
| A04-18 | Exact `A04-01..17` integration on Leo P original data | Only native-band diagnostic receipts; science HOLD until Q0–Q6 |

### H. Development sequencing and explicit exclusions

Implement **04.0 + 04.1 first** as a small read-only,
tested science-core module in the existing
`src/lsbpipe/photometry/` package and existing
`tests/` directory. Then 04.2, 04.3 and
GUI 04.4, followed by 04.5/04.6/04.7
only when their prerequisites are actually
available. The `exact` optional dependency
needs a clear installation contract; do not
silently downgrade fractional integration.

For the first Leo P pilot, report only
band-r *original-coadd*, signed, partially
observed direct aperture fluxes over
several predeclared anchored geometry
hypotheses. No new brick, forced manual
sky selection, new branch, new project
directory, aperture-size optimization to
force a plateau, fake sky noise
calibration, star-flux reconstruction,
automatic sky subtraction or
unqualified total magnitude is permitted.

**Pending researcher review:** this is a
finished scientific design/audit document,
not implementation approval for a
science-quality total. The first code
milestone should start only from the
specific typed aperture contract above;
all status flags and source provenance
must be test-covered before the UI
may display integrated magnitudes.


---

## 2026-11-11 — First Stage 04 science-core implementation

**Scope: implemented numerical engine only, GUI and real Leo P run
NOT yet claimed.** The preceding specification and Q0–Q6 gate
remain authoritative.

Added `src/lsbpipe/photometry/aperture_flux.py`, with a public
import through the existing `lsbpipe.photometry` package:
`ApertureGeometryICRS`, `ApertureSkyScenario`,
`ApertureProvenance`, `PixelEllipse`,
`ApertureFluxRecord`,
`project_icrs_ellipse_to_native`, and
`measure_native_aperture`.

**Input contract.** `measure_native_aperture` takes a validated
immutable `NativeBandFrame` from the existing native quartet
loader, a project- and source-bound Scene/catalogue fingerprint,
the observed source-pixel-availability mask AFTER any explicit
Scene restore operations, and the separate manual-exclusion mask.
The existing `NativeBandFrame.valid` captures
image/INVVAR/NEXP/MASKBITS hard validity. No native-hard
restoration is allowed; no sky/target-protection geometry is
interpreted as a photometric exclusion.
The current `02 Scene & Masks` GUI does not yet provide a
full science-trusted membership receipt. This function
therefore reports **diagnostic valid-area flux only**.

**Projection.** The ICRS geometry stores
(ra,dec,a in arcsec,q=b/a,PA east-of-north). At the
target ICRS centre, Astropy `SkyCoord.spherical_offsets_by`
constructs local +east and +north one-arcsec celestial
offsets. The original WCS maps both onto the native
pixel basis. SVD maps the on-sky ellipse to a
locally affine pixel ellipse with the actual axis lengths,
pixel handedness and photutils CCW angle; NO assumption
that celestial PA and pixel theta coincide.
Sixteen on-sky perimeter checkpoints must agree
with local affine pixel predictions within **0.15
native pixel** (engineering tolerance, not
an external calibration). If this fails
or WCS has unsupported SIP-like distortion,
fail closed rather than advertise an exact
celestial boundary.

**Pixel integration.** Photutils
`EllipticalAperture.to_mask(method="exact")`
provides the overlap fraction of each native
pixel. Only an aperture-local cutout is
materialized; no resampling or interpolation
of source FITS occurs. The product distinguishes
analytic ideal area, in-frame exact-overlap area,
valid unmasked area, outside-frame area and
disjoint/priority-counted hard, source and
manual-exclusion areas. Hard-invalid pixels
are further partitioned exclusively into
nonfinite IMAGE, nonpositive INVVAR, invalid NEXP,
saturation, allmask, configurable contextual
flags, user-exclusion and an explicit
other-invalid category. The raw frame validity
policy is not silently strengthened or
weakened by sky fitting.

The measured finite-aperture flux is signed,
with exact fractional native pixel weights,
and the only default residual sky is
`zero_additional_correction` on the
already sky-subtracted DR11 data.
Optional project-scoped constant and tilted
planes are **conditional non-adopted
sensitivity scenarios**, whose `b_0,b_x,b_y`
coefficients and reference pixel must be
passed explicitly. Rejection and mask
areas are never used as an implicit
corrective denominator. A native circle
is an ICRS ellipse with `q=1`.

**Formal uncertainty, explicit limits.**
The diagonal-only uncertainty is
`sqrt(sum(w_geom^2 / INVVAR))` on exactly
the usable source pixels. This is NOT
calibrated total uncertainty for correlated
coadd pixels. Positive signed integrated
flux permits a **valid-area-only** AB
magnitude with zero point 22.5; negative/zero
flux remains signed and has no magnitude.
Mean surface intensity uses the local
projection-plane pixel area and is
marked approximate. Eight pixel-centre
on-sky sectors report diagnostic
in-frame/valid area but their angular
sector assignment is not exact
wedge–pixel intersection.

**Lineage and immutability.** Records include
native dataset ID, source/project/draft and
Tractor fingerprints, WCS grid fingerprint,
chosen ICRS/native geometry, all non-adopted
sky parameters, overlap method and a
SHA-256 recipe fingerprint. When the
existing validated quartet's FOUR source
FITS hashes are missing, the record
explicitly flags that it cannot
validate the FITS identity. QA fields
are frozen against in-place mutation
and `to_dict` creates a detached
RFC-8259-compatible JSON payload.
The `recipe_sha256` is an input/recipe
identity **not** an independent proof
of source authenticity or total flux.

**Deliberately not implemented yet:**
- No Stage 04 GUI, no aperture-geometry discovery
  or fully automated radial grid, no computed
  curves of growth, no independent galaxy
  source-membership decisions.
- No numerical masked-area reconstruction,
  corrected total, outer extrapolation,
  Sérsic or PSF-deconvolved model.
- No empirical coadd covariance or
  aperture-matched blank-sky calibration;
  no Science `PASS`, no all-band color
  or `TOTAL_VALIDATED` product.
- No project FITS writes or downloads.

**Regression tests:** `tests/test_aperture_flux.py`
exercises exact uniform/rotated/subpixel ellipse
area and signed flux, mutually exclusive
mask accounting, hard NEXP/IVAR defects,
negative/zero flux magnitudes, constant
and plane background algebra, WCS reflection/
rotation/RA wrap, off-frame no-flux
semantics, strict provenance/geometry and
deeply immutable QA JSON.

**Next ordered step after verified CI:**
build 04.2 fixed anchored ICRS radial grids,
direct cumulative and annular flux sums with
shared exact masks; only then stage 04.4 GUI
and later correlated-noise/membership QA.


---

## 2026-11-11 — Stage 04.2 fixed native growth and first GUI bridge

**Status:** numerical and GUI diagnostic implementation; not adopted
science and not yet executed against Leo P native data.

The new existing-package module
`src/lsbpipe/photometry/aperture_growth.py` implements
`GrowthHypothesis`, `predeclared_growth_radii`,
`geometry_hypotheses`, `measure_native_annulus` and
`measure_growth_diagnostic`. All operate on original,
unchanged band-local `NativeBandFrame` images and strictly
same-grid provisional masks. Fixed project ICRS anchor;
no auto-recentering on stellar peaks.

The sampling reference radius is only the area-equivalent
radius of provisional target protection, R=sqrt(N/pi)*local
native pixel scale. It is not a measured galaxy isophote.
When no usable protection is present, the grid uses a
clearly flagged native-field scale. Ten fixed semimajor axes
are [0.12,0.20,0.30,0.42,0.56,0.72,0.90,1.12,1.40,1.75] R.
The selection scale is capped at 85% of the shortest
native frame-edge margin, with exact area coverage retained
at every radius. A cap never proves complete coverage.

Initial geometry hypotheses are (1) project-anchor circle,
(2) explicit selected SGA q/PA hypothesis where available;
otherwise (2–3) fixed q=0.7 hypotheses at PA 0/90 degrees.
The ICRS center is invariant across all rings. SGA center
and an offset Scene isophote do not silently move the
apertures. All hypotheses are labelled nonvalidated.

Independent photutils `EllipticalAnnulus` exact fractional
native pixel integrals yield annular signed flux, usable
area, diagonal-only IVAR error, linear surface intensity
and positive-only AB surface brightness. The difference
between two measured cumulative fluxes is independently
compared with the intervening annular direct sum; it is
not assumed identical without numerical QA. Negative
annular flux remains negative and the CoG is never forced
monotonic. Masked-area reconstruction and asymptotic
extrapolation are absent.

New `src/lsbpipe/gui/aperture_workspace.py` is the fourth
real tab inside the existing `PhotometryWorkspace`.
It uses in-memory Scene original native IMAGE/INVVAR/NEXP/
MASKBITS, matching current preliminary Tractor catalogue,
revisioned ICRS mask draft, effective restored-source mask
and native hard exclusions, with the project manifest's
release/region/brick ID. Stage 03's constant-plane
preferred fit is **never** automatically adopted.
All first direct curves operate on the original DR11
coadd with ZERO additional sky correction.

The `Build aperture curves` action explicitly starts the
worker; simply visiting the tab must not read a FITS
or run any measurement. The original-source lookup never
downloads other bricks or stellar fields. A read-only
multi-geometry growth plot, radius/coverage/flux readout
and complete `Copy aperture QA` expose provenance and
limitations. Project, band, Scene revision, catalog source
and selected SGA shape changes invalidate a stale QA.
Resume stage `aperture` is supported by the existing
project-scoped navigation pointer. Offline GUI screenshots
include the Stage 04 no-native-FITS fail-closed state.

**Scientific HOLD:** source-FITS checksums are not
independently revalidated from in-memory Scene cutouts,
source membership and masked galaxy flux are unresolved,
and correlated coadd noise, sky-systematic uncertainty,
PSF wings and total magnitude are not calibrated.
Do not use available-area AB magnitudes as total
galaxy magnitudes. No observed Leo P Stage 04 flux
is claimed until the researcher runs the new stage
on its pre-existing DR11 brick.

Next: Stage 04.3 large-aperture empirical noise,
sky/geometry/source-membership sensitivity,
PSF-aware injection, then separately gated total
and multi-filter products. No new branch or
project directory is required.


---

## 2026-11-11 — Stage 04.3a conditional systematics and blank covariance

**Scope:** exploratory, project-bound *diagnostic* implemented in
`src/lsbpipe/photometry/aperture_systematics.py`; no total photometry,
noise calibration claim, sky adoption or change to existing FITS / masks.

The existing Stage 04 worker now attaches `systematics_04_3` to the
source- and Scene-revision-bound JSON. The aperture GUI shows the results
under a separate `04.3 · Conditional uncertainty diagnostics` heading
and includes full arrays in `Copy aperture QA`. It runs only on the
existing explicitly loaded native optical coadd and only when the user
requests aperture curves. Previously recorded 04.2 signed fluxes and
formally diagonal `INVVAR` errors retain their original definitions.

### 04.3a.1 Conditional sky-response derivatives

For a frozen native aperture and effective valid-pixel weight
`w_p=G_p V_p`, with reference origin at the native pixel
coordinates of its project anchor, the flux response to a *hypothetical*
residual plane `delta B_p=b0+bx(x_p-x0)+by(y_p-y0)` is

```text
dF/db0 = -sum(w)
dF/dbx = -sum(w * (x - x0))
dF/dby = -sum(w * (y - y0)).
```

These derivatives are calculated from the same exact fractional
aperture weights and mask policy as the baseline flux. They permit
direct, reproducible comparison against finite differences to the
existing `ApertureSkyScenario`. **They are not background error bars**:
a defensible distribution or bound for `b0,bx,by` has not yet been
obtained. In particular, do not set `sigma_sky` from
`RMS/sqrt(N_pixel)` or use the current Stage 03 diagnostic fit as
an automatically adopted correction.

### 04.3a.2 Mask and geometry hypotheses

At the largest common semimajor axis, the code computes a new
*counterfactual* native pixel sum with only provisional catalogue-source
exclusions disabled. The immutable `NativeBandFrame.valid` rejection
and manual exclusions are retained. The delta measures the consequence
of a different policy; **it cannot be identified with lost galaxy light**,
because it may include contaminating sources or target members.

Likewise, the observed outer-aperture flux span across the existing
predeclared circle / SGA q-PA / exploratory q-PA geometries compares
different sky footprints. It is explicitly **not a 1-sigma geometry
uncertainty**, nor can it be added in quadrature to sky/noise errors.

### 04.3a.3 Native-coadd blank annuli

`estimate_blank_annulus_covariance` creates a deterministic grid of
candidate non-overlapping native-pixel placements, with the same exact
`EllipticalAperture` / `EllipticalAnnulus` shapes as the project
apertures. Candidate centres are separated by more than the outer
diameter and cannot intersect the nominal project ellipse; any area
covered by hard-invalid pixels, current automated source exclusions or
explicit manual exclusions rejects a candidate. No pixels are
resampled, no masked pixels are filled, and the source arrays are copied
at entry.

For each accepted placement the full signed vector of annular sums is
stored in memory. Where at least `max(12,2N_rings+2)` placements are
available (with a deterministic maximum of 48), a sample covariance
matrix between **disjoint annular sums** is returned. It captures
some spatially correlated coadd noise and residual background
structure without assuming that coadd `INVVAR` is its full covariance.
Nested CoG covariance can subsequently be constructed from this matrix
via `C_CoG = T C_annulus T^T`; that transformation is **not yet an
adopted science uncertainty**.

If too few non-overlapping clean positions exist, the module returns
`INSUFFICIENT_INDEPENDENT_BLANK_POSITIONS`, a null covariance and the
explicit count — **never a zero variance or falsely precise CoG error**.
Constant or degenerate residuals likewise do not become zero-valued
science errors. Geometric exclusion of the nominal ellipse cannot rule
out diffuse galaxy light outside it, Galactic cirrus, undetected
background objects or large-scale gradients. These are limitations even
when the sample covariance is nondegenerate. The result therefore
retains `covariance_calibrated_for_science=false` and
`adopted_for_science=false`.

### Remaining 04.3 science gates

- Test sky-response uncertainties against independently constrained
  background models, including low-frequency spatial modes.
- Validate blank-position independence, sample size, local stationarity,
  and empirical covariance through repeated injections and held-out
  regions, then compare to per-pixel `INVVAR` predictions.
- Compute robust annular / cumulative joint covariance and propagate it
  through fits without double-counting one sky nuisance parameter.
- Review ambiguous source memberships (particularly the Leo P Gaia
  foreground candidates) and quantify model-conditioned masked flux.
- Verify source FITS checksums and run the native Leo P `r` measurement.
  Until then: **SCIENCE HOLD**.


---

## 2026-11-12 — One-click project-bound QA bundle

The Photometry Workspace title now exposes **Copy QA** (object name
`photometry_copy_project_qa`). Its output uses schema
`desidera.gui.project-diagnostic-bundle/1` and is JSON suitable for a
single ChatGPT message or a UTF-8 text/JSON attachment.

**Contract.** The button does not start background workers, fetch files,
read original FITS, create directories, change the Scene draft, adopt a
fit, or restore previously invalidated results. It snapshots **one
galaxy project**: current loaded Scene band plus all four independent PSF-0
band panels. Unavailable stages are reported as `NOT_RUN`,
`BLOCKED_BY_SCENE_INPUTS`, `WORKER_IN_PROGRESS`, or
`STALE_RESULT_WITHHELD`. Failure of one stage's report collection is
isolated as `COLLECTION_FAILED` with its exception class and message;
this never becomes an empty successful report.

Bundle sections:

* Runtime/package versions and best-effort Git checkout SHA (if a loose
  local Git branch pointer is present; absence is explicit).
* Project UUID, local saved dataset triplet, anchor, current GUI stage,
  navigation and loaded-band identities, draft and catalogue fingerprints,
  local product-discovery state already shown in the GUI.
* **01 PSF-0:** independent g/r/i/z current settings, source and
  recipe provenance, review/stars, science rejection/warnings, and
  non-adopted QA. Per-band errors do not block other bands.
* **02 Scene & Masks:** native FITS/WCS metadata, draft and dirty state,
  Tractor/SGA identification, source registry, candidate decisions, mask
  population counts, isophote sensitivity, full policy ablation
  *per scenario*, and leave-one-out statistics. Repeated long per-radius
  policy traces are stored as compact typed rows with every radius,
  coverage, sectors-positive and measurement status retained; complete
  octant arrays remain available in the original `Copy Scene QA`.
* **03 Background:** native coadd residual-sky diagnostics, candidate
  fit/holdout results, grid perturbations, injections, explicit
  non-adoption status. Stale results are withheld by identity check.
* **04 Aperture:** finite-area flux/area and invalidity breakdown for
  every provisional geometry/radius, signed independent annuli,
  closure QA, sky/mask/geometry sensitivity, blank annulus and CoG
  covariance diagnostics. Repeated source-ID fields are lifted to their
  containing records; per-aperture octant QA is retained.
* Last 20 operation-journal **attempts** where present (bounded, read-only,
  not proof of a completed scientific calculation), and deterministic
  `diagnostic_findings` triage codes such as absent sky QA, insufficient
  blank covariance or strong mask-policy dependence.

The bundle intentionally excludes raw FITS, masks, stamps, CCD pixels,
and any non-trivial numerical image grid. It uses strict JSON serialization
(`allow_nan=False`); NaN/Infinity values are represented as explicitly
flagged `diagnostic_nonfinite_float` objects rather than silent numeric
substitutions. Unrepresentable or broken stages are surfaced as
`COLLECTION_FAILED` with an error message. Existing
individual **Copy Scene/Background/Aperture QA** buttons and explicit PSF
export remain available for deeper source-specific investigations.

**Scientific boundary:** `science_hold=true`,
`scientific_adoption=false`, even if every GUI diagnostic is present.
Missing-band metadata, no source FITS hash, failed source membership,
unvalidated sky model/PSF and uncalibrated coadd covariance must not be
silently interpreted as a passed scientific acceptance gate.


---

## 2026-11-12 — Real Leo P r audit and follow-up diagnostic contracts

The first all-stage, native-`r` `Copy QA` from Leo P shows that
**astrometrically supported foreground sources inside protected target
light** are a first-order photometric limitation. The existing automatic
source screening correctly prevents blanket mask removal of protected
pixels; however, direct finite-aperture flux includes those unmasked
stellar cores. The project report must make this consequence unmistakable.

### Source evidence × circular annular flux (read-only)

The all-stage bundle now includes
`stages.foreground_aperture_overlap`. It joins **current matching**
Scene inventory records with the native, neutral-circle Stage 04 growth
curve. Only Tractor `PSF` sources with Gaia DR3 reference IDs, finite
positive Tractor model-band fluxes, and
`abs(parallax_significance) >= 5` are listed. The source's spherical
angular offset from the saved project target selects the enclosing
neutral circular ring. For each such source, retain its Tractor
OBJID/Gaia ID, parallax, offset, protected-area flag, band model flux,
measured signed ring flux, and ratio of model flux to signed ring flux.

The **ratio is a heuristic red flag, NOT the flux fraction attributable
to that source**. The catalog model and observed rings have different
footprints, the sky/background was not calibrated, and the band PSF
wings have not been scientifically validated. No star is automatically
masked or subtracted. Multiple sources and galaxy light may contribute
to the same ring. If a protected source's model/ring ratio exceeds 0.25,
the report adds the explicit triage
`PROTECTED_BRIGHT_GAIA_SOURCE_MAY_DOMINATE_APERTURE_RING`.
This is conditional on same-band, same-source provenance; absent or
stale inputs fail closed.

For the initial Leo P `r` QA, OBJID 2035 is at 19.87 arcsec with
Gaia parallax 1.088 mas (30.7-sigma relative to reported uncertainty).
Its Tractor `r` model flux is ~1078.75 nMgy, while the fixed
18--25.2 arcsec neutral ring has signed flux ~1126.17 nMgy.
The new report identifies this pairing automatically, without
asserting that their near coincidence establishes a valid
decontaminated Leo P magnitude.

### Stage 03 sky candidate × exact Stage 04 response (read-only)

The all-stage bundle additionally evaluates each *fitted but
unadopted* Stage 03 constant/plane model through the exact Stage 04
sky-response derivatives, transforming the plane's reference origin to
the aperture's target anchor. Outputs are stored as
`stages.conditional_sky_aperture_response` with explicit
`correction_applied=false`,
`sky_uncertainty_inferred=false`,
`science_adopted=false`. These are hypothetical flux shifts if that
model were subtracted, NOT measured sky bias or a calibrated error
bar. Project, band, draft, catalogue and native source identities must
match. Missing or stale states remain explicit.

### Smaller-scale native-coadd blank diagnostics

A single 105-arcsec diameter-scale blank annulus rarely fits in
the 4-arcmin Leo P cutout while remaining independent of the target.
Stage 04.3 now tests a bounded, predeclared decreasing set of radial
prefixes if the full-scale diagnostic is unavailable. **The target
exclusion for every such probe uses the full 105-arcsec galaxy aperture,
not the smaller probe radius.** A conservative centre-distance veto
excludes overlap between the full target ellipse and smaller blank
ellipse. Candidate blank placements must still be disjoint, and
must remain entirely on unmasked, valid native pixels.

The result `smaller_scale_blank_diagnostics` reports the tested
radii, number of accepted placements and sample covariance, if any,
with `must_not_replace_full_scale_covariance=true`. No smaller-scale
covariance is allowed to stand in for the full CoG covariance;
`covariance_calibrated_for_science` remains false. No PSF wings, target
outskirts, spatial nonstationarity, or historical coadd sky subtraction
are certified by this exercise.

**Next scientific actions**: perform the existing GUI
`Probe selected Gaia star · saved PSF-0` on OBJID 2035 and validate
a foreground star wing model; independently test the unadopted residual
sky planes in repeated blank regions; develop membership-aware, PSF
conditioned star handling before interpreting any apparent CoG plateau.
Scientific hold remains mandatory.


---

## 2026-11-12 — End-to-end workflow and scientific gate audit

**Scope:** native DESI Legacy DR11-South Leo P project, Stage01 PSF-0,
Stage02 Scene/Masks, Stage03 Background, Stage04 finite aperture,
and the all-stage Copy QA collector. The baseline real-`r` bundle at
source revision `11513b7a` is an experimental diagnostic and remains
unadopted. All statements below refer to this scope, not to independent
raw-exposure reprocessing or a validated `g/r/i/z` catalog.

### Correct operator order and distinct meanings

1. **Project identity:** open the saved Leo P project; inspect the DR11
   source triplet and target ICRS anchor. GUI file-presence discovery
   is not FITS/WCS or SHA verification. The selected manifest source
   band displayed in the source table may be `g` while Stage02--04
   deliberately operate on a separately loaded, verified native `r`.
2. **PSF-0:** a per-band, independent, read-only/saved *finite-kernel*
   diagnostic. A saved per-band kernel is a prerequisite for the
   explicit Stage02 bright-star PSF probe but **not** a prerequisite
   for Stage03 sky candidate diagnostics or the Stage04 zero-additional-
   correction flux calculation. Review source hashes, star counts,
   holdout spatial coverage and finite kernel support; never assume
   PSF wings are measured. The real Leo P report had saved `r/i/z`
   results and `g: NOT_RUN`; “build available” does not guarantee `griz`.
3. **Stage02 Scene:** choose the exact same native `r` band explicitly,
   load its four FITS/WCS layers and prepare a fresh, same-band
   catalogue mask. Target protection is a source-screening safeguard,
   **not** a photometry exclusion. Baseline experimental Gaia/red
   checkboxes are OFF. Only now run a conditional isophote trial,
   40--80 arcsec source-protection radius sweep or 70-arcsec policy
   ablation. Those experiment radii refer to *distinct purposes* and
   must never be treated as a measured physical galaxy extent.
4. **Stage02 source review:** inspect bright Gaia foreground candidates
   within protected target area, especially OBJID 2035. The saved PSF
   foreground probe must be run **after** final isophote proposal and
   all desired source-mask experiments: re-preparing masks, loading
   a different band or clearing/rebuilding an isophote can invalidate
   its result. Do not derive a star subtraction or science mask directly
   from a Tractor `PSF` model flux or the finite-kernel probe.
5. **Stage03:** fit the zero, constant and tilted-plane *residual sky
   candidates* from the unchanged, source-bound native scene; verify
   independent tile spatial placement as far as the same coadd allows.
   Avoid presenting a held-out same-coadd pixel as independent exposure
   validation, or a 25%-RMS model-selection rule as a significance test.
   The fitter does not adopt or subtract any of its candidates.
6. **Stage04:** measure fixed project-anchored geometric hypotheses
   on original native, already upstream-sky-subtracted coadd pixels.
   The current Stage04 uses **zero additional sky correction** and
   does **not** consume a Stage03 fitted background or Stage01 PSF
   correction. Stage04's exact fractional per-pixel overlap validates
   geometric integration; independent ring-vs-cumulative closure is
   internal numerical consistency, not a check of foreground source
   membership, absolute sky level or total flux.
7. **Copy QA:** copy once while still on the loaded `r` scene and
   after completed workers; its complete in-memory content may be
   several MB and is best shared as a JSON attachment. It marks
   missing/blocked/stale computations, but should not be interpreted
   as an automatic run of all photometric prerequisites. Switching
   to another Scene band invalidates that band's native masks, sky and
   aperture results; make a new report after explicitly reloading
   and re-running for `i`.
8. **Any science use:** STILL BLOCKED pending PSF-wing/foreground
   treatment, source membership model, residual-sky uncertainty and
   subtraction policy, multi-scale correlated-noise covariance,
   flux under masked regions, full input file hashes and
   matched-footprint/matched-PSF multi-band photometry.

### Audit findings and dispositions

* **FIXED (code-level):** Stage04 incorrectly read
  `SceneEditor.sga_use_prior` (target protection) to enable an
  SGA geometry hypothesis. It now reads `sga_use_seed` (shape seed).
  Synthetic GUI regression covers both mutually independent controls.
* **FIXED (sequence/GUI):** leaving PSF-0 after inspecting another
  band could silently reselect that band in an already loaded Scene
  and invalidate masks/sky/aperture. Existing loaded Scene is now
  preserved; choose the Scene band explicitly. PSF still initializes
  an empty Scene's band on first entry. A GUI transition regression
  covers both paths.
* **FIXED (QA evidence):** a negative but formally significant Gaia
  parallax is no longer described as positive foreground-distance
  evidence. The new annular bright-star QA requires a positive
  significant parallax. Sources with proper-motion-only evidence are
  not exhaustively enumerated by that subsection; the complete
  Tractor source inventory remains available.
* **INTENTIONAL HOLD (not a code bug):** Stage01 and Stage03 are not
  currently science-adopted dependencies of Stage04. Cross-stage
  sky-response/foreground-overlap reports are conditional consistency
  checks and apply no correction.
* **SCIENTIFIC BLOCKER:** sky training is spatially nonuniform (24
  training vs 7 holdout tiles on the Leo P r cutout; training RMS
  exceeds holdout by >5x). Model-selection differences in nMgy/pixel
  are coherent offsets that scale with the large source-unmasked
  aperture area; diagonal IVAR alone understates uncertainty.
* **SCIENTIFIC BLOCKER:** OBJID 2035 is a high-significance Gaia
  foreground candidate at about 19.87 arcsec. Its Tractor r-band
  model flux (~1078.75 nMgy) is of the same order as the directly
  observed 18--25.2 arcsec ring flux (~1126.17 nMgy). Equality of
  these values is *not* an independently fitted star contamination
  estimate and does not validate an adjusted galaxy flux.
* **SCIENTIFIC BLOCKER:** full-scale blank-aperture covariance is
  unavailable in the current image cutout. Shorter-scale background
  probes are diagnostics at the reported smaller radii, not covariance
  for the 105-arcsec integrated aperture.
* **REPRODUCIBILITY LIMIT:** Stage01 saved PSF includes source layer
  hashes; the Stage02--04 in-memory native cutout still reports no
  verified full source FITS SHA256. The source table is only a
  last-discovery path-presence report, not a hash or WCS cross-check.

**Science release rule:** no total Leo P magnitude, uncertainty,
color, `SBF` measurement or adopted model subtraction may follow from
these diagnostics. A green CI run checks programmed contracts, not
photometric validity on the sky.


### Additional coordinate-frame audit: selected isophote centre versus project anchor

An operator may pick a centre for Stage02 *isophote detection* that is
different from the saved target ICRS anchor. The Stage02 source inventory
can then contain angular `offset_arcsec` relative to this picked
detection centre, whereas Stage04's fixed neutral circular apertures
always use the **saved project target anchor**. Reusing a Stage02
inventory offset to assign a source to a Stage04 annulus is therefore
incorrect and could yield a false contamination diagnosis.

**Fixed:** the all-stage foreground/annulus QA now independently
computes the source great-circle separation from the exact
`Stage04.target_anchor_icrs` using its catalogued ICRS RA/Dec. The
Stage02 inventory's original offset is retained *separately* and
labelled potentially different-anchor. Missing fixed target centre
fails closed. An explicit regression supplies contradictory centres and
verifies that the project-anchored radial bin is chosen.


---

## 2026-11-12 — Deep follow-on code audit (Stage01–04)

This follow-on review inspected the actual backends and regression contracts,
not only the displayed pilot QA. Distinguish **confirmed code-level
defects** below from **still-unvalidated scientific assumptions**. All
fixes remain in Draft PR #84; no source FITS, existing user mask drafts,
or `main` branch are modified.

### Confirmed defects corrected

**Stage01 PSF target geometry / lineage (P1).** The GUI backend first
verified that a Leo P project centre lay within 45 arcsec of the canonical
Leo P coordinate, but then supplied the *canonical coordinate* to the
PSF-star target-exclusion step. Stage02/03/04 use the saved project
centre. A user moving the project centre by, say, 20 arcsec without
changing the brick would have obtained inconsistent star-selection
protection between stages. PSF build now uses the saved
`project.field.ra_deg, dec_deg`; the canonical coordinate remains only
the Leo P identity gate. Regression explicitly creates a shifted but
eligible Leo P project and checks the actual PSF fit's input anchor.

**Stage01 saved PSF stale-target guard (P1).** Existing PSF checkpoint
restoration verified the FITS quartet, Tractor catalogue, kernel checksum,
source brick/band and recipe, but did not verify that the current saved
project *target centre* was the one used to train the PSF. A moved
project could therefore restore a result protected around a different
centre. Explicit save and read-only restore now require the original
`target_icrs_deg` to agree with the current project field within
`1e-9 deg` per coordinate (wrapped RA). A mismatched saved result is
withheld, recorded in restore warnings, not deleted. Synthetic FITS
regression covers a successfully saved PSF followed by project-centre
editing, source-hash invariance, and rejection at both boundaries.

**Stage02 circular Tractor mask WCS (P1 for anisotropic fields).** The
point-source and degenerate extended-source circular footprints had used
the *geometric mean* pixel scale, even when the native WCS is
anisotropic/sheared. Elliptic catalogue footprints already used the
local full east/north pixel Jacobian. Consequently the circular
footprint could leave contaminating pixels unmasked on the fine-scale
axis. Circular footprints with a celestial WCS now rasterize the
**on-sky circle** in the local east/north tangent plane and bound its
native pixel box with the minimum singular value. WCS-less legacy
diagnostics retain the explicitly approximate scalar fallback.
Synthetic anisotropic-WCS regression exercises the previously missed
region. This fixes geometric coverage, **not** membership, appropriate
PSF wing-mask radii or missing galaxy-light reconstruction.

**Stage03 sky spatial QA tile coordinates (P1 for shifted grids).**
The selector saves native `tile_bounds_native_xyxy` for each sampled
spatial tile. The residual audit incorrectly regenerated bounds
`tile_index * tile_edge`, valid only for an unshifted grid. Since fixed
grid-offset perturbations are supported, this could sample other pixels
or throw a false `Sky sample tile geometry/provenance changed`
exception when used on shifted selections. The QA now reads saved
actual bounds and verifies they lie inside the raster. Edge-tile
counting uses actual dimensions. Tests cover a nonzero offset,
separate train/holdout medians and mismatch rejection.

### Scientific and operational limits, not silently fixed

* **PSF-0 is currently a Leo P-only pilot.** `_leo_p_project`
  explicitly rejects other galaxy names and projects >45 arcsec from the
  hardcoded Leo P centre. Unlike Stage02–04, it is *not a general
  per-object service*. Its finite kernel cannot establish far stellar
  wings or spatial PSF variation. Generalization requires its own
  explicit project-target identity and acceptance workflow.
* **Different target protections have distinct purposes:** PSF-0's
  stellar-training exclusion radius is user-adjustable and need not
  equal Stage02's source-screening protection or Stage04's last
  measured aperture. A larger PSF radius cannot be interpreted as a
  measured photometric galaxy radius.
* **Stage02 catalogue-source footprints are heuristic and source-type
  dependent.** They are neither fitted Tractor model residuals nor
  all-flux wings. Even with geometrically correct footprints, their
  science use needs contamination/galaxy-light-loss experiments and
  high-confidence foreground/target membership separation.
* **Stage03 holdout is only spatially disjoint within one resampled
  coadd.** It is not an independent exposure, covariance model, nor
  proof that its sky tiles are galaxy-free. The nominal same-coadd
  RMS-based constant/plane preference is not a science decision.
  Grid perturbations remain uncalibrated sensitivity experiments.
* **Stage04 zero-additional-correction is intentionally not coupled
  to Stage03's sky fit**, and is independent of whether Stage01 has
  a reconstructed PSF. Valid-area flux with some missing pixels is not
  a full galaxy flux. Its formal diagonal IVAR sigma excludes
  correlated coadd noise and covariance from foreground subtraction,
  masked-pixel inference and sky estimation.
* **Cross-band `r-i` cannot be adopted** until identical apertures,
  compatible PSF response, consistent target/contaminant policy and
  jointly calibrated covariance/sky corrections have been demonstrated.
  Current radius grids depend on the provisional target-protection
  pixel count and finite cutout extent; separate bands can therefore
  differ in the actual sampled radii.
* **Native input identity still differs between PSF0 and Stage02–04:**
  Stage01 hashes all four original native FITS input layers; Stage02
  checks native FITS/WCS cutout metadata but does not currently expose
  complete full-file SHA256 provenance to Stage04. Do not interpret
  the all-stage QA's `source` label as equivalent to a verified
  full-input content digest.
* **The imaging PSF and sky systematics are coupled.** A PSF-star
  kernel fitted with a contaminated local annulus cannot automatically
  be treated as independent evidence for large-radius galaxy
  photometry or foreground-star subtraction. Dedicated wing,
  injection and joint star/galaxy/sky validation are science gates.

### High-value next independent tests

1. On the same native brick and frozen Scene, compare modelled stellar
   subtraction versus conservative masked-pixel photometry for bright
   Gaia foreground candidates (especially Leo P OBJID 2035), with
   heldout stellar wings and galaxy-light-loss injection tests.
2. Acquire sufficiently wide native context or independent overlapping
   sky fields to estimate large-scale background structure and the
   full-aperture coadd covariance. Do not bootstrap 105-arcsec covariance
   from 10–30-arcsec probes without externally validated scaling.
3. Run a common-aperture, common-geometry, PSF-conditioned
   `griz` synthetic challenge (known injected galaxy flux, realistic
   contaminants, gradients and masks) before reporting AB totals,
   colours or SBF distances.
4. Introduce an explicit source-manifest hash gate for Stage02–04 and
   per-band run provenance before durable science-result exports.

**Scientific acceptance status remains HOLD.** A passing test suite
validates the currently encoded contracts and tested synthetic cases;
it is not validation against independent external galaxy photometry.


---

## 2026-11-12 — Original-FITS provenance, unit, and cross-stage lineage gates

A further repository audit identified an actionable gap between the
full-brick verified PSF-0 lineage and the cutout-based
Scene/Background/Aperture pipeline.

### New enforced contract

1. On explicit Stage02 `Load original FITS`, the worker resolves the
   manifest's exact `(release, region, brick, band)` quartet from its
   existing data root. It records pre-load file sizes and nanosecond
   modification times; loads the original cutout using the existing
   strict WCS-aligned, no-padding renderer; and then hashes each of
   IMAGE, INVVAR, NEXP and MASKBITS as original **whole files**, not
   displayed cutout arrays. A file replacement during the load causes
   the result to be withheld. Files are only read, never modified.

2. The *same original IMAGE science HDU* that supplied the pixels must
   have a nonempty recognized nanomaggy `BUNIT`; an absent or explicitly
   incompatible IMAGE BUNIT fails closed before photometry. An
   explicitly incompatible INVVAR BUNIT is also refused. A missing
   INVVAR BUNIT retains the specifically labeled
   `INVVAR_UNITS_FROM_SURVEY_DOCUMENTATION_ASSUMPTION`, **not** a
   measured unit confirmation. The accepted units match the
   established native-PSF frame loader contract.

3. Stage02's source mask preparation checks the loaded native
   `source_quartet_key` against the selected project manifest before
   reading a catalogue; its output carries the same machine key.
   Stage03 refuses mismatched native-image/Tractor-mask identities.
   Stage04 refuses mismatches with the project source triplet, native
   Scene identity, or prepared masks. The GUI's human-readable source
   label is not trusted as the primary identity.

4. Stage03 and Stage04 expose the verified original SHA-256 fingerprints,
   input unit evidence, and source key when present. Stage04's
   exact-aperture recipe hash includes the original file hashes, making
   changes in source bytes alter the derivative diagnostic identity even
   with the same pixels, masked geometry, and object name. An
   inconsistent or unverified hash claim fails closed. Synthetic
   callers without full source file hashes remain explicitly
   `NOT_VERIFIED`; they do not receive science approval.

5. The unified `Copy QA` additionally emits
   `stages.cross_stage_source_lineage`, checking that the current
   background and aperture outputs were produced using matching
   original file hashes, native quartet identity, and IMAGE unit
   evidence. It marks missing or inconsistent evidence rather than
   retrospectively declaring data accepted.

### Limits and performance

* Hashes establish **local byte identity at the time Scene was loaded**.
  They are not independent archive authenticity and do not certify
  upstream DR11 calibration. A file modified *after* loading cannot
  retroactively alter already copied in-memory pixels; a new Scene
  load is required to use it.
* Full-file hashing entails additional read I/O on explicit Stage02
  loading. It is never triggered by merely opening a tab, creating a
  `Copy QA` report, fitting the background or measuring apertures.
* FITS header and SHA validation do **not** validate astrophysical
  sky subtraction, PSF wings, coadd covariance, masked light, or
  galaxy membership. `SCIENCE HOLD` remains unchanged.
* The current Scene viewer retains the existing detailed native
  image/WCS compatibility gate. Header-unit checks act as an
  additional science-readiness veto and are not a replacement for
  full independent survey archive verification.

### Added synthetic regressions

* Four FITS headers with explicit BUNIT and local SHA-256 in a synthetic
  project, followed by changed IMAGE units and mid-load source-file
  replacement; ensure missing/invalid provenance never reaches GUI.
* Exact-aperture SHA receipt changes when a verified original image
  hash changes without any modification to the pixel integration.
* Explicit mismatch between loaded native source key, current project
  catalogue key and Stage04 source triplet is refused.
* Stage03 rejects mismatched catalogue identities and unverified
  asserted source hashes.
* A one-click QA source-lineage comparison detects a changed aperture
  source fingerprint while preserving other stages and a global
  `SCIENCE HOLD` flag.

This is a **code and synthetic-regression audit**. Real Leo P FITS
have not been reloaded in this environment; the new header/hashing
contract still needs a local GUI acceptance run on the scientist's
original DR11 quartet.


---

## 2026-11-12 — Project-anchor invalidation and blank-covariance budget audit

### Confirmed project lifecycle defect: same UUID, moved target

The Stage02 GUI previously memoized project/source identity from
`(project path, UUID, selected quartet, data root)`. Editing the
project's saved ICRS target centre or field size while retaining those
four attributes then rebinding the project could **keep old in-memory
native pixels and source masks but update the GUI target centre**.
Stage04 uses the target centre to project its exact apertures, making
the mixed old-pixel/new-centre state scientifically dangerous even if
`Copy QA` reports a matching project UUID.

**Fixed:** the GUI-bound source key now additionally contains the exact
saved project RA, Dec and field width. Rebinding with changed field
geometry invalidates the loaded Scene, prepared masks, isophote
diagnostics and downstream results. Each real Scene loader output also
carries `loaded_target_icrs` and the source quartet key. The main
thread rejects asynchronous worker results for a different target
centre or brick (including cases where the project UUID and displayed
band agree). The Stage03 snapshot refuses a mismatch between the
current project centre, loaded Scene origin and prepared mask quartet;
Stage04 independently refuses mismatch between its aperture anchor
and the loaded Scene centre. The unified QA now compares the *fixed
project centre* in Scene, Background and Aperture, independently of
an optional operator-picked isophote centre.

The tests deliberately move the target without changing the UUID,
present an obsolete worker result, swap source brick IDs and corrupt
mask lineage. The affected work is blocked or invalidated, never
silently reused. All these controls are non-destructive.

### Confirmed uncertainty-QA sampling-budget defect

The empirical blank-annulus covariance routine allowed up to 40 radial
bins and required at least `2 * n_bins + 2` independent blank
placements, but its fixed default ceiling was 48 positions. At 31
bins, for example, 64 placements were required while at most 48 could
be sampled, guaranteeing a false `INSUFFICIENT_INDEPENDENT_BLANK_POSITIONS`
conclusion regardless of available image area.

**Fixed:** the default ceiling is now adaptively increased to at least
the declared requirement, within the bounded limit of 256 placements.
An explicitly requested lower `max_placements` remains respected and
produces a separate
`EXPLICIT_PLACEMENT_CAP_BELOW_REQUIRED_COVARIANCE_SAMPLE`
status rather than being misreported as lack of clean sky. The report
records `n_max_placements` and whether the ceiling was explicit.
Regression tests cover 31-ring vectors and the explicit impossible
cap. These are *diagnostic* independent placements on a single
correlated coadd, not a calibrated full-galaxy flux covariance.

### Still open

The user must still validate new source/anchor and FITS-header guards
against real original DR11 Leo P files. A synthetic pass is not a
measurement of PSF far wings, correlated noise covariance, bias from
upstream sky subtraction, masked galaxy-light restoration, or source
membership. The full science release remains blocked.


---

## 2026-11-12 — Deep audit: sky-grid references, contamination vetoes and full-source comparisons

This pass extends the verified scientific GUI contract. It does not
alter original FITS, persistent user masks or the `main` branch.

### Stage03 offset-grid control bug

The Stage03 spatial tile QA was already corrected to use the recorded
`tile_bounds_native_xyxy`. The *separate zero-additional-sky
baseline* still reconstructed a tile footprint as
`tile_index * tile_edge`. These are different native pixel regions
when the intentionally predeclared sky grid is offset. That could
produce an inconsistent sample count or a misleading zero-baseline
tile median and thus undermine comparisons with fitted constant or
planar sky models.

**Fixed:** the zero-baseline now measures each frozen held-out
footprint using its true recorded native bounds, validates its
in-image rectangle and exact sample count, and explicitly fails
closed on mismatches. Regression supplies a deliberately displaced
tile and known signed native pixel residuals. The zero baseline still
does **not** establish a correct absolute sky level or authorize
an additional subtraction.

### Stage03 target-protection precondition

With zero protected target pixels, the original
`_automatic_scene_sky_pixels` function could sample the target as
putative sky if other masks happened to allow it. It could then fit
a background against galaxy light.

**Fixed:** a missing or wholly empty target-protection footprint
raises an explicit background selection failure. The original
science arrays remain unchanged; no candidate is adopted.
This is a *necessary* precondition, not a sufficient test that
undetected galaxy outskirts are excluded.

### Stage03 Gaia angular veto with anisotropic WCS

The independent bright-Gaia-star sky veto originally interpreted
arcsecond exclusion radii as circles in pixel coordinates using
a geometric-mean pixel scale. For anisotropic/sheared WCS, this
under-excluded the fine-angular-scale native direction. The same
scale was also used for pixel-space dilations around the target
and catalogue sources, making them too narrow along a fine axis.

**Fixed:** each Gaia bright-star veto is now a local tangent-plane
sky circle rasterized through the full 2×2 native WCS Jacobian.
Distance-transform dilations conservatively use the smallest
angular scale per native pixel, never the geometric mean. The
explicit report retains both scales and identifies the veto as
an *unvalidated* PSF-wing-radius heuristic. Synthetic anisotropic
and degenerate-WCS regressions exercise both footprints. This
does not measure true bright-star wing extents or certify sky
contamination.

### Cross-stage comparison gates

Earlier Copy QA versions computed a numerical response of Stage04
aperture flux to Stage03 unadopted sky fits when project UUID,
band, draft and catalogue hashes matched. This was insufficient
if an older native FITS/WCS, copied catalogue identity, or moved
project centre survived in a numerical output.

**Fixed:** the conditional sky/aperture comparison requires the
*current Scene*, identical machine source quartet, full original
local SHA-256 quartet, native WCS/shape grid fingerprint, IMAGE
nanomaggy unit evidence and fixed target ICRS centre in Stage03
and Stage04. The underlying stages must explicitly report
Scene-load original-file verification, and their aperture target
anchor must agree with the loaded Scene. Missing or contradictory
evidence yields an explicit provenance-mismatch state, not a
numerical flux response.

The independent `cross_stage_source_lineage` report likewise
validates the presence and hexadecimal syntax of each original
SHA-256 and compares the native grid fingerprint, source quartet,
unit evidence and fixed project anchor across current results.
WCS fingerprint is captured at Scene load without rereading
original FITS when Copy QA is pressed.

### Remaining scientific holds

All per-band aperture fluxes are still *partial direct signed
measurements on unmasked native coadd pixels*. PSF-0 is not yet
a validated galaxy/foreground deblender; fitting a background
candidate does not approve subtraction; neither masked-galaxy
light nor large-scale correlated coadd covariance is recovered.
No total Leo P magnitude, colour, asymptotic result or SBF
measurement is released from these QA-only controls.

A green CI run is evidence for code-level contracts in synthetic
and GUI cases. Repeat source/header/PSF/sky controls on the
researcher's local original DR11 Leo P FITS and compare with
independent astrophysical truth before changing `SCIENCE HOLD`.


---

## 2026-11-12 — Deep audit: empirical blank covariance sampling fairness

The Stage04.3 empirical covariance diagnostic originally enumerated
candidate empty apertures in increasing native image rows and columns,
accepting the first `max_placements` clean positions. Thus even a wide
source-free original coadd could yield all 48 accepted blanks close to
one side of the image. A sample covariance computed from that cluster
cannot represent broad-scale field variation. This was a **sampling
design defect**, not an error in the matrix formula
`C_cumulative = T C_annulus T.T`.

**Corrected**: the same predeclared, nonoverlapping native-pixel lattice
is divided into four fixed image quadrants. Within each quadrant its
candidates are permuted with an explicitly recorded deterministic
seed, independently of coadd fluxes; the search cycles through
quadrants before reaching the existing bound on accepted positions or
evaluations. The report now carries `candidate_grid_nodes_by_quadrant`,
`accepted_centres_by_quadrant`,
`n_quadrants_with_accepted_blanks`,
`candidate_order_independent_of_image_brightness`, and the new
placement-policy identifier V2. An actually restricted field still
reports insufficient spatial support; a one-quadrant or strongly
one-sided accepted sample emits explicit warnings. Changing pixel
intensity alone never changes the selected candidate positions if
validity/masking stays fixed.

Dedicated regression checks require all four quadrants to be sampled
on an ample clean image, verify reproducibility and pixel-intensity
independence of the candidate selection, confirm physical separation
between accepted apertures, and flag a contrived image with clean
sky confined to one quadrant.

**Scientific limit**: deterministic quadrants improve selection
representativeness but cannot prove that blank positions are
independent noise realizations: they share the same resampled coadd,
possibly common background structures and unmasked LSB signal. The
calculated matrix remains a **diagnostic only**. For Leo P a 105-arcsec
target aperture in a roughly 4-arcmin cutout still has too few clean,
nonoverlapping full-size blanks; finer-scale tests MUST NOT be
extrapolated to the 105-arcsec flux without external calibration.
The numeric covariance from V2 is generally *different* from V1
by design; this is not a new flux calibration or science acceptance.


---

## 2026-11-12 — Deep projection-plane and source-mask recipe audit

### Confirmed high-declination Stage02 TAN-WCS defect

The provisional Tractor catalogue mask rasterizer previously used
`WCS.celestial.pixel_scale_matrix` (the FITS CD/PC projection-plane
matrix, in degrees per native pixel), then multiplied its **first
row by `cos(CRVAL2)`** before interpreting the result as east/north
angles. In a TAN projection this is a second, erroneous declination
correction: the CD matrix is already in intermediate tangent-plane
degrees, NOT raw coordinate differences in RA and Dec.

This error stretched catalogued PSF and extended-source masks along
the first projection direction by a factor approximately
`1 / cos(CRVAL2)` and also underestimated local sky-plane distances
for the opt-in isolated-red-extended-source classifier. At a field
declination of 75 degrees this factor is ~3.86. The originally
tested Leo P geometry at +18 degrees suffers a smaller (~5%)
but still genuine difference in source-mask footprint.

**Fixed:** `conservative_catalogue_mask` now uses the 2×2 native
CD/PC intermediate-plane matrix directly in arcsec/pixel, with no
extra cosine. All affected circular and elliptical mask footprints
and WCS-based nearest-neighbour classifications share this corrected
local affine geometry. A synthetic high-declination FITS TAN-WCS
regression uses actual `SkyCoord` great-circle separation to show
that a point 12 true arcsec from a 9.5-arcsec provisional star-mask
centre must not be masked.

The prepared Scene source mask now explicitly reports the
`PROVISIONAL_SCENE_EFFECTIVE_AUTO_PLUS_EXCLUDE_V2_TAN_WCS` policy
and `source_mask_native_geometry` descriptor. Stage04 inherits
the actual Scene mask-policy ID into `ApertureProvenance` and
each exact-aperture recipe SHA-256. Even if an old and new policy
happen to produce identical Boolean arrays, their diagnostic
recipes differ. Preexisting synthetic callers without the
version field are explicitly labelled legacy rather than
misrepresented as fresh V2.

**Scope:** this corrects TAN projection-plane geometry near the
source selected native brick. It does not claim accuracy of a
single constant Jacobian across very wide or nonlinear celestial
fields, does not validate Tractor model sizes, Gaia membership,
stellar PSF wings, or missing target-light correction. Reprepare
Scene masks after updating the local GUI before comparing future
Leo P radial fluxes. Different values produced by the corrected
mask geometry must be treated as a new diagnostic recipe, not
astrophysical variability.

### Independent Stage04 sampling defect

The empirical blank-annulus covariance used to accept candidate
positions in monotonically increasing image-row order until the
placement cap. It could therefore fill its entire budget near one
edge even when all four field quadrants contained clean sky.
Its formal sample matrix was mathematically defined but spatially
unrepresentative.

The default candidate order is now deterministic, flux-independent,
and interleaved by native image quadrant, with a recorded seed.
The underlying separation lattice is unchanged. QA reports accepted
counts in all four quadrants and explicitly flags one-quadrant or
strongly one-sided samples. A synthetic source-free field verifies
four-quadrant coverage and pixel-flux independence; another synthetic
field permits blank positions only in one quadrant and checks
warnings. This remains a same-coadd diagnostic covariance, not a
scientifically adopted error model.

**Acceptance gate:** these regressions establish code-level geometry
and reproducibility. No new published Leo P magnitudes, colours,
SBF distances, or adopted noise covariance are inferred; all
real-data scientific results remain `SCIENCE HOLD`.


### Mask geometry provenance must be a cross-stage identity, not just a text label

The Stage02 raster algorithm V2 is now recorded directly in each
current prepared `Scene.contaminant_screening.source_mask_policy_id`.
Stage03 retains this same ID in its residual-sky result; Stage04
inherits it into both the top-level `source_policy_id` and every
`ApertureProvenance`/recipe SHA. `Copy QA` verifies that these agree
for current stages, and withholds any conditional Stage03-sky
counterfactual if a legacy V1 mask is mixed with corrected TAN V2
masks, even when the quartet of original FITS hashes, WCS, target
centre and catalogue SHA happen to match. Regression demonstrates
that inconsistent policies generate a distinct failure state, and
that no numerical hypothetical correction is emitted.

The version is scientific-input provenance: matching versions does
not imply valid foreground membership, accurate outer PSF wings,
or calibrated masked-light reconstruction. `SCIENCE HOLD` remains.


---

## 2026-11-12 — Deep audit: bounded blank-covariance search and exact mask lineage

### Stage04.3 candidate-list memory scaling

The previous V2 stratified sampling fix correctly avoided the row-major
selection bias, but it constructed **all native lattice point tuples**
for every quadrant before applying the bounded scientific evaluation
budget. In a sufficiently large image, millions of candidate centres
could be allocated solely to select at most 256 independent blanks.
This broke the stated resource-bounded design of the aperture QA
worker; CPU/memory costs scaled with possible grid centres rather
than the intended numerical evaluation budget.

**Fix (V3)**: the identical staggered, geometry-only spacing lattice
is represented by up to two arithmetic-progression blocks per native
grid row. The four quadrant totals are computed without materializing
centre coordinates. A recorded deterministic seeded, coprime-stride
index permutation enumerates each quadrant's nodes without duplicates
or a full permutation array, and the consumer stops at the existing
`max_placements` / `candidate_budget`. This reduces candidate
storage from O(number of grid nodes) to O(number of lattice rows),
while preserving bounded candidate evaluations, separation
requirements and geometry-only selection. The placement policy
identifier is versioned as
`STRATIFIED_QUADRANT_BOUNDED_MEMORY_NATIVE_RINGS_V3`; fitted blank
covariance results may change because sample ordering has changed
and must never be treated as adopted science uncertainties.

A deterministic synthetic regression uses a grid representing **over
10 million** possible lattice centres and checks the first 80
generated positions, uniqueness, native bounds and repeatability.
No full survey-science pixel array is allocated for this regression.
The existing 4-quadrant coverage, intensity-independence,
separation and one-quadrant warning tests remain applicable.

### Exact pixel-level source-mask consistency across stages

The draft SHA, catalogue SHA, WCS fingerprint and policy version
identify recipes and inputs but previously did not bind the **actual
rasterized Boolean pixels** consumed by Stage03 and Stage04. An
in-memory mask replacement or a code-path discrepancy could therefore
pass the identity comparison while sky and aperture diagnostics used
different exclusion footprints.

**Fix**: a shared canonical `native_boolean_mask_sha256` hashes
2D Boolean source and manual-exclude masks using explicit native
shape and C-ordered pixel bytes. Stage02 records both receipts
immediately after preparing the real source and manual masks. Stage03
and Stage04 recompute them from their evaluated working masks, and
fail closed before any sky fit or aperture flux if a known receipt
differs. Both results report these exact digests. The one-click
`Copy QA` now tests the same receipts across Scene, Background and
Aperture, and blocks the hypothetical sky/aperture flux response if
the actual raster footprints differ. Native invalidity remains
immutable and is not treated as manually restorable.

Tests cover one deliberately flipped automatic mask pixel, a wrong
manual-mask receipt, invariance to memory strides, and a deliberate
cross-stage source-mask mismatch despite matching project/catalogue
identifiers. Legacy synthetic paths without a prepared mask receipt
remain explicitly outside pixel-level verification; their absence
is not evidence of agreement.

### Scientific constraints unchanged

These changes harden computational reproducibility and resource use.
They do not independently validate true source membership, PSF wings,
spatially correlated noise or sky–galaxy degeneracy. Blank aperture
covariance is still a preliminary **same-coadd** diagnostic, and
direct valid-area aperture fluxes are not total galaxy fluxes.
`SCIENCE HOLD` remains required for Leo P and subsequent objects.


### Stage04 blank covariance: invalid NaNs outside the annular support

Further inspection found that the exact fractional ring integrator used
`sum(weights * image)` across the **full rectangular aperture-mask
bounding box**, after checking that no masked pixel had significant
positive fractional area. This nevertheless multiplies `0 * NaN`
for an invalid original FITS pixel *outside* the ring but *inside*
the rectangular box, which is NaN in IEEE-754 NumPy. A geometrically
valid blank therefore could be converted to a nonfinite sample and
cause the entire estimated matrix to be withheld. Because native
DESI cutouts may contain hard-invalid pixels, the old arithmetic
was not equivalent to summing on the exact accepted aperture support.

**Fixed:** a blank position is now rejected for **any** positive
fractional geometric overlap with an invalid, source-masked or
manually masked pixel. The actual flux sum uses ONLY positive exact
weights on validated native samples, never the surrounding rectangle.
No invalid native pixel is filled, zeroed as data, interpolated or
assigned fictitious variance. A synthetic 321x321 original frame
with deliberately inserted NaNs at zero-weight corner pixels checks
that a finite empirical annulus covariance remains measurable without
altering any native input data. This is a numerical and masking
correctness fix, not a calibration of correlated coadd noise.


### Same-UUID target-protection radius change invalidation

A further GUI lifecycle audit found that the same-object cache/source
identity already included ICRS RA/Dec, saved field width, source
quartet and data root, but omitted
`project.target_radius_arcsec`. The latter changes the native target
exclusion footprint and therefore the source-overlap decisions in
Stage02. Reopening a modified passport with the same UUID and
unchanged sky coordinates could previously leave already prepared
masks and related in-memory diagnostics associated with an obsolete
protection radius.

The Scene binding key now contains the saved target-protection radius
in addition to the prior dimensions. Changing it invalidates the
loaded Scene and dependent masks/isophotes/Background/Aperture, while
rebinding an unchanged passport remains idempotent. A dedicated
synthetic Qt regression tests a radius-only edit with identical UUID,
ICRS target, FITS catalogue and field width. This is conservative:
the original FITS may need to be reloaded even though their bytes have
not changed, but old mask decisions cannot silently survive.

This radius is a provisional source-screening VETO, not an adopted
outer isophote or proof of galaxy membership.


### Reverting a passport change must not resurrect a stale mask draft

The same project can be opened under multiple scientific source-key
identities (different native brick/field centre/target-protection
radius). The mask draft, however, belongs to the **object**, not the
previous source-key snapshot. Before this audit, a GUI cache entry
from the original radius could override a newer, saved mask draft
when the researcher changed the radius back. That silently replaced
user edits in memory with an earlier working revision, even though
the on-disk project held the more recent draft.

On a *same-root, same-UUID* source-identity change, the editor now
flushes its current autosave, discards obsolete cache entries for
that object, invalidates dependent native-scene products, and reads
the latest persisted draft. Cache restoration between genuinely
different object projects remains available. The Qt regression
changes the target radius, saves a newer exclusion edit, reverts the
radius and checks that the new edit survives rather than resurrecting
the original cached draft. This is a data-integrity fix, separate
from the source-mask geometry verification and from science adoption.


### Autosave conflict gate on source switching / project unbinding

A follow-up review found a conflict-resolution hazard in that cache
repair: `flush_autosave` deliberately uses optimistic HEAD revision
checking and preserves unsaved Scene edits on failure, but a subsequent
source switch or project unbind would still discard the only live
in-memory copy. Additionally, `bind_project` updated its target ICRS
before confirming that the previous draft had been saved, which could
make a rejected transition display an inconsistent target.

**Fixed:** all source-identity transitions now check whether a dirty
draft actually became saved after `flush_autosave`; an unresolved
concurrent revision conflict blocks the transition while keeping
the original loaded Scene, bound target, source key and edit draft.
The new target is not committed to GUI state until after this guard.
Project unbinding applies the same safeguard. The Qt regression
creates a conflicting externally updated `HEAD.json`, attempts both
a radius-only rebind and unbind, and checks that neither loses the
unsaved local edit nor overwrites the external saved edit.

This protects the in-session editor operation; it is **not** a claim
that the application or operating system can persist unsaved edits
after a forced process termination. User intervention is still
required to resolve an explicit optimistic-concurrency conflict.


### Stage 04.2 → 04.3: frozen-sky and pixel-evidence integrity (2026-11-12)

A conditional catalogue-mask ablation is valid only when the **same**
residual-sky scenario is used for both measurements. Previously Stage
04.3 recomputed its unmasked counterfactual with the default zero sky,
irrespective of the background scenario used for the supplied Stage
04.2 growth curve. With any nonzero constant or tilted-plane scenario,
the reported mask-only flux difference therefore also contained a
background change. The current GUI uses zero additional correction,
but the public scientific API allowed conditional scenarios.

Stage 04.3 now accepts an explicit frozen `ApertureSkyScenario` (defaulting
to zero for the unchanged GUI path), records it in the diagnostics and
uses it for both masked and unmasked measurements. All Stage 04.2
hypotheses must have consistent project/dataset/WCS/source-policy,
draft/catalogue fingerprints, target anchor, radial grid and sky
evidence. Their deterministic receipt digests are checked. The
outermost direct aperture is re-evaluated for every hypothesis against
the current native frame and masks; recipe digest, valid area, signed
flux and formal diagonal sigma must agree. This also catches
same-area mask substitutions or changed image/variance pixels that
the prior **area-only** check could miss. A changed background scenario
must be supplied explicitly; implicit changes fail closed.

Empirical blank covariance remains a diagnostic of the **uncorrected**
native coadd even if a conditional nonzero sky scenario is supplied;
a warning makes this scope explicit. These changes do not establish
science-grade sky, masked-flux recovery, calibrated covariance or
validated total photometry. The Leo P original FITS are still required
for end-to-end injection-recovery and foreground PSF validation.


### Native raster byte identity across Stage 04.2 → 04.3 (2026-11-12)

A deterministic curve receipt and agreement of the outer-aperture
flux/area/diagonal variance do **not** certify the native input arrays:
swapping two equal-weight pixels can preserve every such scalar while
changing IMAGE; compensating changes in inverse variance can likewise
preserve the formal error. Native WCS and source-catalogue hashes are
separate kinds of provenance.

`NativeBandFrame.native_input_sha256` is a lazy, per-immutable-frame
SHA-256 over the exact in-memory IMAGE, INVVAR, NEXP, MASKBITS and
effective hard-validity arrays, including dtype, shape, dataset identity
and a domain/version prefix. It hashes contiguous read-only buffers
without allocating a second full raster. Stage 04.2 now includes this
value in its deterministic receipt, and Stage 04.3 verifies it before
any comparison; Stage 04's GUI diagnostic also exposes the value.
Synthetic tests exchange source pixels and variance pixels while
preserving their previous outer-aperture integrals and verify that
the new byte-fingerprint safeguard refuses the stale input.

This is a **local computational identity**, not a claim of archived FITS
authenticity, correct WCS, scientifically valid sky/PSF, or calibrated
total photometry. FITS file hashes, WCS fingerprints and native-array
fingerprints have deliberately distinct scopes.


### Experimental PSF checkpoint restore: symlink-safe history (2026-11-12)

The PSF-0 checkpoint reader already excluded symlinked historical
entries in its restore loop, but it **sorted entries by Path.stat()**
before reaching that check. `stat()` follows a symlink; a broken
historical link raised an uncaught error and prevented the GUI from
restoring any otherwise legitimate PSF checkpoint. The sort now uses
`lstat()` (the directory entry itself), while the existing explicit
symlink rejection and warning remain unchanged. This changes only
safe access to saved diagnostic evidence, not the PSF, original FITS,
star selection, kernel or science-adoption policy. A regression
constructs a broken history link where supported by the test platform.


### Stage 03: optimizer convergence and complete bright-star veto (2026-11-12)

The robust Huber IRLS residual-sky fitter formerly returned its last
coefficients with `FITTED` even if the `max_iterations` cap was hit
without numerical convergence. The fitting contract now distinguishes
`NONCONVERGED` from `ILL_CONDITIONED`, keeps coefficients unavailable
and refuses downstream model evaluation/sky comparison until a
predeclared convergence criterion is actually met. No silent
relaxation of tolerances or post-hoc model selection occurs.

The automatic sky selector also previously took any dictionary-valued
`source_inventory` as a sufficient bright-Gaia veto source. In
practice SC1 explicitly bounds its source-level inventory and reports
`truncated`, `stored_count`, `total_relevant_sources` and
`screening_evidence_incomplete`. If one such list is truncated or
inconsistent, a bright source can be absent from the veto and a biased
sky still be fitted. Stage 03 now requires the existing source
inventory completeness/count receipts and returns
`SOURCE_INVENTORY_INCOMPLETE_FOR_AUTOMATIC_SKY` with **no fitted sky**
when they are missing/inconsistent. This is an input-integrity gate,
not a claim that a complete Tractor list has calibrated PSF wings or
proves galaxy-free sky pixels. Synthetic regression cases exercise
the IRLS iteration limit and each incomplete-inventory condition.


### Stage 03 sky optimizer / holdout setting types (2026-11-12)

`BackgroundSettings` previously compared settings numerically but
did not validate their scalar types: Boolean `max_iterations=True`
was interpreted as a one-iteration budget, Boolean holdout threshold
`True` became 1 nMgy/pixel, and fractional/Boolean pixel counts
could pass despite having no meaningful integral sample budget.
These choices can change whether a sky candidate is reported as
converged/holdout-qualified and must not silently map user input into
a different scientific algorithm.

Integer sample limits/iteration budgets now require actual
integer types (without Boolean coercion), and finite float-valued
controls reject strings and Booleans. Accepted NumPy scalars are
canonicalized into Python ints/floats for stable recording. The
fit/holdout algorithms and default settings are unchanged; synthetic
regressions cover malformed parameters and valid NumPy scalar inputs.
This does not turn an uncalibrated sky candidate into an adopted one.


### Stage 04.2 exact annulus: invalid overlap fractions fail closed (2026-11-12)

The cumulative-aperture integrator and blank-annulus covariance
estimator already enforced the exact `photutils` pixel-coverage
contract: finite weights in [0,1], with only ~1e-10 roundoff
allowed. The independently integrated growth-ring routine instead
used `np.clip` **before** checking for nonfinite values, silently
accepting arbitrary negative or >1 fractional overlaps. If such
weights arise from a dependency regression, the signed annular
flux, intensity and variance could be modified without an error,
and a plausible growth curve could be saved.

The ring integrator now checks finite/bounded exact pixel weights
**before** numerical-tolerance clipping, consistent with all other
native aperture implementations. A regression injects negative,
supra-unity and NaN overlap fractions through the `photutils`
interface and requires explicit rejection. This does not replace
correlated-noise covariance, actual isophotal fitting or scientific
sky/PSF validation.


### PSF core and holdout input consistency (2026-11-12)

The finite-native-core evaluator could previously report a measured
kernel FWHM while storing a NaN, infinite, nonpositive or Boolean
finder-FWHM prior. This allowed nonsensical comparative evidence
and potentially a non-JSON-serializable NaN in PSF QA. It now
returns INVALID_DETECTION_FWHM_PRIOR before measuring the core.

The held-out star evaluator for a *delivered* normalized nonnegative
PSF kernel previously allowed even-sized kernel arrays. This could
silently introduce a half-pixel reference-centre ambiguity despite
the upstream native PSF contract requiring odd dimensions of at
least 15 pixels. It now checks this same kernel-shape contract
before fitting holdout amplitudes, preserving identical radial
comparison mathematics for valid kernels. Neither correction
promotes experimental stellar PSFs to scientific adoption.


### Stage 04.3 numerical integrity of exact overlap fractions (2026-11-12)

The Stage 04.2 direct and annular flux integrators reject nonfinite
and nonphysical exact fractional pixel overlaps, but the Stage 04.3
residual-sky response calculation still used unconditional
`np.clip(weights, 0, 1)`. Corrupt negative or supra-unit weights
could therefore silently change the measured native-area response
to constant/tilted sky, despite a previously validated 04.2 curve.
Stage 04.3 now checks finite fractional weight bounds **before**
applying only roundoff-scale clipping, matching the integrators.
A regression injects bad photutils weights only after the independently
verified outer 04.2 aperture and requires a specific refusal before
reporting a sky derivative. No sky uncertainty or science covariance
is inferred by this numerical check.


### Stage 04.2 annulus input validation before off-frame return (2026-11-12)

The independent-ring integrator checked its native-grid source and
manual masks **only after** `get_overlap_slices`. When an annulus
fell wholly outside the original FITS, it returned an underconstrained
record without ever checking whether the supplied masks were Boolean
or had the correct source grid. This could conceal a stale Scene
state in an otherwise valid-looking growth QA, even though the
full-aperture integrator always validates masks on entry.

Ring measurements now validate both same-grid masks before
projecting/rasterizing and before any no-overlap return. The enclosing
growth run already snapshots editable masks once per curve: individual
rings reuse those validated arrays rather than making two costly
full-native-grid copies for every radius.
The positive inner radius, frozen growth hypothesis q/PA and radial
grid are also required to use genuine finite numeric scalars, not
Booleans or numeric-looking strings. Accepted NumPy scalars are
canonicalized to Python floats, consistent with the full-aperture
geometry contract. Synthetic tests check off-frame annuli,
mismatched and non-Boolean masks and malformed geometry controls.
No masked flux is reconstructed or declared scientifically adopted.


### Methodological correction: Stage 04.3 sky-gradient pivot (2026-11-13)

**Confirmed scientific formula defect.** Stage 04.1 computes conditional
sky subtraction using
`B(x,y)=b0+bx*(x-sky.reference_xy_pix[0])+by*(y-sky.reference_xy_pix[1])`.
Stage 04.3's published `dflux_dbx` and `dflux_dby` instead
computed `-sum(w*(x-aperture_center_x))` and similarly in y,
regardless of the frozen sky scenario's actual origin. This is a
different parameterization of a plane and is only equivalent when
the reference pixels coincide. For an asymmetric masked aperture,
nonzero slopes and a distant pivot, the old Stage04.3 derivatives
misrepresented the sensitivity to the very parameters of the
saved Stage04.2 sky scenario.

Stage04.3 now reports `reference_native_xy` equal to the **actual
frozen sky-scenario pivot**, computes `dflux_dbx_native_pixel3` and
`dflux_dby_native_pixel3` with that pivot, and retains separately
named aperture-centred reparameterizations in the QA. The direct
aperture operator, science HOLD and sky selection are unchanged.
A new scientific regression perturbs all three plane coefficients
around a deliberately off-centre pivot and compares the predicted
signed flux change with a fully re-integrated native exact aperture.
The previous test was insensitive because it deliberately aligned
its scenario reference with the reported target centre.


### Stage 04.3: native exact per-radius sky Jacobian (2026-11-13)

Following the confirmed sky-pivot correction, Stage 04.3 now
recomputes the **valid-area first moments** for *every* frozen
Stage 04.2 cumulative radius under the same WCS, masks and sky
reference, not only for the largest aperture. Each row records
`[-sum(w), -sum(w*(x-xref)), -sum(w*(y-yref))]`, conditional
on fixed geometry and validity. A second table records the exact
**disjoint ring** derivatives as differences between successive
cumulative rows, without treating the resulting noise as
independent. Every radial native valid-area measurement is
checked against its original frozen Stage 04.2 aperture area
before accepting any sensitivity result. Invalid exact fractional
weights are rejected before integration at every radius.

The new `radial_sky_response` payload has a dedicated versioned
schema, records parameter order and pivot, explicitly leaves sky
parameter covariance/uncertainty **absent**, and never adopts
a sky model. The original largest-radius `sky_response` fields
remain available and are now derived from the same native moment
operator. Synthetic regressions independently reintegrate each
aperture after changing all three sky parameters using an
off-centre pivot and asymmetric source/manual masks, compare
annular Jacobians with direct cumulative-flux differences, and
check entirely excluded/underconstrained native apertures.
This improves conditional diagnosis, **not** full noise covariance
or total-galaxy flux calibration.


### Stage 03: predeclared 2-D sky pixel sampling comparison (2026-11-13)

The automatic sky selector previously limited dense tiles to 384
valid native pixels using a regular `linspace` over the **row-major
index** of eligible positions. That is deterministic and brightness-
independent, but may alias with coadd/CCD sampling patterns and does
not explicitly balance a tile's two-dimensional spatial support.
This is an *observed algorithmic limitation*, not proof that the
existing Leo P fit is biased.

Stage 03 now offers an alternative, explicitly versioned
`SPATIAL_8X8_PROPORTIONAL_V2` diagnostic subsampling policy:
split each eligible tile into 8×8 native-pixel subcells, assign
the 384-sample budget proportionally to the number of eligible
pixels (largest-remainder allocation), and choose native indices
with a fixed per-tile NumPy seed independent of IMAGE and sky RMS.
The default `ROW_MAJOR_UNIFORM_V1` and its pixels remain **exactly
unchanged**. The new method is included as a sixth predeclared
sky-grid sensitivity experiment, without post-fitting model
selection or science adoption. The selected strategy, seed and
subcell geometry are recorded in the selection/experiment receipts.

Regressions check reproducibility, exact 384-pixel quotas, full
64-cell coverage on an unmasked dense tile, no fabricated pixels
on an irregularly masked tile, and unchanged legacy sampling.
A better real-galaxy sky result cannot be claimed without
independent blocked-out tests; the comparison remains diagnostic.


### Stage 04.3: radial source-mask ablation without reconstructed flux (2026-11-13)

The earlier source-mask counterfactual compared **only** the largest
aperture with every provisional catalogue exclusion disabled.
A single outer flux difference can hide whether the mask removes
central high-SB target knots, off-centre contaminants or diffuse
outer light, and it is not a membership classifier.

The current operator now measures, at each existing exact radial
footprint, the additional *original* valid pixels exposed if
only provisional source exclusion is disabled while all hard
and user-manual exclusions remain fixed. It records the newly
available area, **signed** native flux under the **same frozen sky
plane** and the corresponding conditional relaxed cumulative flux.
Adjacent-ring deltas follow from cumulative differences without
inventing independent errors. The last-radius result is checked
against a separate full native re-integration with the source
mask disabled; a disagreement fails closed.

This is a read-only diagnostic treatment contrast, not source
membership, PSF-model subtraction, masked-galaxy-light restoration
or a reason to add this light to a total. Regressions include
overlapping hard/manual masks, an off-centre tilted sky and
explicit negative source-mask flux. The per-source astrophysical
classification still requires catalog/PSF and injected truth.


### Stage 03: source-inventory row evidence is part of the bright-star sky veto (2026-11-13)

Automatic sky selection previously verified that SC1 inventory
`stored_count` matched the actual number of entries, but
subsequently assumed all rows were mappings and silently skipped
invalid `gaia_g_mag` or pixel coordinates. A malformed row could
either crash on `entry.get` or, more seriously, omit a potentially
bright foreground star from the halo veto while leaving the
inventory-complete flag **True**.

Stage 03 now checks each inventory entry's mapping type and
finite original pixel centre. Any invalid/malformed Gaia G
magnitude causes a countable incomplete-source status; a Gaia DR3
reference ID without a known G magnitude is likewise explicitly
unverifiable as a bright-star veto and therefore blocks automatic
sky candidate fitting. Unknown magnitudes for records *without*
a Gaia DR3 reference remain missing, rather than being fabricated.
The invalid-row and unknown-G counts are recorded and
`SOURCE_INVENTORY_INCOMPLETE_FOR_AUTOMATIC_SKY` is returned before
any fitted sky is reported.

This cannot prove that all uncatalogued bright stars have been
identified; the veto radii and true wing extents remain
uncalibrated. Synthetic regressions cover nonmapping entries,
invalid/missing positions, Boolean/string G magnitudes and Gaia
identities without usable brightness evidence.


### Stage 04.3: use the same conditional sky plane for blank-ring controls (2026-11-13)

The earlier blank covariance measured signed flux on the **raw**
native coadd even when the frozen Stage 04.2 aperture had already
subtracted a nonzero conditional constant/tilted plane. Thus the
reported same-field covariance referred to a *different flux
operator* and could include deterministic fitted sky gradients that
the target measurement had removed.

The versioned `blank-annulus-covariance/2` diagnostic accepts the
same immutable `ApertureSkyScenario` as Stage04.2 and integrates
`sum(w*(I-B))` on the original valid, unmasked blank pixels
without modifying the image. Full-scale and reduced-scale blank
probes inherit that **same scenario**, while the default zero
additional correction remains equivalent to the original raw
coadd diagnostic. The QA records the full sky scenario and
explicitly states that neither it nor the same-coadd covariance
has become scientifically calibrated/adopted.

A controlled synthetic tilted-plane field tests the reduction in
blank covariance when the correct known plane is subtracted.
An independent large-aperture fallback test verifies that all
smaller-scale blank placements use the same frozen sky
parameters. This corrects a conditional-method mismatch; it is
**not** a sky posterior, off-diagonal CCD covariance reconstruction,
or proof that the plane should have been subtracted in science.


### Stage 03: explicit zero-inclusive native sky-model RMS comparison (2026-11-13)

The former `provisional_automatic_model` selected a **numerical**
winner between the constant and tilted-plane fits only, and even
the original QA noted that zero additional sky correction was
absent from that heuristic. The coadd has **already** been sky
subtracted upstream; a two-fit "winner" can have higher error
on exactly the same heldout pixels than doing nothing.

Stage 03 now additionally records a versioned *three-way descriptive*
heldout RMS comparison for the zero-correction baseline, constant
and tilted-plane models, all evaluated on identical fixed native
samples. A tie prefers the simpler/no-subtraction label only
for display. Existing two-fit experiment naming and baseline
model selection are preserved for longitudinal comparisons,
but all new fields explicitly say that the numerical minimum
does **not** choose physical sky, quantify missing galaxy light
or authorize subtraction. A synthetic case with biased training
pixels but pristine held-out pixels confirms that the zero
reference can correctly rank first.

This is a research QA transparency improvement, not a claimed
scientific sky adoption or validated survey background.


### Stage 04.3: independently verify *all* frozen radial fluxes and formal variances (2026-11-13)

The stage-04.3 receipt verifier originally recomputed native flux
and diagonal formal sigma only for the **outermost** Stage04.2
cumulative aperture. The new radial Jacobian introduced a second
exact native fractional overlap for each *inner* radius, but
checked only its covered **area**. An internally re-receipted stale
inner cumulative flux or formal sigma could enter the radial
mask ablation or other diagnostics while all saved outermost
values and WCS/mask geometry remained unchanged.

The native radial loop now integrates `sum(w*(image-sky))` and
`sqrt(sum(w²/INVVAR))` for every frozen radius from its
already available **original native** valid-pixel weights.
Both results are compared against each saved Stage04.2 record,
independently of the saved whole-curve SHA. Underconstrained bins
must remain null rather than fabricated as zero. The check
reuses radial weights, does not reread FITS or multiply the
number of full-size mask snapshots, and does not modify the
science-state or sky/PSF/foreground assumptions.

A synthetic regression deliberately changes an *internal*
cumulative flux or sigma, recomputes the outer growth-receipt
SHA and keeps the final aperture untouched. Both tampered
inputs are now rejected; an entirely masked curve still
passes as underconstrained. This fixes a numerical QA
consistency failure, **not** missing galaxy light or
off-diagonal noise calibration.


### Stage 04.2: enforce exact annulus/CoG flux closure before reporting a measured curve (2026-11-13)

Stage04.2 already computed the flux difference between two
independent numerical methods: a directly integrated native
elliptical annulus and the difference of its enclosing cumulative
apertures. Until this change, it merely **recorded** an arbitrary
nonzero residual, yet labelled the growth curve diagnostically
measured even if the two operators could disagree catastrophically.
This is a confirmed numerical consistency gap, not a model
selection question.

The growth operator now refuses a comparable annulus if
`abs(F_outer-F_inner-F_independent_ring)` exceeds a fixed
roundoff allowance `1e-6 + 1e-9*max(1,abs(F_outer),abs(F_inner),abs(F_ring))`
in native nMgy, with explicit radial index. A fully excluded
ring remains **not measured** (null), but when both cumulative
endpoints exist, their flux difference is still checked against
zero usable *increment*; no synthetic annular flux is assigned.
The result includes the number of closures actually checked
and explicitly disclaims any interpretation as calibrated
sky/PSF/covariance or astrophysical recovery.

A synthetic regression alters an otherwise valid **independent
ring** by 0.01 nMgy while leaving the two direct cumulative
measurements unchanged; Stage04.2 must now fail closed. Valid
masked/underconstrained rings continue to be distinguished
from observed zero flux. This does not impose a monotonic
CoG or exclude real negative annular flux.


### Stage 03: paired noise-weighted versus pixel-weighted residual-sky controls (2026-11-13)

Current Stage03 robust fits use diagonal original native IVAR
in their regression weights, but the frozen held-out residual
RMS candidate ranking intentionally assigns **equal weight**
to every geometrically selected pixel. Heteroscedastic coadd
noise can then move a numerical ranking simply by changing
the relative fraction of noisy and quiet regions. Neither
objective alone represents a calibrated physical sky likelihood:
coadd correlations, systematic outer-galaxy light and PSF
wings are still unresolved.

Stage03 now reports a paired, **descriptive**
`sqrt(sum(ivar*(I-B)^2)/sum(ivar))` for zero additional
correction, constant and plane, over **exactly the same
original native holdout pixels** as the unweighted RMS.
It publishes whether the two rankings disagree and records
explicitly that the weighted statistic is not a chi-square,
correlated-pixel covariance or science-selection criterion.
The existing fit, historical candidate heuristic, masks,
pixel sample and all sky-science holds remain unchanged.
A synthetic heteroscedastic-INVVAR regression reconstructs
each held-out model/residual on original native coordinates
and independently checks the weighted RMS.


### Stage 04.3: blank covariance must sample the same subpixel aperture phase (2026-11-13)

The previous blank-ring controller used whole-pixel integer lattice
centres irrespective of the actual WCS-projected galaxy centre's
*fractional native pixel phase*. Although ellipse axes and PA were
held fixed, an aperture shifted fractionally has different boundary
overlap coefficients. For diagonal original pixel variances, even
a perfectly uniform-noise field has a distinct
`sum(w_p² sigma_p²)` for different subpixel phases. Thus the
background test and the science aperture used different linear
operators, a confirmed methodological inconsistency in a
same-shape blank covariance comparison.

`estimate_blank_annulus_covariance` now translates the original
deterministic, quadrant-stratified native candidate lattice by
the **exact fractional part** of the target's projected centre,
not by fractional resampling of the image. Every candidate
centre is a whole-native-pixel translation of the target
centre. Original masks, native pixel intensities, full-target
exclusion and pairwise blank separations are unchanged; the
two-dimensional native phase, unchanged selection properties
and non-adoption are recorded in schema
`desidera.photometry.blank-annulus-covariance/3`.

A synthetic fractional-WCS test verifies each accepted blank's
native-pixel phase and compares `sum(w²)` in each exact
translated ellipse to the **actual target** aperture. This
removes one known mismatch but does not correct spatial
variability of pixel IVAR, field-scale correlated noise,
sky/foreground contamination, or galaxy light missing from
survey coadds. The covariance remains diagnostic only.


### Stage 04.3: do not mistake blank-field variance for target-pixel IVAR (2026-11-13)

Matching all blank placements to the target's fractional native
pixel phase removes a *geometric* disparity, but the coadd
INVVAR/NEXP field can vary by orders of magnitude across the
galaxy and blank sky regions. The empirical sample variance of
noise from high-coverage blanks cannot be treated as the
galaxy's noise variance in a lower-coverage target aperture.
The same issue arises when target pixels are excluded by a
manual or provisional source mask: blanks still use complete
rings while the measured target flux uses only its valid area.

The blank covariance QA now separately integrates the expected
**original native diagonal-only** formal `sum(w²/INVVAR)` in
each target ring *after its original source/manual/hard masks*,
and for every unmasked complete blank ring. It reports the
target variance, the distribution of blank formal variances
and their target-to-median ratios by radial index. Missing
target ring flux remains null, and no invented zero or
uncertainty correction is made. The formal ratios are explicitly
**not** full-coadd correlation corrections, covariance
normalization factors, sky-error likelihoods or evidence of
galaxy completeness. The existing empirical covariance is
left entirely unscaled.

Synthetic regression uses a known strongly heterogeneous IVAR
field (noise-heavy around the target, quiet in blank fields)
and independently integrates the smallest-ring exact native
weights; a second test checks a fully masked target inner
ring stays null despite valid off-target control rings.
