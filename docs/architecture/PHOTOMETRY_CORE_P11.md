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
