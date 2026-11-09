# P1.2a — Scene & Geometry: sky-coordinate mask editing and initial ellipse seeds

> **P1.3 clarification:** the provisional Scene_v0/user-mask draft is
> permissible before PSF reconstruction, but final PSF-aware Scene_v1
> masks follow a validated per-band PSF and preliminary star-cutout sky.
> See `PSF_SCENE_SKY_ORDER_P13.md`. An initial ellipse seed is not a
> scientifically accepted background or isophote fit.

Date: 2026-11-03. Parent: P1.1 Photometry Workspace (draft PR #39).
**Status:** first scientific primitives + synthetic tests. **Not yet** an
interactive image editor, calibrated source segmentation, final sky model,
iterated ellipse fitter or adopted photometric science.

## 1. Researcher UX: one edit, all compatible bands

### Screen flow

Galaxy Projects -> Galaxy Workspace -> Photometry -> Scene & Geometry.

Proposed layout at 1920×1080:
- Left: task navigation and a searchable ordered list of semantic mask
  layers, with counts, provenance, visibility, lock/disable and scope.
- Center (largest area): local FITS canvas with zoom/pan, ICRS cursor
  readout, comparison between g/r/i/z and WCS-registered vector overlays.
- Right: selected object/region inspector, geometry radius/shape,
  application scope, mask/photometry purpose, edit-history controls,
  initial ellipse controls and preflight blockers.
- Bottom: state of each filter, affected pixel/sky area, validation,
  stale dependencies and explicit save/apply audit log.

For <=1366×768, collapse left workflow into a narrow nav and place
inspector in a scrollable panel/drawer. Do not let masks occupy the
entire fit viewer or hide primary task actions.

### Direct manipulation contract

1. The researcher sees **the original local FITS**, never a remote
   Legacy Viewer preview as the pixel source for editing.
2. Tool modes: select/inspect; circle (single click + radius or
   center/radius drag); ellipse; polygon (click vertices, double-click
   to finish, ESC to cancel); freehand brush with **arcsec** radius
   sampled into sky-coordinate strokes; AI/segmentation candidates
   are disabled until explicitly generated, classified and reviewed.
3. Hotkeys (to be implemented and tested): V select, B brush,
   P polygon, C circle, E ellipse, Space held pan,
   Enter commit polygon, Esc cancel unfinished action,
   Ctrl+Z/Shift+Ctrl+Z undo/redo **draft edits**, Delete removes
   selected edit. Edit deletion is different from science-pixel restoration.
4. Each edit has purpose EXCLUDE, RESTORE_AUTOMATIC, PROTECT_TARGET,
   SKY_TRAIN or SKY_HOLDOUT. Explicit RESTORE_AUTOMATIC affects only
   editable automatic contamination masks, never original hard rejection.
5. Default *scope* = **All compatible bands**. The user may select
   one named band or turn off the effect of one shared region in a
   specific band by creating an explicit, logged band-specific policy.
6. Region list supports enable/disable, rename, purpose, band-scope,
   provenance, ordering/grouping, isolate/highlight and delete.
   The scientific application order cannot be changed by a visually
   convenient list reorder.
7. Save is **explicit**. Preview/undo does not touch FITS or adopted
   runs. Each save creates an immutable version; old versions remain
   available. Staleness is propagated to derived sky/geometry/photometry.
8. Show g/r/z side-by-side and an RGB view where available, with linked
   sky cursor/zoom. A shared vector region overlays each native WCS
   after fresh reprojection; overlays never reuse a different band's
   pixel mask or claim that pixel grids coincide.
9. Report per-band warnings if a region exits the footprint, crosses
   masked holes, has less than the requested useful fraction, or has
   variable PSF wings so a fixed sky radius is insufficient. Do not
   fabricate measurements for unavailable i-band data.

## 2. Layer semantics: a scientifically meaningful mask stack

Native and instrument provenance remains **immutable**:

- HARD: nonfinite IMAGE, bad IVAR, absent coverage, instrument defects,
  native saturation policy, missing/inconsistent WCS. No manual tool
  can restore information missing from the measurement.
- AUTOMATIC (editable): classified star cores/wings, neighboring
  galaxies and candidate artifacts. Store original proposals,
  segmentation provenance and confidence separately.
- MANUAL EXCLUDE: overrides automatic restoration; linked to region ID
  and researcher action.
- RESTORE AUTOMATIC: only revokes an *editable automatic* exclusion
  within a shape; hard data remain hard. To reverse a manual exclude,
  remove/disable that exclusion itself.
- PROTECT TARGET: defines intended galaxy footprint for sky fitting
  and source separation; it never cancels saturation or a real
  foreground star on top of the galaxy. The galaxy's own star-forming
  knots should normally be retained as target light, not auto-rejected.
- SKY_TRAIN / SKY_HOLDOUT: geometrically disjoint regions for background
  fitting and independent testing; contamination and target protection
  are removed from eligible sky sample pixels.

For an individual band b and pixel p:
```text
bad_b(p) = hard_b(p)
         OR manual_exclude_b(p)
         OR [automatic_contamination_b(p) AND NOT restore_auto_b(p)]
```
All operations must separately report an exclusion-reason census.
A mask removes observations from the estimator; it **does not set their
fluxes to zero, invent a replacement galaxy value or silently inpaint**.

When deriving radial profiles, account for the fractional masked
aperture, exact/controlled fractional pixel coverage and missing
azimuthal sectors. Severe missing coverage must trigger QA rather
than an unqualified area-scaled extrapolation.

### Cross-band mask sharing

Persist **world-space vector geometry** (ICRS center/vertices in
degrees, angular sizes in arcsec, PA measured North through East).
Rasterize it afresh for each source WCS. This avoids requiring users
to draw the same ellipse/polygon for g, r, z.

Automatic masks need **band-dependent** generation/dilation:
bright-star wings, PSF size, saturation and the optical appearance
of knots can vary strongly with filter. A region shared in sky
coordinates should not force identical pixel counts or blanket
masking of differently detected sources.

The initial implementation uses circles, sky-ellipses and sky polygons
with small-field projection (<0.5 degree extent); pixel center tests
or explicitly selected supersampled coverage fractions. Curved
boundaries, high distortion, wide-field polygons and precise fractional
area uncertainty require further validation. These are independent
of the older band/grid-specific `SceneState`/boolean `MaskLayer`;
adapt a draft only after WCS validation, never by copying its raster.

### Persistent revisions

- A galaxy project UUID owns `analyses/mask-drafts/`.
- `draft-<sha256>.json` stores immutable vector regions, scopes and
  purposes; `HEAD.json` references the current **draft**.
- Saving uses an atomic HEAD swap and expected-digest validation to
  reject concurrent overwritten edits.
- Drafts never imply scientific acceptance. A future immutable
  photometry run must freeze a specific draft digest and the original
  FITS/WCS fingerprints; later edits create a different analysis branch.
- A malformed/cross-project draft is an error and may not silently
  reset to empty.

## 3. Initial photometric ellipse: AUTO by default, editable always

Do not use one arbitrary threshold image to define a final galaxy
isophote. For LSB sources, the center, ellipticity, PA and extent are
coupled to background estimation, contamination and PSF.

Proposed full chain:

1. Validate FITS/WCS, pick an explicit detection band or a joint
   matched/PSF-aware detection image. A composite detection image is a
   **derived diagnostic**, never a replacement for original fluxes.
2. Use literature/project ICRS position as a **prior starting center**.
   Segment connected diffuse signal at several independently calibrated
   spatial scales/thresholds, excluding hard/contaminant masks and
   protecting legitimate galaxy knots.
3. Determine initial shape from weighted second moments of a connected
   region near the target, not from the brightest source elsewhere.
   On irregular systems multiple physical components may be blended:
   flag ambiguous components and excessive centroid shifts.
4. On each plausible intermediate isophote, refine center, PA and q by
   robust harmonic analysis (Jedrzejewski / `photutils.isophote.Ellipse`);
   assess angular completeness, radial gradient, masks, fits and stop
   codes. The library requires a sensible starting ellipse; its automatic
   center finder alone does not discover the target.
5. Combine stable inner/intermediate geometry with sky/PSF sensitivity.
   Evaluate multiple thresholds, sky models and bands; do not force the
   geometry of small star-forming knots upon the diffuse stellar body.
6. Present candidate ellipse and suggested radial validity interval.
   User may drag center, major/minor axes and orientation, or lock
   parameters individually; manual geometry also retains ICRS/scope
   provenance and its reason.
7. Use fixed geometry to establish robust cross-band colors; later allow
   free center/PA/q as a function of radius if the data support it.

**Implemented in this increment:** `propose_isophotal_seed` uses an
upstream background-subtracted image and provided RMS, normalized
masked Gaussian smoothing, a nearby connected component, weighted
second moments, and an explicit not-adopted state. It does NOT do
multi-threshold consensus, sky estimation, significance calibration
after smoothing, iterative photutils isophotes or final aperture sums.
`seed_as_sky_geometry` converts pixel ellipses through the local WCS
tangent mapping; it is not a wide-field nonlinear distortion solver.

A rejected/non-detection case should offer researcher manual geometry
without silently replacing it with the image center or another star.

## 4. Scientific tests required before acceptance

- Masks: rotated/mirrored WCS, RA crossing 0, high/low Dec, unequal
  pixel scales, overlapping bricks, missing filters, region outside
  footprint, clipping and partial pixels.
- Safety: negative pixels remain untouched; restoring automatic masks
  never changes native hard exclusions; color-specific masks remain
  isolated; target protection never reinstates bad pixels.
- Editing: click/drag correctness at zoom factors and HiDPI, reset pan,
  hit-testing, accidental edits, keyboard cancel, undo/redo, reload and
  concurrent draft saves.
- Automatic candidates: injection of masked stars/galaxies at varying
  depth and PSF; benchmark over- and under-masking, including clumpy
  dwarf galaxies and nuclear stars.
- Ellipse: synthetic Sérsic and exponential images, varying PA/q,
  high gradients, low S/N, star-forming asymmetries, off-center
  catalog positions, bright neighbors, and poorly constrained sky.
- Propagation: changing masks must invalidate dependent sky/isophotes/
  photometry, not cached survey discovery or unchanged raw files.
- Comparison: AutoProf and photutils on matched synthetic and real
  fields, with measured bias and interval coverage.

No scientific claim of `VALID` until real-FITS tests, manual visual
review, sensitivity to sky/PSF and covariance calibration pass.

## 5. Relevant documentation

- photutils isophote: https://docs.astropy.org/projects/photutils/en/stable/api/photutils.isophote.Ellipse.html
- photutils segmentation: https://docs.astropy.org/projects/photutils/en/stable/api/photutils.segmentation.detect_sources.html
- Astropy regions: https://docs.astropy.org/projects/regions/en/stable/api/regions.EllipsePixelRegion.html
- Native source contract: `docs/architecture/PHOTOMETRY_INPUT_CONTRACT_P05.md`
- Existing scene mask semantics: `docs/SCENE_MASKS_V03.md`
- Photometry architecture: `docs/architecture/PHOTOMETRY_CORE_P11.md`

## 6. Next executable slice

P1.2b: integrate a Qt manual mask canvas driven by **loaded and
source-verified local FITS**, with world-space circle/polygon/brush
tools, temporary preview, per-band linked overlays, explicit revisions
and replay tests. Do not make the sketching tool active for a stale
FITS cutout or a remotely displayed preview.

P1.3: implement residual-sky candidates, holdout regions and a
threshold-consensus/iterative isophote initializer, then compare
against photutils and AutoProf before any photometric adoption.


## 2026-11-09 · SC1 preliminary isophote suggestion in GUI

The object-bound **Scene & Masks** editor now offers a bounded
`Find connected isophote` diagnostic for native original Legacy IMAGE and
its independent native hard-validity mask. A project ICRS position anchors
the component, and a researcher may alternatively click a diffuse-light
centre. Settings are **local smoothed-background contrast** (not a
calibrated sigma level) and Gaussian smoothing in native pixels.

- A clipped local background location and empirical smoothed-background
  dispersion are measured in a distant annulus, independent of any later
  adopted sky model. The original IMAGE is never overwritten.
- A 4-connected threshold component close to the selected centre is
  identified inside a bounded search aperture. Failure/no detection and
  aperture-edge truncation are reported rather than inventing an object.
- The actual connected-component **outer polygon** is traced on native
  pixel-cell edges, softly simplified to at most 192 sky vertices and
  transformed by original WCS into ICRS. Its two-pixel display margin,
  local-sky estimate, threshold and smoothing are diagnostic metadata.
- An independent **moment ellipse** is available from the same connected
  component, with bright-peak weights clipped to reduce knot domination,
  centre anchored to the selected galaxy coordinate, and local WCS
  tangent mapping into sky PA/axis ratio.
- Both alternatives are **preview-only** until the user clicks
  `Add proposed region`, which creates a revisioned, non-adopted
  `protect_target` ICRS mask; Undo and the existing immutable draft
  revision/optimistic concurrency mechanism remain active. `Discard`
  writes nothing. The same sky geometry can be reprojected in another
  band; any actual physical differences between the bands still require
  review.
- Native FITS cutout loading, catalogue queries, adopted PSF/sky
  dependencies and saved science calculations are not triggered by
  suggestion. Display stretch may be switched between faint/standard
  asinh and faint linear without changing the scientific arrays.

**Science limitations:** no validated spatial sky model, PSF wing
subtraction, external Gaia/Tractor contaminant removal, correlated-noise
false-positive calibration, multi-scale consensus, isophote fitting
(Jedrzejewski), partial-pixel segmentation uncertainty or final photometry.
Local annulus clipping may be biased by neighboring sources or galaxy
halos. Even if a contour appears visually convincing on Leo P, it is
an exploratory user-reviewed **mask boundary**, not a measured or
science-adopted galaxy isophote. Do not quote its `a`, `q`, `PA` or
footprint as calibrated morphology without independent validation.


## 2026-11-09 · Correction: source screening PRECEDES isophotal suggestion

The exploratory `Find connected isophote` operation now depends on an
explicitly prepared **catalogue-assisted per-band neighbour-mask preview**,
rather than running unmasked on the original image. The preparation stage
`1 · Prepare neighbour masks`:

- Loads the already-local, exact-brick DR11 South Tractor table via the
  existing validated reader (same `release/region/brick/band`, full
  local catalogue SHA-256); it neither downloads nor fabricates catalogue
  entries. Other releases and missing catalogues fail CLOSED.
- Retains immutable native hard invalidity from IMAGE/INVVAR/NEXP/MASKBITS.
  BRIGHT, MEDIUM, GALAXY and BAILOUT are counted as **diagnostic cautions**,
  not blindly turned into blanket hard exclusions. In particular GALAXY
  may contain the target itself.
- Creates conservative **circular screening envelopes** around catalogued
  PSF and extended neighbours with usable band flux, plus explicitly drawn
  manual `exclude` areas. Envelopes approximate a locality to omit; they
  are NOT rendered Tractor profiles, PSF-convolved subtraction, or an
  optimal wing mask.
- Protects a central ICRS circle of the project's configured target radius
  plus existing manual `protect_target` shapes. Any automatic source
  envelope *touching* protection is skipped entirely and counted as an
  ambiguity for manual review. Manually drawn `exclude` may override
  protection; `restore_auto` can revoke only editable auto envelopes.
- Shows the proposed native-pixel contamination mask as a separate
  colored overlay, with masked, skipped-near-target and invalid-pixel
  counts. No FITS is modified, no saved scene draft or PSF is adopted.

Only after this stage does `2 · Find masked isophote` become available.
The connected-component detector and provisional sky estimation operate
on the native IMAGE with hard and proposed source envelopes excluded.
Mask/source/band/anchor/draft identity is checked, and editing a region
or changing the band/anchor invalidates the preparation. The eventual
isophote candidate is still only a cyan preview until `Add proposed region`.

The prepared mask is **pre-PSF and not scientifically accepted**. In
particular the catalogued PSF and galaxy-model estimates have not been
converted into per-source pixel model images, and bright-star wing
leakage and blends overlapping target protection remain unresolved.
A later full Tractor-model/empirical-PSF-aware source model (not the
all-source coadd model image blindly subtracted from data) is required
before calibrated outer-isophote extraction or scientific adoption.
The manual drafts remain useful for resolving explicitly flagged blends.


## 2026-11-09 · SC1 morphology and stability acceptance corrections

Real Leo P `i`-band masking revealed a 47,407-pixel connected area,
`q≈0.60`, and a large ≈59.4-native-pixel centroid displacement
from the selected object ICRS. The northern lobe cannot be promoted to
a measured isophote on one threshold alone.

**Source-type-aware contamination:** where DR11 Tractor supplies
`SHAPE_R`, `SHAPE_E1/E2` for an extended `EXP/DEV/SER` source, a
source-specific *elliptical screening envelope* is oriented using the
catalogue ellipticity convention and the native celestial WCS pixel
Jacobian. See [DR11 Tractor catalogue shapes](https://www.legacysurvey.org/dr11/catalogs/):
`epsilon=(a-b)/(a+b)*exp(2i*phi)`, and the astronomical
north-through-east `PA=(180-phi)%180`. `PSF`, `REX` and objects
without usable shape data remain circular or explicitly use a round
fallback. The radial extent remains deliberately heuristic, NOT an
actual Tractor-rendered pixel model or a calibrated PSF-wing mask.
Source envelopes intersecting the protected galaxy are not masked;
their catalogue positions receive amber **review** rings. These
unresolved neighbours may still contaminate the proposed boundary.

**Object-bound contour stability:** the connected-source detector
can optionally apply a 2-native-pixel morphological opening to remove
thin bridges (toggle OFF to inspect possible real filaments). A separate
diagnostic repeats the local search at neighboring thresholds and
smoothing values, reporting four trials, area ratio, axis-ratio spread
and offset of the light centroid from the selected ICRS anchor.
The actual *native-pixel component intersection-over-union (IoU)*
against the main threshold is also required: matching area and q
alone cannot hide a spatially displaced lobe. A minimum overlap
below 0.50 is explicitly flagged. Failed trials remain failures
in QA, not cherry-picked away.
A strong offset or unstable contour receives
`LOW_CONFIDENCE_REQUIRES_EXPLICIT_OVERRIDE`, is drawn orange, and
cannot be added even as a DRAFT unless the user checks the explicit
manual-review override. This is a **risk flag**, not a p-value or
scientific measurement of galaxy morphology.

**GUI review:** the image remains read-only; the researcher can
toggle the purple neighbour-screening layer, adjust its opacity,
toggle the unresolved-identity amber rings, inspect the concise QA
summary or expand the evidence/trial details, and compare thin-bridge
opening ON/OFF. An automatically derived isophote is **band-local**
by default; `Apply this band-derived shape to all bands` is separate,
off by default and requires explicit action. Ordinary manual ICRS
regions retain their existing all-band capability. All accepted
regions remain versioned, reversible, unadopted drafts.

**Outstanding scientific gates:** an accurate 2D individual-source
Tractor rendering with observational PSF / coadd weighting;
position-dependent empirical-PSF wings; galaxy-neighbour deblending;
residual sky spatial structure and independent sky holdouts;
cross-band detection consistency with different morphology;
robust physical surface-brightness isophotal fitting and injected
recovery tests. No test in this stage establishes those claims.
