# P1.2a — Scene & Geometry: sky-coordinate mask editing and initial ellipse seeds

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
