# P1.2c: Smart Mask — researcher-seeded source-aware contamination masking

> **P1.3 clarification:** this connected-threshold candidate is
> **PROVISIONAL_PRE_PSF**, not the final PSF-aware contaminant mask.
> PSF-informed wing/leakage decisions require a separately validated
> PSF built from selected stars with independently estimated
> star-cutout backgrounds. See `PSF_SCENE_SKY_ORDER_P13.md`.

Status: **experimental numerical proposal only**. Follows stacked PR #40
(mask/geometry drafts) and PR #39 (Photometry Workspace). There is no
image-click UI, full sky estimation, empirical PSF fitting, or automatic
acceptance in this change.

## The requested workflow

1. Scientist reviews local calibrated/native FITS, opens **Smart Mask**.
2. Click a foreground star or draw an approximate circle; UI freezes the
   **ICRS anchor**, source ID, and angular radius, not the pixel coordinates.
3. On an independently processed optical band, a preview candidate is
   computed from source-bearing pixel intensity after a supplied residual
   sky and RMS map, within the user-defined bounded search area.
4. UI (future) shows outlines for three alternative detection thresholds
   and optional estimated wing margins. User adjusts boundary and compares
   original g/r/i/z (when available). No hidden cross-band copying.
5. Per-band candidates can be accepted one by one or reviewed as a group.
   Acceptance *only* creates new non-destructive, versioned scene layers;
   it never changes original FITS, MASKBITS or prior analysis runs.
6. Edited mask version gets a reproducible dependency ID. Any subsequent
   background, geometry, photometry or SBF execution freezes exact
   mask and sky revisions and refuses to silently reuse stale artifacts.

An automatic system must NEVER search for "the best mask" solely by
optimizing remaining image smoothness: overmasking real LSB galaxy light
could appear numerically preferable.

## Tiered source-specific strategies

### Tier A: Connected isophote around a user anchor (pilot implemented)

Works best for **moderately bright, approximately isolated stars**:
- Use original-band WCS to find user ICRS anchor.
- Restrict search to a small sky circle, so the seed cannot jump to a
  bright but distant galactic nucleus or another star.
- Normalize Gaussian filtering by valid-pixel support, without writing
  to input pixels or filling invalid originals.
- Use a residual-sky-subtracted contrast divided by supplied local RMS;
  connect contiguous thresholded pixels to peak near the user click.
- Refuse a component occupying most of the search circle.
- Offer user-selected morphological dilation as a preliminary faint-wing
  margin; always mark extended PSF wing treatment as *unvalidated*.
- Flag component reaching user search limit, native missing pixels,
  possible asymmetric blends, multiple prominent peaks and overlap
  with protected galaxy light.
- Return immutable, per-band candidate raster and explicit flags.
  `to_reviewed_scene_layer(confirmed=True)` requires researcher action.

It is **NOT** a rigorous per-pixel-significance map: source coadds
contain correlated noise; Gaussian smoothing changes noise variance.
Thresholds are heuristics until local noise correlation is calibrated.
Candidate area and dilation sensitivity must enter the systematic QA.

### Tier B: PSF-informed masking (planned; not implemented)

For bright/saturated stars or significant extended wings:
- Identify a usable star PSF in relevant passband and spatial location,
  including empirical bright-star wing profile where possible.
- Use `I_star(r) = F_star P_b(r) + B_b` only where that model is
  justified; require evidence of clipping, saturation and detector
  artifacts. Never fit saturated cores with unsaturated-Poisson
  assumptions.
- Choose contour by predicted leakage into galaxy flux/sky estimates,
  *not solely* by FWHM, source intensity or a fixed n-sigma radius.
- Optionally offer "Model & subtract" as a **separate explicitly
  validated operation**. A subtraction residual still enters the
  photometry/sky uncertainty budget.
- Keep diffraction rays, ghosts, bleed trails and other anisotropic
  structures as separate masks/models, where supported.

### Tier C: Blended or galaxy-overlapping contaminants (planned)

When star and galaxy light overlap, ordinary segmentation cannot know
which pixels belong to each physical source. Offer:
- Manual bound/protected target region.
- Simultaneous galaxy + star PSF + sky fit, then inspect residuals.
- A comparison of mask-only vs PSF-subtraction estimators.
- Explicit flags if source flux and diffuse target light remain
  degenerate. Never silently make the galaxy symmetric by masking a
  real asymmetric feature.

## Shared sky anchor, band-specific mask: precise semantics

A star at a fixed ICRS position has the same celestial anchor in
all bands, but different signal-to-noise, PSF, color-dependent wings,
saturation, defects and local contaminants. Therefore **one click
means one source identity, NOT an identical pixel mask**:

```text
source_id + RA,Dec + user search region
    -> WCS_g -> detect / propose M_g
    -> WCS_r -> detect / propose M_r
    -> WCS_z -> detect / propose M_z
    -> optional WCS_i only where a real dataset exists
```

If no source is detected in a band, report NO_DETECTION and offer a
reviewable transferred *sky geometry* from the reference band; do not
fabricate a detection or a PSF. Choose a conservative cross-band union
only with an explicit policy and report changed valid area per band.
Recompute even when FITS array shapes happen to match.

## Quality and experimental validation

For a candidate C, inspect at least:
- Fraction of source flux predicted to leak outside C under a **tested**
  PSF (undefined before PSF modeling).
- Coverage of independently labelled contaminant truth masks on
  injected artificial stars; false exclusion of galaxy light.
- Stability as threshold, sky, PSF, dilation and band change.
- Fraction and radial location of galaxy target-protection overlap.
- Effects on recovered galaxy surface-brightness profile, total
  magnitude, half-light radius and structural fit, including
  corrupted-mask stress tests.
- Precision of ICRS-to-native-pixel mapping at field edge, RA=0 and
  rotated/mirrored WCS.

No trained confidence/probability or "best mask" score is assigned
without calibration against astrophysical ground truth.

## Source-specific background/skymask iteration

Background estimation requires masking sources, while a reliable
star mask requires a reasonable background. Use a controlled loop:
preliminary background -> provisional seeded masks ->
background refit with holdout -> review candidate stability.
Stop/reject on strong changes, not just optimizer convergence.
The scientific background module P1.3 is still required.

## Software contract in this PR

`src/lsbpipe/photometry/smart_masks.py` provides:
- `StarMaskSeed` — user-requested source identity, ICRS anchor and search circle.
- `propose_seeded_star_mask(NativeBandFrame, StarMaskSeed, ...)` —
  local native-band segmentation candidate requiring supplied residual
  sky, local noise RMS and optional galaxy-protection mask.
- `SmartStarMaskCandidate` — immutable pixel raster, explicit
  non-adopted state, diagnostics and confirmation-gated adaptation to
  existing band/grid-specific `SceneState`.

Implementation deliberately limits the supported search circle to
120 arcsec; the local affine WCS approximation requires real-object
acceptance before extending to larger fields. Current output is
**not** a SkyMaskDraft shape: irregular segmentation is a derived
per-native-WCS raster with its own provenance. The user-shared
`StarMaskSeed` and future per-band segmentation references must be
persisted together in a separately versioned scene manifest, rather
than silently polygonizing/flattening their boundaries.

## References

- Photutils source detection:
  https://docs.astropy.org/projects/photutils/en/stable/api/photutils.segmentation.detect_sources.html
- Photutils multilevel deblend:
  https://docs.astropy.org/projects/photutils/en/stable/api/photutils.segmentation.SourceFinder.html
- Photutils source-mask dilation:
  https://docs.astropy.org/projects/photutils/en/stable/user_guide/segmentation.html
- Photutils PSF modeling:
  https://docs.astropy.org/projects/photutils/en/stable/api/photutils.psf.PSFPhotometry.html

## Next gated implementation

GUI: new local-FITS mask canvas with click/draw ICRS anchor, pending
preview, threshold slider, comparison of per-band candidates,
accept/reject, mask history and keyboard undo/redo. Require offscreen
capture and real Windows WCS click interaction tests.

Science: PSF-informed candidate boundaries and empirical blank-sky/
injection validation, followed by explicit background-model iterations.
