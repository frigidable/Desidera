# Leo P / Scene & Masks — end-to-end logic audit (2026-11-10)

**Status: SCIENCE HOLD / architecture review.** No production photometry
adopted; PR #84 must stay Draft. This report checks the code on branch
`feature/scene-wcs-editor-sc1` (through `8f2962c`) against the full
`desidera.gui.scene-complete-qa/1` export from Leo P on 2026-11-10
17:00:23 UTC, native DR11 South `1554p180` in `g`. The snapshot
predates the recent Gaia-counterpart label fix. It contains metadata,
not original native image arrays. This is a **code and reported-metadata
audit**, not an independent pixel-level remeasurement.

## What is internally consistent

- Tractor screening: `5363 = 493 masked + 114 overlaps + 3330
  outside_cutout + 1426 unusable`. Categories are sequential code
  branches, not a census of all objects with centres outside the cutout.
- The 785-entry image inventory: `467 PROPOSED_EXCLUDE + 114
  PROTECTED_OVERLAP + 204 NOT_SCREENED_OR_UNFITTED`.
  The 26-entry difference between masked sources (493) and
  in-frame proposed-exclude entries (467) is caused by source *centres*
  outside the image with masks that enter the field; not a count bug.
- There are 27 positive G3 reference IDs among 785 entries, 5 among
  114 protected overlaps. The old `near_target_sources` falsely
  labelled all 114 as Gaia counterpart review due to `gaia_g_mag=0`
  placeholders. Fixed in `scene_contaminants.py` by commit
  `ffbe35e`, regression tests `ce4a8d0`; verify new GUI export.
- There are 42 `NOT_SCREENED_OR_UNFITTED` centres within the
  protected circle in addition to the 114 screening overlaps.
- `physical_membership=UNDETERMINED` and
  `photometry_action=NO_NEW_ACTION_READ_ONLY` on all SA1 entries
  is correct. `scene_draft.edits=[]`, `dirty=false` in this report,
  so *no feedback-induced saved edit has been demonstrated*.

## Confirmed architectural and methodological blockers

| ID | Severity | Evidence | Required change |
| --- | --- | --- | --- |
| L01 | **P0** | `SceneEditor._accept_isophote` persists an isophote candidate as `SkyMaskEdit(purpose='protect_target')`; `prepare_scene_contamination` then uses that edit to enlarge `protected` for the next run. Thus the result can change its own future selection policy. | Separate `isophote_proposal` (observation) from `source_policy` (decision); prohibit implicit feedback, test rerun invariance. |
| L02 | **P0** | All geometric source envelopes touching the provisional 60-arcsec protection circle are withheld from masking, including G3 stars. `propose_scene_isophote` calls the extractor directly on `~(native_hard | mask)` without validated foreground PSF subtraction or wings. | Source-level foreground/target/mixed hypotheses plus validated PSF modelling and uncertainty. |
| L03 | **P0** | The radial extent `63.415 arcsec` is estimated from the very same masked image whose source treatment switches near 60 arcsec. This is flagged in the report but not measured with alternate radii. | Freeze alternative mask scenarios, compare with identical image/sky, do not claim physical edge from coincidence. |
| L04 | P1 | SGA-radial GUI comparison counts only earlier masked/overlap lists: within its proposed ellipse there are 1 masked, 83 unresolved and also 36 unscreened centres, i.e. **120** Tractor centres, not 84 total. | Count all source-inventory classes for *each* candidate aperture; distinguish centre counts from footprint intersections. |
| L05 | P1 | Four nearby threshold/smoothing trials use one foreground mask and a constant local sky; they do not vary PSF, source membership, coverage or sky model. | Add controlled scenario matrix and held-out tests, not additional ad hoc confidence thresholds. |
| L06 | P1 | Isophote extractor performs thresholded Gaussian smoothing and 2-pixel binary opening. It does not estimate true outer surface brightness; it can link unmodelled foreground stars, nebular clumps and residual background. | Keep diagnostic segmentation separate from isophotal radial profiles/total-galaxy flux. |
| L07 | P1 | The GUI says 4-arcmin `view`, but `_cutout_pixels` deliberately requests a 1.60× padded square trimmed at native FITS boundaries. The reported 1466×1324 pixels at 0.262 arcsec/pixel cover ~**6.40×5.78 arcmin**; field core is 4 arcmin. | Show actual native image extent and separate 4-arcmin science field outline; report clipping and sky footprint. |
| L08 | P1 | Export has Tractor catalogue SHA and settings but no per-file IMAGE/INVVAR/NEXP/MASKBITS content digests; the `dataset_source` label is not a file identity. | Add immutable native quartet hashes, sky and policy fingerprints, dependency invalidation. |
| L09 | P1 | `sga_informed_radial_seed` uses an empirical, sky-subtracted smoothed annular contrast without calibrated correlated noise, plus fixed SGA centre/q/PA. | Validate independent sky/PSF/geometry alternatives, no claim of independent SGA radius. |
| L10 | P2 | `sga_radial.bins[].median_smoothed_rms` is actually median smoothed **flux divided by** smoothed RMS, a dimensionless contrast. | Rename field with an explicit schema version / compatibility path. |
| L11 | P2 | Old SC1 architecture documentation hard-codes 125 Leo P protected overlaps; current snapshot has 114. | Treat object counts as fingerprint- and scenario-specific, not stable documentation assertions. |

The masked field has 400,042 candidate-excluded pixels (20.61%),
164,759 protected pixels (8.49%), and 385 native hard-invalid
pixels. These are neither independently unbiased galaxy-light estimates
nor complete source-removal models.

## Isophote numerical cross-check

Reference settings: `g`, contrast 2.20, Gaussian smooth 4.0 pixels,
bridge opening 2 pixels. Four detected components have respectively
45,974 / 56,373 / 37,976 / 62,439 pixels; axis ratios
0.541 / 0.559 / 0.365 / 0.551; centroid offsets
39.9 / 50.2 / 33.1 / 49.3 px. Area ratio 1.644, min IoU 0.735,
and `LOW_CONFIDENCE_REQUIRES_EXPLICIT_OVERRIDE`. Successful detection
is **not** scientific validation.

The moment-based *area-equivalent* ellipse is
`a=43.09 arcsec, q=0.541, PA=162.37 deg`; it has a different
definition from the independently labelled (but **not statistically
independent**) SGA-guided radial candidate `a=63.415 arcsec`.

The SGA2025 **PARENT Viewer** supplies only the initial geometry
`a=14.29 arcsec, q=0.532, PA=150.46 deg`, with mutable source flags;
the ~63.4-arcsec extent is inferred from DESI pixels. At 63.42
arcsec the annulus has 95.9% coverage, 6/8 positive sectors,
median *dimensionless* smoothed contrast 0.812; at 65.04 arcsec
the contrast is 0.589, below the hard-coded 0.65 support cutoff.
The proximity to the 60-arcsec mask-policy boundary demands a
controlled test. It does **not** prove the boundary caused the result.

The direct local polygon contains ~97 Tractor centres (including
three G3 IDs; notable `OBJID=2035`, `1870`, `1876`). The SGA-guided
ellipse has ~120 centres, including all five G3 IDs in the protected
overlap inventory. The data do not yet distinguish their PSF-contaminant
flux from target light.

McQuinn et al. (2013) report for Leo P a resolved-stellar structural
scale of `a=70 ± 5 arcsec, ellipticity=0.52, PA=335 deg (modulo 180
→ 155 deg)`. This is useful external context, **not** a validation
of the Desidera surface-brightness segmentation or total photometry.

## Stage-gated repair and test order

1. **P0** Add an independent, versioned proposal storage type. A
   `candidate-isophote` must never be silently persisted as an active
   `protect_target` mask. Introduce a negative regression: identical
   inputs and source policy must produce identical outputs after merely
   saving an exploratory isophote; explicit source-policy changes
   must be separate and audited.
2. **P0** Implement source-level evidence review with real G3 and
   foreground PSF-wing uncertainty. Do not infer physical origin
   from Tractor TYPE, MASKBITS, SGA or proximity. Retain all 785
   entries and 42 faint protected centres for QA.
3. **P0** Run a *frozen-input* protected-radius and neighbour-policy
   sweep (e.g. 40/50/60/70/80 arcsec). Compare contour, SGA radial
   support and sector-wise background without changing the saved Leo P
   project. The sweep must demonstrate whether the 63.4-arcsec
   cutoff moves with source-policy boundary.
4. **P1** Implement 2D spatial sky and multiple sky train/holdout
   regions with complete valid-pixel and azimuthal coverage; forward
   model or independently mask and audit foreground stars.
5. **P1** Explicitly expose padded native WCS footprint, centre and
   4-arcmin science square; calculate a separate science-region
   coverage ledger and per-band content hashes.
6. **P1** Add independent native g/r/i/z detections, SA2
   cross-associations and alternatives for physical memberships.
   Tractor griz model fluxes are *not* equivalent to independent
   detections or PSF-matched aperture colours.
7. **P1** Calibrate methodological bias and coverage with
   injection/recovery across irregular resolved dwarfs, star-forming
   knots, varying sky, saturated wings and footprint gaps.

Do not use a new target radius, hand-selected star classification, or
the 70-arcsec literature scale to tune the current thresholds and
then call the same Leo P run an independent validation.

## References

- Desidera `SCENE_GEOMETRY_P12.md`, `PHOTOMETRY_QUALITY_VALIDATION.md`
- https://www.legacysurvey.org/dr11/description/
- https://www.legacysurvey.org/dr11/catalogs/
- https://www.legacysurvey.org/dr11/bitmasks/
- McQuinn et al. 2013: https://doi.org/10.1088/0004-6256/146/6/145

**Disposition:** Leave PR #84 in Draft. Do not accept an outer
isophote, a source classification, a multiband colour or a total
magnitude from this diagnostic run.
