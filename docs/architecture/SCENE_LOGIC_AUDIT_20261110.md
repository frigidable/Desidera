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


---

## Second-pass forensic audit — 2026-11-10, exact JSON + DR11 bit definitions

**This addendum supersedes any inconsistent terminology, severity or
counts in the first pass above.** Evidence: exact 1,895,875-character
\`desidera.gui.scene-complete-qa/1\` JSON exported at
\`2026-11-10T17:00:23.589755+00:00\`, all 785 source rows,
complete stored contour vertices and trial metadata, checked against
the exact branch implementation and the published DR11 bit tables.
Original per-pixel FITS arrays and an independent Gaia query were
not available for remeasurement. New code commits below describe
engineering repairs; the input JSON **predates** them.

### F1. Specific corrections to the previous audit

| Earlier wording or inference | Verified correction |
|---|---|
| All 156 flagged sources are "inside the 60-inch circle". | 102 \`PROTECTED_OVERLAP\` centres and 42 \`NOT_SCREENED_OR_UNFITTED\` centres are physically **inside** the provisional circle: **144 centres**. Another **12** screened sources have centres outside but geometric source envelopes touching the circle. Thus 156 **inside or intersecting**, not 156 inside. |
| The candidate isophote "automatically" feeds a new protection mask. | \`_accept_isophote\` needs an explicit **Add proposed region** user action and, for low confidence, a separate explicit override. It saves a \`protect_target\` *draft*, not approved photometry. **The feedback path exists if the draft is subsequently used**, but it has not occurred in the supplied Leo P run (\`edits=[]\`). Engineering/design P1 risk; P0 if such a feedback-derived result were adopted without separate validation. |
| The provisional sky annulus might be outside the loaded FITS in unknown amount. | Derive its actual geometry: cutout \`1324×1466\`; target center approximately \`(732.23,733.03)\` pixels in origin-zero coordinates, using TAN/WCS with **PC×CDELT**, not CDELT alone. Sky radii \`423.68–622.28\` pixels, approximately \`111.0–163.0\` arcsec. **1.28%** of the ideal annulus is cut by image boundaries. Of 644,248 geometrically available ring pixels only **487,180** are used: **75.62%** of available annular pixels. Remaining ~24.38% is excluded by native validity / provisional masking, whose spatial selection is unvalidated. |
| The FITS cutout must exactly equal a 4-arcminute view. | Code intentionally uses \`_cutout_pixels(...padding=1.60)\` to include the sky surroundings; native cutout spans approximately **6.40×5.78 arcmin**. This is **documented figure-padding design**, not an unintended WCS or extraction failure. It remains a P2 GUI labelling/provenance issue because science field and display/sky ring footprint should be marked separately. |
| \`OBJID=2035\` is saturated in G. | \`MASKBITS=2064 = MEDIUM (bit 11) + SATUR_Z (bit 4)\`, **not SATUR_G**. Its G-band saturation must be checked in G native pixels; the published catalogue alone does not establish it. Its positive Gaia astrometry and PSF-wing risk remain. |
| SGA source row \`OBJID=2059\` can be treated as an independently fitted Sérsic component. | \`FITBITS=530 = FIT_BACKGROUND + FROZEN + LARGEGALAXY\`: it is **frozen**, and SGA context is not statistically independent morphology evidence. \`OBJID=2076\` is also **FROZEN** (\`FITBITS=18\`). |
| A science- or code-level "P0" is the same category. | Distinguish **SCIENCE HOLD** (no adopted total/outer radius with existing evidence) from demonstrated correctness defects in preview-only SC1. The existing preview does not claim scientific approval, and no actual adopted total was shown. |

### F2. Independently recomputed ledger

- Image: \`1324×1466=1,940,984\` native pixels.
  Candidate-excluded = **400,042**; native hard-invalid = **385**;
  usable = **1,540,653**. By inclusion-exclusion, **96** native
  hard-invalid pixels also lie under a candidate exclusion. These
  counts are internally consistent and do not quantify flux loss.
- Screening of **5,363** local Tractor rows:
  \`493 masked + 114 protected overlaps + 3,330 outside screen + 1,426
  rejected/unfitted = 5,363\`. The 3,330 "outside" label is only
  the **post-eligibility geometry branch**, not all off-frame sources.
- Old SA1 record: **785 in-frame source centres**, partitioned
  \`467 proposed-exclude + 114 overlap + 204 not-screened/unfitted\`.
  Of 114 overlaps only **5** have positive \`G3/REF_ID\`; **27**
  \`G3\` references exist over all 785 rows. \`gaia_g_mag=0\` wrongly
  labelled 114/114 as Gaia in the original JSON; fixed earlier but
  must be reconfirmed on a fresh input.
- **NEW CONFIRMED TRACEABILITY DEFECT:** all 493 proposed masking
  source IDs are available in \`stats.masked_catalog_centres\`, but
  **26** have centre coordinates outside the native image and were
  omitted from the 785-entry SA1 register despite affecting its
  pixels. Independently verified all 26 have off-frame centres.
  Repaired on this branch via adding a centre-in-native-frame
  provenance bit and separate counts. Expected next export **811
  relevant source rows = 785 native-centre rows + 26 off-frame mask
  contributors**; projected pixel mask must remain identical. Legacy
  field \`total_in_native_frame\` remains 785. New branch changes must
  pass all tests and live Leo P re-export before claiming complete QA.
- The proposed SGA-guided ellipse has **120** centre-in-ellipse Tractor
  records, remeasured from source RA/Dec at the reported centre, q and
  PA: 83 overlap, 36 unfit/unscreened, 1 masked; 5 G3 matches.
  Direct polygon has **97** centres: 71 overlap and 26 unfit, 3 G3.
  These positional numbers are not area-overlap counts or flux
  decompositions.

### F3. Unreported DR11 fitting-history evidence — software defect corrected

The \`scene_source_inventory.FITBIT_NAMES\` dictionary previously
decoded only bits 0–7; the official DR11 table also defines
\`TYCHO2\` (8), \`LARGEGALAXY\` (9), \`WALKER\` (10), \`RUNNER\` (11),
\`GAIA_POINTSOURCE\` (12) and \`ITERATIVE\` (13). Leo P contains
**1 / 4 / 1 / 23 / 40** sources with bits 9/10/11/12/13 respectively.
\`MASKBIT_NAMES\` also omitted published bits 0, 8, 9, 19.

- \`OBJID=2059\`: \`SER\`, \`REF_CAT=L4\`, \`FITBITS=530\`
  (\`FIT_BACKGROUND|FROZEN|LARGEGALAXY\`), source-level SGA bit 12.
  Does not justify treating its \`shape_r=16.60''\` as an independent
  photometric fit or validating the Leo P target population.
- \`OBJID=2076\`: central \`REX\` at 0.073 arcsec, \`FITBITS=18\`
  (\`FIT_BACKGROUND|FROZEN\`). Not an independent target core fit.
- \`OBJID=1762\`: \`SER\`, \`FITBITS=11268\`
  (\`HIT_RADIUS_LIMIT|WALKER|RUNNER|ITERATIVE\`), high-priority fit-QA
  failure candidate; cannot infer membership.
- \`OBJID=2035\`: \`FITBITS=4290\`
  (\`FIT_BACKGROUND|MEDIUM|GAIA|GAIA_POINTSOURCE\`), \`MASKBITS=2064\`
  (\`SATUR_Z|MEDIUM\`). Gaia identity and high-significance astrometry
  justify foreground review but not a numerical G-band wing correction.

The complete published DR11 dictionaries and a synthetic regression
for these real-row examples were added to the existing branch. This
restores evidence visibility only: **no physical source classes or
exclusion masks are changed**.

### F4. Distinguish pixel component, filled outline and area ellipse

The reference trial \`component_indices\` has **45,974** supported
pixels. Independently applying the shoelace formula to the **192
stored native contour vertices** gives a filled outer-polygon area
of approximately **49,560.5 native pixels**, **7.80% greater** than
the thresholded component's pixel count. The three alternative trials
show filled-polygon/threshold-pixel area ratios approximately
**1.092 / 1.078 / 1.100**. The area-equivalent moment ellipse is
**defined** to have area 45,974 pixels, not 49,560.
Also, `propose_local_isophote` computes q/PA from the component's
**weighted barycentre** but constructs `EllipseSeed(x_pix=x0,y_pix=y0)`
at the original target anchor, *not* at that barycentre. For the
reference trial, the centroid displacement is 39.9 px = ~10.5 arcsec;
for the sensitivity trials it reaches 50.2 px = ~13.2 arcsec.
Thus the displayed "moments ellipse" is a **target-anchored,
area-equivalent surrogate** with an independently measured orientation,
not a self-consistent freely fitted moment ellipse. It should never
be labelled an independently fitted isophote.

The discrepancy can result from intentional 2-pixel contour dilation,
smoothing and holes/gaps inside the outer ring. The stored contour
**is not the same object** as the connected-pixel segmentation:
using it as \`protect_target\` also changes treatment of all interior
pixels, including possible masked source holes. This strengthens the
need to separate diagnostic contour selection, sky/PSF source policy
and accepted measurement aperture. It does **not** by itself quantify
any photometric bias or prove all extra area is contaminated.

### F5. Sky / outer-radius precision and what remains unknown

Current native \`g\` sky estimate in the source JSON:
\`local_sky≈2.14355×10⁻⁵ nMgy/pixel\`,
\`raw_sky_rms≈2.85187×10⁻³ nMgy/pixel\`,
\`smoothed_sky_rms≈2.71787×10⁻⁴ nMgy/pixel\`.
These values are conditional on a masked, spatially incomplete
annulus and an already upstream sky-subtracted coadd.
The 24.38% annular validity exclusion is not uniformly distributed
by assumption; its sectors need checking from the original arrays.

At outer annuli:
\`a_outer=60.163''\` (contrast 0.970, 7 positive sectors),
\`61.789''\` (0.828, 8),
\`63.415''\` (0.812, 6),
\`65.041''\` (0.589, 5, **fails fixed median threshold .65**),
and \`66.667''\` (0.559, 6, fails).
Coverage remains ~95.7–98.1% across these rings.
**This specific support failure is driven by the code's median contrast
rule, not a documented collapse of radial coverage.**
Whether the contrast drop is astrophysical, affected by the
60-arcsec protection boundary, or due to residual sky/foreground
light CANNOT be decided without native pixels and frozen mask/PSF/sky
scenario reruns. The existing QA flag is appropriately cautious.

### F6. Revised prioritization and still-unverified gates

**Demonstrated repaired engineering defects (pending live/regression
confirmation):** the Gaia zero-placeholder source classification,
incomplete DR11 bit-name decoding, and absent per-source IDs for 26
off-frame mask contributors. They must not silently modify original
image pixels or proposed mask geometry.

**P1 current architectural hazard:** the user's manual adoption of a
low-confidence contour stores it as \`protect_target\`, which changes
future source eligibility. This is *an explicit GUI action*, not
unobserved automated self-feedback. No such edit existed in this JSON.
The safe scientific design is to store an \`isophote_proposal\` separately
and require a distinct versioned source-treatment policy change.

**Science HOLD even when all GUI tests pass:** no validated
sky/PSF-wing decomposition, SGA-independent outer extent, source
membership, missing-light reconstruction, selection-aware uncertainty
or independent injection/recovery. The Q0–Q6 scientific adoption
contract remains in force.

**Still unavailable for independent verification:** per-pixel
IMAGE/INVVAR/NEXP/MASKBITS FITS quartet; local Tractor FITS bytes;
alternative mask radius sweeps; PSF wing residuals for G3 stars;
Gaia proper-motion/quality crosschecks at their observing epochs;
multi-band image detections and calibrated flux/colour recovery.
The provided JSON permits accurate code and *metadata* checks but
cannot certify an astrophysical outer isophote.

### Sources and reproducibility

- Native original JSON fingerprint:
  \`8e8e3f75af868f3ef7cf4c4c9d7ca28210cbca87ef27c98351a6ab237010562c\`
  (**Tractor catalogue hash**, not native IMAGE SHA).
- Official \`MASKBITS\` and \`FITBITS\`:
  https://www.legacysurvey.org/dr11/bitmasks/
- Official Tractor catalogue flux and shape semantics:
  https://www.legacysurvey.org/dr11/catalogs/
- Upstream spline sky subtraction / image-stacking:
  https://www.legacysurvey.org/dr11/description/
- External Leo P stellar geometry comparison:
  https://doi.org/10.1088/0004-6256/146/6/145

**Disposition after second pass:** branch remains Draft, with
scientific measurements blocked. Next experiment should use a newly
exported complete QA report from the amended build, verifying
**785 native centres + 26 off-frame contributors**, **493 masked / 114
overlaps**, new high FITBITS decoding, and unchanged raw pixel
stats. The fact that synthetic tests pass is not enough for Leo P
science acceptance.


---

## Third-pass live GUI verification — 2026-11-10 19:07 UTC

**Input:** full 2,042,300-byte JSON exported from Leo P via
\`Copy Scene QA\` at \`2026-11-10T19:07:11.196149+00:00\`;
file SHA-256
\`22e7c576a11139d21bf08977a6c1802a87228068093f522ffcdac8e70cd402a7\`.
The report is **not** a FITS-file digest. Same object/band/Tractor hash
\`8e8e3f75af868f3ef7cf4c4c9d7ca28210cbca87ef27c98351a6ab237010562c\`,
project \`Leo_P\`, DR11 South brick \`1554p180\`, native G cutout
\`1324×1466\`, mask draft \`edits=[]\` and \`dirty=false\`.
This is the first supplied **live GUI re-export** following the
source-evidence completeness and bit-decoding corrections.

### Newly confirmed passes

1. **Exact source identity ledger:** 811 distinct source IDs with no
   export truncation; **785 centres within native FITS** and **26
   off-frame contributors**. Class counts (native centres)
   \`114 protected / 467 proposed exclude / 204 unfit\`;
   class counts (all relevant sources)
   \`114 protected / 493 proposed exclude / 204 unfit\`.
   In particular, the 26 exterior sources are now visible and
   individually traceable; the native pixel mask is unchanged.
2. **Classification fix:** \`near_target_sources\` now contains
   **5 GAIA_COUNTERPART_REVIEW + 109 UNRESOLVED_TARGET_BLEND**.
   The 811 catalogue entries contain 29 positive G3 references:
   **27 in-frame and two exterior**; 5 are protected overlaps.
   All entries still have \`physical_membership=UNDETERMINED\` and
   \`photometry_action=NO_NEW_ACTION_READ_ONLY\`.
3. **DR11 FITBITS/MASKBITS:** independent decoding of every source's
   integer bitfields against the official DR11 tables matches its
   exported bit-name arrays **811/811**, with **zero unknown high
   bits**. All-bit counts (including exterior contributors):
   \`ITERATIVE=42\`, \`GAIA_POINTSOURCE=24\`, \`WALKER=5\`,
   \`RUNNER=1\`, \`LARGEGALAXY=1\`. Counts differ slightly from the
   earlier in-frame-only 785-row snapshot because the new export
   includes the 26 exterior sources.
4. **Mask and identity arithmetic:** \`5363 = 493+114+3330+1426\`;
   native total \`1,940,984\` px, proposed excluded \`400,042\`,
   native hard-invalid \`385\`, usable \`1,540,653\`. Inclusion/
   exclusion implies 96 pixels both hard-invalid and proposed
   excluded. Protection circle \`164,759\` px agrees with
   \`pi*(60/0.262)^2=164,759.24\` to integer rounding.
5. **Complete provenance consistency within this run:** the
   catalogue SHA matches in \`contaminant_screening\`,
   \`source_inventory\` and \`isophote_proposal\`; all three mask-draft
   digests match the saved draft. SGA lookup and selected prior share
   the same source SHA; SGA **geometry enabled**, but **SGA automatic
   protection disabled**. No worker was running.
6. **Internal coordinate consistency:** using the exported TAN/FK5
   header (including PC×CDELT) to recompute the pixel→world mapping
   independently for all 811 source rows and 192 contour vertices
   reproduces the exported positions at ~\`1e-8 arcsec\` numerical
   tolerance. This validates **internal WCS bookkeeping**, not
   independent calibration against external astrometry.
7. **Contour arithmetic:** the reference component contains
   **45,974 pixels**; shoelace area of the 192-vertex outer contour
   is **49,560.505 px**, ratio **1.0780116**. The area-equivalent
   target-centred ellipse has calculated area **45,973.822 px**,
   consistent with its name but **not** with the filled polygon.
   Distinguish all three area definitions.
8. **Independent spatial source counts:** direct native polygon
   encloses **97 catalogue centres** =
   \`71 protected + 26 unscreened\` (G3 IDs 1870/1876/2035).
   SGA-guided ellipse encloses **120 centres** =
   \`83 protected + 36 unscreened + 1 proposed masked\`
   (G3 IDs 1763/1870/1876/2034/2035). The existing GUI summary
   reports only \`83+1\`, i.e. omits **36** unscreened centres
   from the sentence, though they exist in the full Source evidence.
9. **Sensitivity arithmetic:** four detected contours have areas
   \`45,974 / 56,373 / 37,976 / 62,439\`;
   independently recomputed area ratio \`1.644170\`,
   q spread \`0.194475\`, min recorded IoU \`0.734775\`,
   maximum centroid drift \`50.211px\`.
   The appropriate output remains
   \`LOW_CONFIDENCE_REQUIRES_EXPLICIT_OVERRIDE\`.
   IoU itself cannot be independently recomputed from this export:
   component pixel *indices* are deliberately omitted.

### Interpretation that remains blocked

- **Gaia foreground:** sources 2035, 1870 and 1876 are G3 matches
  with high astrometric significance within the **direct** polygon;
  these three alone contribute ~652.38 nMgy in **Tractor G-band model
  fluxes at their centres**, among ~787.22 nMgy obtained by naively
  summing the *model flux values* of 97 centred catalogue components.
  The resulting ~83% is **NOT the contaminant fraction of observed
  galaxy flux**: fits overlap, profiles extend beyond polygon edges,
  and Tractor's dwarf deblending is unreliable here. It is only a
  strong diagnostic signal demanding PSF/wing modelling before
  integrating observed pixels.
- **Non-independent radial edge:** supported SGA-guided radial
  radius \`63.415''\` is the last accepted annulus with median
  smoothed flux/RMS contrast \`0.812\`; the next at \`65.041''\`
  has \`0.589\` vs an empirical threshold \`0.65\`, while coverage
  remains \`~96%\`. The 60-arcsec source-protection boundary may
  influence contrast, but causal influence is **not demonstrated**.
  \`sga_radial.masked_fraction\` is ~\`0.77%\`, which denotes
  *provisional masking*, not freedom from foreground-star flux.
- **Unchanged feedback hazard:** \`Add proposed region\` can, on
  explicit manual confirmation, store a segmented contour as a
  future \`protect_target\` geometry. No such action occurred
  in this snapshot. The discrepancy between source pixels,
  filled contour and target-centred area ellipse reinforces the
  requirement to split \`proposal\`, \`policy\` and \`measurement\`.
- **Sky / image science:** original IMAGE/INVVAR/NEXP/MASKBITS
  pixels, PSF models, forced source-render images, independent
  residual sky measurements and full pixel noise covariance are
  **absent** from this JSON. No physical isophotal radius,
  integrated flux, decontaminated colour or SBF distance can be
  validated from the metadata alone.

**Decision:** source evidence QA and bit mapping pass within this
single-band, single-object, current-snapshot scope.
**Scene/Mask and outer-isophote science adoption remain on HOLD.**
Do not tune the detection threshold to reproduce the reference
70-arcsec stellar extent. Next controlled tests should hold the
native FITS, SGA prior and source catalog fixed while varying
protection policy, PSF-wing treatment and sky definitions, and
verify source-policy fingerprints/independent holdout regions.

This record is intended as a reproducible **metadata-level science
audit**, not a claim of end-to-end FITS-level validation.


---

## Fourth-pass protection-radius sweep — Leo P (2026-11-10 21:26:22 UTC)

**Input:** actual 2,064,599-byte \`Copy Scene QA\` JSON exported from
Leo P, SHA-256
\`b2b50d1825075490bef95a1efde214b25b5e2299c4dd7b01363733bf88c6abcd\`.
Read-only five-scenario \`desidera.gui.protection-radius-sweep/1\`,
DR11 south brick \`1554p180\`, \`g\` native IMAGE \`1324×1466\`,
unchanged in-memory image, Tractor catalogue SHA-256
\`8e8e3f75af868f3ef7cf4c4c9d7ca28210cbca87ef27c98351a6ab237010562c\`,
SGA2025-PARENT geometry ID 5065058 / SHA-256
\`4ed9faa3b3ea039ef714d766e657698e0c12246dff831b7d7899f812973506e1\`.
**SGA geometry ON, SGA core source protection OFF.**
Contrast 2.2, smoothing 4.0 px, binary opening 2 px.
Draft \`revision=2\`, \`edits=[]\`, \`dirty=false\` and digest
\`2727f8959a9bcd15338f5a3c6abbf12eafe569e34450a14c4c59643cfd502619\`.
No science adoption or project-file writes observed.

### Reconciled results

| Protected circle | Proposed masked IDs | Protected overlap IDs | Excluded pixels | Connected component px | q from moments | SGA-guided radial endpoint (arcsec) | Component IoU vs 60″ |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 40″ | 540 | 67 | 432426 | 39028 | 0.692313 | 53.6591 | 0.848593 |
| 50″ | 517 | 90 | 417589 | 43397 | 0.583408 | 63.4153 | 0.943947 |
| **60″** | **493** | **114** | **400042** | **45974** | **0.541076** | **63.4153** | **1.0** |
| 70″ | 478 | 129 | 388230 | 46102 | 0.540712 | 63.4153 | 0.996704 |
| 80″ | 459 | 148 | 371674 | 46102 | 0.540712 | 63.4153 | 0.996704 |

Independent integer-ledger checks passed:
- Each mask-preparation scenario preserves 607 screened eligible
  sources \`masked+protected\`. Newly excluded IDs relative to
  60″ are respectively 47/24/0/0/0; newly protected-overlap IDs
  0/0/0/15/34. ID changes are nested, unique and agree with
  baseline in-frame source screening. 811 relevant source entries
  and identical catalogue digest are reported in all five cases.
- XOR mask differences relative to 60″ are
  \`32384,17547,0,11812,28368\` px. Because proposed exclusion
  masks nest with radius, the exported mask Jaccards are exactly
  the ratio of the smaller mask area to the larger area in every
  case; no anomalous mask topology or pixel-accounting mismatch.
- Native hard-invalid pixel overlap with excluded pixels is
  **96** at radii 40/50/60″ and **0** at 70/80″; this explains
  the small difference from simply subtracting total hard pixels
  and excluded pixels from the total image area.
- Relative to 60″, connected-component pixel area changes by
  **−15.109%, −5.605%, 0%, +0.278%, +0.278%**. Solving the
  recorded IoU against integer component sizes yields these
  intersections with the 60″ component:
  \`39020,43397,45974,45962,45962\` pixels respectively.
  At 70″, even though the masks monotonically *shrink*, the new
  thresholded region adds 140 px and loses 12 px. This is
  consistent with recomputed normalized smoothing, thresholding,
  and morphological opening; connected components need not
  themselves be strictly nested.
- The 60″ base row reproduces the separately saved
  \`isophote_proposal\` exactly for component size, q, centroid,
  sky, smoothed RMS and SGA radius. All 5 scenarios report
  \`LOW_CONFIDENCE_REQUIRES_EXPLICIT_OVERRIDE\`.

### Physical interpretation, source-level controls

The simplistic hypothesis that the reported 63.415″ SGA-guided
endpoint is directly or automatically equal to the protection
radius **fails**: it remains exactly 63.415″ for protection
50/60/70/80″. But masks matter substantially at radius 40″,
when endpoint drops by **9.7562 arcsec = exactly 6 radial
bins** (bin step ~1.6260″), area drops ~15.1%, and q shifts
from .541 to .692. These tests demonstrate **conditional mask
sensitivity** but do not identify the causal pixels or establish
an astrophysical outer edge.

The 40→50″ change simultaneously restores protection to 23
source identities; the two particularly diagnostic Gaia G3
sources are:
- \`OBJID=1763\`: g Tractor model flux **17.152 nMgy**,
  circular offset 48.317″, SGA-coordinate elliptical radius
  **58.62″**, strong Gaia proper-motion signal.
- \`OBJID=2034\`: g Tractor model flux **1.404 nMgy**,
  circular offset 53.506″, SGA ellipse radius **58.55″**.
- \`OBJID=1762\`: a nearby \`SER\` component **0.73″**
  from 1763, also changes category at 40→50″, with
  FITBITS \`HIT_RADIUS_LIMIT|WALKER|RUNNER|ITERATIVE\`.
  This is a fitting/deblend-review target, **not** proof that
  either source belongs to the galaxy.

At 70″, newly protected Gaia source \`OBJID=2074\`
(g~6.155 nMgy). At 80″ two more are additionally protected,
\`OBJID=1796\` (g~34.975 nMgy) and \`OBJID=2226\`
(g~24.945 nMgy). Expanded protection can preserve **known
foreground-star candidates** rather than remove contamination.
Conversely, the bright central Gaia sources
\`2035\` (g~607.404 nMgy, offset ~19.9″),
\`1876\` (g~39.101 nMgy, ~29.1″) and
\`1870\` (g~5.876 nMgy, ~26.6″) remain protected in
**all five scenarios**. They are a common-mode potential bias,
not tested by this radius sweep.

In all five runs, the reported robust local sky
\`2.1435495000332594e-05 nMgy/pixel\`, smoothed sky RMS
\`0.00027178690896814656 nMgy/pixel\` and number of sky
pixels \`487180\` are **bitwise/numerically identical**.
Thus no observed change in the reported *scalar* sky summary
explains the inter-scenario differences. Sky modelling/systematic
bias and source PSF wings were **not** tested independently.

The reference SGA radial annulus at 63.415″ has median smoothed
contrast 0.812 with coverage .959; the next annulus at 65.041″
has contrast .589, below the fixed acceptance cutoff .65, with
coverage .965. The sweep exports only one detected endpoint per
radius, not the **full radial-bin sequence in each scenario**.
Therefore the 40″ failure mechanism cannot be diagnosed from
this JSON alone. It would be premature to attribute the
53.66→63.42″ transition to the two named Gaia sources.

### Next controlled test, in strict priority order

1. Keep the input IMAGE, local sky, SGA geometry and catalogue
   frozen, but **export all SGA annular bins (coverage, eight
   sector values, median contrast, support reasons) for each
   protection scenario**. This will identify where 40″ loses
   continuity relative to 50″.
2. Perform source-level leave-one-out / small-group controlled
   tests on IDs \`1762,1763,2034\` with explicit mask-policy
   alternatives. The baseline 60″ core Gaia set
   \`2035,1876,1870\` must separately receive PSF-wing models
   or clearly flagged exclusion experiments. Do not use
   Tractor \`TYPE\`, Gaia ID or SGA context to decide target
   membership automatically.
3. Independently vary local/2D sky models, sky holdout sectors,
   effective coverage and PSF-wing treatment. Then perform
   star/galaxy injection–recovery on real coadds, without tuning
   a threshold to reproduce the LBT resolved-star 70″ scale.

**Disposition:** software experiment and internal arithmetic pass;
**photometric science HOLD continues**. Radius 63.4″ is **robust
to these particular protective-radius variants for 50–80″**
but is neither a validated total/limiting isophote nor free
of star-wing, deblending and residual-sky biases.
