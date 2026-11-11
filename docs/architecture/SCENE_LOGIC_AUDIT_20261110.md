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


---

## Fifth-pass source-policy separation and legacy isophote repair — 2026-11-11

**Basis:** direct reinspection of \`conservative_catalogue_mask\`,
\`prepare_scene_contamination\`, \`build_scene_source_inventory\`,
\`SkyMaskEdit\` / \`evaluate_mask_draft\`, the GUI worker and acceptance
paths, original DR11 documentation, and the 2026-11-10
21:26 UTC complete Leo P QA report. This section documents actual code
changes in branch \`feature/scene-wcs-editor-sc1\`, not an adopted
foreground model.

### Source policy: what was broken, and what remains unresolved

The old screening rule
\`if np.any(disk & protected): near+=1; continue\`
vetoes the **entire proposed source footprint** on even one
protected pixel. The radius encodes no physical source membership and
cannot distinguish a real galaxy H II region from a Gaia foreground
star or background system. The apparent stability of the 63.4-arcsec
radial endpoint under 50–80 arcsec radius sweeps does not test the
five Gaia astrometric foreground candidates common to the protected
region (OBJID 1763, 1870, 1876, 2034, 2035). This behaviour is
retained **as explicitly unvalidated screening**, pending a proper
PSF/model-based source policy. It must not be described as foreground
decontamination or science-ready masking.

Independent source evidence is now reported with
\`foreground_evidence_assessment\`,
\`proposed_source_treatment\` and
\`proposed_treatment_is_active=False\`, while
\`physical_membership=UNDETERMINED\` and
\`photometry_action=NO_NEW_ACTION_READ_ONLY\` remain unchanged.
The proposed foreground review requires an identified G3/REF_ID and
a positive parallax significance or absolute PM-component
significance >=5 from inherited Tractor reference astrometry,
provided there is no strongly negative parallax quality alert.
This is **not a posterior probability** or a substitute for
Gaia astrometric-solution quality verification. The protection
worker records explicit IDs/counts of such sources currently left
unmasked, and the GUI exposes a dedicated filter. No new pixel
exclusion or native FITS modification is permitted by these labels.

**Negative-parallax safety rule:** an absolute-significance calculation
by itself would misclassify a strongly negative parallax as foreground
distance evidence. Such sources now receive
`GAIA_NEGATIVE_PARALLAX_QUALITY_REVIEW` and
`GAIA_ASTROMETRY_QUALITY_ALERT`, never an automatic foreground
proposal, even if the nominal proper-motion component also has high
significance. All treatment remains read-only. See Rybizki et al.
(2022, MNRAS 510, 2597), https://doi.org/10.1093/mnras/stab3588,
on spurious Gaia astrometry and negative-parallax quality warnings.

From the supplied unmodified Leo P baseline, the flagged internal
source IDs are **1763, 1870, 1876, 2034, 2035 (5)**.
They require external Gaia quality cross-checks and validated
wavelength-dependent PSF wing treatment. In particular,
\`OBJID=2035\` has large, nonzero reference parallax and a bright
Tractor G model but its MASKBITS do not prove G-band saturation.

### Confirmed separate P0 source of feedback, now mitigated

Before this change, an explicit \`Add proposed region\` click
persisted tentative \`contour\`, moment-ellipse or \`sga_seed\` as
\`SkyMaskEdit(purpose="protect_target")\`. This activated additional
protection in *subsequent* catalogue masking and introduced an
unjustified feedback loop. While the 2026-11-10 Leo P baseline
contained no saved edits, previously saved project drafts elsewhere
could contain this harmful state.

Corrected implementation:

- New purpose \`isophote_candidate\` in the existing mask-draft
  schema. Newly retained isophotes become **non-masking** scene
  annotations, never a new active protection region. They remain
  revisioned, opt-in, per-band scoped by default, redrawable and
  explicitly shareable between bands.
- Legacy \`protect_target\` regions generated automatically by SC1
  are detected by \`sc1_iso_*\` ID and/or provenance markers
  \`AUTO_ISOPHOTE_PRELIMINARY\` or
  \`SGA_INITIAL_GEOMETRY_DRAFT\`, and are **interpreted as inactive
  isophote candidates** by both the Scene preparation path and the
  generic mask-draft evaluator. They are not mutated or automatically
  resaved; previously recorded project digest/revision remain intact.
  The GUI labels them \`legacy isophote · INACTIVE\`.
- Explicit user-authored \`protect_target\` edits retain their previous
  masking-policy effect. Explicit user-authored \`exclude\` keeps
  priority over target protection even inside the circle. Original
  native-hard pixels remain inviolable.

Regression tests cover serialization, band scoping, manually authored
protection, foreground manual exclusion inside the core, previously
saved auto-isophote drafts, truncated old provenance comments and GUI
reopen/rendering.

**Scope boundary:** this corrects *proposal-to-policy feedback*, not
the original source-screening veto itself. No Gaia candidate is
automatically masked or PSF-subtracted. No source treatment can become
a science-adopted mask before calibrated star models and explicit
uncertainty-aware decision provenance exist. The science HOLD remains.

### Next source-policy architecture (pending, NOT implemented)

Separate three immutable layers per native dataset/band:
1. \`SourceEvidence\`: detections, full Tractor metadata, verified
   Gaia identity/epoch/quality and SGA contextual geometry; preserves
   every candidate including below screening threshold.
2. \`SourceTreatmentDecision\`: researcher-verified hypothesis
   (galaxy component / foreground / background / ambiguous),
   source identity, confidence, evidence references, explicit
   policy action and per-band scope, independent of any circular
   target-protection geometry.
3. \`TreatmentRealization\`: rendered PSF/Tractor model with uncertainty,
   missing-wing bounds or reviewed geometric exclusion, fractional
   coverage and source-specific residual QA. A classification alone
   never automatically becomes a pixel mask.

A validated foreground star should be eligible for subtraction or
exclusion even within the provisional target circle. Ambiguous
overlaps must be propagated as alternative treatment scenarios rather
than silently included or excluded. Star-forming knots and genuine
resolved stellar light must not be removed merely because the Tractor
reported multiple PSF/REX components.

### Required revalidation before any scientific adoption

- Existing Leo P \`edits=[]\`: new build must preserve the
  baseline \`493 masked / 114 overlaps / 400042 excluded pixels\`,
  original source-catalogue SHA and five-radius sweep as before.
- Synthetic old \`sc1_iso_*\` accepted drafts: must not alter
  \`protected\`, \`auto_mask\` or resulting contours on reprocessing.
- New \`isophote_candidate\` after Save → Close → Reopen: must not
  alter source policy or undo native-hard flags.
- Manual \`protect_target\`, \`exclude\`, \`restore_auto\` must
  keep their documented priorities and band scopes.
- Gaia/source review flags must remain unadopted evidence and never
  silently cause new mask bits.
- Future independent PSF-wing and spatial-sky systematics must be
  tested with native FITS and injection/recovery; JSON diagnostics
  alone cannot supply this validation.


---

## Sixth-pass experimental inner Gaia masks — Leo P (2026-11-11)

**User-observed defect:** a post-fix screenshot and actual \`Copy Scene QA\`
generated \`2026-11-10T21:55:57Z\` still displayed no
purple/orange masks over the five foreground-star *candidates* in the
target core. The report explicitly contained
\`PROTECTED_GAIA_ASTROMETRY_STARS_NOT_REMOVED\` and
\`foreground_astrometry_protected_review_ids\` =
\`1763, 1870, 1876, 2034, 2035\`.
Baseline, unedited Leo P source stats were exactly 493 proposed
masked and 114 protected overlaps with 400042 masked pixels.
**This was expected from the previous code; earlier changes added
read-only evidence, not an internal masking engine.**

### New, deliberately opt-in preview mode

The SC1 GUI now exposes **Include inner Gaia candidates · QA only**
(default OFF). On the next **Prepare neighbour masks**, the native
same-brick DR11 Tractor catalogue is screened by the same source-type
/ positive SNR / flux and native FITS WCS rules as before, with one
additional *conditional* eligibility rule for sources whose proposed
footprints overlap target protection:

- Tractor PSF source, genuine \`REF_CAT=G3\` and positive \`REF_ID\`;
- nominal Gaia parallax positive and >=5-sigma OR one absolute
  proper-motion component >=5-sigma, with no significantly
  negative-parallax quality alert;
- \`gaia_duplicated_source\` false;
- the source's existing geometric PSF envelope actually intersects
  the current target-protected pixel set;
- source was otherwise eligible for catalogue screening.

These five *existing* Leo P Tractor/Gaia records satisfy that
provisional eligibility criterion on the supplied snapshot, so a
current-band preparation with the checkbox enabled should yield
**5 orange source-footprint previews**: 498 total masked candidate
sources and 109 remaining unresolved overlaps (subject to the same
files and untouched project inputs as the baseline). Exact pixel count
and isophote changes **must be measured**, not inferred from these
identity counts.

The native \`gaia_inner_preview_mask\` is a separate Boolean layer;
preview pixels are overdrawn in orange over ordinary purple
catalogue candidate exclusions. The source inventory uses
\`GAIA_FOREGROUND_MASK_PREVIEW\` and the diagnostic-activity flag while
keeping \`physical_membership=UNDETERMINED\`,
\`photometry_action=NO_NEW_ACTION_READ_ONLY\` and
\`proposed_treatment_is_active=False\`. The layer is included in
the current provisional isophote calculations and optional radius
sweep, and the toggle state, per-source IDs and pixel counts enter
\`Copy Scene QA\`. Turning the checkbox OFF and re-preparing returns to
the exact baseline; neither state is adopted or stored as a project
mask decision.

### Scientific limitations, not resolved by this mode

- DR11 Tractor circular magnitude-dependent source envelopes
  **are not a PSF model, subtraction, wing prediction or calibrated
  profile mask**. The loss of actual galaxy light beneath a foreground
  star is not quantified by zeroing these pixels.
- Gaia astrometric solution fidelity, EDR3 reference epoch and
  reliability are not independently verified. Even extremely large
  formal significance need not imply a valid clean deblend or PSF.
- DR11 \`PSF\` is a fitted *component*, not a unique physical
  star; close pairs and shredding require independent inspection.
- Unknown H II complexes, overlapping background galaxies and
  non-Gaia contaminants are not automatically classified/masked.
- The previous blanket protected-area veto persists for **other**
  sources. The mode is a controlled A/B diagnostic, not the final
  source-treatment architecture.

**Critical safety:** default OFF, no new directories/branches,
native FITS immutable, manual \`exclude\` / \`restore_auto\` edits retain
their usual priority, old preliminary isophotes remain inactive,
and the SCIENCE HOLD is unchanged.

### Reproduction and acceptance targets

For the same frozen Leo P G-band dataset (DR11 south \`1554p180\`,
source catalogue SHA-256
\`8e8e3f75af868f3ef7cf4c4c9d7ca28210cbca87ef27c98351a6ab237010562c\`):

1. With checkbox **OFF**, run Prepare → expect **493 / 114** and no
   orange pixels; save a \`Copy Scene QA\` snapshot.
2. Enable checkbox and run Prepare again → expect **5 orange**
   previews, **498 / 109**, same catalogue identity and protected
   region. Inspect source IDs and any overlap with real target light.
3. Recompute four isophote trials and the SGA radial endpoint,
   comparing native component IoU, q/PA, sky, masked pixels and
   leakage warnings against baseline. No new values are pre-claimed.
4. Toggle OFF and re-prepare: the exact baseline source count and
   pixel mask must be recovered; original source FITS and project
   draft digest/revision remain invariant.

The resulting radial change, if any, is **sensitivity to this
particular geometric masking policy**, not validated decontamination.


---

## Seventh-pass Leo P independent MUSE contamination crossmatch (2026-11-11)

Reference: Evans et al. (2019), A&A 622 A129, DOI https://doi.org/10.1051/0004-6361/201834145; Tables 1-3. User QA snapshot 2026-11-10T22:19:04Z, DR11 South 1554p180 G, Tractor SHA-256 8e8e3f75af868f3ef7cf4c4c9d7ca28210cbca87ef27c98351a6ab237010562c. All results READ-ONLY.

- 19/20 spectroscopically confirmed MUSE background galaxy positions match within 1.5 arcsec to 17 distinct Tractor OBJIDs, all unmasked: 1793, 1794, 1799, 1804, 1827, 1839, 1841, 1845, 1853, 1888, 1976, 1978, 2010, 2084, 2110, 2139, 2151. 17/20 lie within 1 arcsec (16 distinct OBJIDs); MSB 172 has nearest model at 5.94 arcsec, not a credible association. Positional matches must be quality-reviewed before physical association.
- MUSE galaxy pairs MSB 19/37 -> Tractor OBJID 1794 and MSB 26/30 -> OBJID 1845 are many-to-one blends, not two accurately modelled profiles.
- Three seemingly foreground stars from MUSE Table 2 associate with OBJID 1870 (0.17 arcsec), 1876 (0.10 arcsec), 1900 (0.38 arcsec). Only the first two are covered by the opt-in Gaia preview. OBJID 1900 is a high-priority non-Gaia foreground candidate.
- Protect known target light: 16/32 Table 1 candidate/confirmed hot-star positions match 10 distinct Tractor components at 1.5 arcsec; all 5 Table 2 cool/AGB/candidate sources match distinct protected components. Confirmed carbon stars MSB 10 -> 1939 and MSB 8 -> 1991 have model g-r approx. 1.92 and 1.95. LP26 O/HII source matches protected OBJID 2076 within 0.41 arcsec. Color or Tractor TYPE must not trigger blanket exclusion.
- Within 60 arcsec: 144 centres = 5 Gaia-preview masks, 97 protected overlaps, 42 unfitted/unscreened; 139 unmasked. Within the current SGA-supported ellipse: 120 centres = 5 Gaia, 1 standard proposed exclusion, 78 screened protected overlaps, 36 previously omitted unfitted/unscreened; 114 unmasked. These are not all contaminants.
- Current SC1 isophote remains LOW CONFIDENCE; selected 2.20/4px component decreased from 45,974 to 34,057 pixels after opt-in Gaia previews (-25.9 percent), so masking induces substantial sensitivity and must be validated with identical inputs, PSF/sky/residual checks.
- Implemented in the existing Draft PR #84: read-only Source Evidence filters TARGET_PROTECTION_ALL, TARGET_UNSCREENED, and full SGA unmasked counts. Do not automatically mask non-Gaia candidates; introduce a general provenance-rich literature crossmatch and research-confirmed source-mask proposal workflow, with many-to-one blend handling and HST/MUSE validation.

No new branch/directory, FITS modification, active mask promotion or SCIENCE adoption.


---

## Eighth-pass scalable conservative source-policy pilot — 2026-11-11

**Scientific objective:** enable unattended high-precision contamination selection when thousands of low-surface-brightness galaxies lack object-specific MUSE/HST/literature analysis. NO per-object literature or network lookup is required by this algorithm. The current state is a conservative **QA-only candidate exclusion**, NOT a calibrated classifier or accepted galaxy photometry.

Policy `red_extended_v1` is implemented as pure deterministic source evidence in the existing `scene_source_inventory.py` and consumed by `conservative_catalogue_mask` within `scene_contaminants.py`. The `red_extended_preview=True` option is explicit, default False, and source-/band-bound. It evaluates native Tractor sources that intersect the target-protection area and that already meet same-band positive flux/SNR screening. This is additional to, independent of, the Gaia astrometry mask preview.

**High-specificity conjunctive criteria (all must pass):**

- Native TYPE in EXP/DEV/SER; never infer physical identity from a single Tractor type. Do not classify PSF or REX models as secure interlopers by color alone.
- No Gaia G2/G3 or SGA L4 reference, no Gaia duplicate, and both FITBITS and MASKBITS exactly zero. Reject strong catalog warning / forced fit regions.
- Native Tractor shape radius between 0.65 and 4.0 arcsec; strongly extended model but not the huge host target component.
- Positive independently fitted g/r/z Tractor fluxes and inverse variances, with SNR g >= 10, r >= 10, z >= 7 (diagonal-only, not PSF-matched); model g-r in [0.80,2.50] and r-z in [0.55,2.00]. Colors are not dust-corrected; a red stellar population alone is not contamination proof.
- At least 2.5 arcsec from the nearest other native Tractor source, measured with a one-time KD-tree; at least max(12 arcsec, 0.2 times target protection radius) from the galaxy anchor. Both are conservative blend/core vetoes, not membership radii.

The source policy returns deterministic evidence metrics (`g_minus_r`, `r_minus_z`, size, target offset, nearest neighbor and three band SNRs), and the source registry marks `RED_EXTENDED_MASK_PREVIEW` with `MULTIBAND_RED_EXTENDED_ISOLATED_HYPOTHESIS`, `AUTO_CLASSIFIER_PURITY_NOT_VALIDATED`. `physical_membership=UNDETERMINED` and `photometry_action=NO_NEW_ACTION_READ_ONLY` remain unchanged. The independent pixel Boolean preview layer and color-coded cyan overlay are used in provisional isophote A/B and radius-sweep, honor existing manual `restore_auto` operations, do not mutate native FITS or project edits, and reset on dataset/source changes.

**Frozen Leo P QA-snapshot offline predicate check:** 6 selected OBJIDs 1674, 1794, 1804, 2151, 2185, 2290. Three sources (1794, 1804, 2151) independently match spectroscopically confirmed background galaxies from Evans et al. (2019, A&A 622 A129); 3 are still unconfirmed. ZERO selected among the earlier listed known member/hot-star, carbon-star, SGA-core and H II Tractor control IDs. This is **not** a measured 100% purity; there is no representative galaxy validation set. Recall of known spectroscopic background galaxies is deliberately low, especially for blue, PSF and REX-type or low-SNR systems, to protect true galaxy light.

Foreground star 1900 without Gaia information is deliberately **not** automatically selected on red color: known carbon stars Leo P (OBJIDs 1939/1991) also have red Tractor SEDs. This is a demonstrated precision/recall tradeoff, not a bug to fix with unconditional PSF red masking.

**Next release gate for actual mass science:** benchmark across diverse dwarf morphologies, environments and distances, with independent resolved/spectroscopic truth sets; measure per-source precision and recall by projected radius, color and morphology; quantify target-galaxy flux loss, radial-profile bias and false-positive source masking with star/galaxy injection-recovery; check spatial sky, PSF wing residuals and per-band aperture consistency. Revise threshold policy only on training objects and verify on HELD-OUT objects. Until this succeeds, keep opt-in QA mode and SCIENCE HOLD. Do not claim all secure contaminants found, do not auto-adopt any mask or use color as a posterior membership probability.


---

## Ninth-pass verified branch state and WCS-safe red-source QA — 2026-11-11

**Verified before changes:** Draft PR #84 against \`main\`, same
\`feature/scene-wcs-editor-sc1\` branch and HEAD
\`5b9f6e665f955335f804c847134ec87aba5ed27c\`.
Both workflows **SUCCESS**, including Ubuntu + Windows × Python
3.10/3.12, GUI-import, and offline Windows PySide6 screenshot.
This checks software regression, not source-classification purity.

**Independent historical Leo P metadata re-evaluation:** loaded the
unchanged uploaded \`Copy Scene QA\` from
\`2026-11-10T22:19:04.846138Z\`, using all 811 records and the
FITS-native WCS header
(\`PC1_1=-7.2777777777778e-05 deg/pix\`,
\`PC2_2=+7.2777777777778e-05 deg/pix\`,
\`CRVAL2=18 deg\`). Angular scales in the east/north
tangent plane are approximately \`0.249177″/px\` in native x,
\`0.262000″/px\` in native y; a scalar \`0.262″/px\` is
not exact for isolation or target-offset vetoes.

Recomputed the exact \`red_extended_v1\` TYPE, bitfield, size,
g/r/z Tractor SNR, g−r, r−z, centre-distance and nearest-source
cuts from the exported catalogue with a sky-plane KD-tree.
The **six** red/extended preview candidates are unchanged under
the corrected WCS-neighbour test:

| OBJID | Tractor TYPE | g−r | r−z | Nearest neighbour (arcsec) |
| --- | --- | ---: | ---: | ---: |
| 1674 | EXP | 0.93 | 0.87 | 3.41 |
| 1794 | EXP | 1.13 | 0.67 | 3.57 |
| 1804 | DEV | 1.62 | 0.84 | 3.57 |
| 2151 | SER | 1.72 | 1.08 | 7.98 |
| 2185 | SER | 0.93 | 0.64 | 16.86 |
| 2290 | DEV | 1.67 | 1.21 | 3.39 |

Independently matched MUSE background galaxies among these six:
1794, 1804, 2151. Three other IDs do **not** have such
independent confirmation. Of 17 distinct Tractor components
associated with 19/20 published MUSE background positions within
1.5″, 3 are recovered (**17.6% reference-set recall**), by design
low for purity-oriented screening. None of 15 distinct native
Tractor components positionally associated with published hot
stellar and cool/AGB target controls were selected (**0/15
observed false positives in this nonrepresentative control set**).
The three unverified selections prevent any claim of 100% purity.
No true population-level rates, photometric light-loss correction
or calibrated probability can be inferred from this single dwarf.

**Corrective engineering in same PR:**

1. Replace scalar pixel-scale neighbour separation with the
   native WCS PC/CD east/north tangent-plane transformation and
   a KD-tree. Use the same WCS for galaxy-centre offset. Red
   screening without valid native WCS now fails closed.
2. Missing native FITBITS/MASKBITS and NaN nearest-neighbour
   distances now fail the candidate gate. If *no* other source
   exists (\`nearest=+inf\`), source evidence explicitly records
   \`no_catalogue_neighbour\` and a JSON-safe null nearest value.
3. Include \`red_extended_preview\` and classifier version
   \`red_extended_v1\` in radius-sweep \`fixed_settings\`.
   Previously the worker could calculate red-mask scenarios but
   the GUI would reject their result as stale, because the
   policy switch was missing from this input fingerprint.
4. Regression coverage: anisotropic WCS vs scalar-distance false
   isolation, JSON-safe evidence, missing/NaN-fit veto, full
   radius-sweep GUI acceptance and unchanged drafts/native arrays.

**Limitations:** the 811 QA records are sufficient for the
specified inner-field pilot, but the original complete 5363-row
Tractor FITS and native IMAGE/INVVAR pixel arrays are not present
in this export. The positional filtering and WCS metadata are
independently recalculated, **not** the physical mask pixel count or
PSF/sky-corrected isophotes. The next on-device Leo P A/B QA
must verify that the original full Tractor catalogue yields the
same six identities under the enabled policy and record the
actual pixel masks, per-band leakage, sky and contour changes.

**Science status remains HOLD.** The red/extended classifier is
a conservative, opt-in **non-adopted diagnostic**, not a
calibrated unattended mass-photometry mask policy.


---

## Tenth-pass actual Leo P GUI source-rejection audit (2026-11-11)

User supplied a new real GUI screenshot and **complete** `Copy Scene QA` JSON generated `2026-11-10T23:26:04.463963Z` for native DR11 South `1554p180` g-band. Catalog SHA-256 remains `8e8e3f75af868f3ef7cf4c4c9d7ca28210cbca87ef27c98351a6ab237010562c`. Both experimental modes active. This is genuine on-device execution, not an inferred run.

Confirmed: 504 candidate-masked sources (493 ordinary + 5 Gaia + 6 red extended), 103 screened protected overlaps, 204 unscreened/unfitted catalogue sources (42 within the circular target protection), 421494 masked candidate pixels, 13728 Gaia-preview pixels, 7767 red-preview pixels. Native hard-invalid 385. Native low-confidence isophote 34054 px at contrast 2.20, smoothing 4 pix; q=0.36 and centre drift 64.1 pix; SGA radial endpoint ~60.2 arcsec close to the 60-arcsec protection radius, so no scientific outer-isophote adoption.

Inside the circular protection region 135 source centres remain unmasked (93 screened overlaps + 42 unscreened/unfitted); 10 more screened protected-overlap source footprints enter the protection region from outside it. Within the SGA-proposal ellipse 110 source centres remain unmasked (76 screened overlaps + 34 unscreened). **These are not 135/110 real contaminants**; the target galaxy produces many Tractor stellar/nebular components.

External MUSE control crossmatch is *held for validation* and not used by the runtime mask policy. Of 17 Tractor OBJIDs associated with 19/20 published background spectra, three are selected by `red_extended_v1` (1794, 1804, 2151). Unmasked background examples: 1793 (EXP, g-r ~0.75, all other strict gates pass), 1845 (SER, blue), 1853 (PSF), 1976 (DEV shape_r ~11.4 arcsec and blue). Likely non-Gaia foreground star OBJID 1900 remains unmasked, but target carbon stars 1939 and 1991 have comparable red point-like fits; unconditional red-PSF masking would incur false positives.

**Technical follow-up:** added a bounded, per-source `red_extended_v1` first-failure audit for *all current-band catalogue-screened protected overlaps*, including non-selected cases. Added separate refusal metadata for unscreened/unfitted source centres within target protection. The exported source inventory retains policy diagnostics with each source: `red_extended_rule_evidence`, `red_extended_first_failure`, `red_extended_near_miss`. GUI `Source evidence` now filters `RED_POLICY_UNMASKED` and `RED_POLICY_NEAR_MISS`, showing rejections without inventing membership. A near miss means *all preceding policy gates passed*, but g-r is between .70 and .80 and r-z is in-range; it is **not** automatically excluded. Based on the frozen Leo P JSON, OBJID 1793 is the sole near miss. Existing red/Gaia/ordinary geometric masks, project draft, original FITS and science hold remain unchanged.

**Why this matters:** weakening fixed color or morphology thresholds to chase individual spectroscopic known contaminants would overfit one dwarf. Instead, first quantify rejection frequencies and inspect the scatter for multiple independent dwarfs; use held-out MUSE/HST cases to evaluate source false positives and aperture/sky flux loss, before migrating any candidate tier into unattended science masks. Additional measurable native Tractor fields such as model-fit DCHISQ, RCHISQ, FRACFLUX/FRACMASKED, and pixel-level segmentation residuals should be tested in a separate, independently validated classifier; the current source registry does not yet certify their discriminating power.

Regressions: new tests exercise exact red rejection trace for EXP near-miss, red target PSF, low-band-SNR unscreened source, positive red mask and GUI filters. No additional branches or source directories created.


---

## Eleventh-pass Leo P native protection-radius sweep (2026-11-11)

The user's complete GUI `Copy Scene QA` generated `2026-11-10T23:51:27.622925Z` reports a **successful read-only 5/5 experiment** in native g-band DR11 South `1554p180`, with SGA parent geometry ID 5065058, contrast 2.2, Gaussian smoothing 4 px, 2-pixel bridge guard, optional inner Gaia and red-extended candidates enabled. Same native input catalog SHA-256 `8e8e3f75af868f3ef7cf4c4c9d7ca28210cbca87ef27c98351a6ab237010562c`. All project edits remain empty and unadopted; the native FITS and project target radius were not modified.

| Protection radius (arcsec) | SGA-geometry radial estimate (arcsec) | Masked source footprints | Candidate-excluded pixels | Connected component pixels | q | Component centroid shift (pixels) | Component IoU relative to 60 arcsec |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 40 | 50.4070 | 546 | 447671 | 27108 | 0.47178 | 42.87 | 0.79561 |
| 50 | 53.6591 | 526 | 436118 | 31477 | 0.39167 | 56.03 | 0.92433 |
| 60 | 60.1632 | 504 | 421494 | 34054 | 0.36042 | 64.08 | 1.0 |
| 70 | 61.7893 | 490 | 410377 | 34182 | 0.36029 | 64.14 | 0.99555 |
| 80 | 61.7893 | 473 | 397881 | 34182 | 0.36029 | 64.14 | 0.99555 |

**Interpretation:** the preliminary radius increases with the *provisional protection size* for 40–60 arcsec, so those scenarios demonstrate substantial mask dependence; the plateau in the scalar radial metric at 70–80 arcsec is encouraging but **not** a validated galaxy boundary or proof of identical pixel-wise components. Compared to 60 arcsec, 40 arcsec adds 42 source candidates to the excluded class, including known spectroscopic background target OBJID 1793; the same policy also risks suppressing real resolved target stellar light. Between 60 and 70 arcsec, 14 catalogue source identities become protected overlaps, and the native connected-component IoU is 0.99555, but the contour remains unanchored and low-confidence. The empirical constant local sky level `2.1435495e-05 nMgy/pixel`, smoothed RMS `2.717869e-4 nMgy/pixel`, and 487180 training pixels were *numerically unchanged* across these five scenarios, despite the implementation refitting sky per scenario. No independent spatial sky or PSF-wing experiment was performed.

All five scenarios report `LOW_CONFIDENCE_REQUIRES_EXPLICIT_OVERRIDE`, with centroid shifts from 42.9 to 64.1 px, spatial-trial minimum IoUs roughly 0.72–0.77 and contour-area sensitivity ratios 1.45–1.57. The slope of the radial result is not a surface-brightness profile and its convergence is not science adoption. Geometry in SGA parent viewer is a provisional prior, not an independently calibrated 26 mag/arcsec^2 isophote.

**Radial bin resolution (checked against the implementation):** the radial endpoint is selected from discrete annular outer edges, with effective step `max(max(3, 1.5*smooth_sigma_pix)*pixel_scale, max_a/64)`. In the actual Leo P run, the consecutive sweep deltas are exactly `3.252066311, 6.504132622, 1.626033156, 0` arcsec, i.e. **2, 4, 1, 0 bins** at a radial sampling of **1.626033156 arcsec**. The 70/80 equality therefore only establishes the **same radial bin**, not a better-than-arcsecond measured edge. A new `radial_step_arcsec` provenance field is saved on the SGA radial seed and `sga_radial_bin_step_arcsec` on every sweep scenario; these are sampling resolution metadata, **not an error bar**. Do not publish the unrounded radial extent as a physical radius.

**Instrumented follow-up:** previous-radius pairwise diagnostics were added to `run_scene_protection_radius_sweep` (next source revision), including pixel-wise component IoU, mask Jaccard, XOR pixel count, radial-endpoint delta and source-disposition transitions for each adjacent scenario. This directly tests 70-vs-80 pixel identity on the *next* actual application run; equality of summary values alone is insufficient. Regression tests include a failed contour and an unchanged contour despite differing masks, with no science or file writes. In a future physically meaningful photometric stage, compare sky models and PSF wings independently; test a fixed sky across mask scenarios and model source flux retention before choosing a protection aperture for unattended processing.


---

## Twelfth-pass real Leo P adjacent-radius pixel comparison (2026-11-11)

User reran real GUI `Copy Scene QA` at `2026-11-11T00:10:46.868295Z` (native DR11 South brick `1554p180`, g band, catalog SHA-256 `8e8e3f75af868f3ef7cf4c4c9d7ca28210cbca87ef27c98351a6ab237010562c`). The new adjacent-radius instrumentation has now been independently exercised on real Leo P inputs. All 5/5 scenarios completed; all remain `ISOPHOTE_PROPOSED_NOT_ADOPTED` and `LOW_CONFIDENCE_REQUIRES_EXPLICIT_OVERRIDE`.

| Protection change | Mask XOR pixels | Mask Jaccard | Connected component IoU | SGA radial endpoint delta |
|---|---:|---:|---:|---:|
| 40→50 arcsec | 11553 | 0.974193 | 0.860727 | +3.252 arcsec (2 bins) |
| 50→60 arcsec | 14624 | 0.966468 | 0.924326 | +6.504 arcsec (4 bins) |
| 60→70 arcsec | 11117 | 0.973625 | 0.995555 | +1.626 arcsec (1 bin) |
| 70→80 arcsec | **12496** | **0.969550** | **1.000000** | **0 arcsec (0 bins)** |

**Confirmed pixel-wise plateau:** for the last pair the same **34,182 native connected-component pixels** were selected (IoU=1.0) despite 12,496 changed contamination-mask pixels and 17 source identities newly categorised as protected overlaps. Both have moment q=0.3602860567 and 64.1351-pixel centroid offset; the SGA-guided radial endpoint is 61.7893 arcsec and the explicitly exported radial annulus step is 1.62603315556 arcsec. This is evidence of segmentation invariance to this *particular* mask enlargement, not a validated isophote or evidence that masks themselves are unchanged. The 60→70 component IoU was 0.99555, contrasting with marked changes at 40→60. The SGA guiding geometry remains parent prior q=0.53187, PA=150.46 degrees; the moment q=0.36029 is a separate masked-contour diagnostic, not an independently fitted isophote shape.

The experiment retains Gaia interior preview **ON**, strict red extended preview **ON**, SGA geometry **ON**, SGA protection **OFF**, contrast 2.2, smoothing 4 pixels, bridge guard 2 px, with no target edits and no science adoption. All scenario local skies remain exactly `2.1435495e-05 nMgy/pixel` and smoothed RMS `0.0002717869 nMgy/pixel`; these are provisional empirical estimates, not a new sky validation. Four threshold/smoothing trials remain sensitive (minimum IoU ≈0.723 at 70/80; area ratio ≈1.57) and centroid drift >64 px carries `UNANCHORED_CENTROID_DRIFT`.

**Next discriminating experiment:** hold g-band image, parent SGA geometry, threshold, smoothing, centre and protection-radius sweep constant, but run the same diagnostic with the experimental red-extended masks OFF (Gaia ON), then with BOTH experimental masks OFF. Compare the **adjacent-radius 70→80 component IoU**, radial endpoints, retained-source IDs and source/galaxy flux-loss risks separately. A later independent PSF-wing and spatial-sky comparison is necessary before choosing an unattended mask policy. Continue SCIENCE HOLD.


---

## Thirteenth-pass Leo P foreground-mask policy ablation (2026-11-11)

Two full, real on-device g-band QA exports were supplied after a previously verified both-preview-ON sweep: `2026-11-11T00:17:56.857826Z` (Gaia interior ON / red extended OFF) and `2026-11-11T00:18:46.524464Z` (both OFF); comparison run `2026-11-11T00:10:46.868295Z` had both ON. All share project `Leo_P` ID `0876bbc4-7405-481b-88a8-0e3486c7be98`, identical source SHA-256 `8e8e3f75af868f3ef7cf4c4c9d7ca28210cbca87ef27c98351a6ab237010562c`, draft digest `2727f8959a9bcd15338f5a3c6abbf12eafe569e34450a14c4c59643cfd502619`, SGA parent reference 5065058 and SHA, native g DR11 South brick 1554p180, contrast=2.2, smoothing=4 pixels, bridge guard=2, SGA geometry used and core protection disabled. Each completed all five protection radii (40, 50, 60, 70, 80 arcsec), with zero adopted mask edits.

| Radius | Experimental previews | Excluded catalogue pixels | Connected component px | Moment q | Centroid drift px | Radial bin endpoint arcsec | Minimum threshold/smoothing IoU |
|---:|---|---:|---:|---:|---:|---:|---:|
| 60 | Gaia+red | 421494 | 34054 | .360421 | 64.075 | 60.163 | .72566 |
| 60 | Gaia only | 413770 | 34057 | .360422 | 64.077 | 63.415 | .72411 |
| 60 | Neither | 400042 | 45974 | .541076 | 39.859 | 63.415 | .73477 |
| 70 | Gaia+red | 410377 | 34182 | .360286 | 64.135 | 61.789 | .72348 |
| 70 | Gaia only | 402626 | 34185 | .360287 | 64.135 | 63.415 | .72194 |
| 70 | Neither | 388230 | 46102 | .540712 | 39.894 | 63.415 | .73309 |
| 80 | Gaia+red | 397881 | 34182 | .360286 | 64.135 | 61.789 | .72348 |
| 80 | Gaia only | 390130 | 34185 | .360287 | 64.135 | 63.415 | .72194 |
| 80 | Neither | 371674 | 46102 | .540712 | 39.894 | 63.415 | .73309 |

**Key result (fixed-mask-policy stability is NOT cross-policy stability):** for every 40–80 arcsec sweep radius the Gaia-preview-enabled connected component is exactly **11,917 pixels smaller** than the corresponding no-preview component; at radius 70/80 this is 25.85% of the no-preview area, accompanied by moment q changing from .5407 to .3603 and centroid offset increasing 39.89→64.14 pixels. At 70 arcsec, even the theoretically *maximum* possible pixel IoU between no-preview (46,102 px) and both-preview (34,182 px) components is `34182/46102 = 0.741443`, and the actual IoU may be smaller: pixel arrays are intentionally omitted from QA. **Do not claim** precisely 11,917 galaxy-light pixels were removed or that the mask pixels equal the component-area difference; there are foreground sources and topology effects. The red-extended preview changes component *area* only **3 pixels** per radius after Gaia masks, yet changes the *SGA radial endpoint* by two 1.626 arcsec bins at 60 arcsec, one bin at 70/80 arcsec. This independently demonstrates important radial-diagnostic sensitivity despite nearly unchanged component area.

All three policies preserve **within-policy** pixelwise component IoU=1.0 for 70→80 arcsec even while mask pixels continue to change: 12,496 XOR pixels for Gaia+red, 12,496 for Gaia only, and 16,556 for neither. Hence the previously established radial-protection plateau is reproducible *within each mask policy* but cannot validate the mask policy or a scientific galaxy edge. Sky diagnostics (local sky 2.1435495e-05 nMgy/pixel and smoothed RMS 2.717869e-4) are unchanged, with no spatial background or independent PSF model; minimum threshold/smoothing IoUs remain 0.722–0.733 at r=70/80. All candidates remain LOW_CONFIDENCE, read-only and unadopted.

Five Gaia interior candidate OBJIDs are `1763,1870,1876,2034,2035`. Source inventory records significant Gaia astrometry: e.g. OBJID 1763 parallax significance 5.44 and |PMdec| significance 96.98; 1870 parallax 6.38 and |PMra| 147.40; 1876 parallax 5.24 and |PMdec| 77.64; 2034 proper motion |PMra| significance 11.95 even though parallax insignificant; **2035** G≈15.03, parallax significance ≈30.7, bright g-band Tractor model flux ≈607.4 nMgy, angular distance from galaxy centre ≈19.87 arcsec and preview radius ≈52.7 native pixels. These are credible **foreground-source** indicators, *not proof that the whole extended mask can be removed without losing diffuse Leo P light*. Do not infer real-galaxy flux loss from pixel area alone.

**Next useful implementation/test:** obtain exact *cross-policy* pixel overlap by running alternative policies on a fixed original native FITS and fixed radius (70 arcsec), including per-source Gaia mask **leave-one-out ablations** (starting with bright OBJID 2035) and a source/target-light flux budget. Maintain separate contamination removal (including PSF-wing models), masked-target-light recovery and astrophysical isophote fitting; validate with injection/recovery and independent background/PSF assessments before unattended mass photometry. Preserve `SCIENCE HOLD`.


---

## Fourteenth-pass fixed-radius mask-policy ablation QA (2026-11-11)

**Motivation:** Leo P g-band exhibits a 34,185 vs 46,102 connected-component pixel difference at 70" when Gaia-preview footprints are enabled versus disabled. The same component can be unchanged between protection radii 70" and 80" **within** each policy, without showing that the policy, galaxy-light retention or radial photometry is valid.

**Implementation:** Added run_scene_mask_policy_ablation in the existing src/lsbpipe/gui/scene_editor.py, with the GUI button **Mask policy ablation · fixed 70″**, a status line and a new mask_policy_ablation section in the existing Copy Scene QA JSON export. The worker always holds fixed: project, original native FITS, WCS, band, draft digest, centre, SGA priors, image contrast, smoothing, bridge guard and source-protection radius (70 arcsec). Four simultaneous cases are: (1) baseline normal Tractor masks, (2) baseline + experimental Gaia interior masks, (3) baseline + red extended experimental masks, (4) baseline + Gaia + red. The worker then derives all native currently eligible Gaia preview OBJIDs from its actual source inventory and runs a **Gaia-star leave-one-out** for every candidate with the other sources retained under case (4). No hard-coded Leo P IDs, external per-object literature or network requests.

The optional gaia_inner_skip_objids parameter is only valid for explicit Gaia previews, with a strict bound of <=12 unique source identities in the ablation. A left-out star reverts to the original protected-overlap policy for that source only, not to an adopted astrophysical membership. If the inventory is incomplete, catalog SHA changes, or too many/ambiguous Gaia candidates exist, the worker fails closed. All four core policies plus every individual intervention share the original image and native pixel grid.

The output gives **exact** Boolean mask XOR, Jaccard, connected-component pixel IoU/intersection and mask overlap against the same baseline component; source-disposition changes, component area, moment q and centroid drift; SGA radial-bin endpoint and optional local sky/stability diagnostics. Failed contours remain explicit result rows. Actual mask/contour arrays stay in worker memory and are NOT exported. A mask overlapping the no-experimental-mask component is not automatically the amount of *galaxy flux* lost: those pixels may contain unrelated star/foreground emission.

This is READ-ONLY and NOT SCIENCE VALIDATED: no native FITS, project, mask draft or photometric product is changed, even in a failed scenario. Sky is provisionally refitted per policy; no independent sky model, PSF-wing model, source subtraction, physical membership inference, calibrated galaxy-light recovery, isophotal aperture measurement or flux uncertainty estimate is made.

**Regressions:** source-specific Gaia skip preserves other source masks; four-policy scenario masks and two synthetic Gaia interventions check direct IoUs; one missing contour remains in the result; Qt button launches asynchronous worker; Copy Scene QA reports the result; control changes invalidate the diagnostic; existing files and drafts are unchanged. All in the existing Draft PR #84 branch. No new branches or directories.

**Expected next user run:** open Leo_P, load current original g FITS, choose SGA2025-PARENT Leo P initial geometry, SGA core protection OFF, contrast=2.2, smoothing=4 pixels, bridge guard ON. Press **Mask policy ablation · fixed 70″** and then **Copy Scene QA**. Previous real radius=70 Gaia-only scenario had six Gaia candidate centres, so 4+6=10 scenarios are expected *if the same catalogue/inputs are used*. This must be verified from the new exported report. Current opt-in preview checkbox states do not affect which four policies the experiment tests.

**Before production unattended aperture photometry:** independently validate PSF/wing handling and galaxy light hidden behind contaminants; spatial sky and correlated-noise uncertainty; reproducible multi-band common ICRS aperture and valid-pixel fraction; injection/recovery and held-out source truth tests; explicit QC/uncertainty flags to skip failed objects rather than accepting unstable photometry. The current SC1 SCIENCE HOLD remains.
