# Desidera performance audit — v0.2.1.dev2 (2026-11-13)

Status: **development candidate in existing Draft PR #84**,
`feature/scene-wcs-editor-sc1`. No merge, new branch or directory,
and no change to original science FITS. No new astronomical result is
accepted or calibrated by a speed optimization.

## Evidence classification

**CODE-PROVEN:** mathematical operator identities, unnecessary
allocations/algorithmic asymptotics, number of avoided operations,
input immutability, the exact source functions changed.

**CI-VERIFIED:** assertions on synthetic arrays against independent
Photutils reference operations and GitHub Actions test execution
(only after the exact-HEAD CI has completed).

**UNMEASURED:** absolute time improvement for real Leo P FITS,
DR11 brick-level throughput, Windows GUI frame latency, peak process
RSS on the researcher's machine, PSF wall-time, network/disk effects,
and end-to-end Stage01–04 performance. These require the
hardware-specific observational benchmark, not inferred timings.

## Source-observed performance inventory

| Priority | Location and evidence | Mechanism | Action / status |
| --- | --- | --- | --- |
| P0 | `photometry/aperture_systematics.py`: exact blank-ring loop | Every accepted or rejected candidate rerasterized every fractional Photutils ring, despite native integer translations after M23 phase matching | **FIXED:** reuse one exact ring mask per radius, translate only its bounding box in native pixels; independent Photutils and edge-slice tests |
| P0 | `photometry/native_fits.py`: immutable loader | It previously built an intermediate `NativeBandFrame` before hashing and a second complete frame after hashing, duplicating full-grid conversions, validity layers and WCS checks | **FIXED:** exactly one frame construction after FITS headers and full-file SHA validation, final path/stat guard after materialization |
| P0 | `photometry/aperture_growth.py` / `aperture_flux.py` | Growth code copied Scene source/manual masks once but each cumulative aperture recopied both entire native arrays | **FIXED:** freeze two masks once in immutable byte-backed buffers; `_mask_for` reuses only provably bytes-rooted arrays, still snapshots all ordinary GUI/mutable inputs |
| P1 | `background.py`: M33 spatial convex hull | Full cloud of train pixels was copied for Qhull | Already **FIXED** in scientific M33: exact min/max x per occupied image row; retained and tested |
| P1 | `photometry/aperture_growth.py`, `aperture_flux.py` | Repeated `frame.wcs` defensive deep copies and 16-point sky projection/validation at each radial geometry/hypothesis | **OPEN:** assess immutable WCS/projected geometry cache keyed by complete frame WCS, center, radii, q, PA and code version; do not drop existing affine-error gate |
| P1 | `gui/scene_editor.py`: `load_scene_original_cutout` | Full source FITS SHA passes and defensive raster copies every load; additional FITS header opens | **OPEN:** benchmark header/cache/IO policy with exact content-change detection; do not replace full SHA with size+mtime as supposed authenticity |
| P1 | `gui/scene_editor.py`: `_rerender_image` | Whole-array display stretch and QPixmap are recomputed on display-mode changes | **OPEN:** immutable display-only raster cache with bounded bytes and LRU; keep source science data separate |
| P1 | `photometry/scene_contaminants.py` | Repeated per-source local cutout rasterization; variants repeat source mask generation; dense foreground may create high cumulative footprint area | **OPEN:** instrument source count, per-source patch area, variants and wall time; don't cache across changed WCS/Tractor/target protection |
| P1 | `background.py`: `_automatic_scene_sky_pixels` | Native 2-D distance transforms, candidate tiles and repeated fixed-grid sensitivity variants require full-grid passes | **OPEN:** compare masks/tile geometry reuse per frozen source and draft ID; never reuse under changed source or manual masks |
| P1 | `photometry/psf_empirical.py` | Iterative ePSF building, per-star cutouts and multiple spatial/holdout QA modes | **OPEN:** separate default PSF-core profile from explicit advanced diagnostics and instrument durations; cannot discard QA as an 'optimization' |
| P2 | `photometry/contracts.py`: immutable copy policy | Three float64 layers + one uint32 mask + multiple Boolean reason layers are held per native grid; high transient memory pressure | **OPEN:** peak RSS and per-stage allocation profile before contemplating less precision; preserve native nMgy/IVAR fidelity |
| P2 | `gui/photometry_workspace.py`, `gui/scene_editor.py` | Rebinding, JSON QA serialization, scene overlays and workspace redraw can block UI for large objects | **OPEN:** trace main-thread blocked time and p50/p95 UI event latency; retain worker generation/identity checks |
| P2 | `surveys/dataset_integrity.py`, catalogue and project readers | Deep FITS validation is inherently I/O-heavy; inventory and reuse/reopen may repeat metadata traversal | **OPEN:** instrument file accesses and root sizes; prefer provenance-aware incremental *read* caching over skipping safety checks |

These are source-identified mechanisms and prioritization judgments,
not demonstrated relative runtime contributions on a real galaxy.
The initial v0.2.1.dev0 bottlenecks were chosen because the scientific
operator can be preserved exactly without broad architectural change.

## Optimization P0.1 — exact same-phase ring-template translation

Let `w_i(p; x,y)` denote Photutils exact native overlap weights for
ring `i`, located at `(x,y)`. If the blank centre is a native integer
translation `(k,l)` from the fractional-phase-matched target centre,
the pixel-space area overlap satisfies

`w_i(p; x+k,y+l) = w_i(p-(k,l); x,y)`.

The science image is never interpolated or shifted. The new implementation
constructs each target ring's exact Photutils `ApertureMask` just once
and translates its bounding-box slices by integer offsets for all
eligible/attempted blank centres. It crops those slices against the
unchanged original FITS grid and retains the original rejection rule:
every positive-weight intersection with invalid/source/manual pixels
invalidates a blank. All conditional sky coefficients and per-pixel
IVAR are still evaluated at the blank's actual native coordinates.

Reference rasterizations for `R` rings and `B` attempted candidates:
up to `R*(1+B)`. Revised exact template rasterizations: `R`,
plus `O(R*B)` cheap bbox and weighted-sum operations that are
necessary regardless. Computational *operator creation* is reduced;
full blank signal/variance summation remains real work. The report
exposes `n_exact_ring_templates_rasterized` and
`n_integer_translated_ring_evaluations`.

Independent tests compare signed blank covariance/means with separately
rasterized original Photutils apertures under a spatially tilted frozen
sky and heteroscedastic IVAR. An additional test checks border-clipped
and off-frame slice identities. No fixed performance-speed threshold
is used as a CI gate.

## Optimization P0.2 — one full NativeBandFrame build

Former behavior:
- `_validate_open_hdus` verifies all four FITS sources.
- construct a `NativeBandFrame` with empty local SHA metadata,
  materializing immutability and validity;
- close FITS, calculate verified full source SHA-256;
- build **another** `NativeBandFrame` with the same native pixel data
  just to update `source_evidence`.

Revised behavior:
- validate HDUs and BUNIT while FITS are read-only;
- hash all four original full files and check initial size/mtime;
- construct **one** detached, immutable `NativeBandFrame` with final
  SHA-256 evidence while source handles are open;
- recheck safe paths and initial file metadata *after* array copying,
  then close all handles.

This preserves the original BUNIT, SHA, HDU, WCS, immutable-array,
cooperative-cancellation and no-write guarantees. SHA-256 is not
skipped, weakened to mtime or claimed to be archive authenticity.
Any external change during frame materialization with detectable
metadata difference is rejected. A synthetic regression instruments
`NativeBandFrame` constructor count (exactly one) and triggers a
post-materialization source mutation to prove fail-closed detection.

## Optimization P0.3 — one immutable mask snapshot per growth run

At the beginning of Stage04.2, source and manual masks now receive one
`bytes`-rooted 2D native Boolean snapshot each. The normal
`measure_native_aperture` API still independently checks dimensions
and datatype. The private `_mask_for` has a fast path **only** if
the root backing buffer is truly immutable `bytes` or a readonly
memoryview of `bytes`. Ordinary writeable NumPy arrays, and read-only
views ultimately backed by writable arrays, are copied as before.

For `R` cumulative apertures, repeated complete array copies change
from two per radius to two for the entire growth curve. The
independently integrated annuli, apertures, mask SHA receipts and
measured signed flux retain their same definitions. Tests assert
object identity across all radii and verify that `setflags(write=True)`
cannot change the frozen inputs.

## Reproducible synthetic benchmark

```powershell
python -m pip install -e ".[exact,test]"
python scripts/benchmark_desidera_performance.py --size 512 --repeats 3
```

The script generates **only synthetic native arrays**, performs
four separately timed reference-vs-optimized kernel comparisons and
prints a JSON report:
1. two immutable native-frame materializations vs one;
2. ten full-mask snapshots vs byte-backed reusable snapshots;
3. each blank-ring exact Photutils rasterization vs a translated
   single exact ring template per radius;
4. three dense exponential-profile sky-injection grids vs broadcasted
   coordinate offsets and one reusable elliptical-radius operator.

The third kernel independently checks exact-overlap integrated area
before reporting times. Ratios are specific to one machine/run.
No I/O, real FITS, network, GUI, PSF or background-selection throughput
is measured; the result is NOT a performance SLA and no observed
speedup is asserted before executing it. CI includes a low-resolution
smoke run but has **no flaky wall-time acceptance threshold**.

Suggested future measurement protocol: record machine/CPU/RAM,
Python/NumPy/Astropy/Photutils version, cold/warm file cache,
native FITS sizes, real Leo P g/r/i/z timing per Stage01–04,
`time.perf_counter()` median of 5 repetitions, peak RSS
(`psutil.Process(...).memory_info().rss`), source bytes read,
candidate attempts, user-perceived first-image time and Qt main-thread
p50/p95 responsiveness. Run the same identical native data and masks
before and after any change, and verify fluxes/variances/QA receipt
semantics independently.

## Remaining engineering risks, in order

1. **Peak RAM / 16 GB-class workstation:** native quartets + immutable
   byte-backed masks + source/sky/PSF work buffers + GUI RGB/pixmaps
   can coexist. A 2048-square grid has 4,194,304 pixels; just three
   `float64` arrays plus one `uint32` array occupy approximately
   112 MiB before masks, source evidence, sky variants, Qt display
   and temporary arrays. Concurrency across bands risks excessive
   process RSS; measure rather than guessing a safe band count.
2. **Overhead of repeated WCS/ellipse projection:** safe projected
   ellipse caching may be useful with many frozen radii and hypotheses,
   but the current 0.15-native-pixel affine error test must still
   apply to every relevant geometry and original WCS version.
3. **FITS hash/I/O duplication across different entrypoints:**
   may dominate cold starts with larger bricks. A trusted cache needs
   a security and provenance policy that cannot silently treat a
   modified FITS as the same immutable observation.
4. **Scene redraw and source footprint throughput:** profile
   mask rendering, Qt pixmap conversion, red-extended/Gaia
   ablations and their GUI event scheduling separately from science
   CPU kernels.
5. **Advanced PSF diagnostic cost:** instrument ePSF iterations and
   audit modes, but never turn an unrun or timed-out safety check
   into a passed science holdout.

## Release policy

`0.2.1.dev1` is an **installable development candidate** on Draft
PR #84, NOT a Git tag, merged main release, verified performance SLA
or validated astronomy science release. The same on-disk object-project
formats, legacy `lsbpipe` API, preserved original FITS and scientific
SCIENCE HOLD continue. The final change status is determined only by
exact HEAD GitHub Actions and dedicated observational benchmarks.


## Continued performance audit — v0.2.1.dev1, P35–P38 (2026-11-13)

This iteration **continues** the pre-existing source audit rather than
assuming that earlier benchmark runs represent the full scientific
pipeline. Four additional exact-operator improvements were made in
the same Draft PR #84. None changes source-mask policy, sky-model
acceptance, fitted noise covariance, PSF scientific validation,
original FITS, scene-draft lineage, or the science HOLD.

### P35 — repeated 2×2 native WCS SVD per Tractor row — FIXED

`photometry/scene_contaminants.conservative_catalogue_mask` derives
a native-space outer bounding radius from the smallest singular value
of one local tangent-plane WCS Jacobian:

`r_bbox_pix = source_major_axis_arcsec / sigma_min(J_WCS)`.

The Jacobian was computed once, yet its `np.linalg.svd(...)` was
unnecessarily repeated inside the per-Tractor-source loop. For
`N` valid catalogue candidates this caused up to `N` identical
decompositions. The SVD now runs **once for the entire fixed-WCS
mask build**, ahead of candidate iteration; per-source angular radii,
TAN shape orientation, pixel footprints and protection-overlap
decisions remain unchanged.

A synthetic test instruments SVD calls for 72 sources under
anisotropic/sheared WCS, asserts **one decomposition**, and preserves
the existing no-WCS legacy fallback. The scientific limitation is
unchanged: a single local Jacobian is only a TAN-local affine
approximation, not a measured extended PSF or true celestial
isophotal fit.

### P36 — double factorization inside each sky IRLS step — FIXED

`background._weighted_lstsq` called
`np.linalg.svd(weighted_design, compute_uv=False)` to determine
the condition number, then `np.linalg.lstsq(weighted_design, weighted_y)`
to solve the **same weighted design**. The latter already returns
the singular-value spectrum and matrix rank, so performing an
additional decomposition repeats O(n_samples × p²) work at each
robust Huber/IRLS iteration.

The implementation now obtains coefficients, rank and singular values
from the same `lstsq` call. The original
`condition > max_design_condition`, nonfinite-parameter and rank
guards remain. Regressions assert coefficients and condition numbers
against the independent SVD reference; they deliberately disable
Python-level `np.linalg.svd` in the tested function and exercise
a collinear and an explicit condition-threshold failure.

This does **not** reduce sky training samples or relax rejection
of ill-conditioned tilted planes; it only avoids a duplicate
decomposition. The impact on total Stage03 wall-time has not yet
been directly profiled on a real coadd.

### P37 — redundant full-frame sky-injection coordinate arrays — FIXED

`background._audit_current_project_sky_injection_recovery`
formerly created two full `float64` coordinate grids
`yy, xx = np.indices(shape)`, one circular-radius grid and,
for each of the **three predeclared** exponential-profile scale
lengths, recalculated the same elliptical-radius grid:

`r_ell = hypot(x-x_anchor, (y-y_anchor)/0.65)`.

The method now uses 1-D native x/y offset arrays with NumPy
broadcasting, creates the circular selection grid once, and reuses
a single elliptical-radius grid for all three fixed scale lengths.
The selection remains on original valid, unmasked native pixels and
each injected image remains a separate `float64` copy; no model
selection is changed.

At 2048×2048 pixels, the **two removed persistent coordinate
arrays alone** would occupy 64 MiB of `float64` storage
(2×2048²×8 bytes), before any other arrays or temporary
broadcast operations. This is a calculation from array shapes,
**not** an empirically observed reduction in process peak RSS.

An independent test reconstructs the earlier full-density
`np.indices`/ellipse expression and checks the exact integrated
synthetic flux over the complete native grid and surviving aperture
for three profile scales. The extended benchmark measures both
methods on the same synthetic grid without setting absolute timing
requirements.

### P38 — Stage04.3 blank-scale probes copied the same masks repeatedly — FIXED

`aperture_systematics._scene_masks` defensively copied full-frame
source and manual exclusion masks at **every** call to
`estimate_blank_annulus_covariance`, even when invoked by
`diagnose_aperture_systematics` on an already-snapshotted input.
The fallback attempt at smaller annular scales could thus allocate
another pair of full-grid Booleans for each predeclared scale.

`_scene_masks` now snapshots into the existing immutable,
`bytes`-rooted Boolean representation from Stage04.2.
`_snapshot_native_boolean_mask` reuses that buffer only after
walking the complete NumPy `base` chain and confirming that
the ultimate immutable root really is read-only `bytes`.
A NumPy view with `writeable=False` but ultimately mutable
parent still undergoes a defensive copy.

Regressions check object identity over eight nested probes,
mutation-isolation from a GUI-owned writable array, failure of
`setflags(write=True)` on byte-owned masks, and **bitwise-identical**
blank covariance/QA for mutable versus pre-frozen source masks.
This is a local memory optimization, not a workaround for source
identity and mask provenance gates.

### Updated synthetic kernel benchmark

The existing command remains:

```powershell
python scripts/benchmark_desidera_performance.py --size 512 --repeats 3
```

The benchmark now reports four independent comparisons, adding
`sky_frozen_injection_dense_xy_vs_broadcast_reusable_distance`.
The added comparison verifies integrated synthetic flux independently
before reporting two wall times, the formula-based byte count for
the omitted dense coordinate grids, and the number of elliptical
distance-map recomputations avoided.

**No end-to-end observational timing was executed here.** Neither
an apparent synthetic speedup nor a green CI establishes faster
real-world photometry under high catalogue density, PSF iterations,
FITS I/O, GUI event processing or large image cutouts.

### Remaining performance roadmap after P38

1. **Measure real process RSS and cold/warm FITS time:** isolate file
   decompression, SHA-256, FITS decoding, immutable frame
   construction, Qt RGB images, Stage01 empirical PSF, Stage02
   masks, Stage03 variants and Stage04 radial+blank diagnostics.
   Record Python, Astropy, Photutils, NumPy and OS versions.
2. **Repeated Stage04 WCS transformations:** `frame.wcs` returns a
   defensive deep copy, as required for immutability. Most radial
   measurements repeat 16-boundary-point affine-accuracy checks.
   A bounded cache is a potential future improvement, but must be
   keyed by genuine immutable source WCS/geometry/version evidence
   and NEVER skip a changed-radius affine-error limit.
3. **Whole-image distance-transform reuse:** Stage03 repeats target
   and source guard constructions for six predeclared grid/guard
   sensitivity trials. A protected, immutable distance-field cache
   could reduce these passes, but the masks and guards differ across
   trials and stale-cache risk is scientifically significant.
4. **Dense Tractor footprint rasterization:** repeatedly allocating
   large source patches can dominate time for crowded fields.
   Candidate batching or footprint reuse must retain exact
   target-overlap, SGA/Gaia protection policies and source identity.
5. **GUI responsiveness:** Scene Editor source registry rebuilds a
   QListWidget for every filter/search change; display image
   stretches can also hold up the main Qt thread. Consider a
   debounced, bounded list model and cached display-only renditions,
   but test selection stability and stale-generation protections.
6. **Whole pipeline vs kernels:** independently establish p50/p95
   user-visible latencies, stage throughput and failure rates over
   a fixed real-data corpus; no observational performance SLA is
   currently justified.

**Development status: v0.2.1.dev1 candidate, Draft PR #84.**
No new branches/directories. Keep the original FITS and `main`
unaltered, preserve all science safeguards and SCIENCE HOLD.


## Third performance audit — v0.2.1.dev2, P39–P42 (2026-11-13)

This further investigation follows P35–P38 and the earlier P0
optimization cycle, all on the **same existing feature branch and
Draft PR #84**. It adds four additional optimizations without
replacing original FITS, skipping any scientific gate, modifying
coadd pixel values or silently changing source/sky decisions.

### P39 — repeat Stage03 full-image Euclidean distance transforms across fixed grid trials — FIXED

**Source:** `background._automatic_scene_sky_pixels` is run for
the initial automatic sky selection and six predeclared fixed
sensitivity trials. Every call previously recalculated the full
Euclidean distance transform of native **target** protection and
the combined provisional-source/unresolved-footprint veto. Five
of the six variants change only tile edges, tile offsets or the
sampling policy. The other two increase one guard radius each:
`expanded_target_1p5` and `expanded_sources_1p5`.

**Correction:** `diagnose_current_project_residual_sky` creates a
private dictionary, and passes it to the initial sky selector and
six trials. Only the **first, baseline** target and catalogue
guard-distance threshold masks are retained. At every possible
reuse, both the radius and **all original native Boolean mask
pixels** are compared to the cached source mask snapshot.
If either differs, the EDT is recomputed; altered-radius trials
never evict the baseline. There is no module-global cache, user
project cache, cross-band reuse or lookup based only on mtime,
input object identity, sky flux, or number of masked pixels.

With both protected/source masks nonempty and seven selections,
full-grid EDT calls decrease from `7×2=14` to
`2 baseline + 2 expanded guards = 4`. Boolean AND/OR,
source-inventory screening, tile selection, same-coadd holdout,
fit-quality gates and original sky values are still recomputed
for each fixed trial.

**Memory tradeoff:** the cache holds two Boolean source-mask
snapshots and two Boolean threshold-result rasters. At
2048×2048 this is **approximately 16 MiB** of additional
byte-packed-in-NumPy Boolean `itemsize=1` array storage,
before other temporary arrays. Unlike caching float64 EDT
distance fields, this memory use remains bounded and the
cache is discarded when the Stage03 call completes. It
is not a measured reduction of process peak RSS.

**Tests:** compare every cached train/holdout tile mask and the
entire selection QA dictionary against a no-cache run for
all six variants plus baseline; instrument EDT calls to
require exactly four. Deliberately mutate a source-exclusion
mask and a target-protection mask after cache creation;
the changed regions **must not** reuse stale guards. This
is geometry preservation, not evidence of truly empty sky.

### P40 — repeated serialization of an unchanged native WCS header per radial aperture — FIXED

**Source:** `measure_native_aperture` previously computed
`grid_fingerprint(wcs, shape)` for every radius, and
`_verify_frozen_growth_curves` recomputed it again.
This calls `WCS.to_header(relax=True)`, sorts header cards,
constructs a JSON payload and hashes it. The frame already
owns its native celestial WCS and shape independently of
external callers.

**Correction:** `NativeBandFrame.native_grid_fingerprint`
calculates the **same exact** `lsbpipe.scene.grid_fingerprint`
once and memoizes the result in that frame. Stage04
aperture receipts and Stage04.3 frozen-curve verification use
the cached value. The public `frame.wcs` property **still
returns a deep copy**; no caller obtains a mutable pointer
to the private source WCS through the cache. Changing
a separate WCS object does not change the frame fingerprint.

The expensive boundary-based WCS **affine projection quality
check is still performed for every ellipse** and must not be
confused with the cached geometry *identity*. Tests check
exact SHA/string equality against an independently
serialized WCS, one fingerprint serialization across a
ten-radius full growth curve, per-frame (not global)
cache boundaries, and isolation from external WCS edits.
Existing receipt provenance and scientific status remain unchanged.

### P41 — repeated catalogue-search string assembly and intermediate Qt paints — FIXED

**Source:** `gui/scene_editor.SceneEditor._render_source_registry`
joined and case-folded ten fields of Tractor QA text for
**every** record, even with an empty search query. It
then cleared and rebuilt a `QListWidget` item by item,
allowing unnecessary viewport updates for inventories up
to the existing safety cap.

**Correction:** build text only when a real nonempty
search term exists; use Qt uniform item sizes for homogeneous
one-line rows, suppress viewport updates during the batch
rebuild and restore them before status/selection handling.
All filters still select the same source entries, search
remains synchronous and case-insensitive, and no source
evidence or QA state is lost. Existing offline Qt regression
covers all filter/search modes; the test now also verifies
uniform-size mode and repaint re-enablement.

This is a modest **main-thread efficiency** improvement,
not a measured large-catalogue responsiveness guarantee.
For very large source lists, a proper virtualized
`QAbstractListModel`, search-token index and controlled
debouncing may be warranted, but require identity/selection
and keyboard-accessibility tests.

### P42 — full 2D coordinate-stack allocations in Stage04 sector QA — FIXED

**Source:** `measure_native_aperture` constructs **eight
descriptive pixel-centre sectors** using dense 2D native x/y
coordinate grids plus a stacked `2 × n_pixels` matrix to
multiply by the local inverse east/north tangent Jacobian.
That was avoidable temporary allocation in every radial
aperture, especially for larger ellipses.

**Correction:** a private
`_sector_ids_from_broadcast_native_centres` computes 1-D
native x/y offsets, broadcasts their two affine combinations
directly into east and north tangent coordinates, and uses
the **same** `atan2(east, north)`, modulo 2π and floor
into 45-degree labels. It preserves WCS handedness,
the same local tangent basis, the same eight sector
area accounting, and original pixel-centre semantics.
No modification is made to **exact fractional-pixel
flux integration** or to the fact that sector boundaries
remain descriptive rather than true exact angular sectors.

A deterministic seeded regression compares the new
native sector labels **elementwise** against the previous
dense matrix expression on multiple rectangular grids
and shear/rotation combinations. The enlarged benchmark
also checks this identity before reporting elapsed times.

### New performance evidence and verification boundaries

The existing offline benchmark remains the documented
command, now with **seven** reference-vs-optimized
computational kernel comparisons:

```powershell
python scripts/benchmark_desidera_performance.py --size 512 --repeats 3
```

Its new Stage03 variant recomputes baseline and fixed tile
controls with and without guard caching and verifies
deterministic selected-tile counts, guard radii and quadrant
coverage. Its new WCS comparison checks the exact original
header identity. Its native-sector comparison checks identical
eight-sector labels.

A median kernel-time ratio is **not** an end-to-end Leo P
speed measurement or a measurement of Windows GUI latency.
Neither total elapsed time nor resident memory has been
measured for the actual DR11 quartet on the researcher's
machine by this code audit. The benchmark comparisons
are deliberately non-gating to avoid flaky cross-platform
CI timing conditions.

### Remaining prioritized bottlenecks and risks

1. **Original-FITS cold-start I/O:** FITS decompression,
   all four source hashes, redundant opening paths,
   sharing of verified immutable pixel buffers and
   GUI raster construction remain plausible major
   contributors; no shortcut may replace full content
   lineage with `mtime` as if it were immutable identity.
2. **Native RAM and multi-band concurrency:** three
   `float64` full-band layers, a uint32 MASKBITS raster,
   validity/reason masks, sky sensitivity variants,
   ePSF kernel/stars, Qt RGB/pixmap buffers and parallel
   background workers can overlap in memory. Record
   **real peak RSS**, not just summed array sizes.
3. **Stage04 scientific WCS projections:** caching the
   **fingerprint** does not remove the repeated 16-sky-
   boundary-point accuracy tests. Potential reuse of
   radial WCS basis needs per-radius affine-error checks
   and exact saved recipe equivalence before adoption.
4. **Stage03 still-repeated bright-Gaia rasterization
   and Scene/QA object construction:** no source-
   inventory or star-wing safety checks were skipped;
   any future reuse must include verified identities
   and protect the source-screening invariants.
5. **Stage01 empirical PSF and extended diagnostics:**
   profile ePSF iterations, star stamps/holdouts and
   optional wing probes on real original coadds,
   keeping QA completeness and scientific HOLD explicit.
6. **Qt large-catalogue/scene redraw:** evaluate a
   virtualized model, bounded preview-image cache and
   WCS-consistent overlay redraw only with real
   responsiveness and memory metrics.
7. **True performance benchmark:** a pinned real-data
   Leo P / dense-Tractor / large-extended-galaxy
   corpus with identical science receipt identities,
   CPU wall time, peak RSS, I/O byte counts, worker
   concurrency, GUI first-image time, and p50/p95
   event latency across at least cold/warm runs.

**v0.2.1.dev2 is a Draft-PR development revision,
not a tagged scientific release.** The original
NativeBandFrame arrays, FITS data, `main` branch
and SCIENCE HOLD are unchanged.
