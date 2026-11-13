# Desidera performance audit — v0.2.1.dev0 (2026-11-13)

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
The three corrected bottlenecks were chosen because the scientific
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
three separately timed reference-vs-optimized kernel comparisons and
prints a JSON report:
1. two immutable native-frame materializations vs one;
2. ten full-mask snapshots vs byte-backed reusable snapshots;
3. each blank-ring exact Photutils rasterization vs a translated
   single exact ring template per radius.

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

`0.2.1.dev0` is an **installable development candidate** on Draft
PR #84, NOT a Git tag, merged main release, verified performance SLA
or validated astronomy science release. The same on-disk object-project
formats, legacy `lsbpipe` API, preserved original FITS and scientific
SCIENCE HOLD continue. The final change status is determined only by
exact HEAD GitHub Actions and dedicated observational benchmarks.
