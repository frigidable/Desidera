# P0.5b Data Manager — offline GUI evidence and scientific scope

Status: Draft implementation based on P0.5a; **not** a photometry-quality or mosaic release.

## Entry point and persistence

Open an existing galaxy via Project Hub and use the dedicated **Data Manager · local sources** action on the left. The current Remote Sky and FITS Explorer tabs are preserved. The first visit to an object with a saved Legacy FITS root and no dataset index registers that local inventory **once** in a sidecar at `datasets/index.json`. Subsequent visits read the same record and never re-download or move original images.

The registration action is metadata-only. A user can select an additional pre-existing Legacy root. Each external source is a distinct source ID. Project identity, ICRS center, legacy last-selected quartet, original FITS, old QA histories and existing scientific reports remain intact.

## User interface

- Left: source groups by *actual* provider/release/processing region and underlying original assets by brick, band, product kind. Statuses distinguish file presence, SHA-256 and FITS verification. No fake DR10/DR11 for UGC 2162.
- Right: the object's sky field, target aperture and sky-ring geometry. When requested, this adds each available image asset's **native FITS WCS header** footprint, with east-left/north-up orientation and source selection emphasis. It never invents polygons from a brick filename.
- Header actions: **Index local FITS folder**, **Inspect WCS footprints**, **Verify selected file SHA-256**, **Back to FITS Explorer**.
- Footprint computation reads headers only; malformed/missing FITS are counted as rejections and unprocessed assets as explicit display-cap omissions. The diagnostic map does not claim to represent valid-pixel coverage, weight maps, deep noise, PSF, sky gradients or scientific mosaic coverage.

## Quality boundaries

The screen must not classify a dataset as scientifically suitable from WCS/FITS presence alone. Exposure independence and calibration suitability stay UNKNOWN/UNASSESSED. There is no automatic selection, coaddition, PSF matching, artificial survey availability or SBF certification in P0.5b.

Future phases: mask/invvar/nexp-weighted coverage, exact segmentation and sky-region checks; immutable analysis configurations and cross-release comparisons; then guarded Mosaic Builder with exposure-lineage/PSF/sky/covariance QA. The P0.5b GUI lays out where those metrics will appear, but does not manufacture them.

## Validation

Unit tests use artificial WCS headers to verify RA orientation, geometry, real native raster dimensions, explicit rejected inputs and preservation of FITS bytes. Qt tests cover navigation, source cards, file status and checksum action. Windows CI screenshots at 1366x768 and 1600x900 use synthetic external FITS only and must never be represented as real UGC 2162 data.

On-device acceptance: import existing `LSBpipeData/UGC_2162`, confirm stable UUID and no redownload; inspect genuine DR9 North/South WCS boundaries, verify single source SHA and restore FITS Explorer. Compare map footprint visually to real image corners; numerical photometric acceptance is a separate QA task.
