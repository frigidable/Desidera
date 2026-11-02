# LSBpipe P0.5 — Dataset Manager architecture and scientific contract

Status: **foundation; not yet a scientifically valid mosaic, joint coadd, full photometric workflow or GUI feature**.

## 1. Core invariants

- One galaxy has one immutable project UUID, any number of observational data sources and any number of analysis configurations.
- A source is a release/reduction/camera/provenance grouping. An asset is a specific original file. A collection groups assets for a purpose. Results belong to immutable analysis runs, never directly to the galaxy or a mutable GUI selection.
- Presence, content SHA-256, FITS/WCS validation, science suitability, observation independence and adoption are independent evidence states. Missing evidence remains UNKNOWN.
- Original FITS are never modified or silently moved, redownloaded or combined.
- A newer release is not necessarily deeper, independent, or more suitable for diffuse structures.
- Overlapping bricks of the same release are not independent measurements. Different releases may re-use exposures. No default cross-release or cross-region inverse-variance sum.
- A reproducible result must identify exact inputs, revisions, masks, PSF, sky method, calibration, software and scientifically relevant settings.

## 2. Data hierarchy and persistence

GalaxyProject -> DataSource -> DataAsset -> DatasetCollection -> AnalysisConfiguration -> AnalysisRun -> ResultArtifact.

Existing lsbpipe-object.json and lsbpipe-project.json remain untouched. Existing selected_dataset is a last-view hint, NOT the primary science selection. Introduce optional object-bound sidecar datasets/index.json, schema lsbpipe.dataset-index/1, with:

- sources: stable source_id; provider, release, region, import root, calibration status and exposure lineage;
- assets: stable asset_id, source_id, original path, brick, band, product type, size/mtime, optional SHA-256, FITS validity, observation lineage;
- collections: future brick mosaic groups, filter sets and observing-epoch groups;
- analysis_configurations: future immutable analysis snapshots;
- runs: future run receipts and data products.

Proposed directories (no forced relocation):

    UGC_2162/
      lsbpipe-object.json              # existing identity
      lsbpipe-project.json             # optional historical metadata
      dr9/ or data/dr9/                # ORIGINAL files remain in place
      datasets/index.json             # NEW optional registry
      dataset-audits/<run>/audit.json  # existing immutable QA receipts
      processed/{mosaics,psf,masks}/<key>/
      analyses/<config-id>/<run-id>/
      results/<run-id>/

An absent sidecar on a legacy project means not yet indexed, not zero observations. External roots remain linked, not copied. A later relink process must preserve asset IDs and hashes without trusting filesystem paths alone.

## 3. Scientific dataset composition

### Same-release bricks

Coverage joining is a geometry problem first. WCS, flux units and pixel areas must be checked. In overlap, do not combine nominal inverse variances as if the repeated source observations were independent. For a prototype use a deterministic single-contributor footprint selection, retain per-pixel source provenance and test seams, PSF, sky residuals, noise correlation and masked-region recovery. Science mosaic status remains blocked until independent validation.

### Multiple Legacy releases

Primary: select a candidate for a named analysis goal and save that choice.

Compare: analyze DR9/DR10/DR11 in separate configurations; compare profiles, sky levels, calibration and morphology, keeping their shared-exposure limitations explicit.

Joint: fail closed until actual exposure/CCD IDs, processing lineage, calibration, PSF, sky, filter system and covariance have been demonstrated compatible. Same camera or same pixel WCS is insufficient evidence.

DR11 southern data include new DECam imaging and potentially i band. DR11 northern data consist of reprocessed existing DR9 contributions; the presence of i and actual available bricks must be confirmed per location, release and region. UGC 2162 must not be assumed to have DR10/DR11 coverage merely from release existence.

### Cross-survey data

HSC, CFHT, SDSS and DECam are different data sources. Photometric flux units, bandpasses, zeropoints, PSFs and background treatments must be explicit. Color gradients require PSF matching or forward convolution; morphology fits may retain native PSFs.

### SBF

Legacy coadds use Lanczos-3 and have correlated pixel noise. SBF eligibility requires its own power-spectrum noise model, PSF power spectrum, masking and calibrator validity. Do not infer suitability from existence of IMAGE, INVVAR, NEXP and MASKBITS alone. Single exposures or specially validated homogeneous regions may be preferable.

## 4. Dataset quality evaluation

Quality must be goal-specific (outer profile, Sérsic fit, color gradient, SBF). Audit coverage of galaxy center and broad clean sky annulus, background bias and large-scale residuals, saturated/bright-star masks, extended-source noise, PSF spatial variation, calibration and exposure overlap. Point-source limiting magnitude and INVVAR alone are insufficient to judge LSB performance. Missing diagnostic data -> INSUFFICIENT_EVIDENCE, not a favorable invented score.

Any recommended primary source requires a traceable metric breakdown, uncertainty, vetoes, quality-check version, and manual override recorded in the analysis configuration. An injection-recovery benchmark of synthetic diffuse objects against varying sky errors should precede calibrated automatic ranking.

## 5. Computation DAG and reproducibility

Data acquisition -> inventory -> SHA / FITS / WCS / coverage audit -> optional PSF/mosaic -> masks / background model -> radial profile and fits -> color analysis -> optional SBF -> result comparisons.

An immutable run should bind project ID, configuration ID/revision, exact source asset IDs + content SHA, processing version, relevant parameter dict, dependencies, method versions, random seeds (where applicable), code revision and output checksums. Do not reuse cache when an input has unknown or changed SHA.

Changing only a sky model should not redownload original data or rebuild unrelated products. On source modification, invalidate only downstream artifacts while preserving historical run records. A prospective content-addressed key helper is included in the P0 foundation; a full task scheduler is not.

## 6. GUI / UX plan

Project Hub -> Data Manager -> Data Inspection -> Analysis Configurations -> Results.

Source cards grouped by survey, release and processing region display actual local availability, status and verification receipts. Coverage Map overlays exact FITS WCS footprints, object/sky annuli, masks, exposure coverage, weak-weight regions and boundaries; the map cannot invent footprint information from a brick filename.

Separate actions: import/relink; verify checksum; validate FITS/WCS; inspect; configure primary/compare; propose quality-based source; optionally build a vetted mosaic. Joint remains disabled until the required physical independence checks pass.

## 7. Roadmap

P0.5a: new object-bound multi-source registry, explicit checksums, import existing Legacy/external FITS without movement, unlimited metadata discovery, conservative composition policy, cache-key primitive. No photometry or mosaic.

P0.5b: Dataset Manager GUI, multi-brick WCS coverage map, detailed QA statuses, arbitrary source selection.

P0.5c: immutable configuration and run records, dependency DAG and selective cache invalidation.

P0.5d: validated Mosaic Builder, true PSF/background/seam QA, injection-recovery research validation and independent release comparison.

P0.5e: external provider adapters and SBF covariance-aware data eligibility.

Acceptance gates: preserve UGC 2162 original FITS and project UUID; import existing directory without redownload; test >64 source sets, masked shared products, source changes, cache refusal, no automatic joint; test actual DR9 North/South data with a user-reviewed Windows GUI. Synthetic tests do not certify scientific suitability.

## Official technical references

- https://www.legacysurvey.org/dr11/description/
- https://www.legacysurvey.org/dr11/files/
- https://www.legacysurvey.org/dr10/description/
- https://www.legacysurvey.org/dr10/files/

These describe 3600x3600 coadds at 0.262 arcsec/pixel, ~130-pixel overlap, Lanczos-3, sky subtraction, band-dependent calibration and exposure provenance (survey-ccds). They explicitly caution against treating standard coadds as precision-ready measurements.
