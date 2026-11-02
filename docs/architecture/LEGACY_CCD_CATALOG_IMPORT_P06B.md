# DESIDERA P0.6b+ — offline original Legacy Survey brick CCD importer

Status: Draft; scientifically **NOT CERTIFIED**. This work extends the
metadata-only Observation Registry; not an image-processing or photometric
algorithm.

## Authoritative file contracts

- DR9: https://www.legacysurvey.org/dr9/files/ — Image Stacks:
  `<region>/coadd/<AAA>/<brick>/legacysurvey-<brick>-ccds.fits`
- DR11: https://www.legacysurvey.org/dr11/files/ — same per-brick CCD table,
  and potential `i` band in the southern survey data.
- DR10: https://www.legacysurvey.org/dr10/files/

These per-brick files list the CCDs associated with building a brick. Their
catalogue layout and string widths can change across cameras/releases;
this adapter addresses columns by *case-insensitive name*, not fixed offsets.
Documented minimum columns: CAMERA, EXPNUM, CCDNAME, FILTER, BRICK_X0/X1,
BRICK_Y0/Y1. The last four describe the CCD's overlap in brick coordinates.
Positive overlap is required. Additional columns remain untouched.

This API reads a local FITS BINTABLE only. It **does not** parse the global
`survey-ccds-<camera>-<release>` summaries, infer release coverage from a
filename, or fetch a file automatically.

## Conservative acceptance rules

1. User selects an existing Legacy **coadd IMAGE** asset whose byte SHA-256
   has been verified in Dataset Manager. Missing/changed coadd evidence is
   rejected. FITS/WCS validity is a separate status.
2. The original per-brick CCD filename must match the exact selected brick.
   Filename/selection alone **cannot prove** a table belongs to the selected
   release or region; this remains a `CATALOGUE_ASSERTION_NOT_VALIDATED`.
3. A single, nonempty FITS BINTABLE is required. Case-insensitive mandatory
   columns, valid integer EXPNUM, ASCII identifiers, known cameras
   (`decam`, `90prime`, `mosaic`), valid filters (`g,r,i,z`),
   nonzero footprint intersection, and no duplicate CCD rows.
4. Import selects CCDs only for the asset's band; the same per-brick original
   catalogue may serve separate g/r/i/z original coadd IMAGE assets.
5. SHA-256 of the original catalogue, its location, size, mtime, table row
   counts and selected band are recorded with each specific coadd asset.
   A checksum is **content integrity**, not calibration validity.
6. Processing is offline, fail-closed and atomic. Invalid tables change no
   existing lineage. Identical imports are idempotent. A different catalogue
   SHA for an already-linked asset is rejected pending explicit review.
7. The lineage report may flag repeated CCDs or shared exposure groups
   across DR9/DR10/DR11. It **never** certifies disjoint observations,
   independent noise, optical calibration, PSF/sky consistency, or
   suitability for a scientific coadd/SBF.
8. No original FITS is moved or edited; no network request or photometric
   work occurs. Existing project UUID, legacy manifests, QA artifacts, and
   source-asset identifiers remain unchanged.

## User workflow

Within an open galaxy project:

- Data Manager → Index local FITS folder; select an original Legacy
  `IMAGE` (not INVVAR/NEXP/MASKBITS).
- Verify selected file SHA-256 if it is still NOT_HASHED.
- Select the same original IMAGE, choose **Index brick CCDs…** and select
  its local `legacysurvey-<brick>-ccds.fits`.
- Review the CCD catalogue assertion count in the source detail panel.
  The status remains unverified until a separate catalogue association
  audit and physical/exposure completeness checks exist.

Direct Python API:

    from lsbpipe.surveys.legacy_ccd_catalog import (
        inspect_brick_ccd_catalog, import_brick_ccd_catalog)
    preview = inspect_brick_ccd_catalog(
        project_root, coadd_asset_id=asset_id, catalogue_file=local_file)
    receipt = import_brick_ccd_catalog(
        project_root, coadd_asset_id=asset_id, catalogue_file=local_file)

The preview is read-only and does not create a lineage sidecar.

## Known scientific limitations and future acceptance

- A complete physical-exposure independence audit requires verified
  release-specific CCD catalogue provenance, completeness against actual
  imaging contributions, and a canonical cross-release identifier map.
- The official integer EXPNUM does not contain leading zero formatting;
  it is converted to a decimal string without invented zeros.
- Rejection of multiple BINTABLE HDUs is deliberate: the program must not
  silently choose an arbitrary data table.
- A reported per-brick CCD relationship is not a per-pixel noise
  covariance map and should not determine mosaic weights.
- Independent scientific adoption still requires original coadd FITS/WCS
  QA, optical calibration, PSF, sky bias, correlated noise and injection
  recovery. For SBF, a separate noise power spectrum and valid distance
  calibration are mandatory.
- Follow-ups: file provenance/release path checks, diagnostic transparency,
  richer CCD metadata and reporting, conflict-resolution workflows,
  versioned catalogue supersession and validated cross-release comparison.

No official remote CCD files have yet been checked against this local parser
in CI. Synthetic real FITS BINTABLEs exercise its schema and error paths;
a release-specific real-data acceptance campaign is required before claiming
production scientific catalogue compatibility.
