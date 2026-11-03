# Desidera P0.6c: original Legacy CCD catalogue acceptance

**Status: LIMITED REAL-ARCHIVE PARSER ACCEPTANCE PASSED; full scientific lineage/photometry NOT CERTIFIED.**
This package adds evidence collection and strict schema testing, not a scientific result.

## Archive contract and selected samples

Official format documentation:
- [DR9](https://www.legacysurvey.org/dr9/files/)
- [DR11](https://www.legacysurvey.org/dr11/files/)

The test uses these **exact original files** from the official NERSC Legacy archive:

1. `dr9/south/coadd/155/1554p180/legacysurvey-1554p180-ccds.fits`
2. `dr11/south/coadd/155/1554p180/legacysurvey-1554p180-ccds.fits`
3. `dr9/north/coadd/210/2100p542/legacysurvey-2100p542-ccds.fits`

They are initial nominated samples, not yet authenticated to exist. A missing sample must fail the run. Replacement must be an explicit reviewed code change. The paths are source assertions; a matching path alone does not independently authenticate the release or region.

## Run

```bash
python -m pip install -e ".[astro,test]"
python -m pytest -q tests/test_real_legacy_ccd_acceptance.py
python scripts/accept_legacy_ccd_archives.py --report-dir real-ccd-evidence
```

`--skip-download` uses the original saved local FITS. The workflow `.github/workflows/legacy-ccd-real-archive.yml` runs the checks for the PR and uploads JSON/Markdown evidence to GitHub Actions even if a test fails.

## Requirements for a valid limited acceptance

- Every nominated file must be downloaded as original FITS bytes; no synthetic fallback.
- Each file must pass Desidera's strict BINTABLE parser and report content SHA-256, size, column formats, counts of CCDs and distinct exposures.
- An **independent Astropy-only reader** must reproduce exactly all row identities (camera, exposure, detector, filter, row location).
- The DR9 South/DR11 South shared brick comparison must report both exact shared CCD identities and shared exposures.
- No missing or unreadable sample can be treated as a pass. Results must carry `PASS_LIMITED_PARSER_ACCEPTANCE` or `INCOMPLETE_OR_FAILED`.

## Deliberately outside certification

No real coadd IMAGE files are fetched in this stage. This does **not** certify the real end-to-end Dataset Manager integration, pixel-level CCD completeness, release authenticity beyond the archive path, calibration, WCS, PSF, noise, independent measurements, photometry, coaddition, or SBF. Shared CCD detection is positive evidence of reuse; no shared entries **never** proves independence. All scientific adoption flags remain false.

After obtaining successful real-data reports, review actual hashes and column formats, expand instrument coverage if northern `2100p542` does not include both `90prime` and `mosaic`, and separately validate project binding using original coadd IMAGE data.

## Verified real-archive evidence — 2026-11-03 UTC

The [CI run 41546160591](https://github.com/frigidable/Desidera/actions/runs/41546160591) completed successfully. Its downloadable artifact contains the complete JSON and Markdown reports. All three original FITS tables were extracted by the Desidera importer and an independently implemented Astropy reader, with exact row-identity agreement.

| Archive source | Original SHA-256 | CCDs | Exposure groups | Cameras |
|---|---|---:|---:|---|
| DR9 South 1554p180 | `e0d1733171ab464c423ee6687a56a952ecb604f23709bf30995cdb55b6077f78` | 39 | 9 | decam: 39 |
| DR11 South 1554p180 | `8085050046a370eb461896db0b08905195e30cb3cbac93fa8953a1712c3ef7fc` | 39 | 9 | decam: 39 |
| DR9 North 2100p542 | `db3eeb740f6c463b37329ce1588ad44178989696090e386254c025f8390e2bed` | 33 | See JSON report | 90prime: 22, mosaic: 11 |

DR9 South and DR11 South share **39/39 CCDs from 9 shared exposure groups**, despite different file SHA-256 values. The importer therefore must compare physical CCD identity and must not treat different release files as independent exposures.

The original south tables contain `BRICK_X0=-10`, `BRICK_Y0=-8`, and high edges up to `BRICK_X1=3612` or `BRICK_Y1=3607`. P0.6c fixes the parser to accept these valid signed extents while rejecting zero intersection with the standard 3600×3600 brick area.

The three samples validate format compatibility for the observed camera types and these fields only. No inference of pixel-scale observational completeness, release authentication beyond archive URL, or SBF/photometry readiness follows.
