# Desidera P0.6c: original Legacy CCD catalogue acceptance

**Status: acceptance not yet certified until a real-archive workflow completes.**
This package adds evidence collection and strict schema testing, not a scientific result.

## Archive contract and selected samples

Official format documentation:
- [DR9](https://www.legacysurvey.org/dr9/files/)
- [DR11](https://www.legacysurvey.org/dr11/files/)

The test uses these **exact original files** from the official NERSC Legacy archive:

1. `dr9/south/coadd/155/1554p180/legacysurvey-1554p180-ccds.fits`
2. `dr11/south/coadd/155/1554p180/legacysurvey-1554p180-ccds.fits`
3. `dr9/north/coadd/210/2104p543/legacysurvey-2104p543-ccds.fits`

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

After obtaining successful real-data reports, review actual hashes and column formats, expand instrument coverage if northern `2104p543` does not include both `90prime` and `mosaic`, and separately validate project binding using original coadd IMAGE data.
