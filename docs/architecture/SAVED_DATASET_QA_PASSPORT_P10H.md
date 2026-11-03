# P10h — Saved Dataset QA Evidence in the Scientific Galaxy Passport

P10h integrates an existing, explicitly saved Dataset QA report into Galaxy Projects. It makes no network requests, changes no FITS, does no recalibration, sky modelling, photometry or science adoption, and uses the same application/virtual environment in the existing Desidera directory.

## Existing audit provenance

Survey Explorer's Dataset QA already produces an immutable JSON file at `dataset-audits/<run_id>/audit.json`. Its reference is saved in `ObjectProject.last_dataset_audit`. The report identifies project UUID, source data_root, target/field, creation time, a dataset list keyed by release/region/brick/band, FITS/WCS validation status, original-FITS SHA-256 values and file size/mtime. Its scientific status is `DIAGNOSTIC_NOT_ADOPTED`.

The P10h passive passport reads *only* a saved receipt explicitly referenced by the selected project. Its relative path must have the original three-component `dataset-audits/<run_id>/audit.json` layout. Symlinks, missing or oversized receipts, incompatible schemas, another project UUID, mismatched source root, target coordinates or survey selection fail closed. The selected band's record must be unique.

When the selected record says `FITS_WCS_VALIDATED_DIAGNOSTIC`, the passport also requires all four original products IMAGE/INVVAR/NEXP/MASKBITS with exactly the expected relative paths, recorded SHA-256 strings and unchanged current **size and mtime_ns**. This is a bounded metadata-level comparison, not live SHA-256 verification or pixel inspection. It can miss in-place changes where size and mtime are restored. Consequently the UI uses the deliberate wording:

> FITS/WCS checked in saved audit · metadata unchanged.

It **never** says that currently opened files have been independently rehashed or scientifically validated now. Previous FITS/WCS diagnostics do not certify sky/PSF modelling, correlated noise, photometric eligibility, multi-brick completeness or independence of exposures. Validity of the local JSON receipt itself is not cryptographically authenticated.

The compact library list shows `FITS 4/4 · prior WCS QA` only if that bounded check succeeds. Missing or untrustworthy evidence becomes `QA needs rerun`; a previous failed audit is labelled a **historical** issue, not a current verified failure. Without any saved audit, the card says `No saved FITS/WCS audit`.

## UI

The Observations and Saved Dataset QA sections provide independent claims: file presence vs historical diagnostic. The audit receipt's relative path and timestamp live inside the existing expanded details, never in an image thumbnail. The contextual primary action offers review or rerun, but still invokes the original safe **Open project** entrypoint. P10g's extra NEXT STEP line is hidden because it duplicated this button's wording. Catalogue identity, RGB preview and comparison semantics are unchanged.

## Acceptance and scientific limits

Unit tests cover: absence of a report, selection and project UUID mismatch, source-root mismatch, changed file stats, corrupt/unsafe receipt paths, symlinks, invalid report schema/status, a prior invalid-FITS report, and the explicit limitation that unchanged size/mtime do **not** establish unchanged bytes. Windows Qt tests exercise all three user states: pending, historical FITS/WCS validation and stale evidence, including the removal of duplicate NEXT STEP text.

There is still no independent current SHA-256 check at Galaxy Projects startup. For authoritative data revalidation, use the original Dataset QA workflow inside the galaxy workspace; then revisit the passport. Original FITS are read-only throughout.

Existing Windows environment:
`C:\Users\frigi\Documents\PhD\Desidera\.venv`.

Suggested focused tests:
```powershell
Set-Location "C:\Users\frigi\Documents\PhD\Desidera"
$py = ".\.venv\Scripts\python.exe"
& $py -m pytest -q tests/test_galaxy_passport_audit.py tests/gui/test_saved_dataset_audit_passport_gui.py
& $py -m lsbpipe.gui
```
