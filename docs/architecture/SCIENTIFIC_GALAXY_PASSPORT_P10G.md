# P10g — Scientific Galaxy Passport (pre-photometry)

This is a presentation and evidence-readiness upgrade. It does not introduce background estimation, photometric measurements, catalogue network queries or automatic literature adoption.

## Confirmed identity and literature

Only receipt-validated, explicitly confirmed SIMBAD and NED identities may be displayed as associated objects. The card explains why a researcher might link both catalogues: aliases, cross-identification, bibliography and future comparison against published values. SIMBAD and NED can disagree; a positional difference is descriptive and does not automatically invalidate either association.

The existing NED ObjectLookup read model provides an object-type code and optionally a heliocentric redshift with reference code. The card labels these as unadopted external catalogue statements, never as a verified Desidera measurement or distance. Missing values are simply omitted. Invalid evidence remains fail-closed with an explicit diagnostic. Confirmation never changes the immutable project identity or its accepted ICRS position.

## Read-only observation readiness

The new galaxy_passport_readiness module checks selected dataset and original IMAGE, INVVAR, NEXP and MASKBITS presence for the exact current release, region, brick and filter. It checks at most eight deterministic file paths in the local data root and optional data subfolder. No FITS data or headers are opened; no checksum is calculated, directory tree recursively scanned, or source file written.

States are TARGET_NOT_SET, DATASET_NOT_SELECTED, DATA_ROOT_NOT_LINKED, DATA_ROOT_MISSING, NO_MATCHING_ORIGINAL_FITS, INCOMPLETE_ORIGINAL_FITS, FOUR_ORIGINAL_FITS_PRESENT_UNVERIFIED and DATASET_INVALID. The contextual action continues to open the existing galaxy workspace; it never downloads data on the Projects screen.

Four present files mean only that four nonempty paths exist. They do not certify FITS/WCS integrity, checksum validity, PSF, sky, science eligibility, independence or coverage. Existing last_science and last_dataset_audit fields are not promoted because their current selection and freshness have not been independently validated here. A later stage will interpret verified dataset-audit receipts.

## UX details

The card uses the top square image/identity layout already present in P10e. When catalogues have not been linked, it explains their practical value and points to the existing Catalogues menu. After NED confirmation, source-labelled type/redshift and positional comparison are shown. Observational readiness and Next step form separate sections. Long provenance continues to live behind Catalogue & project details.

For projects without original IMAGE FITS or a dataset selection, the thumbnail shows the neutral text No image yet. Exceptions, filenames and stack traces never replace the image. Manual local-preview refresh is disabled if there is no matching local source, and cached remote RGB remains a decorative image unrelated to FITS readiness.

## Windows acceptance

Use the existing repository only: C:\Users\frigi\Documents\PhD\Desidera and its existing .venv.

Run pytest for tests/test_galaxy_passport_readiness.py, tests/gui/test_scientific_galaxy_passport_gui.py and tests/gui/test_project_hub_galaxy_passport_gui.py, then start the GUI.

Check Leo P: source-labelled NED values, confirmed SIMBAD/NED links, comparison offset and explicit FITS status. Check dw 1343+58 without a selected dataset: No image yet and an action to find observations, without a FileNotFoundError message. Original files, object manifests and all confirmation receipts must remain unchanged.

The next phase may add evidence-verified dataset audits and richer literature measurements, but not photometric sky modelling in this package.
