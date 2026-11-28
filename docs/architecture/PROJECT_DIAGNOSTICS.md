# Project diagnostic runner (experimental)

Run from the existing checkout, with the same GUI/astro dependencies used by
Scene. No window or network access is needed. Only saved mask edits are used.

```powershell
python -m desidera diagnose-project --project "C:\Users\frigi\DESIDERAProjects\NGC2146sat" --bands g r --objid 360
```

The explicit saturated foreground candidate is checked in each native band.
The same worker functions as Scene and Sky are invoked. By default the
saved project field/target anchor and the default source-mask policy are
used, not unsaved UI choices or an SGA prior. Use --width, --anchor,
--gaia-inner-preview or --red-extended-preview to explicitly change them.

The existing project analyses folder holds one cache ZIP per band per
heavy stage (SH1 including GEO-1, and SKY-2), and project-diagnostics-qa.zip.
Each archive is atomically replaced, contains numeric NPY arrays plus JSON
(no pickle), and checks cache recipe, source hashes, grids, masks and code
fingerprints. Original FITS, project manifest and mask drafts are read only.
FITS and mask preparation still run on every invocation to verify identity.
Sky-only module edits invalidate SKY but not SH1; other edits conservatively
invalidate both. --force-stage sky or sh1 or all forces recomputation
of the named stage and downstream stages. Caches are not science adoption.
The report includes both bands, pair diagnostics and PNGs, no replacement FITS.
A failed run does not publish a new final report; check its timestamp.
No fresh PSF0, injection recovery, final photometry or GUI cache restoration
is implemented in this first runner.
