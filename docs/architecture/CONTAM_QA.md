# CONTAM-QA / COG acceptance (experimental, not adopted)

Read-only offline triage of an existing OUTER-PHOT g/r report:
```powershell
python scripts/contaminant_audit.py "C:\Users\frigi\DESIDERAProjects\NGC2146sat\analyses\project-diagnostics-qa.zip" --output "C:\Users\frigi\DESIDERAProjects\NGC2146sat\analyses\contam-qa.csv"
```
No rerun, cache change, network, source subtraction, or mask changes. The
CSV is a catalogue-centre inventory, not an assessment of all extended
source footprints, uncatalogued sources, or valid pixels. 120 arcsec is a
review radius, not a statement of measured coverage. Full source IDs are
merged across bands. The gnomonic ellipse uses the OUTER reference
geometry. Catalogue flux is not aperture contamination. Gaia assessment
is only a candidate flag. Negative catalogue flux is retained; undefined
colour remains blank. Physical membership and treatment remain unresolved.

QA(6) evidence: 94 entries in the union of the a<=120 arcsec centre
selection and protected-overlap entries; 12 centres inside 40 arcsec,
12 at 40–65 arcsec, 70 at 65–120 arcsec. No automatic foreground
classification within 65 arcsec. Review OBJIDs 394, 424, 365, 388,
322 jointly with the diffuse component; do not independently subtract
their Tractor models. A Gaia association alone is not membership evidence.

Next numerical stage: local PSF and joint diffuse + compact models on
original validated numerical pixels with hard-mask and WCS identity
checks; retain keep / model / mask hypotheses, verify injection recovery,
report signed residuals and flux changes on matched support, and rerun
geometry and sky only if the scene changed. No new treatment adopted
from this report.

COG acceptance: measure signed annular sums with common support,
sector coverage and covariance; use independent sky controls outside
the tested aperture, not a background tuned to flatten the COG.
Require calibrated multi-ring equivalence to a small residual growth
tolerance, not simply failure to reject zero slope. Test false plateaus
in injection experiments, retain nonmonotonic curves, and distinguish
observed aperture flux from extrapolated Sérsic totals. Extend only
where native coverage and sky control support allow. Existing OUTER
results do not establish a plateau or calibrated uncertainty.

## Numeric snapshot for the next review

Run only when a numerical source decomposition or sky/COG investigation
is needed. No diagnose-project run is required:

```powershell
python scripts/export_contam_inputs.py --project "C:\Users\frigi\DESIDERAProjects\NGC2146sat"
```

Produces analyses/contam-inputs.zip from existing verified SH1/SKY caches
and the same saved PSF0 lineage as the latest QA. Each array is stored as
a non-pickle NPY with a SHA256 in the diagnostic-cache manifest. The
bundle retains the expanded original numerical image, IVAR, separate
hard/source/protection/sky-veto masks, reference halo residual, sky
training support, exact expanded WCS and source positions. No replacement
of missing data or change to existing caches. It is a snapshot, not a
fresh verification of the source FITS. Missing or inconsistent caches,
PSF or WCS block export rather than trigger a long rerun. Existing
export files are not overwritten.
