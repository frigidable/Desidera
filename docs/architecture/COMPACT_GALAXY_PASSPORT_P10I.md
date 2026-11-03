# P10i — Compact scientific passport and source/QA disclosure

This increment is a **presentation-only refinement** of the existing
P10g/P10h read models. It does not calculate photometry, modify FITS or
analysis configurations, change user-confirmed SIMBAD/NED identities, or
reinterpret saved dataset-QA evidence.

## Default Galaxy Projects view

The current image-led layout is retained: square Legacy RGB thumbnail
on the left, galaxy name/ICRS/confirmed SIMBAD and NED hyperlinks on
the right. Underneath are short sections:

- **Published values · NED:** compact, source-identified object type and
  heliocentric redshift, if present. Both are explicitly labelled
  published catalogue values, not adopted DESIDERA results. The reference
  code is still accessible in the expanded original confirmation evidence.
  Redshift is not presented as a distance.
- **Observation data:** selected dataset and four-original-file
  presence; existence is not silently promoted to FITS/WCS integrity.
- **FITS/WCS check:** a concise label such as *Prior FITS/WCS diagnostic ·
  metadata unchanged*, *Saved audit needs rerun*, or
  *No saved FITS/WCS audit*. Each has an explanatory tooltip.

Long sky/PSF/correlated-coadd-noise limitations, exact hash caveats,
SIMBAD–NED positional differences and literature refcodes live in the
existing disclosure, now called **Sources & QA details**. The read-only
reference labels are moved, not removed. If catalogue evidence is
corrupted, the existing fail-closed warning expands the details
automatically.

**Open project** is always the primary CTA for accessible projects.
A brief suggested action (inspect FITS, review saved QA, rerun QA, obtain
source files) is exposed through its tooltip. The action continues to
open the existing workspace with unchanged project ID and no automatic
science adoption. Missing project folders still use *Relink project…*.

The two native Qt popup-menu buttons are now simply *Catalogues* and
*Project tools*. Their native menu indicators remain; redundant manually
typed glyphs were removed. All actions remain available.

## Acceptance contracts

1. Leo P shows two confirmed hyperlinks, a compact source-labelled type
   and heliocentric z, DR11 selected dataset, present FITS count and a
   saved QA status. Its literature refcode, 6″ SIMBAD/NED separation,
   original-file caveat and saved QA hash limitations are visible after
   **Sources & QA details** is expanded, not in the default view.
2. dw 1343+58 has a normal *No image yet / no dataset selected* state,
   a brief reason to connect catalogues, and *Open project*, with the
   recommended *find data* step in the button's tooltip. It must not
   show a technical exception inside the image.
3. A saved diagnostic with unchanged file size/mtime remains **historical**;
   it does not establish unchanged bytes or scientific acceptance.
   Changing input file metadata still prompts a new dataset QA run.
4. No network activity is added on card rendering or selection beyond the
   existing explicit P10e RGB policy. No source file, receipt, accepted
   ICRS coordinate or scientific result is rewritten.

Update and test in the existing Windows installation and .venv:
```powershell
Set-Location "C:\Users\frigi\Documents\PhD\Desidera"
if (git status --porcelain) { throw "Local edits present" }
git pull --ff-only origin main
$py = ".\.venv\Scripts\python.exe"
& $py -m pytest -q tests/gui/test_scientific_galaxy_passport_gui.py tests/gui/test_saved_dataset_audit_passport_gui.py tests/gui/test_project_hub_galaxy_passport_gui.py
& $py -m lsbpipe.gui
```

The green/black block at the upper-right of the user screenshot has
not been attributed to this widget. It requires a separate Qt/GPU
rendering reproduction; this change does not guess at its cause.
