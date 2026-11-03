# P10c — local original-FITS thumbnails in Galaxy Projects

**Status:** cache/display feature under CI. No astrophysical photometry, sky fitting, image resampling for science, PSF modeling, or source adoption.

## Product contract

Galaxy Projects now uses a **lazy selected-card preview**: if the selected
object project has a recorded DR9/DR10/DR11 dataset, a local Legacy
`image-<band>.fits.fz` file, and no compatible cached preview, its
thumbnail is generated in a single background worker. This reads a single
original IMAGE FITS for display only. Other projects in the library are
never bulk-decoded, and project library startup does not contact a remote
service.

- Preferred visual band: **r** if available, otherwise the selected band
  and then a bounded list of available g/i/z originals in the *same*
  release/region/brick. No cross-release substitution.
- Centre: the project's **accepted** target ICRS, not a SIMBAD or NED
  coordinate. TAN celestial WCS is used to locate the target.
- Field: the project's saved search square, with a 1800-native-pixel
  cap to bound decompression/crop memory. A single-brick edge may trim
  the display image; the cached receipt labels any partial/limited display.
  No missing pixels are zero-padded or silently invented.
- Rendering: existing finite-pixel area-mean display reducer, existing
  asinh display stretch, grayscale 8-bit PNG with FITS y-axis inverted for
  display. It is **not** flux-preserving calibrated photometry, a PSF-matched
  mosaic, or a colour composite.
- Cache: `<project>/.desidera/card-preview.png` and
  `card-preview.json` with project UUID, selected release/region/brick,
  image band, source path, source SHA-256 at generation, size/mtime, target
  field, rendering version, PNG SHA-256, display trimming and explicit
  `DISPLAY_ONLY_NOT_SCIENTIFICALLY_ADOPTED` status.
- Reopening Galaxy Projects reads only small local cache metadata and PNG.
  A source size/mtime or selected brick/field change invalidates the cache;
  original FITS **is not rehashed during startup**. PNG corruption is
  detected by its SHA-256. Source hash recorded at generation does not
  prove NERSC authenticity or complete FITS QA.
- A manual **Regenerate local FITS preview** button is available to
  explicitly replace a disposable cached preview. One background task
  at a time. Failure appears in the right-hand card and does not prevent
  opening the galaxy project.
- Legacy standalone `card-preview.png` without a receipt remains
  displayable, marked as having **unverified provenance** for compatibility.

## Safety and regression gates

- Never change or move an original FITS, object manifest, source catalogue,
  analysed data, mask, or QA report; only cache files under `.desidera`
  are written.
- No GUI blocking, auto-download, login, network request, or background
  science run. On failed FITS reads leave existing cache provenance
  unpromoted; require a manual retry after failure during the same session.
- Synthetic compressed-FITS tests prove WCS-target centring, output size,
  PNG readability, automatic cache reuse, source staleness, partial-brick
  trimming, and unchanged original inputs.
- Qt test proves actual Project Hub lazily creates its first thumbnail
  and does not decode FITS on subsequent application launch.
- This change does **not** promote dataset integrity status.
  The existing `FITS linked · integrity not checked` label is deliberately
  retained until a separate review of saved Dataset Manager QA receipts.

## Local Windows acceptance

After CI/merge, in the **existing** repository/venv:

```powershell
Set-Location "C:\Users\frigi\Documents\PhD\Desidera"
git status --short
git pull --ff-only origin main
$python = ".\.venv\Scripts\python.exe"
& $python -m pytest -q tests/test_card_previews.py tests/gui/test_card_preview_autogeneration_gui.py
& $python -m lsbpipe.gui
```

Select **Leo P** in Galaxy Projects. It should show a r-band grayscale
preview generated from its already downloaded DR11 original FITS. The
project folder contains a small `.desidera/card-preview.*` pair. Close
and restart the app: the cached preview appears without a FITS reload.
Confirm that the accepted RA/Dec, SIMBAD and NED receipts, original FITS,
source file size/mtime, scientific QA and configuration snapshots remain
unchanged. The manual button can refresh when needed.

**Remaining next work:** unified provenance-aware QA summaries in the
passport, plus a separately labelled RGB colour-rendering mode. Do not
implement background estimation or photometric extraction here.
