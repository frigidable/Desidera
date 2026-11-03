# P10d — square Galaxy Passport and official Legacy Surveys RGB thumbnails

**Scope:** presentation/UI, visual network cutouts and recoverable cache only.
No photometry, background determination, PSF inference, distance estimate,
automatic catalogue association, FITS replacement or science QA approval.

## Official service

The team-maintained Legacy Surveys JPEG cutout service is documented by
DR9/DR10/DR11 release pages:

- https://www.legacysurvey.org/dr9/description/
- https://www.legacysurvey.org/dr10/description/
- https://www.legacysurvey.org/dr11/description/

The colour image is obtained from
`https://www.legacysurvey.org/viewer/jpeg-cutout` with explicit
`ra`, `dec`, `layer`, `pixscale` and `size` parameters. It is not a
WebEngine screenshot or an original coadd FITS. Pixel dimensions are a
**square 384 × 384** (below the service's documented 512-pixel limit).
The angular field is the *project's* saved search-field side, with
`pixscale = width_arcmin * 60 / 384`. No match is made to NED or SIMBAD
coordinates unless the user explicitly changes the project centre.

Validated release/region mapping:

| Object project dataset | JPEG cutout layer |
|---|---|
| DR9 South | `ls-dr9-south` |
| DR9 North | `ls-dr9-north` |
| DR10 South | `ls-dr10-south` |
| DR11 South | `ls-dr11-south` |

DR10 and DR11 north use DR9 north imaging in the official viewer documentation.
P10d does **not** silently present it as a release-independent DR10/11
north-specific RGB source. For unmatched combinations, only the existing
local FITS thumbnail remains available.

## UI

The right-hand `GalaxyDetailCard` has a *single horizontal content region*:
a fixed-size **square visual panel at left** and a scrollable scientific
passport at right. Images use 184 × 184 logical Qt pixels, or 124 × 124
for narrow inspector widths. Catalogue identity, coordinates, dataset, FITS
status and project history stay separate from the image. The *Open project*
action remains visible at the bottom; no window redesign or second library.

The existing `Regenerate local FITS preview` button is retained.
A new `Fetch Legacy Surveys colour preview` button explicitly initiates
a network call for the **selected** galaxy only. The colour fetch is
bounded to one background worker and does not block GUI events. Opening,
selecting, filtering or re-opening Galaxy Projects **never requests RGB
from the network**.

The rendering priority is:
1. Valid cached remote RGB (labelled `Legacy Surveys RGB · display only`).
2. Valid local original-FITS grayscale cache.
3. No preview.

A missing service, HTTP error, invalid MIME, bad JPEG, unavailable Internet,
invalid region/release or stale receipt does **not** discard the local FITS
preview. Qt also attempts a local fallback if a JPEG passes header validation
but cannot be decoded for display. Remote cutouts do not change source
coordinates or original FITS.

## Separated provenance and cache

Original-FITS display remains stored as
`.desidera/card-preview.png` plus `.desidera/card-preview.json`.

Remote display is stored separately as
`.desidera/card-preview-remote-rgb.jpg` plus
`.desidera/card-preview-remote-rgb.json`.

The remote receipt records project UUID; saved target RA, Dec and
field width; selected data release, region and brick; exact allowed
viewer layer, JPEG size and pixel scale; the requested and final URL;
server MIME; acquisition time; SHA-256 of the saved JPEG; and
`DISPLAY_ONLY_REMOTE_RGB_NOT_SCIENTIFICALLY_ADOPTED`. Stale or corrupted
remote caches do not count as valid. The source hash establishes local
bytes at the time of receipt, **not independent archive integrity or
astrophysical fidelity**.

Requests are restricted to official HTTPS `www.legacysurvey.org`.
Unapproved redirects fail; the content type must be JPEG, body is capped
at 4 MiB, Content-Length mismatches fail, JPEG dimensions must be 384 ×
384. A cache refresh writes only disposable thumbnail/receipt files
atomically; existing galaxy manifest, science runs, catalogue evidence
and original FITS remain untouched.

## Windows acceptance

Use the existing repository and existing virtual environment only:

```powershell
Set-Location "C:\Users\frigi\Documents\PhD\Desidera"
if (git status --porcelain) { throw "Local edits present" }
git pull --ff-only origin main
$python = ".\.venv\Scripts\python.exe"
& $python -m pytest -q tests/test_remote_card_previews.py tests/gui/test_legacy_rgb_card_gui.py tests/gui/test_galaxy_cards_gui.py
& $python -m lsbpipe.gui
```

Open Leo P. Its currently cached grayscale thumbnail should now be a
square on the left; the SIMBAD/NED passport should appear to its right.
Click **Fetch Legacy Surveys colour preview**. Successful live
Internet access should give an RGB image from
`ls-dr11-south` in that same square. Restart offline: cached RGB
persists without a network request. Try deleting *only* the remote
cache in a disposable test project (not your scientific source files);
the grayscale preview should be shown again. Real Windows network
compatibility must be accepted separately from mocked CI.

## Deferred scope

Updating the observational QA status in Galaxy Projects is a separate
P10e change. A viewer JPEG is decorative and cannot be used as source
photometry, original archive evidence, an SBF image, or a replacement for
a photometrically calibrated coadd.
