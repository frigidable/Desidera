# P10e — automatic Legacy colour previews and uncluttered galaxy passport

## User-visible behavior

The previously working remote RGB cache remains unchanged. This step makes
it automatic **only for the currently selected galaxy**. After 650 ms of
stable selection, the Project Hub requests a colour image using the same
official Legacy Surveys cutout endpoint; rapidly browsing the list cancels
the pending request. One background fetch is allowed at a time. There is
no bulk prefetch for the galaxy library.

A valid RGB cache is used immediately without network access, otherwise
the existing cached local FITS monochrome image is displayed until colour
arrives. On failure, the local preview is retained (no science data changes);
the failed project/request key is not automatically retried again during
the same application session. **Project tools → Refresh colour preview**
allows explicit retry. An optional user-level setting
\`DESIDERA_DISABLE_AUTOMATIC_RGB=1\` suppresses automatic network activity
while leaving explicit refresh in the menu.

The cache is deliberately tied to the project ID + exact URL, including
dataset/release/layer, search field and project ICRS; selection switches
cannot attach the wrong downloaded preview. It remains display-only, and
does not modify original FITS, adopted ICRS, source measurements or
SIMBAD/NED identity receipts. It is not a substitute for science QA.

## Compact card design

The card now has two stages of text flow. A square image at left shares
its *hero* with the object's name, coordinates and confirmed primary
catalogue names at right. Further fields span the **full inspector width
under the square image** rather than continuing indefinitely in a narrow
column. This is intentional responsive wrapping; it is not a fake
HTML floating image nor a second static text pane. The image has no
caption underneath. Provenance is available through a tooltip and
existing cache receipts.

The compact overview shows dataset and local FITS status. Detailed
catalogue coordinates, SIMBAD/NED offsets and literature redshift
with refcodes, saved project path and last-opened timestamp appear only in
**Catalogue & project details**, a user-expandable section. Invalid
evidence remains explicit and automatically expands the diagnostic view
instead of hiding a provenance problem.

The **Open project** button remains the primary action. A two-button row
provides **Catalogues** (SIMBAD/NED) and **Project tools**
(local original-FITS preview refresh, colour refresh, immutable analysis
configurations and removal from the *recent* list). No feature is deleted:
the old full-width controls are retained for compatibility but no longer
occupy the right panel.

## Scientific and network boundaries

The colour image is visual only, never background-subtracted/scientifically
adopted, and is never used as a photometric input. Any ambiguity between
SIMBAD and NED remains visibly sourced and unadopted. The remote
downloader's HTTPS host allowlist, JPEG MIME/size checks, hash/receipt
validation and finite 20-second timeout remain unchanged.

## Windows acceptance — current repository and .venv only

After full CI and merge:

\`\`\`powershell
Set-Location "C:\Users\frigi\Documents\PhD\Desidera"
if (git status --porcelain) { throw "Local changes exist" }
git pull --ff-only origin main
$python = ".\.venv\Scripts\python.exe"
& $python -m pytest -q tests/gui/test_automatic_rgb_passport_gui.py tests/gui/test_legacy_rgb_card_gui.py tests/gui/test_galaxy_cards_gui.py
& $python -m lsbpipe.gui
\`\`\`

On **Leo P**, the colour image should load automatically after selection
if not cached. There is no Fetch button; the image has no caption.
SIMBAD/NED names are beside it, detailed measurements and evidence
are disclosed on request. Switching rapidly across several galaxies
should not fetch every candidate. Disconnect the network and reopen a
project with cached RGB: it must still display the existing colour image.
For a project without a cached RGB and with no network, cached local FITS
remains the fallback, and the Project Hub stays responsive.
