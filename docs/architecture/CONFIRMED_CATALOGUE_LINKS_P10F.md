# P10f — confirmed SIMBAD/NED object hyperlinks

Galaxy Projects and the SIMBAD/NED confirmation dialogs now show an explicit
clickable link beside each **user-confirmed** catalogue identity. Opening a link
uses the platform's default web browser; there is no request while a project
is opened, refreshed, or displayed.

Official HTTPS navigation targets:

- SIMBAD: \`/simbad/sim-id?Ident=<confirmed main_id>\` on
  \`simbad.cds.unistra.fr\`, or the official \`simbad.harvard.edu\`
  mirror when the source confirmation records that endpoint.
- NED: \`https://ned.ipac.caltech.edu/byname?objname=<confirmed main_id>\`.

Links are generated at display time from **receipt-validated object identities**
provided by \`load_galaxy_knowledge\` or directly from the validated
SIMBAD/NED confirmation receipts. A stored candidate or nearby source, a
free-text project alias, or a name-only query cannot create a trusted link.
Any unreadable/corrupt evidence causes the existing scientific evidence
warning to appear, with **no link**.

Only these allowlisted HTTPS hosts/paths can be constructed. Object names
are URL-encoded and HTML-escaped before being passed to Qt's rich-text
\`QLabel\`. The link is user-activated and never triggers background catalogue
lookups. The link destination is a *name-resolution viewer page*, not a
cryptographically persistent object record ID or a new identity assertion;
its displayed values may change as the database evolves.

No project UUID, accepted position, analysis configuration, FITS source,
derived photometry, external evidence receipt or science run is rewritten.

## Manual verification (Windows)

Use the same repository and existing virtual environment:

\`\`\`powershell
Set-Location "C:\Users\frigi\Documents\PhD\Desidera"
if (git status --porcelain) { throw "Local changes exist" }
git pull --ff-only origin main
$python = ".\.venv\Scripts\python.exe"
& $python -m pytest -q tests/test_catalogue_links.py tests/gui/test_confirmed_catalogue_links_gui.py tests/gui/test_knowledge_passport_gui.py tests/gui/test_ned_knowledge_gui.py tests/gui/test_galaxy_knowledge_gui.py
& $python -m lsbpipe.gui
\`\`\`

Choose Leo P. SIMBAD's confirmed **NAME Leo P** and NED's confirmed
**PGC1 5065058 NED001** should each be clickable in the primary passport.
Click each to open the corresponding named object page. Open each catalogue
dialog to check the same link under its saved confirmation. Disconnect the
network: both links should still render (navigation obviously needs access).
Unconfirmed candidates should not appear as hyperlinks.
