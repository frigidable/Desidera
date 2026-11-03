# P10b — NED independent Galaxy Knowledge provider

**Scope:** manually initiated NED ObjectLookup identity discovery, immutable
NED evidence, explicit confirmation and cross-provider catalogue passport.
This does **not** implement scientific photometry, sky estimation, SBF,
literature aggregation or automatic distance adoption.

## Verified external contract

The official NED ObjectLookup service is documented at:
https://ned.ipac.caltech.edu/Documents/Guides/Interface/ObjectLookup

- Endpoint: `https://ned.ipac.caltech.edu/srs/ObjectLookup` (HTTPS POST).
- Body: form-encoded `json={"name":{"v":"Leo P"},"aliases":{"v":true}}`.
- Response: JSON; `StatusCode=100`, major `Version=2`, and
  `ResultCode=0,1,2,3` for non-interpreted/ambiguous/unknown/known names.
- `Preferred.Name`, `Preferred.Position.RA/Dec`, `Preferred.ObjType`,
  `Preferred.Redshift`, refcodes, and `Aliases` are separate NED statements.
- Only `ResultCode=3` with finite valid positions can be confirmed; even
  then a user must make the association explicitly.
- Unknown API major versions, unknown statuses, invalid positions,
  malformed schema or excessive response sizes fail closed. No automatic
  fallback to the deprecated `/cgi-bin/objsearch` services.

NED warns that its **legacy APIs no longer include data ingested after
January 2026**. This increment uses the documented JSON service, not the
legacy HTML/CGI API. Later position searches may use the published
`/tap/sync` ADQL endpoint, after independently validating its table schema:
https://ned.ipac.caltech.edu/Documents/Guides/Interface/TAP

## Local architecture

- `src/lsbpipe/surveys/ned_knowledge.py`: opt-in bounded ObjectLookup
  adapter, normalization, source-response preservation, independent immutable
  SHA-256 NED snapshots and explicit identity confirmations.
- `src/lsbpipe/gui/ned_knowledge_dialog.py`: offline-first review,
  manual network query, ambiguity display, confirmation guard.
- `src/lsbpipe/surveys/knowledge_summary.py`: independently read NED and
  SIMBAD receipts, compute positions relative to the **project** ICRS and
  calculate descriptive SIMBAD–NED great-circle offsets.
- `src/lsbpipe/gui/project_hub.py`: `Galaxy Knowledge · NED…` button
  beside SIMBAD. Confirmed NED name joins local project search keys.
- `src/lsbpipe/gui/galaxy_cards.py`: labelled catalogue identity and object
  type, catalogue coordinates, redshift tagged as **unadopted**, separate
  source evidence hashes and positional offset.

NED and SIMBAD are **different catalogues**, not interchangeable mirrors.
They can disagree in preferred names, positions, object types and redshifts.
We display both; no automatic reconciliation, redshift averaging, identity
auto-matching, or changes to user target coordinates or source FITS.

NED `Preferred.Redshift` comes from selected NED literature and is
displayed, if available, with its provenance reference code and explicitly
unadopted status. It is not a distance, not a measured velocity in every
frame, and never a substitute for TRGB/SBF measurements.

### Windows test (one established repository and one virtual environment)

After CI and merge:

```powershell
Set-Location "C:\Users\frigi\Documents\PhD\Desidera"
git status --short
git pull --ff-only origin main
$python = ".\.venv\Scripts\python.exe"
& $python -m pytest -q tests/test_ned_knowledge.py tests/gui/test_ned_knowledge_gui.py
& $python -m lsbpipe.gui
```

No new folders, worktrees or Python environments. Stop before `git pull`
if there are uncommitted changes.

Open **Galaxy Projects → Leo P → Galaxy Knowledge · NED…**. The dialog
must open offline without querying NED. Click **Search NED ObjectLookup…**.
Review the returned name, type, position and bibliography. If the result is
known and unambiguous, explicitly confirm. Reopen the program offline:
SIMBAD and NED must be shown separately in the card. Try searching by the
confirmed NED preferred name. Check that original FITS hashes, project
coordinates, project aliases and science reports are unchanged.

NED availability, up-to-date live API response compatibility and network
configuration on the Windows host are **not** established by offline CI.
If no result is returned, inspect the error text; do not assume NED and
SIMBAD use equivalent naming conventions or coverage.

## Future increments

- NED cone/TAP candidate discovery and explicit conflict resolution.
- Extended redshift/distance bibliographic measurements with method and
  uncertainty provenance (never auto-adopt or blend).
- Local FITS thumbnail generation and coherent QA cards.
- Photometric core remains separate and deferred.
