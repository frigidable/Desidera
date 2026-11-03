# P10a — source-neutral Galaxy Knowledge read model and SIMBAD passport

**Status: implementation under acceptance testing.** This is a UI/identity
improvement, **not** photometry, an approved scientific measurement, or NED
integration. It is based on the current merged `main`.

## Project boundaries

- The primary research entity is the immutable object-project UUID; project
  name/ICRS target coordinate/aliases remain under user control.
- A *catalogue identity* is a source-specific external assertion. It is stored
  separately from the project manifest in immutable, SHA-256 addressed
  `knowledge/simbad/snapshot-*.json` and
  `knowledge/simbad/confirmation-*.json` receipts.
- The only association displayed as **confirmed** is an explicit
  user-confirmed candidate, validated against the saved original query.
  Unconfirmed SIMBAD cone neighbours, even if closer to the target center,
  must never become galaxy aliases or searchable identities.
- P10a adds a source-neutral `CatalogueIdentity` /
  `GalaxyKnowledgeSummary` read model
  (`src/lsbpipe/surveys/knowledge_summary.py`). Its first concrete provider
  is SIMBAD. **NED is not implemented** and must use independent
  query/identity/evidence receipts when added.
- Opening or searching local projects does no remote calls, FITS reading,
  science run, image download, or project manifest modification. The card
  reads only lightweight saved JSON receipts.
- The catalogued ICRS coordinates and project target ICRS are displayed
  separately. Their great-circle separation is recomputed **relative to the
  currently saved target** (including RA wrap-around); no coordinates are
  automatically adopted or reconciled.
- The panel identifies the source, accepted name, separation, and a
  shortened confirmation receipt SHA-256. A corrupt/wrong-project receipt
  produces `EVIDENCE_ERROR` on the card, not a fabricated empty/confirmed
  state.
- The full original SHA-256 receipt remains available on disk; no file is
  renamed or rewritten.

## Workflow

1. Galaxy Projects → select a project → **Galaxy Knowledge · SIMBAD…**.
2. Manually search, review alternatives, and explicitly confirm the correct
   object. Querying and confirming preserves historical evidence.
3. Close the dialog. Project Hub rebuilds the read model and its search index.
   The right-hand card shows `SIMBAD: <main_id> · user-confirmed association`.
4. Search for the *confirmed* external identifier to find the project. Names
   from unconfirmed cone search results are **not** indexed.
5. Close/reopen DESIDERA: the user-confirmed identity is visible offline.
   No science result, project coordinates, original image or existing
   photometry is modified.

The source-neutral model keeps a separate provider/receipt identity and
coordinate-bibcode record so NED can later coexist without being
misrepresented as a SIMBAD mirror. A future read model must explicitly flag
cross-provider conflicts and source-specific measurements with units,
methods and bibliography. No automatic adoption or averaging is permitted.

## Windows local pilot (one existing repository and one existing .venv)

```powershell
Set-Location "C:\Users\frigi\Documents\PhD\Desidera"
git status --short
git pull --ff-only origin main
$python = ".\.venv\Scripts\python.exe"
& $python -m pytest -q tests/test_knowledge_summary.py tests/gui/test_knowledge_passport_gui.py
& $python -m lsbpipe.gui
```

Do not run `git pull` with uncommitted changes; do not create a second
worktree or Python environment. Run the commands after PR merge.

**Leo P acceptance:** the already confirmed `NAME Leo P` must appear in
its Galaxy Projects card, while its initial target ICRS stays at the saved
RA/Dec; SIMBAD star candidates remain invisible to the project's confirmed
alias search. The source remains readable offline. If a receipt is corrupt,
the card reports an evidence error.

## Next stages and limitations

- P10b: NED identity adapter through the supported new API, with
  independent evidence, literature-source metadata and conflict-aware joins.
- P10c: local FITS thumbnail generation and extensible project cards.
- P10d: unified audit summaries across Dataset Manager, QA and Galaxy
  Projects with no false science certification.
- P10a does not yet import SIMBAD alternative identifiers or object type,
  nor NED redshift/distance measurements. The raw candidate
  `matched_id` is retained as evidence but is **not** automatically adopted.
- An SHA-256 receipt establishes local content integrity, not remote source
  authenticity or astrometric truth. Timestamp/freshness does not imply
  scientific applicability.
