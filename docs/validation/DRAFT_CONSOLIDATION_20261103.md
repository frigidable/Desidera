# Historical PR consolidation audit — 2026-11-03

**ARCHIVAL RECORD:** This describes an earlier pre-merge checkpoint and
its original release gates, not a pending instruction to merge those
historical PRs today. Re-check the actual `main` ancestry and open PR
bases before taking any integration action. The scientific caveats below
remain relevant.

**Status: integration candidate; do not claim a science-validated release merely because CI passes.**

## Source branches and ancestry

- Main development chain: open drafts **#3–#21, #23, #25, #26 and #27** (stacked by Git history; #24 was already integrated into #22).
- Experimental scene/mask/background/photometry sensitivity branch: **#2**.
- Compact right-side Galaxy Passport integrated with Project Hub: **#22**.
- Integration baseline began from `feature/p09a-simbad-mirror-failover`, exact SHA `66e27ac62c5f62cf162695cf7f21f6ecd2d6f251`.
- The #2 branch was incorporated as a true two-parent Git merge commit `1d4c118ccdab864a10ae2f03a186269c3b34db71`, with upstream head `ae0e9dff16a6b2a9cced3e5eac59e75e9dd021b9`.
- The #22 branch was incorporated as a true two-parent Git merge commit `5f04cd304d37c9a1f0b8f3c523532a2e2b1363a7`, with upstream head `4b15f8a748b9975033f6e5944d445e23bd85c69b`.
- Both merges preserve the full historical branch ancestry. No force-push, rebase, history rewrite or file deletion was used.
- `AGENTS.md` was reconciled by retaining the later DESIDERA GUI/science contract and adding the early scene/background scientific invariants.
- Overlapping `src/lsbpipe/gui/project_hub.py` was resolved in favor of the compact read-only GalaxyDetailCard while retaining visible SIMBAD and immutable-analysis actions and legacy widget API compatibility.
- A dedicated integration regression lives in `tests/gui/test_consolidated_project_hub_gui.py`.

## Release gates

1. Current combined commit must pass the full scientific suite with Python 3.10 and 3.12 on Ubuntu and Windows, plus Windows GUI import.
2. Windows Qt GUI test suite and offline screenshot generation must pass, with no regressions to Project Hub navigation or original FITS safety.
3. Integration must compare cleanly as a descendant of the historical `main`. All merge commit parent SHAs must be preserved and reachable.
4. Before touching `main`, require stable green CI on the integration PR head, exact-head merge protection, and no unresolved review blockers.
5. After the default branch is updated, verify the final `main` SHA, its CI run(s), and that all original draft branch heads are reachable from the final main tree.

## Scientific limitations that survive merging

A merged codebase is **not** a certified pipeline for science photometry, SBF, independent-exposure stacking, residual-background correction, published stellar populations, or catalog distances. PR #2 explicitly describes experimental residual sky and sensitivity diagnostics; P0.6c proves limited real-table reader compatibility for selected DR9/DR11 CCD catalogues, not pixel-level completeness. SIMBAD and Harvard failover were tested using offline mocks and require manual network validation on the user's Windows PC. Immutable configurations preserve *recorded* hashes but do not rehash original inputs when frozen.

No old scientific measurements are automatically promoted. The original FITS are not changed by merging code. Existing galaxy project folders and local session state are outside Git and should be backed up prior to the first local GUI run with the merged branch.

## Post-merge pilot

On Windows, save any existing local changes, then fetch and update a clean worktree from `main`. Run:

```powershell
python -m pip install -e ".[exact,gui,knowledge,test]" pytest-qt
python -m pytest -q
python -m pytest -q tests/gui
python -m lsbpipe.gui
```

Use a disposable test galaxy project first, checking Project Hub card selection, manual SIMBAD lookup, configuration freeze, the Dataset Manager, and actual original-FITS QA without any automatic scientific adoption.

The existence of these notes is not a substitute for the recorded GitHub Actions results.
