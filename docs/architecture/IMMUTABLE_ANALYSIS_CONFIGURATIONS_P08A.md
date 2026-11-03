# DESIDERA P0.8a — immutable analysis setup snapshots

**Status:** Draft, setup metadata only, NOT a scientific run or photometric acceptance.

## Data model and behavior

The Project Hub provides **Analysis configurations…**, a dialog to review saved configurations and **Freeze current setup…** after giving it a descriptive label.

A frozen configuration contains:
- Immutable project UUID and target ICRS position (when specified).
- Selected release / region / brick / band (when selected).
- Target aperture and sky-annulus radii.
- The dataset registry's selected original assets, their source IDs, kinds, and **historically recorded** SHA-256 verification states.
- A content digest of the setup plus the snapshot's own SHA-256 identity.

Every receipt lives at `<galaxy>/analyses/configurations/config-<sha256>.json`. Writing uses exclusive creation and never overwrites a prior configuration. The existing object manifest, dataset index, CCD lineage sidecar, science QA runs and original FITS are unchanged.

A read verifies the full content hash, stored setup hash, and immutable project UUID. The review action compares current **project and dataset registry metadata** to the stored setup. It returns `SETUP_METADATA_MATCHES_INPUTS_NOT_REHASHED` or `PROJECT_OR_DATASET_METADATA_CHANGED`. **It does not rehash original FITS, prove exposure independence, or mark scientific results valid.**

When no dataset or SHA-verified assets exist, the stored record explicitly states `SETUP_ONLY_NO_DATASET`, `NO_REGISTERED_ORIGINAL_INPUTS`, or `PARTIAL_OR_UNVERIFIED_INPUT_EVIDENCE`. Registered asset hashes may have been verified previously; they are never represented as verified during the freeze.

**Compatibility note:** P0.8a setup snapshots are not yet bound into `lsbpipe.science-run-contract/1` receipts. That requires a future explicit run-contract integration with code revision, verified input materialization and dependency records, and tests on actual science runs.

## Windows manual acceptance

From a clean checkout in PowerShell 7:

```powershell
cd C:\Users\frigi\Documents\PhD\Desidera
git fetch origin
git switch --track origin/feature/p08a-immutable-analysis-configurations
python -m pip install -e ".[exact,gui,knowledge,test]" pytest-qt
python -m pytest -q tests/test_galaxy_knowledge.py tests/test_analysis_configurations.py
python -m pytest -q tests/gui/test_galaxy_knowledge_gui.py tests/gui/test_analysis_configurations_gui.py
python -m lsbpipe.gui
```

Adapt the folder to your actual repository location. Do not switch branches with uncommitted local changes; commit or stash first.

**GUI pilot:**
1. Open any existing Galaxy Project. Without opening its science workspace, choose **Analysis configurations…** and **Freeze current setup…**; give the snapshot a short label.
2. Close and reopen the dialog; the entry persists and can be inspected offline.
3. Change an aperture radius or dataset selection through the existing GUI, save the project, and reopen the configuration dialog. The original snapshot remains unchanged, but comparison must report a different setup.
4. Open **Galaxy Knowledge · SIMBAD…**. Search on demand, compare candidates, confirm one; verify the project target coordinates have not moved. Opening either dialog offline must never fetch remote resources.
5. For data-backed configurations, first index and verify the original FITS using Dataset Manager. The snapshot then preserves the recorded SHA-256 receipts; this is not a fresh file-integrity check.

Inspect under `analyses/configurations/` and `knowledge/simbad/`. No folders containing FITS should be moved or rewritten by these operations.

## Acceptance gates

- Backend and Qt offline tests on Windows/Ubuntu; no network in CI.
- Content tampering and wrong-project records fail closed.
- Original project metadata and FITS remain byte-identical after creating, reading and comparing configurations.
- Frozen configs never trigger photometry or silently choose a dataset.
- Manual Windows acceptance on actual existing project is still required.

This is a bounded foundation for fully reproducible configurations, not a finished physical inference engine.
