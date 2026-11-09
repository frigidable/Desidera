# Desidera — integrated repository checkpoint (2026-11-09)

This file is a source-first, immutable-history-oriented checkpoint after consolidating PRs #39–#43. Its **code integration base** is merge commit `071d5c8ee89e06565f2377e9bcbec3e1b10f2c92`. The present documentation/regression-maintenance commit descends from that merge; Git itself records its final SHA. An engineering milestone must never be described as science validation.

## Integrated PR ancestry and exact merge SHAs

| PR | Topic | Merge commit |
| --- | --- | --- |
| [#39](https://github.com/frigidable/Desidera/pull/39) | Object-centred Photometry Workspace foundation | `754746b98472ec6dc82e88cd69320143615e490a` |
| [#40](https://github.com/frigidable/Desidera/pull/40) | Shared ICRS scene/mask geometry | `f2460d14340683bb2b345e11f51f737892b93562` |
| [#41](https://github.com/frigidable/Desidera/pull/41) | Researcher-seeded Smart Mask | `b1021a0982ef70821234d53deba3917ecdb4d126` |
| [#42](https://github.com/frigidable/Desidera/pull/42) | PSF-first scene / sky workflow dependencies | `600ec06576e3e76398421c283ca98cda6d4adf76` |
| [#43](https://github.com/frigidable/Desidera/pull/43) | Leo P empirical PSF, scientifically non-adopted P16–P20 diagnostics and simplified researcher UI | `071d5c8ee89e06565f2377e9bcbec3e1b10f2c92` |

All five PRs are **merged**, their original feature tips are ancestors of `main`, and there were **zero open pull requests** at this integration checkpoint. Merge commits preserve the original scientific commit history and QA notes.

## User-facing state (only actual implemented operations)

- Galaxy Project Hub: create/open/relink object projects and resume from the saved local source selection, preserving original project UUID and manifests.
- Source management: inspect/select/download original native Legacy Survey products by explicit user action; previews are separate from calibrated FITS.
- Photometry Workspace: select original native optical band, inspect source provenance, build empirical PSF, inspect the 2D kernel and observed/model radial profiles, review candidate stars, export explicit diagnostic data and QA, or run the available multiband batch sequentially.
- GUI cleanup completed in PR #43: removed eight **non-executable roadmap cards** and duplicative methodology panels; original six exploratory PSF recipes still exist for legacy QA but their controls are disclosed only in Advanced diagnostics. Standard defaults to Gaia/Tractor strict selection with a registered FWHM-scaled recipe. Missing Gaia/Tractor does **not** silently fall back to unverified DAO-only selection. One compatible band with native FITS is selected on entry regardless of the optional catalog's availability.
- API and disk compatibility preserved: `src/lsbpipe/` remains the scientific implementation; `src/desidera/` remains the public facade and launchers. No bulk rename and no changes to original local FITS, immutable source hashes, catalog identities, manual star-review receipts, or historical PSF QA records.

## Scientific closure versus remaining HOLD gates

The empirical PSF0 *engineering study* for `dr11/south/1554p180/g` has a frozen 81×81 native ePSF core, 10 original training stars, 3 core holdouts and 48 original catalog-only outer-profile diagnostic identities. The corrected P18 native WCS and spatially stratified P20 /2 background controls are recorded in [PSF_CORE_WINGS_VALIDATION.md](../qa/PSF_CORE_WINGS_VALIDATION.md). FWHM ≈5.137 native pixels; fully enclosed radial kernel support 40 native pixels (~7.79 FWHM); outer conditional residuals cannot uniquely identify a physically flux-normalized wide halo.

**`engineering_stage_finished=true` does not mean `science_adopted=true`.** The following remain **FALSE/HOLD** for the pilot: `psf_wings_validated`, `science_adopted`, spatial representativity, calibrated native-coadd noise covariance, stellar contamination separation, and wide-angle scattered-light correction. Galaxy injection/recovery is deferred to the separate photometry stage, by explicit research scope. No other brick was introduced to qualify the Leo P PSF.

The `g`-band pilot is not evidence that `r/i/z` models have been calibrated; the registered PSF0 implementation remains object-specific rather than a released all-galaxy general-science pipeline.

## Validation and evidence

Exact PR #43 head: `d8088a70c1585effc64e8b267ca33fb159862dc1`.

- [Scientific/Survey Explorer workflow run 42617056781](https://github.com/frigidable/Desidera/actions/runs/42617056781): **success** on Ubuntu and Windows / Python 3.10 and 3.12, plus GUI-import checks. Ubuntu 3.12: **721 passed, 42 skipped, 1 warning**.
- [Windows offline GUI workflow run 42617056697](https://github.com/frigidable/Desidera/actions/runs/42617056697): **success**, **736 passed, 3 skipped, 1 warning** in the scientific/UI mixed checks, and **201 GUI tests passed**; offline screenshots generated as CI artifacts. Screenshot creation is not equivalent to independent visual approval, and no scientific measurements are established from synthetic screenshots.
- The subsequent checkpoint/docs regression commit should be rechecked on its own GitHub `main` run when available. Do not reuse earlier run IDs as proof of that later SHA.
- No locally mounted real FITS were altered during these GitHub-only operations.

## Remote branch triage and deletion boundary

At this checkpoint, the GitHub repository still exposed **45 remote non-`main` branches**; this count is *not* the number of outstanding changes.

- **41 feature/fix/history branch tips**, including heads of PRs #39–#43, are now ancestry-merged into `main`. They are safe candidates for deletion **only after** verifying `ahead_by=0` on the latest remote main. Ref deletion is a separate operation, not implied by PR merging.
- **Four branch histories** still formally diverge at this checkpoint: `feature/galaxy-passport-projecthub-integration-20261103`, `maintenance/current-docs-prune-20261107`, `maintenance/repository-hygiene-20261107`, and `maintenance/source-first-audit-20261107`. Their ahead counts are 3, 1, 1 and 2 respectively. Recheck the live comparison before deleting any of these refs. The older passport branch contains an obsolete Project Hub implementation compared to `main`; merging it verbatim would regress the working application.
- Historic maintenance branch content has been reconciled: current README and docs index already incorporate the important source-first changes, current ignore rules are a superset of the retired hygiene branch, and two useful portable regression files (`tests/test_repository_docs.py` and `tests/test_desidera_entrypoints.py`) are incorporated by this checkpoint.
- Do not merge obsolete UI/doc snapshots only to make the branch counter reach zero. The remaining diverged refs should be explicitly retired/deleted after verifying that their unique changes are superseded. This session's GitHub connector can merge/update files, but exposes **no delete-branch action**; no remote refs were deleted or force-moved, and no local branch or directory was created.

## Next development baseline

Use the integrated `main` as the one source of truth; keep user-facing operations short and scientific QA accessible through reports rather than visible pilot controls. Preserve original FITS, project identity, scientific statuses and the complete P16–P20 QA log; never enable empirical PSF wing scattering correction by housekeeping. Later work on generalized stellar PSF, multi-object photometry, sky/covariance and galaxy injection-recovery requires separate scientific validation.
