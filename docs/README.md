# DESIDERA documentation and evidence index

**Source-first reading rule (2026-11-07):** Phase-specific documents in this
repository contain proposals, implementation notes and dated QA evidence.
Some were written for older branches or describe work subsequently changed.
They must not be treated as a maintained inventory of current capabilities.

For current functionality, start with `src/`, `pyproject.toml`, corresponding
`tests/` and the GitHub Actions result for the exact commit. `AGENTS.md`
defines editing and scientific-safety constraints. Older reports remain useful
as historical evidence of the inputs and tests they actually describe, not
as evidence of acceptance for a later code version.

The `main` branch represents integrated *code*, not necessarily accepted
astronomical measurements; stacked draft PRs are not part of that baseline.

## Start here

| Question | Source |
| --- | --- |
| How do I install, launch or test the current default branch? | [Repository README](../README.md) |
| What may change without science review? | [AGENTS.md](../AGENTS.md) |
| What are the current implementation priorities? | [GitHub roadmap — issue #44](https://github.com/frigidable/Desidera/issues/44) and its scoped work-package issues |
| Which changes are still experimental/unmerged? | [Open pull requests](https://github.com/frigidable/Desidera/pulls) — inspect each PR's **base branch** |
| What has actually passed in CI? | [GitHub Actions runs](https://github.com/frigidable/Desidera/actions) for the **exact commit** |
| What is the distinction between historical proposals and implemented functions? | This index, the corresponding source/tests, and the latest PR/issue evidence |

The active research gates are tracked in [PSF / Leo P #59](https://github.com/frigidable/Desidera/issues/59),
[measured photometry #62](https://github.com/frigidable/Desidera/issues/62),
[uncertainty/coverage #66](https://github.com/frigidable/Desidera/issues/66) and
[quantity-specific science release #68](https://github.com/frigidable/Desidera/issues/68).
A passing synthetic or offline Qt test does not satisfy these scientific
acceptance criteria by itself.

## Architecture and scientific invariants

- [Accepted prototype decisions](ARCHITECTURE_DECISIONS.md): signed flux,
  aperture geometry, mask/coverage distinctions and observational error limits.
- [Technical design v0.2](TECHNICAL_DESIGN_V02.md): early FITS/WCS integration.
- [Object-centric Project Hub](architecture/OBJECT_CENTRIC_PROJECT_HUB_DESIGN.md):
  project identity and resumed, portable data workflows.
- [Native FITS input contract](architecture/PHOTOMETRY_INPUT_CONTRACT_P05.md):
  band-specific data, WCS, validity and units.
- [Dataset Manager](architecture/DATASET_MANAGER_P050.md) and
  [dataset integrity](architecture/DATASET_INTEGRITY_P02.md):
  observations, local files, lineage and checks.
- [Science run engine](architecture/SCIENCE_RUN_ENGINE_P03.md) and
  [immutable analysis configurations](architecture/IMMUTABLE_ANALYSIS_CONFIGURATIONS_P08A.md):
  execution and reproducibility contracts.
- [GUI worker lifecycle](architecture/GUI_WORKER_LIFECYCLE_P04.md):
  asynchronous Qt operations without stale state.
- [Galaxy Knowledge / SIMBAD](architecture/GALAXY_KNOWLEDGE_SIMBAD_P07A.md),
  [NED lookup](architecture/NED_OBJECTLOOKUP_P10B.md) and
  [catalogue read model](architecture/GALAXY_KNOWLEDGE_READ_MODEL_P10A.md):
  external evidence, confirmation and offline handling.

The other phase-specific notes in [`architecture/`](architecture/) record
implementation decisions for project consistency, observation selection,
pixel-level coverage, remote failover, galaxy cards, catalogue links and
passport presentations. A versioned phase label is **not** a release or a
science-quality endorsement.

## Current GUI/UX engineering references

- [Product brief](ux/PRODUCT_BRIEF.md), [design system](ux/DESIGN_SYSTEM.md),
  [screen states](ux/SCREEN_STATES.md) and [quality gates](ux/QUALITY_GATES.md):
  read together before changing Qt screens.
- [Autonomy playbook](ux/AUTONOMY_PLAYBOOK.md): branch discipline, test evidence
  and review requirements.
- [Progressive research workflow](ux/PROGRESSIVE_RESEARCH_WORKFLOW_P043.md),
  [Science QA Workbench](ux/SCIENCE_QA_WORKBENCH.md) and
  [Survey Explorer v0.4](SURVEY_EXPLORER_V04.md): task-specific behavior.
- [Redesign progress](ux/REDESIGN_PROGRESS.md) and
  [GUI next gates](ux/GUI_NEXT_GATES_AUDIT.md): dated audit evidence, not
  necessarily unresolved current defects.
- [`ux/`](ux/) contains additional screen- and milestone-specific records,
  including dataset, mask, FITS display and card-preview pilots.

## Observational evidence and validation

- [Validation report v0.2](VALIDATION_REPORT_V02.md),
  [earlier report](VALIDATION_REPORT.md): dated, environment-specific
  test records.
- [Real-FITS DR9 audit](REAL_FITS_AUDIT_DR9_0401P012.md),
  [real-FITS pilot](REAL_FITS_PILOT_P3.md) and
  [real Legacy CCD archive acceptance](validation/REAL_LEGACY_CCD_ACCEPTANCE_P06C.md):
  inspect the exact inputs and acceptance scope before citing any result.
- [Background](BACKGROUND_V03.md), [scene masks](SCENE_MASKS_V03.md),
  [sensitivity](SENSITIVITY_V03.md) and [coverage](COVERAGE_P2.md):
  original diagnostic/experimental contracts and limits.
- [Draft consolidation audit](validation/DRAFT_CONSOLIDATION_20261103.md):
  **historical integration checkpoint**, not a current merge instruction.

## Code and test locations

| Path | Purpose |
| --- | --- |
| `src/desidera/` | Distribution name, entry points and compatibility-facing launchers |
| `src/lsbpipe/` | Scientific core and original compatibility package |
| `src/lsbpipe/surveys/` | Survey discovery, datasets, projects, lineage and knowledge services |
| `src/lsbpipe/photometry/` | Photometry contracts and native FITS adapters |
| `src/lsbpipe/gui/` | PySide6 views, controllers and workflows |
| `tests/`, `tests/gui/` | Scientific and offline GUI regression tests |
| `scripts/` | Synthetic fixtures, GUI captures and explicit archive acceptance |
| `.github/workflows/` | Cross-platform science tests, Windows GUI captures, optional archive acceptance |

From the existing checkout, the broad local regression commands are:

```powershell
python -m pip install -e ".[exact,gui,knowledge,test]" pytest-qt
python -m pytest -q tests --ignore=tests/gui
python -m pytest -q tests/gui
```

A passing run confirms only the cases exercised. The observational
acceptance gates in the roadmap remain separate.

## Editing policy

Keep existing document paths stable when possible: PR descriptions, issues,
and historical research notes link to them. Prefer correcting dated status
claims or adding links here over bulk renames, duplicated issue trackers,
new folder hierarchies or deleting historical scientific evidence. Update this
index when a genuinely new canonical document is introduced. Keep on-disk
galaxy project data, original FITS and generated QA products outside Git.
