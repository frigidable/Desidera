# DESIDERA — documentation and evidence

This index distinguishes **operational instructions**, **engineering/scientific contracts** and **dated historical evidence**. Many documents were written for the former LSBpipe name or feature branches that have since been integrated. Their old headings, checklists and PR references are not a live implementation inventory.

For current capabilities and runnable entry points, use the [top-level README](../README.md), `src/`, `pyproject.toml`, tests and [exact-commit Actions results](https://github.com/frigidable/Desidera/actions). A passing test validates its tested behavior, not adoption of an astronomical measurement.

## Living guidance

- [Contribution and scientific-safety rules](../AGENTS.md): accepted development boundaries, original FITS and scientific provenance.
- [Product brief](ux/PRODUCT_BRIEF.md): current object-centric researcher workflow and explicitly unadopted science.
- [Contribution workflow](ux/AUTONOMY_PLAYBOOK.md): selecting the correct branch, recording CI evidence and obtaining review.
- [UX design system](ux/DESIGN_SYSTEM.md), [screen states](ux/SCREEN_STATES.md) and [quality gates](ux/QUALITY_GATES.md): design objectives and acceptance criteria; screenshot outcomes require actual evidence.
- [Project roadmap #44](https://github.com/frigidable/Desidera/issues/44): current work packages and release/science gates. [PRs](https://github.com/frigidable/Desidera/pulls) show work that has not entered the default branch.

## Current integration checkpoint

- [Repository checkpoint, 2026-11-09](validation/REPOSITORY_CHECKPOINT_2026-11-09.md): merged PRs #39–#43, CI and science statuses, user-workflow simplification, and historical branch triage.
- [Empirical PSF core and broad-wing audit](qa/PSF_CORE_WINGS_VALIDATION.md): the preserved Leo P DR11 South g-band P16–P20 analysis, measured radius limitations and unresolved wide-angle wings.

## Scientific and implementation contracts

These explain constraints and intended behavior; re-check details against current source before relying on an old phase note.

- [Aperture geometry, signed signal and mask/coverage conventions](ARCHITECTURE_DECISIONS.md)
- [Native FITS/WCS input contract](architecture/PHOTOMETRY_INPUT_CONTRACT_P05.md)
- [Galaxy project identity and resume behavior](architecture/OBJECT_CENTRIC_PROJECT_HUB_DESIGN.md)
- [Dataset Manager](architecture/DATASET_MANAGER_P050.md), [integrity](architecture/DATASET_INTEGRITY_P02.md) and [project consistency](architecture/PROJECT_CONSISTENCY_P01.md)
- [Science run engine](architecture/SCIENCE_RUN_ENGINE_P03.md) and [immutable configurations](architecture/IMMUTABLE_ANALYSIS_CONFIGURATIONS_P08A.md)
- [GUI worker lifecycle](architecture/GUI_WORKER_LIFECYCLE_P04.md)
- [SIMBAD](architecture/GALAXY_KNOWLEDGE_SIMBAD_P07A.md), [NED](architecture/NED_OBJECTLOOKUP_P10B.md) and [catalogue read model](architecture/GALAXY_KNOWLEDGE_READ_MODEL_P10A.md)

Other `docs/architecture/` files contain feature-specific engineering records. They should not be read as proof that a feature is accepted science.

## Historical evidence — not current requirements

- [Original prototype v0.1](VALIDATION_REPORT.md), [v0.2 validation](VALIDATION_REPORT_V02.md), [technical design](TECHNICAL_DESIGN_V02.md): dated results and assumptions; old test counts are not present-day CI counts.
- [Real DR9 FITS audit](REAL_FITS_AUDIT_DR9_0401P012.md), [early FITS pilot](REAL_FITS_PILOT_P3.md), [real CCD-table acceptance](validation/REAL_LEGACY_CCD_ACCEPTANCE_P06C.md): inspect source hashes, exact samples and limitations. Real-data parser evidence is not photometric certification.
- [Scene masks](SCENE_MASKS_V03.md), [background](BACKGROUND_V03.md), [sensitivity](SENSITIVITY_V03.md), [coverage](COVERAGE_P2.md): experimental research contracts, not adopted photometry.
- [Survey Explorer v0.4](SURVEY_EXPLORER_V04.md), [redesign progress](ux/REDESIGN_PROGRESS.md), [past GUI defect audit](ux/GUI_NEXT_GATES_AUDIT.md): older UI iterations and observations; some recorded problems have since been fixed.
- [Draft consolidation audit](validation/DRAFT_CONSOLIDATION_20261103.md): pre-merge snapshot, preserved for its ancestry and scientific caveats rather than its old integration instructions.

The standalone `docs/ux/IMPLEMENTATION_PLAN.md` was a superseded plan tied to an older development branch and is retired from the working tree. Git history retains its content. Dated reports that record actual measurements, scientific caveats or experiment provenance remain in place.

## Source layout and reproducibility

| Location | Role |
| --- | --- |
| `src/desidera/` | Public package facade and current CLI/GUI launchers |
| `src/lsbpipe/` | Preserved numerical implementation and legacy imports |
| `src/lsbpipe/surveys/` | Datasets, FITS sources, projects, catalogue and lineage services |
| `src/lsbpipe/photometry/` | Experimental photometry data contracts |
| `src/lsbpipe/gui/` | PySide6 project, preview and scientific QA interface |
| `tests/`, `tests/gui/` | Scientific and offline GUI regression checks |
| `scripts/` | Synthetic fixtures, CI screenshots and explicitly invoked real-archive audits |

The scripts are not packaged research data. The original scientific FITS, local projects, screenshots and generated outputs are excluded from Git. No file deletion from this repository should be interpreted as permission to delete a local galaxy project.

For changes to science output, the scientific acceptance issues — [PSF #59](https://github.com/frigidable/Desidera/issues/59), [photometry #62](https://github.com/frigidable/Desidera/issues/62), [uncertainties #66](https://github.com/frigidable/Desidera/issues/66), [release #68](https://github.com/frigidable/Desidera/issues/68) — remain separate from software CI.
