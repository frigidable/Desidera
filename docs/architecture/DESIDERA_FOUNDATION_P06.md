# DESIDERA P0.6 — scientific foundation before photometry

Status: **metadata-only initial implementation, Draft**. No code in this branch
claims measured surface brightness, validated CCD calibration, sky reconstruction
or physical source independence.

## Why this stage comes first

The project materials identify four common architectural subsystems:
Observation Registry, Reconstruction Engine, Uncertainty Engine and Scientific
Eligibility Engine. For every future measurement, we need an auditable input
identity and explicit statement of what the data can and cannot establish.
The literature audit also warns that coadd-versus-CCD agreement is not an
independent validation when both derive from the same exposures, and that
preprocessing may have removed diffuse flux before downstream modeling.

The development sequence must **not** imply that these engines are already
scientifically implemented. This PR makes the first narrow, verifiable
improvement in the Observation Registry and fixes interfaces for later work.

## Invariants (hard acceptance)

1. A galaxy has one stable project UUID; original FITS, old object passport
   (`lsbpipe-object.json`), `datasets/index.json`, previously saved diagnostics
   and archived runs are preserved.
2. A `DataSource` is a reduction/provider/release grouping. A `DataAsset`
   is an original file. A CCD observation is an independently identified
   physical detector exposure, not a coadd file, source group, brick or release.
3. Equal `archive + instrument + exposure_id + detector_id` means the
   *supplied catalogue assertions* refer to the same CCD, regardless of
   DR9/DR10/DR11. The detector-independent exposure group is kept separately.
   Exact identifiers are retained, including leading zeroes.
4. WCS coincidence, brick geometry, identical photometric bands and absence
   of a registered CCD match are never evidence of independent exposures.
5. A recorded catalogue row is only an *assertion* until its catalog, identity
   mapping and relation to the actual coadd have been independently verified.
   The P0.6 API deliberately cannot promote assertions to certified evidence.
6. Missing lineage means `UNKNOWN`, not "no exposure" or "independent".
   Same CCD or same exposure group is sufficient to flag possible reuse, but
   apparently disjoint identifiers do NOT prove independent noise/systematics.
7. Joint science coaddition remains disabled under
   `dataset_registry.compatibility_policy` regardless of this metadata.
8. All registry edits are offline, atomic, additive and idempotent.
   Invalid batches must not partially mutate the sidecar. No FITS downloads,
   conversion, science output or project reorganization occurs.

## First implementation: observation lineage sidecar

`src/lsbpipe/surveys/observation_lineage.py` adds:

- canonical `CCDIdentity` with stable content-derived CCD and exposure-group
  identifiers, independent of reduction/release;
- `datasets/observation-lineage.json`, schema
  `desidera.observation-lineage/1`, separate from the existing
  `lsbpipe.dataset-index/1`;
- `register_ccd_assertions(project_root, coadd_asset_id=...,
  catalogue_reference=..., records=[...])`;
- `compare_coadd_lineage(project_root, asset_a, asset_b)`;
- structural/foreign-key validation against actual registered native
  optical coadd IMAGE assets, no acceptance of invented missing assets.

Inputs are **explicit metadata**, *not* a claim that an official
`legacysurvey-<brick>-ccds.fits` parser has been written. Each record has
`archive`, `instrument`, exact string `exposure_id`, `detector_id`, exact
coadd `band`, and unique `record_locator` in `catalogue_reference`.
All contributions record `CATALOGUE_ASSERTION_NOT_VALIDATED`.
Source catalogue file hashing and official metadata inspection are next-step
work; arbitrary caller-supplied catalogue metadata must not become scientific
acceptance by merely passing a matching identifier.

`compare_coadd_lineage` returns one of:

| State | Conservative meaning |
|---|---|
| `SAME_ASSET` | Same file asset requested twice |
| `SHARED_CCD_IDS_ASSERTED` | Catalogue assertions point to a common detector image |
| `SHARED_EXPOSURE_IDS_ASSERTED` | Different CCDs share an exposure-group identifier |
| `UNKNOWN_LINEAGE_NOT_INDEXED` | One or both assets have no registered CCD identity |
| `NO_SHARED_REGISTERED_IDS_NOT_CERTIFIED` | No listed IDs match; completeness and independence remain unknown |

Every comparison returns `independence_certified=false` and
`joint_scientific_combination_allowed=false`. **This conservatism is
intentional**, not a limitation to be bypassed by adding a GUI toggle.

## How this fits the existing architecture

`src/lsbpipe/surveys/dataset_registry.py` remains the owner of data-source and
asset IDs, paths and checksum states. It must not be silently migrated. The
lineage sidecar uses its asset IDs and project UUID as foreign keys.

`src/lsbpipe/photometry/contracts.py` remains the read-only native coadd input
interface with PSF/sky/covariance evidence. It does not need a new photometry
algorithm merely to ingest CCD metadata.

The existing P0.3 `run_contracts.py` already has immutable run identifiers,
fingerprints and a limited DAG. Before the science engine it should acquire
versioned `AnalysisConfiguration` and `ScientificResult` contracts, retaining
distinct provenance and science-eligibility states.

The unmerged `feature/scene-background-engine` includes existing
experimental scene/mask/background/sensitivity components. Integrate only
after API/data-shape reconciliation with `NativeBandFrame`, independent
tests, and explicit non-certification of its uncalibrated background QA.

## Ordered follow-up increments, each requiring its own acceptance

| Phase | Code deliverable | Key acceptance |
|---|---|---|
| P0.6a | Schema/version/backcompat and branch integration | Old project and FITS bytes unchanged; full CI |
| **P0.6b (this PR)** | CCD/exposure identity and unverified relations | Same CCD across releases recognized; no false independence |
| P0.6b+ | Real offline Legacy brick-CCDs adapter and catalogue hash | Link only rows grounded in checked files, detect incomplete/ambiguous maps |
| P0.6c | SIMBAD/NED object identity + local catalogue cache | User confirmation for ambiguous matches, no automatic overwrite of target |
| P0.6d | Immutable analysis configuration, result type, selective DAG invalidation | Background-only change cannot redownload source FITS |
| P0.6e | Goal-based eligibility declarations and evidence states | Unknown remains unknown; refusal is explainable and reproducible |
| P0.6f | Bandpass/calibration metadata, WISE *catalogue* assets, scene adapters | WISE forced flux never treated as validated native optical FITS |

P0.6 should establish **interfaces** for reconstruction, identifiability and
uncertainty, not numerical confidence bounds before injection-recovery.
The benchmark plan must include tests of sky-removal stages, field-scale
surface-brightness limit `mu_lim(L)`, PSF wings, and exact catalogue/erratum
versions. No astrophysical accuracy claim follows from synthetically passing
metadata-only tests.

## Baselines and integration

The GUI that was manually verified free of top-left ghost controls is from
`fix/webengine-render-diagnostics-p056`, SHA
`7d42b0fa0ea4d434dc5acebb36077c011438d8c4`. This new branch starts
at exactly that commit. PR #19 and its ancestor Draft PRs are *not* merged to
main; do not silently merge the separate scene-background branch or rewrite
project metadata. The final integration sequence must first review ancestor
PRs and compatibility, then rebase/merge with full scientific + GUI CI.

**No existing Python distribution namespace, import path, on-disk schema,
original file or GUI workflow is renamed by P0.6b.**
