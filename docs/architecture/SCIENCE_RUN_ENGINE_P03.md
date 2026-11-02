# P0.3 — Project-scoped reproducible science runs

**Implementation:** Draft PR #6, `feature/science-qa-gui`.
**Scope:** deterministic run identity, trustworthy completed history,
explicit source and parameter staleness, limited dependency DAG.
Not a photometric computation engine, background model, or scheduler.

## Data and scientific identity

The project UUID in `lsbpipe-object.json` is the immutable owner.
Each GUI science action freezes `ScienceContext` with target ICRS,
release/region/brick/band, circular QA radii, full source-root paths
and `project_id`. The worker must recreate the exact same frozen
context. It rejects mismatches before publishing a result.

For object-bound workflows, run outputs live at
`<object project>/qa-runs/<run_id>/` regardless of whether original FITS
are under the project or in an external linked datastore.
Existing standalone/legacy runs retain their
`<original-data-parent>/qa-runs` layout. Previous receipts stay
readable. No old original or QA file is automatically overwritten.

Current receipts retain backward-compatible schema
`lsbpipe.gui-science-qa/1`, and add a `run_contract` field with schema
`lsbpipe.science-run-contract/1`. The contract contains:
- UUIDv4 `project_id` (nullable only for legacy/standalone workflows).
- Explicit release/region/brick/band dataset IDs; paired comparisons
  record distinct North and South identifiers.
- `parameter_sha256` (canonical JSON, SHA-256).
- `input_set_sha256` of the *content* fingerprints for every actual
  original IMAGE/INVVAR/NEXP/MASKBITS product. Exact file paths and
  filesystem mtimes are **not** hashed into this semantic signature;
  relocating a project alone does not invalidate its science.
- Git source revision if available in checked-out source or a trusted
  build variable; `null` when unavailable (not a fabricated commit).
- Explicit `dependencies` referencing immutable prior run IDs.
- Overall signature combining project, parameters, original content,
  source revision, and dependencies.
- Unchanged `DIAGNOSTIC_NOT_ADOPTED` science-status label.

`square` and `compare` use the target field but not circular object/
sky radii in their scientific parameters: changing a radius that their
algorithms do not use must not produce a false parameter-staleness
signal. `census`, `masks`, `figure` include those radii.

## Explicit review states

The check is read-only. It verifies completion manifests and the SHA-256
of saved artifacts. It also validates transitive referenced run
dependencies when present. Results are classified as:

| Status | Meaning |
| --- | --- |
| `RECORDED_NOT_REVALIDATED` | Receipt/artifacts valid, but original FITS have not been rehashed against current local inputs. |
| `CURRENT_INPUTS_AND_PARAMETERS` | Caller explicitly supplied current input hashes and selected science parameters; both match recorded signatures. Does **not** certify photometry or source integrity upstream. |
| `STALE_PARAMETERS` | Relevant algorithmic parameters changed. |
| `STALE_INPUTS` | The separately rehashed original FITS differ from recorded input hashes. |
| `PROJECT_MISMATCH` | Receipt belongs to another UUID. |
| `DEPENDENCY_UNAVAILABLE` | A parent run is absent, altered or linked to a different project, or the graph is cyclic. |
| `ARTIFACT_MISSING` / `ARTIFACT_MODIFIED` | A previously published output is absent or has a different SHA-256. |
| `LEGACY_UNVERIFIED` | A historic v1 manifest lacks the P0.3 run contract. |
| `INVALID_RECEIPT` | Corrupt, malformed, incomplete or unexpected receipt. |

**Recorded is not current**: reopening or listing results does not
rehash multi-GB science inputs by default. Checking inputs is an
explicit future user operation; rerunning the science QA produces a
new completed record. Results with another action or current
selection must never be promoted to the "latest accepted science".

## Dependencies and failures

A new run may declare earlier immutable run IDs, with no cross-project
references, missing parents, modified artifacts, cyclic graphs or
traversal paths permitted. Each dependency is rechecked before and
after the dependent computation. The DAG order routine independently
detects nonexistent nodes and cycles. A parent changing later
invalidates the child's review status.

Runs are staged in hidden `.partial` directories; only completed
results are atomically published. Historical incomplete directories
remain outside completed-history lookup and are never assumed
scientifically successful. Formerly published runs are never overwritten.

## UI

Science QA now contains **Review saved scientific runs…**. It
opens a read-only table after inspecting completed receipts in the
background, showing action, dataset, timestamp and review state.
Output manifests can be opened for inspection; no history action
automatically accepts or reruns a science measurement. The default
state of a well-formed run without new input verification is explicitly
`RECORDED_NOT_REVALIDATED`.

## Future gates

- Source rehash and version comparison as an explicit opt-in UI command
  with a progress indicator and reliable cached integrity metadata.
- Full persisted run DAG across data-integrity audits, sky models,
  mask generations, PSF models, photometry and uncertainties.
- Immutable source-revision capture for packaged releases with no Git
  working tree; environment pinning and dependency lockfiles.
- Cancellation/cooperative process lifetime and stable resumable jobs.
- Real original-FITS Windows pilot for UGC 2162, beyond synthetic CI.

This stage deliberately leaves numerical photometric processing out.
