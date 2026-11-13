# Desidera infrastructure and architecture audit — 2026-11-13

**Original audit:** source-backed implementation of the
`frigidable/Desidera` repository, originally in
`feature/scene-wcs-editor-sc1`, PR #84. Source baseline
`6036485323dfb16af7b0d82753b7093306719e8e`
(`0.2.1.dev3`), corrected engineering version `0.2.1.dev4`.
PR #84 was subsequently **merged into main** on 2026-11-13
as squash commit `07f139b69614ee3cc1a1780bbfa9f2fa54452e1f`.
Later corrections are developed in separate narrowly scoped PRs.

**Scope:** application boundaries, installable packaging, GitHub Actions
coverage and trust, GUI/science dependency separation, project-state
storage, filesystem safety, offline/online execution, long-term
maintainability and scientific provenance. This is NOT a scientific
validation, full supply-chain penetration test, performance benchmark,
or human GUI screenshot sign-off.

**Source of evidence:** exact repository tree, `pyproject.toml`,
the three existing workflows in `.github/workflows/`, the
`src/lsbpipe` and `src/desidera` packages, relevant unit/GUI tests,
and prior documented science/performance contracts. Claims about
remaining risks are distinguished from reproduced defects. Current
GitHub branch protection, external permissions, real workstation
concurrency and real-data acceptance are not established here.

## 1. Actual architecture and boundary contract

```text
desidera  (distribution and modern command facade)
  ├── desidera.cli / desidera.gui  → existing lsbpipe launchers
  └── lsbpipe                         legacy stable scientific API
       ├── geometry, profiles, survey, scene, background
       ├── surveys/
       │    ├── object_projects, local_project, project registry
       │    ├── dataset_registry, dataset_integrity
       │    ├── operation_journal, task_control, I/O and catalogues
       │    └── other survey and observation-lineage adapters
       ├── photometry/
       │    ├── native_fits, contracts
       │    ├── PSF/wing diagnostics, scene/mask and sky representations
       │    ├── aperture_flux, aperture_growth, aperture_systematics
       │    └── workflow_dependencies, provenance and QA
       └── gui/
            ├── project_hub, survey_explorer, dataset_manager
            ├── photometry_workspace, scene_editor, sky_workspace
            └── PSF panels, GUI orchestrators and view models
```

Desired dependency direction: **GUI → project/survey adapters →
validated native scientific data and pure algorithms**. No PySide6
or `lsbpipe.gui` import in the scientific/survey core.
UI-only preview pixels are not original FITS. Historical `lsbpipe`
imports, project schema IDs and original-FITS paths remain supported.
The new `tests/test_infrastructure_contracts.py` parses the core
Python AST to enforce the Qt-boundary rule in ordinary headless CI.

The project and scientific-result persistence roles remain distinct:

| Evidence | Existing storage | Claim |
| --- | --- | --- |
| Object passport | `lsbpipe-object.json` | Project UUID, target, selected dataset and links; **not photometric truth** |
| Legacy session | `lsbpipe-project.json` | Historical metadata preserved, not implicitly rewritten |
| Recent projects | User-local `projects.json` | UI registry, not an authoritative science record |
| Dataset registry | `datasets/index.json` | Imported assets and checksum/compatibility declarations; separate FITS QA |
| UI resume | `photometry-workspace-state.json` | Active tab/band; not a science receipt |
| Operation journal | `operations/*.json` | Attempted/completed UI jobs; no automatic science adoption |
| Native FITS | Original Legacy quartet | IMAGE/INVVAR/NEXP/MASKBITS; source data must be immutable |
| Scientific QA | Versioned diagnostic receipts and run artifacts | Conditional measured quantities and lineage; still not adopted calibration |

All of these contracts must remain valid before the data model is
renamed or a generalized executor replaces the existing controllers.

## 2. Confirmed findings and landed corrections

### A01 — GUI CI omitted major core changes — **FIXED / P0**

**Evidence:** `.github/workflows/ui-baseline.yml` executed the
ENTIRE scientific regression suite, GUI regressions and screenshots,
but its `pull_request.paths` listed selected subpackages/tests
rather than all Python sources. A change to `src/lsbpipe/background.py`,
`src/lsbpipe/scene.py` or a non-whitelisted scientific test did not
necessarily queue the Windows full-suite GUI workflow.

**Correction:** `ui-baseline.yml` and `survey-tests.yml` both
cover `src/**`, `tests/**`, `scripts/**`, `pyproject.toml`,
`README.md`, `AGENTS.md`, **all** `docs/**` (including
`docs/README.md`) and `.github/workflows/**`. The additional
documentation/governance triggers deliberately ensure that the
*final documentation-only HEAD commit* gets a new full CI verdict.
This costs CI minutes; a later change to event filtering must
maintain this exact-HEAD claim explicitly.

**Regression:** an offline infrastructure contract test verifies
coverage of the broad source/test and architecture directories.
The test itself does not simulate GitHub event routing; the
actual final workflows must complete successfully.

### A02 — unbounded concurrent GUI screenshot runs — **FIXED / P0**

**Evidence:** the scientific matrix already had per-PR
`cancel-in-progress: true`, but `ui-baseline.yml` did not.
Successive commits queued multiple Windows scientific+GUI+
screenshot workflows, consuming runner time and generating
stale evidence.

**Correction:** per-PR/job-family `concurrency` now cancels
superseded GUI runs. Science matrix jobs have
`timeout-minutes: 20`, GUI import has
`timeout-minutes: 10`, and the wheel smoke job 15.
The separate UI workflow retains its existing 20-minute bound.
A cancelled intermediate workflow is **not** a failure of
the FINAL SHA; final HEAD needs a completed successful verdict.

### A03 — archived real FITS acceptance downloaded on PR — **FIXED / P0**

**Evidence:** `legacy-ccd-real-archive.yml` has a
`pull_request` trigger and its original
`accept_legacy_ccd_archives.py` execution fetched real
public Legacy CCD archive data without an interactive opt-in.
External network availability, large downloads and upstream
archive changes can make a PR run costly or nondeterministic.
This conflicts with the requirement that original-data retrieval
be a deliberate, explicitly requested action.

**Correction:** on PR, the workflow still runs the deterministic
**offline** `tests/test_real_legacy_ccd_acceptance.py`.
The real archive download, evidence printer and evidence upload
are gated on `github.event_name == 'workflow_dispatch'`.
That manual action still performs the genuine archive checks,
and its evidence is preserved even on failure. No official
remote data are claimed downloaded in this audit.

### A04 — floating third-party Action references — **FIXED / P1**

All three GitHub Actions files used moving major tags:
`checkout@v4`, `setup-python@v5`, `upload-artifact@v4`.
The live GitHub public tag references were checked on
2026-11-13 and the resolved **40-character commit SHAs**
pinned in workflows (with human-readable tag/date comments):

- `actions/checkout`: `11d5960a326750d5838078e36cf38b85af677262`
- `actions/setup-python`: `a26af69be951a213d495a4c3e4e4022e16d87065`
- `actions/upload-artifact`: `ea165f8d65b6e75b540449e92b4886f43607fa02`

A source-level regression rejects short/floating `uses`
references in these workflows. This reduces mutable-tag drift,
**not** broader PyPI/build-backend/supplier compromise risk.
Action commits need periodic security review and planned updates.

### A05 — editable-only packaging validation — **FIXED / P1**

**Evidence:** the original CI installed `pip install -e ...`
and tested imported sources in the checkout, but did not
build/install the distribution wheel. The existence of
`pyproject.toml` scripts did not prove that a noneditable
installation includes both `desidera` and `lsbpipe`
packages with all five historical/new CLI entry points.

**Correction:** `survey-tests.yml` gains an independent
Ubuntu/Python 3.12 `package-wheel` job. It builds a real
wheel, installs that wheel with declared runtime dependencies
and executes a metadata/import/entrypoint check from the
runner TEMP directory, outside the source tree.
It verifies exactly the public, backward-compatible APIs
without launching optional GUI or claiming science success.
`pip check` is also run in the scientific matrix,
GUI workflow, import smoke and wheel installation.

**Limit:** wheel smoke is not yet a full clean-machine
Windows/macOS desktop installer, sdist acceptance,
standalone binary or security code-signing system.
No release tag is created by this work.

### A06 — implicit GUI test dependency — **FIXED / P1**

**Evidence:** `ui-baseline.yml` installed `pytest-qt`
as an untracked extra command-line package instead of a
declared project optional dependency; README repeated the
manual workaround.

**Correction:** `pyproject.toml` now declares
`test-gui = ["pytest-qt>=4.4"]`, GUI CI installs
`.[exact,gui,test,test-gui]` and the README documents
the same public setup. The infrastructure test checks
the *installed* distribution's optional dependency metadata.

### A07 — portable project passport could follow symlink escape — **FIXED / P0**

**Evidence:** `read_object_project` already rejected lexical
`..\` traversal in relative `data_root` paths but
joined `folder / relative_path` without checking whether
intermediate symlinks/junctions resolve OUTSIDE the project.
Unlike `datasets/index.json` and the UI resume pointer,
`lsbpipe-object.json` and the user-local `projects.json`
were also not explicitly refused when symlinked.
This could import metadata from unexpected locations or
misrepresent a supposedly portable in-project data path.

**Correction:**
- refuse symlinked object-passport path on both read and save;
- refuse a symlinked recent-project registry before reading;
- require `Path(relative_data_root).resolve()` to remain inside
  the resolved project root, including Windows junction behavior;
- retain explicit absolute external data-root references and
  legitimate relative `.` and ordinary `data/` references.

These are **metadata-path validations**, not an assertion that
arbitrary disk-content mutation between checks is impossible.
On platforms without symlink permissions, those specific
negative tests are intentionally skipped.

### A08 — two reads within one project-library update — **PARTIALLY FIXED / P1**

**Evidence:** `register_object_project` initially loaded the
registry's `projects` list, transformed it, and loaded
the full payload **again** for dismissals and publication.
If the registry changed between those reads, the procedure
could mix a stale `projects` list with newer other state.

**Correction:** use a single validated `_registry_payload()`
snapshot per registration. A negative regression instruments
and demands exactly one read. Existing atomic
`os.replace` file publication is retained.

**Historical limitation (fixed separately in follow-up branch):**
at the time of PR #84, two application processes could still read
the same registry, perform independent modifications and overwrite
one another. The new `fix/project-library-concurrency` change
introduces a stable OS-level advisory lock covering read–modify–write
for this **one user-local projects registry**; read its separate PR
and exact-commit CI before treating that follow-up as accepted.
Dataset indices, project passports and GUI resume/operation journals
are **not** automatically transactional merely because this registry
uses a lock.

### A09 — untested erosion of GUI/science separation — **GUARDED / P1**

The scientific-core architecture already intended a headless
library. A new AST-level regression scans actual
`src/lsbpipe/photometry/**/*.py`, `surveys/**/*.py`,
`background.py`, `scene.py`, `survey.py` and
`mask_engine.py`. It rejects direct imports of
`PySide6`, `PyQt5/6`, `qtpy`, `lsbpipe.gui` and
`desidera.gui`.

This is a **static boundary** for direct Python imports,
not a complete transitive runtime import graph,
plugin-entrypoint audit or a proof that no subprocess
could invoke a GUI. The existing `desidera` facade
and `lsbpipe` imports intentionally remain stable.

## 3. Open architecture and infrastructure risks

The following are **identified areas requiring follow-up**.
No completion, quantified speedup or science acceptance is
claimed merely from this review.

| Priority | Risk/area | Impact and proposed next acceptance criterion |
| --- | --- | --- |
| **P0** | Cross-process project/registry/analysis-index race | `_atomic_json` ensures per-file atomic publication but is **not** a transaction across two GUI windows; add an advisory cross-platform lock (or verifiable compare-and-swap) around the whole read–modify–write, conflict rejection, then multiprocessing race/injected crash tests |
| **P0** | Real-CCD and original-FITS provenance | Manual archive workflow is separated from PR CI, but true native input and independent-exposure acceptance still need actual source-hash and astronomical QA; no diagnostic may become adopted merely by passing infrastructure checks |
| **P1** | Large GUI controllers | Source-tree inspection finds `scene_editor.py` ~202 KB, `survey_explorer.py` ~140 KB and `psf0_panel.py` ~115 KB. Large multi-responsibility Qt classes can obscure worker lifetime, source identity and scene invalidation. Refactor incrementally behind unchanged public classes with state-machine/keyboard/resize screenshots and no new application directories |
| **P1** | Scientific module complexity | `psf_empirical.py` ~96 KB, `psf_wing_profile.py` ~72 KB, and `aperture_systematics.py` ~59 KB; extracting pure computations and data validation could improve test isolation, but must not alter published QA schemas or numerical integration |
| **P1** | Data/analysis schema migration governance | Coexisting legacy `lsbpipe-project.json`, object-passport, dataset index, resume pointer and multiple science receipts require versioned migrations, canonical ownership and explicit rollback/dry-run tests; never silently rewrite old originals |
| **P1** | Python/astronomy dependency drift | `numpy>=1.24`, `scipy>=1.10`, `astropy>=6.0`, `photutils>=2.0` and Qt lower bounds are intentionally broad; future major releases may alter photometry or ABI. Current CI tests only Python 3.10/3.12 on Ubuntu/Windows; add a reviewed constraints/lock policy and record actual installed versions per reproducible science run |
| **P1** | End-to-end native FITS resource budgets | Large coadds, PSF stamps, sky masks and Qt display buffers share memory. Current kernel tests do not bound process peak RSS, concurrency, cancellation latency or response to disk-full; a representative cold/warm observational benchmark remains required |
| **P1** | Supply-chain and release attestation | GitHub actions now pinned, but PyPI dependencies are not fully locked and there is no signed/reproducible official release artifact. A later release policy should include dependency-review, wheel/sdist provenance, reviewed tags and checksums |
| **P1** | Distribution governance / license | Repository root has no explicit `LICENSE` file in the inspected tree. Do not select an open-source legal license without the owner's decision. Public visibility is not a license grant. Formal changelog/release notes and ownership review would help publication |
| **P2** | Operation journal growth | `interrupted_operations` bounds returned results to 100 but still scans the journal directory; a long-running project could accrue thousands of job files. Add non-destructive indexed/bounded-read behavior before introducing automatic cleanup |
| **P2** | GUI source inventories and rendering | Qt `QListWidget` still rebuilds on filters; no measurements of large-catalogue UI p95, screenshot visual approval or widget accessibility are recorded |
| **P2** | PR reviewability / historical architecture drift | PR #84 contains many science/GUI/performance changes under an old SC1 title. A reviewer should use feature-specific commits and current audit indexes; new branch proliferation is not a substitute for explicit subfeature/QA sign-off. Dated architectural documents are not automatically current contracts |
| **P2** | Platform coverage | Windows+Ubuntu scientific tests and Windows Qt are tested. macOS and official Windows installation outside CI remain unverified; native survey FITS and OS-specific filesystem behavior require targeted validation |
| **P2** | Structural static checks | AST guard prevents direct Qt imports, not dynamic plugin import, forbidden Qt usage through subprocesses or hidden circular dependencies. Future tooling could include a module dependency graph and packaging-only import profiling |

## 4. Release and acceptance gates

**Engineering acceptance for this specific change** requires:

1. Every workflow file parses and GitHub accepts its triggers,
   pinned action SHAs and job permissions.
2. Exact-HEAD `survey-tests.yml`: all four combinations
   Windows/Ubuntu × Python 3.10/3.12 pass, GUI import job passes,
   and the **new installed wheel** job passes.
3. Exact-HEAD `ui-baseline.yml`: scientific suite, GUI suite,
   screenshots and artifact upload complete successfully.
4. Explicit `test-gui` installation and source-level
   infrastructure tests pass. Symlink/relative-traversal tests
   pass where filesystem permissions allow, and leave all
   external-test original bytes unchanged.
5. The archive acceptance PR path runs **only offline tests**.
   No real original FITS are downloaded in routine PR CI.
6. `main`, original native FITS and existing scientific
   receipts are not overwritten; PR remains Draft.
7. Green CI remains **engineering regression evidence**:
   no science adoption of sky, PSF wings, full photometric
   magnitude, coadd covariance, SBF distance or galaxy
   source membership.

**Observed CI evidence for implementation commit
`f4189b5ae6b53fc4070c0cb1e391be4f4f34676e`:**

- Scientific matrix on Ubuntu/Windows × Python 3.10/3.12,
  GUI-import and the **new noneditable wheel job**: [all
  successful](https://github.com/frigidable/Desidera/actions/runs/43503182675).
- Windows full scientific/GUI workflow: [successful](https://github.com/frigidable/Desidera/actions/runs/43503182660):
  **963 scientific tests passed**, 3 skipped and one expected
  `NoDetectionsWarning`; **264 GUI tests passed**; 31 general
  and 25 photometry screenshots produced (not manually
  reviewed for visual acceptance).
- Real CCD archive acceptance workflow: [successful offline
  validation](https://github.com/frigidable/Desidera/actions/runs/43503182650).
  The actual remote FITS download, archive printout and
  artifact upload were explicitly **SKIPPED** on the PR event.

These results establish correctness only for that exact
implementation commit, not the eventual distribution release.
The final PR description and GitHub checks identify the
latest HEAD, including any subsequent documentation-only
commits. In no case does green CI certify astrophysical
sky/PSF/SBF science adoption.

## 5. Practical prioritized next phase

Resolve cross-process registry and scene-state writes first,
with one narrowly scoped lock/transaction contract, a
Windows-compatible implementation, crash-injection tests and
zero mutation of original FITS. In parallel, instrument user-
visible GUI/worker performance and capture actual real-data
stage memory budgets. Only then plan module decomposition
behind existing interfaces, reviewed schema migrations,
dependency constraints and eventual license/release policy.

Do **not** automatically restructure the large existing GUI
classes, rename `lsbpipe` APIs, start another feature branch,
merge to `main`, download additional telescope FITS, or
change scientific calibration from this infrastructure audit.


## 6. Follow-up: object project library race — separate focused PR

The original A08 registry **single-snapshot** repair prevented
intra-call mixing of payload generations but not cross-process
lost updates. In `fix/project-library-concurrency`, a
`project_locking.exclusive_metadata_lock` context manager now
holds one stable advisory file lock through the entire
user-local project-library transaction:

- Linux/Unix uses `fcntl.flock(LOCK_EX | LOCK_NB)`;
  Windows uses `msvcrt.locking(LK_NBLCK, 1)`, with a bounded
  retry timeout. Both are OS locks released after process death.
- The `projects.json.lock` sidecar remains on disk to prevent
  accidentally splitting the lock domain between waiting processes.
  Its bytes are a non-scientific sentinel, **not** a project index.
- `create_object_project` protects duplicate-name checks,
  passport creation and first registration under one library
  transaction, using the inner already-locked registry writer
  to avoid nested non-reentrant OS locking.
- `register_object_project` and `forget_recent_project` now
  serialize reading, validation, modification and atomic replacement.
  Reads by `recent_projects` remain lock-free and safe against
  *incomplete JSON publication* because the writer still uses
  same-directory `os.replace`. A reader may see the prior or
  new complete generation.
- Reject symlinked/nonregular lock paths and an inode/path swap
  observed while acquiring the lock. No source FITS are opened,
  renamed, downloaded or otherwise modified by the mechanism.
- Regression tests invoke **actual different Python processes**,
  force a paused stale reader while another writer tries to update,
  assert that both project UUIDs survive, test lock-timeout and
  OS cleanup on abrupt `os._exit`, and reject a symlink lock.

**Limits:** the protocol coordinates only cooperating local
Desidera processes on filesystems supporting the stated OS locks.
It is not a distributed/SMB/NFS concurrency proof. Multiple
metadata files still lack a unified project-level transaction;
a process crash after writing an object passport but before
registering it may require ordinary explicit project import.
The dataset/analysis registry and operation journals have their
own separate update paths and remain subject to a further
targeted audit. The scientific HOLD and existing data schemas
are unchanged.


## 7. Dataset-index transaction follow-up — separate focused PR

This audit identified a second real **lost update** channel in
`surveys/dataset_registry.py`. The index is one atomic
`datasets/index.json` sidecar per galaxy; an atomic JSON replace
only prevents partial files, **not** two writers reading the same
prior generation and both publishing different derived snapshots.

**Affected paths before the repair:** `import_legacy_fits`,
`register_external_fits`, and `verify_asset_checksum`.
The last one was particularly problematic: it read a whole
registry generation, hashed potentially large source FITS files,
then saved the **original stale registry**, erasing any new
asset registrations made while the hash was running.

The narrowly scoped `fix/dataset-index-transactions` branch:
- reuses PR #85's cross-platform OS advisory lock in a stable
  project-local `datasets/index.json.lock` sidecar;
- serializes the full read–modify–atomic-replace transaction for
  FITS discovery imports and user-selected external FITS assets;
- performs full source-file SHA-256 **outside** the index lock,
  then acquires it, **re-reads the current index**, and verifies
  asset ID, source, path reference, size, mtime, product type,
  band, import identity, previous checksum and current on-disk
  file metadata before committing the original SHA-256 receipt;
- rejects stale or changed metadata rather than publishing a
  valid hash over another file's current identity;
- preserves any concurrently imported unrelated sources and
  assets; retains the previous `CHECKSUM_VERIFIED` and
  FITS-unvalidated semantics rather than asserting scientific
  image validity; and
- supplies actual separate-process races (concurrent external
  imports, import during a deliberately slow SHA-256), a
  concurrent metadata-mutation refusal and a symlinked-lock
  negative test.

**Scope limits:** neither project-passport nor GUI resume/journal
transactions are covered by this lock. Third-party software
that disregards the sidecar lock can still rewrite
`index.json`; no general distributed/SMB/NFS protocol is claimed.
SHA validation detects ordinary timestamp/size changes around
the read/commit but is not a proof against malicious filesystems
that replace content while preserving those metadata fields.
No FITS is modified, merged, repixelized or scientifically
adopted. Acceptance requires exact-HEAD Windows+Ubuntu CI and
the deliberately interleaved process tests before merge.


## 8. Legacy local-session optimistic concurrency — 2026-11-14 follow-up

**Confirmed before correction:** the historical
`lsbpipe-project.json` session retained the selected dataset,
photometric-screening radii, the target field and a pointer to the
last GUI QA artifact. `save_local_project` checked only the
legacy schema then atomically replaced the whole JSON. Two separate
Desidera instances could both open a valid original-FITS project;
a later save of one stale in-memory `LocalProject` silently
overwrote the other instance's dataset choice, QA pointer or
screening parameters. This was independent of the already-fixed
`projects.json` library, `datasets/index.json` and
`operations/*.json` transaction locks.

**Correction (focused PR `fix/legacy-project-session-integrity`):**

- A private `LocalProject.metadata_sha256` / `metadata_observed`
  pair records the SHA-256 of **the exact historical JSON bytes**
  read at `open_local_project`. These private fields are never
  written into the historical project schema and are **not**
  an original-FITS checksum or scientific QA certificate.
- `save_local_project` holds one stable cross-platform
  `lsbpipe-project.json.lock` while reading the current
  file, verifying both its schema and the originally observed
  exact revision, and atomically publishing a replacement.
  If another process changed the file, the stale writer
  raises an actionable `RuntimeError` instructing a
  *reopen* rather than silently discarding another window's
  changes. A new session observed as initially absent may
  create its JSON once. A directly constructed, never-opened
  `LocalProject` cannot overwrite an existing session.
- After successful publication, the in-memory object receives
  the new content hash and may be saved again without a reopen.
  This protects legitimate repeated edits by the same window.
- Legacy `lsbpipe-project.json`, `download-report.json`
  and historical QA manifest readers now reject metadata
  *file symlinks* rather than following externally redirected
  JSON. The user-local `last-project.json` hint also declines
  to read a symlink or a pointer under a symlinked settings
  directory; its writer rejects symlinked/nonregular paths.
- Regression checks two independently opened projects,
  stale-save rejection without overwriting the successful
  edit, same-instance repeated save, exact byte isolation of
  placeholder original FITS, and a deterministic **two-process**
  stale-session interleaving. Negative symlink tests verify
  that external metadata and the legacy last-session pointer
  remain unchanged. These tests run offline.

**Limitations:** a cross-file distributed transaction still
does not exist. The legacy user-local `last-project.json`
pointer is deliberately an independent last-active-session
hint and can lag behind a successfully written project JSON if
the process exits between those operations. No data-loss
recovery or automatic merge of concurrent edits is promised;
users must reopen the project when a revision conflict is
detected. The advisory lock only coordinates cooperating
local processes on a filesystem with correct OS lock
semantics; third-party writers that ignore it are not
excluded, although changed file bytes are normally detected
before a cooperative save. No original native FITS are
rewritten, and `SCIENCE HOLD` remains.


## 9. Galaxy object-passport stale-state protection — 2026-11-14 follow-up

**Affected source:** `surveys/object_projects.py` writes the authoritative
object-centric `lsbpipe-object.json`, which contains the immutable
galaxy project UUID, display name, target field, dataset selection,
aliases, confirmed/unconfirmed catalogue associations, latest QA
references and screening radii. PR #85 fixed the *user-local
project-library index*, but it did not serialize or version this
distinct per-galaxy passport. The previously atomic, same-directory
`os.replace` could still silently publish stale in-memory fields
from a second GUI instance with the **same valid UUID**.

**Focused correction:** `fix/object-passport-stale-writes` adds
private, nonserialized `ObjectProject.manifest_sha256` and
`manifest_observed` revision fields. The read path computes SHA-256
of the **same file bytes** from which the passport was decoded.
Before any save, the writer acquires a stable per-galaxy
`lsbpipe-object.json.lock`, rechecks the on-disk schema, UUID and
exact observed byte hash, then atomically publishes the new record.
It updates its own revision token on success, allowing repeated saves
from the same GUI instance. A stale second writer or a caller that
constructed a record without ever observing the existing passport
must explicitly **reopen** rather than overwrite the first writer.

**Regression requirements:** independent Python processes open the
same passport before either writes, and only the first committed
writer's target/project fields survive. A second test checks that
stale catalogue metadata cannot erase changed selected datasets or
QA pointers. A negative test refuses a reconstructed record that
knows the correct UUID but lacks an authenticated manifest revision.
No original FITS pixels, project schema keys, absolute/relative FITS
links or scientific QA acceptance flags are modified.

**Limits and follow-up:** read revision tokens are optimistic
concurrency checks against cooperative local writers, not
cryptographic guarantees of a trusted machine, transactional
multi-file snapshots or independent astronomy validation. GUI
resume pointers still intentionally use last-writer-wins semantics
for active-tab selection and need no lost-update transaction
unless their role changes. The highest-priority remaining
cross-file issue is crash recovery when *several* different
metadata stores must be published together (passport, dataset
index, journal and GUI resume). A future design should specify a
recoverable operation protocol with a single unambiguous owner
for each source-of-truth field instead of adding hidden locks
around the entire application. `SCIENCE HOLD` persists.
