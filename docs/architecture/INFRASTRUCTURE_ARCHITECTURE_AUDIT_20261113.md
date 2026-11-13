# Desidera infrastructure and architecture audit — 2026-11-13

**Status:** source-backed implementation audit of the existing
`frigidable/Desidera` repository, branch
`feature/scene-wcs-editor-sc1`, Draft PR #84. Source baseline
`6036485323dfb16af7b0d82753b7093306719e8e`
(`0.2.1.dev3`); corrected candidate `0.2.1.dev4`.

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
`README.md`, `AGENTS.md`, `docs/architecture/**`,
`docs/ux/**` and `.github/workflows/**`. The additional
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

**NOT FIXED:** two application processes can still simultaneously
read the same registry, perform independent modifications and
replace it in arbitrary order. A multi-process transaction would
need one explicit lock protocol **around read–modify–write**,
collision handling, bounded retries and crash recovery. A single
snapshot merely avoids inconsistency *within one operation*.

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
