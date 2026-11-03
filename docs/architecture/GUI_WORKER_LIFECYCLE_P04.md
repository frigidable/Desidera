# P0.4 — Qt GUI lifecycle and worker hardening

Status: incremental hardening in Draft PR #6, not completion of all GUI
refactoring or a Windows real-FITS acceptance run.

## Independent task-execution boundary

The former inline `TaskWorker` and `_TaskDispatch` have moved from
`survey_explorer.py` to `gui/task_execution.py`. Scientific numerical
functions do not import Qt. The work function receives a progress callback
and is run on a QThread. The GUI-owned `TaskDispatch` has decorated
`@Slot` methods for progress, success, error and stop; **never** connect a
plain Python lambda that manipulates widgets directly to worker-thread
signals. Signals use explicit `Qt.QueuedConnection` to reach GUI-owned
QObject methods.

Each task receives a `TaskToken` with a thread-safe stop event. A token
is tied to the owning `SurveyExplorer` and is checked before invoking
GUI callbacks. Superseded/closed-window callbacks may not change the
new project's GUI. Only a single job is allowed per workspace.

## Request stop

The workspace contains **Request stop** beside the progress bar. Clicking
sets a cooperative request; it does **not** call `QThread.terminate` or
kill a native thread. A stop is observed at the next wrapped progress
checkpoint and again when the work function returns. Therefore:
- Search/download/FITS work with frequent progress callbacks can finish
  earlier, but operations blocking on network or a FITS read may still
  finish their current internal step before observing the stop.
- Scientific functions may publish their atomic, immutable QA receipt
  before their top-level call returns. **Request stop does not roll back
  results or promise no files were written.**
- A late success after the stop request is discarded from the GUI. It
  is not presented as scientifically accepted, and previously completed
  diagnostic reports remain inspectable from the run history.
- A cancelled run emits no ordinary failure popup; actual worker failures
  still use the GUI-thread error pathway. One job may be started after the
  previous QThread has actually finished and the UI is idle again.
- Close/switch remains disabled while a QThread is still alive. There is
  no forced exit and no dead worker deletion.

Cancellation is currently **best-effort cooperative**, not a proof that
all numerical algorithms can interrupt at fine granularity. Full
inner-loop cancellation requires explicit checkpoints in each future
algorithm (sky models, fitting iterations, source detection, etc.).

## FITS rendering lifetime

FITS Explorer uses a bounded raster cache and background NumPy conversion
with generation tokens. It now exposes an idempotent
`shutdown_renderer()`, which rejects future `set_dataset`/render
requests and stops its QTimer before disposing its
`ThreadPoolExecutor`. The owning `SurveyExplorer.closeEvent` calls
this explicitly on accepted close. Merely closing a QMainWindow does not
guarantee every child QWidget receives its own `closeEvent`.

Running native NumPy work may continue briefly after
`shutdown(wait=False)`, but its result cannot update a deleted Qt
widget because timer polling and the renderer generation have been
invalidated. This shutdown path does not silently downsample original
pixels or alter FITS.

## Regression evidence

Real Qt GUI tests check:
- Work executed on QThread; progress/completion/error/stop callbacks
  executed via QObject slots on the GUI thread.
- Stop before/during work; cancellation at checkpoints.
- A non-cooperative work function completing after stop must not
  publish a stale success, and the next task can be started normally.
- Closing/returning to the Project Hub while a worker runs is blocked.
- Project switching after safe completion tears down the child FITS
  display thread and opens an independent fresh renderer.
- Renderer shutdown is idempotent, stops timer polling and refuses reuse.
- Existing local science QA, Project Hub, FITS raster/cutout regression
  and scientific test matrix remain unchanged.

## Remaining P0.4 work

- Refactor the 2000+ line `SurveyExplorer` into targeted subcontrollers
  with a documented state-transition graph; preserve the tested CLI/GUI API
  and do not change numerical photometry.
- Add proper progress/cancel checkpoints inside costly FITS operations,
  download streaming, and validation, without corrupting atomic receipts.
- Add memory/RSS and QWebEngine repeated-open/close stress tests on actual
  Windows high-DPI systems and on-device original DR9 FITS.
- Audit cross-project WebEngine JavaScript callbacks, retained Python
  references and shutdown after unexpected exceptions.
- Implement persistent crash-recovery markers for incomplete jobs,
  while never interpreting a `.partial` directory as completed science.
- Extend source-wide static checks, GUI automated accessibility checks and
  packaging across supported Python/Qt versions.

Do not describe P0.4 as a completed high-stress manual Windows acceptance
until the original-FITS test suite is actually run on the user's machine.

## P0.4.1 — hardening follow-up

The previous single `SurveyExplorer._start` implementation was replaced with
thin compatibility wrappers and a dedicated
`gui.task_execution.WorkspaceTaskController`. The controller owns thread
creation, thread shutdown, start/stop controls, GUI dispatch, and the
object-project operation journal. Legacy test hooks
`_start(action, finished, on_error)`, `_request_stop()` and `_idle()`
remain compatible.

A subtle asynchronous race is now explicitly guarded: delivery of
`worker.completed/failed/cancelled` and delivery of `QThread.finished` are
separate queued events. The controller and `TaskDispatch` hold a
**two-event barrier**. The active token, worker, journal and QThread wrapper
are not released or deleted until *both* worker terminal callback and
thread exit have been observed on the Qt GUI thread. Tests exercise both
delivery orders and ensure the end callback runs exactly once.

The science and survey libraries now share a GUI-independent
`surveys.task_control.OperationCancelled` and optional
`cancel_check` callback. `TaskWorker` makes a checkpoint available to
the work closure independently of GUI progress emission. Long operations
that support fine-grained checkpoints include:
- Survey HTTP transfer: before starting a file, between 1-MiB chunks,
  during SHA-256 of an existing file, and before atomic rename.
- Multi-file downloads: between products and before publishing the batch
  receipt. Cancellation does not turn into an ordinary per-file error.
  Complete files remain; `.part` files are cleaned. A cancelled batch
  is not falsely recorded as complete.
- Dataset integrity: between 256-row raster slabs, between brick footprint
  evaluations, during original FITS hashing and before audit publication.
  An explicit stop is propagated rather than classified as `INVALID_FITS`.
- HTTP HEAD verification: before each product and after a finished request.
  A blocked HTTP request still runs until timeout or underlying return.

### Recoverable operation journal

Every new operation started inside an object project creates a separate,
atomic `operations/<run-id>.json` record with
`lsbpipe.operation-journal/1`, UUID, action and `RUNNING` state.
Requesting stop updates it to `STOP_REQUESTED`. After terminal callback
delivery and QThread exit the state becomes `COMPLETED`, `CANCELLED`,
`FAILED` or `IGNORED`. Journal states are **not scientific results**:
all records carry `NOT_SCIENTIFIC_EVIDENCE`.

On reopening the project, the hub checks for journals still in
`RUNNING` or `STOP_REQUESTED` state. It shows a nonmodal warning that
the previous process may have crashed. There is no blind replay, deletion,
source FITS mutation or automatic acceptance of potentially existing QA
output. New jobs get independent journal records; uncertain old records
remain available for a later explicit recovery audit.

### Acceptance boundaries

CI exercises interrupted and normally completed journals, foreign UUIDs,
cancelled downloads, unchanged original FITS, real Qt signals, both queued
terminal/exit orders and project reopen with a warning. Production-grade
cancellation of numerical kernels that lack checkpoints, process-level
crash injection, high-DPI Windows stress tests, multi-gigabyte I/O
benchmarks and real UGC 2162 FITS validation remain future gates.

**P0.4 is incrementally hardened; it is not yet an independently certified
full GUI refactor.** Additional controllers for target, survey and science
workflows will follow after the worker lifecycle is stable on CI.


## P0.4.2 — controller boundary refactor

The former `SurveyExplorer` now delegates three coherent responsibilities
while preserving its established public methods and Qt signal wiring:

- `gui.target_controller.TargetController`: strict ICRS `SkyField` reading,
  invalid-entry feedback, coordinate-pair clipboard insertion and conservative
  North/South hosted-preview region. `TargetTextField` moved alongside it,
  while remaining re-exported from `gui.survey_explorer`.
- `gui.survey_data_controller.SurveyDataController`: source-product table,
  file filters, local-versus-remote evidence, complete quartet selection and
  checked-path selection. Visible-only filtering must NEVER mutate hidden
  selected products. Choosing a quartet explicitly selects exactly
  IMAGE, INVVAR, NEXP and MASKBITS, clearing unrelated hidden selections.
- `gui.science_context_controller.ScienceContextController`: frozen
  `ScienceContext` capture with project UUID and checks for stale selection,
  parameters and FITS raster identity. No file writes or science computations.

`SurveyExplorer` keeps compatibility forwarding methods: `_field`,
`_validated_field`, `_target_edited`, `paste_coordinate_pair`,
`_preferred_viewer_region`, `_fill_files`, `_filter_products`,
`_selected`, `select_science_quartet`, `_current_science_context`,
`_science_ticket_current` and `_qa_context_current`.

This is functional separation, not numerical science or visual redesign.
Tests inspecting old source method bodies now inspect the correct
controller implementation. Real Qt tests also verify the behavioral
invariants and reopen several isolated galaxy workspaces.

Multiple successful open/close cycles are NOT an RSS/VRAM/QtWebEngine
leak acceptance test; real Windows high-DPI, original large DR9 FITS
and forced process termination remain outstanding. Additional GUI
workflow extraction and deterministic science API contracts are next.
