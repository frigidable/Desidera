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
