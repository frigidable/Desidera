# LSBpipe — progressive UX contract (P0.4.3)

This document defines the ordinary researcher journey and mandatory
safety invariants. It replaces a module-centric surface with a
**three-stage workflow**, without removing advanced scientific tools.

## Core researcher journey

| Stage | Researcher sees | Primary action | Details not shown until needed |
| --- | --- | --- | --- |
| 01 Target | ICRS RA, Dec, field width, Find data, sky preview | Find data | Multiple DRs, clipboard, external viewer, project diagnostics |
| 02 Data | Intersecting bricks, selected file count, download | Explicitly choose one science quartet, then Download | HTTP HEAD, raw file/type filters, manual bulk selection, aperture/coverage settings, manifest export |
| 03 Inspect | FITS Explorer, quick Open FITS, read-only audit | Open FITS Explorer | Extra mask statistics, regional comparisons, PNG export, full reports, technical provenance |

Projects with an existing local FITS inventory enter Stage 03 directly.
A new object without data enters Stage 01; opening its workspace must
NOT automatically show an empty product list on the right. Discovery of
bricks enters Stage 02. Nothing is adopted as science merely because
the UI advances.

The stage caption is a *non-interactive status cue*, not an additional
series of buttons. It always tells the user the next meaningful action.
When the selected brick is missing or the target is edited after a
search, the workflow falls back to the appropriate earlier stage.
It must not auto-select an alternate brick or band.

## Progressive disclosure

- Ordinary workflow emphasizes a single research task at a time.
  Scientific diagnostics are optional, not first-class navigation.
- Target-specific secondary controls are under **More target options**.
- Archive HTTP verification, bulk selection, survey filters, coverage
  ring radii and manifest save are under **Advanced file options**.
- The science panel has two prominent actions: integrity audit of
  existing originals and **Open FITS Explorer** for a selected brick.
  More algorithms are under **Advanced diagnostics**.
- Verbose output, immutable run receipts and history are under
  **Report details & history**. A short outcome label remains visible.
- The FITS graphics viewport dominates its screen; contrast and layer
  controls are in a narrow sidebar, masking and provenance in their
  own disclosures. Four-layer QA is shown only on request.

Technical provenance and QA fields **must remain readable**. Visual
disclosure cannot be used to drop physical units, suppress diagnostic
status, invent missing coverage or change data-masking policy.

## Safe-by-default science and transfer

The archive manifest initially selects **zero** files. A survey
search may offer many bricks/products, but cannot implicitly schedule
large downloads. The researcher must select a brick and science
quartet, or explicitly choose individual entries. The Download button
is disabled until a nonempty selection exists. Existing local project
inventories retain their recorded on-disk selection.

The UI must preserve separate North/South families, original FITS
integrity, frozen science context and cancellation state.
Background progress/Stop controls appear only while a task exists;
the underlying task journal and safe-checkpoint cancellation remain
unchanged.

## Visual language and accessibility

Use short, verb-led action labels; one heading hierarchy; consistent
spacing; OS palette instead of forced light/dark colors. Scientific
metadata belongs in accessible tooltips and disclosure panes, not
multi-line fixed headers over images. A hidden control must never be
needed merely to understand an invalid target or a failed operation:
alerts and the workflow cue remain visible.

The Project Hub remains the one landing screen. Opening an object
must not trigger downloads or scientifically accept an old QA result.
File reopening is offline. A separate explicit user action can
launch the hosted Legacy Viewer.

## Acceptance criteria

Automated tests should demonstrate:
- clean Stage 01 with no empty brick grid, no idle progress bar and no
  unrequested archive panel;
- Stage 02 only after an actual catalog result; zero default-selected
  archive downloads; selected quartet gives exactly four files;
- Stage 03 for actual local FITS inventory, and a visible quick Open
  FITS action when a brick/band is selected;
- extra target/file/QA/report controls only when opened;
- data edits invalidate stale FITS/QA without silently moving the
  project's ICRS position or switching dataset;
- Science QA and GUI rendering keep their original units and
  diagnostics; original FITS never change.

## Follow-up design gates

This implementation is a first major structural simplification, not
the final visual-design acceptance. A later stage needs:
1. Real Windows screenshot walkthrough at 1366×768, 1920×1080 and
   2560×1440, including 125–150% system DPI.
2. Keyboard focus, Tab order, screen-reader names and error-state QA.
3. Real-world UGC 2162 offline FITS flow and new-object discovery flow.
4. A dedicated file-download wizard or compact single-brick picker if
   researchers still find the product table overwhelming.
5. Tiled FITS rendering and memory profiling for large cutouts.
6. A photometry execution screen only after scientifically validated
   P0.5+ algorithms exist; no dead-end fake photometry buttons.

Avoid equating passing synthetic CI with positive real-user UX evidence.

## Post-capture refinement: keep the actual FITS image dominant

A four-state Windows offline screenshot review (initial, archive found,
science quartet selected, FITS inspection; 1366x768) exposed a spatial
problem even after the first simplification: left navigation + right
inspector + the FITS sidebar left too little width for the actual
galaxy image.

Two refinements follow from that evidence:

- The ordinary Stage 02 archive inspector has **no default 88-row
  spreadsheet**. It shows a clear four-FITS preset and Download action.
  A researcher opens **Browse / customize files** to inspect the complete
  product table, HTTP evidence, file filters and manual selection. The
  default input selection is empty, never a mass download.
- Entering **FITS Explorer** or interactive **Four-layer QA** hides the
  right-hand archive/science inspector. It can be reopened explicitly,
  but should not automatically squeeze the image every time a local
  dataset is loaded. Ordinary FITS viewing therefore has the left
  project/navigation context and a wide central inspection area.

The right inspector retains a readable minimum width when explicitly
opened. Captured comparison figures are offscreen Qt layout evidence
only, not proof of ergonomic behavior on real Windows displays.
