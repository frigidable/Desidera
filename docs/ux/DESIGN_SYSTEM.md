# LSBpipe Scientific Workbench — design system v0.1

Status: historical v0.1 design proposal, **not an audit of the current UI**.
Apply its visual principles to new work only after inspecting current Qt
widgets and screenshots; do not treat proposed components or color tokens
as implemented or visually approved.
Implemented Qt widgets must be measured against actual screenshot evidence
before this becomes a visual-regression baseline.

## Windows GUI failure analysis — QML Text property (15 November 2026)

Exact-head Windows baseline `43745395770` proved 988 science tests passed,
but timed out within the GUI suite after reaching the Project Hub tests.
The `if: always()` screenshot stage supplied an actionable QML load error:
`project_hub_preview.qml:510: Cannot assign to non-existent property
"selectByMouse"`. Three expanded-science-source `Text` items
incorrectly used a `TextField`/`TextEdit`-only property. Fixed the
three declarations before continuing the Observatory Studio Workspace
acceptance; retain `selectByMouse` on the actual search `TextField`.
Add an early isolated QML hub preflight and verbose GUI test names with
pytest faulthandler timeouts so a QML component construction failure
cannot silently absorb the entire 20-minute Windows screenshot budget.
This was a **GUI** defect, not proof of a PSF or WCS failure.

---

## Galaxy Workspace navigation and science layout — desktop pilot

The real Leo P Windows session showed THREE overlapping navigational layers:
sidebar Projects/Sources/Photometry, a second outer `Remote sky / FITS
Explorer / Photometry` tab bar, and a Photometry-specific Workspace/Sources
toolbar. This makes it hard to understand which workspace owns the active
scientific state.

The first Observatory Studio Galaxy Workspace increment keeps the established
`SurveyExplorer.preview_tabs` as a **nonvisual page host**, preserving its
object-bound project identity, FITS viewer and existing signals/resume. Its
visible tab bar is always hidden; a single explicit **Sky / FITS /
Photometry** switcher lives under Projects/Sources in the persistent sidebar.
Native Data Manager remains a guarded Sources action; FITS comparison
workspaces remain accessible through their original buttons. Selected-state
styling is driven by the actual current page, not a guessed readiness state.
In embedded mode only, `PhotometryWorkspace` hides redundant Workspace and
Sources buttons, recovers header height and retains their signals for
standalone/compatibility behavior. The four PSF/Scene/Sky/Apertures tabs
are a DIFFERENT (scientific) level and deliberately remain visible.

At desktop widths, `PSF0Panel` now places **model + real radial comparison**
on the left and **quality + existing controls** in one inspector on the right.
On narrower viewports it reflows vertically without computing a PSF,
changing a source band, changing starred exclusions or promoting any result
to scientific acceptance. The no-model/no-profile placeholders are still
explicitly empty. Scene Editor preserves the original QGraphicsView/FITS
canvas and ICRS draft semantics, but has a named contextual inspector
with bounded width and unchanged Draw/Isophote/SGA task routing.
Actual first embedded 1366×768 offscreen captures showed the Scene input
toolbar needlessly split into two rows: band and cutout-width selectors
expanded to hundreds of pixels. Bound the band selector to 88px and native
view-width input to 176px, align them to the left, and use one toolbar row
at viewports ≥920px while preserving the original 760px/high-font wrap.
This recovers image-canvas height without changing stretch, pixels, or
coordinates. A PSF-0 GUI check initially expected an empty band to show
model/profile cards; that contradicts the existing fail-closed policy.
The test now verifies those cards are absent when no original FITS/PSF
exists and that the QA/setup inspector remains accessible.

Test this as layout-only: no new project folders, no backend imports from
QML, no remote WebEngine loading during offline capture, and no implicit
read/write of native FITS. CI captures the real embedded Galaxy Workspace
at 1366×768 and 1600×900, in both PSF and Scene states. The user's
physical Windows test remains the source of truth for native DPI,
focus and real science-image presentation.

---

## Real Windows Project Hub usability review (15 November 2026)

Real on-device screenshots from a five-galaxy library (Leo P, dw
1343+58, KDG 61, KKR 25, UGC 2162) exposed defects not reproduced by
the earlier three-object CI fixture:

- With the source evidence expanded, the disclosure button lived inside
  the scrollable evidence column. Qt Quick scrolled to the focused button,
  clipping the beginning of the passport and forcing the user to read a
  continuous dense column with weak separation between provider receipts.
- At a 1680×900 desktop viewport, the fifth galaxy's list card was partly
  hidden by the pinned New/Import toolbar.
- The 76px cached RGB preview was too small relative to a wide empty
  inspector, obscuring its value as a visual identification aid.

Fix: keep the evidence disclosure control **outside the scroll area**,
as a permanently visible command above the recommended next action.
Switch between a concise status summary and a **read-only
source-documentation mode** with separated catalogue, FITS/WCS QA and
project-provenance text sections. Avoid height-coupled nested cards in a
Qt Quick ScrollView: the first attempt made the offscreen GUI smoke step
unusually slow, so the source layout was simplified before acceptance.
Reset evidence scroll to the top on each mode switch; clicking sources
must never modify the selected galaxy UUID.
Use read-only, plain-text, line-wrapped source data with 13px reading
type rather than the tiny dense receipt block. Standard QML `Text` does
not support `selectByMouse` on this PySide6 runner: unlike `TextField`,
it must not declare that property. Direct text selection/copy is a later
design requirement, to implement with a supported read-only `TextEdit`
or a dedicated Copy provenance control, tested on Windows. Enlarge the cached
display-only thumbnail when available, but keep the ordinary placeholder
compact when no image is cached. Reduce the 5-row catalogue card height
to 84px and spacing to 6px; short windows still scroll normally.

These changes affect presentation only. They do not reinterpret the saved
FITS/WCS diagnostic as new content-hash verification or imply adopted
surface photometry. Add matched light/dark offscreen captures of *both*
summary and expanded documentation, and GUI tests exercising disclosure,
filtering, the selected project UUID and on-disk passport immutability.

---

## Qt Quick scientific-evidence and actions pilot

The optional Qt Quick inspector now exposes the same evidence model as
`GalaxyDetailCard`: confirmed SIMBAD/NED associations (never unverified
cone candidates), published NED descriptors marked unadopted, existing
`qa_status` and receipt-derived `audit.explanation`, the actual project
root, aliases, and the readiness-derived next action. The QML expanded
details are independent of whether the *Widgets* disclosure is expanded;
both read the same immutable `GalaxyCardSnapshot`. Text from catalogues
and paths must be rendered with `Text.PlainText`, not interpreted as QML
rich text. A thumbnail, when available, is a previously cached display
preview; the original science FITS image is never passed to QML.

**Existing user actions are reused exactly** through the native Project
Hub: SIMBAD identity dialog, NED evidence dialog, analysis configurations,
local-FITS thumbnail generation, explicit remote-RGB request, and guarded
removal from the recent-project list. A whitelist in the bridge checks each
existing `QAction.isEnabled()`; QML cannot enable an unavailable action,
bypass project UUID guards or launch networking from merely viewing a card.
The remote-RGB item is explicitly labelled as a network request. During a
native modal action the experimental window yields focus to Widgets, then
resynchronizes the same model. Creating or importing a project retains
the established native dialogs and the one-galaxy project model.

Missing or corrupt receipts fail closed, never as a successful association
or science validation. A historical FITS/WCS audit describes unchanged
file *metadata*, not a fresh content-hash verification, PSF readiness or
scientific acceptance. Tests exercise unverified identities, disabled
actions, original UUID persistence, audit receipt provenance and detection
of changed FITS file metadata. QML remains **experimental** until native
Windows GPU, accessibility, input/navigation and functional parity are
verified. The currently working Widgets interface remains authoritative.

The first Windows GUI rerun found a test-harness error: children of a
QML `ScrollView` can have a visual `QQuickItem.childItems()` parent chain
different from their `QObject.findChild()` chain. The regression now checks
the mounted ScrollView, the actual evidence projections and the rendered
screenshot; its internal QQuickItem hierarchy is an implementation detail,
not a stable test contract. Visual inspection also caught a black scrollbar
in both themes. The pilot initializes the Qt Quick `Basic` control style only
once per process (before its first QML view) and draws a theme-colored
ScrollBar explicitly, without touching native Qt Widgets/FITS styles.

---

## Qt Quick Project Hub A/B experiment (15 November 2026)

The first **measured implementation comparison** is opt-in: use
**Display settings → Compare Qt Quick Project Hub (experimental)…** to
open a separate `QQuickView` window. Normal startup and scientific work
continue to use the existing Qt Widgets Project Hub and QGraphicsView FITS
renderer. This experiment is not the application's default GUI.

Implementation: `project_hub_qml.py` exposes a `QAbstractListModel` which
projects the exact **visible items and UUIDs** of the already-existing
Widgets library. It does not create a second registry or calculate data
readiness. `project_hub_preview.qml` renders the same name, status,
coordinates, selected galaxy and readiness-derived next action as the Qt
Widgets card. Search, selection and open/new/import commands call the
established Project Hub functions (including strict relink and project
identity checks). Missing data are labelled as missing; SCIENCE HOLD is
always explicit. There are no fabricated thumbnails, magnitudes or results.

This is a **standalone QQuickView**, not a QQuickWidget embedded beside
the FITS canvas. The pilot therefore avoids changing render ownership of
scientific widgets and allows a meaningful comparison before considering
a mixed-window shell. The `.qml` resource is included in the Python wheel.
The preview closes when an object workspace opens; it has no independent
session state and performs no background work.

Visual review of the initial actual QML renders exposed two concrete
problems: the placeholder orbit glyph was unavailable in the runner's font,
and resetting the QML list to apply a theme accidentally selected a different
galaxy. Both are release-blocking. The QML pilot now draws vector orbit marks
from Qt Quick shapes, and its list mirrors the existing Hub selection without
writing it during model resets. Arrow navigation is an explicit user action
which alone changes the selected UUID. A regression test and the screenshot
script assert the exact project ID survives light/dark transitions.

Qt CI uses the software scene-graph adaptation in offscreen Windows, and
captures matched **1366×768, three-object** Widgets and Qt Quick light/dark
examples. These are comparative screenshots, not proof of real Windows
125–200% DPI, Fluent material, GPU performance, WebEngine interoperability
or keyboard/screen-reader acceptance. Treat the prototype as a candidate,
not a decision.

Before making Qt Quick the default, evaluate with real on-device evidence:
- Visual hierarchy and readable metadata at 1366×768, 1600×900, 1920×1080,
  with large fonts; list/inspector search, selection and missing project.
- Actual Qt keyboard access, screen reader semantics, focus, dialog/return
  navigation, and error states; no regression of native Widgets actions.
- Repeated cold and warm launch timings, peak and steady-state RAM, and
  rendering latency on supported Windows hardware and without the remote
  viewer. The synthetic CI timings must not be treated as a benchmark.
- Safe handling of 100+ galaxies, exact UUID preservation and corrupt
  index recovery. Never infer scientific adoption from file presence.
- Only then consider the architectural cost of embedding QML and any extra
  Qt Quick dependencies. Maintain one testable Python project/readiness
  domain model and no second astronomical pipeline.

---

## Project Hub catalogue pilot — visual architecture

The first Observatory Studio migration increments the **existing** offline
Project Hub rather than adding a second launcher or adopting QML before
a benchmark. Each row is painted by `GalaxyProjectRowDelegate`, a lightweight
`QStyledItemDelegate` over the same `QListWidget` model. The row's native
`DisplayRole` text and immutable `UserRole` project record are preserved
for keyboard access, assistive technology, search, and object-safe relinking.
The additional `UserRole + 2` carries only the `(name, status, coordinates)`
presentation tuple. There is no FITS decoding, network request or science
acceptance in a delegate `paint()` operation.

A list row displays three unambiguous lines: **galaxy name**, **workflow
state**, and **ICRS coordinates or a missing-data explanation**. Missing
projects use a warning-color status with ordinary explanatory text, while
available projects use the interface accent; neither means validated science.
A two-pixel focus border remains visible independently of selection color.
The library header shows the number of saved galaxies and filtered matches.
Empty and no-match states are distinguished and retain keyboard controls.
The selected galaxy passport also reveals the *recommended next step* from
`snapshot.readiness.next_action_label`; it cannot label a dataset or science
measurement as validated on file presence alone. An empty selection hides
the readiness callout.

Three-project synthetic screenshots at 1366×768 and 1600×900 complement
the established empty and single-object baselines. Keep this phase entirely
in Qt Widgets. Evaluate a QML alternative only if a measured prototype
improves layout, accessibility and maintenance without regressing real FITS
interaction. No new design-specific project registry or directory is needed.

---
## Observatory implementation (15 November 2026)

The earlier v0.1 color table below is a **historical proposal**, not the
active implementation. The current Qt Widgets palette is defined in
`src/lsbpipe/gui/visual_theme.py`. Project Hub now exposes:
**Display settings → Interface appearance → Follow system / Light / Dark**.
The exclusive choice is stored in Qt `QSettings` as `appearance/theme`,
outside the galaxy project and all scientific FITS/QA products.
Invalid stored values revert to Follow system.

| Active token | Light | Dark |
| --- | --- | --- |
| Background | #F2F5FA | #0C1424 |
| Surface | #FFFFFF | #172238 |
| Raised panel | #E9F0FA | #1D2C46 |
| Primary text | #1A2B45 | #E9F1FD |
| Secondary text | #4D617A | #AABAD1 |
| Primary action | #215CAA | #65B8FA |
| Focus outline | #117F8C | #74D8D1 |
| Warning text | #865313 | #EEC47A |
| Error text | #A33848 | #F3A0A5 |
| Border | #CED9E8 | #344862 |

Changing appearance refreshes open Qt workbench surfaces, including
the selected galaxy card, semantic actions, vector icons and Scene task
selection. It never fetches data, opens FITS, runs a worker or changes
scientific adoption status. The four colored Photometry stage icons
identify *methods*, not evidence quality. FITS pixels, display colormaps
and remote QtWebEngine content remain outside Observatory QSS.

GUI CI now exercises the explicit theme selection and captures separate
light/dark Project Hub screenshots at 1366×768. These are controlled
offscreen Qt appearances, **not** tests of actual dark OS settings,
Windows high contrast, NVDA, monitor DPI or authentic survey observations.
Those E3/E4 acceptance gates remain open. Scientific HOLD is unchanged.

---
## Visual intent
A precise scientific workstation rather than a marketing dashboard.
The FITS canvas and evidence occupy primary visual attention. Quiet,
consistent panels carry controls and metadata. Avoid gradients as a
substitute for hierarchy, unnecessary cards, and oversaturated badges.

## Color tokens (hex)
| Role | Dark | Light |
| --- | --- | --- |
| Window | #111923 | #F6F8FC |
| Panel | #182331 | #FFFFFF |
| Elevated | #223143 | #E9EEF5 |
| Primary text | #E8F1FA | #172232 |
| Secondary text | #AAB9C8 | #415569 |
| Interactive accent | #59B7FF | #1C5D99 |
| Warning text | #F8C678 | #8C4F00 |
| Error text | #FF9E9E | #9B2F3A |
| Success text | #81D4A3 | #176641 |
| Border | #334558 | #CAD4DF |

The primary/secondary/accent pairs above have been numerically checked
against their primary panel/window backgrounds for >= 4.5:1 text contrast.
This does not certify every hover, disabled, overlay or state combination:
verify actual rendered foreground/background pairs before shipping.
Never rely on hue alone for scientific data categories or warnings.

## Typography and geometry
- Font: system UI (Segoe UI on Windows), fallback Noto Sans / sans-serif.
  Preserve astronomical symbols and Greek characters correctly.
- Base body: 13–14 px at 100%; dense table: 12–13 px; title: 18–20 px;
  smaller metadata: >= 12 px where practicable.
- Spacing: 4, 8, 12, 16, 24, 32 px; radii 4/8 px; consistent borders.
- Controls: nominal >= 32 px high (36 px for primary actions);
  do not shrink interactive hit regions to fit a crowded sidebar.
- Use a single primary action per state; secondary actions subdued but
  discoverable. Labels describe actions, not internal method names.
- Support resizable windows, high DPI, font scaling and keyboard focus.
  Never encode pixel fixed heights as the only way to display controls.
- Splitter state must be adjustable and recoverable. Scrollable inspectors
  are preferable to hidden controls on a short screen.

## Information density and progressive disclosure

- Default workflow surfaces should display the **current object, relevant
  science/data state, one primary action and the main visual**. Metadata, log
  receipts and alternative methods remain accessible from their contextual
  screen instead of appearing as permanent instructions.
- A compact form uses short *visible* labels, with explicit units, complete
  accessible names and contextual help in tooltips. Scientific status such as
  `DIAGNOSTIC_NOT_ADOPTED` may be shortened visually only when the non-adopted
  meaning remains legible on the active screen.
- Avoid cascading disclosure controls within another disclosure. Prefer one
  named optional-settings section per task; once expanded, show its actual
  inputs directly instead of requiring a second reveal.
- Never fill an empty workspace with large panels of disabled batch buttons
  or repeated placeholders. A concise no-data outcome and an actionable
  source/workspace navigation should take precedence. Researchers can switch
  to Advanced mode to review the complete scientific controls.
- Distinguish **separate workflow stages** (e.g. Scene/Background/Photometry)
  from redundant nested tabs. Stage navigation remains visible and stable;
  additional filters/bands are context selectors, not new primary screens.
- Preserve full error messages, provenance, FITS raw values, masking
  semantics, units and QA reports through explicit accessible detail actions.
  Never drop a scientific caveat merely to meet a visual density target.
- At 1366×768 inspect both default and advanced screen states. Count visible
  disclosures and duplicated paragraphs, check that keyboard focus finds
  the advanced evidence and that the image/plot still dominates.

## Information hierarchy
Default layout (only when viewport permits):
1. Compact top row: project/target, release/band, data source, task state.
2. Left navigation: Survey, Data, Scene, Background, Photometry, Results
   according to actual capability.
3. Center: largest continuous surface for hosted preview/local FITS
   visualization; clear preview/FITS status and coordinate scale.
4. Right contextual inspector for current task and relevant options.
5. Collapsible bottom run/status evidence log (do not force modal dialogs
   for routine non-blocking messages).

At ~1366x768, keep task hierarchy usable via scroll/stacking; do not
assume all panes can remain full width. Avoid nested permanent toolbars.

## Reusable Qt components (target contracts)
- `TargetHeader`: ICRS, selected release, selected band, source provenance.
- `WorkspaceNavigation`: capabilities, not speculative runnable modes.
- `ViewerSurface`: hosted preview vs local FITS explicitly separated.
- `ContextInspector`: only properties relevant to the current operation.
- `EvidenceStatus`: status + source + uncertainty + recommended action.
- `ScienceFileTable`: product identity, release/region/brick/band,
  evidence, known or unknown size, selection persistence.
- `TaskProgress`: quantitative progress only when measured; safe retries.
- `QualityReport`: diagnostic vs adopted status and source evidence.
Do not create components solely because this list exists. Extract them
when a real screen has a working scenario and accompanying tests.

## Scientific visualization contracts
- A celestial image must expose FITS source, band, coordinate convention,
  scale and displayed stretch without changing underlying pixel data.
- Image, IVAR, NEXP, MASKBITS, derived model and residual are distinct
  layers with separate provenance. Displayed previews are not analysis data.
- Overlay geometry must preserve WCS alignment; changing image source
  invalidates stale overlays unless reprojected/verified.
- Distinguish valid, missing and masked pixels; negative measured flux
  remains negative. Never clip negative data as a scientific transform.
- Plots label axes with physically correct units, uncertainty conventions,
  and clear data/model/residual encoding; default to accessibility-friendly
  patterns, line styles and symbols.

## UX language
- "Preview (remote)" / "FITS (local)" instead of ambiguous "Image".
- "Candidate file" / "HTTP checked" / "Downloaded" are different stages.
- "Coverage diagnostic — not adopted" for local coverage screening.
- Prefer explicit, repairable error messages to catch-all failure dialogs.

## Influence and licensing
Use official Qt accessibility principles, napari scientific-layer
interactions, and established research-software UX patterns as conceptual
references; do not copy third-party proprietary skins or assets.
Qt/PySide6 remains the first implementation target.
