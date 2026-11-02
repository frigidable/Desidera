# Galaxy Projects — compact right-hand object passport (2026-11-03)

## Decision

Keep the current DESIDERA Galaxy Projects screen's structure:

- **Left:** the existing searchable list of recent galaxies, not a grid.
- **Right:** exactly one compact galaxy passport for the selected list item.
- **Default splitter:** enlarge the previous approximately 40%-wide right pane
  by roughly 30% of its *own* width: target about 52% right / 48% left.
  The divider stays draggable; no fixed-width layout or nested gallery.
- At small window heights or scaled fonts, only the passport's detail content
  scrolls; the object name and primary Open/Relink action remain accessible.

The preview intentionally retains the project's original context rather than
displaying invented distances, velocities, photometry or QA metrics.

## Safety and parallel development

Branch: `feature/galaxy-cards-isolated-20261103`.
Base: `fix/webengine-render-diagnostics-p056` at
`7d42b0fa0ea4d434dc5acebb36077c011438d8c4` (#19).

This is an **isolated draft**. No edits to existing
`src/lsbpipe/gui/project_hub.py`, project registry, saved manifests, original
FITS, science runs or Dataset Manager. In particular, none of the changes
touches work on the CCD observation-lineage branch (#20). Integration into
the production launcher is separately reviewed after the ongoing GUI work.

## Implemented widgets

- `GalaxyDetailCard` — reusable PySide6 right-side inspector. Displays:
  - Object name and ICRS coordinates;
  - Optional *locally cached* thumbnail if it exists; otherwise an honest
    empty image state;
  - Saved dataset (release, region, brick, band);
  - Folder availability and FITS link state, **without pretending QA passed**;
  - Aliases, last-opened value and project path;
  - Missing-project/relink label and a primary Open/Relink signal.
- `GalaxyCardsPreview` — opt-in standalone two-pane harness with list search,
  selection synchronization, 48/52 resizable splitter and the single passport.
  It is not the production ProjectHub.
- `galaxy_card_model.py` — read-only adapter of the existing registry to
  stable card snapshots. It never reads full FITS or writes manifests.

No network activity, catalogue resolution, downloads, photometry or automatic
thumbnail generation are initiated by opening this view. The preview reads only
the optional dedicated `<root>/.desidera/card-preview.png` that already exists.

## Standalone inspection

With DESIDERA installed with GUI dependencies:

```powershell
python -m lsbpipe.gui.galaxy_cards
```

This displays the user's saved galaxy list and a right-side passport. Clicking
Open/Relink in the standalone preview explains that the real launcher must be
used. In the eventual production integration, the same `openRequested`
project-ID signal must invoke the existing `ProjectHub` open/relink code.

## Validation and integration gates

- `python -m pytest tests/test_galaxy_card_model.py`
- `QT_QPA_PLATFORM=offscreen python -m pytest tests/gui/test_galaxy_cards_gui.py`
- Existing ProjectHub tests, regression of no-network launch and no-unintended
  science/QA adoption.
- Visual acceptance screenshots at 1100x720, 1366x768, and Windows 125–200%
  DPI; readable at 760px width and with long project paths.
- Confirm the card fits *within* the existing right pane; no full-width new page
  and no secondary grid of cards.

After the active Project Hub branch settles, integration changes may replace
its right-panel label stack with `GalaxyDetailCard` while keeping the left list,
existing names/roles and application startup. Preserve selection identity using
the immutable project ID. Existing Remove from Recent stays the existing
ProjectHub action; it is not replicated as a file-deleting action in the card.

Do not merge this draft until CI and Windows GUI acceptance are verified.
