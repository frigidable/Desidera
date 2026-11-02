# Galaxy Cards P0 — isolated implementation slice (2026-11-03)

## Scope and branch safety

Branch: `feature/galaxy-cards-isolated-20261103`.
Base: `fix/webengine-render-diagnostics-p056` at
`7d42b0fa0ea4d434dc5acebb36077c011438d8c4` (Draft PR #19).

The current Project Hub, its main-list widget, the project registry, science
modules, manifests, FITS data and existing tests are **unchanged**.
This initial slice adds an opt-in component and two new test modules only.
It is intentionally not merged and does not depend on the observation-lineage
work in PR #20. Request integration after rebasing on the reviewed GUI tip.

## Present behavior

- Read existing `recent_projects()` rows and saved object metadata without
  touching `last_opened_at`, files or scientific reports.
- Display one card per registry entry: identity, ICRS coordinates, selected
  dataset, local-path state and a conservative FITS-availability description.
- Search name, aliases, project path, identifier, coordinates and dataset.
- Responsive grid, selected state, keyboard Enter to request opening and
  Space to select; double click emits an open request for host integration.
- Missing roots, old local projects and invalid manifests get explicit,
  non-fatal labels.
- Optional preview is loaded only if a separate
  `<project-root>/.desidera/card-preview.png` already exists. There is no
  internet access, FITS decoding, background work or thumbnail generation.
  Empty previews are deliberately labelled rather than substituted by
  unrelated catalogue imagery.

Smoke-preview command from installed project with GUI dependencies:

```powershell
python -m lsbpipe.gui.galaxy_cards
```

This preview **does not replace** the production `desidera-gui` entry point.

## Next integration step (separate, explicitly reviewed change)

Once the active Project Hub branch is accepted, add an optional
`Cards / List` switch to `ProjectHub`, keeping the existing list as an
available view for dense selections and accessibility. Connect the gallery's
`openRequested(project_id)` and `selectionChanged(project_id)` to the
existing project lookup by registry identity; never use card display names as
file or object identifiers. Preserve selected identity across views and
filters, and reuse the existing relink/open/remove behavior without
modifying the project manifest merely by selecting a card.

Separate follow-up: a bounded, explicitly invoked thumbnail generator that
uses valid local FITS with verified band/WCS association. Its RGB renders must
be visually labelled as display products, never evidence of photometric
measurement or FITS QA. Cache writes should be atomic and located outside
immutable observations; no automatic image download on the Hub screen.

## Validation

- `python -m pytest tests/test_galaxy_card_model.py`
- `python -m pytest tests/gui/test_galaxy_cards_gui.py` with PySide6, pytest-qt,
  and `QT_QPA_PLATFORM=offscreen` in headless CI.
- Real Windows high-DPI and keyboard checks remain required before merging.
- Existing Project Hub GUI regression tests must still pass after integration.

Current statuses intentionally distinguish `FITS linked · integrity not
checked` from verified data, and catalogue/project identification from
scientific result adoption.
