# UI/UX quality gates and evidence requirements

These rules are **stage-aware**. Initial legacy captures are baselines
for investigation, not an endorsement of the legacy UI.

## Mandatory gates
**G0 Science isolation — blocking:** no accidental modifications to core
photometry, WCS, FITS validity, mask policies, transfer semantics or
adoption status. Existing regression tests must run; any failure described
and resolved before release.

**G1 Deterministic GUI execution — blocking once GUI is modified:**
construct and dispose the GUI on a Windows runner without network,
capture a valid image, and publish machine-readable layout metadata.
A screenshot without successful capture checks is not passing evidence.

**G2 Scenario correctness — blocking for altered workflows:** verify
real user-level state transitions and negative paths, including
file-selection semantics, missing data and error recovery.
No fake measurements or previews labeled as source FITS.

**G3 Layout usability — blocking before marking UI-ready:** inspect
screenshots at required viewports and check essential controls remain
reachable. No clipped primary controls, occluded dialogs, inaccessible
actions, unreadable labels or overflow that prevents task completion.

**G4 Accessibility — blocking before marking UI-ready:** keyboard
focus order, accessible labels, non-color-only statuses, useful
contrast (WCAG AA 4.5:1 for ordinary text where applicable), DPI/font
scaling, and notification accessibility.

**G5 Science communication — blocking:** data release, band,
remote/local source, completeness/evidence and diagnostic/adoption
qualification must be unambiguous. Do not score this away.

**G6 Visual consistency — quality target:** token adherence, coherent
density, typography, alignment and iconography. Record actual
screenshots and independently review them before acceptance.

## Cognitive-load acceptance criteria

A visual review of a GUI change must record (a) the primary next action,
(b) count and depth of nested disclosures on the active workflow, (c) whether
disabled/no-data panels dominate the viewport, (d) duplicate explanation
lines and (e) whether warnings, masks and original FITS provenance remain
reachable **without color identification**. A new accordion is not a default
fix for a crowded screen: first consider removing duplication, shortening
labels while preserving accessible details, or regrouping controls.

The researcher must be able to operate the default task without expanding
multiple nested panels. Advanced/research-specific controls may be
separately available but must retain their exact scientific meaning and
functional behavior. Evaluate this in the 1366×768 and high-font fixture,
in addition to the established screenshot and regression requirements.
Static source review alone never certifies this gate.

## Evidence levels
- E0: source reviewed. No visual claims.
- E1: local/CI tests executed; no screenshot examined.
- E2: screenshot artifact successfully captured and inspected.
- E3: scenario interaction exercised on real GUI; visual and functional
  evidence inspected.
- E4: separate real-Windows/on-device test including hosted Viewer and
  actual FITS path, documented limitations and researcher acceptance.
G0/G1 passing does not imply E2/E3/E4.

## Scoring (human-assisted only)
Candidate composite weights:
science clarity 25, workflow 20, visualization 20, interaction 15,
visual consistency 10, accessibility 10 (total 100).
Suggested target >=85 and **zero blocking G0–G5 issues**, but the
composite is not a validated UX metric. Never auto-assign an 85/100
score from code appearance or screenshot existence. Per-category
score requires written evidence, reviewer and linked screenshots.

## Automated initial P0 gate
The baseline CI is intentionally narrower: it runs a deterministic
offline application capture and checks PNG/JSON presence/dimensions.
It does NOT claim visual quality, online map tile rendering,
science validation or complete GUI acceptance. Visual hard gates
will be introduced only after observed baselines and fixtures exist.
