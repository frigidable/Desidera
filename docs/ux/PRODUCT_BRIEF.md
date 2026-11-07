# DESIDERA — object-centred scientific workbench

**Living product requirements (2026-11-07), not a science-validation certificate.**
This replaces the 2026-10-30 Survey Explorer-only planning snapshot.
Actual behavior is determined by the code and tests for the relevant commit.

## Researcher and primary workflow

The primary researcher studies galaxies, including low-surface-brightness systems, using original DESI Legacy Surveys FITS and independently sourced astronomical context. DESIDERA should keep an individual object's project, data, QA and analysis history together while preventing accidental scientific adoption.

1. **Projects.** Start in the Galaxy Project Hub. Create a project with an optional ICRS position, select an existing galaxy, import an older local project or relink a moved folder. Keep project UUIDs and previously downloaded science FITS.
2. **Context.** For the selected galaxy, inspect saved source-labelled facts, local observation readiness and confirmed SIMBAD/NED identity receipts. Remote catalogue searches require user actions and candidate confirmation; nothing automatically changes adopted coordinates or science results.
3. **Data.** Discover survey regions/bricks and explicitly select FITS products. Distinguish catalogue coverage, URL/HEAD evidence and verified local data. Reopening must not redownload original FITS.
4. **Inspection.** Browse original FITS, WCS positions, the source quartet (IMAGE, INVVAR, NEXP, MASKBITS), pixel validity and diagnostic quality maps. A Legacy Viewer or RGB preview is display-only and must never be presented as original calibrated input.
5. **Analysis and QA.** Offer implemented measurement/diagnostic workflows behind their actual prerequisites, preserving signed flux, explicit sky assumptions, source hashes, science-context identity and run provenance. An experimental result is not adopted by default.

## Current engineering scope

The Project Hub, survey discovery/download workflow, FITS views, dataset services, galaxy information receipts and diagnostic QA have Python implementations and regression tests in `src/lsbpipe/`. `src/desidera/` provides the current public identity and launchers. The fixed-ellipse numerical core is implemented; scene-mask and candidate-background/sensitivity code exists as experimental research tooling.

Do not confuse **code present** with **science method validated**. High-confidence PSF/galaxy separation, calibrated sky estimation, cross-band total-flux inference, accepted SBF distances and publication-grade uncertainty propagation are not automatically supplied by the presence of prototype modules. Keep their work packages separate and inspect exact-commit science evidence.

## Interaction and scientific principles

- Give the selected galaxy, release, band, source, operation state and evidence class a stable visible identity.
- Keep the FITS science image prominent. Preserve access to per-pixel values, masks, units and celestial coordinates, rather than relying on display stretch or RGB color.
- Make new-project creation and existing-project recovery visible before entering the Explorer.
- Keep preview, candidate, URL probe, downloaded, inspected, diagnostic and science-adopted statuses distinct.
- Never substitute zero or an extrapolated value for missing/unusable pixels without an explicit science-reviewed model.
- Do not infer independence of DR9/DR11 CCD exposures merely from release identifiers.
- Keep source FITS and historical project identifiers immutable; make diagnostics separate, identifiable outputs.
- Provide meaningful error, offline, loading and cancelled states, without claiming unsupported progress or features.
- Any network fetch of remote science data requires the documented user workflow. A selected-object RGB preview can use the separately documented automatic display cache, which the user may disable.

## Scope control and release gates

Current priorities and ownership belong to [GitHub roadmap #44](https://github.com/frigidable/Desidera/issues/44). Software CI, offline Qt screenshots, independent real-FITS verification, on-device GUI review and quantity-specific scientific acceptance are distinct gates. Never use an obsolete phase plan, experimental branch name or screenshot alone to mark a new capability science-ready.

GUI-only changes must not silently modify WCS, FITS inputs, photometric mathematics, transfer behavior, masks, uncertainty assumptions or results adoption.
