# DESIDERA

**Experimental galaxy-analysis workbench** for DESI Legacy Surveys FITS data, local galaxy projects, exploratory photometry and reproducible scientific QA.

**Scientific status:** research prototype. A feature being implemented, a successful test or a saved diagnostic does **not** certify a calibrated astronomical measurement. No PSF-based or SBF distance measurement, background estimate, reconstructed total flux or catalogue association is adopted automatically.

**Performance development candidate:** `v0.2.1.dev2` is in existing Draft PR #84 (not tagged or merged). It retains previous exact-operator improvements to native frame loading, radial masks, blank-aperture geometry, Tractor WCS screening and sky fitting. The new audit adds reuse of exactly verified Stage03 distance-transform guards during fixed grid trials, memoization of immutable Stage04 WCS/header fingerprints, memory-conscious broadcasting of native sector coordinates, and more efficient Qt Tractor catalogue-list updates. These are engineering optimizations; absolute real-FITS speedup and peak process RAM have not been benchmarked. See [the performance audit](docs/architecture/PERFORMANCE_AUDIT_20261113.md). Reproduce seven synthetic computational-kernel comparisons offline with `python scripts/benchmark_desidera_performance.py --size 512 --repeats 3`. Scientific adoption remains on HOLD.

## Quick start

Requires Python 3.10+ and an editable checkout:

```powershell
python -m pip install -e ".[exact,gui,knowledge,test]" pytest-qt
python -m desidera.gui
```

The desktop application opens the **Galaxy Project Hub**, where a researcher can create an object-centred project, reopen an existing galaxy folder and continue from its locally saved state. Photometry Workspace runs **core + holdout empirical PSF QA by default**; the historically intensive same-brick P15–P20 analysis is available only through an explicit checkbox in **Advanced diagnostics**. The complete scientific pilot remains archived; an unrun audit is never reported as passed. The current project-creation form is a **single-page dialog**, not the older three-step wizard described in historical design notes.

The main source tree remains `src/lsbpipe/` for compatibility with existing research projects and imports. `src/desidera/` contains public package/launcher compatibility code; the installable distribution is called `desidera`. Both `python -m desidera.gui` and `desidera-gui` use the same Project Hub launcher. Historical `lsbpipe` imports, `python -m lsbpipe.gui` and on-disk project identifiers are intentionally preserved.

## Implemented engineering capabilities (not science release claims)

| Area | Current scope |
| --- | --- |
| **Galaxy projects** | Local create/import/relink/reopen operations, recent-project registry and preserved older manifests |
| **Survey navigation** | Explicit DR9/DR10/DR11 search/brick selection, candidate source URLs and download management |
| **Original FITS** | IMAGE, INVVAR, NEXP and MASKBITS dataset inspection, WCS-aware views, pixel-level coverage and diagnostic masks |
| **Prototype photometry** | Local fixed-ellipse integration, signed flux/valid-area accounting and `measure-fits` output |
| **Empirical PSF (Leo P pilot)** | Original DR11 South native coadd, Gaia/Tractor star selection, finite 81×81 g-band ePSF, held-out-star and native blank-sky QA; wide wings and science use **remain HOLD** |
| **Research QA** | Diagnostic scenes, residual-background candidates, sensitivity checks, run/lineage evidence and immutable configurations |
| **Catalogues / previews** | Opt-in SIMBAD/NED searches with user-confirmed identities, local FITS thumbnails and display-only Legacy RGB previews |

These capabilities have different levels of validation. Read their code and tests before treating a particular mode or science quantity as accepted. The formal research backlog and adoption gates are tracked in [issue #44](https://github.com/frigidable/Desidera/issues/44) and its work packages; [unmerged PRs](https://github.com/frigidable/Desidera/pulls) may contain additional prototypes.

### Current integration checkpoint

The complete scientific/UI implementation from PRs #39–#43 is now in `main`. The user-facing photometry workspace displays the real PSF workflow rather than eight non-executable planning cards, and keeps historical PSF recipes inside advanced diagnostics. The registered PSF-0 experiment remains **Leo P only**, with a real archived **g-band brick 1554p180** pilot; the existence of g/r/i/z tabs is not evidence of per-band science validation. **No empirical wing PSF or scattered-light correction is science-adopted.**

For exact merge SHAs, QA scope, CI evidence and classification of historical remote branches, read the [2026-11-09 integration checkpoint](docs/validation/REPOSITORY_CHECKPOINT_2026-11-09.md). The detailed [PSF core/wings audit](docs/qa/PSF_CORE_WINGS_VALIDATION.md) remains the scientific source of truth.

### Network and data safety

Creating or reopening a project must not redownload its science FITS. Archive downloads require the documented download actions. A display-only RGB preview **may be requested automatically for the selected galaxy** when it is absent from cache; disable this with `DESIDERA_DISABLE_AUTOMATIC_RGB=1`. Remote Viewer tiles and RGB previews are not calibrated science inputs.

Original FITS and existing project manifests must remain intact. A saved candidate, completed download or visually plausible image is not proof of science-grade coverage, WCS compatibility, background quality or measurement adoption.

## Command-line tools

```powershell
python -m desidera --help
python -m desidera demo --output demo_output
python -m lsbpipe.surveys --help
```

The `desidera` CLI exposes the legacy reproducible synthetic demonstration and `measure-fits` for pre-existing, appropriately matched Legacy coadd FITS. To generate a **synthetic test triplet**, not observational data:

```powershell
python -m pip install -e ".[exact,demo,test]"
python scripts/make_mock_legacy_fits.py --output mock_fits --size 128
python -m desidera measure-fits --image mock_fits/mock-image-r.fits --ivar mock_fits/mock-invvar-r.fits --maskbits mock_fits/mock-maskbits.fits --release dr9 --band r --ra 96.0619 --dec 77.9684 --q 0.63 --pa 127 --max-radius 8 --step 1 --background 0.0003 --backend exact --output mock_output
```

The `--background` input is **caller-supplied residual sky in nMgy/pixel**, not an automatically estimated sky model. Sky PA (`--pa`) is measured from north through east; radial distances are arcseconds. Outputs describe measured **valid-pixel flux**, not reconstructed whole-galaxy flux. The diagonal inverse-variance errors do not capture all correlated-coadd uncertainties. Read [photometry input contracts](docs/architecture/PHOTOMETRY_INPUT_CONTRACT_P05.md) before applying this to real telescope data.

## Testing and contribution

```powershell
python -m pip install -e ".[exact,gui,knowledge,test]" pytest-qt
python -m pytest -q tests --ignore=tests/gui
python -m pytest -q tests/gui
```

These are local commands, not a claim that they were executed on a given machine. Consult [GitHub Actions](https://github.com/frigidable/Desidera/actions) for exact-commit CI results. Offscreen Windows Qt captures do not replace manual Windows/Legacy Viewer/FITS inspection, and successful tests do not imply scientific acceptance.

See [contribution and scientific-safety rules](AGENTS.md), the [documentation and evidence index](docs/README.md), and the [current roadmap](https://github.com/frigidable/Desidera/issues/44). Historical phase reports remain accessible under `docs/`, but they describe their dated branch/commit, not necessarily the present `main`.

## Scientific references

- [DESI Legacy Surveys DR9 data format](https://www.legacysurvey.org/dr9/files/) and [bitmask semantics](https://www.legacysurvey.org/dr9/bitmasks/)
- [DESI Legacy Surveys DR11 data format](https://www.legacysurvey.org/dr11/files/) and [bitmask semantics](https://www.legacysurvey.org/dr11/bitmasks/)
- [Astropy Cutout2D](https://docs.astropy.org/en/stable/api/astropy.nddata.utils.Cutout2D.html) and [Photutils aperture methods](https://photutils.readthedocs.io/en/stable/user_guide/aperture.html)
