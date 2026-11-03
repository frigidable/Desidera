# PSF-0 first observational pilot: Leo P

Date: 2026-11-03. Status: **numerical, experimental pipeline implementation**.
This document defines executable input and output contracts; it is not
a report of a successful analysis on Leo P FITS. Synthetic tests and the
Windows/Linux CI may run; user-local native FITS must be verified separately.

## The scientific pilot is Leo P, not NGC2146sat

Leo P is the **first mandatory real-object control**. NGC2146sat is a
secondary methodological comparison because a researcher-built empirical
PSF already exists there; it is not the first pilot and must not replace
Leo P in script defaults, test documentation, or acceptance gates.

Working saved Desidera project center:
- ICRS (155.438012 deg, +18.088025 deg)
- Selected observational set per the user's current GUI: DR11 South,
  brick 1554p180, r-band (4/4 local FITS displayed in the earlier
  project screenshot; the repository does not contain these local FITS).
- Additional g/z are **candidates only**, verified by file existence
  and then independently through the strict original-FITS adapter;
  never create a missing band, download it implicitly or claim coverage.
- Literature positions may differ by tens of arcseconds depending on
  adopted center/definition. The script checks the project center
  against the recorded working center within 45 arcsec; this does NOT
  validate physical membership, morphological center or photometric
  ellipse. Leo P has resolved stars and irregular star formation,
  so target masking must not use the bright knots as an automatic
  empirical PSF training sample.

The Leo P source identity is checked *before* any numerical pilot runs.
A generous initial target protection area (at least 90 arcsec radius,
configurable and requiring researcher review) is used to avoid selecting
Leo P stars as foreground PSF trainers. This preliminary circle is not
a scientifically validated target extent.

## Numerical core now provided

`src/lsbpipe/photometry/psf_empirical.py`:
1. Requires an existing `NativeBandFrame` from the strict quartet
   adapter. Never accepts a remote Viewer RGB image or cached screen
   resolution in place of photometry.
2. Generates conservative target protection around the explicit
   project ICRS center using local WCS (small-field approximation).
3. Computes initial source-detection sky and RMS from valid pixels
   outside protection; this is **not final galaxy sky**.
4. Uses `photutils.DAOStarFinder` for candidate positions with a
   user-controlled approximate FWHM and absolute detection threshold.
   Detection does not by itself certify star/galaxy classification.
5. Applies saturation/ALLMASK overrides, native validity,
   sufficient unmasked cutout fraction, crowding, star core contrast,
   sky-annulus quality and local sky-gradient rejection.
6. Fits local planar sky to each star's annular pixels using robust
   iterations. This is deliberately separate from galaxy residual
   sky, and may still absorb wide PSF wings; that limitation is flagged.
7. Builds native-grid empirical core ePSF using photutils
   `EPSFStar`/`EPSFBuilder` from independently sky-subtracted
   cutouts, explicitly weighted by native diagonal INVVAR.
   The configured native sampling is 1. Oversampling 2/4 is **not**
   silently claimed: adaptive oversampling and its star-count tests
   remain a follow-up requirement.
8. Holds out spatially sorted sources from model construction and
   refits only their amplitude to the empirical kernel; reports
   median core absolute residual relative to peak. This is a
   **diagnostic** rather than calibrated chi-square or independent
   high-precision PSF error budget.
9. Keeps training and held-out stars with sky positions and local
   sky metrics; rejects a result with insufficient sample or severe
   negative fitted-PSF ringing.
10. Emits `EXPERIMENTAL_PSF_NOT_ADOPTED` and warning flags for
    unvalidated wings, local spatial representativity, correlated
    coadd noise, unverified Gaia identity, preliminary sky and
    potential PSF reconstruction nonconvergence.
11. If and only if input quartets carry the strict loader's local
    SHA-256 evidence, creates a dataset/band/instrument/selection/sky/
    recipe-bound `PSFLineageKey`. A digest is **not** a quality
    score or a proof of remote archive provenance.
12. On explicit `--save-diagnostic` only, creates a versioned
    `analyses/psf-pilots/psf0-<hash>/` with an original-normalized
    native kernel FITS, QA JSON and caution. It refuses overwrite
    and never touches original FITS. It does not insert the output
    into the scientific accepted-results registry.

## Primary workflow: existing Windows GUI and galaxy projects

The only supported day-to-day developer and researcher path is the existing
`C:\\Users\\frigi\\Documents\\PhD\\Desidera` checkout and its
existing `.venv`. Do **not** create a separate worktree, venv, or copy
of the project. The existing Leo P project is opened via Galaxy Projects.

1. Run the established `lsbpipe.gui` application from the same
   checkout, with the existing virtual environment.
2. In **Galaxy Projects**, open the already saved **Leo P**.
3. Enter **Galaxy Workspace → Photometry**.
4. In **04 · Empirical PSF**, choose an actually present band
   (e.g., `r`, even if the saved Galaxy Workspace default is `g`),
   review the initial guessed detection width and the explicit
   conservative target protection radius.
5. Click **Build PSF preview**. Existing background task controller
   validates real original FITS/WCS/units and computes a diagnostic
   ePSF without writing an output FITS or modifying the project selection.
   GUI stays responsive, with existing Task Stop semantics; the
   photutils fitting phase itself is not presently interruptible.
6. If successful, inspect training vs independent holdout stellar
   positions, reconstructed native-grid PSF image, residual diagnostic
   and all scientific flags. The result is **experimental**, not
   adopted science.
7. Optionally click **Save diagnostic FITS + report**. An independent
   GUI task rehashes all four original files to detect changes,
   then publishes an immutable experimental derivative within
   the SAME Leo P project directory. Previously saved outputs are
   never overwritten.
8. If the model is scientifically rejected (e.g., strong negative
   ePSF ringing), inspect the rejected star sample, the negative
   mass fraction and background/selection flags directly in the GUI;
   no invalid PSF kernel is accepted or exported. No automatic
   weakening of QA thresholds is performed.

The user's first actual run on 2026-11-03 found Leo P, with the saved
selected `DR11 South / 1554p180 / g` band. The requested **r**
original quartet exists independently. That real-r run returned
`PSF0_FAILED_NOT_ADOPTED` due to
`Strong negative EPSF ringing; inspect star selection/sky`.
It is therefore **not** a validated observational PSF.

**Real g-band GUI evidence (researcher screenshots, 2026-11-03):** all
four FITS paths were found; the initial empirical reconstruction used
74 training stars and 19 independent holdout stars, with 28 isolation,
5 pixel-quality/target and 1 image-edge rejection. The signed
ePSF had 4.09% negative absolute mass against an **experimental,
uncalibrated 2.50% screening threshold** and was correctly rejected.
These observations are not a measured FWHM, a validated PSF or evidence
of unusable stellar data; they motivate radial, per-star and
solver-version diagnostics.

The new GUI examines **unmodified signed ePSF values**, separates
negative mass into core/intermediate/outer zones defined using the
*assumed detection FWHM* (not a measured FWHM), records the effective
photutils version/solver defaults, shows a cyan-positive/red-negative
diagnostic image, and offers **Copy QA report** without writing files
or adopting a science result. Do not lower the 2.5% provisional
cutoff or adopt the clipped kernel simply to pass this diagnostic.

Photutils 3.1 changed default ePSF smoothing and center-fit boxes;
at a native-grid assumed FWHM of roughly 4.8 px its auto smoothing
may apply no smoothing, as documented by photutils. This is a
**hypothesis requiring real g/r star-cutout verification**, not a
demonstrated cause of the observed 4.09% negative mass.
See https://docs.astropy.org/projects/photutils/en/stable/whats_new/3.1.html
and https://docs.astropy.org/projects/photutils/en/stable/user_guide/epsf_building.html .

After the stacked PRs are reviewed and available in the existing
checkout, start the GUI with the existing interpreter:

```powershell
Set-Location "C:\\Users\\frigi\\Documents\\PhD\\Desidera"
$py = ".\\.venv\\Scripts\\python.exe"
& $py -m lsbpipe.gui
```

If `photutils` is not installed in the **existing** `.venv`, install
the project's current optional science dependencies there; do not
create a new `.venv`. The code on the current unmerged PR branch is
not obtained merely by running `git pull origin main`.

## Optional legacy CLI diagnostic (not required for GUI work)

From the existing Desidera repository, with dependencies installed:

```powershell
Set-Location "C:\Users\frigi\Documents\PhD\Desidera"
git status -sb
# After approved PRs are merged, sync main with --ff-only.

$py = ".\.venv\Scripts\python.exe"
& $py scripts/psf0_leo_p.py --project-root "C:\path\to\Leo_P" --bands r
# Read-only by default: produces JSON stdout but writes no derivative FITS.

# After reviewing the pilot:
& $py scripts/psf0_leo_p.py --project-root "C:\path\to\Leo_P" --bands r g z --save-diagnostic
```

Do **not** enter a guessed Leo P folder path. The actual project
directory must be selected from Galaxy Projects / its saved manifest.
The script will not download missing bands and will stop without
promoting synthetic or source-unverified data.

CLI flag `--fwhm-pix` sets an explicit *detection guess*, not a
measured PSF width. The source-image scale in arcsec/pixel must be
read from WCS; never assume 0.262 exactly on arbitrary CCD files.

## Controlled post-failure diagnosis: Leo P g-band (PSF-0)

The actual GUI screenshot provides **4.09% signed negative absolute ePSF
mass entirely at radii >2 × the guessed 4.8 px FWHM**, with **zero**
reported contribution from the two inner zones. The source star cutout
background was estimated in an annulus at 10–15 pixels. Since the
reported outer zone begins at 9.6 pixels, overlap makes outer-wing
sky subtraction a **high-priority hypothesis**, not an established
cause. The solver also reports `converged=False` after 8 iterations
in photutils 3.1.1. Those are separate failure indicators.

The GUI now provides four reproducible diagnostic experiments:

| Experiment | Change | Scientific question |
|---|---|---|
| `baseline` | 31px stamp, sky 10–15px, isolation 32px, 8 iter | Reproduce the original 4.09% g result |
| `outer_sky` | 51px stamp, sky 19–24px, isolation 55px, 8 iter | Are the outer negative wings sensitive to a more distant local sky? |
| `wing_smoothing_off` | 31px, default sky/stars, disable photutils final wing smoothing | Is the photutils 3.1 final wing smoothing contributing to the signed outer residual? |
| `longer_convergence` | 31px, default sky, 16 iter | Does the ePSF approach centroid/shape convergence with more iterations? |

Each run is still **non-adopted** and uses the existing Leo P FITS
without downloads, new folders/environments or permanent project-band
changes. Diagnostic results and cumulative same-band attempt
summaries remain visible in the Photometry Workspace; **Copy QA
report** includes the current detailed result and recent attempt
summary. Neither the comparison history nor a lower negative
fraction certifies the model.

The `outer_sky` mode alters the cutout, sky annulus **and isolation
footprint** together. Different stars may be included; the outcome is
therefore not a clean single-factor causal test. Check stellar
overlap, number of accepted/heldout stars and their ICRS positions
before interpreting its effect on negative PSF mass.

Note: in photutils 3.1, `smoothing_kernel='auto'` may return
`smoothing_kernel=None` for well-behaved native-grid PSF widths,
while **independent final wing smoothing may still be enabled**.
The controls report these two mechanisms separately. Official API:
https://docs.astropy.org/projects/photutils/en/stable/api/photutils.psf.EPSFBuilder.html .

The diagnostic computes additional **independent signed-kernel
holdout residuals** before any negative-value clipping, in the same
core/intermediate/outer radial zones. For each held-out star, fit
only the positive source amplitude in the core while keeping the
unmodified signed candidate ePSF; evaluate signed and absolute
residuals divided by the observed peak outside the fitted core.
Then summarize the medians across held-out stars per zone. This
guards against an apparently cleaner-looking PSF which predicts
real star wings less accurately. These diagnostics are not a
noise-whitened chi-square or a calibrated PSF-leakage budget.

Interpret candidate changes in this order: (a) sufficient and
comparable clean stars, (b) convergence fraction, (c) signed mass
distribution AND independently held-out *outer-zone residuals*,
(d) PSF spatial/wide-angle representativity, (e) estimated
galaxy-profile bias. A variant that merely passes the provisional
2.5% mass threshold **must not** be automatically adopted or
ranked above a lower-bias alternative.

## Blocking issues before publishing Leo P science

These items are required, even after successful synthetic unit tests:

- Run on the actual Leo P **original** DR11 South r-band quartet,
  record source hashes, image size and which FITS paths exist.
- Inspect all detected stellar candidates (training, holdout and
  rejected) on native FITS. Gaia evidence may be used for
  cross-identification; this is not yet automated.
- Assess which field is large enough: a small galaxy cutout may
  contain too few isolated stars for an empirical PSF. In that case,
  do not silently shrink the training/holdout thresholds—use a
  verified whole-brick reference or explicitly documented survey PSF.
- Check galaxy protection conservatively against real Leo P diffuse
  light and stellar population; don't mask real galaxy members as
  external PSF stars.
- Validate ePSF centroid, subpixel phase coverage, field variation,
  FWHM/ellipticity/radial encircled energy, 2D residuals and heldout
  stars. Current core result alone does not do those tests.
- Build extended wings with more suitable bright stars, distinguish
  the **remaining** DECam halo from survey-subtracted light.
  Do not model/subtract the full halo a second time.
- Test PSF choice on injected extended/diffuse Leo P-like galaxies,
  colors and cross-band PSF matching; propagate systematics to the
  scientific result. Only after those gates may a PSF be adopted.
- Validate 2/4x ePSF oversampling with star count and subpixel phases
  if any band proves insufficiently sampled. Do not force it by default.
- Any SBF analysis requires its own spatial-frequency/noise-mask QA.

## Scientific references

- Photutils ePSF building, which requires background-subtracted,
  carefully selected, isolated stars:
  https://docs.astropy.org/projects/photutils/en/stable/user_guide/epsf_building.html
- EPSFBuilder API and convergence diagnostics:
  https://docs.astropy.org/projects/photutils/en/stable/api/photutils.psf.EPSFBuilder.html
- Legacy DR11 per-CCD PSFEx and DECam extended halo treatment:
  https://www.legacysurvey.org/dr11/description/
- Leo P early physical/structural literature:
  https://doi.org/10.1088/0004-6256/146/6/145

**No observational Leo P PSF measurement has been produced as part
of creating this PR.** Its success is evaluated only after executing
the script against the user's real, locally stored original FITS.
