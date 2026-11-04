# PSF-0: FWHM-scaled local-star sky geometry (experimental)

**Scope:** Leo P-first DR11 south / g,r,i,z native coadd PSF-0. The original FITS are read-only. This change is a diagnostic **recipe change**, not a statement that PSF wings/sky are solved, and not a scientifically adopted PSF. PSFSIZE values are survey-weighted seeing estimates, **not** empirical ePSF measurements.

## Diagnosed problem

For Leo P DR11 `g`, the original native coadd `psfsize-g` map gives `1.721 arcsec` at `0.2620 arcsec/pix`, a detection FWHM prior of `6.57 pix`. The former fixed Extended sky annulus `19–24 pix` samples only **2.89–3.65 × FWHM**, where stellar wings may bias the fitted local sky plane. Sigma clipping does not eliminate smooth wing light. This invalidates any automatic claim that a low stellar residual confirms correct background removal.

Official photutils guidance: https://docs.astropy.org/projects/photutils/en/stable/user_guide/epsf_building.html . DR11 seeing/FWHM: https://www.legacysurvey.org/dr11/files/ .

## Active, exactly three, deliberately non-adopted recipes

| Mode | Cutout / annulus / isolation | Interpretation |
|---|---|---|
| `baseline` | fixed 31 pix / 10–15 px / 32 px | Retained *as an explicitly wing-sensitive reference only*, NOT correct sky |
| `outer_sky` | adaptive `r_inner=ceil(4.5 FWHM)`, `r_outer=ceil(5.75 FWHM)`, `size=max(51,2*(r_outer+2)+1)` odd, `isolation=size+4` | Conservative expanded local sky sampling, **not proof of sky free from wings** |
| `outer_sky_robust` | same geometry, plus ≤2 passes of *guarded* train-only residual screening | Screen is still disabled automatically when it tracks star amplitude, and holdouts are not used to choose rejections |

At `FWHM=6.57 pix`, both extended modes have `81 × 81` cutouts, `30–38 pix` sky annulus, `85 px` minimum neighbour isolation; ratios **4.57–5.78 × FWHM**. At `FWHM=4.8 pix`, extended geometry is `61 × 61` and `22–28 pix`. Both examples are based on a *guessed seeing prior*, not verified physical ePSF.

If an adequate annulus would require a cutout larger than the PSF-0 safety limit of `91 × 91`, **raise an actionable error** rather than silently moving the sky closer to the star. A baseline reference remains possible but cannot be presented as a validated replacement. The enlarged cutout and isolation constraints may reject many stars; the number of accepted train/holdout sources must be evaluated on real Leo P FITS, not simulated as a pass.

## Scientific restrictions and next validation

1. Per-star 2D signed residual, core/outer profiles, radial trend of sky annulus and octant spread remain **descriptive**, not covariance-calibrated chi-square.
2. `EPSFBuilder` requires background-subtracted input. The Desidera robust planar sky fit is upstream of the ePSF, but it is not a joint star-wing + residual sky model.
3. A `>=4.5 FWHM` annulus may still contain the wings of bright sources, especially if nearby sources or wide-angle scatter are present.
4. Changing seeing or cutout size changes the detection/selection footprint. Cross-recipe comparisons require matched Gaia source IDs, consistent input hashes, and a scientifically justified common population.
5. The batch dispatcher uses distinct `PSFSIZE_{G,R,I,Z}` priors, except explicit per-band manual overrides. A manual `star-XXXX` exclusion from another FWHM is **rejected**, not silently transferred.
6. The detailed `griz` debug JSON includes exact effective recipe + SHA256 and original source hashes for rejected and completed outputs. Strict project UUID, release/brick, and optical-band checks prevent cross-band displays.
7. Future science QA: injection–recovery under a simulated smooth wing, blend, bright/saturated star and sky-gradient grid; multi-annulus background sensitivity; stratified heldout star comparisons (magnitude, colour, position), photutils subpixel-phase diagnostics, and covariance-aware goodness of fit. None are completed here.

**Scientific status:** `EXPERIMENTAL_NOT_ADOPTED` / `DIAGNOSTIC_NOT_ADOPTED`. Input source files have not been modified.


## Leo P griz evidence (2026-11-04) and two scientifically distinct catalog policies

The real DR11 south brick `1554p180` four-band QA was inspected on the local
native FITS, with independently verified Tractor source SHA-256 and per-band
PSFSIZE priors. Under the FWHM-scaled `outer_sky` recipe:

| Band | Prior FWHM (native px) | Strict Gaia/Tractor train+holdout | Signed negative mass | ePSFBuilder converged |
|---|---:|---:|---:|---|
| g | 6.57 | selection underflow: 10 accepted, minimum 13 | not fitted | not fitted |
| r | 4.88 | 16+4 | 1.235% | yes |
| i | 5.12 | 10+3 | 2.133% | **no** (8 iterations) |
| z | 3.99 | 39+10 | 1.430% | yes |

This is observational diagnostic evidence, not a calibrated model adoption.
The g run examined 969 DAO candidates, of which 692 were first rejected by
`NOT_GAIA_DR3`; this category includes **unverified sources of all types**,
not automatically usable stars. Other leading g reasons were 109 strict
Tractor neighbours and 109 missing centroid crossmatches. The i model is
not converged even though its final fitted-star centroid convergence
fraction is 100%: the Photutils build stopped at the configured maximum,
so the complete convergence status must be used. The negative-mass figures
are measured on **signed pre-clipping kernels**, not photometric flux losses.

Two noninterchangeable candidate-training policies are now provided:

- `gaia_strict` (unchanged): original native DR11 Tractor `TYPE=PSF`
  and `REF_CAT=G3`, valid Gaia ID, masks, fitted flux, isolation, shape
  and local sky. Both training and holdout stars have Gaia-confirmed IDs.
- `tractor_psf_validated` (new **EXPERIMENTAL TRAINING POLICY**): same
  original DR11 Tractor brick and the same conservative source/cutout
  gates, with the relative core-width reference measured exclusively
  from Gaia DR3-confirmed stars, but `TYPE=PSF` stars lacking Gaia IDs
  are eligible **only for training**. The expanded mode rejects the
  sample if fewer than five Gaia reference cores establish that prior. Gaia-confirmed stars are independently reserved for the
  holdout **before** building or iteratively screening the training PSF.
  It fails closed if there are fewer than the required Gaia holdouts,
  even if many non-Gaia Tractor stars are available.

Neither `TYPE=PSF` nor a fitted-flux S/N from a coadd is an independent
spectroscopic or high-resolution proof of stellar morphology. Gaia
holdouts are independent of the *ePSF fitting and train-only trimming*, but
the initial candidate preselection (e.g. masks, core moment width) still
uses their cutout pixels. Consequently, Gaia-only holdout residuals are
not an unbiased validation of the preselection classifier. Different
training policies may yield different PSF star samples. Compare
source-ID overlap and report brightness, position, and noise dependence
before attributing residual changes to any individual gate.

A per-candidate `candidate_screening_audit` in the full griz QA JSON now
records each inspected DAO detection's position, peak, nearest DAO
neighbour, catalog-match evidence (including ref-cat, optional Gaia ID,
Tractor type, fitted flux/S/N, mask/fit bits and nearest catalog neighbour),
first refusal reason, passed gates, local-sky QA and shape QA where they
were actually reached, and its pre-model train/holdout role. **Missing
metrics mean that the earlier failure prevented that test**; they are
not imputed. The audit reports `n_screened_candidates` and
`n_unexamined_candidates`, so DAO-only diagnostic caps cannot be
mistaken for complete census. No full native pixel arrays are exported.

The training-only policy is a research control, not a recommendation to
replace Gaia-strict as the default, and it cannot be science-adopted.
The native FITS/Tractor sources remain read-only. Further acceptance
requires explicit per-band injection/recovery, sky-wing and spatial
PSF validation and covariance-aware heldout residual assessment.


## Independent real-image follow-up: Tractor-vetted pilot (2026-11-04 15:09 UTC filename)

The uncommitted local researcher export
`Desidera_LeoP_PSF0_griz_debug_20261104_150909.json` uses the exact
DR11 south `1554p180` Tractor SHA-256
`8e8e3f75af868f3ef7cf4c4c9d7ca28210cbca87ef27c98351a6ab237010562c`.
These figures are **observed diagnostics only**, not clean PSF measurements:

| Band | Total DAO examined | Train (Gaia/Tractor-only) | Gaia-only holdout | Signed negative absolute mass | Solver/QA outcome |
|---|---:|---:|---:|---:|---|
| g | 972 | 11 candidates before split (10 Gaia, 1 Tractor-only) | no model | n/a | Insufficient total 11 vs 13 required |
| r | 882 | 16/10 | 4 | 1.707% | Converged, still NOT_ADOPTED |
| i | 1427 | 10/12 | 3 | 4.825% | Not converged after 8 iterations, blocked |
| z | 675 | 36/39 | 10 | 2.719% | Converged but 2.5% experimental screen failed |

In r/i/z the **training** median native positive-moment core sigmas for
Gaia and Tractor-only are respectively
2.30/3.15, 2.34/3.18, and 2.19/2.68 pixels. The corresponding
median native **peak** fluxes in nMgy/pix are
2.879/0.094, 2.654/0.106, and 2.970/0.214. Comparing unmatched
medians is **not a causal test for galaxy contamination**:
positive-clipped moments are S/N dependent. For example, in the
z 0.2–0.5 nMgy/pix peak bin, already observed median widths
were 2.46 pixels for six Gaia-associated objects and 2.58 for
19 Tractor-only objects. The other bands have poor faint-end Gaia
coverage and hence insufficient brightness-matched controls.

The negative power remains concentrated in the model's **outer**
radii: i contributes about 4.70% of total unsigned mass beyond
2 assumed FWHM and z about 2.72% beyond 2 assumed FWHM. Negative
absolute mass is an uncalibrated shape heuristic, NOT an error rate
or negative physical light. Photutils `converged` is the actual
builder-level Boolean, not just the fraction of recentered stars;
i had a final converged fraction of 1.0 but `converged=False`
because an additional center-shift criterion was unmet.

**Instrumentation repair in this PR:** preserve exact candidate and
training-cap roles in the full QA JSON (including underflow and
cap-unused rows), record median peak flux, native peak/RMS and
positive-core moment width separately for Gaia and Tractor-only,
and provide fixed, **descriptive** peak-flux strata with sample-size
gates on the interpretation. No new automated classifier, PSF
acceptance criterion or loosened selection threshold is introduced.

**Next physical validation:** match Gaia holdout `REF_ID` and source
hashes between strict and pilot runs, perform matched-brightness
train-set comparisons where actual overlap exists, test cutout
background/ring sensitivity, and then investigate PSF spatial
variation and broad wings on real data. Do not use more Tractor-only
stars as evidence of better ePSF by itself.


## Approved PSF-page visualization contract (2026-11-04)

The approved GUI prototype at the project level requires a **persistent**
per-selected-band scientific view: a native PSF model image, a radial
profile/residual diagram and a quality/status summary. These three science
surfaces must never be hidden behind the detailed-QA disclosure. The left
Galaxy Workspace navigation and independent g/r/i/z states remain intact;
per-band configuration and griz batch/export follow the science results.

Data binding is deliberately stricter than a cosmetic illustration:

- **Successful experimental fit:** the 2D model is the *actual in-memory*
  normalized, nonnegative \`PSF0Result.kernel\`. It remains labelled as
  provisional/non-adopted; its signed pre-clipping diagnostic is preserved
  separately in QA. The radial diagram is generated from the real,
  in-memory, local-sky-subtracted native cutouts and the *signed aligned*
  model from the existing star-review worker snapshot. Each training star
  uses its measured native peak as a display normalization and identical
  1-pixel annular masks in observed and model arrays. Each displayed radius
  pools valid stars using the median, for observed intensity, matched-model
  intensity and **observed−model signed residual**.
- **Signed-kernel quality rejection:** an actual computed but
  scientifically rejected signed model may be shown in diagnostic colors
  with explicit REJECTED/non-adopted status. Radial plots are permitted
  only if the worker also supplied the actual paired cutouts/model.
- **Selection underflow / no computed fit / missing inputs:** model and
  radial cards show explicit unavailable explanations. The GUI does **not**
  invent a PSF model, normalize an unavailable model, fabricate residual
  curves or reuse a prior-band model. The complete selection failure
  remains accessible in per-band full QA.
- The radial diagram is **TRAINING-IN-SAMPLE**, not an independent Gaia
  holdout, calibrated chi-square, full PSF-wing validation, or spatial
  PSF constraint. Only stars still marked \`train\` are aggregated;
  \`excluded\` and \`auto_rejected\` cutouts cannot silently contribute.
  The original independent heldout metrics and immutable source
  provenance remain accessible through the same QA surface.

Implementation reuses \`psf_star_review.py\` rather than introducing an
additional image pipeline or directories. All images are rendered in the Qt
GUI from previously computed arrays, without network access, FITS writes,
or changes to any scientific fitting/selection thresholds.
