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
