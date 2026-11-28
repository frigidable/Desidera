# OUTER-PHOT / INJ-3 (conditional)
Use the same diagnose-project invocation. Extended profiles on the existing
expanded native field, common g/r photometric support and fixed GEO-1
geometry. No new sky-only exclusions applied to the target flux mask.
Signed sums, nominal area coverage (including off-frame area), eight
elliptical sectors and vertical-strip sensitivity. Contiguous support
screen: ring/cumulative coverage >=70%, every sector >=50%; two final
signed growth increments <=2% is descriptive only, not a calibrated plateau.
Reference radial halo and matched hypothesis refits with 75/90 arcsec
protection; additional sky guards 60/75/90 arcsec. Sky candidates may
only be evaluated in an enclosed circular region containing the ellipse.
Report 40/50/60/70 arcsec finite-core Sersic fits as support permits.
Reference reporting convention: constant g, plane r, guard 75 arcsec,
40 arcsec outer/3 arcsec inner. Not a model adoption or statistical error.
INJ-3 uses three target-like profiles n=.4/.5/.7, Re=22/26/29 arcsec,
F40 ratios .85/1/1.15. Coverage and conservative training-hull preflight,
separate 40 and 60 arcsec experiments, up to two and one sites respectively.
Actual refits may still block. Frozen masks, radial halo only, finite-core
PSF, no extra Poisson noise. Sensitivity table separates model/window/guard
ranges and injection error ranges; no quadrature or calibrated uncertainty.

Native/expanded reference overlap is checked on identical observed pixels.
The reference 75 arcsec guard retains all radius windows; other guard and
halo trials use 40 arcsec and the largest supported predeclared radius.
For each fit, outer rings outside the fitted interval are predicted on
observed support up to the screened limiting radius. These are withheld
from that profile fit, not necessarily independent of upstream masking/sky.
Display of interpolated sky corrections is restricted to the verified disk.


## Scenario-local support correction (NGC2146sat, QA 2026-11-27)

The previous global intersection included the protect-90 halo, whose
g-band lower training radius is 37.967561 arcsec from OBJID 360. Native
and expanded reference finite-pixel masks reconstructed from the report
PNGs agree inside the 40 arcsec ellipse (40924 pixels). Applying that
radial exclusion reproduces all 4448 missing pixels and the 36476-pixel
intersection exactly. No hard/source mask or sector threshold is relaxed.

Use one common band mask per halo scenario; preserve the SH1-C ensemble
radial baseline. The reference is no longer vetoed by alternative halo
model domains. Report per-scenario hashes, radial aperture counts,
nonfinite reference-pixel counts, and the old global intersection loss.
A refit with missing sectors remains blocked. Compare refits with the
reference at the same sky guard on their pairwise intersection, including
signed ring differences and reference fits where the same support passes
the screen. Sensitivity ranges exclude different full-crop mask hashes.
No halo extrapolation or calibrated uncertainty is introduced.
