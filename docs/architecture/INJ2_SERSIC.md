# INJ-2 and conditional Sersic extrapolation

Batch diagnose-project now verifies stored native PSF0 receipts against
the original FITS quartet, then caches INJ-2 separately. Missing PSF
evidence blocks that band's intrinsic Sersic fit and the paired PSF
experiment. Finite signed core kernels are normalized, not clipped;
their wings/spatial representativity are not validated.

The real target is fitted at fixed common GEO-1 reference geometry in
signed linear annular means, with finite-core forward convolution on
the same observed pixels. Compare windows 30/35/40 arcsec, inner 3/5
arcsec and all supported SKY-CLEAN hypotheses. Report analytic Sersic
total, enclosed observed and predicted flux separately, n, Re, mu0,
tail fraction and degeneracy/bound flags. No calibrated uncertainties
or empirical saturation claim.

Cross-convolve g by r's PSF and r by g's PSF on strict common pixel
support. Exclude any output affected by a masked pixel or by the edge
of the supported background domain. No interpolation. Record both
matched and native flux on identical output pixels, coverage and colour
only when coverage >=75%. This is a finite-core experiment, not proof
of true PSF matching.

INJ-2 uses at most three sites separated by 90 arcsec and at least
100 arcsec from the real target. Site selection uses coverage alone,
minimum 75% in central/cumulative/outer ellipse regions of every shape.
Inject three Sersic profiles (Re=23,30,38 arcsec; n=.55,.85,1.15) with
full convolved F40 equal to .8,1,1.2 times the target's observed F40
under the plane scenario (fallback zero if plane blocked). Same shapes
for both bands; band amplitudes differ. Refit original radial halo
hypothesis on control and injection, with extra 60 arcsec protection
around the fake galaxy and cleaned sky-only vetoes. No re-selection
of radial hypothesis; no new directional spike fit. Compare on same
supported pixels, refit additional sky and geometry, fit Sersic with
recovered geometry. Record absolute and differential errors. Do not
silently relax coverage, interpolate missing pixels or adopt totals.
No extra Poisson noise or source detection/mask redrafting after
injection. This remains a conditional coadd-level experiment.

No GUI changes or additional user actions. Inspect inj2-qa.json in the
existing unified archive. Real FITS and saved masks remain read-only.

The SH1 refit executes only the explicitly matched hypothesis, with the
same model code and gates as the full ensemble. A regression test verifies
that its candidate equals the corresponding full-ensemble candidate.
The last two signed curve-of-growth increments are screened against 2% of
the outer measured flux; this is descriptive only and cannot establish a
total or a noise-calibrated limiting radius. Large measured/model aperture
disagreement is flagged independently of optimizer convergence.
