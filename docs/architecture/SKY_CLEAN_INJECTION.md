# SKY-CLEAN and INJ-1 (experimental)
The existing diagnose-project command now runs a downstream validation stage
and caches it in analyses/project-diagnostics-validation.zip. No GUI actions
or new science FITS. Old Scene and SH1/SKY-2 results remain available.

Sky-only screening uses a 20 arcsec guard around each current-band saturated
Tractor centre (no flux/IVAR gate), and positive compact excess on a normalized
Gaussian high-pass (1.5 and 10 pixels, 8 robust sigma, >=3 connected pixels,
4 arcsec growth). Negative components are counted as a control. This is a
conditional heuristic, not a source catalogue or PSF wing measurement. No
photometric mask or saved draft changes. No network.

Compare before/after source screening with a common g/r mask, same tile
and aperture definitions, zero/constant/plane/quadratic sky. Before vs after
changes pixel support; don't interpret as paired predictive improvement.
A two-sided mean-minus-median block stress and empty-aperture tails are
reported, not used to adopt or calibrate a model.

Up to two spatially separated injection sites are selected by coverage only,
at least 100 arcsec from the real target, with >=75% coverage within 40 arcsec.
At each site inject two observed exponential profiles (Re=15,25 arcsec,
mu0=24.5,25.5 mag/arcsec2; q=.56,.70; PA=44.5,75 deg). They are sampled
on the native grid and extend into the sky training domain. Use the same
frozen mask for both control and injection and exclude 45 arcsec around
the injected centre from sky training. Refit sky and GEO-1, record signed
aperture recovery on true fixed observed support and absolute geometry errors.
Absolute recovery includes pre-existing field structure, differential
recovery measures the sky-fit transfer on the same pixels. Do not interpret
the nuisance fit half-intensity radius as effective radius. No upstream SH1
refit, intrinsic PSF convolution, source detection after injection, injected
Poisson noise or total-flux/uncertainty adoption. If no suitable site exists,
report it rather than relax coverage. Results remain NOT_ADOPTED.
