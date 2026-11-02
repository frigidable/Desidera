# P0.5d — Native pixel coverage diagnostics

Status: DIAGNOSTIC_NOT_ADOPTED. This phase does not create a photometric mosaic, rank source quality, or modify original FITS.

## Method

The Data Manager offers Audit pixel coverage, an explicit offline worker task. It groups local IMAGE, INVVAR, NEXP by identical source UUID, brick and filter, and uses the MASKBITS belonging to the same source and brick. No substitution between North and South, regardless of matching WCS.

Each complete quartet passes the existing original-FITS HDU, WCS, unit and pixel-shape consistency checks. Pixel counts are accumulated in bounded row stripes, using pixel-center ICRS separations from the saved object center for two regions: the target disk and the project's sky annulus.

For each region separately, the diagnostic reports the proportions with finite IMAGE, positive finite INVVAR, positive finite NEXP, and all three jointly. It also reports the observed fractions of NPRIMARY, BRIGHT, MEDIUM, BAILOUT, band-specific SATUR and ALLMASK bits separately. These flags are NOT all automatically classified as rejected pixels.

Within-image numerical fractions use only region pixels contained in that original FITS. To prevent confusing these with *full sky* coverage, a separate approximate footprint fraction is reported from the known sky area and local WCS pixel area. This is explicitly approximate, not an exact WCS-polygon integral.

Incomplete inputs remain INCOMPLETE_FILES (for example DR9 South g/INVVAR is not registered). Invalid FITS, changed original metadata and inconsistent WCS become INVALID_INPUT. The diagnostic neither reuses results of another source nor invents a percentage for missing files. At most 16 quartets and 20 million native pixels per quartet are inspected in one call; omitted counts are explicit.

The geometric map can now show all source outlines or one selected source alone, without moving sky coordinates to make coincident outlines look different.

## Scientific limitations

A numerical weight/exposure fraction is NOT LSB depth, clean background coverage, noise power spectrum, PSF matching, validated mask policy or an SBF-eligible image. Legacy coadds are resampled using Lanczos-3 and contain correlated noise. The inverse-variance map does not supply full covariance. Different release/region coadds are never treated as independent exposures by default. Source SHA receipts are distinct from the file size/mtime preflight performed here.

Official references:
https://www.legacysurvey.org/dr9/files/
https://www.legacysurvey.org/dr9/bitmasks/

## Next methodological gate

Before an automatic source recommendation or validated Mosaic Builder:
- explicit, versioned, science-specific masking semantics;
- valid-pixel maps and local blank-sky contamination;
- PSF spatial variation and extended wings;
- correlated noise and residual-background large-scale power;
- injection/recovery of diffuse artificial galaxies across source variants.
