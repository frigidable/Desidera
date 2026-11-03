# P0.5 — Scientific photometry input contract (native DESI Legacy coadds)

Status: code + synthetic acceptance tests in Draft PR #6. This is a
**validation interface only**. It contains no actual sky modeling, flux
integration, PSF fitting, correlated-noise calibration, isophote solver,
image matching or construction of science mosaics. The existence of
`NativeBandFrame` is **not** an assertion that a dwarf galaxy is
photometrically measurable.

## Authoritative data and conventions

Official data descriptions:

- [Legacy Surveys DR9 files](https://www.legacysurvey.org/dr9/files/):
  IMAGE is a coadded **nanomaggy per original pixel**, INVVAR is its
  diagonal inverse variance in **nanomaggy^-2**, NEXP is the number of
  individual contributing exposures. IMAGE is a Lanczos-3 resampled
  coadd, distinct from Tractor single-exposure data.
- [DR9 optical MASKBITS](https://www.legacysurvey.org/dr9/bitmasks/):
  optical mask lives in HDU1; a signed int16 MASKBITS represents a bit
  field, *not* signed numerical values. HDU2 and later may be WISE maps.
- [Astropy proj_plane_pixel_area](https://docs.astropy.org/en/stable/api/astropy.wcs.utils.proj_plane_pixel_area.html):
  a WCS pixel area computed using the pixel-scale matrix gives local
  *projection-plane area*, not exact spherical area in a distorted WCS.

Arrays are always `(y, x)` in NumPy, sky coordinates are ICRS RA/Dec in
degrees, FITS pixel coordinates are 0-based when passed to WCS methods.
Stored image values are **not divided by pixel area** on ingestion.
A later surface-brightness method must explicitly use an area model
to convert nMgy/pixel to nMgy/arcsec², then apply AB calibration with
`m_AB = 22.5 - 2.5 log10(F_nMgy)` **only for a positive integrated
flux**. Do not compute `log10` of negative background-subtracted pixels.

No extinction, k-correction, zeropoint adjustment, convolution or
sky subtraction is performed by P0.5. The upstream Legacy coadd has
already undergone its survey background handling; *residual* sky is
**unknown until independently modeled**.

## NativeBandFrame contract

`lsbpipe.photometry.NativeBandFrame` stores a single original
`release/region/brick/band`. P0.5 supports separately identified
optical g/r/z and region-available i only; North does not claim i-band
where unavailable. Every frame contains:

| Field | Numerical semantics |
| --- | --- |
| `image_nMgy[y,x]` | Original integrated flux, nMgy/native pixel |
| `invvar_nMgy2[y,x]` | Original inverse variance, 1/nMgy², **diagonal evidence only** |
| `nexp[y,x]` | Exposure-count diagnostic (never a surrogate variance map) |
| `maskbits[y,x]` | Unsigned interpreted optical bit field; signed int16 bit 15 retained |
| `wcs` | Detached celestial WCS copy, checked on image-grid round trips |
| `valid[y,x]` | Explicit computed boolean exclusion mask; original pixel values unchanged |
| `mask_reasons` | Separate per-bit and numeric exclusion diagnostics |
| `source_evidence` | Original input paths, HDUs, units, locally computed hashes when requested |
| `local_pixel_area_arcsec2` | Approximate local projection-plane area at WCS reference |
| `science_state` | `NATIVE_VALIDATED_NOT_PHOTOMETRY_READY` |

The frame takes ownership of *immutable byte-backed copies* of the
four arrays and mask reasons. Mutating a caller's original array or a
returned WCS copy cannot silently change the persisted scientific
snapshot. An explicit conversion/adapter stage must produce any
resampled or sky-subtracted derivatives.

The field does **not** pretend that native coadd INVVAR yields the
correct variance of an aperture flux: adjacent resampled pixels are
correlated. Covariance is explicitly
`DIAGONAL_IVAR_CORRELATIONS_UNKNOWN`, not a diagonal-only assumption.
Subsequent background estimation, blank-aperture uncertainty and
covariance propagation must be implemented independently and validated
with injections or appropriate residual measurements.

### Mask policy

`MaskPolicy` defaults to excluding nonfinite IMAGE, nonpositive or
nonfinite INVVAR, zero/nonfinite NEXP, and band-specific `SATUR_X`
and `ALLMASK_X`. The BRIGHT, GALAXY and NPRIMARY bits are separately
reported but **not globally rejected by default**: these flags
may indicate contaminating stars, regions to be treated by a
neighbor/star model, or a primary-brick boundary requiring overlap
management. Alternate explicit policies can reject each of them.
The user's additional exclusions remain a distinct named mask
`USER_EXCLUSION`. A mask bit is not a PSF model or proof of
contamination-free surface photometry.

### Geometry and pixel area

Inputs must have 2-dimensional celestial WCS, finite invertible local
pixel scale matrix and sub-millipixel round trips on test locations
(corners and center). Different source HDUs must match the same original
image grid to 0.001 pixel using the existing P0.2 five-sample test.
Distorted WCS that cannot be characterized by a local affine area
fails closed pending an explicit per-pixel solid-angle model.

Even for an ordinary TAN coadd, a scalar computed from the WCS Jacobian
is a **local approximation**, not exact per-pixel spherical area.
Before publication-quality surface-brightness integrations, introduce
a per-pixel or rigorously controlled area map where needed and
quantify any residual effect.

## Ingestion adapter

`load_native_legacy_quartet(data_root, release, region, brick, band)`
reads original local FITS **without archive access**. It reuses
existing P0.2 science HDU identification, WCS-grid consistency and
optional FITS CHECKSUM/DATASUM tests. Original optical MASKBITS is
always HDU1, never WISE extensions. Paths are constructed by the
established Legacy product naming rules, not user-provided relative
path traversal. Files are read-only, and input metadata is
rechecked before/after ingestion; local SHA-256 may be calculated.

- An explicitly present recognized IMAGE `BUNIT` is accepted.
- A missing IMAGE `BUNIT` **fails by default**. An explicit
  `require_explicit_bunit=False` opt-in adopts the survey
  documentation's nMgy/pixel convention and records
  `LEGACY_SURVEY_DOCUMENTATION_ASSUMPTION` as provenance.
- An explicit INVVAR BUNIT incompatible with the known inverse
  nanomaggy units fails; a missing one receives the equivalent
  recorded assumption. No unit conversion is silently guessed.
- Source hashes are **local fingerprints** for reproducibility, not
  authenticated comparisons against the remote archive. Opting out
  of hashing is explicitly recorded.
- Cancellation is cooperative; no original FITS is changed.
- Loading a successful quartet is not a scientific-quality verdict.

## Preparation and multi-brick isolation

`PhotometryPreparation` preserves a native frame plus
`SkyDescription`, `PSFDescription` and optional neighbor-model raster.
All sky/neighbor raster units are nMgy per original pixel, with
matching shapes. Sky defaults to `NOT_ESTIMATED`, *never implicitly
zero*. PSF defaults to `UNKNOWN`, *never an assumed Gaussian*.
When supplied, empirical PSF kernel must have odd dimensions,
nonnegative finite entries and unit normalization; the presence
of an empirical kernel does not certify spatial PSF uniformity.
`pending` explicitly enumerates unmodeled sky, PSF, correlated
noise and the absent photometry algorithm.

`NativeFrameCollection` preserves each
`(release, region, band)` group separately, and reports
`NATIVE_UNMOSAICKED_NO_REPROJECTION`. It refuses
`require_single_native()` if there are several bricks. North/South
may be **compared**, but never automatically coadded as statistically
independent measurements. No interpolation, variance propagation,
effective PSF transport, overlap selection or flux-conserving
mosaic is implied by the collection.

## Future scientific gates

P0.5 defines interfaces, not a certified full data model. Subsequent
scientific stages must separately implement and validate:

1. Sky and contaminant models with confidence regions and
   independent residual/blank-aperture diagnostics.
2. WCS-aware masked apertures (including fractional pixel overlap),
   area weighting, corrected integrated fluxes, growth curves,
   Sérsic fits and uncertainty budgets.
3. Spatially varying empirical PSFs, common-band convolution when
   necessary, and star-halo/diffraction-ray handling.
4. Coadd noise correlation/covariance calibration; simultaneous
   sky/PSF uncertainty and controlled negative-flux treatment.
5. A formal multi-brick footprint and flux-conserving mosaic API
   with per-output-pixel coverage, weights, mask transport,
   provenance and covariance/PSF propagation.
6. Propagation of input hashes and preparation signatures into
   the P0.3 scientific execution graph.
7. Source-level photometric calibration, Galactic extinction and
   other astrophysical corrections as distinct reversible
   transformations.
8. Real FITS acceptance on UGC 2162 and several unrelated LSB
   galaxies; synthetic CI alone is insufficient.

**No P0.5 input should yet be labelled SCIENCE_READY.**
