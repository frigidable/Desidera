"""Photometric input contracts and experimental geometry/mask drafts."""
from .contracts import (
    NativeBandFrame, NativeFrameCollection, MaskPolicy, PSFDescription,
    SkyDescription, PhotometryPreparation, CONTRACT_SCHEMA,
)
from .native_fits import load_native_legacy_quartet

__all__ = [
    "NativeBandFrame", "NativeFrameCollection", "MaskPolicy",
    "PSFDescription", "SkyDescription", "PhotometryPreparation",
    "CONTRACT_SCHEMA", "load_native_legacy_quartet",
]

# Experimental preparations are not adopted photometric results.
from .sky_masks import (
    SkyMaskShape, SkyMaskEdit, SkyMaskDraft,
    rasterize_sky_shape, evaluate_mask_draft,
)
from .mask_drafts import read_active_mask_draft, save_mask_draft
from .ellipse_seed import (
    EllipseSeed, propose_isophotal_seed, seed_as_sky_geometry,
)
__all__ += [
    "SkyMaskShape", "SkyMaskEdit", "SkyMaskDraft",
    "rasterize_sky_shape", "evaluate_mask_draft",
    "read_active_mask_draft", "save_mask_draft",
    "EllipseSeed", "propose_isophotal_seed", "seed_as_sky_geometry",
]

# Human-seeded stellar segmentation proposals; no automatic science adoption.
from .smart_masks import (
    StarMaskSeed, SmartStarMaskCandidate, propose_seeded_star_mask,
    SMART_MASK_STATUS,
)
__all__ += [
    "StarMaskSeed", "SmartStarMaskCandidate", "propose_seeded_star_mask",
    "SMART_MASK_STATUS",
]

# Planning-only workflow DAG and observation-bound PSF cache identity.
from .workflow_dependencies import (
    STAGE_ORDER, DEPENDENCIES, invalidated_descendants,
    require_stage_dependencies, PSFLineageKey,
    inspect_empirical_psf_coverage, require_empirical_psf_coverage,
)
__all__ += [
    "STAGE_ORDER", "DEPENDENCIES", "invalidated_descendants",
    "require_stage_dependencies", "PSFLineageKey",
    "inspect_empirical_psf_coverage", "require_empirical_psf_coverage",
]

# Experimental per-band native-coadd empirical PSF: holdout QA, never adopted.
from .psf_empirical import (
    PSF0Config, PSF0Star, PSF0Result, build_empirical_coadd_psf,
)
from .psf_export import export_psf0_diagnostic
__all__ += [
    "PSF0Config", "PSF0Star", "PSF0Result",
    "build_empirical_coadd_psf", "export_psf0_diagnostic",
]


# Stage 04.0–04.1 read-only finite signed native-aperture measurements.
# photutils/astropy remain lazy numerical requirements: GUI imports do not
# acquire a global FITS/photometry dependency or execute measurements.
from .aperture_flux import (
    ApertureGeometryICRS, ApertureSkyScenario, ApertureProvenance,
    ApertureFluxRecord, PixelEllipse,
    project_icrs_ellipse_to_native, measure_native_aperture,
)
__all__ += [
    "ApertureGeometryICRS", "ApertureSkyScenario", "ApertureProvenance",
    "ApertureFluxRecord", "PixelEllipse",
    "project_icrs_ellipse_to_native", "measure_native_aperture",
]


# Stage 04.2 fixed-ellipse curves of growth and exact annular diagnostics.
from .aperture_growth import (
    GrowthHypothesis, AnnulusFluxRecord,
    geometry_hypotheses, predeclared_growth_radii,
    measure_native_annulus, measure_growth_diagnostic,
)
__all__ += [
    "GrowthHypothesis", "AnnulusFluxRecord",
    "geometry_hypotheses", "predeclared_growth_radii",
    "measure_native_annulus", "measure_growth_diagnostic",
]

# Stage 04.3: exploratory sky/mask/geometry sensitivity and native-coadd
# blank-aperture covariance. No science error-bar adoption.
from .aperture_systematics import (
    estimate_blank_annulus_covariance, diagnose_aperture_systematics,
)
__all__ += [
    "estimate_blank_annulus_covariance", "diagnose_aperture_systematics",
]
