"""No PSF0 heuristic may silently become a validated scientific limit."""
import numpy as np
import pytest

from lsbpipe.photometry.psf_parameter_audit import (
    audit_psf0_parameters, audit_star_sky_ring_sensitivity,
    audit_gaia_catalog_seeds_vs_dao,
)


def test_all_existing_psf0_thresholds_explicitly_classified():
    pytest.importorskip("astropy")
    from lsbpipe.photometry.psf_empirical import PSF0Config
    from lsbpipe.photometry.psf_gui_backend import psf0_experiment_recipe
    cfg = psf0_experiment_recipe(6.57, "outer_sky_isolation_footprint")
    report = audit_psf0_parameters(cfg)
    assert set(report["parameters"]) == set(vars(cfg))
    assert report["parameters"]["sky_annulus_inner_pix"]["basis"] == "HEURISTIC_GEOMETRY"
    assert report["parameters"]["guessed_fwhm_pix"]["basis"] == "OBSERVATIONAL_PRIOR_OR_MANUAL"
    assert report["geometry"]["ring_circles_disjoint_at_isolation_distance"] is True
    assert report["geometry"]["square_stamps_can_overlap_at_isolation_distance"] is True
    assert report["automatic_parameter_science_certification"] is False
    assert report["parameters"]["min_holdout_stars"]["calibrated_from_original_fits"] is False


def test_sky_ring_probe_detects_injected_wings_without_autotuning():
    yy, xx = np.indices((81, 81), dtype=float)
    rr = np.hypot(xx-40, yy-40)
    # This synthetic broad halo contaminates nearer sky probes.  It is not
    # intended as a realistic Leo P or DECam PSF.
    original = 4. + 20*np.exp(-.5*(rr/12.)**2)
    original_before = original.copy()
    valid = np.ones_like(original, dtype=bool)
    output = audit_star_sky_ring_sensitivity(
        original, valid, center_xy=(40., 40.),
        inner_radius=30., outer_radius=38.,
        reference_sky_nMgy_pix=4., reference_rms_nMgy_pix=.1,
        guessed_fwhm_pix=6.57)
    assert output["status"].startswith("ORIGINAL_PIXELS")
    assert output["max_absolute_inner_ring_sky_change_in_local_rms"] > 1.
    assert output["changes_psf_selection_or_sky_model"] is False
    assert output["outer_ring_beyond_cutout_not_tested"] is True
    np.testing.assert_array_equal(original, original_before)


def test_sky_ring_probe_stable_for_exact_planar_sky():
    yy, xx = np.indices((81,81),dtype=float)
    # Tiny deterministic nonzero fluctuations keep the robust RMS defined.
    img = (3. + .01*(xx-40) - .02*(yy-40)
           + .0002*np.sin(.2*xx)*np.cos(.35*yy))
    result = audit_star_sky_ring_sensitivity(
        img, np.ones_like(img,dtype=bool),
        center_xy=(40.,40.), inner_radius=30, outer_radius=38,
        reference_sky_nMgy_pix=3., reference_rms_nMgy_pix=.1,
        guessed_fwhm_pix=6.57)
    assert len(result["probes"]) == 2
    assert all(item["status"]=="FIT_FOR_SENSITIVITY_ONLY" for item in result["probes"])
    assert result["max_absolute_inner_ring_sky_change_in_local_rms"] < .02


def test_sky_ring_probes_refuse_invalid_input():
    with pytest.raises(ValueError):
        audit_star_sky_ring_sensitivity(
            np.ones((20,20)), np.ones((10,10),dtype=bool),
            center_xy=(10,10),inner_radius=4,outer_radius=9,
            reference_sky_nMgy_pix=1.,reference_rms_nMgy_pix=1.,
            guessed_fwhm_pix=4.)


class FakeMatcher:
    dataset_id = "dr11/south/1554p180/g"
    match_radius_pix = 2.5
    rows = (
        {"ref_cat":"G3","ref_id":1,"type":"PSF"},
        {"ref_cat":"G3","ref_id":2,"type":"PSF"},
        {"ref_cat":"G3","ref_id":3,"type":"PSF"},
        {"ref_cat":"G3","ref_id":4,"type":"PSF"},
        {"ref_cat":"G3","ref_id":4,"type":"PSF"},
    )
    pixel_positions = np.array([[50.,50.],[130.,50.],[210.,50.],
                                [270.,50.],[275.,50.]])
    def match(self,x,y,*,band,neighbor_search_radius_pix,require_gaia):
        if x == 210.:
            return None, "TRACTOR_NEIGHBOR_WITHIN_PSF_CUTOUT"
        return {"gaia_source_id":int(x)}, None


def test_gaia_only_catalog_preflight_does_not_promote_or_remove_stars():
    m = FakeMatcher()
    result = audit_gaia_catalog_seeds_vs_dao(
        m, [(50.7,50.2)], frame_shape=(350,350),
        protected=np.zeros((350,350),dtype=bool),band="g",
        cutout_size_pix=81,guessed_fwhm_pix=6.57)
    assert result["n_catalog_gaia_type_psf_rows"] == 5
    assert result["n_catalog_only_preflight_pass"] == 3
    assert result["n_near_dao_centroid_not_proven_associated"] == 1
    assert result["n_catalog_preflight_without_near_dao"] == 2
    assert result["preflight_rejected_reasons"]["DUPLICATE_GAIA_SOURCE_ID"] == 1
    assert result["preflight_rejected_reasons"]["TRACTOR_NEIGHBOR_WITHIN_PSF_CUTOUT"] == 1
    assert result["no_changes_to_stellar_sample_or_holdouts"] is True
    assert result["science_adopted"] is False


def test_catalog_preflight_does_not_confuse_galaxy_with_valid_seed():
    m = FakeMatcher()
    mask=np.zeros((350,350),dtype=bool)
    mask[50,130]=True
    result = audit_gaia_catalog_seeds_vs_dao(
        m, np.empty((0,2)),frame_shape=mask.shape,protected=mask,
        band="g",cutout_size_pix=81,guessed_fwhm_pix=6.57)
    assert result["preflight_rejected_reasons"]["GALAXY_PROTECTED_AT_CATALOG_SEED"]==1
    assert result["n_near_dao_centroid_not_proven_associated"]==0
