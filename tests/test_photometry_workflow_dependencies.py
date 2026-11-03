"""Contract tests for PSF/scene/sky stage order and selective invalidation."""
from dataclasses import replace
import pytest

from lsbpipe.photometry.workflow_dependencies import (
    DEPENDENCIES, STAGE_ORDER, PSFLineageKey,
    invalidated_descendants, require_stage_dependencies,
)


def _digest(character: str) -> str:
    return character * 64


def _psf(**changes) -> PSFLineageKey:
    record = dict(
        observation_id="dr9/north/1554p180/g",
        instrument="90prime", band="g",
        source_sha256=_digest("a"),
        star_selection_sha256=_digest("b"),
        cutout_sky_sha256=_digest("c"),
        recipe_sha256=_digest("d"),
        image_kind="coadd", model_scope="global")
    record.update(changes)
    return PSFLineageKey(**record)


def test_scientific_order_prepares_stars_and_sky_before_psf():
    assert STAGE_ORDER.index("scene_initialization") < STAGE_ORDER.index("preliminary_sky")
    assert STAGE_ORDER.index("preliminary_sky") < STAGE_ORDER.index("psf_reconstruction")
    assert STAGE_ORDER.index("psf_reconstruction") < STAGE_ORDER.index("final_scene_masks")
    assert STAGE_ORDER.index("final_scene_masks") < STAGE_ORDER.index("sky_reconstruction")
    assert STAGE_ORDER.index("sky_reconstruction") < STAGE_ORDER.index("galaxy_photometry")
    assert "sky_reconstruction" not in DEPENDENCIES["psf_reconstruction"]
    assert "psf_cutout_sky" in DEPENDENCIES["psf_reconstruction"]
    assert require_stage_dependencies(
        "psf_reconstruction", ("data_preparation", "psf_star_selection")
    ) == ("psf_cutout_sky",)


def test_manual_mask_edit_invalidates_later_products_but_not_psf():
    invalidated = invalidated_descendants(("manual_mask_draft",))
    assert invalidated == (
        "final_scene_masks", "sky_reconstruction",
        "refinement_assessment", "galaxy_photometry")
    assert "psf_reconstruction" not in invalidated
    assert "psf_star_selection" not in invalidated


def test_revised_galaxy_sky_does_not_force_rebuild_of_star_psf():
    invalidated = invalidated_descendants(("sky_reconstruction",))
    assert invalidated == ("refinement_assessment", "galaxy_photometry")
    assert "psf_reconstruction" not in invalidated
    assert "psf_cutout_sky" not in invalidated


def test_new_psf_star_sample_or_cutout_sky_invalidates_psf_dependents():
    for name in ("psf_star_selection", "psf_cutout_sky", "psf_reconstruction"):
        invalidated = invalidated_descendants((name,))
        assert "final_scene_masks" in invalidated
        assert "sky_reconstruction" in invalidated
        assert "galaxy_photometry" in invalidated
        if name != "psf_reconstruction":
            assert "psf_reconstruction" in invalidated


def test_changes_to_original_data_invalidate_every_dependent_artifact():
    invalidated = invalidated_descendants(("data_preparation",))
    assert "psf_reconstruction" in invalidated
    assert "galaxy_photometry" in invalidated
    assert "manual_mask_draft" not in invalidated
    assert invalidated_descendants(()) == ()
    with pytest.raises((ValueError, TypeError)):
        invalidated_descendants(("made_up_stage",))
    with pytest.raises(TypeError):
        invalidated_descendants("scene_initialization")


def test_global_psf_key_reused_by_same_dataset_across_projects():
    original = _psf()
    another_project_same_observation = _psf()
    assert original.digest == another_project_same_observation.digest
    assert original.digest != _psf(band="r").digest
    assert original.digest != _psf(instrument="decam").digest
    assert original.digest != _psf(source_sha256=_digest("e")).digest
    assert original.digest != _psf(image_kind="ccd").digest
    assert original.digest != _psf(cutout_sky_sha256=_digest("f")).digest
    assert original.digest != _psf(star_selection_sha256=_digest("f")).digest
    assert original.digest != _psf(model_scope="spatial").digest
    assert replace(original, observation_id="dr11/north/1554p180/g").digest != original.digest


def test_local_psf_tied_to_exact_evaluation_position():
    with pytest.raises(ValueError, match="needs its ICRS"):
        _psf(model_scope="local")
    local = _psf(model_scope="local", reference_icrs_deg=(155.438012, 18.088025))
    assert local.digest != _psf(
        model_scope="local", reference_icrs_deg=(155.438112, 18.088025)).digest
    with pytest.raises(ValueError, match="keyed by dataset"):
        _psf(model_scope="global", reference_icrs_deg=(155., 18.))


def test_psf_lineage_rejects_missing_or_invalid_provenance():
    with pytest.raises(ValueError, match="SHA-256"):
        _psf(source_sha256="unchecked")
    with pytest.raises(ValueError):
        _psf(observation_id="")
    with pytest.raises(ValueError):
        _psf(model_scope="unverified")
    with pytest.raises(ValueError):
        _psf(model_scope="local", reference_icrs_deg=(400., 18.))
    with pytest.raises(ValueError):
        _psf(band="W1")
