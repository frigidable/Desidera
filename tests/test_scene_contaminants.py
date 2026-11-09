"""Scientific safety checks for source-bound provisional Scene masking."""
from __future__ import annotations

from types import SimpleNamespace
import numpy as np
import pytest
from astropy.wcs import WCS

from lsbpipe.photometry.scene_contaminants import (
    conservative_catalogue_mask, prepare_scene_contamination,
)
from lsbpipe.surveys.legacy import SkyField
from lsbpipe.surveys.object_projects import create_object_project, save_object_project


def _wcs():
    w = WCS(naxis=2)
    w.wcs.crpix = [121.,121.]
    w.wcs.crval = [155.,18.]
    w.wcs.cd = np.diag([-1.,1.])/3600.
    w.wcs.ctype = ["RA---TAN","DEC--TAN"]
    return w


def _row(type_="PSF", flux=100., shape_r=None):
    return {"type":type_, "flux":flux, "flux_ivar":5.,
            "shape_r_arcsec":shape_r}


def test_catalogue_mask_excludes_outer_stars_but_not_target_or_bad_data():
    protected = np.zeros((241,241),bool)
    yy, xx = np.indices(protected.shape)
    protected[(xx-120.)**2+(yy-120.)**2 < 25.**2] = True
    positions = np.array([[120., 120.], [181.,112.], [55.,73.],
                          [500., 400.]])
    rows = (_row(flux=150.),_row(flux=1000.),
            _row("EXP",flux=80.,shape_r=6.),_row())
    mask, counts = conservative_catalogue_mask(
        protected.shape, positions, rows, arcsec_per_pixel=1.,
        protected=protected)
    assert counts["masked_catalog_sources"] == 2
    assert counts["skipped_target_overlap"] == 1
    assert counts["outside_cutout"] == 1
    assert mask[112,181] and mask[73,55]
    assert not mask[120,120]
    assert not np.any(mask & protected)
    assert not mask.flags.writeable


def test_full_preparation_uses_exact_catalogue_identity_and_preserves_originals(
        monkeypatch, tmp_path):
    import lsbpipe.photometry.scene_contaminants as module
    monkeypatch.setenv("LSBPIPE_SESSION_DIR",str(tmp_path/"session"))
    project = create_object_project(
        "Leo P",tmp_path/"LeoP",field=SkyField(155.,18.,4.))
    project.data_root = project.root/"data"
    project.data_root.mkdir()
    project.selected_dataset = ("dr11","south","1554p180","g")
    project.target_radius_arcsec = 20.
    save_object_project(project)
    yy,xx = np.indices((241,241))
    original = (np.sin(xx/11.)+.5*np.cos(yy/8.))
    image = original.copy()
    hard = np.zeros(image.shape,bool)
    hard[111,181] = True
    bits = np.zeros(image.shape,np.uint32)
    bits[120,120] = 1 << 12  # GALAXY may flag target itself.
    bits[110,180] = 1 << 1   # BRIGHT warning is not unconditional exclusion.
    dataset = {
        "project_id":project.project_id,
        "source":"DR11 south / 1554p180 / g",
        "band":"g", "wcs":_wcs(),
        "pixel_shape":image.shape,
        "native_hard":hard, "data":{"image":image, "maskbits":bits},
    }
    from lsbpipe.photometry.sky_masks import SkyMaskDraft
    draft = SkyMaskDraft(project.project_id)
    catalogue = SimpleNamespace(
        source_path=str(project.data_root/"tractor-original.fits"),
        source_sha256="f"*64,
        pixel_positions=np.array([[120.,120.],[181.,112.]]),
        rows=(_row(),_row(flux=1000.)))
    calls = []
    def fake_read(root, **kwargs):
        calls.append((root,kwargs["release"],kwargs["region"],
                      kwargs["brick"],kwargs["frame"].dataset_id))
        return catalogue
    monkeypatch.setattr(module,"read_native_tractor_gaia",fake_read)
    result=prepare_scene_contamination(dataset, project, draft,
                                       (155.,18.))
    assert calls == [(project.data_root,"dr11","south","1554p180",
                      "dr11/south/1554p180/g")]
    assert result["stats"]["masked_catalog_sources"] == 1
    assert result["stats"]["skipped_target_overlap"] == 1
    assert result["mask"][112,181] and not result["mask"][120,120]
    assert result["maskbits_caution_counts"]["GALAXY"] == 1
    assert result["maskbits_caution_counts"]["BRIGHT"] == 1
    assert result["draft_digest"] == draft.digest
    assert result["catalog_sha256"] == "f"*64
    assert not result["mask"].flags.writeable
    assert np.array_equal(image,original)
    assert np.array_equal(dataset["data"]["maskbits"],bits)


def test_missing_catalogue_and_unsupported_release_fail_closed(
        monkeypatch,tmp_path):
    import lsbpipe.photometry.scene_contaminants as module
    monkeypatch.setenv("LSBPIPE_SESSION_DIR",str(tmp_path/"session"))
    project=create_object_project(
        "Unknown LSB",tmp_path/"unknown",field=SkyField(155.,18.,4.))
    project.data_root=project.root/"data"
    project.data_root.mkdir()
    project.selected_dataset=("dr11","south","1554p180","r")
    save_object_project(project)
    pixels=np.zeros((241,241))
    dataset={"project_id":project.project_id,"band":"r",
             "source":"original","wcs":_wcs(),
             "pixel_shape":pixels.shape,
             "native_hard":np.zeros_like(pixels,bool),
             "data":{"image":pixels,"maskbits":np.zeros_like(pixels,np.uint32)}}
    from lsbpipe.photometry.sky_masks import SkyMaskDraft
    with pytest.raises(FileNotFoundError):
        prepare_scene_contamination(
            dataset,project,SkyMaskDraft(project.project_id),(155.,18.))
    project.selected_dataset=("dr9","north","1554p180","r")
    with pytest.raises(ValueError,match="DR11 South only"):
        prepare_scene_contamination(
            dataset,project,SkyMaskDraft(project.project_id),(155.,18.))


def test_tractor_ellipticity_creates_native_wcs_ellipses_not_round_disks():
    """DR11 Tractor epsilon shape is a sky-plane PA, not raw pixel angle."""
    ny=nx=241
    protected=np.zeros((ny,nx),bool)
    catalog=[{"type":"EXP","flux":100., "flux_ivar":5.,
              "shape_r_arcsec":6., "shape_e1":.5, "shape_e2":0.,
              "objid":23}]
    mask,stats=conservative_catalogue_mask(
        (ny,nx),np.array([[120.,120.]]),catalog,
        arcsec_per_pixel=1.,protected=protected,wcs=_wcs())
    assert stats["elliptical_catalog_sources"]==1
    assert mask[120,120]
    assert mask[137,120]  # phi=0 -> PA=0, major axis along north
    assert not mask[120,137]  # extended source is no longer a round mask
    assert stats["masked_catalog_sources"]==1


def test_tractor_missing_or_unphysical_shape_falls_back_without_guessing_pa():
    protected=np.zeros((241,241),bool)
    rows=[{**_row("SER",flux=100.,shape_r=6.),
           "shape_e1":None,"shape_e2":None}]
    mask,stats=conservative_catalogue_mask(
        (241,241),np.array([[110.,110.]]),rows,
        arcsec_per_pixel=1.,protected=protected,wcs=_wcs())
    assert stats["elliptical_catalog_sources"]==0
    assert stats["missing_shape_fallbacks"]==1
    assert mask[110,127] and mask[127,110]


def test_catalogue_target_overlap_records_source_for_manual_review():
    protected=np.zeros((241,241),bool)
    protected[110:131,110:131]=True
    rows=[{**_row(), "objid":418}]
    mask,stats=conservative_catalogue_mask(
        (241,241),np.array([[136.,120.]]),rows,
        arcsec_per_pixel=1.,protected=protected)
    assert stats["skipped_target_overlap"]==1
    assert not mask.any()
    assert stats["near_target_sources"][0]["objid"]==418
