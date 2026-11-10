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
    assert result["target_protection_radius_arcsec"] == pytest.approx(20.)
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


def test_explicit_sga_core_preserves_candidate_shreds_outside_old_circle(
        monkeypatch, tmp_path):
    """SGA central size is protection evidence, NOT the LSB outer limit."""
    import lsbpipe.photometry.scene_contaminants as module
    from lsbpipe.photometry.sky_masks import SkyMaskDraft
    from lsbpipe.photometry.sga_priors import SGAGeometryPrior
    monkeypatch.setenv("LSBPIPE_SESSION_DIR",str(tmp_path/"session"))
    project=create_object_project(
        "Leo P",tmp_path/"LeoP_SGA",field=SkyField(155.,18.,4.))
    project.data_root=project.root/"data"
    project.data_root.mkdir()
    project.selected_dataset=("dr11","south","1554p180","i")
    project.target_radius_arcsec=16.
    save_object_project(project)
    image=np.zeros((241,241))
    dataset={"project_id":project.project_id,"band":"i",
             "source":"immutable-original","wcs":_wcs(),
             "pixel_shape":image.shape,
             "native_hard":np.zeros_like(image,bool),
             "data":{"image":image,
                     "maskbits":np.zeros_like(image,dtype=np.uint32)}}
    catalog=SimpleNamespace(
        source_path="native-tractor",source_sha256="a"*64,
        pixel_positions=np.array([[160.,120.],[209.,120.]]),
        rows=({**_row("PSF",flux=100.),"objid":31},
              {**_row("PSF",flux=100.),"objid":32}))
    monkeypatch.setattr(module,"read_native_tractor_gaia",
                        lambda *args,**kwargs:catalog)
    prior=SGAGeometryPrior(
        release="SGA2020",sga_id=71,name="Leo P",
        ra_deg=155.,dec_deg=18.,offset_arcsec=0.,
        d26_arcmin=2.,d26_ref="SB26",ba=.7,pa_deg=0.,
        ellipsebit=0,ellipsemode=0,flags=(),
        source_sha256="b"*64,source_path="SGA2020.fits")
    draft=SkyMaskDraft(project.project_id)
    plain=module.prepare_scene_contamination(
        dataset,project,draft,(155.,18.))
    assert plain["stats"]["masked_catalog_sources"]==2
    protected=module.prepare_scene_contamination(
        dataset,project,draft,(155.,18.),sga_prior=prior)
    assert protected["stats"]["masked_catalog_sources"]==1
    assert protected["stats"]["possible_sga_shreds"]==1
    assert protected["stats"]["skipped_target_overlap"]==1
    row=protected["stats"]["near_target_sources"][0]
    assert row["classification"]=="POSSIBLE_GALAXY_SHRED"
    assert row["objid"]==31
    assert not protected["mask"][120,160]
    assert protected["mask"][120,209]
    assert protected["sga_prior_id"]==71
    assert protected["sga_prior_sha256"]=="b"*64
    assert "SGA_CORE_ONLY_NOT_OUTER_LSB_BOUNDARY" in protected["flags"]
    assert np.array_equal(image,dataset["data"]["image"])
    assert not protected["mask"].flags.writeable


def test_sga_neighbour_with_gaia_counterpart_stays_unresolved_not_auto_removed():
    protected=np.zeros((241,241),bool)
    sga_core=np.zeros_like(protected)
    sga_core[100:140,100:140]=True
    protected |= sga_core
    rows=[{**_row("PSF",flux=30.), "objid":37,
           "gaia_g_mag":17.4,"ref_cat":"G2","ref_id":456}]
    mask,stats=conservative_catalogue_mask(
        (241,241),np.array([[125.,120.]]),rows,
        protected=protected,arcsec_per_pixel=1.,sga_core=sga_core)
    assert not mask.any()
    assert stats["possible_sga_shreds"]==0
    assert stats["near_target_sources"][0]["classification"]==(
        "GAIA_COUNTERPART_REVIEW")


def test_source_inventory_does_not_promote_gaia_psf_or_shreds_to_membership():
    from lsbpipe.photometry.scene_source_inventory import (
        build_scene_source_inventory, SOURCE_EVIDENCE_SCHEMA,
    )
    shape=(241,241)
    protected=np.zeros(shape,bool)
    protected[90:151,90:151]=True
    positions=np.array([[120.,120.],[180.,120.],[55.,50.],
                        [110.,130.]])
    records=(
        # Target-overlapping: Gaia match with significant, but potentially
        # unreliable astrometric evidence. DO NOT automatically mask.
        {**_row(),"objid":23,"ref_cat":"G3","ref_id":999,
         "parallax_mas":1.5,"parallax_ivar":36.,
         "pmra_masyr":6.,"pmra_ivar":4.},
        # Candidate mask outside target; not yet a physical foreground star.
        {**_row(),"objid":24,"ref_cat":"", "ref_id":0},
        # Retained Gaia source without ANY fitted band flux; old SC1
        # screening silently ignores it. Must remain in source inventory.
        {**_row(flux=25.),"objid":25,"ref_cat":"G3","ref_id":123456,
         "flux_ivar":0.,"gaia_g_mag":17.},
        # Forced PSF in local resolved-dwarf / SGA context: could be
        # galaxy light; neither MASKBITS nor FITBITS proves membership.
        {**_row(),"objid":26,"ref_cat":"L4","ref_id":77,
         "fitbits":1<<0,"maskbits":(1<<17)|(1<<12)},
    )
    proposed,screen=conservative_catalogue_mask(
        shape,positions,records,arcsec_per_pixel=1.,
        protected=protected)
    assert screen["skipped_target_overlap"]==2
    assert screen["masked_catalog_sources"]==1
    assert screen["unusable_catalog_sources"]==1
    assert not proposed[120,120] and proposed[120,180]
    native=SimpleNamespace(rows=records,pixel_positions=positions,
                           source_path="native-tractor.fits",
                           source_sha256="f"*64,
                           dataset_id="dr11/south/1554p180/i")
    inventory=build_scene_source_inventory(
        native,pixel_shape=shape,protected=protected,
        screen_stats=screen,release="dr11",region="south",
        brick="1554p180",band="i")
    assert inventory["schema"]==SOURCE_EVIDENCE_SCHEMA
    assert inventory["catalog_sha256"]=="f"*64
    assert inventory["total_in_native_frame"]==4
    by_id={row["objid"]:row for row in inventory["entries"]}
    assert set(by_id)=={23,24,25,26}
    assert by_id[23]["screening_disposition"]=="PROTECTED_OVERLAP"
    assert "GAIA_ASTROMETRIC_SIGNAL_FOREGROUND_REVIEW" in by_id[23][
        "evidence_flags"]
    assert by_id[23]["parallax_significance"]==pytest.approx(9.)
    assert by_id[23]["gaia_dr3_source_id"]==999
    assert by_id[24]["screening_disposition"]=="PROPOSED_EXCLUDE"
    assert "GEOMETRIC_MASK_ONLY_NO_TRACTOR_MODEL" in by_id[24][
        "evidence_flags"]
    assert by_id[25]["screening_disposition"]=="NOT_SCREENED_OR_UNFITTED"
    assert "RETAINED_GAIA_SOURCE_POSSIBLE" in by_id[25][
        "evidence_flags"]
    assert "NO_FITTED_TRACTOR_BAND_FLUX" in by_id[25]["evidence_flags"]
    assert by_id[26]["screening_disposition"]=="PROTECTED_OVERLAP"
    assert "PSF_TYPE_NOT_STELLAR_MEMBERSHIP_PROOF" in by_id[26][
        "evidence_flags"]
    assert "RESOLVED_DWARF_CONTEXT_NOT_MEMBERSHIP_PROOF" in by_id[26][
        "evidence_flags"]
    assert by_id[26]["reference_catalog"]=="L4"
    assert all(row["physical_membership"]=="UNDETERMINED" for row in by_id.values())
    assert all(row["photometry_action"]=="NO_NEW_ACTION_READ_ONLY"
               for row in by_id.values())
    assert not proposed[120,120]
    assert inventory["truncated"]==0


def test_source_inventory_explicit_truncation_is_not_silent():
    from lsbpipe.photometry.scene_source_inventory import (
        build_scene_source_inventory,
    )
    points=np.array([[30.,30.],[31.,30.],[32.,30.]])
    rows=tuple({**_row(),"objid":i} for i in range(3))
    catalog=SimpleNamespace(rows=rows,pixel_positions=points,
                            source_path="local",source_sha256="a"*64)
    result=build_scene_source_inventory(
        catalog,pixel_shape=(61,61),
        protected=np.zeros((61,61),bool),screen_stats={},
        release="dr11",region="south",brick="1554p180",
        band="g",max_records=2)
    assert result["total_in_native_frame"]==3
    assert result["stored_count"]==2
    assert result["truncated"]==1
    assert "INVENTORY_TRUNCATED" in result["cautions"]
