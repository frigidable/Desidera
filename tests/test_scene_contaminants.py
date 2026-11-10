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


def test_zero_gaia_magnitude_placeholders_do_not_create_fake_gaia_matches():
    """DR11 non-Gaia rows commonly have gaia_g_mag=0.0.

    A null/zero magnitude must not turn every protected blend and every
    masked-centre record into a "Gaia counterpart". Preserve real G2
    and G3 catalogue references without claiming stellar membership.
    """
    protected = np.zeros((241, 241), dtype=bool)
    protected[100:140, 100:140] = True
    rows = (
        {**_row(), "objid": 11, "ref_cat": "", "ref_id": 0,
         "gaia_g_mag": 0.0},
        {**_row(), "objid": 12, "ref_cat": "G2", "ref_id": 789,
         "gaia_g_mag": 17.4},
        {**_row(), "objid": 13, "ref_cat": "", "ref_id": 0,
         "gaia_g_mag": 0.0},
        {**_row(), "objid": 14, "ref_cat": "G3", "ref_id": 456,
         "gaia_g_mag": 18.0},
    )
    positions = np.array([[116., 120.], [125., 125.],
                          [180., 120.], [210., 120.]])
    mask, stats = conservative_catalogue_mask(
        protected.shape, positions, rows,
        arcsec_per_pixel=1., protected=protected)
    assert stats["skipped_target_overlap"] == 2
    assert stats["masked_catalog_sources"] == 2
    overlap = {v["objid"]:v for v in stats["near_target_sources"]}
    assert overlap[11]["classification"] == "UNRESOLVED_TARGET_BLEND"
    assert overlap[12]["classification"] == "GAIA_COUNTERPART_REVIEW"
    centres = {v["objid"]:v for v in stats["masked_catalog_centres"]}
    assert centres[13]["gaia_counterpart"] is False
    assert centres[14]["gaia_counterpart"] is True
    assert not mask[120,116] and not mask[125,125]
    assert mask[120,180] and mask[120,210]


def test_gaia_reference_requires_valid_catalogue_identifier():
    """A Gaia magnitude alone is insufficient to claim matched identity."""
    protected = np.zeros((241, 241), dtype=bool)
    protected[105:135, 105:135] = True
    rows = (
        {**_row(), "objid": 21, "ref_cat": "G3", "ref_id": 0,
         "gaia_g_mag": 17.0},
        {**_row(), "objid": 22, "ref_cat": "G3", "ref_id": None,
         "gaia_g_mag": 16.0},
        {**_row(), "objid": 23, "ref_cat": "", "ref_id": 99,
         "gaia_g_mag": 17.0},
    )
    positions = np.array([[112., 113.], [121., 113.], [130., 115.]])
    _, stats = conservative_catalogue_mask(
        protected.shape, positions, rows,
        arcsec_per_pixel=1., protected=protected)
    assert stats["skipped_target_overlap"] == 3
    assert all(v["classification"] == "UNRESOLVED_TARGET_BLEND"
               for v in stats["near_target_sources"])


def test_source_registry_decodes_all_dr11_bit_names():
    """DR11 flags 8..13 must not disappear from the evidence browser.

    In Leo P, FITBITS ITERATIVE, GAIA_POINTSOURCE, WALKER and RUNNER
    occur even though the old registry only decoded FITBITS 0..7.
    These tags describe fitting history, never a stellar/galaxy truth.
    """
    from lsbpipe.photometry.scene_source_inventory import (
        build_scene_source_inventory, FITBIT_NAMES, MASKBIT_NAMES,
    )
    for bit, name in {
        8: "TYCHO2", 9: "LARGEGALAXY", 10: "WALKER",
        11: "RUNNER", 12: "GAIA_POINTSOURCE", 13: "ITERATIVE",
    }.items():
        assert FITBIT_NAMES[bit] == name
    for bit, name in {
        0: "NPRIMARY", 8: "WISEM1", 9: "WISEM2",
        19: "WISE_GAIA",
    }.items():
        assert MASKBIT_NAMES[bit] == name
    rows = (
        {**_row("SER"), "objid": 2059,
         "ref_cat": "L4", "ref_id": 5053883,
         "fitbits": (1<<9) | (1<<4) | (1<<1),
         "maskbits": (1<<12)},
        {**_row("SER"), "objid": 1762,
         "fitbits": (1<<2) | (1<<10) | (1<<11) | (1<<13),
         "maskbits": (1<<0) | (1<<19)},
        {**_row("PSF"), "objid": 2035,
         "ref_cat": "G3", "ref_id": 623919492384922624,
         "fitbits": (1<<7) | (1<<12)},
    )
    points = np.array([[12., 12.], [20., 20.], [28., 28.]])
    catalogue = SimpleNamespace(
        rows=rows, pixel_positions=points,
        source_path="original-tractor.fits", source_sha256="a"*64)
    result = build_scene_source_inventory(
        catalogue, pixel_shape=(41, 41),
        protected=np.zeros((41, 41), dtype=bool),
        screen_stats={}, release="dr11", region="south",
        brick="1554p180", band="g")
    by_id = {entry["objid"]:entry for entry in result["entries"]}
    assert {"LARGEGALAXY", "FROZEN", "FIT_BACKGROUND"}.issubset(
        set(by_id[2059]["fitbits_names"]))
    assert {"HIT_RADIUS_LIMIT", "WALKER", "RUNNER", "ITERATIVE"}.issubset(
        set(by_id[1762]["fitbits_names"]))
    assert {"NPRIMARY", "WISE_GAIA"}.issubset(
        set(by_id[1762]["maskbits_names"]))
    assert "GAIA_POINTSOURCE" in by_id[2035]["fitbits_names"]
    assert all(v["physical_membership"] == "UNDETERMINED"
               for v in by_id.values())


def test_off_frame_mask_contributors_keep_original_source_identity():
    """Sources outside the cutout can mask *in-frame* pixels.

    The old registry omitted these accepted screen sources, leaving
    excluded pixels with no inspectable Tractor identity. Preserve the
    2D native-centre count and report contributing off-frame centres.
    """
    from lsbpipe.photometry.scene_source_inventory import (
        build_scene_source_inventory,
    )
    protected = np.zeros((61, 61), dtype=bool)
    protected[22:39, 22:39] = True
    rows = (
        {**_row(flux=100.), "objid": 1},  # protected source
        {**_row(flux=100.), "objid": 2},  # source just outside right edge
        {**_row(flux=100.), "objid": 3},  # ordinary in-frame exclusion
        {**_row(flux=100.), "objid": 4},  # too distant to affect cutout
    )
    points = np.array([[30., 30.], [62., 30.],
                       [10., 10.], [100., 40.]])
    mask, stats = conservative_catalogue_mask(
        protected.shape, points, rows,
        arcsec_per_pixel=1., protected=protected)
    assert stats["masked_catalog_sources"] == 2
    assert stats["skipped_target_overlap"] == 1
    assert stats["outside_cutout"] == 1
    assert mask[30, 60] and mask[10, 10]
    catalog = SimpleNamespace(
        rows=rows, pixel_positions=points,
        source_path="local-original-tractor.fits",
        source_sha256="d"*64)
    registry = build_scene_source_inventory(
        catalog, pixel_shape=protected.shape,
        protected=protected, screen_stats=stats,
        release="dr11", region="south", brick="1554p180", band="g")
    assert registry["total_in_native_frame"] == 2
    assert registry["contributing_off_frame"] == 1
    assert registry["total_relevant_sources"] == 3
    assert registry["stored_count"] == 3
    assert registry["screening_counts_all_native"] == {
        "PROTECTED_OVERLAP": 1, "PROPOSED_EXCLUDE": 1,
        "NOT_SCREENED_OR_UNFITTED": 0,
    }
    assert registry["screening_counts_all_relevant"] == {
        "PROTECTED_OVERLAP": 1, "PROPOSED_EXCLUDE": 2,
        "NOT_SCREENED_OR_UNFITTED": 0,
    }
    assert "OFF_FRAME_CONTRIBUTING_MASK_SOURCES" in registry["cautions"]
    by_id = {row["objid"]:row for row in registry["entries"]}
    assert set(by_id) == {1, 2, 3}
    assert by_id[2]["centre_in_native_frame"] is False
    assert by_id[2]["screening_disposition"] == "PROPOSED_EXCLUDE"
    assert all(row["physical_membership"] == "UNDETERMINED"
               for row in by_id.values())



def test_gaia_astrometry_review_is_not_automatic_mask_and_legacy_iso_is_inactive(
        tmp_path, monkeypatch):
    """A saved exploratory outline cannot suppress neighbour screening.

    Gaia foreground *evidence* is shown, but no stellar PSF is silently
    removed from the protected target. A separate deliberate manual
    exclusion still works, including within target protection.
    """
    import lsbpipe.photometry.scene_contaminants as module
    from lsbpipe.photometry.sky_masks import (
        SkyMaskDraft, SkyMaskEdit, SkyMaskShape,
    )
    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path/"session"))
    project = create_object_project(
        "Leo P", tmp_path/"LeoP_Legacy",
        field=SkyField(155., 18., 4.))
    project.data_root = project.root/"data"
    project.data_root.mkdir()
    project.selected_dataset = ("dr11", "south", "1554p180", "g")
    project.target_radius_arcsec = 20.
    save_object_project(project)
    image = np.zeros((241, 241), dtype=float)
    dataset = {
        "project_id": project.project_id, "band": "g",
        "source": "native-test", "wcs": _wcs(),
        "pixel_shape": image.shape, "native_hard": np.zeros_like(image, bool),
        "data": {"image": image,
                 "maskbits": np.zeros_like(image, dtype=np.uint32)},
    }
    rows = (
        {**_row(flux=600.), "objid": 2035,
         "ref_cat": "G3", "ref_id": 623919492384922624,
         "parallax_mas": 1.1, "parallax_ivar": 900.,
         "pmra_masyr": -1.5, "pmra_ivar": 400.},
        {**_row(flux=100.), "objid": 1674,
         "ref_cat": "", "ref_id": 0},
    )
    catalogue = SimpleNamespace(
        source_path=str(project.data_root/"tractor.fits"),
        source_sha256="a"*64,
        pixel_positions=np.array([[120., 120.], [180., 120.]]),
        rows=rows)
    monkeypatch.setattr(module, "read_native_tractor_gaia",
                        lambda *args, **kwargs: catalogue)
    empty = SkyMaskDraft(project.project_id)
    legacy = SkyMaskEdit(
        "sc1_iso_older", "protect_target",
        SkyMaskShape("circle", (155., 18.), radius_arcsec=90.),
        comment="AUTO_ISOPHOTE_PRELIMINARY kind=contour")
    assert legacy.applies_to("g")
    previously_saved = empty.with_edit(legacy)
    baseline = module.prepare_scene_contamination(
        dataset, project, empty, (155., 18.))
    recovered = module.prepare_scene_contamination(
        dataset, project, previously_saved, (155., 18.))
    assert baseline["stats"]["masked_catalog_sources"] == 1
    assert recovered["stats"]["masked_catalog_sources"] == 1
    assert baseline["stats"]["skipped_target_overlap"] == 1
    assert np.array_equal(recovered["mask"], baseline["mask"])
    assert np.array_equal(recovered["protected"], baseline["protected"])
    new_proposal = SkyMaskEdit(
        "sc1_iso_new_proposal", "isophote_candidate",
        SkyMaskShape("circle", (155.,18.), radius_arcsec=90.),
        band="g")
    annotated = module.prepare_scene_contamination(
        dataset, project, empty.with_edit(new_proposal), (155.,18.))
    assert np.array_equal(annotated["mask"], baseline["mask"])
    assert np.array_equal(annotated["protected"], baseline["protected"])
    assert annotated["stats"]["masked_catalog_sources"] == 1
    assert annotated["stats"]["skipped_target_overlap"] == 1
    assert annotated["stats"]["legacy_isophote_protect_edits_ignored"] == 0
    assert recovered["stats"]["legacy_isophote_protect_edits_ignored"] == 1
    assert "LEGACY_ISOPHOTE_PROTECT_EDITS_INACTIVE" in recovered["flags"]
    assert "PROTECTED_GAIA_ASTROMETRY_STARS_NOT_REMOVED" in recovered["flags"]
    assert recovered["stats"]["protected_gaia_astrometry_review"] == 1
    assert recovered["foreground_astrometry_protected_review_ids"] == (
        "dr11/south/1554p180/OBJID=2035",)
    by_id = {e["objid"]: e for e in recovered["source_inventory"]["entries"]}
    assert by_id[2035]["foreground_evidence_assessment"] == (
        "GAIA_ASTROMETRY_FOREGROUND_CANDIDATE_REVIEW")
    assert by_id[2035]["proposed_source_treatment"] == (
        "REVIEW_FOREGROUND_PSF_OR_MANUAL_EXCLUSION")
    assert by_id[2035]["proposed_treatment_is_active"] is False
    assert by_id[2035]["physical_membership"] == "UNDETERMINED"
    assert by_id[2035]["photometry_action"] == "NO_NEW_ACTION_READ_ONLY"
    assert not recovered["mask"][120, 120]
    assert recovered["mask"][120, 180]
    manual = SkyMaskEdit(
        "user_exclude_g3", "exclude",
        SkyMaskShape("circle", (155.,18.), radius_arcsec=3.),
        band="g", comment="Independent foreground-star QA")
    manual_result = module.prepare_scene_contamination(
        dataset, project, previously_saved.with_edit(manual), (155.,18.))
    assert manual_result["mask"][120,120]
    assert manual_result["mask"][120,180]
    assert manual_result["stats"]["protected_gaia_astrometry_review"] == 1
    assert np.array_equal(image, dataset["data"]["image"])



def test_negative_gaia_parallax_alert_is_not_foreground_evidence():
    """Large negative parallax can indicate an unreliable Gaia solution.

    Even nominally significant component proper motion cannot make a
    strongly negative astrometric solution a validated foreground star.
    The QA evidence remains read-only in both cases.
    """
    from lsbpipe.photometry.scene_source_inventory import (
        build_scene_source_inventory,
    )
    records = (
        {**_row(), "objid": 901, "ref_cat": "G3", "ref_id": 123,
         "parallax_mas": -4., "parallax_ivar": 25.,
         "pmra_masyr": 10., "pmra_ivar": 25.},
        {**_row(), "objid": 902, "ref_cat": "G3", "ref_id": 124,
         "parallax_mas": 1., "parallax_ivar": 36.,
         "pmra_masyr": 0., "pmra_ivar": 1.},
    )
    catalog = SimpleNamespace(
        rows=records, pixel_positions=np.array([[12.,12.],[24.,24.]]),
        source_path="original-tractor.fits", source_sha256="b"*64)
    protected = np.ones((41,41), dtype=bool)
    e = build_scene_source_inventory(
        catalog, pixel_shape=protected.shape, protected=protected,
        screen_stats={}, release="dr11", region="south",
        brick="1554p180", band="g")["entries"]
    bad, good = e
    assert bad["objid"] == 901 and good["objid"] == 902
    assert bad["parallax_significance"] == pytest.approx(20.)
    assert bad["foreground_evidence_assessment"] == (
        "GAIA_ASTROMETRY_QUALITY_ALERT")
    assert "GAIA_NEGATIVE_PARALLAX_QUALITY_REVIEW" in bad["evidence_flags"]
    assert "GAIA_ASTROMETRIC_SIGNAL_FOREGROUND_REVIEW" not in bad["evidence_flags"]
    assert bad["proposed_source_treatment"] == "VERIFY_GAIA_ASTROMETRY_QUALITY"
    assert good["foreground_evidence_assessment"] == (
        "GAIA_ASTROMETRY_FOREGROUND_CANDIDATE_REVIEW")
    assert "GAIA_ASTROMETRIC_SIGNAL_FOREGROUND_REVIEW" in good["evidence_flags"]
    assert all(v["physical_membership"] == "UNDETERMINED" for v in e)
    assert all(v["photometry_action"] == "NO_NEW_ACTION_READ_ONLY" for v in e)
    assert all(v["proposed_treatment_is_active"] is False for v in e)



def test_internal_gaia_mask_preview_does_not_touch_non_gaia_or_negative_parallax():
    """Only opt-in candidate PSF footprints enter target-protected pixels.

    A galaxy knot and a doubtful negative-parallax Gaia source remain
    unresolved; brightness/type alone are never membership decisions.
    """
    protected = np.zeros((241, 241), bool)
    yy, xx = np.indices(protected.shape)
    protected[(xx-120.)**2 + (yy-120.)**2 < 40.**2] = True
    positions = np.array([[120.,120.], [160.,120.],
                          [120.,160.], [214.,120.]])
    rows = (
        {**_row("PSF", flux=600.), "objid": 2035,
         "ref_cat": "G3", "ref_id": 12345,
         "parallax_mas": 1.1, "parallax_ivar": 900.,
         "pmra_masyr": 1., "pmra_ivar": 900.},
        {**_row("REX", flux=100., shape_r=2.), "objid": 2076,
         "ref_cat": "L4", "ref_id": 12},
        {**_row("PSF", flux=100.), "objid": 901,
         "ref_cat": "G3", "ref_id": 987,
         "parallax_mas": -4., "parallax_ivar": 25.,
         "pmra_masyr": 10., "pmra_ivar": 25.},
        {**_row("PSF", flux=100.), "objid": 1000},
    )
    baseline, old = conservative_catalogue_mask(
        protected.shape, positions, rows, arcsec_per_pixel=1.,
        protected=protected)
    preview, new = conservative_catalogue_mask(
        protected.shape, positions, rows, arcsec_per_pixel=1.,
        protected=protected, gaia_inner_preview=True)
    assert old["masked_catalog_sources"] == 1
    assert old["skipped_target_overlap"] == 3
    assert new["masked_catalog_sources"] == 2
    assert new["skipped_target_overlap"] == 2
    assert new["gaia_inner_preview_sources"] == 1
    assert [e["objid"] for e in new["gaia_inner_preview_centres"]] == [2035]
    assert not baseline[120,120] and preview[120,120]
    assert not baseline[120,160] and not preview[120,160]
    assert not baseline[160,120] and not preview[160,120]
    assert baseline[120,214] and preview[120,214]
    assert np.array_equal(
        preview, baseline | new["gaia_inner_preview_mask"])
    assert not new["gaia_inner_preview_mask"][120,214]
    assert not preview.flags.writeable
    assert not new["gaia_inner_preview_mask"].flags.writeable


def test_opt_in_internal_gaia_mask_does_not_change_source_membership(
        monkeypatch, tmp_path):
    """Prepare returns a separate auditable layer, not an adopted science mask."""
    import lsbpipe.photometry.scene_contaminants as module
    from lsbpipe.photometry.sky_masks import SkyMaskDraft
    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path/"session"))
    project = create_object_project(
        "Leo P", tmp_path/"Leo_P_preview",
        field=SkyField(155.,18.,4.))
    project.data_root = project.root/"data"
    project.data_root.mkdir()
    project.selected_dataset = ("dr11","south","1554p180","g")
    project.target_radius_arcsec = 20.
    save_object_project(project)
    image = np.ones((241,241))
    ds = {"project_id":project.project_id, "band":"g",
          "source":"native-test", "wcs":_wcs(),
          "pixel_shape":image.shape,
          "native_hard":np.zeros(image.shape,bool),
          "data":{"image":image, "maskbits":np.zeros(image.shape,np.uint32)}}
    rows = (
        {**_row(flux=600.), "objid":2035,"ref_cat":"G3","ref_id":12345,
         "parallax_mas":1.1,"parallax_ivar":900.},
        {**_row(flux=100.),"objid":900,"ref_cat":"","ref_id":0},
    )
    catalog = SimpleNamespace(
        rows=rows,pixel_positions=np.array([[120.,120.],[205.,120.]]),
        source_path="native-tractor.fits",source_sha256="b"*64)
    monkeypatch.setattr(module,"read_native_tractor_gaia",
                        lambda *args,**kwargs:catalog)
    draft=SkyMaskDraft(project.project_id)
    base=module.prepare_scene_contamination(
        ds,project,draft,(155.,18.))
    new=module.prepare_scene_contamination(
        ds,project,draft,(155.,18.),gaia_inner_preview=True)
    assert not base["gaia_inner_preview_enabled"]
    assert new["gaia_inner_preview_enabled"]
    assert base["stats"]["masked_catalog_sources"]==1
    assert new["stats"]["masked_catalog_sources"]==2
    assert base["stats"]["skipped_target_overlap"]==1
    assert new["stats"]["skipped_target_overlap"]==0
    assert new["stats"]["gaia_inner_preview_sources"]==1
    assert new["stats"]["gaia_inner_preview_pixels"]>0
    assert np.array_equal(new["mask"], base["mask"] |
                          new["gaia_inner_preview_mask"])
    assert not base["mask"][120,120] and new["mask"][120,120]
    assert new["catalog_sha256"]==base["catalog_sha256"]=="b"*64
    assert new["draft_digest"]==base["draft_digest"]==draft.digest
    assert "EXPLICIT_GAIA_INTERIOR_MASK_PREVIEW_ONLY" in new["flags"]
    assert "GAIA_WING_AND_GALAXY_FLUX_LOSS_NOT_QUANTIFIED" in new["flags"]
    e = {x["objid"]:x for x in new["source_inventory"]["entries"]}
    assert e[2035]["screening_disposition"]=="GAIA_FOREGROUND_MASK_PREVIEW"
    assert e[2035]["diagnostic_inner_gaia_mask_preview_active"] is True
    assert e[2035]["photometry_action"]=="NO_NEW_ACTION_READ_ONLY"
    assert e[2035]["physical_membership"]=="UNDETERMINED"
    assert e[900]["screening_disposition"]=="PROPOSED_EXCLUDE"
    assert new["source_inventory"]["screening_counts_all_relevant"][
        "GAIA_FOREGROUND_MASK_PREVIEW"]==1
    assert not new["mask"].flags.writeable
    assert not new["gaia_inner_preview_mask"].flags.writeable
    assert np.array_equal(image,ds["data"]["image"])
    assert draft.edits==()



def _red_extended_row(objid, *, source_type="EXP", red=True,
                      blue=False, shape_r=1.1, fitbits=0, maskbits=0,
                      ref_cat="", ref_id=0):
    """Synthetic multi-band Tractor source; NOT a real observation."""
    g, r, z = (1., 3., 6.) if red else (1., 1.1, 1.25)
    if blue:
        g, r, z = 4., 1.5, 1.
    return {
        **_row(source_type, flux=g, shape_r=shape_r),
        "objid":objid, "flux_ivar":100.,
        "ref_cat":ref_cat, "ref_id":ref_id,
        "fitbits":fitbits, "maskbits":maskbits,
        "model_fluxes_nMgy":{"g":g,"r":r,"z":z},
        "model_flux_ivars":{"g":400.,"r":400.,"z":400.},
    }


def test_red_extended_v1_multicue_policy_vetoes_member_like_sources():
    from lsbpipe.photometry.scene_source_inventory import (
        conservative_red_extended_candidate,
    )
    extra={"offset_arcsec":32.,"neighbour_arcsec":5.,
           "target_radius_arcsec":60.}
    candidate=_red_extended_row(100)
    passes,evidence=conservative_red_extended_candidate(candidate,**extra)
    assert passes and evidence["policy"]=="red_extended_v1"
    assert evidence["passes"] is True
    assert evidence["g_minus_r"] >= .8
    # Genuine dwarf carbon stars may be very red, but PSF-only
    # photometry and colors are NOT proof of background membership.
    assert not conservative_red_extended_candidate(
        _red_extended_row(101,source_type="PSF"),**extra)[0]
    # Bluer star-forming knots and spurious SGA Tractor shreds stay protected.
    assert not conservative_red_extended_candidate(
        _red_extended_row(102,blue=True),**extra)[0]
    assert not conservative_red_extended_candidate(
        _red_extended_row(103,ref_cat="L4",ref_id=15),**extra)[0]
    assert not conservative_red_extended_candidate(
        _red_extended_row(104,maskbits=1<<12),**extra)[0]
    assert not conservative_red_extended_candidate(
        _red_extended_row(105,fitbits=18),**extra)[0]
    assert not conservative_red_extended_candidate(
        _red_extended_row(106,shape_r=.4),**extra)[0]
    assert not conservative_red_extended_candidate(
        _red_extended_row(107),**{**extra,"offset_arcsec":4.})[0]
    assert not conservative_red_extended_candidate(
        _red_extended_row(108),**{**extra,"neighbour_arcsec":1.})[0]
    assert not conservative_red_extended_candidate(
        {**candidate,"model_flux_ivars":{"g":0.,"r":400.,"z":400.}},
        **extra)[0]
    assert not conservative_red_extended_candidate(
        {**candidate,"model_fluxes_nMgy":{"g":1.,"r":3.,"z":None}},
        **extra)[0]
    assert not conservative_red_extended_candidate(
        {**candidate,"gaia_duplicated_source":True},**extra)[0]
    assert not conservative_red_extended_candidate(
        {**candidate,"maskbits":None},**extra)[0]
    assert not conservative_red_extended_candidate(
        candidate,**{**extra,"neighbour_arcsec":float("nan")})[0]


def test_red_extended_preview_only_masks_isolated_red_galaxies_not_blue_knots():
    """An all-catalogue non-Gaia procedure independent of external papers."""
    center=(120.,120.)
    protected=np.zeros((241,241),bool)
    yy,xx=np.indices(protected.shape)
    protected[(xx-120.)**2+(yy-120.)**2<50.**2]=True
    rows=(
        _red_extended_row(1794),  # red+extended+isolated
        _red_extended_row(1946,blue=True), # blue galaxy star cluster
        _red_extended_row(1939,source_type="PSF"), # red target carbon star
        _red_extended_row(2059,ref_cat="L4",ref_id=1), # target component
        _red_extended_row(2048,fitbits=18), # source forced by GALAXY
        _red_extended_row(1888,source_type="PSF",red=False),
    )
    pos=np.array([[160.,120.],[130.,125.],[115.,135.],
                  [115.,110.],[135.,115.],[110.,150.]])
    baseline,legacy=conservative_catalogue_mask(
        protected.shape,pos,rows,arcsec_per_pixel=1.,
        protected=protected,wcs=_wcs())
    preview,stats=conservative_catalogue_mask(
        protected.shape,pos,rows,arcsec_per_pixel=1.,
        protected=protected,wcs=_wcs(),
        red_extended_preview=True,anchor_pixel=center,
        target_radius_arcsec=50.)
    assert legacy["masked_catalog_sources"]==0
    assert stats["red_extended_preview_sources"]==1
    assert stats["masked_catalog_sources"]==1
    assert stats["skipped_target_overlap"]==5
    assert [v["objid"] for v in stats["red_extended_preview_centres"]]==[1794]
    assert stats["red_extended_preview_centres"][0]["evidence"]["policy"]==(
        "red_extended_v1")
    assert preview[120,160] and not baseline[120,160]
    assert not preview[125,130]
    assert not preview[135,115]
    assert np.array_equal(preview, stats["red_extended_preview_mask"])
    assert not preview.flags.writeable
    assert not stats["red_extended_preview_mask"].flags.writeable
    # Cannot enable isolated color classifier without target geometry.
    with pytest.raises(ValueError,match="current target geometry"):
        conservative_catalogue_mask(
            protected.shape,pos,rows,arcsec_per_pixel=1.,
            protected=protected,red_extended_preview=True)


def test_red_extended_mask_preview_keeps_origins_and_source_roles_separate(
        monkeypatch,tmp_path):
    import lsbpipe.photometry.scene_contaminants as module
    from lsbpipe.photometry.sky_masks import SkyMaskDraft
    monkeypatch.setenv("LSBPIPE_SESSION_DIR",str(tmp_path/"session"))
    project=create_object_project(
        "Leo P",tmp_path/"LeoP_RedPreview",field=SkyField(155.,18.,4.))
    project.data_root=project.root/"data"
    project.data_root.mkdir()
    project.selected_dataset=("dr11","south","1554p180","g")
    project.target_radius_arcsec=50.
    save_object_project(project)
    image=np.ones((241,241),dtype=float)
    data={"project_id":project.project_id,"band":"g",
          "source":"native-test","wcs":_wcs(),"pixel_shape":image.shape,
          "native_hard":np.zeros_like(image,bool),
          "data":{"image":image,"maskbits":np.zeros_like(image,np.uint32)}}
    rows=(_red_extended_row(1794),_red_extended_row(1946,blue=True))
    cat=SimpleNamespace(rows=rows,
        pixel_positions=np.array([[160.,120.],[130.,130.]]),
        source_path="native-tractor.fits",source_sha256="a"*64)
    monkeypatch.setattr(module,"read_native_tractor_gaia",
                        lambda *args,**kwargs:cat)
    draft=SkyMaskDraft(project.project_id)
    old=module.prepare_scene_contamination(
        data,project,draft,(155.,18.))
    new=module.prepare_scene_contamination(
        data,project,draft,(155.,18.),red_extended_preview=True)
    assert old["stats"]["masked_catalog_sources"]==0
    assert new["stats"]["masked_catalog_sources"]==1
    assert new["stats"]["red_extended_preview_sources"]==1
    assert new["red_extended_preview_enabled"] is True
    assert new["gaia_inner_preview_enabled"] is False
    assert new["stats"]["red_extended_preview_pixels"]>0
    assert new["mask"][120,160] and not old["mask"][120,160]
    assert not new["mask"][130,130]
    assert np.array_equal(new["mask"],new["red_extended_preview_mask"])
    assert "EXPLICIT_RED_EXTENDED_MASK_PREVIEW_ONLY" in new["flags"]
    assert "COLOR_MORPHOLOGY_CLASSIFIER_PURITY_NOT_VALIDATED" in new["flags"]
    e={v["objid"]:v for v in new["source_inventory"]["entries"]}
    assert e[1794]["screening_disposition"]=="RED_EXTENDED_MASK_PREVIEW"
    assert e[1794]["diagnostic_red_extended_mask_preview_active"] is True
    assert e[1946]["screening_disposition"]=="PROTECTED_OVERLAP"
    assert e[1946]["diagnostic_red_extended_mask_preview_active"] is False
    assert all(v["physical_membership"]=="UNDETERMINED" for v in e.values())
    assert all(v["photometry_action"]=="NO_NEW_ACTION_READ_ONLY"
               for v in e.values())
    assert all(v["proposed_treatment_is_active"] is False
               for v in e.values())
    assert draft.edits==()
    assert np.array_equal(image,data["data"]["image"])



def test_red_extended_isolation_uses_wcs_not_scalar_pixel_scale():
    """A nearby interloper must not appear isolated in an anisotropic WCS.

    Here the neighbour is five native pixels away, but just one
    arcsecond in the actual world tangent plane. The former scalar
    pixel-scale method would incorrectly pass the 2.5-arcsec veto.
    """
    wcs = _wcs()
    wcs.wcs.cd = np.diag([-1., .2])/3600.
    protected = np.ones((241,241), dtype=bool)
    candidate = _red_extended_row(1794)
    neighbour = _row("PSF", flux=50.)
    positions = np.array([[160.,120.], [160.,125.]])
    mask, stats = conservative_catalogue_mask(
        protected.shape, positions, (candidate,neighbour),
        arcsec_per_pixel=1., protected=protected, wcs=wcs,
        red_extended_preview=True,anchor_pixel=(120.,120.),
        target_radius_arcsec=60.)
    assert stats["red_extended_preview_sources"] == 0
    assert stats["skipped_target_overlap"] == 2
    assert not mask.any()
    # Removing the close neighbour makes exactly the same red object
    # pass, with a sky-plane offset accounting for cos(dec) and scale.
    alone, counts = conservative_catalogue_mask(
        protected.shape, positions[:1], (candidate,),
        arcsec_per_pixel=1., protected=protected, wcs=wcs,
        red_extended_preview=True,anchor_pixel=(120.,120.),
        target_radius_arcsec=60.)
    assert counts["red_extended_preview_sources"] == 1
    evidence = counts["red_extended_preview_centres"][0]["evidence"]
    assert evidence["distance_geometry"] == (
        "NATIVE_WCS_LOCAL_EAST_NORTH_ARCSEC")
    assert 35. < evidence["offset_arcsec"] < 40.
    assert evidence["nearest_arcsec"] is None
    assert evidence["no_catalogue_neighbour"] is True
    assert alone[120,160]
    # Fail closed if caller attempts to classify using no native WCS.
    with pytest.raises(ValueError,match="requires native celestial WCS"):
        conservative_catalogue_mask(
            protected.shape, positions[:1], (candidate,),
            arcsec_per_pixel=1., protected=protected,
            red_extended_preview=True,anchor_pixel=(120.,120.),
            target_radius_arcsec=60.)


def test_red_extended_policy_records_evidence_for_selected_source():
    """All decisions are reversible and metrics survive Source evidence."""
    import lsbpipe.photometry.scene_source_inventory as si
    protected = np.ones((241,241), dtype=bool)
    candidate = _red_extended_row(2290)
    points = np.array([[170.,120.]])
    mask, stats = conservative_catalogue_mask(
        protected.shape, points, (candidate,),
        arcsec_per_pixel=1., protected=protected, wcs=_wcs(),
        red_extended_preview=True,anchor_pixel=(120.,120.),
        target_radius_arcsec=60.)
    cat = SimpleNamespace(
        pixel_positions=points, rows=(candidate,),
        source_path="native-tractor.fits",source_sha256="c"*64)
    inv = si.build_scene_source_inventory(
        cat,pixel_shape=protected.shape,protected=protected,
        screen_stats=stats,release="dr11",region="south",
        brick="1554p180",band="g",anchor_icrs=(155.,18.))
    source = inv["entries"][0]
    assert source["screening_disposition"] == "RED_EXTENDED_MASK_PREVIEW"
    evidence = source["red_extended_rule_evidence"]
    assert evidence["policy"]=="red_extended_v1"
    assert evidence["passes"] is True
    assert evidence["distance_geometry"] == (
        "NATIVE_WCS_LOCAL_EAST_NORTH_ARCSEC")
    assert evidence["snr_g"] >= 10
    assert source["physical_membership"] == "UNDETERMINED"
    assert source["photometry_action"] == "NO_NEW_ACTION_READ_ONLY"
    assert source["proposed_treatment_is_active"] is False
    assert "MULTIBAND_RED_EXTENDED_ISOLATED_HYPOTHESIS" in (
        source["evidence_flags"])
    assert int(mask.sum()) > 0
