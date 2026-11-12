"""SC1 native WCS/ICRS scene editor: functional offline Qt regressions."""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pytest
from astropy.wcs import WCS

pytest.importorskip("pytestqt")
pytest.importorskip("PySide6")

from lsbpipe.gui.photometry_workspace import PhotometryWorkspace
from lsbpipe.gui.scene_editor import (
    SceneEditor, scene_shape_from_world_points, scene_shape_outline,
)
from lsbpipe.photometry.mask_drafts import (
    read_active_mask_draft, save_mask_draft,
)
from lsbpipe.photometry.sky_masks import (
    SkyMaskDraft, SkyMaskEdit, SkyMaskShape, evaluate_mask_draft,
)
from lsbpipe.surveys.legacy import SkyField
from lsbpipe.surveys.object_projects import (
    create_object_project, save_object_project)


def _native_wcs(*, flipped=False):
    w = WCS(naxis=2)
    w.wcs.crpix = [49.,49.]
    w.wcs.crval = [359.999, 60.]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.cd = np.array(
        [[.26, .055], [-.040, -.26]] if flipped else
        [[-.26, .055], [.040, .26]], dtype=float)/3600.
    return w


def _data(project_id, *, band="g", wcs=None):
    size=96
    yy,xx = np.indices((size,size))
    image = (yy+.25*xx)*1e-4
    hard = np.zeros((size,size),dtype=bool)
    hard[48,48] = True
    rgba=np.zeros((size,size,4),dtype=np.uint8)
    rgba[hard] = (238,89,104,160)
    return {
        "status":"NATIVE_SCENE_CUTOUT_NOT_SCIENCE_ADOPTED",
        "data":{"image":image,"invvar":np.ones_like(image),
                "nexp":np.ones_like(image),
                "maskbits":np.zeros_like(image,dtype=np.uint32)},
        "wcs":wcs or _native_wcs(),
        "pixel_shape":(size,size),
        "rgb":np.full((size,size,3),132,dtype=np.uint8),
        "hard_rgba":rgba,
        "native_hard":hard,
        "band":band,
        "source":f"DR11 south · 1554p600 · {band}",
        "display_width_arcmin":6.,
        "project_field_width_arcmin":6.,
        "project_id":project_id,
        "hdu_indices":{"image":1,"invvar":1,"nexp":1,"maskbits":1},
    }


def _click_at_native(editor,x,y):
    # The same FITS/Qt vertical inversion used by the real image viewer.
    editor._on_click(x+.5,editor._loaded["pixel_shape"][0]-y-.5)


def _make_editor(qtbot,tmp_path,monkeypatch):
    monkeypatch.setenv("LSBPIPE_SESSION_DIR",str(tmp_path/"session"))
    p=create_object_project("Leo P",tmp_path/"Leo_P",
                            field=SkyField(359.999,60.,6.))
    w=SceneEditor()
    qtbot.addWidget(w)
    w.resize(1100,700)
    w.show()
    w.bind_project(p.root)
    assert w._draft.project_id==p.project_id
    assert not (p.root/"analyses").exists()
    w._accept_dataset(_data(p.project_id))
    return w,p


def test_same_uuid_project_centre_rebind_invalidates_loaded_scene(
        qtbot,tmp_path,monkeypatch):
    """Moving a project target must discard previous native pixels/masks."""
    scene,project=_make_editor(qtbot,tmp_path,monkeypatch)
    original_key=scene._source_key
    assert scene._loaded is not None
    assert scene._project_target == (359.999, 60.)
    project.field=SkyField(359.9995, 60.00025, 5.)
    save_object_project(project)
    scene.bind_project(project.root)
    assert scene._project_id==project.project_id
    assert scene._project_target==(359.9995,60.00025)
    assert scene._source_key!=original_key
    assert scene._loaded is None
    assert scene._contamination is None
    assert scene._isophote_proposal is None
    assert scene._foreground_psf_result is None
    assert scene._draft is not None
    # Even an asynchronous load result from the old passport cannot
    # be displayed as if it used the new target anchor.
    old=_data(project.project_id)
    old["loaded_target_icrs"]=[359.999,60.]
    with pytest.raises(ValueError,match="target ICRS differs"):
        scene._accept_dataset(old)
    assert scene._loaded is None


def test_same_uuid_target_protection_radius_change_invalidates_scene(
        qtbot,tmp_path,monkeypatch):
    """Source masks depend on the target protection radius, not just RA/Dec."""
    scene,project=_make_editor(qtbot,tmp_path,monkeypatch)
    previous_key=scene._source_key
    previous_target=scene._project_target
    original_radius=float(project.target_radius_arcsec)
    assert scene._loaded is not None
    assert scene._draft is not None
    project.target_radius_arcsec=float(project.target_radius_arcsec)+12.
    save_object_project(project)
    scene.bind_project(project.root)
    assert scene._source_key!=previous_key
    assert scene._project_target==previous_target
    assert scene._loaded is None
    assert scene._contamination is None
    assert scene._isophote_proposal is None
    assert scene._foreground_psf_result is None
    assert scene._draft is not None
    # Rebinding the unchanged passport is idempotent.
    new_key=scene._source_key
    scene.bind_project(project.root)
    assert scene._source_key==new_key
    # Changing the radius back must NOT resurrect an old draft from the
    # previous source-key cache after a newer edit was saved meanwhile.
    newer=scene.draft.with_edit(SkyMaskEdit(
        "newer_object_scope_edit", "exclude",
        SkyMaskShape("circle", (359.999,60.), radius_arcsec=4.)))
    newer_sha=save_mask_draft(
        project.root,newer,expected_digest=scene._saved_digest)
    project.target_radius_arcsec=original_radius
    save_object_project(project)
    scene.bind_project(project.root)
    assert scene._source_key==previous_key
    assert scene.draft.digest==newer_sha
    assert len(scene.draft.edits)==1
    assert scene.draft.edits[0].edit_id=="newer_object_scope_edit"
    assert scene._loaded is None


def test_failed_autosave_blocks_same_object_source_switch_without_losing_edits(
        qtbot,tmp_path,monkeypatch):
    """A stale external HEAD must never erase the live unsaved draft."""
    scene,project=_make_editor(qtbot,tmp_path,monkeypatch)
    original_key=scene._source_key
    original_scene=scene._loaded
    original_target=scene._project_target
    unsaved=scene.draft.with_edit(SkyMaskEdit(
        "live_not_saved", "exclude",
        SkyMaskShape("circle",(359.999,60.),radius_arcsec=3.)))
    scene._draft=unsaved
    external=SkyMaskDraft(project.project_id).with_edit(SkyMaskEdit(
        "other_window", "protect_target",
        SkyMaskShape("circle",(359.999,60.),radius_arcsec=5.)))
    external_sha=save_mask_draft(
        project.root,external,expected_digest=scene._saved_digest)
    project.target_radius_arcsec+=12.
    save_object_project(project)
    scene.bind_project(project.root)
    assert scene._source_key==original_key
    assert scene._loaded is original_scene
    assert scene._project_target==original_target
    assert scene._draft==unsaved
    assert scene._draft.digest!=scene._saved_digest
    assert "SOURCE SWITCH BLOCKED" in scene.status.text()
    scene.bind_project(None)
    assert "PROJECT UNBIND BLOCKED" in scene.status.text()
    assert scene._source_key==original_key
    assert scene._loaded is original_scene
    assert scene._draft==unsaved
    saved,sha=read_active_mask_draft(project.root)
    assert sha==external_sha
    assert saved==external


def test_scene_rejects_cross_brick_worker_result_without_mutating_gui(
        qtbot,tmp_path,monkeypatch):
    scene,project=_make_editor(qtbot,tmp_path,monkeypatch)
    old_dataset=scene._loaded
    dataset=_data(project.project_id)
    dataset["source_quartet_key"]=["dr11","south","1555p600","g"]
    with pytest.raises(ValueError,match="quartet differs"):
        scene._accept_dataset(dataset)
    assert scene._loaded is old_dataset


def test_scene_is_real_second_stage_without_auto_loading_fits(
        qtbot,tmp_path,monkeypatch):
    monkeypatch.setenv("LSBPIPE_SESSION_DIR",str(tmp_path/"session"))
    p=create_object_project("Leo P",tmp_path/"LeoP",
                            field=SkyField(359.999,60.,6.))
    import lsbpipe.gui.scene_editor as module
    calls=[]
    monkeypatch.setattr(
        module,"load_scene_original_cutout",
        lambda *args,**kwargs:calls.append((args,kwargs)))
    view=PhotometryWorkspace()
    qtbot.addWidget(view)
    view.show()
    view.bind_project(p)
    assert view.stage_tabs.count()==4
    assert view.stage_tabs.tabText(0)=="1 · PSF-0"
    assert view.stage_tabs.tabText(1)=="2 · Scene & Masks"
    assert view.stage_tabs.tabText(2)=="3 · Background"
    assert view.stage_tabs.tabText(3)=="4 · Aperture Photometry"
    assert view.aperture_stage._result is None
    assert view.sky_stage._result is None
    assert view.scene_editor._project_id==p.project_id
    view.stage_tabs.setCurrentWidget(view.scene_editor)
    assert calls==[]
    assert view.scene_editor._loaded is None
    assert not (p.root/"analyses").exists()
    assert view.scene_editor.draft.edits==()
    assert "NOT validated" in view.scene_editor.warning.text()


def test_scene_loader_records_source_quartet_hashes_and_detects_race(
        tmp_path, monkeypatch):
    """A loaded Scene carries exact original file hashes, not filenames."""
    import lsbpipe.gui.scene_editor as module
    from lsbpipe.surveys.dataset_integrity import _expected_paths
    from lsbpipe.surveys.object_projects import save_object_project
    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path/"session"))
    project=create_object_project(
        "Leo P", tmp_path/"Leo_P",
        field=SkyField(359.999, 60., 6.))
    project.selected_dataset=("dr11","south","1554p600","r")
    project.data_root=tmp_path/"originals"
    project.data_root.mkdir()
    save_object_project(project)
    paths=_expected_paths(project.data_root, project.selected_dataset)
    from astropy.io import fits
    pixels=np.arange(96*96,dtype=float).reshape(96,96)*1e-6
    source_arrays={
        "image":pixels.astype("f4"),
        "invvar":np.ones(pixels.shape,dtype="f4"),
        "nexp":np.ones(pixels.shape,dtype="i2"),
        "maskbits":np.zeros(pixels.shape,dtype="i2"),
    }
    for kind,path in paths.items():
        path.parent.mkdir(parents=True,exist_ok=True)
        header=_native_wcs().to_header()
        if kind=="image":
            header["BUNIT"]="nanomaggies"
        if kind=="invvar":
            header["BUNIT"]="1/nanomaggies^2"
        fits.PrimaryHDU(data=source_arrays[kind],header=header).writeto(path)
    originals={name:file.read_bytes() for name,file in paths.items()}
    data={"image":pixels,"invvar":np.ones_like(pixels),
          "nexp":np.ones_like(pixels),
          "maskbits":np.zeros_like(pixels,dtype=np.int16)}
    def fake_read(*, directory, field, hit, band):
        assert band=="r"
        return {
            "data_root":str(project.data_root),
            "data":data, "wcs":_native_wcs(),
            "pixel_shape":pixels.shape,
            "source":"DR11 south · 1554p600 · r",
            "hdu_indices":{kind:0 for kind in paths},
        }
    monkeypatch.setattr(module, "prepare_fits_explorer",fake_read)
    loaded=module.load_scene_original_cutout(str(project.root),"r",4.)
    assert loaded["original_fits_provenance_status"] == (
        "FULL_ORIGINAL_FITS_SHA256_VERIFIED_AT_SCENE_LOAD")
    assert loaded["original_fits_sha256_is_archive_authenticity"] is False
    assert loaded["image_unit_evidence"] == (
        "FITS_IMAGE_BUNIT_NANOMAGGIES_CONFIRMED")
    assert loaded["invvar_unit_evidence"] == "FITS_INVVAR_BUNIT_CONFIRMED"
    assert loaded["original_fits_sha256_local"] == {
        kind:hashlib.sha256(payload).hexdigest()
        for kind,payload in originals.items()}
    assert all(file.read_bytes()==originals[kind]
               for kind,file in paths.items())
    with fits.open(paths["image"],mode="update") as hdus:
        hdus[0].header["BUNIT"]="ADU"
    with pytest.raises(ValueError, match="not nanomaggies"):
        module.load_scene_original_cutout(str(project.root),"r",4.)
    with fits.open(paths["image"],mode="update") as hdus:
        del hdus[0].header["BUNIT"]
    with pytest.raises(ValueError, match="IMAGE BUNIT missing"):
        module.load_scene_original_cutout(str(project.root),"r",4.)
    with fits.open(paths["image"],mode="update") as hdus:
        hdus[0].header["BUNIT"]="nanomaggies"
    with fits.open(paths["invvar"],mode="update") as hdus:
        del hdus[0].header["BUNIT"]
    assumed=module.load_scene_original_cutout(str(project.root),"r",4.)
    assert assumed["invvar_unit_evidence"] == (
        "INVVAR_UNITS_FROM_SURVEY_DOCUMENTATION_ASSUMPTION")
    def racing_reader(**kwargs):
        receipt=fake_read(**kwargs)
        paths["image"].write_bytes(b"changed-original-during-load")
        return receipt
    monkeypatch.setattr(module,"prepare_fits_explorer",racing_reader)
    with pytest.raises(RuntimeError,match="changed during Scene load"):
        module.load_scene_original_cutout(str(project.root),"r",4.)


def test_icrs_circle_ellipse_and_wraparound_on_native_rotated_wcs():
    w=_native_wcs(flipped=True)
    # Each shape is constructed from *world* positions sampled via the
    # actual WCS, never by interpreting RA delta as an east angle.
    c=(float(x)%360. for x in w.all_pix2world(48.,48.,0))
    center=tuple(c)
    endpoint=tuple(map(float,w.all_pix2world(62.,48.,0)))
    endpoint=(endpoint[0]%360.,endpoint[1])
    circle=scene_shape_from_world_points(
        "circle",(center,endpoint))
    ellipse=scene_shape_from_world_points(
        "ellipse",(center,endpoint),axis_ratio=.43)
    assert circle.radius_arcsec==pytest.approx(ellipse.a_arcsec)
    assert ellipse.q==.43
    assert 0<=ellipse.pa_deg<180
    for shape in (circle,ellipse):
        outline=scene_shape_outline(shape,w,96)
        assert outline.boundingRect().width()>5.
        assert outline.boundingRect().height()>5.
        assert np.isfinite(outline.boundingRect().center().x())
    vertices=[
        tuple(map(float,w.all_pix2world(x,y,0)))
        for x,y in [(34.,34.),(65.,34.),(60.,65.),(34.,65.)]]
    vertices=[(ra%360.,dec) for ra,dec in vertices]
    polygon=scene_shape_from_world_points("polygon",tuple(vertices))
    assert len(polygon.vertices)==4
    assert scene_shape_outline(polygon,w,96).boundingRect().width()>10.


def test_live_scene_edits_undo_redo_review_save_and_reopen(
        qtbot,tmp_path,monkeypatch):
    w,p=_make_editor(qtbot,tmp_path,monkeypatch)
    native=tmp_path/"unmodified_original.fits"
    native.write_bytes(b"DO NOT WRITE NATIVE FITS")
    before=hashlib.sha256(native.read_bytes()).hexdigest()
    w.mode.setCurrentIndex(w.mode.findData("ellipse"))
    w.axis_ratio.setValue(.55)
    _click_at_native(w,48,48)
    _click_at_native(w,70,53)
    assert len(w.draft.edits)==1
    ellipse=w.draft.edits[0]
    assert ellipse.purpose=="protect_target"
    assert ellipse.shape.kind=="ellipse"
    assert ellipse.shape.q==.55
    assert len(w._paths)==1
    w.mode.setCurrentIndex(w.mode.findData("polygon"))
    w.purpose.setCurrentIndex(w.purpose.findData("exclude"))
    w.scope.setCurrentIndex(w.scope.findData("band"))
    for x,y in [(23,22),(34,22),(36,38),(23,37)]:
        _click_at_native(w,x,y)
    w.finish_polygon()
    assert len(w.draft.edits)==2
    assert w.draft.edits[-1].band=="g"
    assert w.draft.edits[-1].shape.kind=="polygon"
    w.undo()
    assert len(w.draft.edits)==1
    w.redo()
    assert len(w.draft.edits)==2
    w.edits.setCurrentRow(0)
    w.remove_selected()
    assert len(w.draft.edits)==1
    w.undo()
    assert len(w.draft.edits)==2

    w.save_draft()
    loaded,digest=read_active_mask_draft(p.root)
    assert loaded==w.draft
    assert digest==w._saved_digest
    assert "Not science adopted" in w.status.text()
    oldfile=p.root/"analyses"/"mask-drafts"/f"draft-{digest}.json"
    assert oldfile.exists()
    assert hashlib.sha256(native.read_bytes()).hexdigest()==before

    other=SceneEditor()
    qtbot.addWidget(other)
    other.bind_project(p.root)
    assert other.draft==loaded
    assert other._loaded is None  # loading saved scene never opens FITS
    other._accept_dataset(_data(p.project_id))
    assert len(other._paths)==2
    # Hard invalid pixels remain prohibited even in protected target.
    qa=evaluate_mask_draft(
        other.draft,other._loaded["wcs"],(96,96),"g",
        native_hard=other._loaded["native_hard"])
    assert qa["native_hard"][48,48]
    assert qa["excluded"][48,48]
    assert qa["status"]=="MASK_DRAFT_NOT_ADOPTED"
    assert not (p.root/"qa-runs").exists()


def test_band_local_native_overlays_and_shared_icrs_geometry(
        qtbot,tmp_path,monkeypatch):
    w,p=_make_editor(qtbot,tmp_path,monkeypatch)
    w.mode.setCurrentIndex(w.mode.findData("circle"))
    w.scope.setCurrentIndex(w.scope.findData("all"))
    _click_at_native(w,48,48)
    _click_at_native(w,65,48)
    circle=w.draft.edits[0]
    assert circle.band is None
    gpath=scene_shape_outline(circle.shape,_native_wcs(),96)
    w.set_band("r")
    assert w.draft.edits==(circle,)
    assert w._loaded is None
    assert len(w._paths)==0
    w._accept_dataset(_data(p.project_id,band="r",wcs=_native_wcs(flipped=True)))
    assert len(w._paths)==1
    rpath=scene_shape_outline(circle.shape,_native_wcs(flipped=True),96)
    assert rpath.boundingRect().width()>0
    assert abs(gpath.boundingRect().center().x() -
               rpath.boundingRect().center().x())<2
    assert "R band" not in w.status.text() or "original" in w.status.text()
    w.set_band("g")
    assert w._loaded is None
    assert w.draft.edits==(circle,)


def test_optimistic_concurrency_refuses_to_overwrite_other_scene_draft(
        qtbot,tmp_path,monkeypatch):
    w,p=_make_editor(qtbot,tmp_path,monkeypatch)
    w.mode.setCurrentIndex(w.mode.findData("circle"))
    _click_at_native(w,48,48)
    _click_at_native(w,58,48)
    w.save_draft()
    saved,sha=read_active_mask_draft(p.root)
    external=saved.with_edit(SkyMaskEdit(
        "other_researcher", "exclude",
        SkyMaskShape("circle",(359.999,60.),radius_arcsec=3.)))
    other_sha=save_mask_draft(p.root,external,expected_digest=sha)
    _click_at_native(w,50,40)
    _click_at_native(w,66,41)
    w.save_draft()
    assert "NOT saved" in w.status.text()
    current,digest=read_active_mask_draft(p.root)
    assert digest==other_sha
    assert current==external
    assert len(w.draft.edits)==2


def test_explicit_asynchronous_load_drops_stale_band_result(
        qtbot,tmp_path,monkeypatch):
    w,p=_make_editor(qtbot,tmp_path,monkeypatch)
    import lsbpipe.gui.scene_editor as module
    from threading import Event
    started=Event()
    release=Event()
    def slow_loader(*args,**kwargs):
        started.set()
        release.wait(2)
        return _data(p.project_id,band=args[1])
    monkeypatch.setattr(module,"load_scene_original_cutout",slow_loader)
    w.load_original()
    assert started.wait(2)
    w.set_band("r")
    release.set()
    qtbot.waitUntil(lambda: w._future is None, timeout=4000)
    assert w._loaded is None
    assert w.draft is not None
    assert w.load_btn.isEnabled()
    w.load_original()
    qtbot.waitUntil(lambda: w._loaded is not None, timeout=4000)
    assert w._loaded["band"]=="r"
    w.shutdown()



def test_redraw_selected_region_preview_and_identity_undo(
        qtbot,tmp_path,monkeypatch):
    editor,project=_make_editor(qtbot,tmp_path,monkeypatch)
    editor.mode.setCurrentIndex(editor.mode.findData("polygon"))
    for px,py in [(20,20),(38,20),(38,39),(20,39)]:
        _click_at_native(editor,px,py)
        assert len(editor._pending_scene)>0
        assert editor._preview_path is not None
        assert editor._preview_path.scene() is editor.graphics
    editor.finish_polygon()
    assert editor._preview_path is None
    assert editor._pending_world==[]
    original=editor.draft.edits[0]
    original_geometry=original.shape.vertices
    editor.edits.setCurrentRow(0)
    assert editor.edit_btn.isEnabled()
    editor.edit_btn.click()
    assert editor._replace_index==0
    for px,py in [(26,28),(54,28),(54,53),(26,53)]:
        _click_at_native(editor,px,py)
    editor.finish_polygon()
    assert len(editor.draft.edits)==1
    edited=editor.draft.edits[0]
    assert edited.edit_id==original.edit_id
    assert edited.shape.vertices!=original_geometry
    assert editor._preview_path is None
    editor.undo()
    assert editor.draft.edits==(original,)
    editor.redo()
    assert editor.draft.edits==(edited,)
    editor.save_draft()
    reopened,sha=read_active_mask_draft(project.root)
    assert reopened.edits==(edited,)
    assert sha==editor._saved_digest



def test_scene_circle_outline_uses_true_spherical_bearing_and_arcsec_radius():
    """First circular outline vertex is at north, not at a misread PA."""
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    w=_native_wcs(flipped=True)
    center=(359.999,60.)
    radius=8.
    shape=SkyMaskShape("circle",center,radius_arcsec=radius)
    path=scene_shape_outline(shape,w,96)
    north=SkyCoord(center[0]*u.deg,center[1]*u.deg).spherical_offsets_by(
        0*u.deg,radius*u.arcsec)
    expected_x,expected_y=w.all_world2pix(north.ra.deg,north.dec.deg,0)
    first=path.elementAt(0)
    assert first.x==pytest.approx(expected_x+.5,abs=1e-5)
    assert first.y==pytest.approx(96.-expected_y-.5,abs=1e-5)
    center_x,center_y=w.all_world2pix(*center,0)
    assert path.boundingRect().center().x()==pytest.approx(
        center_x+.5,abs=.1)
    assert path.boundingRect().center().y()==pytest.approx(
        96.-center_y-.5,abs=.1)
    assert path.boundingRect().width()>40.
    assert path.boundingRect().height()>40.


def test_preliminary_isophote_preview_is_not_saved_until_explicit_accept(
        qtbot, tmp_path, monkeypatch):
    editor, project = _make_editor(qtbot, tmp_path, monkeypatch)
    data = _data(project.project_id)
    yy, xx = np.indices((96,96))
    rng = np.random.default_rng(71)
    # Native coadd in nMgy/pix with a resolved, faint elongated source
    # at the exact ICRS centre plus unrelated compact sources.
    ang = np.deg2rad(25.)
    major = np.cos(ang)*(xx-48.) + np.sin(ang)*(yy-48.)
    minor = -np.sin(ang)*(xx-48.) + np.cos(ang)*(yy-48.)
    data["data"]["image"] = (
        .7 + rng.normal(0., .18, size=xx.shape)
        + 1.7*np.exp(-.5*((major/9.)**2+(minor/5.)**2)))
    editor._accept_dataset(data)
    original = data["data"]["image"].copy()
    assert editor.iso_prepare.isEnabled()
    assert not editor.iso_generate.isEnabled()  # no raw-image fallback
    assert editor.draft.edits == ()

    import lsbpipe.gui.scene_editor as module
    from types import SimpleNamespace
    source_masks = []

    def fake_prepare(dataset, root, draft, centre):
        mask = np.zeros(dataset["pixel_shape"], bool)
        mask[5:12, 78:85] = True
        source_masks.append((dataset["band"], root))
        return {
            "status":"CATALOGUE_CONTAMINATION_MASK_NOT_ADOPTED",
            "mask":mask, "auto_mask":mask.copy(),
            "protected":np.zeros_like(mask),
            "hard":dataset["native_hard"].copy(),
            "project_id":project.project_id,
            "draft_digest":draft.digest,
            "dataset_source":dataset["source"], "band":dataset["band"],
            "anchor_icrs":centre, "catalog_sha256":"a"*64,
            "stats":{"masked_catalog_sources":1,
                     "skipped_target_overlap":0,
                     "excluded_pixels":49, "native_hard_pixels":1},
            "flags":("PSF_WINGS_NOT_VALIDATED",),
        }

    monkeypatch.setattr(module, "prepare_scene_neighbour_masks",
                        fake_prepare)
    editor.iso_prepare.click()
    qtbot.waitUntil(
        lambda: editor._future is None
        and editor._contamination is not None, timeout=10000)
    assert source_masks == [("g", str(project.root))]
    assert editor.iso_generate.isEnabled()
    assert editor._contamination_overlay is not None
    assert editor.draft.edits == ()

    editor.iso_smooth.setValue(2.)
    editor.iso_threshold.setValue(2.5)
    editor.iso_generate.click()
    qtbot.waitUntil(
        lambda: editor._future is None
        and editor._isophote_proposal is not None, timeout=10000)
    assert editor._isophote_preview is not None
    assert editor.iso_summary.text()
    assert "area ratio" in editor.iso_summary.text()
    assert editor.iso_candidate.count() >= 2
    assert editor.iso_candidate.currentData() == 0
    reference = editor._selected_isophote()["contour"]
    editor.iso_candidate.setCurrentIndex(1)
    alternative = editor._selected_isophote()["contour"]
    assert alternative != reference
    assert editor._isophote_preview is not None
    assert editor.draft.edits == ()
    editor.iso_candidate.setCurrentIndex(0)
    if not editor.iso_accept.isEnabled():
        assert editor._isophote_proposal["sensitivity"][
            "confidence"].startswith("LOW_CONFIDENCE")
        editor._accept_isophote()  # direct method cannot bypass the gate
        assert editor.draft.edits == ()
        editor.iso_low_conf_override.setChecked(True)
    assert editor.iso_accept.isEnabled()
    assert editor.draft.edits == ()
    assert not (project.root/"analyses").exists()
    assert "NOT a calibrated outer isophote" in editor.iso_details.text()
    assert editor._isophote_proposal["contour"].kind == "polygon"
    assert editor._isophote_proposal["ellipse"].kind == "ellipse"
    editor.iso_shape.setCurrentIndex(editor.iso_shape.findData("ellipse"))
    assert editor._isophote_preview is not None
    editor.iso_discard.click()
    assert editor._isophote_proposal is None
    assert editor.draft.edits == ()
    assert not (project.root/"analyses").exists()

    editor.iso_generate.click()
    qtbot.waitUntil(
        lambda: editor._future is None
        and editor._isophote_proposal is not None, timeout=10000)
    editor.iso_shape.setCurrentIndex(editor.iso_shape.findData("contour"))
    assert editor.iso_candidate.count() >= 2
    editor.iso_candidate.setCurrentIndex(1)
    chosen = editor._selected_isophote()["contour"]
    chosen_trial = editor.iso_candidate.currentData()
    if not editor.iso_accept.isEnabled():
        editor.iso_low_conf_override.setChecked(True)
    editor.iso_accept.click()
    assert len(editor.draft.edits) == 1
    edit = editor.draft.edits[0]
    assert edit.shape.kind == "polygon"
    assert edit.shape == chosen
    assert f"trial={chosen_trial}" in edit.comment
    assert edit.purpose == "isophote_candidate"
    assert edit.band == "g"  # native-band default for source-derived contours
    assert editor.purpose.findData("isophote_candidate") >= 0
    assert "no mask" in editor.iso_accept.text().lower()
    assert "AUTO_ISOPHOTE_PRELIMINARY" in edit.comment
    assert np.array_equal(data["data"]["image"],original)
    assert editor._isophote_proposal is None
    editor.save_draft()
    saved, _ = read_active_mask_draft(project.root)
    assert saved.edits == (edit,)
    editor.undo()
    assert editor.draft.edits == ()
    editor.redo()
    assert editor.draft.edits == (edit,)
    editor.set_band("r")
    editor._accept_dataset(_data(project.project_id, band="r",
                                 wcs=_native_wcs(flipped=True)))
    assert len(editor._paths) == 0  # never leak g segmentation to r implicitly
    assert editor._contamination is None
    assert not editor.iso_generate.isEnabled()


def test_isophote_worker_result_is_discarded_after_band_switch(
        qtbot, tmp_path, monkeypatch):
    editor, project = _make_editor(qtbot, tmp_path, monkeypatch)
    import lsbpipe.gui.scene_editor as module
    from threading import Event
    started, release = Event(), Event()

    def delayed(*args):
        started.set()
        release.wait(2)
        return {"status":"LOCAL_ISOPHOTE_PROPOSAL_NOT_ADOPTED",
                "band":"g", "project_id":project.project_id,
                "source":"old source"}
    monkeypatch.setattr(module, "propose_scene_isophote", delayed)
    def prepare(dataset, root, draft, center):
        mask = np.zeros(dataset["pixel_shape"],bool)
        return {
            "status":"CATALOGUE_CONTAMINATION_MASK_NOT_ADOPTED",
            "mask":mask, "auto_mask":mask,
            "project_id":project.project_id,
            "dataset_source":dataset["source"], "band":dataset["band"],
            "draft_digest":draft.digest, "anchor_icrs":center,
            "catalog_sha256":"f"*64,
            "stats":{"masked_catalog_sources":0,
                     "skipped_target_overlap":0,
                     "excluded_pixels":0, "native_hard_pixels":1},
            "flags":(),
        }
    monkeypatch.setattr(module, "prepare_scene_neighbour_masks",prepare)
    editor.iso_prepare.click()
    qtbot.waitUntil(
        lambda: editor._future is None and editor._contamination is not None,
        timeout=4000)
    editor.iso_generate.click()
    assert started.wait(2.)
    editor.set_band("r")
    release.set()
    qtbot.waitUntil(lambda: editor._future is None, timeout=4000)
    assert editor._isophote_proposal is None
    assert editor.draft.edits == ()


def test_isophote_share_to_all_bands_requires_explicit_opt_in(qtbot,tmp_path,
                                                                monkeypatch):
    editor, project=_make_editor(qtbot,tmp_path,monkeypatch)
    from lsbpipe.photometry.sky_masks import SkyMaskShape
    from types import SimpleNamespace
    assert not editor.iso_share_bands.isChecked()
    shape=SkyMaskShape("ellipse",center=(359.999,60.),a_arcsec=8.,
                       q=.7,pa_deg=20.)
    editor._contamination={
        "catalog_sha256":"d"*64,"draft_digest":editor.draft.digest}
    editor._isophote_proposal={
        "status":"LOCAL_ISOPHOTE_PROPOSAL_NOT_ADOPTED",
        "band":"g","project_id":project.project_id,
        "source":editor._loaded["source"],
        "catalog_sha256":"d"*64,"draft_digest":editor.draft.digest,
        "ellipse":shape,"contour":shape,
        "bridge_guard_pix":2,
        "sensitivity":{"confidence":"REVIEWABLE_PRELIMINARY","flags":()},
        "proposal":SimpleNamespace(
            threshold_contrast=2.2,smoothing_pix=4.,
            local_sky=0.,smoothed_sky_rms=.1)}
    editor.iso_shape.setCurrentIndex(editor.iso_shape.findData("ellipse"))
    editor._show_isophote_preview()
    assert editor.iso_accept.isEnabled()
    editor.iso_share_bands.setChecked(True)
    editor.iso_accept.click()
    assert editor.draft.edits[0].band is None
    assert "scope=all_optical_EXPLICIT" in editor.draft.edits[0].comment
    editor.set_band("r")
    editor._accept_dataset(_data(project.project_id,band="r",
                                 wcs=_native_wcs(flipped=True)))
    assert len(editor._paths)==1


def test_mask_opacity_and_conflict_review_markers(qtbot,tmp_path,monkeypatch):
    editor,project=_make_editor(qtbot,tmp_path,monkeypatch)
    data=_data(project.project_id)
    visible=np.zeros(data["pixel_shape"],bool)
    visible[20:30,20:30]=True
    result={
       "status":"CATALOGUE_CONTAMINATION_MASK_NOT_ADOPTED",
       "mask":visible,"auto_mask":visible.copy(),
       "project_id":project.project_id,"dataset_source":data["source"],
       "band":"g","draft_digest":editor.draft.digest,
       "anchor_icrs":editor._project_target,"catalog_sha256":"b"*64,
       "stats":{"masked_catalog_sources":1,
                "elliptical_catalog_sources":1,
                "skipped_target_overlap":1,
                "near_target_sources":(
                    {"x":48.,"y":48.,"objid":77,
                     "type":"PSF","radius_pix":9.},),
                "excluded_pixels":100,"native_hard_pixels":1},
       "flags":("PSF_WINGS_NOT_VALIDATED",)}
    editor._accept_contamination_result(result)
    assert len(editor._conflict_paths)==1
    assert editor._contamination_overlay.isVisible()
    assert not editor.iso_show_conflicts.isChecked()
    assert not editor._conflict_paths[0].isVisible()
    editor.iso_show_conflicts.setChecked(True)
    assert editor._conflict_paths[0].isVisible()
    editor.iso_opacity.setValue(20)
    assert editor._contamination_overlay.opacity()==pytest.approx(.2)
    editor.iso_mask_toggle.setChecked(False)
    assert not editor._contamination_overlay.isVisible()
    editor.iso_show_conflicts.setChecked(False)
    assert not editor._conflict_paths[0].isVisible()
    editor._discard_prepared_masks()
    assert len(editor._conflict_paths)==0
    assert not editor.iso_generate.isEnabled()


def test_sga_choice_explicit_core_protection_and_wcs_preview(
        qtbot,tmp_path,monkeypatch):
    from lsbpipe.photometry.sga_priors import SGAGeometryPrior,SGAPriorLookup
    import lsbpipe.gui.scene_editor as module
    editor,project=_make_editor(qtbot,tmp_path,monkeypatch)
    p=SGAGeometryPrior(
        release="SGA2020",sga_id=72,name="Leo P",
        ra_deg=359.999,dec_deg=60.,offset_arcsec=0.,
        d26_arcmin=.26,d26_ref="SB26",ba=.6,pa_deg=45.,
        ellipsebit=0,ellipsemode=0,flags=(),
        source_sha256="f"*64,source_path="localSGA.fits")
    second=SGAGeometryPrior(
        release="SGA2020",sga_id=73,name="Other target",
        ra_deg=359.9992,dec_deg=60.,offset_arcsec=.4,
        d26_arcmin=.10,d26_ref="LEDA",ba=.9,pa_deg=90.,
        ellipsebit=0,ellipsemode=0,
        flags=("D26_NOT_DIRECT_SB26",),
        source_sha256="f"*64,source_path="localSGA.fits")
    result=SGAPriorLookup(
        "SGA2020","localSGA.fits","f"*64,(p,second))
    editor._accept_sga_lookup(result)
    assert editor._sga_candidate is None
    assert editor._sga_item is None
    assert editor.sga_candidates.count()==3
    assert not editor.sga_use_prior.isChecked()
    editor.sga_candidates.setCurrentIndex(1)
    assert editor._sga_candidate==p
    assert editor._sga_item is not None
    assert "CENTRAL PRIOR" in editor.sga_status.text()
    assert not editor.sga_use_prior.isChecked()
    calls=[]
    def fake_prep(*args):
        calls.append(args)
        image=editor._loaded
        zeros=np.zeros(image["pixel_shape"],bool)
        return {
          "status":"CATALOGUE_CONTAMINATION_MASK_NOT_ADOPTED",
          "mask":zeros,"auto_mask":zeros.copy(),
          "project_id":project.project_id,"band":image["band"],
          "dataset_source":image["source"],
          "draft_digest":editor.draft.digest,
          "anchor_icrs":editor._project_target,
          "catalog_sha256":"a"*64,
          "sga_prior_id":p.sga_id,
          "sga_prior_sha256":p.source_sha256,
          "stats":{"masked_catalog_sources":0,
                   "skipped_target_overlap":0,
                   "possible_sga_shreds":0,
                   "excluded_pixels":0,"native_hard_pixels":1},
          "flags":("SGA_CORE_ONLY_NOT_OUTER_LSB_BOUNDARY",)}
    monkeypatch.setattr(module,"prepare_scene_neighbour_masks",fake_prep)
    editor.sga_use_prior.setChecked(True)
    assert editor._contamination is None
    editor.iso_prepare.click()
    qtbot.waitUntil(
        lambda: editor._future is None and editor._contamination is not None,
        timeout=5000)
    assert len(calls)==1 and calls[0][-1] is p
    assert not editor._contamination["mask"].any()
    editor.sga_candidates.setCurrentIndex(2)
    assert editor._sga_candidate is second
    assert not editor.sga_use_prior.isChecked()
    assert editor._contamination is None
    editor.set_band("r")
    editor._accept_dataset(_data(project.project_id,band="r",
                                 wcs=_native_wcs(flipped=True)))
    assert editor._sga_item is not None
    assert editor._sga_candidate is second
    assert not editor.sga_use_prior.isChecked()
    assert editor.draft.edits==()


def test_no_sga_matches_does_not_autocreate_protect_regions(
        qtbot,tmp_path,monkeypatch):
    from lsbpipe.photometry.sga_priors import SGAPriorLookup
    editor,project=_make_editor(qtbot,tmp_path,monkeypatch)
    editor._accept_sga_lookup(
        SGAPriorLookup("SGA2020","local","e"*64,()))
    assert editor.sga_candidates.count()==1
    assert editor._sga_candidate is None
    assert not editor.sga_use_prior.isEnabled()
    assert editor.draft.edits==()
    assert "No SGA match" in editor.sga_status.text()


def test_on_demand_desi_sga2025_never_needs_local_full_catalogue(
        qtbot,tmp_path,monkeypatch):
    from lsbpipe.photometry.sga_priors import (
        SGAGeometryPrior, SGAPriorLookup,
    )
    from lsbpipe.surveys.object_projects import save_object_project
    import lsbpipe.gui.scene_editor as module
    editor,project=_make_editor(qtbot,tmp_path,monkeypatch)
    project.selected_dataset=("dr11","south","1554p600","g")
    save_object_project(project)
    # The public Viewer returns a provisional ellipse: neither verified
    # final SGAID nor the full SGA2025 D26 photometry/fitting flags.
    prior=SGAGeometryPrior(
        release="SGA2025",sga_id=731,name="SGA-731",
        ra_deg=359.999,dec_deg=60.,offset_arcsec=0.,
        d26_arcmin=.30,d26_ref="VIEWER_RADIUS_ARCSEC",
        ba=.7,pa_deg=55.,ellipsebit=0,ellipsemode=0,
        flags=("SGA2025_VIEWER_OVERLAY_NOT_MERGED_CATALOG",),
        source_sha256="d"*64,
        source_path="https://www.legacysurvey.org/viewer/"
                    "sga2025-south/1/cat.json?test=1")
    lookup=SGAPriorLookup(
        "SGA2025",prior.source_path,prior.source_sha256,(prior,))
    calls=[]
    def fake_fetch(coords,*,region,radius_arcmin):
        calls.append((coords,region,radius_arcmin))
        return lookup
    monkeypatch.setattr(module,"fetch_viewer_sga2025_priors",
                        fake_fetch)
    assert editor.sga_find.isEnabled()
    editor.sga_find.click()
    qtbot.waitUntil(
        lambda: editor._future is None and editor._sga_lookup is not None,
        timeout=5000)
    assert calls==[((359.999,60.),"south",3.)]
    assert editor.sga_candidates.count()==2
    assert editor._sga_candidate is None
    assert editor.draft.edits==()
    editor.sga_candidates.setCurrentIndex(1)
    assert editor._sga_candidate is prior
    assert editor._sga_item is not None
    assert "Viewer overlay diameter" in editor.sga_status.text()
    assert "Viewer ref" in editor.sga_status.text()
    assert "D26" not in editor.sga_status.text().split("q=")[0]
    assert not editor.sga_use_prior.isChecked()
    assert editor._contamination is None
    editor.set_band("r")
    editor._accept_dataset(_data(project.project_id,band="r",
                                 wcs=_native_wcs(flipped=True)))
    assert editor._sga_item is not None
    assert editor.draft.edits==()



def test_leo_p_zero_radius_viewer_rejected_then_parent_requested_explicitly(
        qtbot, tmp_path, monkeypatch):
    from lsbpipe.photometry.sga_priors import (
        SGAGeometryPrior, SGAPriorLookup,
    )
    from lsbpipe.surveys.object_projects import save_object_project
    import lsbpipe.gui.scene_editor as module

    editor, project = _make_editor(qtbot, tmp_path, monkeypatch)
    invalid = SGAPriorLookup(
        "SGA2025", "https://www.legacysurvey.org/viewer/"
        "sga2025-south/1/cat.json", "a"*64, (),
        unusable_count=1,
        unusable_examples=("SGA-0 (Viewer ref 0, 0.1 arcsec): "
                           "missing/nonpositive overlay radius",))
    editor._accept_sga_lookup(invalid)
    assert editor.sga_candidates.count() == 1
    assert editor._sga_candidate is None
    assert editor._sga_item is None
    assert not editor.sga_use_prior.isEnabled()
    assert "ZERO usable initial ellipses" in editor.sga_status.text()
    assert "SGA-0" in editor.sga_status.text()
    assert editor.draft.edits == ()

    project.selected_dataset = ("dr11", "south", "1554p600", "g")
    save_object_project(project)
    parent = SGAGeometryPrior(
        release="SGA2025-PARENT", sga_id=42, name="Leo P",
        ra_deg=359.999, dec_deg=60., offset_arcsec=0.,
        d26_arcmin=.42, d26_ref="VIEWER_PARENT_DIAMETER",
        ba=.65, pa_deg=80., ellipsebit=0, ellipsemode=0,
        flags=("PARENT_INITIAL_GEOMETRY_NOT_FITTED_ISOPHOTE",),
        source_sha256="c"*64,
        source_path="https://www.legacysurvey.org/viewer/"
                    "sga2025-parent/1/cat.json")
    lookup = SGAPriorLookup(
        "SGA2025-PARENT", parent.source_path, parent.source_sha256,
        (parent,))
    seen = []
    def fake_fetch(coords, **kwargs):
        seen.append((coords, kwargs))
        return lookup
    monkeypatch.setattr(module, "fetch_viewer_sga2025_priors",
                        fake_fetch)
    editor.sga_parent_find.click()
    qtbot.waitUntil(
        lambda: editor._future is None and editor._sga_lookup is lookup,
        timeout=5000)
    assert len(seen) == 1
    assert seen[0][1]["layer"] == "parent"
    assert seen[0][1]["region"] == "south"
    assert editor.sga_candidates.count() == 2
    editor.sga_candidates.setCurrentIndex(1)
    assert editor._sga_item is not None
    assert "Parent initial diameter" in editor.sga_status.text()
    assert "Parent PGC" in editor.sga_status.text()
    assert not editor.sga_use_prior.isChecked()
    assert editor.draft.edits == ()


def test_sga_informed_radial_initial_ellipse_requires_explicit_review(
        qtbot,tmp_path,monkeypatch):
    from lsbpipe.photometry.sga_priors import SGAGeometryPrior,SGAPriorLookup
    import lsbpipe.gui.scene_editor as module
    editor,project=_make_editor(qtbot,tmp_path,monkeypatch)
    data=_data(project.project_id,band="g")
    n=321
    xx,yy=np.meshgrid(np.arange(n),np.arange(n))
    rng=np.random.default_rng(140)
    image=(1.+rng.normal(0.,.15,(n,n))
           +2.4*np.exp(-.5*(((xx-160.)/29.)**2+
                             ((yy-160.)/19.)**2)))
    data["data"]["image"]=image
    data["data"]["invvar"]=np.ones_like(image)
    data["data"]["nexp"]=np.ones_like(image)
    data["data"]["maskbits"]=np.zeros(image.shape,np.uint32)
    data["pixel_shape"]=image.shape
    data["rgb"]=np.full((*image.shape,3),128,np.uint8)
    data["native_hard"]=np.zeros(image.shape,bool)
    data["hard_rgba"]=np.zeros((*image.shape,4),np.uint8)
    w=_native_wcs()
    w.wcs.crpix=[161.,161.]
    data["wcs"]=w
    editor._accept_dataset(data)
    original=image.copy()
    prior=SGAGeometryPrior(
        release="SGA2025-PARENT",sga_id=5065058,name="Leo P",
        ra_deg=359.999,dec_deg=60.,offset_arcsec=0.,
        d26_arcmin=2.*1.8/60.,
        d26_ref="VIEWER_PARENT_DIAMETER",
        ba=.65,pa_deg=140.,ellipsebit=0,ellipsemode=0,
        flags=("PARENT_INITIAL_GEOMETRY_NOT_FITTED_ISOPHOTE",),
        source_sha256="c"*64,source_path="https://www.legacysurvey.org/"
        "viewer/sga2025-parent/1/cat.json")
    lookup=SGAPriorLookup(
        "SGA2025-PARENT",prior.source_path,prior.source_sha256,(prior,))
    editor._accept_sga_lookup(lookup)
    editor.sga_candidates.setCurrentIndex(1)
    assert editor.sga_use_seed.isChecked()
    assert not editor.sga_use_prior.isChecked()
    assert editor.draft.edits==()

    def fake_prepare(dataset,root,draft,center):
        mask=np.zeros(dataset["pixel_shape"],bool)
        mask[145:156,182:192]=True
        return {
            "status":"CATALOGUE_CONTAMINATION_MASK_NOT_ADOPTED",
            "project_id":project.project_id,"band":dataset["band"],
            "dataset_source":dataset["source"],"draft_digest":draft.digest,
            "anchor_icrs":center,"catalog_sha256":"a"*64,
            "mask":mask,"auto_mask":mask.copy(),
            "sga_prior_id":None,
            "stats":{"masked_catalog_sources":1,
                     "masked_catalog_centres":(
                         {"x":186.,"y":150.,"objid":42,"type":"PSF",
                          "gaia_counterpart":False},),
                     "skipped_target_overlap":1,
                     "near_target_sources":(
                         {"x":160.,"y":160.,"objid":73,
                          "type":"PSF","radius_pix":12.,
                          "classification":"UNRESOLVED_TARGET_BLEND"},),
                     "excluded_pixels":int(mask.sum()),
                     "native_hard_pixels":0},
            "flags":("PSF_WINGS_NOT_VALIDATED",),
        }

    monkeypatch.setattr(module,"prepare_scene_neighbour_masks",fake_prepare)
    editor.iso_prepare.click()
    qtbot.waitUntil(
        lambda: editor._future is None and editor._contamination is not None,
        timeout=10000)
    # Inject inventory-only evidence AFTER the GUI accepts the
    # ordinary minimal test fixture; no source-registry widget contract
    # is mocked by this targeted radial-count regression.
    editor._contamination["source_inventory"]={"entries":(
        {"x":160.,"y":159.,"objid":88,"tractor_type":"PSF",
         "source_id":"dr11/south/test/OBJID=88",
         "screening_disposition":"NOT_SCREENED_OR_UNFITTED"},)}
    editor.iso_smooth.setValue(3.)
    editor.iso_generate.click()
    qtbot.waitUntil(
        lambda: editor._future is None
        and editor._isophote_proposal is not None,timeout=20000)
    proposal=editor._isophote_proposal
    assert proposal["sga_comparison"]["sga_ref"]==5065058
    assert proposal["sga_comparison"]["reference_q"]==pytest.approx(.65)
    assert proposal["sga_comparison"]["reference_pa_deg"]==pytest.approx(140.)
    assert proposal["sga_radial"] is not None
    assert proposal["sga_seed_shape"] is not None
    assert proposal["sga_comparison"][
        "unresolved_catalog_centres_in_initial_ellipse"]==1
    assert proposal["sga_comparison"][
        "unscreened_catalog_centres_in_initial_ellipse"]==1
    assert proposal["sga_comparison"][
        "all_unmasked_catalog_centres_in_initial_ellipse"]==2
    assert "1 protected overlaps + 1 unscreened/unfitted" in (
        editor.iso_sga_summary.text())
    assert proposal["sga_seed_shape"].a_arcsec>prior.core_shape().a_arcsec
    assert "SGA INITIAL:" in editor.iso_sga_summary.text()
    assert "Independent moments:" in editor.iso_sga_summary.text()
    assert "no auto restores or new exclusions" in (
        editor.iso_sga_summary.text())
    review=proposal["sga_comparison"].get("review_candidates",())
    assert len(editor._sga_shred_paths)==len(review)
    assert not editor.iso_show_sga_shreds.isChecked()
    if review:
        assert not editor._sga_shred_paths[0].isVisible()
        editor.iso_show_sga_shreds.setChecked(True)
        assert editor._sga_shred_paths[0].isVisible()
    assert editor.iso_shape.findData("sga_seed")>=0
    editor.iso_shape.setCurrentIndex(
        editor.iso_shape.findData("sga_seed"))
    assert editor._isophote_preview is not None
    assert "SGA INITIAL (not a fitted isophote)" in editor.iso_summary.text()
    assert "q=0.650" in editor.iso_summary.text()
    assert "PA=140.0" in editor.iso_summary.text()
    assert "nan" not in editor.iso_summary.text()
    assert not editor.iso_accept.isEnabled()
    editor._accept_isophote()
    assert editor.draft.edits==()
    editor.iso_low_conf_override.setChecked(True)
    assert editor.iso_accept.isEnabled()
    shape=proposal["sga_seed_shape"]
    editor.iso_accept.click()
    assert editor.draft.edits[0].shape==shape
    assert editor.draft.edits[0].band=="g"
    assert "kind=sga_seed" in editor.draft.edits[0].comment
    assert "SGAref=5065058" in editor.draft.edits[0].comment
    assert "SGAsha256=" + "c"*64 in editor.draft.edits[0].comment
    assert "TractorSHA256=" + "a"*64 in editor.draft.edits[0].comment
    assert "MANUAL_REVIEW_NO_SCIENCE_ADOPTION" in editor.draft.edits[0].comment
    assert np.array_equal(data["data"]["image"],original)
    assert editor._isophote_proposal is None
    assert editor._sga_shred_paths==[]


def test_single_detected_contour_trial_does_not_falsely_claim_stability():
    from lsbpipe.gui.scene_editor import scene_contour_stability_readout
    # The Leo P screenshot showed 1 successful reference trial and
    # 3 failures, giving a trivial 1.00 area ratio and NaN IoU.
    result = {
        "trials": ({"state":"DETECTED"},)
                  + ({"state":"NO_STABLE_COMPONENT"},)*3,
        "area_ratio":1.0,
        "minimum_spatial_iou":None,
    }
    assert scene_contour_stability_readout(result) == (
        1,4,"n/a","n/a")
    result = {
        "trials":({"state":"DETECTED"},{"state":"DETECTED"},
                  {"state":"NO_STABLE_COMPONENT"},
                  {"state":"NO_STABLE_COMPONENT"}),
        "area_ratio":1.5,
        "minimum_spatial_iou":.73,
    }
    assert scene_contour_stability_readout(result) == (
        2,4,"1.50","0.73")


def test_sga_protected_60arcsec_radius_coincidence_is_not_claimed_as_real_edge():
    from lsbpipe.photometry.ellipse_seed import SGAInformedRadialSeed
    from lsbpipe.gui.scene_editor import scene_sga_mask_boundary_warning
    radial=SGAInformedRadialSeed(
        reference_a_arcsec=14.28,detected_a_arcsec=60.2,
        q_sky=.532,pa_deg=150.5,centre_icrs=(359.999,60.),
        bins=(),masked_fraction=0.,
        flags=("RADIAL_EXTENT_NEAR_PROVISIONAL_TARGET_PROTECTION_RADIUS",
               "MASK_BOUNDARY_SENSITIVITY_NOT_TESTED"),
        target_protection_radius_arcsec=60.)
    warning=scene_sga_mask_boundary_warning(radial,60.)
    assert "60.0 arcsec" in warning
    assert "MASK-DEPENDENT" in warning
    assert "does not prove a galaxy boundary" in warning
    assert not scene_sga_mask_boundary_warning(radial,None)
    other=SGAInformedRadialSeed(
        reference_a_arcsec=14.28,detected_a_arcsec=60.2,
        q_sky=.532,pa_deg=150.5,centre_icrs=(359.999,60.),
        bins=(),masked_fraction=0.,flags=())
    assert not scene_sga_mask_boundary_warning(other,60.)


def test_readonly_source_inventory_shows_all_ambiguous_and_unfitted_gaia(
        qtbot, tmp_path, monkeypatch):
    from types import SimpleNamespace
    from PySide6.QtWidgets import QApplication
    from lsbpipe.photometry.scene_source_inventory import (
        build_scene_source_inventory,
    )
    from lsbpipe.photometry.scene_contaminants import (
        conservative_catalogue_mask,
    )
    editor, project = _make_editor(qtbot, tmp_path, monkeypatch)
    data = _data(project.project_id)
    points = np.array([[48.,48.],[75.,75.],[45.,45.]])
    protected = np.zeros(data["pixel_shape"],bool)
    protected[38:59,38:59] = True
    rows = (
        {"objid":40, "type":"PSF","flux":15.,"flux_ivar":12.,
         "ref_cat":"G3","ref_id":4567,"gaia_g_mag":20.,
         "parallax_mas":1.2,"parallax_ivar":100.,
         "pmra_masyr":None,"pmra_ivar":None,
         "maskbits":(1<<17)|(1<<12),"fitbits":1<<0},
        {"objid":41,"type":"EXP","flux":120.,"flux_ivar":8.,
         "ref_cat":"","ref_id":0,"maskbits":0,"fitbits":0,
         "shape_r_arcsec":2.},
        {"objid":42,"type":"PSF","flux":45.,"flux_ivar":0.,
         "ref_cat":"G3","ref_id":9999,"gaia_g_mag":16.5,
         "maskbits":1<<1,"fitbits":1<<4},
    )
    mask, stats = conservative_catalogue_mask(
        data["pixel_shape"],points,rows,
        arcsec_per_pixel=1.,protected=protected)
    assert stats["skipped_target_overlap"]==1
    assert stats["masked_catalog_sources"]==1
    assert stats["unusable_catalog_sources"]==1
    catalog = SimpleNamespace(
        rows=rows,pixel_positions=points,
        source_path="local/tractor-1554p600.fits",
        source_sha256="d"*64,
        dataset_id="dr11/south/1554p600/g")
    evidence = build_scene_source_inventory(
        catalog,pixel_shape=data["pixel_shape"],protected=protected,
        screen_stats=stats,release="dr11",region="south",
        brick="1554p600",band="g")
    result = {
        "status":"CATALOGUE_CONTAMINATION_MASK_NOT_ADOPTED",
        "mask":mask,"auto_mask":mask.copy(),
        "project_id":project.project_id,
        "dataset_source":data["source"],
        "band":"g", "draft_digest":editor.draft.digest,
        "anchor_icrs":editor._project_target,
        "catalog_sha256":"d"*64,
        "source_inventory":evidence,
        "stats":{**stats,
                 "excluded_pixels":int(mask.sum()),
                 "native_hard_pixels":1},
        "flags":("CATALOGUE_GEOMETRY_NOT_TRACTOR_MODEL_RENDER",),
    }
    editor._accept_contamination_result(result)
    assert editor.source_registry_toggle.isEnabled()
    assert editor.source_registry_list.count()==3
    assert "3 shown of 3" in editor.source_registry_status.text()
    assert editor.source_registry_details.toPlainText()
    editor.source_registry_toggle.setChecked(True)
    editor.source_registry_filter.setCurrentIndex(
        editor.source_registry_filter.findData("PROTECTED_OVERLAP"))
    assert editor.source_registry_list.count()==1
    assert "OBJID 40" in editor.source_registry_list.item(0).text()
    assert "UNDETERMINED" in editor.source_registry_details.toPlainText()
    assert "GAIA_ASTROMETRIC_SIGNAL_FOREGROUND_REVIEW" in (
        editor.source_registry_details.toPlainText())
    editor.source_registry_filter.setCurrentIndex(
        editor.source_registry_filter.findData("TARGET_PROTECTION_ALL"))
    assert editor.source_registry_list.count()==2
    editor.source_registry_filter.setCurrentIndex(
        editor.source_registry_filter.findData("TARGET_UNSCREENED"))
    assert editor.source_registry_list.count()==1
    assert "OBJID 42" in editor.source_registry_list.item(0).text()
    editor.source_registry_filter.setCurrentIndex(
        editor.source_registry_filter.findData("NOT_SCREENED_OR_UNFITTED"))
    assert editor.source_registry_list.count()==1
    assert "OBJID 42" in editor.source_registry_list.item(0).text()
    assert "RETAINED_GAIA_SOURCE_POSSIBLE" in (
        editor.source_registry_details.toPlainText())
    assert editor.source_registry_filter.findData("RED_POLICY_UNMASKED") >= 0
    assert editor.source_registry_filter.findData("RED_POLICY_NEAR_MISS") >= 0
    # Deliberately synthetic GUI-only evidence; the physics predicate
    # is independently tested with real-valued multi-band inputs.
    reviewed = evidence["entries"][0]
    reviewed["red_extended_first_failure"]="G_MINUS_R_OUTSIDE_STRICT_RANGE"
    reviewed["red_extended_near_miss"]=True
    reviewed["red_extended_rule_evidence"]={
        "policy":"red_extended_v1", "passes":False,
        "first_failure":"G_MINUS_R_OUTSIDE_STRICT_RANGE",
        "g_minus_r":.75,"r_minus_z":.70}
    editor.source_registry_filter.setCurrentIndex(
        editor.source_registry_filter.findData("RED_POLICY_NEAR_MISS"))
    assert editor.source_registry_list.count()==1
    assert "RED NEAR-MISS" in editor.source_registry_list.item(0).text()
    assert "G_MINUS_R_OUTSIDE_STRICT_RANGE" in (
        editor.source_registry_details.toPlainText())
    editor.source_registry_filter.setCurrentIndex(
        editor.source_registry_filter.findData("RED_POLICY_UNMASKED"))
    assert editor.source_registry_list.count()==1
    editor.source_registry_filter.setCurrentIndex(
        editor.source_registry_filter.findData("NOT_SCREENED_OR_UNFITTED"))
    editor.source_registry_copy.click()
    exported = QApplication.clipboard().text()
    assert "source_catalog_sha256\t" + "d"*64 in exported
    assert "dr11/south/1554p600/OBJID=42" in exported
    assert "NO_NEW_ACTION_READ_ONLY" in exported
    assert "OBJID=40" not in exported
    editor.source_registry_search.setText("nonexistent_source")
    assert editor.source_registry_list.count()==0
    assert not editor.source_registry_copy.isEnabled()
    editor._discard_prepared_masks()
    assert editor._source_registry is None
    assert editor._source_registry_visible==[]
    assert not editor.source_registry_toggle.isEnabled()
    assert not editor.iso_generate.isEnabled()
    assert editor.draft.edits==()
    assert not (project.root/"analyses").exists()



def test_one_click_complete_scene_qa_copies_all_trials_and_catalogue(
        qtbot, tmp_path, monkeypatch):
    """Full clipboard snapshot is read-only, includes failed sensitivity
    trials and full source metadata, but excludes native pixel arrays."""
    import json
    from PySide6.QtWidgets import QApplication
    from lsbpipe.photometry.ellipse_seed import (
        EllipseSeed, LocalIsophotalProposal,
    )

    editor, project = _make_editor(qtbot, tmp_path, monkeypatch)
    record = {
        "status": "CATALOGUE_CONTAMINATION_MASK_NOT_ADOPTED",
        "band": "g",
        "project_id": project.project_id,
        "catalog_sha256": "e"*64,
        "dataset_source": editor._loaded["source"],
        "draft_digest": editor.draft.digest,
        "stats": {
            "masked_catalog_sources": 493,
            "skipped_target_overlap": 114,
            "near_target_sources": ({"objid": 2035, "x": 48., "y": 51.},),
        },
        "flags": ("PSF_WINGS_NOT_VALIDATED",),
        "mask": np.ones((96, 96), dtype=bool),
        "auto_mask": np.ones((96, 96), dtype=bool),
    }
    editor._contamination = record
    editor._source_registry = {
        "schema": "test-source-evidence",
        "catalog_sha256": "e"*64,
        "stored_count": 2,
        "total_in_native_frame": 2,
        "entries": (
            {"objid": 2035, "source_id": "dr11/south/1554p600/OBJID=2035",
             "physical_membership": "UNDETERMINED",
             "evidence_flags": ("GAIA_ASTROMETRIC_SIGNAL_FOREGROUND_REVIEW",)},
            {"objid": 2076, "source_id": "dr11/south/1554p600/OBJID=2076",
             "physical_membership": "UNDETERMINED",
             "evidence_flags": ("SGA_GALAXY_CONTEXT_NOT_MEMBERSHIP_PROOF",)},
        ),
    }
    seed = EllipseSeed(x_pix=48., y_pix=48., sma_pix=13., q=.5,
                       theta_native_deg=30., component_pixels=150,
                       centroid_offset_pix=1.)
    proposal = LocalIsophotalProposal(
        seed=seed, contour_xy=((40., 40.), (60., 40.), (50., 60.)),
        local_sky=-.001, raw_sky_rms=.002, smoothed_sky_rms=.0004,
        sky_pixels=1600, threshold_contrast=2.2, smoothing_pix=4.,
        flags=("BACKGROUND_CONSTANT_NOT_VALIDATED",),
        component_indices=np.arange(10000, dtype=np.int64))
    trials = (
        {"threshold": 2.2, "smooth_pix": 4., "state": "DETECTED",
         "connected_pixels": 150, "overlap_to_reference": 1.},
        {"threshold": 1.7, "smooth_pix": 4., "state": "DETECTED",
         "connected_pixels": 200, "overlap_to_reference": .6},
        {"threshold": 2.85, "smooth_pix": 4., "state": "NO_STABLE_COMPONENT",
         "reason": "Disconnected pixel boundary"},
        {"threshold": 2.2, "smooth_pix": 5.4, "state": "DETECTED",
         "connected_pixels": 180, "overlap_to_reference": .8},
    )
    editor._isophote_proposal = {
        "status": "LOCAL_ISOPHOTE_PROPOSAL_NOT_ADOPTED",
        "proposal": proposal,
        "alternatives": ({"model": proposal}, None, None, None),
        "sensitivity": {"confidence": "LOW_CONFIDENCE", "trials": trials,
                        "minimum_spatial_iou": .6,
                        "flags": ("CONTOUR_FAILS_NEARBY_THRESHOLDS",)},
        "catalog_sha256": "e"*64,
    }
    editor.iso_summary.setText("LOW CONFIDENCE: no science adoption")
    editor.iso_details.setText("Independent PSF wings not validated")
    before_draft = editor.draft
    before_image = editor._loaded["data"]["image"].copy()
    editor.copy_scene_qa_btn.click()
    snapshot = json.loads(QApplication.clipboard().text())
    assert snapshot["schema"] == "desidera.gui.scene-complete-qa/1"
    assert snapshot["science_status"] == "UNADOPTED_DIAGNOSTIC_ONLY"
    assert snapshot["native_fits"]["native_pixel_shape"] == [96, 96]
    assert snapshot["contaminant_screening"]["stats"]["skipped_target_overlap"] == 114
    assert snapshot["contaminant_screening"]["proposed_mask_pixel_count"] == 96*96
    assert len(snapshot["source_inventory"]["entries"]) == 2
    assert snapshot["source_inventory"]["entries"][1]["objid"] == 2076
    assert snapshot["source_inventory"]["all_stored_entries_included"]
    qa = snapshot["isophote_proposal"]["sensitivity"]
    assert len(qa["trials"]) == 4
    assert qa["trials"][2]["state"] == "NO_STABLE_COMPONENT"
    assert "Disconnected" in qa["trials"][2]["reason"]
    assert snapshot["isophote_proposal"]["proposal"]["component_indices"] == {
        "pixel_array_omitted": True, "shape": [10000], "dtype": "int64",
        "element_count": 10000,
    }
    assert snapshot["display_messages"]["isophote_details_all_trials"] == (
        "Independent PSF wings not validated")
    assert editor.draft == before_draft
    assert np.array_equal(editor._loaded["data"]["image"], before_image)
    assert not (project.root/"analyses").exists()


def test_scene_qa_button_without_loaded_fits_stays_readonly(
        qtbot, tmp_path, monkeypatch):
    """Unprepared GUI exports explicit absence, never guesses FITS results."""
    import json
    from PySide6.QtWidgets import QApplication

    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path/"session"))
    project = create_object_project(
        "Leo P", tmp_path/"Leo_P", field=SkyField(359.999, 60., 6.))
    editor = SceneEditor()
    qtbot.addWidget(editor)
    editor.bind_project(project.root)
    assert editor._loaded is None
    editor.copy_scene_qa_btn.click()
    qa = json.loads(QApplication.clipboard().text())
    assert qa["native_fits"] is None
    assert qa["isophote_proposal"] is None
    assert qa["contaminant_screening"] is None
    assert qa["project"]["project_id"] == project.project_id
    assert qa["science_status"] == "UNADOPTED_DIAGNOSTIC_ONLY"
    assert not (project.root/"analyses").exists()



def test_radius_sweep_is_non_destructive_and_tracks_mask_isophote_changes(
        qtbot, tmp_path, monkeypatch):
    """Only target protection radius differs, on one immutable band/scene.

    Source membership and photometric adoption are unchanged. Failure
    scenarios remain in the denominator instead of disappearing.
    """
    import lsbpipe.gui.scene_editor as module
    from lsbpipe.photometry.ellipse_seed import EllipseSeed, LocalIsophotalProposal

    editor, project = _make_editor(qtbot, tmp_path, monkeypatch)
    image_before = editor._loaded["data"]["image"].copy()
    original_radius = project.target_radius_arcsec
    draft_before = editor.draft.digest
    observed = []
    def fake_prepare(dataset, scenario_project, draft, anchor, *,
                     sga_prior=None, **kwargs):
        radius = scenario_project.target_radius_arcsec
        observed.append(("mask",radius,anchor))
        mask = np.zeros(dataset["pixel_shape"],dtype=bool)
        mask[:int(radius/10),:4] = True
        return {
            "status":"CATALOGUE_CONTAMINATION_MASK_NOT_ADOPTED",
            "project_id":project.project_id, "band":"g",
            "dataset_source":dataset["source"],
            "draft_digest":draft.digest,"anchor_icrs":anchor,
            "mask":mask,"catalog_sha256":"a"*64,
            "source_inventory":{
                "total_relevant_sources":3, "truncated":0,
                "screening_evidence_incomplete":False,
                "entries": (
                    {"source_id":"1","screening_disposition":
                     "PROPOSED_EXCLUDE"},
                    {"source_id":"2","screening_disposition":
                     ("PROTECTED_OVERLAP" if radius >= 60.
                      else "PROPOSED_EXCLUDE")},
                    {"source_id":"3","screening_disposition":
                     "NOT_SCREENED_OR_UNFITTED"},
                ),
            },
            "stats":{
                "masked_catalog_sources":(1 if radius>=60. else 2),
                "skipped_target_overlap":(1 if radius>=60. else 0),
                "excluded_pixels":int(mask.sum()),
                "protected_pixels":int(radius*4),
                "usable_pixels":mask.size-int(mask.sum()),
                "native_hard_pixels":1,
            },
        }
    def fake_propose(dataset, anchor, contrast, smooth, prepared,
                     bridge_guard_pix=0, sga_prior=None):
        radius = float(prepared["stats"]["protected_pixels"])/4
        observed.append(("detect",radius,anchor))
        if radius == 80.:
            raise ValueError("No connected emission at 80 arcsec")
        n=int(radius)
        component=np.arange(n,dtype=np.int64)
        seed=EllipseSeed(
            x_pix=48., y_pix=48., sma_pix=12., q=.55,
            theta_native_deg=30., component_pixels=n,
            centroid_offset_pix=float(radius/10))
        model=LocalIsophotalProposal(
            seed=seed,contour_xy=((42.,42.),(56.,42.),(48.,56.)),
            local_sky=float(radius)/1e6,raw_sky_rms=.002,
            smoothed_sky_rms=.0002,sky_pixels=2400,
            threshold_contrast=contrast,smoothing_pix=smooth,
            flags=("BACKGROUND_CONSTANT_NOT_VALIDATED",),
            component_indices=component)
        return {
            "proposal":model,"sga_radial":None,"sga_warning":None,
            "sensitivity":{
                "area_ratio":1.2,"minimum_spatial_iou":.7,
                "confidence":"LOW_CONFIDENCE_REQUIRES_EXPLICIT_OVERRIDE",
                "flags":("UNANCHORED_CENTROID_DRIFT",),
                "trials":(
                    {"threshold":contrast,"smooth_pix":smooth,
                     "connected_pixels":n,"state":"DETECTED"},
                    {"threshold":contrast-.5,"smooth_pix":smooth,
                     "state":"DETECTED"},
                    {"threshold":contrast+.65,"smooth_pix":smooth,
                     "state":"DETECTED"},
                    {"threshold":contrast,"smooth_pix":smooth*1.35,
                     "state":"DETECTED"}),
            },
        }
    monkeypatch.setattr(module,"prepare_scene_contamination",fake_prepare)
    monkeypatch.setattr(module,"propose_scene_isophote",fake_propose)
    result=module.run_scene_protection_radius_sweep(
        editor._loaded, project.root, editor.draft,
        editor._project_target,
        threshold_contrast=2.2, smoothing_pix=4.,
        bridge_guard_pix=2)
    assert result["schema"]=="desidera.gui.protection-radius-sweep/1"
    assert result["status"]=="READ_ONLY_DIAGNOSTIC_NOT_SCIENCE_VALIDATED"
    assert [r["protection_radius_arcsec"] for r in result["scenarios"]] == [
        40.,50.,60.,70.,80.]
    assert len(observed)==10  # 5 preparations, 5 detection attempts
    assert result["catalog_sha256"]=="a"*64
    assert all(r["catalog_sha256"]=="a"*64 for r in result["scenarios"])
    baseline=result["scenarios"][2]
    assert baseline["masked_pixel_change_vs_60arcsec"]==0
    assert baseline["mask_jaccard_vs_60arcsec"]==pytest.approx(1.)
    assert baseline["component_iou_vs_60arcsec"]==pytest.approx(1.)
    assert baseline["component_sky_nMgy_pixel"]==pytest.approx(60e-6)
    assert baseline["new_excluded_ids_vs_60arcsec"]==[]
    assert result["scenarios"][0]["new_excluded_ids_vs_60arcsec"]==["2"]
    assert result["scenarios"][3]["component_iou_vs_60arcsec"] < 1.
    assert result["scenarios"][4]["status"]=="ISOPHOTE_NOT_FOUND"
    assert "80 arcsec" in result["scenarios"][4]["isophote_error"]
    assert result["scenarios"][4]["component_iou_vs_60arcsec"] is None
    # Compare neighbouring scenarios directly; relative-to-60 overlap
    # cannot prove two large-radius components are pixel-identical.
    assert result["scenarios"][0]["previous_radius_arcsec"] is None
    assert result["scenarios"][0]["component_iou_vs_previous_radius"] is None
    assert baseline["previous_radius_arcsec"]==50.
    assert baseline["component_iou_vs_previous_radius"]==pytest.approx(
        50./60.)
    assert baseline["new_protected_overlap_ids_vs_previous_radius"]==["2"]
    assert result["scenarios"][3]["previous_radius_arcsec"]==60.
    assert result["scenarios"][3]["component_iou_vs_previous_radius"]==(
        pytest.approx(60./70.))
    assert result["scenarios"][4]["previous_radius_arcsec"]==70.
    assert result["scenarios"][4]["component_iou_vs_previous_radius"] is None
    assert result["scenarios"][4]["mask_jaccard_vs_previous_radius"] < 1.
    assert result["scenarios"][4]["sga_radius_delta_vs_previous_arcsec"] is None
    assert len(baseline["all_four_trials"])==4
    assert "NO_PSF_WING_MODEL_SCENARIO_EXECUTED" in result["limitations"]
    assert "NO_INDEPENDENT_SPATIAL_SKY_MODEL_SCENARIO_EXECUTED" in result["limitations"]
    assert editor.draft.digest==draft_before
    assert project.target_radius_arcsec==original_radius
    assert np.array_equal(editor._loaded["data"]["image"],image_before)
    assert not (project.root/"analyses").exists()
    editor._accept_radius_sweep_result(result)
    assert editor._sweep_result is result
    assert "4/5 isophotes" in editor.radius_sweep_status.text()
    editor.copy_scene_qa_btn.click()
    from PySide6.QtWidgets import QApplication
    import json
    export=json.loads(QApplication.clipboard().text())
    assert len(export["protection_radius_sweep"]["scenarios"])==5
    assert export["protection_radius_sweep"]["status"]==(
        "READ_ONLY_DIAGNOSTIC_NOT_SCIENCE_VALIDATED")
    editor.iso_threshold.setValue(2.7)
    assert editor._sweep_result is None

    # The red classifier must be part of the immutable input
    # fingerprint. Without the recorded flag, the GUI would reject
    # all valid red-on sweeps as stale.
    # A component can be completely unchanged while neighbouring masks
    # still differ: expose both facts rather than assuming zero mask change.
    def fake_stable_at_80(dataset, anchor, contrast, smooth, prepared,
                          bridge_guard_pix=0, sga_prior=None):
        radius=float(prepared["stats"]["protected_pixels"])/4
        if radius == 80.:
            stable=dict(prepared)
            stable["stats"]=dict(prepared["stats"],
                                 protected_pixels=70.*4)
            return fake_propose(dataset,anchor,contrast,smooth,stable,
                                bridge_guard_pix=bridge_guard_pix,
                                sga_prior=sga_prior)
        return fake_propose(dataset,anchor,contrast,smooth,prepared,
                            bridge_guard_pix=bridge_guard_pix,
                            sga_prior=sga_prior)
    monkeypatch.setattr(module,"propose_scene_isophote",
                        fake_stable_at_80)
    plateau=module.run_scene_protection_radius_sweep(
        editor._loaded,project.root,editor.draft,editor._project_target,
        threshold_contrast=2.7)
    assert plateau["scenarios"][4]["component_iou_vs_previous_radius"]==1.
    assert plateau["scenarios"][4]["mask_jaccard_vs_previous_radius"] < 1.
    assert plateau["scenarios"][4]["component_iou_vs_60arcsec"] < 1.
    assert plateau["scenarios"][4]["masked_pixel_change_vs_previous_radius"]>0
    assert editor.draft.digest==draft_before
    editor.iso_red_extended_preview.setChecked(True)
    extra=module.run_scene_protection_radius_sweep(
        editor._loaded, project.root, editor.draft,
        editor._project_target, threshold_contrast=2.7,
        red_extended_preview=True)
    assert extra["fixed_settings"]["red_extended_preview"] is True
    assert extra["fixed_settings"]["red_extended_policy"]=="red_extended_v1"
    editor._accept_radius_sweep_result(extra)
    assert editor._sweep_result is extra
    assert editor.draft.digest==draft_before
    assert np.array_equal(editor._loaded["data"]["image"],image_before)


def test_radius_sweep_refuses_invalid_radii_without_changing_project(
        qtbot, tmp_path, monkeypatch):
    import lsbpipe.gui.scene_editor as module
    editor, project=_make_editor(qtbot,tmp_path,monkeypatch)
    project_before=project.target_radius_arcsec
    for radii in ((40.,50.,70.,80.),(40.,60.,60.),
                  (40.,60.,100.),(-1.,60.)):
        with pytest.raises(ValueError,match="Sweep radii"):
            module.run_scene_protection_radius_sweep(
                editor._loaded, project.root, editor.draft,
                editor._project_target, radii_arcsec=radii)
    assert project.target_radius_arcsec==project_before
    assert editor.draft.edits==()
    assert not (project.root/"analyses").exists()



def test_saved_legacy_synthetic_isophote_is_clearly_inactive_in_scene_gui(
        qtbot, tmp_path, monkeypatch):
    """A reopened old SC1 tentative contour must not silently protect pixels."""
    from lsbpipe.photometry.sky_masks import SkyMaskEdit, SkyMaskShape
    editor, project = _make_editor(qtbot, tmp_path, monkeypatch)
    legacy = SkyMaskEdit(
        "sc1_iso_reopened_test", "protect_target",
        SkyMaskShape("circle", (359.999, 60.), radius_arcsec=18.),
        band="g",
        comment="AUTO_ISOPHOTE_PRELIMINARY trial=0 manual_override=1")
    editor._push_edits(editor.draft.edits+(legacy,))
    assert editor.draft.edits[0] == legacy
    assert "legacy isophote · INACTIVE" in editor.edits.item(0).text()
    assert len(editor._paths) == 1
    editor.edits.setCurrentRow(0)
    editor.edit_selected()
    assert editor.purpose.currentData() == "isophote_candidate"
    editor._cancel_pending()
    editor.save_draft()
    saved, _ = read_active_mask_draft(project.root)
    assert saved.edits == (legacy,)
    assert editor.draft.edits[0].purpose == "protect_target"
    # Original on-disk revision preserved: only its *interpretation*
    # changes (legacy AUTO_ISOPHOTE no longer active source protection).
    assert "AUTO_ISOPHOTE_PRELIMINARY" in saved.edits[0].comment



def test_opt_in_gaia_inner_mask_preview_recomputes_isophote_inputs(
        qtbot, tmp_path, monkeypatch):
    """Gaia experimental checkbox is explicit, band-safe and reversible."""
    import lsbpipe.gui.scene_editor as module
    editor, project = _make_editor(qtbot, tmp_path, monkeypatch)
    editor._accept_dataset(_data(project.project_id))
    assert not editor.iso_gaia_inner_preview.isChecked()
    requests = []

    def prepare(dataset, root, draft, anchor, **kwargs):
        active = bool(kwargs.get("gaia_inner_preview", False))
        requests.append(active)
        outside = np.zeros(dataset["pixel_shape"], bool)
        outside[8:10,8:10] = True
        inner = np.zeros(dataset["pixel_shape"], bool)
        if active:
            inner[45:48,45:48] = True
        full = outside | inner
        return {
            "status":"CATALOGUE_CONTAMINATION_MASK_NOT_ADOPTED",
            "mask":full, "auto_mask":full.copy(),
            "gaia_inner_preview_mask":inner,
            "gaia_inner_preview_enabled":active,
            "protected":np.zeros_like(full),
            "hard":dataset["native_hard"].copy(),
            "project_id":project.project_id,
            "draft_digest":draft.digest,
            "dataset_source":dataset["source"],
            "band":dataset["band"],
            "anchor_icrs":anchor,
            "catalog_sha256":"b"*64,
            "stats":{
                "masked_catalog_sources":1+int(active),
                "skipped_target_overlap":1-int(active),
                "gaia_inner_preview_sources":int(active),
                "gaia_inner_preview_pixels":int(inner.sum()),
                "protected_gaia_astrometry_review":1,
                "excluded_pixels":int(full.sum()),
                "native_hard_pixels":1,
            },
            "flags":("EXPLICIT_GAIA_INTERIOR_MASK_PREVIEW_ONLY",)
                    if active else ("PSF_WINGS_NOT_VALIDATED",),
        }

    monkeypatch.setattr(module, "prepare_scene_neighbour_masks", prepare)
    editor.iso_prepare.click()
    qtbot.waitUntil(
        lambda: editor._future is None and editor._contamination is not None,
        timeout=10000)
    assert requests == [False]
    assert not editor._contamination["mask"][46,46]
    editor.iso_gaia_inner_preview.setChecked(True)
    assert editor._contamination is None
    assert editor._isophote_proposal is None
    editor.iso_prepare.click()
    qtbot.waitUntil(
        lambda: editor._future is None and editor._contamination is not None,
        timeout=10000)
    assert requests == [False, True]
    assert editor._contamination["mask"][46,46]
    assert editor._contamination["gaia_inner_preview_mask"][46,46]
    assert "1 inner Gaia masks" in editor.iso_mask_status.text()
    assert "ORANGE" in editor.iso_mask_status.text()
    assert editor._complete_scene_qa_snapshot()["controls"][
        "gaia_inner_mask_preview"] is True
    assert editor.draft.edits == ()
    editor.iso_gaia_inner_preview.setChecked(False)
    assert editor._contamination is None
    editor.iso_prepare.click()
    qtbot.waitUntil(
        lambda: editor._future is None and editor._contamination is not None,
        timeout=10000)
    assert requests == [False, True, False]
    assert not editor._contamination["mask"][46,46]
    assert editor.draft.edits == ()
    # Experiment selection MUST NOT leak to another native band.
    editor.iso_gaia_inner_preview.setChecked(True)
    editor.set_band("r")
    assert not editor.iso_gaia_inner_preview.isChecked()
    assert editor._contamination is None



def test_red_extended_qa_preview_is_opt_in_reversible_and_not_persisted(
        qtbot, tmp_path, monkeypatch):
    """Unadopted high-specificity classifier crosses GUI worker safely."""
    import lsbpipe.gui.scene_editor as module
    editor, project = _make_editor(qtbot, tmp_path, monkeypatch)
    assert not editor.iso_red_extended_preview.isChecked()
    assert editor.source_registry_filter.findData(
        "RED_EXTENDED_MASK_PREVIEW") >= 0
    seen = []

    def prepare(dataset, root, draft, anchor, **kwargs):
        red = bool(kwargs.get("red_extended_preview", False))
        gaia = bool(kwargs.get("gaia_inner_preview", False))
        seen.append((red,gaia))
        outside = np.zeros(dataset["pixel_shape"], bool)
        outside[2:6,2:6] = True
        red_mask = np.zeros_like(outside)
        if red:
            red_mask[35:40,35:40] = True
        result_mask = outside | red_mask
        return {
            "status":"CATALOGUE_CONTAMINATION_MASK_NOT_ADOPTED",
            "mask":result_mask,"auto_mask":result_mask.copy(),
            "red_extended_preview_mask":red_mask,
            "red_extended_preview_enabled":red,
            "gaia_inner_preview_enabled":gaia,
            "gaia_inner_preview_mask":np.zeros_like(outside),
            "protected":np.zeros_like(outside),
            "hard":dataset["native_hard"].copy(),
            "project_id":project.project_id,
            "draft_digest":draft.digest,
            "dataset_source":dataset["source"],
            "band":dataset["band"],
            "anchor_icrs":anchor,
            "catalog_sha256":"a"*64,
            "stats":{
                "masked_catalog_sources":1+int(red),
                "skipped_target_overlap":1-int(red),
                "gaia_inner_preview_sources":0,
                "red_extended_preview_sources":int(red),
                "red_extended_preview_pixels":int(red_mask.sum()),
                "excluded_pixels":int(result_mask.sum()),
                "native_hard_pixels":0},
            "flags":("EXPLICIT_RED_EXTENDED_MASK_PREVIEW_ONLY",)
                    if red else ("PSF_WINGS_NOT_VALIDATED",),
        }

    monkeypatch.setattr(module, "prepare_scene_neighbour_masks", prepare)
    editor.iso_prepare.click()
    qtbot.waitUntil(
        lambda:editor._future is None and editor._contamination is not None,
        timeout=10000)
    assert seen==[(False,False)]
    assert not editor._contamination["mask"][37,37]
    editor.iso_red_extended_preview.setChecked(True)
    assert editor._contamination is None
    editor.iso_prepare.click()
    qtbot.waitUntil(
        lambda:editor._future is None and editor._contamination is not None,
        timeout=10000)
    assert seen==[(False,False),(True,False)]
    assert editor._contamination["mask"][37,37]
    assert editor._contamination["red_extended_preview_mask"][37,37]
    assert "1 red extended masks" in editor.iso_mask_status.text()
    assert "CYAN" in editor.iso_mask_status.text()
    assert editor._complete_scene_qa_snapshot()["controls"][
        "red_extended_mask_preview"] is True
    assert editor.draft.edits==()
    editor.iso_red_extended_preview.setChecked(False)
    assert editor._contamination is None
    editor.iso_prepare.click()
    qtbot.waitUntil(
        lambda:editor._future is None and editor._contamination is not None,
        timeout=10000)
    assert seen[-1]==(False,False)
    assert not editor._contamination["mask"][37,37]
    assert editor.draft.edits==()
    editor.iso_red_extended_preview.setChecked(True)
    editor.set_band("r")
    assert not editor.iso_red_extended_preview.isChecked()



def test_policy_ablation_fixed_radius_exact_pixels_and_starwise_leave_one_out(
        qtbot,tmp_path,monkeypatch):
    """Four native-mask policies plus two independent star interventions.

    The GUI accepts the worker result, and a setting change invalidates it.
    No native image, source FITS, mask draft or project setting is modified.
    """
    import lsbpipe.gui.scene_editor as module
    from lsbpipe.photometry.ellipse_seed import EllipseSeed, LocalIsophotalProposal

    editor,project=_make_editor(qtbot,tmp_path,monkeypatch)
    assert not (project.root/"analyses").exists()
    draft_before=editor.draft.digest
    source_before=editor._loaded["data"]["image"].copy()
    project_radius=project.target_radius_arcsec
    calls=[]
    from types import SimpleNamespace
    radial_bins=({
        "outer_arcsec":63.,"sectors_positive":6,
        "sector_contrast_by_octant":(1.2, None, .9, .8, 1.1, .7, .8, .7),
    },)
    radial=SimpleNamespace(
        detected_a_arcsec=63.,radial_step_arcsec=1.5,
        flags=("CONSTANT_LOCAL_SKY_NOT_VALIDATED",),bins=radial_bins)

    def prepare(dataset,scenario_project,draft,anchor,*,sga_prior=None,
                gaia_inner_preview=False,red_extended_preview=False,
                gaia_inner_skip_objids=()):
        assert scenario_project.target_radius_arcsec==70.
        assert draft.digest==draft_before
        gaia=bool(gaia_inner_preview)
        red=bool(red_extended_preview)
        skipped=set(gaia_inner_skip_objids)
        calls.append((gaia,red,tuple(sorted(skipped))))
        mask=np.zeros(dataset["pixel_shape"],dtype=bool)
        protected=np.zeros_like(mask)
        protected[10:80,10:80]=True
        mask[:2,:2]=True
        if gaia and 11 not in skipped:
            mask[10:14,10:14]=True
        if gaia and 22 not in skipped:
            mask[18:23,18:23]=True
        if red:
            mask[29:33,29:33]=True
        dispositions=[]
        for oid in (11,22):
            dispositions.append({
                "source_id":f"dr11/south/1554p180/OBJID={oid}",
                "screening_disposition":(
                    "GAIA_FOREGROUND_MASK_PREVIEW"
                    if gaia and oid not in skipped
                    else "PROTECTED_OVERLAP")})
        dispositions.append({
            "source_id":"dr11/south/1554p180/OBJID=1794",
            "screening_disposition":(
                "RED_EXTENDED_MASK_PREVIEW"
                if red else "PROTECTED_OVERLAP")})
        mask.setflags(write=False)
        return {
            "status":"CATALOGUE_CONTAMINATION_MASK_NOT_ADOPTED",
            "project_id":project.project_id,"band":"g",
            "dataset_source":dataset["source"],"draft_digest":draft.digest,
            "anchor_icrs":anchor,"mask":mask,"protected":protected,
            "catalog_sha256":"a"*64,
            "source_inventory":{
                "total_relevant_sources":3,"truncated":False,
                "screening_evidence_incomplete":False,
                "entries":tuple(dispositions)},
            "stats":{
                "masked_catalog_sources":1+int(red)+
                    (2-len(skipped) if gaia else 0),
                "skipped_target_overlap":3-int(red)-
                    (2-len(skipped) if gaia else 0),
                "gaia_inner_preview_sources":(
                    2-len(skipped) if gaia else 0),
                "gaia_inner_preview_pixels":(
                    (16 if 11 not in skipped else 0)+
                    (25 if 22 not in skipped else 0)
                    if gaia else 0),
                "red_extended_preview_sources":int(red),
                "red_extended_preview_pixels":16 if red else 0,
                "native_hard_pixels":1,
            },
        }

    def propose(dataset,anchor,contrast,smooth,prepared,
                bridge_guard_pix=0,sga_prior=None):
        m=prepared["mask"]
        first=bool(m[11,11])
        second=bool(m[20,20])
        red=bool(m[30,30])
        n=60-8*int(first)-12*int(second)-2*int(red)
        if first and not second and red:
            raise ValueError("Synthetic contour failure for one intervention")
        seed=EllipseSeed(
            x_pix=48.,y_pix=48.,sma_pix=12.,q=.6,
            theta_native_deg=30.,component_pixels=n,
            centroid_offset_pix=float(60-n)/2)
        model=LocalIsophotalProposal(
            seed=seed,contour_xy=((40.,40.),(55.,40.),(49.,55.)),
            local_sky=2.e-5,raw_sky_rms=.002,
            smoothed_sky_rms=.0002,sky_pixels=800,
            threshold_contrast=contrast,smoothing_pix=smooth,
            flags=("BACKGROUND_CONSTANT_NOT_VALIDATED",),
            component_indices=np.arange(n,dtype=np.int64))
        return {"proposal":model,"sga_radial":radial,
                "sensitivity":{
                    "minimum_spatial_iou":.72,
                    "confidence":"LOW_CONFIDENCE_REQUIRES_EXPLICIT_OVERRIDE",
                    "flags":("UNANCHORED_CENTROID_DRIFT",)},
                "sga_warning":None}

    monkeypatch.setattr(module,"prepare_scene_contamination",prepare)
    monkeypatch.setattr(module,"propose_scene_isophote",propose)
    report=module.run_scene_mask_policy_ablation(
        editor._loaded,project.root,editor.draft,editor._project_target)
    assert report["schema"]=="desidera.gui.mask-policy-ablation/1"
    assert report["status"]=="READ_ONLY_DIAGNOSTIC_NOT_SCIENCE_VALIDATED"
    assert report["gaia_candidate_objids"]==[11,22]
    assert report["gaia_star_leave_one_out_count"]==2
    assert len(report["scenarios"])==6
    assert calls==[(False,False,()),(True,False,()),
                   (False,True,()),(True,True,()),
                   (True,True,(11,)),(True,True,(22,))]
    entries={v["policy"]:v for v in report["scenarios"]}
    assert entries["none"]["component_pixels"]==60
    assert entries["gaia"]["component_pixels"]==40
    assert entries["red"]["component_pixels"]==58
    assert entries["both"]["component_pixels"]==38
    assert entries["both"]["sga_radial_bins"]==[{
        "outer_arcsec":63.,"sectors_positive":6,
        "sector_contrast_by_octant":[1.2,None,.9,.8,1.1,.7,.8,.7]}]
    assert entries["both"]["sga_radial_flags"]==[
        "CONSTANT_LOCAL_SKY_NOT_VALIDATED"]
    assert entries["both_without_gaia_22"]["sga_radial_bins"]==[]
    assert entries["both"]["vs_no_experimental_masks"][
        "component_iou"]==pytest.approx(38/60.)
    assert entries["both"]["vs_no_experimental_masks"][
        "mask_xor_pixels"]==57
    # Component-index sampling on the flattened grid is exact; the
    # synthetic baseline component occupies only the first image row.
    assert entries["both"]["masked_pixels_over_baseline_component"]==2
    assert entries["both_without_gaia_11"]["component_pixels"]==46
    assert entries["both_without_gaia_11"][
        "vs_both_experimental_masks"]["component_iou"]==pytest.approx(38/46.)
    assert entries["both_without_gaia_22"]["component_status"]==(
        "ISOPHOTE_NOT_FOUND")
    assert entries["both_without_gaia_22"][
        "vs_both_experimental_masks"]["component_iou"] is None
    assert "Synthetic contour failure" in entries[
        "both_without_gaia_22"]["component_error"]
    assert report["fixed_settings"]["protection_radius_arcsec"]==70.
    assert project.target_radius_arcsec==project_radius
    assert editor.draft.edits==()
    assert editor.draft.digest==draft_before
    assert np.array_equal(editor._loaded["data"]["image"],source_before)
    editor._accept_mask_policy_ablation_result(report)
    assert editor._policy_ablation_result is report
    qa=editor._complete_scene_qa_snapshot()
    # Snapshot normalizes nested diagnostic data for JSON: equivalent
    # content is required, not shared mutable dictionary identity.
    assert qa["mask_policy_ablation"] == report
    assert qa["mask_policy_ablation"] is not report
    assert "READ-ONLY 5/6" in editor.policy_ablation_status.text()
    assert "NOT galaxy-light loss" in editor.policy_ablation_status.text()
    # Regression for a real Leo P user exporting null QA during a running
    # ablation: even a completed Future is not an accepted result until
    # the GUI has dispatched it on the main Qt thread.
    from concurrent.futures import Future
    pending=Future()
    editor._future=pending
    editor._work_kind="policy_ablation"
    editor._refresh_controls()
    assert not editor.copy_scene_qa_btn.isEnabled()
    assert "not been accepted" in editor.copy_scene_qa_btn.toolTip()
    pending.set_result(report)
    editor._refresh_controls()
    assert not editor.copy_scene_qa_btn.isEnabled()
    editor._future=None
    editor._work_kind=None
    editor._refresh_controls()
    assert editor.copy_scene_qa_btn.isEnabled()
    # Exercise the actual Qt-worker dispatch as well as pure backend.
    editor.policy_ablation_btn.click()
    qtbot.waitUntil(
        lambda:editor._future is None and
        editor._policy_ablation_result is not None,
        timeout=15000)
    assert len(calls)==12
    assert "READ-ONLY 5/6" in editor.policy_ablation_status.text()
    # Experimental checkbox state is NOT a policy-ablation input:
    # the worker computes all four policies regardless of the UI toggle.
    editor.iso_gaia_inner_preview.setChecked(True)
    assert editor._policy_ablation_result is None
    assert editor.draft.edits==()
    assert not (project.root/"analyses").exists()
    with pytest.raises(ValueError,match="Ablation radius"):
        module.run_scene_mask_policy_ablation(
            editor._loaded,project.root,editor.draft,editor._project_target,
            protection_radius_arcsec=1e6)



def test_psf0_foreground_probe_is_explicit_read_only_and_scene_qa_bound(
        qtbot,tmp_path,monkeypatch):
    import lsbpipe.gui.scene_editor as module
    editor,project=_make_editor(qtbot,tmp_path,monkeypatch)
    assert editor.psf0_fg_btn.objectName()=="scene_psf0_foreground_probe"
    assert not editor.psf0_fg_btn.isEnabled()
    candidate={
        "source_id":"dr11/south/1554p600/OBJID=2035",
        "objid":2035, "x":48., "y":48., "tractor_type":"PSF",
        "gaia_dr3_source_id":623919492384922624,
        "foreground_evidence_assessment":
            "GAIA_ASTROMETRY_FOREGROUND_CANDIDATE_REVIEW",
    }
    editor._source_registry={
        "catalog_sha256":"a"*64,"entries":[candidate]}
    zero=np.zeros(editor._loaded["pixel_shape"],dtype=bool)
    editor._contamination={
        "catalog_sha256":"a"*64,"draft_digest":editor.draft.digest,
        "mask":zero,"auto_mask":zero.copy()}
    editor._source_registry_visible=[candidate]
    editor.source_registry_list.blockSignals(True)
    editor.source_registry_list.addItem("OBJID 2035")
    editor.source_registry_list.setCurrentRow(0)
    editor.source_registry_list.blockSignals(False)
    calls=[]
    def fake_probe(dataset,root,source,inventory,centre,indices):
        calls.append((dataset["band"],source["objid"],centre,indices))
        return {
            "status":"FINITE_PSF0_LOCAL_PROBE_DIAGNOSTIC_NOT_ADOPTED",
            "project_id":project.project_id,
            "source_id":candidate["source_id"],
            "source_objid":2035,"dataset_source":dataset["source"],
            "band":"g","tractor_catalog_sha256":"a"*64,
            "measured_fwhm_native_pix":5.137,
            "star_to_target_centre_native_pix":76.,
            "finite_kernel_half_width_native_pix":40,
            "galaxy_centre_outside_finite_kernel_square":True,
            "sky_fit_comparison":{
                "constant":{
                    "heldout_weighted_rms_empirical_not_significance":1.2},
                "planar":{
                    "heldout_weighted_rms_empirical_not_significance":1.1}},
            "science_adopted":False,"source_light_subtracted":False,
            "original_pixels_modified":False,
            "project_only_outer_sector_probe":{
                "schema":"desidera.scene.project-only-foreground-sector/1",
                "status":"SAME_PROJECT_DIRECTIONAL_SKY_SENSITIVITY_ONLY",
                "additional_bricks_or_catalogues_used":False,
                "science_adopted":False,
            },
        }
    monkeypatch.setattr(module,"run_scene_psf0_foreground_probe",fake_probe)
    editor._refresh_controls()
    assert editor.psf0_fg_btn.isEnabled()
    editor.psf0_fg_btn.click()
    assert editor._work_kind=="psf_foreground"
    assert not editor.copy_scene_qa_btn.isEnabled()
    qtbot.waitUntil(lambda: editor._future is None,timeout=5000)
    assert calls==[("g",2035,(359.999,60.),None)]
    assert editor._foreground_psf_result["source_objid"]==2035
    assert "UNKNOWN" in editor.psf0_fg_status.text()
    assert "Same-project direction/sky control" in editor.psf0_fg_status.text()
    snap=editor._complete_scene_qa_snapshot()
    assert snap["psf0_foreground_probe"]["source_id"]==candidate["source_id"]
    assert snap["psf0_foreground_probe"]["project_only_outer_sector_probe"][
        "additional_bricks_or_catalogues_used"] is False
    assert snap["psf0_foreground_probe"]["science_adopted"] is False
    editor.set_band("r")
    assert editor._foreground_psf_result is None
    assert not editor.psf0_fg_btn.isEnabled()
