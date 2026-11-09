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
from lsbpipe.surveys.object_projects import create_object_project


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
    assert view.stage_tabs.count()==2
    assert view.stage_tabs.tabText(0)=="1 · PSF-0"
    assert view.stage_tabs.tabText(1)=="2 · Scene & Masks"
    assert view.scene_editor._project_id==p.project_id
    view.stage_tabs.setCurrentWidget(view.scene_editor)
    assert calls==[]
    assert view.scene_editor._loaded is None
    assert not (p.root/"analyses").exists()
    assert view.scene_editor.draft.edits==()
    assert "NOT validated" in view.scene_editor.warning.text()


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
    assert edit.purpose == "protect_target"
    assert edit.band == "g"  # native-band default for source-derived contours
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
