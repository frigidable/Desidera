"""Leo P is the FIRST real-object control; test pilot never invents its data."""
from pathlib import Path

import numpy as np
import pytest
pytest.importorskip("astropy")
pytest.importorskip("photutils")
from astropy.io import fits
from astropy.wcs import WCS

from lsbpipe.surveys.object_projects import (
    create_object_project, save_object_project)
from lsbpipe.surveys.legacy import SkyField
from lsbpipe.surveys.dataset_integrity import _expected_paths
from lsbpipe.photometry.psf_empirical import PSF0Config
from scripts.psf0_leo_p import (
    LEO_P_WORKING_ICRS, inspect_leo_p_inputs, run_leo_p_pilot)


def _project(tmp_path, *, name="Leo P", ra=155.438012, dec=18.088025):
    p = create_object_project(
        name, tmp_path / name.replace(" ","_"),
        field=SkyField(ra, dec, 6.))
    p.selected_dataset = ("dr11", "south", "1554p180", "r")
    p.data_root = tmp_path / "originals"
    p.data_root.mkdir(exist_ok=True)
    save_object_project(p)
    return p


def test_leo_p_pilot_is_read_only_without_existing_fits(tmp_path, monkeypatch):
    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path / "session"))
    p = _project(tmp_path)
    before = (p.root / "lsbpipe-object.json").read_bytes()
    info = inspect_leo_p_inputs(p.root)
    assert info["working_center_icrs_deg"] == list(LEO_P_WORKING_ICRS)
    assert info["dataset"] == {
        "release":"dr11", "region":"south", "brick":"1554p180",
        "selected_band":"r"}
    assert [b["band"] for b in info["bands"]] == ["r","g","z"]
    assert all(b["state"] == "MISSING_OR_UNSAFE_NATIVE_FILES"
               for b in info["bands"])
    assert not (p.root / "analyses").exists()
    assert (p.root / "lsbpipe-object.json").read_bytes() == before
    result = run_leo_p_pilot(p.root)
    assert all(b["status"] == "SKIPPED_MISSING_NATIVE_QUARTET"
               for b in result["band_results"])
    assert not (p.root / "analyses").exists()


def test_leo_p_rejects_other_objects_and_inconsistent_centers(tmp_path, monkeypatch):
    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path / "session"))
    with pytest.raises(ValueError, match="Leo P"):
        inspect_leo_p_inputs(_project(tmp_path, name="UGC 2162").root)
    p = _project(tmp_path / "shifted", dec=18.18)
    with pytest.raises(ValueError, match="differs"):
        inspect_leo_p_inputs(p.root)


def _synthetic_legacy_files(project):
    rng = np.random.default_rng(173)
    shape = (280,280)
    y,x=np.indices(shape,dtype=float)
    image=1.2+rng.normal(0.,.28,shape)
    for j,(cx,cy) in enumerate((cx,cy) for cy in (35,78,122,166,209,248)
                                      for cx in (35,78,122,166,209,248)):
        image += 145.*np.exp(-.5*((x-cx-.12)**2+(y-cy+.29)**2)/2.12**2)
    w = WCS(naxis=2)
    w.wcs.crpix=[140.5,140.5]
    w.wcs.crval=list(LEO_P_WORKING_ICRS)
    w.wcs.cdelt=[-.262/3600., .262/3600.]
    w.wcs.ctype=["RA---TAN","DEC--TAN"]
    paths=_expected_paths(project.data_root, project.selected_dataset)
    for kind,path in paths.items():
        path.parent.mkdir(parents=True,exist_ok=True)
        h=w.to_header()
        if kind=="image":
            h["BUNIT"]="nanomaggies"
            fits.PrimaryHDU(data=image.astype(np.float32),header=h).writeto(path)
        elif kind=="invvar":
            h["BUNIT"]="1/nanomaggies^2"
            fits.PrimaryHDU(data=np.ones(shape,dtype=np.float32)*10.,
                            header=h).writeto(path)
        elif kind=="nexp":
            fits.PrimaryHDU(data=np.ones(shape,dtype=np.int16)*2,
                            header=h).writeto(path)
        else:
            fits.HDUList([fits.PrimaryHDU(),
                          fits.ImageHDU(data=np.zeros(shape,dtype=np.int16),
                                        header=h)]).writeto(path)
    return paths


def test_synthetic_native_legacy_fits_runs_leo_p_pilot_without_false_science(
        tmp_path, monkeypatch):
    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path / "session"))
    p = _project(tmp_path)
    paths = _synthetic_legacy_files(p)
    originals = {k:path.read_bytes() for k,path in paths.items()}
    recipe = PSF0Config(
        max_training_stars=30, min_train_stars=6, min_holdout_stars=2)
    outcome = run_leo_p_pilot(
        p.root, bands=("r", "g"), config=recipe,
        protection_arcsec=12., save_diagnostic=True)
    r,g=outcome["band_results"]
    assert r["status"] == "EXPERIMENTAL_PSF_NOT_ADOPTED", r
    assert r["report"]["lineage_id"]
    assert r["report"]["science_adopted"] is False
    assert g["status"] == "SKIPPED_MISSING_NATIVE_QUARTET"
    folder=Path(r["report"]["diagnostic_saved_to"])
    assert folder.is_dir()
    assert (folder/"psf0_report.json").exists()
    assert (folder/"empirical_psf.fits").exists()
    with fits.open(folder/"empirical_psf.fits") as hdul:
        assert np.isclose(hdul[0].data.sum(),1.,atol=1e-8)
        assert hdul[0].header["PSFVALID"] is False
        assert hdul[0].header["BAND"]=="r"
    assert all(paths[k].read_bytes()==data for k,data in originals.items())
    again=run_leo_p_pilot(
        p.root,bands=("r",), config=recipe,
        protection_arcsec=12.,save_diagnostic=True)
    assert again["band_results"][0]["status"]=="PSF0_FAILED_NOT_ADOPTED"
    assert "already exists" in again["band_results"][0]["reason"]
    assert (p.root/"lsbpipe-object.json").is_file()


def test_gui_preview_refuses_export_after_original_fits_changed(tmp_path, monkeypatch):
    from lsbpipe.photometry.psf_gui_backend import run_gui_psf0, save_gui_psf0

    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path / "session"))
    project = _project(tmp_path)
    paths = _synthetic_legacy_files(project)
    preview = run_gui_psf0(
        project.root, band="r",
        config=PSF0Config(max_training_stars=30,
                          min_train_stars=6, min_holdout_stars=2),
        protection_arcsec=12.)
    assert preview["status"] == "EXPERIMENTAL_PSF_NOT_ADOPTED"
    assert preview["result"].lineage_id is not None
    assert not (project.root/"analyses"/"psf-pilots").exists()
    # Intentional change after successful preview invalidates source SHA.
    with paths["image"].open("ab") as stream:
        stream.write(b"extra bytes after original preview")
    with pytest.raises(ValueError, match="changed since PSF build"):
        save_gui_psf0(preview)
    assert not (project.root/"analyses"/"psf-pilots").exists()


def test_explicit_manual_psf_star_exclusions_persist_without_modifying_fits(
        tmp_path, monkeypatch):
    from dataclasses import replace
    import json
    from lsbpipe.photometry.psf_gui_backend import (
        run_gui_psf0, commit_psf0_manual_review)
    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path/"session"))
    p = _project(tmp_path)
    originals = _synthetic_legacy_files(p)
    contents = {kind:path.read_bytes() for kind,path in originals.items()}
    config = PSF0Config(
        max_training_stars=30, min_train_stars=6, min_holdout_stars=2)
    first = run_gui_psf0(
        p.root, band="r", config=config, protection_arcsec=12.)
    assert first["status"] == "EXPERIMENTAL_PSF_NOT_ADOPTED"
    selected = first["star_review"]["stars"][0]["star_id"]
    reviewed = run_gui_psf0(
        p.root, band="r",
        config=replace(config, excluded_training_star_ids=(selected,)),
        protection_arcsec=12.)
    assert reviewed["status"] == "EXPERIMENTAL_PSF_NOT_ADOPTED"
    assert len(reviewed["result"].training) == len(first["result"].training)-1
    assert {s.star_id for s in reviewed["result"].holdout} == {
        s.star_id for s in first["result"].holdout}
    with pytest.raises(ValueError, match="explicitly applied"):
        commit_psf0_manual_review(reviewed)
    assert not (p.root/"analyses").exists()
    receipt_path = Path(commit_psf0_manual_review(
        reviewed, explicitly_confirmed=True))
    assert receipt_path.is_file()
    receipt = json.loads(receipt_path.read_text("utf-8"))
    assert receipt["status"] == "MANUAL_SELECTION_NOT_PSF_ADOPTION"
    assert receipt["science_adopted"] is False
    assert receipt["excluded_training_star_ids"] == [selected]
    assert receipt["holdout_protected"] is True
    assert all(originals[kind].read_bytes()==data
               for kind, data in contents.items())
    again = commit_psf0_manual_review(reviewed, explicitly_confirmed=True)
    assert again == str(receipt_path)
    assert list(receipt_path.parent.glob("review-*.json")) == [receipt_path]
    recovered = run_gui_psf0(
        p.root, band="r", config=config, protection_arcsec=12.)
    assert recovered["saved_manual_review"]["excluded_training_star_ids"] == (
        selected,)
    assert recovered["saved_manual_review"]["auto_applied"] is False
    assert {star.star_id for star in recovered["result"].training} == {
        star.star_id for star in first["result"].training}
    with originals["image"].open("ab") as changed:
        changed.write(b"external source modification")
    with pytest.raises(ValueError, match="Original FITS changed"):
        commit_psf0_manual_review(reviewed, explicitly_confirmed=True)



def test_gui_r_band_lineage_survives_project_default_band_switch(
        tmp_path, monkeypatch):
    """End-to-end native synthetic FITS: analyzed r != saved default g/i."""
    from lsbpipe.photometry.psf_gui_backend import (
        run_gui_psf0, save_gui_psf0, commit_psf0_manual_review)

    monkeypatch.setenv("LSBPIPE_SESSION_DIR", str(tmp_path/"session"))
    project = _project(tmp_path)
    paths = _synthetic_legacy_files(project)
    originals = {kind: p.read_bytes() for kind, p in paths.items()}
    project.selected_dataset = ("dr11", "south", "1554p180", "g")
    save_object_project(project)
    preview = run_gui_psf0(
        project.root, band="r",
        config=PSF0Config(
            max_training_stars=30, min_train_stars=6, min_holdout_stars=2),
        protection_arcsec=12.)
    assert preview["status"] == "EXPERIMENTAL_PSF_NOT_ADOPTED"
    assert preview["dataset"] == ("dr11", "south", "1554p180", "r")
    assert preview["result"].dataset_id.endswith("/r")

    # Another default filter is just project navigation, not a change to
    # the original r FITS. Both immutable user-initiated derivatives remain
    # tied to the actual native r science source.
    project.selected_dataset = ("dr11", "south", "1554p180", "i")
    save_object_project(project)
    receipt = Path(commit_psf0_manual_review(
        preview, explicitly_confirmed=True))
    assert receipt.is_file()
    output = Path(save_gui_psf0(preview))
    assert (output/"empirical_psf.fits").is_file()
    with fits.open(output/"empirical_psf.fits") as hdul:
        assert hdul[0].header["BAND"] == "r"
        assert hdul[0].header["PSFVALID"] is False
    assert all(paths[k].read_bytes() == original
               for k, original in originals.items())

    project.selected_dataset = ("dr11", "south", "1555p180", "r")
    save_object_project(project)
    with pytest.raises(ValueError, match="selected dataset changed"):
        save_gui_psf0(preview)
