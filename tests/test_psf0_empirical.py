"""PSF-0 synthetic, source-isolated regression: no fabricated Leo P measurement."""
import json
import numpy as np
import pytest

pytest.importorskip("astropy")
pytest.importorskip("photutils")
from astropy.wcs import WCS

from lsbpipe.photometry.contracts import NativeBandFrame
from lsbpipe.photometry.psf_empirical import (
    PSF0Config, STATUS, build_empirical_coadd_psf)
from lsbpipe.photometry.psf_export import export_psf0_diagnostic


LEO_P = (155.438012, 18.088025)


def pilot_wcs(shape=(280,280), *, band="r"):
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [shape[1]/2.+.5, shape[0]/2.+.5]
    wcs.wcs.crval = list(LEO_P)
    wcs.wcs.cd = np.array([[-.262,0],[0,.262]])/3600.
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    return wcs


def fake_frame(*, band="r", nstars="many", amplitude=145., bad=False,
               instrument_region="south"):
    rng = np.random.default_rng(173)
    shape = (280,280)
    yy,xx = np.indices(shape, dtype=float)
    img = 1.2 + 0.00013*xx + rng.normal(0, .28, shape)
    sites = [(x, y) for y in (35,78,122,166,209,248)
             for x in (35,78,122,166,209,248)]
    if nstars == "few":
        sites = sites[:5]
    if nstars == "none":
        sites = []
    for j,(x0,y0) in enumerate(sites):
        # Deliberately noninteger subpixel centroids and variable stellar flux.
        x0 += .12 if j%2 else -.24
        y0 += .29 if j%3 else -.19
        flux = amplitude*(1+0.15*(j%3))
        img += flux*np.exp(-.5*((xx-x0)**2+(yy-y0)**2) / (2.12**2))
    # A galaxy-like feature near the ICRS target must be excluded
    # from stellar PSF training even if its knots are bright.
    img += 21.*np.exp(-.5*((xx-140)**2+(yy-140)**2)/ (8.**2))
    ivar = np.ones(shape)*1./(.28**2)
    if bad:
        ivar[35,35] = 0.
    return NativeBandFrame(
        image_nMgy=img, invvar_nMgy2=ivar,
        nexp=np.ones(shape)*2.,
        maskbits=np.zeros(shape, dtype=np.int32),
        native_wcs=pilot_wcs(shape),
        release="dr11", region=instrument_region,
        brick="1554p180", band=band,
        source_evidence={"synthetic":True,
                         "data_status":"SYNTHETIC_NOT_OBSERVATIONAL_FITS"})


@pytest.fixture(scope="module")
def built():
    return build_empirical_coadd_psf(
        fake_frame(), target_icrs_deg=LEO_P,
        target_protection_arcsec=13.,
        config=PSF0Config(min_neighbor_separation_pix=32., min_train_stars=6,
                          min_holdout_stars=2, max_training_stars=30))


def test_returns_normalized_native_grid_epsf_with_disjoint_holdouts(built):
    result = built
    assert result.status == STATUS
    assert result.kernel.shape == (31,31)
    assert result.kernel.sum() == pytest.approx(1.,abs=1e-10)
    assert result.kernel.min() >= 0
    assert not result.kernel.flags.writeable
    with pytest.raises(ValueError):
        result.kernel.setflags(write=True)
    assert len(result.training) >= 6
    assert len(result.holdout) >= 2
    assert {s.star_id for s in result.training}.isdisjoint(
        {s.star_id for s in result.holdout})
    assert result.metrics["holdout_median_abs_core_residual_over_peak"] < .12
    assert result.lineage_id is None
    assert "SYNTHETIC_TEST_INPUT" in result.flags
    assert "CORE_PSF_ONLY_EXTENDED_WINGS_UNVALIDATED" in result.flags
    assert "SPATIAL_REPRESENTATIVITY_NOT_VALIDATED" in result.flags
    assert result.instrument == "DECAM"


def test_approximate_core_profile_and_no_galaxy_centroids(built):
    kernel = built.kernel
    c = kernel.shape[0]//2
    # A reasonably resolved Gaussian core should dominate outer pixels.
    assert kernel[c,c] > 8.*kernel[c+7,c]
    assert kernel[c,c] > kernel[0,0]
    assert all(np.hypot(star.x_pix-140,star.y_pix-140) > 13.
               for star in built.training+built.holdout)
    # No psf wing accuracy or galaxy measurement is claimed.
    report = built.report()
    assert report["science_adopted"] is False
    assert report["psf_wings_validated"] is False
    assert report["noise_covariance_validated"] is False
    assert report["n_train"] == len(built.training)
    assert report["n_holdout"] == len(built.holdout)
    assert report["kernel_sha256"]
    json.dumps(report, allow_nan=False)


def test_other_band_produces_own_identity_not_cloned_kernel(built):
    g = build_empirical_coadd_psf(
        fake_frame(band="g"), target_icrs_deg=LEO_P,
        target_protection_arcsec=13.)
    assert g.dataset_id.endswith("/g") and built.dataset_id.endswith("/r")
    assert g.lineage_id is None
    assert g.report()["status"] == STATUS


def test_no_data_or_few_sources_blocks_psf():
    for population in ("few", "none"):
        with pytest.raises(ValueError, match="star|stellar|source|candidates|detected"):
            build_empirical_coadd_psf(
                fake_frame(nstars=population), target_icrs_deg=LEO_P,
                target_protection_arcsec=13.)


def test_source_quality_rejects_invalid_and_protected_pixels():
    frame = fake_frame(bad=True)
    result = build_empirical_coadd_psf(
        frame, target_icrs_deg=LEO_P, target_protection_arcsec=13.)
    assert not frame.valid[35,35]
    assert result.status == STATUS
    assert not any(abs(star.x_pix-35)<2 and abs(star.y_pix-35)<2
                   for star in result.training+result.holdout)


def test_synthetic_psf_must_not_be_exported_as_observational_fits(tmp_path,built):
    with pytest.raises(ValueError, match="Explicit researcher"):
        export_psf0_diagnostic(built, tmp_path)
    with pytest.raises(ValueError, match="SHA-256-verified"):
        export_psf0_diagnostic(
            built, tmp_path, explicitly_requested=True)
    assert not list(tmp_path.glob("*"))


def test_invalid_recipe_and_wrong_input_fail_closed():
    with pytest.raises(ValueError, match="oversampling"):
        PSF0Config(psf_oversampling=2)
    with pytest.raises(ValueError, match="Star cutout"):
        PSF0Config(cutout_size_pix=30)
    with pytest.raises(TypeError, match="NativeBandFrame"):
        build_empirical_coadd_psf(
            np.zeros((30,30)), target_icrs_deg=LEO_P,
            target_protection_arcsec=90.)
    with pytest.raises(ValueError, match="Invalid ICRS"):
        build_empirical_coadd_psf(
            fake_frame(), target_icrs_deg=(410.,18.),
            target_protection_arcsec=90.)


def test_rejected_epsf_retains_star_diagnostics_without_export(monkeypatch):
    """A failed ePSF must not be quietly clipped or accepted for Leo P."""
    from types import SimpleNamespace
    import photutils.psf
    from lsbpipe.photometry.psf_empirical import PSF0QualityError

    class RingingBuilder:
        def __init__(self, **kwargs):
            self.shape = kwargs["shape"]

        def __call__(self, stars):
            kernel = np.zeros(self.shape, dtype=float)
            c = self.shape[0] // 2
            kernel[c, c] = 1.
            kernel[c+5, c] = -.2
            return SimpleNamespace(
                epsf=SimpleNamespace(data=kernel), converged=True)

    monkeypatch.setattr(photutils.psf, "EPSFBuilder", RingingBuilder)
    with pytest.raises(PSF0QualityError) as raised:
        build_empirical_coadd_psf(
            fake_frame(), target_icrs_deg=LEO_P,
            target_protection_arcsec=13., config=PSF0Config(
                min_train_stars=6, min_holdout_stars=2,
                max_training_stars=30))
    diagnostic = raised.value.diagnostics
    assert diagnostic["status"] == "PSF0_QUALITY_REJECTED_NOT_ADOPTED"
    assert diagnostic["negative_absolute_mass_fraction"] > .025
    assert diagnostic["n_train"] >= 6
    assert diagnostic["n_holdout"] >= 2
    assert diagnostic["training"]
    assert diagnostic["holdout"]
    assert len(diagnostic["training"]) == diagnostic["n_train"]
    assert "STRONG_NEGATIVE_EPSF_RINGING" == diagnostic["reason_code"]


def test_photutils31_wing_smoothing_ab_is_executed_and_audited():
    import inspect
    from photutils.psf import EPSFBuilder
    from lsbpipe.photometry.psf_empirical import PSF0QualityError
    if "wing_smoothing" not in inspect.signature(EPSFBuilder).parameters:
        pytest.skip("photutils <3.1 lacks native wing_smoothing option")
    config = PSF0Config(
        wing_smoothing=False, min_train_stars=6,
        min_holdout_stars=2, max_training_stars=30)
    try:
        proposal = build_empirical_coadd_psf(
            fake_frame(), target_icrs_deg=LEO_P,
            target_protection_arcsec=13., config=config)
    except PSF0QualityError as exc:
        diagnostic = exc.diagnostics
        assert diagnostic["status"] == "PSF0_QUALITY_REJECTED_NOT_ADOPTED"
        solver = diagnostic["solver_diagnostics"]
        assert diagnostic["psf_recipe"]["wing_smoothing"] is False
    else:
        assert proposal.status == STATUS
        solver = proposal.metrics["solver_diagnostics"]
        assert proposal.recipe.wing_smoothing is False
    assert solver["wing_smoothing_requested"] is False
    assert solver["wing_smoothing_explicitly_set"] is True
    assert solver["max_iterations_requested"] == 8
    assert solver["photutils_version"]
    history = solver["iteration_history_if_reported"]
    assert isinstance(history, list)
    assert history, "photutils 3.1 must expose actual ePSF iteration evidence"
    assert all("iteration" in step for step in history)
    assert len(history) <= 25
    assert "converged_if_reported" in solver


def test_manual_training_star_review_recomputes_without_touching_holdouts():
    from dataclasses import replace
    from lsbpipe.photometry.psf_empirical import PSF0QualityError

    frame = fake_frame()
    base = PSF0Config(
        min_train_stars=6, min_holdout_stars=2,
        max_training_stars=30, min_neighbor_separation_pix=32.)
    def run(cfg):
        try:
            outcome = build_empirical_coadd_psf(
                frame, target_icrs_deg=LEO_P,
                target_protection_arcsec=13., config=cfg)
            return outcome.training, outcome.holdout, outcome.review_snapshot
        except PSF0QualityError as exc:
            from types import SimpleNamespace
            train = [SimpleNamespace(**star) for star in
                     exc.diagnostics["training"]]
            holdout = [SimpleNamespace(**star) for star in
                       exc.diagnostics["holdout"]]
            return train, holdout, exc.review_snapshot

    train0, holdout0, review0 = run(base)
    assert review0["schema"] == "desidera.psf0-star-review-in-memory/1"
    assert len(review0["stars"]) == len(train0)
    assert all(s["cutout"].shape == (31, 31) for s in review0["stars"])
    assert all(s["model"].shape == s["cutout"].shape
               for s in review0["stars"])
    to_remove = train0[0].star_id
    excluded_cfg = replace(base, excluded_training_star_ids=(to_remove,))
    train1, holdout1, review1 = run(excluded_cfg)
    assert len(train1) == len(train0)-1
    assert {s.star_id for s in holdout0} == {s.star_id for s in holdout1}
    assert to_remove not in {s.star_id for s in train1}
    assert {s["star_id"] for s in review1["stars"]} == {
        s["star_id"] for s in review0["stars"]}
    assert next(s for s in review1["stars"]
                if s["star_id"] == to_remove)["role"] == "excluded"
    assert base.digest != excluded_cfg.digest
    with pytest.raises(ValueError, match="training sample|protected holdouts"):
        build_empirical_coadd_psf(
            frame, target_icrs_deg=LEO_P, target_protection_arcsec=13.,
            config=replace(base, excluded_training_star_ids=(
                holdout0[0].star_id,)))


def test_manual_exclusion_recipe_rejects_malformed_or_duplicated_ids():
    from dataclasses import replace
    with pytest.raises(ValueError, match="sorted, unique"):
        PSF0Config(excluded_training_star_ids=("star-0001", "star-0001"))
    with pytest.raises(ValueError, match="sorted, unique"):
        PSF0Config(excluded_training_star_ids=("star-2",))
    with pytest.raises(ValueError, match="sorted, unique"):
        PSF0Config(excluded_training_star_ids=("star-0003", "star-0002"))


def test_gaia_strict_requires_catalog_and_preserves_holdout_when_retraining(built):
    from dataclasses import replace
    from lsbpipe.photometry.psf_tractor_gaia import TractorGaiaMatcher
    from lsbpipe.photometry.psf_empirical import PSF0QualityError

    frame = fake_frame()
    base = PSF0Config(
        min_neighbor_separation_pix=32., min_train_stars=6,
        min_holdout_stars=2,max_training_stars=30,
        star_catalog_policy="gaia_strict")
    with pytest.raises(ValueError,match="GAIA_STRICT_REQUIRES_NATIVE_DR11_TRACTOR"):
        build_empirical_coadd_psf(
            frame,target_icrs_deg=LEO_P,target_protection_arcsec=13.,
            config=base)

    references = list(built.training) + list(built.holdout)
    rows = tuple({
        "ref_cat":"G3","ref_id":int(1000+i),
        "objid":i,"type":"PSF","fitbits":0,
        "maskbits":0,"gaia_duplicated_source":False,
        "flux":4.,"flux_ivar":100.,"gaia_g_mag":18.,
    } for i,_ in enumerate(references))
    positions=np.array([(s.x_pix,s.y_pix) for s in references])
    matcher=TractorGaiaMatcher(
        "synthetic-tractor.fits","c"*64,frame.dataset_id,rows,positions)

    def run(recipe):
        try:
            result=build_empirical_coadd_psf(
                frame,target_icrs_deg=LEO_P,target_protection_arcsec=13.,
                config=recipe,catalog_matcher=matcher)
            return result.training, result.holdout, result.review_snapshot
        except PSF0QualityError as exc:
            from types import SimpleNamespace
            return (
                [SimpleNamespace(**x) for x in exc.diagnostics["training"]],
                [SimpleNamespace(**x) for x in exc.diagnostics["holdout"]],
                exc.review_snapshot)
    train,holdout,snapshot=run(base)
    assert len(train)>=base.min_train_stars
    assert len(holdout)>=base.min_holdout_stars
    assert all(s.gaia_source_id is not None for s in train+holdout)
    assert snapshot["catalog_sha256"]=="c"*64
    assert snapshot["catalog_policy"]=="gaia_strict"
    assert all(item["gaia"]["gaia_source_id"]>0 for item in snapshot["stars"])
    selected=train[0].star_id
    reduced, holdout_again, review_again=run(
        replace(base,excluded_training_star_ids=(selected,)))
    assert selected not in {s.star_id for s in reduced}
    assert {s.star_id for s in holdout}=={
        s.star_id for s in holdout_again}
    assert len(reduced)==len(train)-1
    assert next(item for item in review_again["stars"]
                if item["star_id"]==selected)["role"]=="excluded"


def test_baseline_sky_annulus_has_explicit_star_wing_degeneracy_warning():
    from lsbpipe.photometry.psf_empirical import _star_cutout
    frame=fake_frame()
    mask=np.zeros(frame.shape,dtype=bool)
    src, reason=_star_cutout(
        frame,78.,78.,PSF0Config(),mask)
    assert reason is None
    assert src["sky_qa"]["method"] == "ROBUST_LOCAL_PLANE_ON_STELLAR_ANNULUS"
    assert src["sky_qa"]["wing_sky_degeneracy_warning"] is True
    assert src["sky_qa"]["annulus_inner_in_guessed_fwhm"] == pytest.approx(
        10./4.8)
    assert src["sky_qa"]["sky_uncertainty_covariance_calibrated"] is False
    assert src["sky_qa"]["sky_versus_psf_wings_separable_from_this_cutout"] is False
    assert 0 < src["sky_qa"]["annulus_sigma_clip_retained_fraction"] <= 1


def test_training_outlier_screening_does_not_touch_fixed_holdouts():
    from dataclasses import replace
    from types import SimpleNamespace
    from lsbpipe.photometry.psf_empirical import PSF0QualityError

    frame=fake_frame()
    initial=PSF0Config(
        min_train_stars=6,min_holdout_stars=2,
        max_training_stars=30,min_neighbor_separation_pix=32.)
    def run(cfg):
        try:
            result=build_empirical_coadd_psf(
                frame,target_icrs_deg=LEO_P,
                target_protection_arcsec=13.,config=cfg)
            return result.training,result.holdout,result.review_snapshot
        except PSF0QualityError as exc:
            return (
                [SimpleNamespace(**row) for row in exc.diagnostics["training"]],
                [SimpleNamespace(**row) for row in exc.diagnostics["holdout"]],
                exc.review_snapshot)
    train0,holdout0,_=run(initial)
    train1,holdout1,review=run(
        replace(initial,training_outlier_passes=2))
    assert {s.star_id for s in holdout0} == {s.star_id for s in holdout1}
    assert len(train1) >= 6
    assert review["training_residual_screening_history"]
    assert all(h["uses_holdout_to_select"] is False
               for h in review["training_residual_screening_history"])
    assert all(h["coadd_chi_square_calibrated"] is False
               for h in review["training_residual_screening_history"])
    auto=set(review["automatic_training_rejections"])
    assert auto.isdisjoint({s.star_id for s in holdout0})
    assert auto.isdisjoint({s.star_id for s in train1})
    assert len(auto) <= max(1,int(.2*len(train0)))


def test_brightness_rank_guard_preserves_original_training_and_holdout(
        monkeypatch):
    from dataclasses import replace
    import lsbpipe.photometry.psf_empirical as module

    frame=fake_frame()
    cfg=PSF0Config(
        min_train_stars=6,min_holdout_stars=2,
        max_training_stars=30,min_neighbor_separation_pix=32.)
    def synthetic_brightness_coupled_score(entry, kernel, *,
                                           guessed_fwhm_pix):
        amplitude=float(entry[-1]["peak_nMgy"])
        return {
            "star_id":entry[0],
            "score":1. + 0.001*amplitude + 30.*(amplitude > 175.),
            "amplitude_nMgy_diagnostic":amplitude,
            "metric":"SYNTHETIC_BRIGHTNESS_BIAS_NOT_A_PSF_QUALITY_SCORE",
        }
    monkeypatch.setattr(
        module, "training_star_residual", synthetic_brightness_coupled_score)
    initial=module.build_empirical_coadd_psf(
        frame,target_icrs_deg=LEO_P,
        target_protection_arcsec=13.,config=cfg)
    guarded=module.build_empirical_coadd_psf(
        frame,target_icrs_deg=LEO_P,target_protection_arcsec=13.,
        config=replace(cfg,training_outlier_passes=2))
    evidence=guarded.review_snapshot["training_residual_screening_history"]
    assert evidence and evidence[0]["n_above_cutoff"] > 0
    assert evidence[0]["brightness_score_coupling_qa"]["auto_trim_blocked"]
    assert evidence[0]["proposed_before_brightness_guard"]
    assert evidence[0]["selected_for_exclusion"] == []
    assert evidence[0]["n_removed_this_pass"] == 0
    assert evidence[0]["n_proposed_before_brightness_guard"] > 0
    assert guarded.review_snapshot["automatic_training_rejections"] == ()
    assert len(guarded.training) == len(initial.training)
    assert {s.star_id for s in guarded.holdout} == {
        s.star_id for s in initial.holdout}
    assert {s.star_id for s in guarded.training} == {
        s.star_id for s in initial.training}
    assert np.allclose(guarded.kernel,initial.kernel,atol=1e-12)
    assert "AUTO_TRAINING_TRIM_BLOCKED_BY_SCORE_FLUX_COUPLING" in guarded.flags
