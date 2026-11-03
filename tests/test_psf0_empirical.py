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
