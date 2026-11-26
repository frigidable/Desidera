import numpy as np
import pytest
from astropy.wcs import WCS
from lsbpipe.photometry.stellar_halo_sectors import sector_diagnostics,paired_sector_diagnostics


def fixture():
    w=WCS(naxis=2)
    w.wcs.crpix=[101.,101.]; w.wcs.crval=[96.,78.]
    w.wcs.cdelt=[-.5/3600,.5/3600]
    w.wcs.ctype=["RA---TAN","DEC--TAN"]
    yy,xx=np.indices((201,201))
    a=np.where(yy>=100,2.,-1.)
    a[100,100]=np.nan
    c=dict(radial_knots_arcsec=[0.,100.],halo_at_knots_nMgy_pix=[0.,0.],
        sky_constant_nMgy_pix=0.,sky_x_nMgy_pix_per_arcsec=0.,sky_y_nMgy_pix_per_arcsec=0.)
    p=dict(layers={"original":a.copy(),"residual":a.copy()},
        target_xy_native_pix=[100.,100.],source_xy_native_pix=[100.,20.],
        pixel_scale_arcsec=.5,effective_source_mask_sha256="x")
    return w,p,{"candidates":[c]}


def test_sectors_partition_and_signed_exclusion():
    w,p,c=fixture()
    ra,dec=w.all_pix2world(100.,110.,0)
    e=dict(source_id="test",objid=1,x=100.,y=110.,ra_deg=ra,dec_deg=dec,centre_in_native_frame=True)
    before=p["layers"]["original"].copy()
    d=sector_diagnostics(p,c,w,[e])
    for lo in d["radial_edges_arcsec"][:-1]:
        rows={r["sector"]:r for r in d["rows"] if r["r_inner_arcsec"]==lo}
        assert rows["N"]["used_pixels"]+rows["S"]["used_pixels"]==rows["ALL"]["used_pixels"]
        assert sum(rows[k]["used_pixels"] for k in ("NE","NW","SE","SW"))==rows["ALL"]["used_pixels"]
        assert rows["ALL"]["reference_residual_flux_nMgy"]==pytest.approx(
            rows["N"]["reference_residual_flux_nMgy"]+rows["S"]["reference_residual_flux_nMgy"])
    a=d["individual_source_exclusion_trials"][0]
    assert a["removed_pixels"]>0
    assert a["baseline_reference_flux_nMgy"]-a["excluded_reference_flux_nMgy"]==pytest.approx(a["removed_signed_flux_nMgy"])
    assert np.array_equal(before,p["layers"]["original"],equal_nan=True)


def test_pair_intersection_and_grid_refusal():
    w,p,c=fixture()
    import copy
    base=dict(project_id="p",source_id="s",tractor_catalog_sha256="c",draft_digest="d",
        original_scene_grid_fingerprint="grid",sh1d_policy={"anchor_icrs":[96,78]},
        sh1c_preview=p,original_fits_sha256_local={"image":"hash"},
        sh1_sector_measurement={"sh1b_model_comparison":c})
    g=dict(base,band="g"); r=copy.deepcopy(dict(base,band="r"))
    r["sh1c_preview"]["layers"]["residual"][100,101]=np.nan
    d=paired_sector_diagnostics(g,r,w)
    assert d["common_pixel_count"]==p["layers"]["residual"].size-2
    assert not d["psf_matched"] and not d["science_adopted"]
    for a,b in zip(d["bands"]["g"]["rows"],d["bands"]["r"]["rows"]):
        assert a["used_pixels"]==b["used_pixels"]
        assert a["reference_residual_flux_nMgy"]==b["reference_residual_flux_nMgy"]
    with pytest.raises(ValueError,match="grid_fingerprint"):
        paired_sector_diagnostics(g,dict(r,original_scene_grid_fingerprint="other"),w)
