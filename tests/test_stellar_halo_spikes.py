import numpy as np
import pytest
from lsbpipe.photometry.stellar_halo_spikes import (
    coords,safe_pixels,transverse_rows,mirrored_template,scenario_sums,ensemble_refit_deltas)


def test_away_template_ignores_protected_target():
    shape=(601,601); star=(300.,300.); target=(300.,440.); scale=.25
    t,q=coords(shape,star,target,scale)
    a=np.where(np.abs(q)<2.,.003,0.)+np.zeros(shape)
    a[t>0]+=10.
    z=np.zeros(shape,dtype=bool)
    safe=safe_pixels(a,np.ones(shape),z,z,t>0,z,[],scale)
    rows=transverse_rows(a,safe,t,q,2.)
    assert any(row["template_supported"] for row in rows)
    assert all(row["signed_excess_nMgy_pix"]==pytest.approx(.003)
        for row in rows if row["template_supported"])
    m=mirrored_template(t,q,2.,rows)
    assert m[440,300]==pytest.approx(.003)
    assert m[160,300]==pytest.approx(.003)
    assert m[440,330]==0.
    assert np.isnan(m[300,300])
    assert np.isnan(m[600,300])
    assert np.all(a[t>0]>=10.)


def test_matched_subsets_and_flux_conservation():
    shape=(201,201); star=(100.,30.); target=(100.,100.); scale=.5
    t,q=coords(shape,star,target,scale); strip=np.abs(q)<2.
    a=np.ones(shape); b=a-.1; r=a-.2
    s=np.where(strip,.03,0.)+np.zeros(shape)
    s[100,100]=np.nan; r[100,101]=np.nan
    rows=scenario_sums(a,b,r,s,strip,target,scale)
    row=next(v for v in rows if v["kind"]=="aperture" and v["outer_arcsec"]==20 and v["sector"]=="ALL")
    assert row["off_strip_pixels"]<row["common_pixels"]
    assert row["template_supported_pixels"]==row["common_pixels"]-1
    assert row["refit_template_matched_flux_nMgy"]-row["spike_corrected_flux_nMgy"]==pytest.approx(row["spike_template_flux_nMgy"])
    assert row["baseline_off_strip_flux_nMgy"]-row["refit_off_strip_flux_nMgy"]==pytest.approx(.1*row["off_strip_pixels"])


def test_refit_ensemble_deltas_on_identical_pixels():
    a=np.ones((201,201)); base=a.copy(); base[100,100]=np.nan
    p=dict(layers={"original":a,"residual":base},source_xy_native_pix=[100.,30.],
           target_xy_native_pix=[100.,100.],pixel_scale_arcsec=.5)
    c=dict(tile_size_pix=24,assumed_zero_halo_radius_arcsec=70.,sky_model="constant",
           radial_knots_arcsec=[10.,70.],halo_at_knots_nMgy_pix=[0.,0.],
           sky_constant_nMgy_pix=.1,sky_x_nMgy_pix_per_arcsec=0.,sky_y_nMgy_pix_per_arcsec=0.,
           training_r_range_arcsec=[10.,80.])
    d=dict(c,sky_constant_nMgy_pix=.2,training_r_range_arcsec=[25.,80.])
    rows=ensemble_refit_deltas(p,[c],{(24,70.,"constant"):d},2.)
    for row in rows:
        assert row["used_pixels"]>0
        assert row["common_radial_support_arcsec"]==[25.,80.]
        assert row["refit_delta_range_nMgy"]==pytest.approx([-.1*row["used_pixels"]]*2)
    assert rows[1]["used_pixels"]<rows[0]["used_pixels"]
    assert np.isnan(base[100,100])
