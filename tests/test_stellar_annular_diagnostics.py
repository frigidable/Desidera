"""Offline numerical checks; no synthetic result is science adoption."""
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.signal import fftconvolve

from lsbpipe.photometry.stellar_contamination import (
    RayMask, stellar_exclusion_layers, compact_halo_template,
    fit_independent_sky_halo,
)
from lsbpipe.photometry.annular_inference import (
    annular_support, adjoint_annular_operator, propagate_linear_covariance,
    fit_conditional_gls, conditional_holdout,
)
from lsbpipe.photometry.uncertainty_budget import (
    conditional_magnitude_budget, linear_sigma,
)
from lsbpipe.scene import SceneState, MaskLayer, LayerKind, evaluate_mask


def test_native_stellar_mask_direction_and_no_restore():
    WCS = pytest.importorskip("astropy.wcs").WCS
    w = WCS(naxis=2)
    w.wcs.crpix = [16., 16.]
    w.wcs.crval = [96., 78.]
    w.wcs.cdelt = [-1/3600, 1/3600]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    frame = SimpleNamespace(wcs=w, shape=(31,31), band="r")
    core, ray = stellar_exclusion_layers(
        frame, source_id="star", ra_deg=96., dec_deg=78., core_arcsec=3.,
        rays=(RayMask(0., ((3.,4.), (12.,2.))),), row_chunk=7)
    assert core.data[15,15] and not core.data[25,15]
    assert ray.data[25,15] and not ray.data[5,15]
    assert not ray.data[29,15]
    assert not core.data.flags.writeable
    hard = np.zeros(frame.shape,dtype=bool)
    hard[0,0] = True
    scene = SceneState("test",core.grid_id,"r",frame.shape,(
        MaskLayer("native",LayerKind.SATURATION,hard,core.grid_id,"r"),
        core,ray,
        MaskLayer("protect",LayerKind.TARGET_PROTECT,np.ones(frame.shape,dtype=bool),core.grid_id,"r"),
    ))
    valid=evaluate_mask(scene,np.ones(frame.shape),np.ones(frame.shape),
                        grid_id=core.grid_id,operation="photometry").valid
    assert not valid[0,0] and not valid[25,15] and not valid[15,15]
    assert valid[5,15]


def test_halo_template_core_and_cutoff():
    a=compact_halo_template([10,26,30,69.99,70,100],cutoff_arcsec=70,
                            power=1,pivot_arcsec=30,core_arcsec=26)
    assert np.isnan(a[0]) and a[1]>1
    assert a[2]==pytest.approx(1)
    assert 0<a[3]<1e-6
    assert np.array_equal(a[4:],[0,0])
    with pytest.raises(ValueError):
        compact_halo_template([30],cutoff_arcsec=25,power=0,pivot_arcsec=30,core_arcsec=26)


def test_separate_halo_sky_and_influence():
    xx=np.linspace(-1,1,12)
    x=np.column_stack((np.ones(12),xx))
    h=np.r_[np.zeros(5),1.,2.,3.,4.,1.,2.,3.]
    y=x@np.array([.002,.001])+h*.003
    far=np.arange(12)<5
    near=(np.arange(12)>=5)&(np.arange(12)<9)
    test=~(far|near)
    r=fit_independent_sky_halo(y,x,h,far=far,near=near,holdout=test)
    np.testing.assert_allclose(r["sky_coefficients"],[.002,.001],atol=1e-10)
    assert r["amplitude"]==pytest.approx(.003,abs=1e-10)
    assert r["holdout_rmse"]<1e-10
    assert np.all(r["sky_influence"][:,test]==0)
    assert np.all(r["halo_influence"][test]==0)
    np.testing.assert_allclose(r["halo_influence"]@y,.003,atol=1e-10)
    with pytest.raises(ValueError):
        fit_independent_sky_halo(y,x,h,far=far,near=near|far,holdout=test)
    with pytest.raises(ValueError):
        fit_independent_sky_halo(y,np.ones_like(x),h,far=far,near=near,holdout=test)


def test_coverage_must_be_contiguous_and_missing_sector_fails():
    rings=np.array([[0,0,1,1,2,2]])
    sectors=np.array([[0,1,0,1,0,1]])
    valid=np.array([[True,True,True,False,True,True]])
    r=annular_support(rings,sectors,valid,n_rings=3,n_sectors=2,
                      ring_min=.5,cumulative_min=.5,sector_min=.5)
    assert r["eligible"].tolist()==[True,False,False]
    assert r["coverage_scope"]=="IN_FRAME_ONLY"


def test_psf_adjoint_matches_direct_asymmetric_convolution():
    y,x=np.mgrid[-4:5,-5:6]
    r=np.hypot(x,y)
    nodes=np.linspace(0,8,81)
    profile=np.exp(-nodes/2)
    masks=np.array([r<2,(r>=2)&(r<4)])
    masks[:,4,5]=False
    kernel=np.array([[0,1,0],[1,4,2],[0,3,0]],dtype=float)/11
    a=adjoint_annular_operator(r,nodes,masks,kernel,grid_id="same",mask_grid_id="same")
    im=fftconvolve(np.interp(r,nodes,profile),kernel,mode="same")
    np.testing.assert_allclose(a@profile,[im[m].sum() for m in masks],atol=1e-12)
    with pytest.raises(ValueError):
        adjoint_annular_operator(r,nodes,masks,kernel,grid_id="one",mask_grid_id="two")


def test_covariance_and_gls_holdout():
    c=np.array([[1.,.2,.3],[.2,1.,.4],[.3,.4,1.]])
    l=np.array([[1.,-1.,0.],[0.,0.,1.]])
    np.testing.assert_allclose(propagate_linear_covariance(l,c,np.eye(2)),
                               l@c@l.T+np.eye(2))
    train=np.array([True,True,False])
    test=~train
    d=np.array([2.,2.,3.])
    fit=fit_conditional_gls(d,c,lambda p:np.repeat(p[0],3),[1.],([0.],[5.]),train=train)
    assert fit["success"] and not fit["boundary"]
    np.testing.assert_allclose(fit["parameters"],[2.],atol=1e-8)
    np.testing.assert_allclose(fit["conditional_parameter_covariance"],[[.6]],atol=1e-8)
    result=conditional_holdout(d,fit["prediction"],c,np.ones((3,1)),
                               fit["conditional_parameter_covariance"],train=train,test=test)
    a=c[2:3,:2]@np.linalg.inv(c[:2,:2])
    expected=1-(a@c[:2,2:3]).item()+.6*(1-a.sum())**2
    assert result["predictive_covariance"][0,0]==pytest.approx(expected)
    assert result["chi2"]==pytest.approx(1/expected)
    with pytest.raises(ValueError):
        conditional_holdout(d,fit["prediction"],c,np.ones((3,1)),
                            [[.6]],train=train,test=train)


def test_real_ngc2146sat_conditional_calibration_regression():
    # Numerical excerpt: NGC2146sat_uncertainty_budget.zip / manifest.json.
    # Source core-mask-gls.zip SHA256:
    # b677b20fa4d53beea80bf315fbeceeea1db0e00c36c2c7cc7cb125f684e38cb7
    # Reference core26_compact_refit_50_L15_A50, band order g,r.
    stat=[[.005344788351425037,.001112544539598875],
          [.0011125445395988752,.005068015528628374]]
    b=conditional_magnitude_budget(stat,[3.214,2.165],.009500338882207872,
                                   np.eye(2)*.02**2)
    c=b["conditional_plus_assumed_cal_cov"]
    np.testing.assert_allclose(c,[[.006677118952708755,.0017405768207996368],
                                  [.0017405768207996368,.00589106776534264]],rtol=1e-10)
    assert linear_sigma(c,[1,-1])==pytest.approx(.095327,abs=1e-5)


def test_real_control_comparison_contract_on_synthetic_vectors():
    from lsbpipe.photometry.annular_inference import control_model_comparison
    x=np.arange(4,dtype=float)
    models={
        "constant":(lambda p:np.repeat(p[0],4),[1.],([-10.],[10.])),
        "linear":(lambda p:p[0]+p[1]*x,[1.,0.],([-10.,-10.],[10.,10.])),
    }
    kwargs=dict(null_model="constant",train=np.ones(4,dtype=bool),
                reference_support_id="same",control_support_ids=("same","same"),
                control_ids=("a","b"),overlap=np.eye(2,dtype=bool),
                control_provenance="synthetic test only")
    result=control_model_comparison(2+x,np.repeat(2.,4),[.2*x,2*x],
                                    np.eye(4),models,**kwargs)
    assert result["target"]["delta_chi2"]==pytest.approx(5.)
    assert result["n_equal_or_exceed"]==1
    assert result["calibrated_false_positive_rate"] is None
    kwargs["overlap"]=np.ones((2,2),dtype=bool)
    with pytest.raises(ValueError):
        control_model_comparison(2+x,np.repeat(2.,4),[.2*x,2*x],
                                 np.eye(4),models,**kwargs)
