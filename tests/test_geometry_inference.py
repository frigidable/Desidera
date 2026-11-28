import numpy as np
import pytest
from scipy.signal import fftconvolve

from lsbpipe.photometry.geometry_inference import _radius, fit_block_geometry
from lsbpipe.photometry.radial_models import EllipseNode, NestedEllipse, sersic_intensity


def fixture():
    n,e=np.mgrid[-19:20,-19:20]
    seed=EllipseNode(0.,0.,.66,35.)
    truth=NestedEllipse(EllipseNode(0.,-.1,.6375,33.5),
                        EllipseNode(.3,-.3,.6825,36.5),4.,16.)
    radius=truth.radius(e,n)
    kernels=(np.array([[0,1,0],[1,4,2],[0,3,0]],dtype=float)/11,
             np.array([[1,2,1],[2,4,2],[1,2,1]],dtype=float)/16)
    profiles=np.array([[.05,7.,.7],[.07,6.,.8]])
    data=[fftconvolve(sersic_intensity(radius,*p),k,mode="same")
          for p,k in zip(profiles,kernels)]
    valid=np.ones(e.shape,dtype=bool)
    valid[15:17,16:18]=False
    for image in data:
        image[~valid]=np.nan
    region=np.hypot(e,n)<=17.
    train=np.zeros((13,13),dtype=bool)
    held=np.zeros_like(train)
    train[1:12,1:8]=True
    held[1:12,9]=True
    kwargs=dict(band_order=("g","r"),valid=valid,region=region,
        train_tiles=train,holdout_tiles=held,block_size=3,
        min_block_fraction=.5,buffer_tiles=1,block_scales=[1e-5,1e-5],
        seed=seed,transition_start=4.,transition_end=16.,
        midpoint_bounds=([-.7,-.7,.58,30.],[.7,.7,.75,40.]),
        delta_limits=[.6,.6,.06,5.],delta_scales=[1.,1.,.1,10.],
        profile_initial=profiles,
        profile_bounds=([[.01,3.,.35]]*2,[[.2,12.,1.4]]*2),max_nfev=120)
    return e,n,data,kernels,truth,kwargs


def test_fast_implicit_radius_matches_existing_contract():
    e,n,_,_,truth,_=fixture()
    np.testing.assert_allclose(_radius(truth,e,n),truth.radius(e,n),atol=2e-10,rtol=0)


def test_psf_geometry_recovery_and_frozen_holdout():
    e,n,data,kernels,truth,kwargs=fixture()
    result=fit_block_geometry(e,n,data,kernels,**kwargs)
    a=result["candidates"]
    assert all(r["success"] for r in a.values())
    assert a["variable"]["training_weighted_score"] < a["constant"]["training_weighted_score"]
    assert a["variable"]["holdout_weighted_score"] < a["constant"]["holdout_weighted_score"]
    assert a["variable"]["holdout_mae_by_band"].max()<2e-5
    for endpoint in ("inner","outer"):
        assert a["variable"]["geometry"][endpoint]["q"]==pytest.approx(getattr(truth,endpoint).q,abs=.015)
    assert result["selected_geometry"] is None
    assert not result["covariance_calibrated"]
    assert not np.any(result["train"]&result["holdout"])
    # Inject an arbitrary signal into the held-out tile columns only.
    # Training solutions and regularization must not change.
    changed=[image.copy() for image in data]
    for image in changed:
        image[:,27:30]+=.005
    again=fit_block_geometry(e,n,changed,kernels,**kwargs)
    for mode in a:
        r=again["candidates"][mode]
        np.testing.assert_allclose(r["log_profile_parameters"],a[mode]["log_profile_parameters"],
                                   rtol=0,atol=1e-9)
        assert r["geometry"]==a[mode]["geometry"]
        assert r["holdout_weighted_score"]>a[mode]["holdout_weighted_score"]


def test_reject_non_nested_box_and_holdout_overlap():
    e,n,data,kernels,_,kwargs=fixture()
    with pytest.raises(ValueError,match="nested"):
        fit_block_geometry(e,n,data,kernels,**dict(kwargs,delta_limits=[8.,8.,.3,30.]))
    bad=kwargs["train_tiles"].copy()
    bad[3,8]=True
    with pytest.raises(ValueError,match="buffer"):
        fit_block_geometry(e,n,data,kernels,**dict(kwargs,train_tiles=bad))
    with pytest.raises(ValueError,match="disjoint"):
        fit_block_geometry(e,n,data,kernels,**dict(kwargs,train_tiles=kwargs["holdout_tiles"]))
