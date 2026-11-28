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


def test_native_geometry_identity_and_review(monkeypatch):
    WCS=pytest.importorskip("astropy.wcs").WCS
    from lsbpipe.scene import SceneState, MaskLayer, LayerKind
    from lsbpipe.photometry.contracts import NativeBandFrame, PSFDescription
    from lsbpipe.photometry.native_annular_diagnostic import NativeCorrection
    import lsbpipe.photometry.geometry_inference as gi
    w=WCS(naxis=2)
    w.wcs.crpix=[5.,5.]
    w.wcs.crval=[96.,78.]
    w.wcs.cdelt=[-1/3600,1/3600]
    w.wcs.ctype=["RA---TAN","DEC--TAN"]
    shape=(9,9)
    frames=[];scenes=[];corrections=[];psfs=[];identities=[]
    for i,band in enumerate(("g","r")):
        f=NativeBandFrame(np.full(shape,2.+i),np.ones(shape),np.ones(shape),
                          np.zeros(shape,dtype=np.uint32),w,
                          release="dr9",region="south",brick="0965p780",band=band)
        grid=f.native_grid_fingerprint
        mask=np.zeros(shape,dtype=bool);mask[2+i,2]=True
        s=SceneState("synthetic",grid,band,shape,(
            MaskLayer("mask",LayerKind.INSTRUMENTAL_DEFECT,mask,grid,band),))
        c=NativeCorrection("residual_sky",grid,f.native_input_sha256,
                           "synthetic fixed correction",np.ones(shape))
        k=np.zeros((3,3));k[1,1]=1
        frames.append(f);scenes.append(s);corrections.append((c,))
        psfs.append(PSFDescription("EMPIRICAL_KERNEL",k,{"fixture":"synthetic"}))
        identities.append((f.native_input_sha256,grid))
    def probe(e,n,images,kernels,*,valid,region,band_order,**kw):
        assert band_order==("g","r")
        assert not valid[2,2] and not valid[3,2]
        assert valid.sum()==79
        np.testing.assert_allclose(images[0][valid],1.)
        np.testing.assert_allclose(images[1][valid],2.)
        assert abs(e[4,4])+abs(n[4,4])<1e-7
        return {"status":"DIAGNOSTIC_NOT_ADOPTED"}
    monkeypatch.setattr(gi,"fit_block_geometry",probe)
    kw=dict(anchor_icrs=(96.,78.),corrections_by_band=corrections,
        psfs=psfs,psf_identities=identities,region=np.ones(shape,dtype=bool),
        region_grid_id=grid,train_tiles=np.ones((3,3),dtype=bool),
        holdout_tiles=np.zeros((3,3),dtype=bool))
    result=gi.fit_native_geometry_candidates(frames,scenes,**kw)
    assert not result["geometry_adopted"]
    with pytest.raises(ValueError,match="PSF"):
        gi.fit_native_geometry_candidates(frames,scenes,**dict(kw,psf_identities=identities[::-1]))
    node={"east_arcsec":0.,"north_arcsec":0.,"q":.6,"pa_deg":45.}
    candidate={"success":True,"boundary":False,"geometry":{"inner":node,
        "outer":node,"transition_start":2.,"transition_end":8.}}
    result["candidates"]={"constant":candidate}
    with pytest.raises(ValueError):
        gi.reviewed_geometry_candidate(result,"constant",confirmed=False)
    assert gi.reviewed_geometry_candidate(result,"constant",confirmed=True).inner.q==.6
    candidate["boundary"]=True
    with pytest.raises(ValueError):
        gi.reviewed_geometry_candidate(result,"constant",confirmed=True)
