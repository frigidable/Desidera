import numpy as np
import pytest

from lsbpipe.photometry.contracts import NativeBandFrame, PSFDescription
from lsbpipe.photometry.radial_models import EllipseNode, NestedEllipse, sersic_intensity
from lsbpipe.photometry.native_annular_diagnostic import (
    NativeCorrection, prepare_native_annuli, fit_common_support_sersic,
)
from lsbpipe.scene import SceneState, MaskLayer, LayerKind


def test_native_joint_workflow_on_synthetic_frames():
    WCS=pytest.importorskip("astropy.wcs").WCS
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    w=WCS(naxis=2)
    w.wcs.crpix=[21.,21.]
    w.wcs.crval=[96.,78.]
    w.wcs.cdelt=[-1/3600,1/3600]
    w.wcs.ctype=["RA---TAN","DEC--TAN"]
    yy,xx=np.indices((41,41))
    center=SkyCoord(96*u.deg,78*u.deg,frame="icrs")
    coords=w.pixel_to_world(xx,yy).icrs
    distance=center.separation(coords).arcsec
    pa=center.position_angle(coords).rad
    node=EllipseNode(0,0,.8,0)
    geometry=NestedEllipse(node,node,2,10)
    radius=geometry.radius(distance*np.sin(pa),distance*np.cos(pa))
    truth=[.06,5.,.6,.09,4.5,.7]
    kernel=np.zeros((3,3));kernel[1,1]=1.
    psf=PSFDescription("EMPIRICAL_KERNEL",kernel,{"fixture":"synthetic delta"})
    records=[]
    for i,band in enumerate(("g","r")):
        im=.02+sersic_intensity(radius,*truth[3*i:3*i+3])
        frame=NativeBandFrame(im,np.full(im.shape,1e4),np.ones(im.shape),
                              np.zeros(im.shape,dtype=np.uint32),w,
                              release="dr9",region="south",brick="0965p780",band=band)
        grid=frame.native_grid_fingerprint
        hard=np.zeros(im.shape,dtype=bool);hard[21,21]=True
        scene=SceneState("synthetic",grid,band,frame.shape,(
            MaskLayer("hard",LayerKind.INSTRUMENTAL_DEFECT,hard,grid,band),))
        sky=NativeCorrection("residual_sky",grid,frame.native_input_sha256,
                             "synthetic known sky",np.full(im.shape,.02))
        params=dict(anchor_icrs=(96,78),ring_edges=np.arange(0,14,2.),
                    radial_nodes=np.arange(0,41,.02),psf=psf,
                    psf_input_sha256=frame.native_input_sha256,psf_grid_id=grid,
                    n_sectors=4,ring_min=.5,cumulative_min=.5,sector_min=.4)
        rec=prepare_native_annuli(frame,scene,geometry,corrections=(sky,),**params)
        assert rec["coverage"]["eligible"].all()
        assert rec["status"]=="DIAGNOSTIC_NOT_ADOPTED"
        assert not rec["full_geometry_containment_verified"]
        assert rec["observed_flux_nMgy"].shape==(6,)
        assert rec["operator"].shape==(6,len(params["radial_nodes"]))
        np.testing.assert_allclose(rec["operator"]@sersic_intensity(
            params["radial_nodes"],*truth[3*i:3*i+3]),rec["observed_flux_nMgy"],
            rtol=1e-4,atol=1e-6)
        wrong=NativeCorrection("residual_sky","wrong",frame.native_input_sha256,
                               "wrong grid",np.full(im.shape,.02))
        with pytest.raises(ValueError):
            prepare_native_annuli(frame,scene,geometry,corrections=(wrong,),**params)
        records.append(rec)
    train=np.array([True]*5+[False])
    fit=fit_common_support_sersic(records,np.eye(12)*.0004,band_order=("g","r"),
        covariance_source="synthetic conditional covariance",train_rings=train,test_rings=~train,
        initial=np.log([.05,4.,.7,.08,5.,.6]),
        bounds=(np.log([.001,1.,.3]*2),np.log([1.,20.,2.]*2)))
    assert fit["success"] and not fit["boundary"]
    np.testing.assert_allclose(np.exp(fit["parameters"]),truth,rtol=.005)
    assert fit["holdout"]["chi2"]<.001
    assert fit["model_extrapolations"]["g"]["kind"]=="UNTRUNCATED_MODEL_EXTRAPOLATION"
    assert not fit["covariance_calibrated"]
    with pytest.raises(ValueError):
        fit_common_support_sersic(records,np.eye(12),band_order=("r","g"),
            covariance_source="wrong band order",train_rings=train,test_rings=~train,
            initial=np.log(truth),bounds=(np.log([.001,1.,.3]*2),np.log([1.,20.,2.]*2)))


def test_native_independent_sky_halo_with_buffered_sector():
    WCS=pytest.importorskip("astropy.wcs").WCS
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    from lsbpipe.photometry.stellar_contamination import (
        compact_halo_template, fit_native_independent_sky_halo,
    )
    w=WCS(naxis=2)
    w.wcs.crpix=[41.,41.]
    w.wcs.crval=[96.,78.]
    w.wcs.cdelt=[-1/3600,1/3600]
    w.wcs.ctype=["RA---TAN","DEC--TAN"]
    yy,xx=np.indices((81,81))
    center=SkyCoord(96*u.deg,78*u.deg,frame="icrs")
    coords=w.pixel_to_world(xx,yy).icrs
    r=center.separation(coords).arcsec
    pa=center.position_angle(coords).rad
    sky=.02+.001*r*np.sin(pa)/25-.002*r*np.cos(pa)/25
    h=compact_halo_template(r,cutoff_arcsec=20,power=0,pivot_arcsec=8,core_arcsec=4)
    image=sky+.03*np.nan_to_num(h,nan=0.)
    image[r<4]=99.
    frame=NativeBandFrame(image,np.full(image.shape,1e4),np.ones(image.shape),
                          np.zeros(image.shape,dtype=np.uint32),w,
                          release="dr9",region="south",brick="0965p780",band="r")
    scene=SceneState("synthetic",frame.native_grid_fingerprint,"r",frame.shape,())
    result=fit_native_independent_sky_halo(frame,scene,
        star_icrs=(96.,78.),galaxy_icrs=(96.1,78.),core_arcsec=4,cutoff_arcsec=20,
        power=0,pivot_arcsec=8,galaxy_guard_arcsec=10,far_star_guard_arcsec=25,
        rays=(),block_size=3,minimum_block_fraction=.8,holdout_pa_deg=0,
        holdout_width_deg=90,buffer_deg=10,sky_degree=1)
    assert result["amplitude"]==pytest.approx(.03,abs=1e-9)
    np.testing.assert_allclose(result["sky_coefficients"],[.02,.001,-.002],atol=1e-9)
    assert result["holdout_rmse"]<1e-9
    assert np.isnan(result["halo"].values[40,40])
    assert result["proposed_layers"][0].data[40,40]
    assert not np.any(result["far"]&result["holdout"])
    assert not result["covariance_calibrated"]
