import numpy as np
import pytest
from scipy.spatial.distance import cdist

from lsbpipe.photometry.block_covariance import (
    correction_influence, within_tile_diagonal_remainder,
    conditional_block_covariance,
)


def test_sky_halo_influence_coupling_and_disjoint_training():
    counts=np.array([[3,2,0,0],[0,1,0,0]])
    x=np.array([[5.],[1.]])
    h=np.array([2.,.4])
    b=np.array([[0.,0.,.5,.5]])
    a=np.array([0.,0.,-1.,1.])
    train=np.array([False,False,True,True])
    l=correction_influence(counts,x,h,b,a,training_tiles=train)
    np.testing.assert_allclose(l,counts-x@b-np.outer(h,a))
    values=np.array([.3,.7,1.,1.2])
    np.testing.assert_allclose(l@values,counts@values-x@(b@values)-h*(a@values))
    with pytest.raises(ValueError):
        correction_influence(counts,x,h,b,a,training_tiles=np.array([False,True,True,True]))
    with pytest.raises(ValueError):
        correction_influence(counts,x,h,b,np.ones(4),training_tiles=train)


def test_within_tile_covariance_is_orthogonal_to_tile_sum():
    v=np.array([[1.,3.,0.],[2.,0.,0.]])
    d=within_tile_diagonal_remainder(v)
    np.testing.assert_allclose(d,[[2/3,-2/3],[-2/3,2/3]])
    np.testing.assert_allclose(d@np.ones(2),0,atol=1e-14)
    assert np.linalg.eigvalsh(d).min()>-1e-14


def test_spatial_block_covariance_matches_full_kron():
    xy=np.array([[0.,0.],[3.,0.],[0.,4.]])
    l=[np.array([[1.,2.,-1.],[0.,2.,-2.]]),np.array([[2.,-1.,0.]])]
    rho=np.array([[1.,.4],[.4,1.]])
    scales=np.array([.2,.5])
    rem=[np.eye(2)*.1,np.eye(1)*.2]
    from scipy.linalg import block_diag
    k=.6*np.exp(-cdist(xy,xy)/3)+.4*np.eye(3)
    op=block_diag(*l)
    truth=op@np.kron(rho*np.outer(scales,scales),k)@op.T+block_diag(*rem)
    result=conditional_block_covariance(xy,l,scales,rho,rem,band_order=("g","r"),
                                        length_arcsec=3,mix=.6)
    np.testing.assert_allclose(result["covariance"],truth,atol=1e-14)
    assert not result["covariance_calibrated"]
    for bad in ([[1,2],[2,1]],[[2,0],[0,1]]):
        with pytest.raises(ValueError):
            conditional_block_covariance(xy,l,scales,bad,rem,band_order=("g","r"),
                                         length_arcsec=3,mix=.6)
    with pytest.raises(ValueError):
        conditional_block_covariance(xy,l,scales,rho,rem,band_order=("g","r"),
                                     length_arcsec=3,mix=.6,max_tiles=2)


def test_single_ring_exact_zero_and_paired_controls():
    from lsbpipe.photometry.block_covariance import paired_control_statistics
    assert np.array_equal(within_tile_diagonal_remainder([[.1,.2,3.7]]),[[0.]])
    r=np.array([[-2.,-3.],[-1.,-1.],[0.,.2],[1.,1.2],[2.,4.]])
    result=paired_control_statistics(r,band_order=("g","r"),tile_ids=range(5),
                                     minimum_tiles=4,clip_sigma=3.)
    assert result["n_matched_tiles"]==5
    np.testing.assert_allclose(result["scales"],1.4826*np.median(
        np.abs(r-np.median(r,axis=0)),axis=0))
    np.testing.assert_allclose(np.diag(result["correlation"]),1.)
    assert np.linalg.eigvalsh(result["correlation"]).min()>=0
    with pytest.raises(ValueError):
        paired_control_statistics(r,band_order=("g","r"),tile_ids=[0,0,1,2,3],
                                  minimum_tiles=4,clip_sigma=3.)


def test_noise_split_and_exact_independent_ivar_propagation():
    from lsbpipe.photometry.block_covariance import (
        split_control_block_noise, independent_ivar_covariance,
    )
    rho=np.array([[1.,.3],[.3,1.]])
    result=split_control_block_noise([2.,3.],rho,[[1.,2.],[1.,2.]])
    s=result["coherent_scales"]
    np.testing.assert_allclose(result["coherent_correlation"]*np.outer(s,s),
                               [[3.,1.8],[1.8,7.]])
    with pytest.raises(ValueError):
        split_control_block_noise([1.,1.],[[1.,.9],[.9,1.]],[[.5,.5]])
    q=np.array([[-1.,.5],[2.,1.]])
    v=np.array([.2,.3])
    t=np.array([.1,.4])
    np.testing.assert_allclose(independent_ivar_covariance(v,q,t),
                               np.diag(v)+q@np.diag(t)@q.T)


def test_native_scene_block_covariance_pipeline():
    WCS=pytest.importorskip("astropy.wcs").WCS
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    from lsbpipe.scene import SceneState
    from lsbpipe.photometry.contracts import NativeBandFrame, PSFDescription
    from lsbpipe.photometry.stellar_contamination import (
        compact_halo_template, fit_native_independent_sky_halo,
    )
    from lsbpipe.photometry.radial_models import EllipseNode, NestedEllipse
    from lsbpipe.photometry.native_annular_diagnostic import (
        common_native_validity, prepare_native_annuli,
    )
    from lsbpipe.photometry.block_covariance import native_scene_block_covariance
    w=WCS(naxis=2)
    w.wcs.crpix=[41.,41.]
    w.wcs.crval=[96.,78.]
    w.wcs.cdelt=[-1/3600,1/3600]
    w.wcs.ctype=["RA---TAN","DEC--TAN"]
    yy,xx=np.indices((81,81))
    world=w.pixel_to_world(xx,yy).icrs
    star=SkyCoord(96*u.deg,78*u.deg,frame="icrs")
    gal=w.pixel_to_world(20.,60.).icrs
    gi=(gal.ra.deg,gal.dec.deg)
    rad=star.separation(world).arcsec
    pa=star.position_angle(world).rad
    template=compact_halo_template(rad,cutoff_arcsec=20,power=0,pivot_arcsec=8,core_arcsec=4)
    base=.02+.001*rad*np.sin(pa)/25-.002*rad*np.cos(pa)/25
    rng=np.random.default_rng(1234)
    noise=rng.normal(0,.001,(2,27,27))
    frames=[];scenes=[];fits=[]
    for i,band in enumerate(("g","r")):
        fluct=noise[i,yy//3,xx//3]+.3*noise[0,yy//3,xx//3]
        im=base+(.03+i*.01)*np.nan_to_num(template,nan=0.)+fluct
        im[rad<4]=99.
        f=NativeBandFrame(im,np.full(im.shape,1e8),np.ones(im.shape),
                          np.zeros(im.shape,dtype=np.uint32),w,
                          release="dr9",region="south",brick="0965p780",band=band)
        scene=SceneState("synthetic",f.native_grid_fingerprint,band,f.shape,())
        fit=fit_native_independent_sky_halo(f,scene,star_icrs=(96.,78.),galaxy_icrs=gi,
            core_arcsec=4,cutoff_arcsec=20,power=0,pivot_arcsec=8,
            galaxy_guard_arcsec=11,far_star_guard_arcsec=25,rays=(),block_size=3,
            minimum_block_fraction=.8,holdout_pa_deg=0,holdout_width_deg=90,buffer_deg=10)
        for layer in fit["proposed_layers"]:
            scene=scene.with_layer(layer)
        frames.append(f);scenes.append(scene);fits.append(fit)
    common=common_native_validity(frames,scenes)
    assert not common[40,40]
    kernel=np.zeros((3,3));kernel[1,1]=1
    psf=PSFDescription("EMPIRICAL_KERNEL",kernel,{"fixture":"synthetic"})
    node=EllipseNode(0,0,.9,0)
    geometry=NestedEllipse(node,node,2,10)
    records=[]
    for f,s,h in zip(frames,scenes,fits):
        records.append(prepare_native_annuli(f,s,geometry,
            anchor_icrs=gi,ring_edges=[0,2,4,6],radial_nodes=np.linspace(0,120,601),
            corrections=(h["sky"],h["halo"]),psf=psf,
            psf_input_sha256=f.native_input_sha256,psf_grid_id=f.native_grid_fingerprint,
            n_sectors=4,ring_min=.5,cumulative_min=.5,sector_min=.4,
            common_valid=common,retain_pixel_support=True))
    result=native_scene_block_covariance(frames,records,fits,length_arcsec=5,mix=.5,
                                           minimum_control_tiles=8,clip_sigma=3.)
    assert result["covariance"].shape==(6,6)
    assert result["n_matched_far_controls"]>=8
    assert np.linalg.eigvalsh(result["covariance"]).min()>0
    assert not result["covariance_calibrated"]
    assert result["covariance_model"]=="ADDITIVE_COHERENT_TILE_AND_INDEPENDENT_IVAR"
    bad=dict(records[0],ring_index_sha256="not-the-original")
    with pytest.raises(ValueError):
        native_scene_block_covariance(frames,[bad,records[1]],fits,
            length_arcsec=5,mix=.5,minimum_control_tiles=8,clip_sigma=3.)
