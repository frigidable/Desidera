import numpy as np
import pytest
from scipy.signal import fftconvolve

from lsbpipe.photometry.native_workflow import (
    WorkflowBlocked, prepare_scene_candidates, run_reviewed_annular_workflow,
)


def test_review_gates():
    with pytest.raises(WorkflowBlocked,match="MASK_REVIEW"):
        run_reviewed_annular_workflow(None,masks_confirmed=False,geometry_confirmed=False,
            geometry=None,anchor_icrs=(96,78),annular_settings={},
            covariance_settings={},fit_settings={})


def test_full_native_workflow(tmp_path):
    WCS=pytest.importorskip("astropy.wcs").WCS
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    from lsbpipe.scene import SceneState
    from lsbpipe.photometry.contracts import NativeBandFrame, PSFDescription
    from lsbpipe.photometry.radial_models import EllipseNode, NestedEllipse, sersic_intensity
    from lsbpipe.photometry.stellar_contamination import compact_halo_template
    from lsbpipe.photometry.native_workflow import geometry_candidates
    import zipfile
    import json

    w=WCS(naxis=2)
    w.wcs.crpix=[41.,41.]
    w.wcs.crval=[96.,78.]
    w.wcs.cdelt=[-1/3600,1/3600]
    w.wcs.ctype=["RA---TAN","DEC--TAN"]
    yy,xx=np.indices((81,81))
    world=w.pixel_to_world(xx,yy).icrs
    star=SkyCoord(96*u.deg,78*u.deg,frame="icrs")
    gal=w.pixel_to_world(20.3,60.3).icrs
    anchor=(gal.ra.deg,gal.dec.deg)
    r=star.separation(world).arcsec
    pa=star.position_angle(world).rad
    halo=compact_halo_template(r,cutoff_arcsec=20,power=0,pivot_arcsec=8,core_arcsec=4)
    sky=.02+.001*r*np.sin(pa)/25-.002*r*np.cos(pa)/25
    node=EllipseNode(0.,0.,.9,10.)
    geo=NestedEllipse(node,node,2,10)
    gr=gal.separation(world).arcsec
    gp=gal.position_angle(world).rad
    radius=geo.radius(gr*np.sin(gp),gr*np.cos(gp))
    rng=np.random.default_rng(20261128)
    noise=rng.normal(0,.00008,(2,27,27))
    kernels=[np.array([[0,1,0],[1,4,2],[0,3,0]],float)/11,
             np.array([[1,2,1],[2,4,2],[1,2,1]],float)/16]
    truth=np.array([[.06,3.5,.6],[.08,3.3,.65]])
    frames=[];scenes=[];psfs=[];identities=[];settings=[]
    for i,b in enumerate(("g","r")):
        light=fftconvolve(sersic_intensity(radius,*truth[i]),kernels[i],mode="same")
        fluct=noise[i,yy//3,xx//3]+.25*noise[0,yy//3,xx//3]
        image=sky+(.03+.01*i)*np.nan_to_num(halo,nan=0.)+light+fluct
        image[r<4]=99.
        f=NativeBandFrame(image,np.full(image.shape,1e10),np.ones(image.shape),
                          np.zeros(image.shape,dtype=np.uint32),w,
                          release="dr9",region="south",brick="0965p780",band=b)
        grid=f.native_grid_fingerprint
        frames.append(f);scenes.append(SceneState("synthetic",grid,b,f.shape,()))
        psfs.append(PSFDescription("EMPIRICAL_KERNEL",kernels[i],{"fixture":"synthetic"}))
        identities.append((f.native_input_sha256,grid))
        settings.append(dict(star_icrs=(96.,78.),galaxy_icrs=anchor,
            core_arcsec=4,cutoff_arcsec=20,power=0,pivot_arcsec=8,
            galaxy_guard_arcsec=17,far_star_guard_arcsec=25,rays=(),
            block_size=3,minimum_block_fraction=.8,holdout_pa_deg=0,
            holdout_width_deg=90,buffer_deg=10))
    prepared=prepare_scene_candidates(frames,scenes,halo_settings=settings,
        compact_sources=[None,None],psfs=psfs,psf_identities=identities)
    assert all(len(s.layers)==0 for s in scenes)
    assert all(len(s.layers)==2 for s in prepared.proposed_scenes)
    with pytest.raises(WorkflowBlocked,match="MASK_REVIEW"):
        geometry_candidates(prepared,masks_confirmed=False,anchor_icrs=anchor,
            region=np.ones(frames[0].shape,dtype=bool),region_grid_id=identities[0][1],
            geometry_settings={})
    ann=dict(ring_edges=np.arange(0,14,2.),radial_nodes=np.linspace(0,120,6001),
             n_sectors=4,ring_min=.7,cumulative_min=.7,sector_min=.5)
    cov=dict(length_arcsec=5,mix=.5,minimum_control_tiles=8,clip_sigma=3.)
    tr=np.array([True]*5+[False])
    fs=dict(train_rings=tr,test_rings=~tr,initial=np.log(truth.ravel()),
            bounds=(np.log([.005,1.,.3]*2),np.log([.3,8.,1.5]*2)))
    args=dict(masks_confirmed=True,geometry_confirmed=True,geometry=geo,
              anchor_icrs=anchor,annular_settings=ann,covariance_settings=cov,fit_settings=fs)
    path=tmp_path/"joint.zip"
    with pytest.raises(WorkflowBlocked,match="GEOMETRY_REVIEW"):
        run_reviewed_annular_workflow(prepared,**dict(args,geometry_confirmed=False),output=path)
    assert not path.exists()
    out=run_reviewed_annular_workflow(prepared,**args,output=path,
        asymptotic_settings={"aperture_count":6},publication=True)
    assert "asymptotic_diagnostic" in out["fit"]
    assert out["stage"]=="DIAGNOSTIC_COMPLETE"
    assert out["fit"]["holdout"] is not None
    assert not out["science_adopted"] and not out["covariance_calibrated"]
    np.testing.assert_allclose(np.exp(out["fit"]["parameters"]).reshape(2,3),
                               truth,rtol=.15)
    with zipfile.ZipFile(path) as z:
        assert z.testzip() is None
        fit=json.loads(z.read("fit.json"))
        assert fit["workflow"]["masks_explicitly_reviewed"]
        assert "asymptotic.json" in z.namelist()
        import io
        import hashlib
        manifest=json.loads(z.read("manifest.json"))
        assert manifest["six_figure_diagnostic"]
        names=("01_stellar_halo","02_prepared_images","03_curve_of_growth",
               "04_surface_brightness","05_model_residual","06_radial_model")
        for name in names:
            assert z.read(name+".png").startswith(b"\x89PNG\r\n\x1a\n")
            assert b"<svg" in z.read(name+".svg")
        for name,digest in manifest["files_sha256"].items():
            assert hashlib.sha256(z.read(name)).hexdigest()==digest
        qa=json.loads(z.read("publication_qa.json"))
        assert not qa["masked_pixels_reconstructed"]
        assert max(r["max_model_ring_difference_nMgy"] for r in qa["ring_closure"])<1e-8
        with np.load(io.BytesIO(z.read("publication_arrays.npz")),allow_pickle=False) as a:
            assert np.isnan(a["band0_prepared"][~a["common_valid"]]).all()
            assert np.isnan(a["band0_residual"][~a["common_valid"]]).all()
            np.testing.assert_allclose(a["band0_prepared"]-a["band0_model"],
                                       a["band0_residual"])
            np.testing.assert_allclose(a["band0_observed_ring_flux"],
                                       out["records"][0]["observed_flux_nMgy"])
        assert not fit["workflow"]["full_geometry_containment_verified"]
    assert all(len(s.layers)==0 for s in scenes)
    assert tuple(f.native_input_sha256 for f in frames)==tuple(i[0] for i in identities)
    from lsbpipe.photometry.publication_products import publication_arrays
    from lsbpipe.photometry.native_annular_diagnostic import NativeCorrection
    original=prepared.corrections[0][0]
    wrong=NativeCorrection(original.kind,original.grid_id,original.input_sha256,
                           original.source,original.values+.01)
    with pytest.raises(ValueError,match="corrections differ"):
        publication_arrays(out["records"],out["covariance"]["covariance"],out["fit"],
            frames=prepared.frames,scenes=prepared.proposed_scenes,psfs=prepared.psfs,
            corrections_by_band=((wrong,)+prepared.corrections[0][1:],prepared.corrections[1]))
    # Force an invalid requested support rather than silently dropping rings.
    bad=dict(fs,train_rings=np.ones(6,dtype=bool))
    with pytest.raises(WorkflowBlocked,match="COVERAGE"):
        run_reviewed_annular_workflow(prepared,**dict(args,fit_settings=bad))

    from lsbpipe.photometry.pixel_injection import run_differential_sersic_injection
    profiles=truth.copy()
    profiles[:,0]*=.25
    recovered=run_differential_sersic_injection(frames,scenes,psfs=psfs,
        psf_identities=identities,halo_settings=settings,compact_sources=[None,None],
        injection_geometry=geo,injection_anchor_icrs=anchor,
        injection_profiles=profiles,workflow_settings=args)
    assert not recovered["science_adopted"]
    assert recovered["sky_halo_covariance_refitted"]
    for band in ("g","r"):
        row=recovered["bands"][band]
        assert np.all(row["expected_valid_ring_flux_nMgy"]>0)
        assert np.isfinite(row["ring_flux_bias_nMgy"]).all()
        assert row["fitted_model_total_change_nMgy"]>0
    assert tuple(f.native_input_sha256 for f in frames)==tuple(i[0] for i in identities)
