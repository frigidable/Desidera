import numpy as np
import pytest

from lsbpipe.photometry.pixel_injection import inject_native_sersic
from lsbpipe.photometry.radial_models import EllipseNode, NestedEllipse
from lsbpipe.photometry.contracts import NativeBandFrame, PSFDescription
from lsbpipe.photometry.native_annular_diagnostic import (
    NativeCorrection, prepare_native_annuli,
)
from lsbpipe.photometry.annular_export import export_annular_diagnostic
from lsbpipe.scene import SceneState


def test_native_injection_preserves_parent_and_excludes_synthetic_export(tmp_path):
    WCS=pytest.importorskip("astropy.wcs").WCS
    w=WCS(naxis=2)
    w.wcs.crpix=[21.,21.]
    w.wcs.crval=[96.,78.]
    w.wcs.cdelt=[-.5/3600,.5/3600]
    w.wcs.ctype=["RA---TAN","DEC--TAN"]
    shape=(41,41)
    image=np.full(shape,-.002)
    image[0,0]=np.nan
    extra=np.zeros(shape,dtype=bool);extra[20,21]=True
    frame=NativeBandFrame(image,np.full(shape,1e6),np.ones(shape),
        np.zeros(shape,dtype=np.uint32),w,release="dr9",region="south",
        brick="0965p780",band="r",extra_exclusion=extra)
    kernel=np.array([[1,2,1],[2,4,2],[1,2,1]],float)/16
    psf=PSFDescription("EMPIRICAL_KERNEL",kernel,{"fixture":"synthetic"})
    node=EllipseNode(0.,0.,.8,30.)
    geometry=NestedEllipse(node,node,2.,8.)
    identity=(frame.native_input_sha256,frame.native_grid_fingerprint)
    injected=inject_native_sersic(frame,psf,psf_identity=identity,
        geometry=geometry,anchor_icrs=(96,78),profile=(.03,1.5,.7))
    child=injected.frame
    assert child.native_input_sha256!=identity[0]
    assert frame.native_input_sha256==identity[0]
    assert child.native_grid_fingerprint==identity[1]
    assert np.isnan(child.image_nMgy[0,0])
    np.testing.assert_array_equal(child.valid,frame.valid)
    np.testing.assert_array_equal(child.invvar_nMgy2,frame.invvar_nMgy2)
    np.testing.assert_array_equal(child.nexp,frame.nexp)
    np.testing.assert_array_equal(child.maskbits,frame.maskbits)
    np.testing.assert_array_equal(frame.image_nMgy,image)
    np.testing.assert_allclose(child.image_nMgy[frame.valid]-image[frame.valid],
                               injected.signal_nMgy_per_pixel[frame.valid],atol=1e-15)
    assert not injected.signal_nMgy_per_pixel.flags.writeable
    assert not injected.metadata["additional_noise_simulated"]
    scene=SceneState("synthetic",identity[1],"r",shape,())
    sky=NativeCorrection("residual_sky",identity[1],child.native_input_sha256,
                          "synthetic known residual sky",np.zeros(shape))
    record=prepare_native_annuli(child,scene,geometry,
        anchor_icrs=(96,78),ring_edges=[0,2,4,6],radial_nodes=np.linspace(0,30,1501),
        corrections=(sky,),psf=psf,psf_input_sha256=child.native_input_sha256,
        psf_grid_id=identity[1],n_sectors=4,ring_min=.5,cumulative_min=.5,sector_min=.4)
    assert record["synthetic_injection_parent_sha256"]==identity[0]
    with pytest.raises(ValueError,match="Synthetic injection"):
        export_annular_diagnostic(tmp_path/"observational.zip",[record],np.eye(3),{})
    assert not (tmp_path/"observational.zip").exists()
