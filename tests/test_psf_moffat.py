import numpy as np
import pytest
from lsbpipe.photometry.psf_moffat import pixel_moffat, fit_joint_moffat


def test_joint_moffat_recovers_noiseless_shape_with_distinct_centres_flux_and_sky():
    n=41
    yy,xx=np.indices((n,n),dtype=float)
    shape=np.array([np.log(5.2),np.log(4.5),.12,2.7])
    entries=[]
    for i in range(5):
        cx,cy=20+.13*i,20-.08*i
        dx,dy=.15-.06*i,-.12+.03*i
        x,y=xx-cx,yy-cy
        im=(70+17*i)*pixel_moffat(x-dx,y-dy,shape)+.02+.0003*x-.0007*y
        w=np.ones_like(im)*100
        w[1:3,1:4]=0
        entries.append((str(i),{"image":im,"weights":w,"center_in_cutout":(cx,cy)}))
    k,qa=fit_joint_moffat(entries,5.,n,18.)
    assert k.shape==(n,n) and np.all(k>=0)
    assert k.sum()==pytest.approx(1.)
    assert qa["shape_parameters"]==pytest.approx(shape,abs=.015)
    assert qa["training_only"] and not qa["science_adopted"]
    assert qa["training_stars"][0]["centroid_offset_xy_pix"]==pytest.approx([.15,-.12],abs=.015)
    assert qa["training_stars"][0]["analytic_flux_nMgy_model_dependent"]==pytest.approx(70,rel=.015)


def test_joint_moffat_rejects_unusable_input():
    with pytest.raises(ValueError,match="four training"):
        fit_joint_moffat([],5.,41,18.)
