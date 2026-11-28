import numpy as np
import pytest
from scipy.special import gammaincinv, gammaincc

from lsbpipe.photometry.asymptotic_diagnostic import linearized_parameter_response
from lsbpipe.photometry.radial_models import EllipseNode, NestedEllipse, sersic_model_total


def test_parameter_response_retains_training_measurement_correlation():
    j=np.array([[1.,0.],[1.,1.],[1.,2.],[1.,3.]])
    c=np.array([[2.,.3,.2,.4],[.3,1.,.1,.2],[.2,.1,1.5,.3],[.4,.2,.3,1.2]])
    train=np.array([True,True,True,False])
    a=linearized_parameter_response(j,c,train)
    ci=np.linalg.inv(c[:3,:3])
    expected=np.linalg.solve(j[:3].T@ci@j[:3],j[:3].T@ci)
    np.testing.assert_allclose(a[:,:3],expected)
    np.testing.assert_array_equal(a[:,3],0.)
    np.testing.assert_allclose(a@c@a.T,np.linalg.inv(j[:3].T@ci@j[:3]))
    # Hybrid residual retains the cross term between observations and fit.
    s=np.ones(4)
    g=s-j.sum(axis=0)@a
    assert g@c@g.T==pytest.approx(s@c@s.T+j.sum(axis=0)@a@c@a.T@j.sum(axis=0)
                                -2*s@c@a.T@j.sum(axis=0))
    with pytest.raises(ValueError):
        linearized_parameter_response(np.ones((4,2)),c,train)


def test_intrinsic_tail_is_stable_and_closes_for_variable_geometry():
    node=EllipseNode(0,0,.7,30.)
    g=NestedEllipse(node,node,4.,16.)
    result=sersic_model_total(g,1.,2.,1.,aperture_radius=40.)
    expected=gammaincc(2,gammaincinv(2,.5)*20.)
    assert result["intrinsic_outer_tail"]>0
    assert result["intrinsic_outer_tail"]/result["model_total"]==pytest.approx(expected,rel=1e-10,abs=0)
    g=NestedEllipse(node,EllipseNode(.1,.1,.73,31.),4.,16.)
    for a in (0.,2.,8.,16.,40.):
        r=sersic_model_total(g,1.,5.,.8,aperture_radius=a)
        assert r["intrinsic_enclosed_flux"]+r["intrinsic_outer_tail"]==pytest.approx(r["model_total"],rel=1e-9)
        assert r["intrinsic_outer_tail"]>=0
