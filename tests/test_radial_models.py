import numpy as np
import pytest
from scipy.special import gammaincinv

from lsbpipe.photometry.radial_models import (
    EllipseNode, NestedEllipse, sersic_intensity, sersic_model_total,
    signed_curve_of_growth, conditional_plateau,
)


def test_constant_ellipse_and_sersic_normalization():
    node=EllipseNode(1.,2.,.6,0.)
    g=NestedEllipse(node,node,10.,70.)
    np.testing.assert_allclose(g.radius(np.array([1.,4.,1.]),
                                      np.array([2.,2.,7.])),[0.,5.,5.],atol=1e-12)
    m=sersic_model_total(g,.3,10.,1.)
    b=gammaincinv(2.,.5)
    assert m["model_total"]==pytest.approx(2*np.pi*.6*.3*100/b**2)
    assert m["model_half_light_arcsec"]==pytest.approx(10.)
    assert m["kind"]=="UNTRUNCATED_MODEL_EXTRAPOLATION"


def test_nested_ellipse_inverts_moving_contours():
    g=NestedEllipse(EllipseNode(0,0,.6,44),EllipseNode(.3,.5,.67,48),10,70)
    a=np.array([5.,20.,40.,69.,90.])
    x,y,q,pa=g.at(a)
    phi=.7
    e=x+a*np.cos(phi)*np.sin(pa)+a*q*np.sin(phi)*np.cos(pa)
    n=y+a*np.cos(phi)*np.cos(pa)-a*q*np.sin(phi)*np.sin(pa)
    np.testing.assert_allclose(g.radius(e,n),a,atol=1e-10)
    m=sersic_model_total(g,.3,28.,.6)
    assert m["model_half_light_arcsec"]>0
    assert m["model_half_light_arcsec"]!=pytest.approx(28.)
    with pytest.raises(ValueError):
        NestedEllipse(EllipseNode(0,0,.6,0),EllipseNode(100,0,.6,0),10,70)


def test_signed_growth_and_plateau_gate():
    f=[100.,-1.,2.,-3.]
    c=np.eye(4)*.01
    ok=[True,True,False,True]
    g=signed_curve_of_growth(f,c,ok)
    np.testing.assert_array_equal(g["flux"],[100.,99.,101.,98.])
    assert g["covariance"][1,3]==pytest.approx(.02)
    assert g["eligible"].tolist()==[True,True,False,False]
    p=conditional_plateau(f,c,ok,inner_count=1,outer_count=2,
                          fractional_tolerance=.02,normal_multiplier=1.96)
    assert p["conditional_equivalent"] and not p["adoptable"]
    p=conditional_plateau(f,c,ok,inner_count=1,outer_count=3,
                          fractional_tolerance=.02,normal_multiplier=1.96)
    assert not p["conditional_equivalent"]
