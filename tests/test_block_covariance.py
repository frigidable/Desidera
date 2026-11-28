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
