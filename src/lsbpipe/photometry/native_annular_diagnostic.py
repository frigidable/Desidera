"""Native-frame orchestration of the experimental annular diagnostic.

Explicit raster corrections and PSF identity are required. This function
does not select a background, halo, geometry or covariance automatically.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from hashlib import sha256
import math
import numpy as np

from lsbpipe.scene import SceneState, evaluate_mask, Operation
from .contracts import NativeBandFrame, PSFDescription
from .radial_models import NestedEllipse
from .annular_inference import annular_support, adjoint_annular_operator
from .sky_masks import native_boolean_mask_sha256
from .uncertainty_budget import STATUS


@dataclass(frozen=True)
class NativeCorrection:
    """One explicitly identified correction in nMgy/native pixel."""
    kind: str
    grid_id: str
    input_sha256: str
    source: str
    values: np.ndarray

    def __post_init__(self):
        if self.kind not in ("residual_sky","stellar_halo","compact_sources"):
            raise ValueError("Unknown or overlapping correction kind")
        if not all(isinstance(v,str) and v for v in (self.grid_id,self.input_sha256,self.source)):
            raise ValueError("Correction provenance required")
        v=np.asarray(self.values,dtype=np.float64)
        if v.ndim!=2 or not v.size:
            raise ValueError("Expected a native correction raster")
        object.__setattr__(self,"values",np.frombuffer(v.tobytes(),dtype=v.dtype).reshape(v.shape))

    @property
    def sha256(self):
        return sha256(self.values.tobytes()).hexdigest()


def prepare_native_annuli(frame, scene, geometry, *, anchor_icrs, ring_edges,
                          radial_nodes, corrections, psf, psf_input_sha256,
                          psf_grid_id, n_sectors, ring_min, cumulative_min,
                          sector_min, max_work_pixels=2_000_000, row_chunk=128):
    """Build observed signed fluxes and a finite native PSF response matrix.

    The tangent plane is a small-field sky-offset approximation. Coverage
    fractions refer ONLY to the rasterized in-frame geometry. The caller
    must additionally verify no geometric truncation at frame boundaries.
    No interpolation/reprojection of pixels, no correction of missing light.
    """
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    if not isinstance(frame,NativeBandFrame) or not isinstance(scene,SceneState):
        raise ValueError("Validated native frame and explicit scene required")
    if not isinstance(geometry,NestedEllipse):
        raise ValueError("Explicit nested geometry hypothesis required")
    grid=frame.native_grid_fingerprint
    digest=frame.native_input_sha256
    if scene.band!=frame.band or scene.grid_id!=grid or scene.shape!=frame.shape:
        raise ValueError("Scene identity mismatch")
    if (not isinstance(psf,PSFDescription) or psf.state!="EMPIRICAL_KERNEL" or
            psf_input_sha256!=digest or psf_grid_id!=grid):
        raise ValueError("Explicit same-frame finite empirical PSF required")
    corrections=tuple(corrections)
    if not corrections or any(not isinstance(c,NativeCorrection) for c in corrections):
        raise ValueError("Explicit corrections including residual sky required")
    kinds=[c.kind for c in corrections]
    if "residual_sky" not in kinds or len(set(kinds))!=len(kinds):
        raise ValueError("Require one residual sky and no duplicate correction kinds")
    valid=evaluate_mask(scene,frame.image_nMgy,frame.invvar_nMgy2,grid_id=grid,
                        operation=Operation.PHOTOMETRY,observation_valid=frame.valid).valid
    image=np.array(frame.image_nMgy,copy=True)
    for c in corrections:
        if c.grid_id!=grid or c.input_sha256!=digest or c.values.shape!=frame.shape:
            raise ValueError("Correction input/grid mismatch")
        if not np.isfinite(c.values[valid]).all():
            raise ValueError("Nonfinite correction on usable pixels")
        image[valid]-=c.values[valid]
    if not np.isfinite(image[valid]).all():
        raise ValueError("Invalid corrected flux")
    edges=np.asarray(ring_edges,dtype=float)
    if (edges.ndim!=1 or len(edges)<2 or edges[0]!=0 or
            not np.isfinite(edges).all() or np.any(np.diff(edges)<=0)):
        raise ValueError("Increasing ring edges beginning at zero required")
    ra,dec=map(float,anchor_icrs)
    if not math.isfinite(ra) or not math.isfinite(dec) or not 0<=ra<360 or not -90<=dec<=90:
        raise ValueError("Invalid ICRS anchor")
    if (not isinstance(max_work_pixels,int) or isinstance(max_work_pixels,bool) or
            max_work_pixels<1 or np.prod(frame.shape)>max_work_pixels):
        raise ValueError("Native diagnostic exceeds explicit work-pixel budget; use a validated cutout or raise budget")
    if not isinstance(row_chunk,int) or isinstance(row_chunk,bool) or row_chunk<1:
        raise ValueError("Invalid row chunk")
    center=SkyCoord(ra=ra*u.deg,dec=dec*u.deg,frame="icrs")
    if not isinstance(n_sectors,int) or n_sectors<1:
        raise ValueError("Explicit sector count required")
    radius=np.empty(frame.shape,dtype=float)
    sectors=np.empty(frame.shape,dtype=np.int32)
    wcs=frame.wcs
    for start in range(0,frame.shape[0],row_chunk):
        stop=min(start+row_chunk,frame.shape[0])
        yy,xx=np.mgrid[start:stop,:frame.shape[1]]
        coords=wcs.pixel_to_world(xx,yy).icrs
        sep=center.separation(coords)
        if np.max(sep.deg)>.5:
            raise ValueError("Small-field tangent-plane diagnostic restricted to 0.5 degree")
        pa=center.position_angle(coords).rad
        east,north=sep.arcsec*np.sin(pa),sep.arcsec*np.cos(pa)
        a=geometry.radius(east,north)
        cx,cy,q,theta=geometry.at(a)
        dx,dy=east-cx,north-cy
        phi=np.mod(np.arctan2((dx*np.cos(theta)-dy*np.sin(theta))/q,
                             dx*np.sin(theta)+dy*np.cos(theta)),2*np.pi)
        radius[start:stop]=a
        sectors[start:stop]=np.minimum(n_sectors-1,np.floor(phi/(2*np.pi)*n_sectors).astype(int))
    labels=np.searchsorted(edges,radius,side="right")-1
    nr=len(edges)-1
    labels[labels>=nr]=-1
    support=annular_support(labels,sectors,valid,n_rings=nr,n_sectors=n_sectors,
                            ring_min=ring_min,cumulative_min=cumulative_min,sector_min=sector_min)
    # Keep any rings with no usable pixels explicitly missing (not zero).
    masks=np.array([valid&(labels==j) for j in range(nr)])
    measured=masks.any(axis=(1,2))
    flux=np.full(nr,np.nan)
    formal=np.full(nr,np.nan)
    for j in np.flatnonzero(measured):
        m=masks[j]
        flux[j]=image[m].sum()
        formal[j]=(1/frame.invvar_nMgy2[m]).sum()
    op=None
    if measured.any():
        op=adjoint_annular_operator(radius,radial_nodes,masks[measured],psf.kernel,
                                    grid_id=grid,mask_grid_id=scene.grid_id)
    return {"status":STATUS,"dataset_id":frame.dataset_id,"grid_id":grid,
            "native_input_sha256":digest,"scene_revision":scene.revision,
            "valid_sha256":native_boolean_mask_sha256(valid),
            "corrections":[{"kind":c.kind,"source":c.source,"sha256":c.sha256} for c in corrections],
            "psf_sha256":sha256(np.asarray(psf.kernel,dtype=np.float64).tobytes()).hexdigest(),
            "anchor_icrs":[ra,dec],"geometry":asdict(geometry),"ring_edges_arcsec":edges.copy(),
            "coverage_thresholds":{"n_sectors":n_sectors,"ring":ring_min,"cumulative":cumulative_min,"sector":sector_min},
            "psf_shape":list(psf.kernel.shape),
            "radial_nodes_arcsec":np.asarray(radial_nodes,dtype=float).copy(),
            "coverage":support,"observed_flux_nMgy":flux,
            "formal_diagonal_variance_nMgy2":formal,"operator":op,
            "operator_ring_indices":np.flatnonzero(measured),
            "covariance_calibrated":False,"full_geometry_containment_verified":False}
