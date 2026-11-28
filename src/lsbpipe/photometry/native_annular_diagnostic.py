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
from .radial_models import NestedEllipse, sersic_intensity, sersic_model_total, EllipseNode
from .annular_inference import (annular_support, adjoint_annular_operator,
                                fit_conditional_gls, conditional_holdout)
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
    # This conservative diagnostic never restores a native saturation or
    # ALLMASK flag even if the frame's exploratory policy allowed it.
    hard_valid=frame.valid & ~frame.mask_reasons[f"SATUR_{frame.band.upper()}"]
    hard_valid &= ~frame.mask_reasons[f"ALLMASK_{frame.band.upper()}"]
    valid=evaluate_mask(scene,frame.image_nMgy,frame.invvar_nMgy2,grid_id=grid,
                        operation=Operation.PHOTOMETRY,observation_valid=hard_valid).valid
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
            "pixel_area_arcsec2":frame.local_pixel_area_arcsec2,
            "valid_sha256":native_boolean_mask_sha256(valid),
            "corrections":[{"kind":c.kind,"source":c.source,"sha256":c.sha256} for c in corrections],
            "psf_sha256":sha256(np.asarray(psf.kernel,dtype=np.float64).tobytes()).hexdigest(),
            "anchor_icrs":[ra,dec],"geometry":asdict(geometry),"ring_edges_arcsec":edges.copy(),
            "coverage_thresholds":{"n_sectors":n_sectors,"ring":ring_min,"cumulative":cumulative_min,"sector":sector_min},
            "psf_shape":list(psf.kernel.shape),
            "radial_nodes_arcsec":np.asarray(radial_nodes,dtype=float).copy(),
            "coverage":support,"observed_flux_nMgy":flux,
            "valid_pixel_count":masks.sum(axis=(1,2)),
            "formal_diagonal_variance_nMgy2":formal,"operator":op,
            "operator_ring_indices":np.flatnonzero(measured),
            "covariance_calibrated":False,"full_geometry_containment_verified":False}


def fit_common_support_sersic(records, covariance, *, band_order,
                              covariance_source, train_rings, test_rings,
                              initial, bounds):
    """Joint conditional fit on explicitly common support, never science adoption.

    Records are ordered by band_order, which also specifies covariance block
    order. Each parameter triple is log(I0 per pixel), log(scale arcsec),
    log(Sersic index). Same native WCS and exact same observed support are
    required for this particular joint diagnostic. The caller supplies the
    full conditional covariance, not just diagonal IVAR.
    """
    records=tuple(records)
    bands=tuple(band_order)
    if (not records or len(records)!=len(bands) or len(set(bands))!=len(bands)
            or not isinstance(covariance_source,str) or not covariance_source.strip()):
        raise ValueError("Explicit band order and covariance provenance required")
    first=records[0]
    n=len(first["ring_edges_arcsec"])-1
    for rec,b in zip(records,bands):
        if rec["status"]!=STATUS or rec["dataset_id"].split("/")[-1]!=b:
            raise ValueError("Record status or band order mismatch")
        if (rec["dataset_id"].rsplit("/",1)[0]!=first["dataset_id"].rsplit("/",1)[0]
                or rec["grid_id"]!=first["grid_id"]
                or rec["valid_sha256"]!=first["valid_sha256"]
                or rec["geometry"]!=first["geometry"]
                or rec["anchor_icrs"]!=first["anchor_icrs"]
                or not np.array_equal(rec["ring_edges_arcsec"],first["ring_edges_arcsec"])):
            raise ValueError("Joint fit requires identical native geometry and support")
        if not np.array_equal(rec["operator_ring_indices"],np.arange(n)):
            raise ValueError("Missing ring support; prepare a smaller explicitly declared range")
    tr,te=np.asarray(train_rings),np.asarray(test_rings)
    if (tr.shape!=(n,) or te.shape!=(n,) or tr.dtype.kind!="b" or te.dtype.kind!="b"
            or not tr.any() or np.any(tr&te)):
        raise ValueError("Explicit disjoint ring selections required")
    for rec in records:
        if np.any((tr|te)&~rec["coverage"]["eligible"]):
            raise ValueError("Training/holdout outside contiguous coverage gate")
    p=np.asarray(initial,dtype=float)
    if p.shape!=(3*len(bands),):
        raise ValueError("Three log parameters per ordered band required")
    d=np.concatenate([rec["observed_flux_nMgy"] for rec in records])
    def predictor(z):
        result=[]
        for i,rec in enumerate(records):
            amp,re,idx=np.exp(z[3*i:3*i+3])
            profile=sersic_intensity(rec["radial_nodes_arcsec"],amp,re,idx)
            result.append(rec["operator"]@profile)
        return np.concatenate(result)
    train=np.tile(tr,len(bands))
    test=np.tile(te,len(bands))
    result=fit_conditional_gls(d,covariance,predictor,p,bounds,train=train)
    result["band_order"]=bands
    result["covariance_source"]=covariance_source
    result["train_rings"]=tr.copy()
    result["test_rings"]=te.copy()
    result["record_sha256"]=tuple(annular_record_sha256(rec) for rec in records)
    result["covariance_sha256"]=numeric_sha256(covariance)
    result["model_extrapolations"]=None
    result["holdout"]=None
    result["full_geometry_containment_verified"]=False
    if not result["success"]:
        return result
    z=result["parameters"]
    geo=first["geometry"]
    geometry=NestedEllipse(EllipseNode(**geo["inner"]),EllipseNode(**geo["outer"]),
                           geo["transition_start"],geo["transition_end"])
    models={}
    for i,(rec,b) in enumerate(zip(records,bands)):
        amp,re,idx=np.exp(z[3*i:3*i+3])
        models[b]=sersic_model_total(geometry,amp/rec["pixel_area_arcsec2"],re,idx)
        models[b]["index"]=float(idx)
    result["model_extrapolations"]=models
    if test.any() and result["conditional_parameter_covariance"] is not None:
        eps=1e-5
        jac=np.column_stack([(predictor(z+v*eps)-predictor(z-v*eps))/(2*eps)
                             for v in np.eye(len(z))])
        result["holdout"]=conditional_holdout(
            d,result["prediction"],covariance,jac,
            result["conditional_parameter_covariance"],train=train,test=test)
    return result


def numeric_sha256(values):
    """Canonical float64 array digest including shape."""
    a=np.ascontiguousarray(values,dtype="<f8")
    h=sha256(str(a.shape).encode("ascii")+b":")
    h.update(a.tobytes())
    return h.hexdigest()


def annular_record_sha256(rec):
    """Bind fit/export to the same measured ring arrays and provenance."""
    import json
    fields=("status","dataset_id","grid_id","native_input_sha256",
            "scene_revision","valid_sha256","corrections","psf_sha256",
            "anchor_icrs","geometry","coverage_thresholds","pixel_area_arcsec2")
    meta={k:rec[k] for k in fields}
    h=sha256(json.dumps(meta,sort_keys=True,allow_nan=False).encode("utf-8"))
    for name in ("ring_edges_arcsec","radial_nodes_arcsec","observed_flux_nMgy",
                 "formal_diagonal_variance_nMgy2","operator","operator_ring_indices",
                 "valid_pixel_count"):
        h.update(name.encode("ascii"))
        h.update(numeric_sha256(rec[name]).encode("ascii"))
    for name in ("ring_fraction","cumulative_fraction","sector_fraction","eligible"):
        h.update(name.encode("ascii"))
        h.update(numeric_sha256(rec["coverage"][name]).encode("ascii"))
    return h.hexdigest()
