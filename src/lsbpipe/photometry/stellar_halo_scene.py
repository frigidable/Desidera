"""Source-bound SH1-B2 crop of a verified original native brick.

No padding, reprojection, project mutation or network.  The diagnostic
grid is deliberately separate from the displayed galaxy Scene grid.
"""
from __future__ import annotations
import math
import numpy as np
from astropy.wcs.utils import proj_plane_pixel_scales


def expanded_halo_cutout(frame, scene, source_icrs, *, radius_arcsec=180.):
    from lsbpipe.scene import grid_fingerprint
    from lsbpipe.survey import legacy_validity
    from lsbpipe.photometry.contracts import _immutable
    key=tuple(scene["source_quartet_key"])
    if frame.dataset_id!="/".join(key):
        raise ValueError("SH1-B2 native quartet identity mismatch")
    evidence=frame.source_evidence
    hashes=dict(evidence["local_hashes"])
    if (set(hashes)!={"image","invvar","nexp","maskbits"} or
            hashes!=scene["original_fits_sha256_local"]):
        raise ValueError("SH1-B2 original FITS hashes differ from loaded Scene")
    if dict(evidence["hdu_indices"])!=dict(scene["hdu_indices"]):
        raise ValueError("SH1-B2 original FITS HDUs differ from loaded Scene")
    wcs=frame.wcs
    scale=np.asarray(proj_plane_pixel_scales(wcs.celestial))*3600.
    if (scale.shape!=(2,) or not np.isfinite(scale).all() or min(scale)<=0
            or max(scale)/min(scale)>1.01 or not 115<=radius_arcsec<=180):
        raise ValueError("SH1-B2 invalid native scale or requested radius")
    # Verify that the original display grid maps by an integer translation
    # to the native brick, not merely by matching shape or file hash.
    ny,nx=scene["pixel_shape"]
    points=np.array([[0.,0.],[nx-1.,0.],[0.,ny-1.],[nx-1.,ny-1.],[(nx-1)/2,(ny-1)/2]])
    world=scene["wcs"].all_pix2world(points,0)
    offsets=wcs.all_world2pix(world,0)-points
    translation=np.rint(offsets[0])
    if not np.isfinite(offsets).all() or np.max(np.abs(offsets-translation))>1e-3:
        raise ValueError("SH1-B2 Scene is not an integer native-brick translation")
    sx,sy=map(float,wcs.all_world2pix(*source_icrs,0))
    h,w=frame.shape
    if not (0<=sx<w and 0<=sy<h):
        raise ValueError("SH1-B2 star lies outside original brick")
    half=math.ceil(radius_arcsec/min(scale))+4
    requested=[math.floor(sx)-half,math.floor(sy)-half,
               math.floor(sx)+half+1,math.floor(sy)+half+1]
    x0,y0,x1,y1=max(0,requested[0]),max(0,requested[1]),min(w,requested[2]),min(h,requested[3])
    sl=np.s_[y0:y1,x0:x1]
    native={"image":frame.image_nMgy,"invvar":frame.invvar_nMgy2,
            "nexp":frame.nexp,"maskbits":frame.maskbits}
    # Immutable, detached crop; keep the same Scene validity semantics.
    data={k:_immutable(v[sl],v.dtype) for k,v in native.items()}
    cutwcs=wcs.slice(sl,numpy_order=True)
    shape=data["image"].shape
    valid,_=legacy_validity(data["image"],data["invvar"],data["maskbits"],key[3])
    hard=_immutable(~valid|~np.isfinite(data["nexp"])|(data["nexp"]<=0),np.bool_)
    grid=grid_fingerprint(cutwcs,shape)
    out={"data":data,"wcs":cutwcs,"pixel_shape":shape,"native_hard":hard,
         "band":key[3],"source":scene["source"],"project_id":scene["project_id"],
         "source_quartet_key":list(key),"grid_fingerprint":grid,
         "original_fits_sha256_local":hashes}
    meta={"schema":"desidera.sh1b2-native-brick-window/1",
          "requested_radius_arcsec":radius_arcsec,"full_brick_shape":list(frame.shape),
          "requested_bounds_xyxy_exclusive":requested,
          "actual_bounds_xyxy_exclusive":[x0,y0,x1,y1],
          "clipped_by_brick":requested!=[x0,y0,x1,y1],
          "native_shape":list(shape),"grid_fingerprint":grid,
          "original_scene_grid_fingerprint":scene["grid_fingerprint"],
          "scene_origin_full_brick_xy":translation.astype(int).tolist(),
          "source_xy_full_brick":[sx,sy],
          "edge_distances_arcsec":{"negative_x":float(sx*scale[0]),
                                    "positive_x":float((w-1-sx)*scale[0]),
                                    "negative_y":float(sy*scale[1]),
                                    "positive_y":float((h-1-sy)*scale[1])},
          "wcs_fits_header":dict(cutwcs.to_header()),
          "no_padding":True,"no_reprojection":True,"original_fits_sha256_local":hashes}
    return out,meta
