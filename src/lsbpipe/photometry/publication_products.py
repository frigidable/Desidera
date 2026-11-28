"""Six reproducible native-grid diagnostic figures, not adopted photometry."""
from __future__ import annotations

from hashlib import sha256
import io
import numpy as np
from scipy.signal import fftconvolve

from .contracts import NativeBandFrame, PSFDescription
from .native_annular_diagnostic import (
    NativeCorrection, annular_record_sha256, common_native_validity, numeric_sha256,
)
from .radial_models import NestedEllipse, EllipseNode, sersic_intensity
from .sky_masks import native_boolean_mask_sha256
from .uncertainty_budget import STATUS, _covariance


def publication_arrays(records, covariance, fit, *, frames, scenes,
                       corrections_by_band, psfs, max_work_pixels=2_000_000):
    """Rebuild exactly the signed data and PSF model used by the joint fit.

    No image interpolation over masks. The radial interpolation is the
    same model-only operation used by the adjoint annular operator.
    Returned display arrays use NaN outside common observed support.
    """
    from astropy.coordinates import SkyCoord
    from astropy import units as u

    records,frames,scenes,psfs=map(tuple,(records,frames,scenes,psfs))
    corrections_by_band=tuple(tuple(c) for c in corrections_by_band)
    if (not records or any(len(v)!=len(records) for v in
                          (frames,scenes,psfs,corrections_by_band))
            or fit.get("status")!=STATUS or not fit.get("success")
            or fit.get("boundary") or fit.get("conditional_parameter_covariance") is None):
        raise ValueError("Complete interior joint diagnostic and exact inputs required")
    if any("synthetic_injection_parent_sha256" in r for r in records):
        raise ValueError("Injected records cannot enter observational publication products")
    if (tuple(annular_record_sha256(r) for r in records)!=tuple(fit["record_sha256"])
            or numeric_sha256(covariance)!=fit["covariance_sha256"]):
        raise ValueError("Publication fit and record identity mismatch")
    if any(not isinstance(f,NativeBandFrame) for f in frames):
        raise ValueError("Validated native frames required")
    first=frames[0]
    if (not isinstance(max_work_pixels,int) or max_work_pixels<1
            or np.prod(first.shape)>max_work_pixels):
        raise ValueError("Publication work-pixel budget exceeded")
    valid=common_native_validity(frames,scenes)
    if not valid.any():
        raise ValueError("No common observed support")
    bands=tuple(fit["band_order"])
    if bands!=tuple(f.band for f in frames):
        raise ValueError("Publication band order mismatch")
    nr=len(records[0]["ring_edges_arcsec"])-1
    c=_covariance(covariance,nr*len(bands))
    z=np.asarray(fit["parameters"],dtype=float)
    pred=np.asarray(fit["prediction"],dtype=float)
    if (z.shape!=(3*len(bands),) or pred.shape!=(nr*len(bands),)
            or not np.isfinite(z).all() or not np.isfinite(pred).all()):
        raise ValueError("Invalid joint model parameters")
    g=records[0]["geometry"]
    geo=NestedEllipse(EllipseNode(**g["inner"]),EllipseNode(**g["outer"]),
                      g["transition_start"],g["transition_end"])
    center=SkyCoord(*records[0]["anchor_icrs"],unit=u.deg,frame="icrs")
    yy,xx=np.indices(first.shape)
    sky=first.wcs.pixel_to_world(xx,yy).icrs
    sep=center.separation(sky)
    if not np.isfinite(sep.deg).all() or np.max(sep.deg)>.5:
        raise ValueError("Small-field diagnostic outside supported WCS domain")
    pa=center.position_angle(sky).rad
    radius=geo.radius(sep.arcsec*np.sin(pa),sep.arcsec*np.cos(pa))
    edges=np.asarray(records[0]["ring_edges_arcsec"])
    labels=np.searchsorted(edges,radius,side="right")-1
    labels[labels>=nr]=-1
    ring=np.where(valid,labels,-1).astype(np.int32)
    inside=ring>=0
    if not inside.any():
        raise ValueError("No measured rings")
    arrays={"common_valid":valid,"ring_index":ring,"radius_arcsec":radius,
            "ring_edges_arcsec":edges,"covariance_nMgy2":c,
            "train_rings":np.asarray(fit["train_rings"]),
            "test_rings":np.asarray(fit["test_rings"])}
    rows=[]
    for i,(f,s,p,cs,r) in enumerate(zip(frames,scenes,psfs,corrections_by_band,records)):
        if (r["dataset_id"]!=f.dataset_id or r["native_input_sha256"]!=f.native_input_sha256
                or r["grid_id"]!=f.native_grid_fingerprint or r["scene_revision"]!=s.revision
                or r["valid_sha256"]!=native_boolean_mask_sha256(valid)
                or r["ring_index_sha256"]!=numeric_sha256(ring)
                or r["geometry"]!=g or r["anchor_icrs"]!=records[0]["anchor_icrs"]
                or not np.array_equal(r["ring_edges_arcsec"],edges)):
            raise ValueError("Publication scene or native ring identity mismatch")
        if (not isinstance(p,PSFDescription) or p.state!="EMPIRICAL_KERNEL"
                or sha256(np.asarray(p.kernel,dtype=np.float64).tobytes()).hexdigest()!=r["psf_sha256"]):
            raise ValueError("Publication PSF identity mismatch")
        if any(not isinstance(q,NativeCorrection) for q in cs):
            raise ValueError("Explicit native corrections required")
        if [{"kind":q.kind,"source":q.source,"sha256":q.sha256} for q in cs]!=r["corrections"]:
            raise ValueError("Publication corrections differ from measurement")
        image=np.array(f.image_nMgy,copy=True)
        before=image.copy()
        halo=np.zeros(f.shape)
        has_halo=False
        for q in cs:
            if (q.grid_id!=r["grid_id"] or q.input_sha256!=r["native_input_sha256"]
                    or q.values.shape!=f.shape or not np.isfinite(q.values[valid]).all()):
                raise ValueError("Publication correction support mismatch")
            image[valid]-=q.values[valid]
            if q.kind=="stellar_halo":
                halo=q.values
                has_halo=True
            else:
                before[valid]-=q.values[valid]
        nodes=np.asarray(r["radial_nodes_arcsec"])
        if radius.min()<nodes[0] or radius.max()>nodes[-1]:
            raise ValueError("Model nodes do not cover native grid")
        profile=sersic_intensity(nodes,*np.exp(z[3*i:3*i+3]))
        intrinsic=np.interp(radius,nodes,profile)
        model=fftconvolve(intrinsic,p.kernel,mode="same")
        obs=np.bincount(ring[inside],weights=image[inside],minlength=nr)
        mod=np.bincount(ring[inside],weights=model[inside],minlength=nr)
        sl=slice(i*nr,(i+1)*nr)
        if (not np.allclose(obs,r["observed_flux_nMgy"],rtol=1e-9,atol=1e-8)
                or not np.allclose(mod,pred[sl],rtol=1e-8,atol=1e-8)):
            raise ValueError("Publication raster fails observed/PSF model ring closure")
        prefix=f"band{i}_"
        for key,val in (("before_halo",before),("halo",halo),("prepared",image),
                        ("model",model),("residual",image-model)):
            arrays[prefix+key]=np.where(valid,val,np.nan)
        arrays[prefix+"observed_ring_flux"]=obs
        arrays[prefix+"model_ring_flux"]=mod
        arrays[prefix+"eligible"]=np.asarray(r["coverage"]["eligible"])
        arrays[prefix+"valid_area_arcsec2"]=np.asarray(r["valid_pixel_count"])*r["pixel_area_arcsec2"]
        rows.append({"band":f.band,"dataset_id":f.dataset_id,"has_halo":has_halo,
                     "max_observed_ring_difference_nMgy":float(np.max(np.abs(obs-r["observed_flux_nMgy"]))),
                     "max_model_ring_difference_nMgy":float(np.max(np.abs(mod-pred[sl])))})
    return arrays,rows


def render_publication_products(records, covariance, fit, **context):
    """Six figure pairs plus the exact plotted arrays; no pyplot or GUI."""
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    import json

    a,closure=publication_arrays(records,covariance,fit,**context)
    bands=tuple(fit["band_order"])
    nb=len(bands)
    edges=a["ring_edges_arcsec"]
    nr=len(edges)-1
    x=edges[1:]
    c=a["covariance_nMgy2"]
    valid=a["common_valid"]
    ypix,xpix=np.nonzero(a["ring_index"]>=0)
    pad=max(3,int(.1*max(np.ptp(ypix),np.ptp(xpix))))
    extent=(max(0,int(xpix.min())-pad)-.5,min(valid.shape[1],int(xpix.max())+pad+1)-.5,
            max(0,int(ypix.min())-pad)-.5,min(valid.shape[0],int(ypix.max())+pad+1)-.5)
    files={}
    def figure(name,cols,rows=1,wide=4.0):
        fig=Figure(figsize=(wide*cols,3.3*rows+.6),layout="constrained")
        FigureCanvasAgg(fig)
        fig.suptitle(name+" | DIAGNOSTIC_NOT_ADOPTED",fontsize=10)
        return fig,fig.subplots(rows,cols,squeeze=False)
    def save(fig,name):
        for ext in ("png","svg"):
            b=io.BytesIO()
            fig.savefig(b,format=ext,dpi=160)
            files[f"{name}.{ext}"]=b.getvalue()
        fig.clear()
    def limits(*values,symmetric=False):
        v=np.concatenate([np.asarray(k)[np.isfinite(k)] for k in values])
        if not len(v):
            return -1.,1.
        lo,hi=np.percentile(v,[2,98])
        if symmetric:
            h=max(abs(lo),abs(hi),1e-15)
            return -h,h
        if hi<=lo:
            return lo-max(abs(lo)*.01,1e-15),hi+max(abs(hi)*.01,1e-15)
        return lo,hi
    def map_panel(fig,ax,image,title,lim,*,crop=False,contours=False,cmap="gray"):
        im=ax.imshow(image,origin="lower",interpolation="nearest",cmap=cmap,
                     vmin=lim[0],vmax=lim[1])
        ax.set_title(title,fontsize=9)
        ax.set_xlabel("Native x [pixel]")
        ax.set_ylabel("Native y [pixel]")
        if contours:
            levels=edges[1:][edges[1:]<np.nanmax(a["radius_arcsec"])]
            if len(levels):
                ax.contour(a["radius_arcsec"],levels=levels,colors="tab:orange",linewidths=.55)
        if crop:
            ax.set_xlim(extent[:2]);ax.set_ylim(extent[2:])
        fig.colorbar(im,ax=ax,shrink=.75,label="nMgy / native pixel")
    fig,axes=figure("Stellar contamination",3,nb)
    for i,b in enumerate(bands):
        p=f"band{i}_"
        lim=limits(a[p+"before_halo"],a[p+"prepared"])
        map_panel(fig,axes[i,0],a[p+"before_halo"],b+" | before halo removal",lim)
        map_panel(fig,axes[i,1],a[p+"halo"],b+" | halo model on valid support",
                  limits(a[p+"halo"]))
        map_panel(fig,axes[i,2],a[p+"prepared"],b+" | after correction",lim)
    save(fig,"01_stellar_halo")
    fig,axes=figure("Prepared native observations",nb)
    for i,b in enumerate(bands):
        im=a[f"band{i}_prepared"]
        map_panel(fig,axes[0,i],im,b+" | common support",limits(im),crop=True,contours=True)
    save(fig,"02_prepared_images")
    fig,axes=figure("Signed observed curve of growth",1,wide=7)
    ax=axes[0,0]
    l=np.tril(np.ones((nr,nr)))
    for i,b in enumerate(bands):
        p=f"band{i}_";sl=slice(i*nr,(i+1)*nr)
        obs=l@a[p+"observed_ring_flux"];model=l@a[p+"model_ring_flux"]
        err=np.sqrt(np.maximum(0,np.diag(l@c[sl,sl]@l.T)))
        good=a[p+"eligible"]
        line=ax.errorbar(x[good],obs[good],yerr=err[good],fmt="o",label=b+" observed")
        color=line[0].get_color()
        ax.plot(x,model,"--",color=color,label=b+" PSF model, same support")
        ax.plot(x[~good],obs[~good],"x",color=color)
    ax.set(xlabel="Semimajor radius [arcsec]",ylabel="Cumulative signed flux [nMgy]")
    ax.legend(fontsize=8);save(fig,"03_curve_of_growth")
    fig,axes=figure("Observed surface brightness (positive flux only)",1,wide=7)
    ax=axes[0,0]
    for i,b in enumerate(bands):
        p=f"band{i}_";sl=slice(i*nr,(i+1)*nr)
        area=a[p+"valid_area_arcsec2"]
        flux=a[p+"observed_ring_flux"]
        err=np.sqrt(np.maximum(0,np.diag(c[sl,sl])))
        good=a[p+"eligible"]&(flux>0)&(area>0)
        robust=good&(flux>err)
        color=f"C{i%10}"
        mu=lambda f,ar:22.5-2.5*np.log10(f/ar)
        if robust.any():
            value=mu(flux[robust],area[robust])
            lo=value-mu(flux[robust]+err[robust],area[robust])
            hi=mu(flux[robust]-err[robust],area[robust])-value
            ax.errorbar(x[robust],value,yerr=[lo,hi],fmt="o",color=color,label=b+" positive, finite interval")
        weak=good&~robust
        if weak.any():
            ax.plot(x[weak],mu(flux[weak],area[weak]),"o",mfc="none",color=color,
                    label=b+" positive, lower flux interval crosses zero")
        ax.plot([],[]," ",label=f"{b}: {int((~good).sum())} rings omitted (nonpositive or ineligible)")
    ax.set(xlabel="Semimajor radius [arcsec]",ylabel="AB mag / arcsec² (before extinction)")
    ax.invert_yaxis();ax.legend(fontsize=7);save(fig,"04_surface_brightness")
    fig,axes=figure("Observed image, PSF model and signed residual",3,nb)
    for i,b in enumerate(bands):
        p=f"band{i}_";lim=limits(a[p+"prepared"],a[p+"model"])
        map_panel(fig,axes[i,0],a[p+"prepared"],b+" | observed",lim,crop=True)
        map_panel(fig,axes[i,1],a[p+"model"],b+" | PSF model",lim,crop=True)
        map_panel(fig,axes[i,2],a[p+"residual"],b+" | data minus model",
                  limits(a[p+"residual"],symmetric=True),crop=True,cmap="RdBu_r")
    save(fig,"05_model_residual")
    fig,axes=figure("Signed radial intensity and residual",1,2,wide=7)
    for i,b in enumerate(bands):
        p=f"band{i}_";sl=slice(i*nr,(i+1)*nr)
        area=a[p+"valid_area_arcsec2"]
        obs=a[p+"observed_ring_flux"]/area
        model=a[p+"model_ring_flux"]/area
        err=np.sqrt(np.maximum(0,np.diag(c[sl,sl])))/area
        good=a[p+"eligible"]
        color=f"C{i%10}"
        axes[0,0].errorbar(x[good],obs[good],yerr=err[good],fmt="o",color=color,label=b+" observed")
        axes[0,0].plot(x,model,"--",color=color,label=b+" model")
        axes[1,0].plot(x[good],(obs-model)[good],"o-",color=color,label=b)
        for ax,v in ((axes[0,0],obs),(axes[1,0],obs-model)):
            ax.plot(x[~good],v[~good],"x",color=color)
    axes[0,0].set_ylabel("Mean signed intensity [nMgy / arcsec²]")
    axes[1,0].set(xlabel="Semimajor radius [arcsec]",ylabel="Raw residual [nMgy / arcsec²]")
    axes[1,0].axhline(0,color="gray",linewidth=.8)
    for ax in axes.ravel(): ax.legend(fontsize=8)
    save(fig,"06_radial_model")
    b=io.BytesIO()
    np.savez_compressed(b,**a)
    files["publication_arrays.npz"]=b.getvalue()
    files["publication_qa.json"]=json.dumps({"status":STATUS,"band_order":bands,
        "ring_closure":closure,"record_sha256":fit["record_sha256"],
        "covariance_sha256":fit["covariance_sha256"],
        "masked_pixels_reconstructed":False,"covariance_calibrated":False,
        "full_geometry_containment_verified":False,
        "residual_kind":"RAW_DATA_MINUS_MODEL_NOT_CONDITIONAL_HOLDOUT"},
        indent=2,allow_nan=False).encode("utf-8")
    return files
