"""Explicit nested-ellipse hypotheses, radial models and signed growth QA.

No geometry optimization or selection is implicit. A model total is an
extrapolation, not a measured integral or an accepted asymptotic flux.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import gammaincinv, gammainc, gammaln

from .uncertainty_budget import _covariance, STATUS


@dataclass(frozen=True)
class EllipseNode:
    east_arcsec: float
    north_arcsec: float
    q: float
    pa_deg: float

    def __post_init__(self):
        for name in ("east_arcsec","north_arcsec","q","pa_deg"):
            value=getattr(self,name)
            if isinstance(value,(bool,np.bool_,str)) or not math.isfinite(float(value)):
                raise ValueError("Finite ellipse parameters required")
            object.__setattr__(self,name,float(value))
        if not 0 < self.q <= 1:
            raise ValueError("Invalid ellipse axis ratio")


@dataclass(frozen=True)
class NestedEllipse:
    """Smoothstep transition; angles explicitly unwrapped, east of north."""
    inner: EllipseNode
    outer: EllipseNode
    transition_start: float
    transition_end: float

    def __post_init__(self):
        if not isinstance(self.inner,EllipseNode) or not isinstance(self.outer,EllipseNode):
            raise ValueError("Explicit inner and outer nodes required")
        lo,hi=float(self.transition_start),float(self.transition_end)
        if not math.isfinite(lo) or not math.isfinite(hi) or not 0 <= lo < hi:
            raise ValueError("Invalid transition radii")
        object.__setattr__(self,"transition_start",lo)
        object.__setattr__(self,"transition_end",hi)
        # A conservative sufficient bound: the outward normal speed of
        # a homothetic ellipse is >= q. Subtract maximal center, shape
        # and rotation speeds. No sampled-only non-crossing guarantee.
        a,b=self.inner,self.outer
        dcenter=math.hypot(b.east_arcsec-a.east_arcsec,b.north_arcsec-a.north_arcsec)
        dq=abs(b.q-a.q)
        dpa=abs(math.radians(b.pa_deg-a.pa_deg))
        bound=1.5/(hi-lo)*(dcenter+hi*(dq+dpa))
        if min(a.q,b.q)<=bound:
            raise ValueError("Cannot certify nested contours with conservative derivative bound")

    def at(self,a):
        a=np.asarray(a,dtype=float)
        if not np.isfinite(a).all() or np.any(a<0):
            raise ValueError("Invalid radius")
        z=np.clip((a-self.transition_start)/(self.transition_end-self.transition_start),0,1)
        t=z*z*(3-2*z)
        i,o=self.inner,self.outer
        return tuple(v+t*(w-v) for v,w in zip(
            (i.east_arcsec,i.north_arcsec,i.q,np.deg2rad(i.pa_deg)),
            (o.east_arcsec,o.north_arcsec,o.q,np.deg2rad(o.pa_deg))))

    def area(self,a):
        return np.pi*np.asarray(a,dtype=float)**2*self.at(a)[2]

    def radius(self,east,north):
        """Implicit pixel-center radius on a caller-validated tangent plane."""
        e,n=map(lambda v:np.asarray(v,dtype=float),(east,north))
        if e.shape!=n.shape or not np.isfinite(e).all() or not np.isfinite(n).all():
            raise ValueError("Invalid tangent-plane coordinates")
        def rho(a):
            x,y,q,pa=self.at(a)
            dx,dy=e-x,n-y
            return np.hypot(dx*np.sin(pa)+dy*np.cos(pa),
                            (dx*np.cos(pa)-dy*np.sin(pa))/q)
        low=np.zeros_like(e)
        high=np.maximum(self.transition_end,rho(np.full_like(e,self.transition_end)))+1
        for _ in range(52):
            mid=(low+high)/2
            outside=rho(mid)>mid
            low=np.where(outside,mid,low)
            high=np.where(outside,high,mid)
        return (low+high)/2


def sersic_intensity(radius, central_intensity, scale_radius, index):
    """Central-intensity convention; scale radius need not be true half-light."""
    r=np.asarray(radius,dtype=float)
    i0,re,n=map(float,(central_intensity,scale_radius,index))
    if not np.isfinite(r).all() or np.any(r<0) or not all(math.isfinite(v) and v>0 for v in (i0,re,n)):
        raise ValueError("Invalid positive Sersic parameters")
    return i0*np.exp(-gammaincinv(2*n,.5)*(r/re)**(1/n))


def sersic_model_total(geometry, central_intensity, scale_radius, index):
    """Model-only integration using d(pi q(a) a^2)/da and analytic outer tail.

    central_intensity is per square arcsec, all radii in arcsec. The result
    assumes the entire untruncated model, including unobserved regions.
    """
    if not isinstance(geometry,NestedEllipse):
        raise ValueError("Explicit nested geometry required")
    i0,re,n=map(float,(central_intensity,scale_radius,index))
    sersic_intensity([0],i0,re,n)
    bn=gammaincinv(2*n,.5)
    lo,hi=geometry.transition_start,geometry.transition_end
    qi,qo=geometry.inner.q,geometry.outer.q
    norm=2*np.pi*i0*n*re**2*np.exp(gammaln(2*n)-2*n*np.log(bn))
    if not np.isfinite(norm) or norm<=0:
        raise ValueError("Nonfinite model total")
    def primitive(a):
        return norm*gammainc(2*n,bn*(a/re)**(1/n))
    def integrand(a):
        z=(a-lo)/(hi-lo)
        q=qi+(qo-qi)*z*z*(3-2*z)
        dq=(qo-qi)*6*z*(1-z)/(hi-lo)
        return float(sersic_intensity(a,i0,re,n)*np.pi*(2*a*q+a*a*dq))
    mid=quad(integrand,lo,hi,epsabs=norm*1e-10,epsrel=1e-9)[0]
    total=qi*primitive(lo)+mid+qo*(norm-primitive(hi))
    def enclosed(a):
        if a<=lo:
            return qi*primitive(a)
        if a<hi:
            return qi*primitive(lo)+quad(integrand,lo,a,epsabs=norm*1e-10,epsrel=1e-9)[0]
        return qi*primitive(lo)+mid+qo*(primitive(a)-primitive(hi))
    upper=max(hi,re)
    for _ in range(100):
        if enclosed(upper)>=total/2:
            break
        upper*=2
    else:
        raise ValueError("Cannot bracket model half-light radius")
    half=brentq(lambda a:enclosed(a)-total/2,0.,upper)
    return {"status":STATUS,"kind":"UNTRUNCATED_MODEL_EXTRAPOLATION",
            "model_total":float(total),"model_half_light_arcsec":float(half),
            "profile_scale_arcsec":re}


def signed_curve_of_growth(ring_flux, covariance, eligible):
    """Signed observed annulus sums; never fill holes or enforce monotonicity."""
    f=np.asarray(ring_flux,dtype=float)
    ok=np.asarray(eligible)
    if f.ndim!=1 or not f.size or not np.isfinite(f).all() or ok.shape!=f.shape or ok.dtype.kind!="b":
        raise ValueError("Invalid signed flux/coverage arrays")
    c=_covariance(covariance,len(f))
    l=np.tril(np.ones((len(f),len(f))))
    cc=l@c@l.T
    return {"status":STATUS,"kind":"OBSERVED_VALID_PIXEL_FLUX",
            "flux":np.cumsum(f),"covariance":cc,
            "eligible":np.logical_and.accumulate(ok)}


def conditional_plateau(ring_flux,covariance,eligible,*,inner_count,outer_count,
                        fractional_tolerance,normal_multiplier):
    """Conservative conditional check of equivalence to negligible growth.

    Not an asymptotic total or a calibrated confidence statement. Both
    the tolerance and multiplier must be predeclared; the reference is the
    observed signed inner flux, not an inferred total.
    """
    g=signed_curve_of_growth(ring_flux,covariance,eligible)
    f=g["flux"]; c=g["covariance"]
    i,j=inner_count-1,outer_count-1
    tol,z=float(fractional_tolerance),float(normal_multiplier)
    if not 0<=i<j<len(f) or not 0<tol<1 or not math.isfinite(z) or z<=0 or f[i]<=0:
        raise ValueError("Invalid plateau interval, tolerance or reference flux")
    delta=f[j]-f[i]
    # Test both sides of the equivalence interval including covariance of
    # the random reference flux: +/- delta - tol*Finner <= 0.
    upper_sigma=np.sqrt(max(0.,c[j,j]+(1+tol)**2*c[i,i]-2*(1+tol)*c[i,j]))
    lower_sigma=np.sqrt(max(0.,c[j,j]+(1-tol)**2*c[i,i]-2*(1-tol)*c[i,j]))
    upper=delta-tol*f[i]+z*upper_sigma
    lower=-delta-tol*f[i]+z*lower_sigma
    return {"status":STATUS,"delta_flux":float(delta),"tolerance_flux":float(tol*f[i]),
            "eligible":bool(g["eligible"][j]),"conditional_equivalent":bool(g["eligible"][j] and upper<=0 and lower<=0),
            "adoptable":False}
