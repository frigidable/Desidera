import numpy as np
from lsbpipe.photometry.stellar_halo_spikes import coords
from lsbpipe.photometry.stellar_halo_spike_extent import scan_axis,contiguous_support,compact_template


def test_finite_ray_not_extended_by_isolated_outer_signal():
    shape=(901,901); scale=.25
    t,q=coords(shape,(450.,450.),(450.,650.),scale)
    # away ray only, two inner consecutive bins, isolated far signal.
    ray=(t<=-12)&(t>-32)&(abs(q)<1.5)
    far=(t<=-62)&(t>-67)&(abs(q)<1.5)
    a=np.where(ray|far,.01,0.)
    safe=t<0
    rows=scan_axis(a,safe,t,q,scale)
    support=contiguous_support(rows)
    assert support["detected_outer_edge_arcsec"]==32.
    assert support["endpoint_status"]=="NEXT_BIN_NOT_DETECTED"
    assert any(r["inner_arcsec"]==62. and r["status"]=="DETECTED" for r in rows)
    m=compact_template(t,q,rows,support)
    assert m[530,450]>0
    assert m[710,450]==0
    assert m[530,470]==0


def test_gap_and_null_control_are_not_physical_endpoints():
    shape=(901,901); t,q=coords(shape,(450.,450.),(450.,650.),.25)
    a=np.where((t<=-12)&(t>-42)&(abs(q)<1.5),.01,0.)
    safe=(t<0)&~((t<=-27)&(t>-32))
    rows=scan_axis(a,safe,t,q,.25)
    s=contiguous_support(rows)
    assert s["detected_outer_edge_arcsec"]==27.
    assert s["endpoint_status"]=="COVERAGE_GAP"
    null=scan_axis(np.zeros(shape),t<0,t,q,.25)
    assert contiguous_support(null)["status"]=="NO_CONTIGUOUS_INNER_DETECTION"
