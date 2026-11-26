import numpy as np
import pytest
from lsbpipe.photometry.stellar_halo import measure_saturated_halo_sectors


def test_saturated_sector_profiles_preserve_signed_flux_and_target_veto():
    n=321
    yy,xx=np.indices((n,n))
    sx,sy=160.,110.
    r=np.hypot(xx-sx,yy-sy)
    image=12/(1+r)**1.5-.1
    hard=(r<10)|(abs(yy-sy)<2)
    image[hard]=1e6
    protected=(xx-160)**2+(yy-220)**2<60**2
    image[protected]+=20
    source=np.zeros_like(hard)
    source[30:50,150:170]=True
    veto=np.zeros_like(hard)
    before=image.copy()
    result=measure_saturated_halo_sectors(image,np.ones_like(image),hard,
        source_xy=(sx,sy),target_xy=(160,220),pixel_scale_arcsec=.262,
        source_mask=source,protected=protected,unresolved_sky_veto=veto,
        neighbor_centres=[(160,40)])
    assert result["status"]=="SH1_SECTORS_MEASURED_NOT_ADOPTED"
    assert not result["sky_subtracted"] and not result["halo_fitted"]
    assert not result["compact_psf_used"] and not result["catalog_flux_used"]
    assert np.array_equal(image,before)
    rows=result["rows"]
    facing=next(x for x in rows if x["sector_index"]==0 and x["r_inner_arcsec"]==25.)
    assert facing["descriptive"]["median_nMgy_pix"]>10
    assert facing["outside_protection_and_veto"]["n"]==0
    away=next(x for x in rows if x["sector_index"]==4 and x["r_inner_arcsec"]==16.)
    assert away["outside_protection_and_veto"]["median_nMgy_pix"]<0
    assert away["adequate_for_descriptive_comparison"]
    edge=next(x for x in rows if x["sector_index"]==4 and x["r_inner_arcsec"]==100.)
    assert edge["in_frame_pixels"]==0 and not edge["adequate_for_descriptive_comparison"]


def test_saturated_sector_rejects_unaligned_masks():
    a=np.zeros((10,10))
    with pytest.raises(ValueError,match="aligned"):
        measure_saturated_halo_sectors(a,a,a,source_xy=(1,1),target_xy=(5,5),
            pixel_scale_arcsec=.262,source_mask=a,protected=a[:2],
            unresolved_sky_veto=a)
