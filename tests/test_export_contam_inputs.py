import importlib.util
from pathlib import Path

import numpy as np
import pytest

spec=importlib.util.spec_from_file_location("export_contam_inputs",Path(__file__).resolve().parents[1]/"scripts"/"export_contam_inputs.py")
exporter=importlib.util.module_from_spec(spec)
spec.loader.exec_module(exporter)


def test_extract_masks_and_grid():
    a=np.ones((4,5))
    mask=np.zeros_like(a,dtype=bool)
    inputs=dict(image=a,iv=a,hard=mask,source_mask=mask,protected=mask,veto=mask,
                source_xy=(0,0),target_xy=(1,1),scale=.262,expanded_grid_fingerprint="a")
    probe=dict(sh1c_preview=dict(sky1_inputs=inputs,geo1={"reference":{}}),
               sh1b2_native_window={"wcs_fits_header":"header"})
    s=dict(expanded_grid_fingerprint="a",expanded_shape=a.shape,
           expanded_layers={"SH1C":a},expanded_safe_mask=mask,source_entries=[])
    sky={"expanded_sky":s}
    out=exporter.extract_arrays(probe,sky)
    assert out["arrays"]["reference_halo_residual"].shape==(4,5)
    assert out["arrays"]["hard"].dtype==np.bool_
    s["expanded_grid_fingerprint"]="b"
    with pytest.raises(ValueError,match="grid mismatch"):
        exporter.extract_arrays(probe,sky)


def test_cache_rejects_mismatch(tmp_path):
    path=tmp_path/"cache.zip"
    path.write_bytes(b"wrong")
    with pytest.raises(ValueError,match="does not match"):
        exporter.verified_cache(path,{"archive_sha256":"wrong","key":"x"})
