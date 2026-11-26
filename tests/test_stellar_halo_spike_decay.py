import numpy as np
import pytest
from lsbpipe.photometry.stellar_halo_spikes import coords
from lsbpipe.photometry.stellar_halo_spike_decay import (
    basis,prepare_blocks,evaluate,positive_cv,injection_trials,decay_experiment)


def test_decay_grid_recovers_across_missing_blocks():
    t,q=coords((901,401),(200.,450.),(200.,650.),.25)
    a=.01*basis(-t,q,20.,.786)+.001+.0002*q
    safe=(t<0)&~((t<-27)&(t>-42))
    blocks,coverage=prepare_blocks(a,safe,t,q,.25)
    f=evaluate(blocks)
    assert any(r["status"]=="UNMEASURED" for r in coverage)
    assert f["best"]["length_arcsec"]==20.
    assert f["best"]["sigma_arcsec"]==.786
    assert f["best"]["amplitude_at_17arcsec_nMgy_pix"]==pytest.approx(.01)
    assert positive_cv(f)
    for v in injection_trials(blocks,f):
        assert v["recovered"]["length_arcsec"]==v["injected_length_arcsec"]
        assert v["recovered"]["sigma_arcsec"]==v["injected_sigma_arcsec"]


def test_no_signal_and_protected_target_do_not_generate_model():
    t,q=coords((901,401),(200.,450.),(200.,650.),.25)
    a=np.zeros(t.shape); a[t>0]=100.
    result=decay_experiment(a,t<0,source_xy=(200.,450.),target_xy=(200.,650.),
        scale=.25,preview={},refit_residual=a)
    assert not result["positive_cv"]
    assert not result["layers"]
    assert not result["science_adopted"]
