import numpy as np
import pytest

from lsbpipe.photometry.uncertainty_budget import (
    conditional_magnitude_budget, linear_sigma, scenario_envelope,
)


def test_shared_extinction_and_correlated_zeropoint():
    stat = np.array([[.04**2, .001], [.001, .05**2]])
    zp = np.full((2, 2), .02**2)
    b = conditional_magnitude_budget(stat, [3.214, 2.165], .01, zp)
    assert b['status'] == 'DIAGNOSTIC_NOT_ADOPTED'
    assert linear_sigma(b['assumed_cal_cov'], [1, -1]) == pytest.approx(.01049)
    assert linear_sigma(b['conditional_plus_assumed_cal_cov'], [1, -1])**2 == pytest.approx(.04**2+.05**2-.002+.01049**2)
    assert np.array_equal(stat, b['stat_cov'])


@pytest.mark.parametrize('bad', [[[1, 2], [0, 1]], [[1, 2], [2, 1]], [[1, 0], [0, float('nan')]]])
def test_invalid_covariance(bad):
    with pytest.raises(ValueError):
        conditional_magnitude_budget(bad, [3, 2], .01, np.eye(2))


def test_envelope_explicit_eligibility_and_signed_offsets():
    result = scenario_envelope([16.1, 16.4, 99], 16.2, [True, True, False])
    assert result['n'] == 2
    assert result['delta_min'] == pytest.approx(-.1)
    assert result['delta_max'] == pytest.approx(.2)
    assert result['kind'] == 'TESTED_EXTREMA_NOT_CONFIDENCE_INTERVAL'
    with pytest.raises(ValueError):
        scenario_envelope([1], 1, [False])
    with pytest.raises(ValueError):
        scenario_envelope([1], 1, [1])
