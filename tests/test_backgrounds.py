import numpy as np
import pytest

from desitter_qft.backgrounds import DeSitter


def test_scale_factor_relations():
    dS = DeSitter(0.7)
    t = np.linspace(-1, 5, 50)
    eta = dS.conformal_time(t)
    assert np.all(eta < 0)
    # a(t) = a(eta(t)) : two independent forms of the scale factor agree
    np.testing.assert_allclose(dS.scale_factor_cosmic(t), dS.scale_factor_conformal(eta), rtol=1e-13)
    # a' / a^2 = H in conformal time: da/deta = a^2 H
    h = 1e-6
    da = (dS.scale_factor_conformal(eta + h) - dS.scale_factor_conformal(eta - h)) / (2 * h)
    np.testing.assert_allclose(da / dS.scale_factor_conformal(eta) ** 2, 0.7, rtol=1e-6)
    np.testing.assert_allclose(dS.cosmic_time(eta), t, atol=1e-12)


def test_scalar_curvature_value_and_dimension():
    assert DeSitter(2.0).ricci_scalar == 12 * 4.0
    assert DeSitter(2.0).cosmological_constant == 3 * 4.0


def test_domain_errors():
    with pytest.raises(ValueError):
        DeSitter(-1.0)
    with pytest.raises(ValueError):
        DeSitter(1.0).scale_factor_conformal(0.1)


def test_horizon_exit():
    dS = DeSitter(1.3)
    k = 5.0
    eta = dS.horizon_exit_time(k)
    assert abs(k / dS.scale_factor_conformal(eta) - dS.H) < 1e-12  # physical wavenumber = H
