import numpy as np
import sympy as sp

from desitter_qft.experiments import tensor_experiment
from desitter_qft.tensor_modes import (derive_tensor_equation, dimensionless_tensor_spectrum,
                                       tensor_spectrum_late_time)


def test_tensor_equation_follows_from_linearised_einstein_equations():
    r = derive_tensor_equation()
    assert r["background_residual"] == 0
    assert sp.simplify(r["difference_from_claim"]) == 0


def test_late_time_tensor_spectrum_value():
    assert abs(dimensionless_tensor_spectrum(-1e-8, 1.0, 1.0, 1.0) / tensor_spectrum_late_time() - 1) < 1e-8
    assert abs(tensor_spectrum_late_time(2.0, 3.0) - 2 * 4.0 / (np.pi**2 * 9.0)) < 1e-15


def test_numerical_tensor_spectrum_agrees_with_analytic():
    d = tensor_experiment()
    assert np.max(np.abs(d["Dt2_num"] / d["Dt2_ana"] - 1)) < 1e-6
    # The mode grows ~ 1/|eta| towards eta = -1e-3, so |v||v'| ~ 1e9 and a relative integration error of 1e-12
    # corresponds to an absolute Wronskian residual of order 1e-7..1e-6 (accuracy loss is reported, not rescaled away).
    assert d["wronskian"]["max_residual"] < 1e-5
