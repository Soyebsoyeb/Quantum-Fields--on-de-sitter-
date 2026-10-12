import numpy as np
import pytest

from desitter_qft.analytic_modes import bunch_davies_mode, bunch_davies_mode_derivative
from desitter_qft.initial_conditions import adiabatic_data, exact_bd_data
from desitter_qft.mode_equations import ModeEquation
from desitter_qft.normalization import check_normalization, wronskian, wronskian_residual
from desitter_qft.numerical_solver import solve_mode


@pytest.mark.parametrize("nu", [0.5, 1.0, 1.27, 1.5, 0.8j])
def test_analytic_modes_satisfy_canonical_wronskian(nu):
    eta = np.linspace(-60, -0.01, 300)
    W = wronskian(bunch_davies_mode(eta, 1.0, nu), bunch_davies_mode_derivative(eta, 1.0, nu))
    np.testing.assert_allclose(W, 1j, atol=1e-10)


def test_wronskian_is_conserved_and_detects_wrong_normalisation():
    eta = np.linspace(-30, -0.1, 50)
    v = bunch_davies_mode(eta, 1.0, 1.5)
    dv = bunch_davies_mode_derivative(eta, 1.0, 1.5)
    assert wronskian_residual(v, dv).max() < 1e-12
    assert wronskian_residual(2.0 * v, 2.0 * dv).min() > 2.9   # |4i - i| = 3: normalisation error is flagged


def test_numerical_wronskian_report_and_no_rescaling():
    eq = ModeEquation(1.0, 0.8, 0.0)
    v0, dv0 = exact_bd_data(eq, -60.0)
    s = solve_mode(eq, -60.0, -0.05, v0, dv0, rtol=1e-10, atol=1e-12)
    rep = check_normalization(s.eta, s.v, s.dv, rtol=1e-10, atol=1e-12, method="DOP853")
    assert rep.passes and rep.max_residual < 1e-8
    assert "no rescaling" in rep.note
    assert rep.initial_residual < 1e-12   # exact initial data are canonically normalised


@pytest.mark.parametrize("order", [0, 2, 4])
def test_adiabatic_initial_data_are_normalised_to_stated_order(order):
    eq = ModeEquation(1.0, 0.0, 0.0)
    v0, dv0 = adiabatic_data(eq, -300.0, order)
    # residual of the Wronskian at the initial point decreases with adiabatic order
    res = abs(wronskian(v0, dv0) - 1j)
    assert res < [1e-6, 1e-9, 1e-10][order // 2 if order else 0] * 1e3


def test_solver_rejects_bad_interval():
    eq = ModeEquation(1.0)
    with pytest.raises(ValueError):
        solve_mode(eq, -1.0, 0.5, 1 + 0j, 0j)
