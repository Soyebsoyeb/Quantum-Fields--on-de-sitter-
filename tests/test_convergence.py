import numpy as np

from desitter_qft.analytic_modes import bunch_davies_mode
from desitter_qft.initial_conditions import adiabatic_data, exact_bd_data
from desitter_qft.mode_equations import ModeEquation
from desitter_qft.normalization import wronskian_residual
from desitter_qft.numerical_solver import solve_mode


def test_error_decreases_monotonically_with_tolerance_and_wronskian_follows():
    eq = ModeEquation(1.0, 0.0, 0.0)
    v0, dv0 = exact_bd_data(eq, -50.0)
    eta = np.linspace(-50, -0.1, 200)
    ref = bunch_davies_mode(eta, 1.0, 1.5)
    errs, ws = [], []
    for rt in (1e-5, 1e-7, 1e-9, 1e-11):
        s = solve_mode(eq, -50, -0.1, v0, dv0, rtol=rt, atol=rt * 1e-2, eta_eval=eta)
        errs.append(np.max(np.abs(s.v - ref) / np.abs(ref)))
        ws.append(wronskian_residual(s.v, s.dv).max())
    assert all(a > b for a, b in zip(errs, errs[1:]))
    assert all(a > b for a, b in zip(ws, ws[1:]))
    assert errs[-1] < 1e-9 and errs[0] / errs[-1] > 1e3


def test_adiabatic_vacuum_converges_to_bunch_davies_with_order_and_eta0():
    eq = ModeEquation(1.0, 0.8, 0.0)
    ref = abs(bunch_davies_mode(-0.5, 1.0, eq.nu)) ** 2
    def err(order, eta0):
        v0, dv0 = adiabatic_data(eq, eta0, order)
        s = solve_mode(eq, eta0, -0.5, v0, dv0, rtol=1e-12, atol=1e-14, n_out=2)
        return abs(abs(s.v[-1]) ** 2 / ref - 1)
    # earlier start -> smaller error ; higher order -> smaller error at fixed eta0
    assert err(2, -60) < err(2, -12) and err(0, -60) < err(0, -12)
    assert err(4, -20) < err(2, -20) < err(0, -20)
    assert err(4, -100) < 1e-10


def test_late_time_numerics_do_not_lose_normalisation_for_heavy_field():
    eq = ModeEquation(1.0, 1.7, 0.0)
    v0, dv0 = exact_bd_data(eq, -100.0)
    s = solve_mode(eq, -100.0, -1e-3, v0, dv0, rtol=1e-11, atol=1e-13)
    assert wronskian_residual(s.v, s.dv).max() < 1e-8
