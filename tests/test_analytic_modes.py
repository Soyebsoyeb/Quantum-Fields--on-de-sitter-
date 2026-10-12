import mpmath as mp
import numpy as np
import pytest

from desitter_qft.analytic_modes import (bunch_davies_mode, bunch_davies_mode_derivative, conformal_mode,
                                         minimal_massless_mode)
from desitter_qft.mode_equations import ModeEquation


@pytest.mark.parametrize("m,xi", [(0.0, 1 / 6), (0.0, 0.0), (0.8, 0.0), (1.7, 0.0), (0.3, 0.05)])
def test_hankel_mode_solves_ode_by_substitution(m, xi):
    eq = ModeEquation(1.0, m, xi)
    nu = eq.nu
    f = lambda e: 0.5 * mp.sqrt(mp.pi) * mp.exp(1j * mp.pi * (nu / 2 + 0.25)) * mp.sqrt(-e) * mp.hankel1(nu, -eq.k * e)
    for e in (-4.0, -1.0, -0.2):
        res = mp.diff(f, e, 2) + eq.omega2(e) * f(e)
        assert abs(res) < 1e-9


def test_derivative_matches_finite_difference():
    eta = np.array([-5.0, -1.3, -0.4])
    for nu in (0.5, 1.27, 1.5, 0.6j):
        h = 1e-5
        fd = (bunch_davies_mode(eta + h, 1.0, nu) - bunch_davies_mode(eta - h, 1.0, nu)) / (2 * h)
        np.testing.assert_allclose(bunch_davies_mode_derivative(eta, 1.0, nu), fd, rtol=1e-6, atol=1e-8)


def test_conformal_benchmark_matches_nu_half():
    eta = np.linspace(-10, -0.1, 40)
    v, dv = conformal_mode(eta, 1.4)
    np.testing.assert_allclose(v, bunch_davies_mode(eta, 1.4, 0.5), atol=1e-13)
    np.testing.assert_allclose(dv, bunch_davies_mode_derivative(eta, 1.4, 0.5), atol=1e-12)


def test_minimal_massless_benchmark_matches_nu_three_halves():
    eta = np.linspace(-10, -0.05, 40)
    v, dv = minimal_massless_mode(eta, 0.9)
    np.testing.assert_allclose(v, bunch_davies_mode(eta, 0.9, 1.5), atol=1e-12)
    np.testing.assert_allclose(dv, bunch_davies_mode_derivative(eta, 0.9, 1.5), atol=1e-11)


@pytest.mark.parametrize("nu", [0.3, 1.0, 1.4, 0.7j])
def test_bunch_davies_positive_frequency_early_time(nu):
    """v -> exp(-ik eta)/sqrt(2k) as k|eta| -> infinity (phase included)."""
    eta = -2000.0
    ref = np.exp(-1j * 1.0 * eta) / np.sqrt(2.0)
    assert abs(bunch_davies_mode(eta, 1.0, nu) - ref) < 5e-3 * abs(ref) * (1 + abs(nu))


def test_wrong_xi_fails_substitution_negative_control():
    """The closed-form conformal mode is NOT a solution when xi = 0: the check has teeth."""
    eq = ModeEquation(1.0, 0.0, 0.0)
    v, dv = conformal_mode(-1.0, 1.0)
    d2 = -v  # exp(-ik eta)'' = -k^2 v
    assert abs(eq.residual(-1.0, v, dv, d2)) > 0.1


def test_late_time_massless_power_is_scale_invariant():
    from desitter_qft.analytic_modes import field_mode
    for k in (0.1, 1.0, 10.0):
        eta = -1e-7 / k
        P = k**3 * abs(field_mode(eta, k, 1.5, H=1.0)) ** 2 / (2 * np.pi**2)
        assert abs(P - 1 / (4 * np.pi**2)) < 1e-6 / (4 * np.pi**2) * 10
