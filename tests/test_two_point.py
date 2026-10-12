import numpy as np
import pytest

from desitter_qft.two_point import (coincident_phi2_pointsplit, invariant_ode_residual, invariant_Z,
                                    phi2_cutoff, phi2_cutoff_minimal_massless_analytic,
                                    short_distance_coefficients, wightman_closed_form,
                                    wightman_closed_form_equal_time, wightman_closed_form_u, wightman_mode_sum)
from desitter_qft.infrared import free_massive_stochastic_phi2, stochastic_equilibrium_moments
from desitter_qft.mode_equations import nu as nu_of
import mpmath as mp


def test_conformal_case_matches_flat_space_green_function_unequal_times():
    eta, etap, r, H = -0.7, -0.4, 1.3, 1.0
    exact = H**2 * eta * etap / (4 * np.pi**2 * (r**2 - (eta - etap) ** 2))
    assert abs(wightman_closed_form(invariant_Z(eta, etap, r), 0.5, H) - exact) < 1e-14


@pytest.mark.parametrize("nu", [0.3, 1.1, 1.45, 0.5j])
def test_closed_form_satisfies_invariant_equation(nu):
    for z in (-0.8, 0.0, 0.6):
        res, g = invariant_ode_residual(z, nu)
        assert res / g < 1e-20


@pytest.mark.parametrize("nu", [0.2, 1.0, 1.3])
def test_mode_sum_equals_closed_form(nu):
    for r in (0.1, 1.0, 5.0):
        a, b = wightman_mode_sum(-1.0, r, nu), wightman_closed_form_equal_time(-1.0, r, nu)
        assert abs(a / b - 1) < 1e-9


def test_mode_sum_imaginary_nu():
    a, b = wightman_mode_sum(-1.0, 1.0, 0.8j), wightman_closed_form_equal_time(-1.0, 1.0, 0.8j)
    assert abs(a / b - 1) < 1e-8


def test_scaling_with_eta_de_sitter_invariance_of_equal_time_function():
    """G depends on r/|eta| only through Z, times nothing else: G(eta, r) = G(1, r/|eta|)."""
    for nu in (0.7, 1.4):
        assert abs(wightman_closed_form_equal_time(-3.0, 6.0, nu) - wightman_closed_form_equal_time(-1.0, 2.0, nu)) < 1e-15


def test_hadamard_leading_singularity_is_flat_space_form():
    H, eta = 1.0, -1.0
    a = 1 / (H * abs(eta))
    for nu in (0.5, 1.0, 1.4, 0.7j):
        r = 1e-4
        val = 4 * np.pi**2 * (a * r) ** 2 * wightman_closed_form_equal_time(eta, r, nu, H)
        assert abs(val - 1) < 1e-6


@pytest.mark.parametrize("nu", [0.5, 1.0, 1.4, 0.8j])
def test_short_distance_expansion_coefficients_independent_check(nu):
    """DLMF 15.8.10-derived (A, B, C) reproduce the exact function up to O(u ln u)."""
    A, B, C = short_distance_coefficients(nu)
    for u in (1e-5, 1e-7):
        G = float(wightman_closed_form_u(u, nu))
        assert abs(G - (A / u + B * np.log(u) + C)) < 20 * u * abs(np.log(u)) * (1 + abs(B))


def test_phi2_cutoff_numeric_matches_closed_form_nu_three_halves():
    for kuv in (5.0, 100.0):
        a = phi2_cutoff(-0.5, kuv, 1e-3, 1.5)
        b = phi2_cutoff_minimal_massless_analytic(-0.5, kuv, 1e-3)
        assert abs(a / b - 1) < 1e-9


def test_cutoff_is_not_a_renormalisation_quadratic_divergence():
    a1 = phi2_cutoff(-0.5, 200.0, 1e-3, 1.5)
    a2 = phi2_cutoff(-0.5, 2000.0, 1e-3, 1.5)
    assert abs(a2 / a1 - 100) < 5   # grows ~ Lambda^2: no finite limit


def test_light_field_limit_matches_stochastic_equilibrium_variance():
    """Independent IR benchmark: point-split <phi^2> -> 3H^4/(8 pi^2 m^2) (leading in m/H) and equals the
    Fokker-Planck Gaussian equilibrium.  The finite-part differs by O(H^2): scheme-dependent constants."""
    for m in (0.02, 0.05):
        pe = coincident_phi2_pointsplit(nu_of(m, 0.0).real)
        sy = free_massive_stochastic_phi2(m)
        assert abs(pe / sy - 1) < 20 * m**2
        assert abs(stochastic_equilibrium_moments(lambda p: 0.5 * m**2 * p**2)[2] / sy - 1) < 1e-8


def test_pointsplit_value_depends_on_renormalisation_scale():
    nu = 1.2
    assert abs(coincident_phi2_pointsplit(nu, mu=1.0) - coincident_phi2_pointsplit(nu, mu=2.0)) > 1e-4
    nu_c = 0.5  # conformal massless: B = 0, scale-independent
    assert abs(coincident_phi2_pointsplit(nu_c, mu=1.0) - coincident_phi2_pointsplit(nu_c, mu=2.0)) < 1e-15
