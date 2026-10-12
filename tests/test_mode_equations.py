import numpy as np
import sympy as sp

from desitter_qft import mode_equations as me
from desitter_qft import scalar_action as sa


def test_klein_gordon_from_action():
    assert sp.simplify(sa.euler_lagrange_vs_klein_gordon()[2]) == 0


def test_cosmic_time_kg_has_3H_friction():
    diff, R = sa.cosmic_time_equation_difference()
    assert sp.simplify(diff) == 0 and sp.simplify(R) == 12 * sp.Symbol("H", positive=True) ** 2


def test_conformal_equation_two_independent_derivations():
    d, R, *_ = sa.conformal_time_equation()
    assert sp.simplify(d) == 0
    assert sp.simplify(sa.cosmic_to_conformal_chain_rule()) == 0


def test_rescaled_equation_and_dS_reduction():
    assert sp.simplify(me.derive_rescaled_mode_equation()["difference"]) == 0
    assert sp.simplify(me.de_sitter_reduction()) == 0


def test_conformally_coupled_massless_limit_cancels_time_dependence():
    eta = np.linspace(-10, -0.1, 50)
    w2 = me.omega_squared(eta, k=1.3, m=0.0, xi=1.0 / 6.0, H=2.0)
    np.testing.assert_allclose(w2, 1.3**2, rtol=1e-14)


def test_minimal_massless_effective_potential():
    eta = np.linspace(-10, -0.1, 50)
    np.testing.assert_allclose(me.omega_squared(eta, 1.0, 0.0, 0.0, 1.0), 1.0 - 2.0 / eta**2, rtol=1e-14)


def test_nu_regimes():
    assert me.nu_squared(0.0, 0.0) == 2.25 and me.nu(0.0, 0.0) == 1.5
    assert abs(me.nu(0.0, 1 / 6) - 0.5) < 1e-14
    assert me.nu(2.0, 0.0).imag > 0 and me.nu(2.0, 0.0).real == 0   # principal series
    assert "tachyonic" in me.regime(me.nu_squared(0.0, -0.1))      # nu > 3/2
    assert "imaginary" in me.regime(me.nu_squared(2.0, 0.0))
    assert "massless-minimal" in me.regime(2.25)
