"""Rescaled mode equation  v'' + [k^2 + a^2 m^2 + (6 xi - 1) a''/a] v = 0.

Symbolic derivation (category 2) and numerical coefficient functions (category 3 support).
In de Sitter, a''/a = 2/eta^2, so the mode equation becomes
    v'' + [k^2 - (nu^2 - 1/4)/eta^2] v = 0,    nu^2 = 9/4 - m^2/H^2 - 12 xi.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import sympy as sp


def derive_rescaled_mode_equation():
    """Substitute phi = v/a into  phi'' + 2(a'/a)phi' + k^2 phi + a^2 (m^2 + xi R) phi = 0, R = 6a''/a^3.

    Returns dict with the residual operator and the claimed target; `difference` must be 0.
    """
    eta = sp.Symbol("eta", real=True)
    k, m, xi = sp.symbols("k m xi", real=True)
    a = sp.Function("a", positive=True)(eta)
    v = sp.Function("v")(eta)
    phi = v / a
    R = 6 * sp.diff(a, eta, 2) / a**3
    lhs = sp.diff(phi, eta, 2) + 2 * sp.diff(a, eta) / a * sp.diff(phi, eta) + k**2 * phi + a**2 * (m**2 + xi * R) * phi
    target = (sp.diff(v, eta, 2) + (k**2 + a**2 * m**2 + (6 * xi - 1) * sp.diff(a, eta, 2) / a) * v) / a
    return {"lhs": sp.simplify(lhs), "target": target, "difference": sp.simplify(lhs - target)}


def de_sitter_reduction():
    """Show (6 xi - 1) a''/a + a^2 m^2 = -(nu^2 - 1/4)/eta^2 for a = -1/(H eta)."""
    eta = sp.Symbol("eta", negative=True)
    H, m, xi = sp.symbols("H m xi", positive=True)
    a = -1 / (H * eta)
    nu2 = sp.Rational(9, 4) - m**2 / H**2 - 12 * xi
    expr = a**2 * m**2 + (6 * xi - 1) * sp.diff(a, eta, 2) / a
    return sp.simplify(expr + (nu2 - sp.Rational(1, 4)) / eta**2)


def nu_squared(m: float, xi: float, H: float = 1.0) -> float:
    return 9.0 / 4.0 - (m / H) ** 2 - 12.0 * xi


def nu(m: float, xi: float, H: float = 1.0) -> complex:
    """Principal square root of nu^2: real for nu^2 >= 0, purely imaginary otherwise."""
    n2 = nu_squared(m, xi, H)
    return complex(np.sqrt(n2)) if n2 >= 0 else 1j * np.sqrt(-n2)


def regime(nu2: float, tol: float = 1e-12) -> str:
    """Classify the parameter regime. The labels describe the mode behaviour, not 'allowed/forbidden'."""
    if abs(nu2 - 2.25) < tol:
        return "massless-minimal (nu=3/2): IR-marginal, |v_k|^2 ~ k^-3, no dS-invariant BD state"
    if nu2 > 2.25:
        return "nu>3/2 (tachyonic effective mass): k^2|phi_k|^2 not IR-integrable; BD two-point function IR divergent"
    if abs(nu2) < tol:
        return "nu=0: limiting case; Hankel H_0; complementary/principal boundary"
    if nu2 > 0:
        return "real nu in (0,3/2): complementary series, power-law late-time decay phi ~ eta^{3/2-nu}"
    return "imaginary nu: principal series (heavy), oscillatory decay phi ~ eta^{3/2} cos(mu ln|eta|+...)"


def omega_squared(eta, k, m=0.0, xi=0.0, H=1.0):
    """Time-dependent frequency  w^2 = k^2 + a^2 m^2 + (6 xi - 1) a''/a = k^2 - (nu^2 - 1/4)/eta^2."""
    eta = np.asarray(eta, dtype=float)
    return k**2 - (nu_squared(m, xi, H) - 0.25) / eta**2


@dataclass(frozen=True)
class ModeEquation:
    """v'' + [k^2 - (nu2 - 1/4)/eta^2] v = 0, parameterised by (k, m, xi, H)."""

    k: float
    m: float = 0.0
    xi: float = 0.0
    H: float = 1.0

    @property
    def nu2(self) -> float:
        return nu_squared(self.m, self.xi, self.H)

    @property
    def nu(self) -> complex:
        return nu(self.m, self.xi, self.H)

    def omega2(self, eta):
        return omega_squared(eta, self.k, self.m, self.xi, self.H)

    def residual(self, eta, v, dv, d2v):
        """v'' + w^2 v for externally supplied v, v', v''."""
        return d2v + self.omega2(eta) * v
