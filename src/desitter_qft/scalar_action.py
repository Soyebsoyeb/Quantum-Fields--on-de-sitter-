r"""Symbolic derivation of the Klein-Gordon equation from the action.

S = -1/2 \int d^4x sqrt(-g) [ g^{mn} d_m phi d_n phi + (m^2 + xi R) phi^2 ]

Everything here is a *symbolic verification* (category 2): sympy's Euler-Lagrange
operator is compared with the covariant expression sqrt(-g)(Box - m^2 - xi R) phi.
Boundary terms from integration by parts are assumed to vanish (fields of compact
support or sufficiently fast falloff in the spatial directions and fixed at the
initial/final time slices).
"""
from __future__ import annotations

import sympy as sp
from sympy.calculus.euler import euler_equations

from .backgrounds import frw_conformal_metric, flat_chart_metric_cosmic, ricci_scalar


def covariant_box(g: sp.Matrix, X, phi, sqrt_mg):
    """Box phi = (1/sqrt(-g)) d_m ( sqrt(-g) g^{mn} d_n phi )."""
    ginv = g.inv()
    n = len(X)
    return sp.simplify(
        sum(sp.diff(sqrt_mg * ginv[m, nn] * sp.diff(phi, X[nn]), X[m]) for m in range(n) for nn in range(n)) / sqrt_mg
    )


def lagrangian_density(g, X, phi, sqrt_mg, m, xi, R):
    ginv = g.inv()
    n = len(X)
    kin = sum(ginv[i, j] * sp.diff(phi, X[i]) * sp.diff(phi, X[j]) for i in range(n) for j in range(n))
    return -sp.Rational(1, 2) * sqrt_mg * (kin + (m**2 + xi * R) * phi**2)


def euler_lagrange_vs_klein_gordon(metric=None):
    """Return (EL, KG, difference) for a metric (default: generic conformal FRW).

    EL  : sympy Euler-Lagrange expression dL/dphi - d_m dL/d(d_m phi)
    KG  : sqrt(-g) (Box - m^2 - xi R) phi
    The difference must simplify to zero.
    """
    if metric is None:
        g, X, a = frw_conformal_metric()
        sqrt_mg = a**4  # a > 0 assumed; sqrt(-det g) = a^4
    else:
        g, X, sqrt_mg = metric
    m, xi = sp.symbols("m xi", real=True)
    phi = sp.Function("phi")(*X)
    R = ricci_scalar(g, X)
    L = lagrangian_density(g, X, phi, sqrt_mg, m, xi, R)
    EL = euler_equations(L, phi, list(X))[0].lhs
    KG = sqrt_mg * (covariant_box(g, X, phi, sqrt_mg) - (m**2 + xi * R) * phi)
    return sp.simplify(EL), sp.simplify(KG), sp.simplify(EL - KG)


def conformal_time_equation():
    """phi'' + 2 (a'/a) phi' - nabla^2 phi + a^2 (m^2 + xi R) phi = 0 with R = 6 a''/a^3.

    Returns (equation_lhs, R, phi, a, symbols).  Derived from covariant Box, not recalled.
    """
    g, X, a = frw_conformal_metric()
    eta = X[0]
    m, xi = sp.symbols("m xi", real=True)
    phi = sp.Function("phi")(*X)
    sqrt_mg = a**4
    R = ricci_scalar(g, X)
    box = covariant_box(g, X, phi, sqrt_mg)
    kg = -a**2 * (box - (m**2 + xi * R) * phi)  # multiply by -a^2 so phi'' has unit coefficient
    lap = sum(sp.diff(phi, c, 2) for c in X[1:])
    target = sp.diff(phi, eta, 2) + 2 * sp.diff(a, eta) / a * sp.diff(phi, eta) - lap + a**2 * (m**2 + xi * R) * phi
    return sp.simplify(kg - target), R, phi, a, (m, xi)


def cosmic_time_equation_difference():
    """Check  phi_tt + 3H phi_t - nabla^2 phi / a^2 + (m^2 + xi R) phi  ==  -(Box - m^2 - xi R) phi
    in the cosmic-time chart; returns sympy difference (should be 0).
    """
    g, X = flat_chart_metric_cosmic()
    t = X[0]
    H = sp.Symbol("H", positive=True)
    a = sp.exp(H * t)
    m, xi = sp.symbols("m xi", real=True)
    phi = sp.Function("phi")(*X)
    sqrt_mg = a**3
    R = ricci_scalar(g, X)
    box = covariant_box(g, X, phi, sqrt_mg)
    lap = sum(sp.diff(phi, c, 2) for c in X[1:])
    kg = sp.diff(phi, t, 2) + 3 * H * sp.diff(phi, t) - lap / a**2 + (m**2 + xi * R) * phi
    return sp.simplify(kg + (box - (m**2 + xi * R) * phi)), sp.simplify(R)


def cosmic_to_conformal_chain_rule():
    """Independent derivation of the conformal-time equation by the chain rule.

    Take phi(t) = Phi(eta(t)) with eta = -exp(-H t)/H (spatially homogeneous part; the
    gradient term maps trivially).  Then
        phi_tt + 3H phi_t + (m^2 + xi R) phi  ==  a^{-2} [Phi'' + 2(a'/a) Phi' + a^2 (m^2 + xi R) Phi]
    Returns the simplified difference (should be 0).
    """
    t = sp.Symbol("t", real=True)
    H = sp.Symbol("H", positive=True)
    m, xi = sp.symbols("m xi", real=True)
    eta = sp.Symbol("eta", negative=True)
    Phi = sp.Function("Phi")
    a_eta = -1 / (H * eta)
    eta_t = -sp.exp(-H * t) / H
    phi_t = Phi(eta_t)
    R = 12 * H**2
    lhs = sp.diff(phi_t, t, 2) + 3 * H * sp.diff(phi_t, t) + (m**2 + xi * R) * phi_t
    # conformal operator evaluated on Phi(eta), then substituted eta -> eta(t)
    conf = (sp.diff(Phi(eta), eta, 2) + 2 * sp.diff(a_eta, eta) / a_eta * sp.diff(Phi(eta), eta) + a_eta**2 * (m**2 + xi * R) * Phi(eta)) / a_eta**2
    conf_t = conf.subs(eta, eta_t).doit()
    return sp.simplify((lhs - conf_t).doit())
