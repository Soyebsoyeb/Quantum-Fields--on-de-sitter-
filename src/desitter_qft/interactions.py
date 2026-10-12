r"""Scalar self-interaction  S_int = -\int d^4x sqrt(-g) (lambda/4!) phi^4  (category 1/2).

* `interaction_hamiltonian_density` derives H_int by a Legendre transform (symbolic verification).
* `check_dimensions` verifies that every term of the action carries mass dimension 4 and that [lambda] = 0.
No interaction-picture *state* is assumed here; see in_in.py for the contour prescription.
"""
from __future__ import annotations

import sympy as sp


def interaction_hamiltonian_density():
    """Legendre transform in conformal time for L = a^2/2 (phi'^2 - (grad phi)^2) - a^4 [m^2/2 + xi R/2] phi^2 - a^4 lam/4! phi^4.

    Returns (H_full - H_free) expressed in (phi, pi); expected  a^4 lam phi^4/4!.  Because the interaction
    contains no derivatives, pi = a^2 phi' is unchanged and H_int = -L_int."""
    a, lam, m, xi, R = sp.symbols("a lambda m xi R", positive=True)
    phi, dphi, grad2, pi = sp.symbols("phi dphi grad2 pi", real=True)

    def lagr(l):
        return a**2 / 2 * (dphi**2 - grad2) - a**4 * (m**2 + xi * R) * phi**2 / 2 - a**4 * l / sp.factorial(4) * phi**4

    def ham(l):
        L = lagr(l)
        p = sp.diff(L, dphi)
        sol = sp.solve(sp.Eq(pi, p), dphi)[0]
        return sp.simplify((p * dphi - L).subs(dphi, sol))

    H_int = sp.simplify(ham(lam) - ham(0))
    return H_int, a**4 * lam * phi**4 / 24


def check_dimensions():
    """Mass-dimension bookkeeping with c = hbar = 1.  Returns dict of term -> dimension and a pass flag."""
    d = {"phi": 1, "m": 1, "H": 1, "R": 2, "lambda": 0, "xi": 0, "d4x": -4, "sqrt_g": 0, "d_mu": 1}
    terms = {
        "kinetic  sqrt(-g) (d phi)^2": d["d4x"] + d["sqrt_g"] + 2 * d["d_mu"] + 2 * d["phi"],
        "mass     sqrt(-g) m^2 phi^2": d["d4x"] + 2 * d["m"] + 2 * d["phi"],
        "curv.    sqrt(-g) xi R phi^2": d["d4x"] + d["xi"] + d["R"] + 2 * d["phi"],
        "interact sqrt(-g) lam phi^4": d["d4x"] + d["lambda"] + 4 * d["phi"],
    }
    return terms, all(v == 0 for v in terms.values())
