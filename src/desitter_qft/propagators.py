r"""Mode-space (fixed comoving k) two-point functions and Green functions for chi = a phi.

With v_k a normalised mode (v v*' - v* v' = i) define, for the operator D = d_eta^2 + w_k^2(eta):
    G^+(eta,eta') = v(eta) v*(eta')                      Wightman,   D G^+ = 0
    G^F(eta,eta') = th(eta-eta') v(eta)v*(eta') + th(eta'-eta) v*(eta)v(eta')     D G^F = -i delta
    G^R(eta,eta') = i th(eta-eta') [v(eta)v*(eta') - v*(eta)v(eta')] = -2 th Im[v(eta) v*(eta')]  D G^R = +delta
    G^A(eta,eta') = G^R(eta',eta)                        D G^A = +delta
The delta-function strengths follow from the Wronskian: the derivative jump at eta = eta' is
d_eta G^R|_{eta'+} = i (v' v* - v*' v) = i(-i) = 1.
These are *mode-space* objects (category 1/3); position-space propagators require the k-integral and a
distributional prescription (see two_point.py and docs/propagators.md).
"""
from __future__ import annotations

import numpy as np


def wightman_k(v, eta, etap):
    return v(eta) * np.conj(v(etap))


def feynman_k(v, eta, etap):
    eta, etap = np.broadcast_arrays(np.asarray(eta, float), np.asarray(etap, float))
    return np.where(eta >= etap, v(eta) * np.conj(v(etap)), np.conj(v(eta)) * v(etap))


def retarded_k(v, eta, etap):
    eta, etap = np.broadcast_arrays(np.asarray(eta, float), np.asarray(etap, float))
    return np.where(eta > etap, -2.0 * np.imag(v(eta) * np.conj(v(etap))), 0.0)


def advanced_k(v, eta, etap):
    return retarded_k(v, etap, eta)


def retarded_from_ode(eq, etap, eta_grid, rtol=1e-12, atol=1e-14):
    """Independent construction: solve D y = 0 for eta>eta' with y(eta')=0, y'(eta')=1.
    Uses only the real ODE, no Hankel functions."""
    from .numerical_solver import solve_mode

    eta_grid = np.asarray(eta_grid, float)
    sol = solve_mode(eq, etap, eta_grid[-1], 0j, 1.0 + 0j, rtol=rtol, atol=atol, eta_eval=eta_grid)
    return sol.v.real
