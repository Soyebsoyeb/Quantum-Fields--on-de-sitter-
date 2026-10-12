r"""Tensor perturbations on a fixed de Sitter background (linear order, TT gauge).

Metric  ds^2 = a^2 [ -d eta^2 + (delta_ij + h_ij) dx^i dx^j ],  h_ii = 0 = d_i h_ij.
* `derive_tensor_equation` (symbolic verification): for the plus polarisation h_xx = -h_yy = eps h(eta, z) the
  first-order part of R_xx - Lambda g_xx (Lambda = 3H^2) equals (1/2)[h'' + 2(a'/a) h' - d_z^2 h].
* Quadratic action S2 = (M_p^2/8) \int d eta d^3x a^2 [ (h_ij')^2 - (d_k h_ij)^2 ],  M_p^2 = 1/(8 pi G)
  (standard result, e.g. Maldacena 2003, Mukhanov-Feldman-Brandenberger 1992; NOT re-derived symbolically here:
  the equation of motion above IS derived, the overall normalisation of S2 is quoted).
* Canonical variable per polarisation: h_ij = sum_lambda e^lambda_ij h_lambda, e^l_ij e^l'_ij = 2 delta^{ll'};
  v_lambda = (M_p a / sqrt 2) h_lambda is then a massless, minimally coupled scalar (nu = 3/2) in the sense of
  mode_equations.py, so P_h(k) per polarisation = (2/M_p^2) |phi_k|^2 with phi_k = v/a.
* Late-time dimensionless spectrum  Delta_t^2 = k^3/(2 pi^2) <h_ij h_ij>' = 2 H^2/(pi^2 M_p^2).
This is a fixed-background *linear* perturbation theory. Quantising h_ij as a field on dS is the free
graviton sector; constraints, gauge fixing and interactions are discussed only conceptually (docs/tensor_perturbations.md).
"""
from __future__ import annotations

import numpy as np
import sympy as sp

from .analytic_modes import bunch_davies_mode


def derive_tensor_equation():
    from .backgrounds import ricci_tensor

    eta = sp.Symbol("eta", negative=True)
    x, y, z = sp.symbols("x y z", real=True)
    H = sp.Symbol("H", positive=True)
    eps = sp.Symbol("epsilon")
    a = -1 / (H * eta)
    h = sp.Function("h")(eta, z)
    g = sp.diag(-(a**2), a**2 * (1 + eps * h), a**2 * (1 - eps * h), a**2)
    Ric = ricci_tensor(g, (eta, x, y, z))
    E = Ric - 3 * H**2 * g
    E0 = sp.simplify(E[1, 1].subs(eps, 0))
    E1 = sp.simplify(sp.diff(E[1, 1], eps).subs(eps, 0))
    target = sp.Rational(1, 2) * (sp.diff(h, eta, 2) + 2 * sp.diff(a, eta) / a * sp.diff(h, eta) - sp.diff(h, z, 2))
    return {"background_residual": E0, "first_order": E1, "difference_from_claim": sp.simplify(E1 - target)}


def tensor_power_per_polarisation(eta, k, H=1.0, Mp=1.0):
    """P_h(k, eta) = (2/Mp^2) H^2 eta^2 |v_k|^2 with the nu = 3/2 Bunch-Davies mode."""
    return 2.0 / Mp**2 * H**2 * np.asarray(eta) ** 2 * np.abs(bunch_davies_mode(eta, k, 1.5)) ** 2


def dimensionless_tensor_spectrum(eta, k, H=1.0, Mp=1.0):
    """Delta_t^2(k, eta) = k^3/(2 pi^2) * (sum over 2 polarisations, e_ij e_ij = 2 each) P_h = k^3/(2pi^2) 4 P_h/ 1... see below.
    <h_ij h_ij>' = sum_lambda 2 P_h = 4 P_h  =>  Delta_t^2 = k^3 * 4 P_h / (2 pi^2)."""
    return k**3 * 4.0 * tensor_power_per_polarisation(eta, k, H, Mp) / (2 * np.pi**2)


def tensor_spectrum_late_time(H=1.0, Mp=1.0):
    return 2.0 * H**2 / (np.pi**2 * Mp**2)
