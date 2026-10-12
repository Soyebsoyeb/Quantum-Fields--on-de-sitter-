r"""Infrared structure: coarse-grained variance and the stochastic (Starobinsky-Yokoyama) equilibrium.

* `ir_variance`: <phi^2>_{k_ir < k < 1/|eta|} from the exact modes (numerical experiment).
* `stochastic_equilibrium_moments`: P_eq(phi) ~ exp[-8 pi^2 V(phi)/(3 H^4)], the stationary solution of the
  Fokker-Planck equation of Starobinsky (1986) / Starobinsky-Yokoyama (1994).  This is an *effective
  long-wavelength description*; it is used here as an independent IR benchmark, not as a derivation from the
  full quantum theory.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import quad

from .two_point import phi2_cutoff


def ir_variance(eta, k_ir, nu, H=1.0, x_uv=1.0):
    """Variance of modes with k_ir < k < x_uv/|eta| (super-horizon coarse graining)."""
    return phi2_cutoff(eta, x_uv / abs(eta), k_ir, nu, H)


def stochastic_equilibrium_moments(V, H=1.0, phi_max=None, n_moments=(2,)):
    """<phi^n> in P_eq ~ exp(-8 pi^2 V/(3 H^4)); V is a callable."""
    c = 8 * np.pi**2 / (3 * H**4)
    L = phi_max if phi_max is not None else 250 * H
    w = lambda p: np.exp(-c * V(p))
    Z = quad(w, -L, L, limit=500, points=[0.0])[0]
    return {n: quad(lambda p: p**n * w(p), -L, L, limit=500, points=[0.0])[0] / Z for n in n_moments}


def free_massive_stochastic_phi2(m, H=1.0):
    """Closed form  3 H^4/(8 pi^2 m^2)  for V = m^2 phi^2/2 (Gaussian P_eq)."""
    return 3 * H**4 / (8 * np.pi**2 * m**2)
