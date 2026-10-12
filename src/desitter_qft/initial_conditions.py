r"""Initial data for numerical mode evolution.

Two kinds:
 * `exact_bd_data`      : analytic Bunch-Davies value and derivative at eta0 (isolates solver error).
 * `adiabatic_data`     : WKB (adiabatic) vacuum of order 0, 2, 4 built *without* Hankel functions
                           (tests the physical prescription: the state is positive frequency at early times).
Adiabatic ansatz v = (2W)^{-1/2} exp(-i \int W),  W^2 = w^2 - [ W''/(2W) - 3 W'^2/(4 W^2) ].
We generate W by symbolic iteration of the right-hand side, starting from W = w.
The overall phase of v is a convention (it drops out of every observable); the numerical solution's
phase is fixed by this choice and is reported when compared with the analytic mode.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np
import sympy as sp

from .analytic_modes import bunch_davies_mode, bunch_davies_mode_derivative
from .mode_equations import ModeEquation


def exact_bd_data(eq: ModeEquation, eta0: float):
    return complex(bunch_davies_mode(eta0, eq.k, eq.nu)), complex(bunch_davies_mode_derivative(eta0, eq.k, eq.nu))


@lru_cache(maxsize=8)
def _wkb_expressions(order: int):
    """Return lambdified W(eta; k, c) and W'(eta; k, c) with w^2 = k^2 + c/eta^2, c = -(nu^2 - 1/4)."""
    eta, k, c = sp.symbols("eta k c", real=True)
    w2 = k**2 + c / eta**2
    W = sp.sqrt(w2)
    for _ in range(order // 2):
        W2 = w2 - (sp.diff(W, eta, 2) / (2 * W) - 3 * sp.diff(W, eta) ** 2 / (4 * W**2))
        W = sp.sqrt(W2)
    return sp.lambdify((eta, k, c), W, "numpy"), sp.lambdify((eta, k, c), sp.diff(W, eta), "numpy")


def adiabatic_data(eq: ModeEquation, eta0: float, order: int = 2):
    """(v, v') at eta0 from the adiabatic vacuum of given even order (0, 2, 4)."""
    if order not in (0, 2, 4):
        raise ValueError("order must be 0, 2 or 4")
    Wf, dWf = _wkb_expressions(order)
    c = -(eq.nu2 - 0.25)
    W = complex(Wf(eta0, eq.k, c))
    dW = complex(dWf(eta0, eq.k, c))
    v = 1.0 / np.sqrt(2 * W)
    dv = (-1j * W - dW / (2 * W)) * v
    return complex(v), complex(dv)
