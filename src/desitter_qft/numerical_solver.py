"""Adaptive numerical integration of the rescaled mode equation (category 3: numerical experiment).

State vector y = (Re v, Im v, Re v', Im v').  Default method DOP853 (Hairer-Norsett-Wanner).
Nothing is rescaled after integration; normalisation is audited in normalization.py.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp

from .mode_equations import ModeEquation
from .provenance import environment_info


@dataclass
class ModeSolution:
    eta: np.ndarray
    v: np.ndarray
    dv: np.ndarray
    equation: ModeEquation
    solver: dict
    stats: dict
    environment: dict = field(default_factory=dict)


def solve_mode(
    eq: ModeEquation,
    eta0: float,
    eta1: float,
    v0: complex,
    dv0: complex,
    n_out: int = 400,
    method: str = "DOP853",
    rtol: float = 1e-10,
    atol: float = 1e-12,
    eta_eval=None,
) -> ModeSolution:
    if not (eta0 < eta1 < 0):
        raise ValueError("require eta0 < eta1 < 0")
    nu2m = eq.nu2 - 0.25

    def rhs(eta, y):
        w2 = eq.k**2 - nu2m / eta**2
        return [y[2], y[3], -w2 * y[0], -w2 * y[1]]

    t_eval = np.linspace(eta0, eta1, n_out) if eta_eval is None else np.asarray(eta_eval)
    sol = solve_ivp(rhs, (eta0, eta1), [v0.real, v0.imag, dv0.real, dv0.imag], method=method,
                    rtol=rtol, atol=atol, t_eval=t_eval)
    if not sol.success:
        raise RuntimeError(f"integration failed: {sol.message}")
    y = sol.y
    return ModeSolution(
        eta=sol.t,
        v=y[0] + 1j * y[1],
        dv=y[2] + 1j * y[3],
        equation=eq,
        solver={"method": method, "rtol": rtol, "atol": atol, "eta0": eta0, "eta1": eta1, "precision": "float64"},
        stats={"nfev": int(sol.nfev), "njev": int(sol.njev), "nlu": int(sol.nlu)},
        environment=environment_info(),
    )
