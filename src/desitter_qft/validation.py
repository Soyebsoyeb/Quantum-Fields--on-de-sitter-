"""Programmatic validation suite: every entry is an independent check with a stated tolerance.

Categories follow the repository taxonomy (symbolic verification, numerical experiment,
perturbative approximation, ...).  `run_all()` is exposed as `desitter-qft validate` and its
output is archived in results/validation.json.
"""
from __future__ import annotations

import json
from typing import Callable

import numpy as np
import sympy as sp

from . import backgrounds as bg
from . import mode_equations as me
from . import scalar_action as sa
from .analytic_modes import bunch_davies_mode, bunch_davies_mode_derivative, minimal_massless_mode
from .normalization import wronskian_residual

CHECKS: list = []


def check(name: str, category: str, tol: float):
    def deco(f: Callable[[], float]):
        CHECKS.append((name, category, tol, f))
        return f
    return deco


def _zero(expr) -> float:
    return 0.0 if sp.simplify(expr) == 0 else 1.0


@check("R = 12 H^2 (flat cosmic chart)", "symbolic verification", 0.0)
def _c1():
    g, X = bg.flat_chart_metric_cosmic()
    return _zero(bg.ricci_scalar(g, X) - 12 * sp.Symbol("H", positive=True) ** 2)


@check("R = 12 H^2 (conformal chart)", "symbolic verification", 0.0)
def _c2():
    g, X = bg.flat_chart_metric_conformal()
    return _zero(bg.ricci_scalar(g, X) - 12 * sp.Symbol("H", positive=True) ** 2)


@check("R = 12 H^2 (static patch)", "symbolic verification", 0.0)
def _c3():
    g, X = bg.static_patch_metric()
    return _zero(bg.ricci_scalar(g, X) - 12 * sp.Symbol("H", positive=True) ** 2)


@check("Euler-Lagrange = sqrt(-g)(Box - m^2 - xi R) phi", "symbolic verification", 0.0)
def _c4():
    return _zero(sa.euler_lagrange_vs_klein_gordon()[2])


@check("cosmic-time equation with 3H phi-dot term", "symbolic verification", 0.0)
def _c5():
    return _zero(sa.cosmic_time_equation_difference()[0])


@check("conformal-time equation vs chain rule from cosmic time", "symbolic verification", 0.0)
def _c6():
    return _zero(sa.cosmic_to_conformal_chain_rule())


@check("rescaled mode equation v'' + [k^2 + a^2 m^2 + (6 xi - 1) a''/a] v = 0", "symbolic verification", 0.0)
def _c7():
    return _zero(me.derive_rescaled_mode_equation()["difference"])


@check("de Sitter reduction to -(nu^2 - 1/4)/eta^2", "symbolic verification", 0.0)
def _c8():
    return _zero(me.de_sitter_reduction())


@check("Hankel mode Wronskian |W - i| (nu = 1.27, 0.8i, 3/2)", "analytic result, numerical check", 1e-10)
def _c9():
    eta = np.linspace(-40, -0.01, 400)
    return max(float(wronskian_residual(bunch_davies_mode(eta, 1.0, nu), bunch_davies_mode_derivative(eta, 1.0, nu)).max())
               for nu in (1.2689, 0.8j, 1.5))


@check("Hankel mode satisfies the ODE (mpmath derivatives)", "analytic result, numerical check", 1e-9)
def _c10():
    import mpmath as mp
    nu, k = 1.2689, 1.0
    f = lambda e: 0.5 * mp.sqrt(mp.pi) * mp.exp(1j * mp.pi * (nu / 2 + 0.25)) * mp.sqrt(-e) * mp.hankel1(nu, -k * e)
    return max(float(abs(mp.diff(f, e, 2) + (k**2 - (nu**2 - 0.25) / e**2) * f(e))) for e in (-3.0, -0.7, -0.1))


@check("closed-form nu = 3/2 mode equals Hankel mode", "analytic result, numerical check", 1e-12)
def _c11():
    eta = np.linspace(-20, -0.05, 200)
    v, dv = minimal_massless_mode(eta, 1.0)
    return float(max(np.abs(v - bunch_davies_mode(eta, 1.0, 1.5)).max(),
                     np.abs(dv - bunch_davies_mode_derivative(eta, 1.0, 1.5)).max()))


@check("numerical (DOP853) vs analytic mode, max abs error", "numerical experiment", 1e-7)
def _c12():
    from . import experiments as ex
    return max(r["max_abs_err"] for r in ex.free_modes_experiment(n=300).values())


@check("Wightman mode sum vs closed form (nu = 1, 1.3), max rel diff", "numerical experiment", 1e-9)
def _c13():
    from .two_point import wightman_closed_form_equal_time as cf, wightman_mode_sum as ms
    return max(abs(ms(-1.0, r, nu) - cf(-1.0, r, nu)) / cf(-1.0, r, nu) for nu in (1.0, 1.3) for r in (0.2, 1.0, 4.0))


@check("closed form satisfies the de Sitter-invariant ODE (relative residual)", "symbolic/numerical verification", 1e-20)
def _c14():
    from .two_point import invariant_ode_residual
    worst = 0.0
    for z in (-0.6, 0.2, 0.7):
        for nu in (0.3, 1.1, 0.5j):
            r, g = invariant_ode_residual(z, nu)
            worst = max(worst, r / g)
    return worst


@check("QM in-in O(g) slope vs exact diagonalisation (rel)", "perturbative approximation vs exact", 1e-4)
def _c15():
    from .perturbation_theory import qm_first_order_slope, qm_inin_prediction
    return abs(qm_first_order_slope(1.0) / qm_inin_prediction(1.0, 1.0) - 1)


@check("conformal contact correlator: numeric vs analytic (rel)", "perturbative approximation", 1e-9)
def _c16():
    from .in_in import conformal_contact_correlator as an, contact_correlator as nm
    ks = [1.0, 1.3, 0.7, 1.1]
    return max(abs(nm(e, ks, 0.5) / an(e, ks, 1.0) - 1) for e in (-2.0, -0.5, -0.05))


@check("tensor equation from linearised Einstein equations", "symbolic verification", 0.0)
def _c17():
    from .tensor_modes import derive_tensor_equation
    r = derive_tensor_equation()
    return _zero(r["difference_from_claim"]) + _zero(r["background_residual"])


def run_all():
    out = []
    for name, cat, tol, f in CHECKS:
        try:
            val = float(f())
            ok = val <= tol
        except Exception as e:  # report failures instead of hiding them
            val, ok = float("nan"), False
            name = f"{name} [exception: {type(e).__name__}: {e}]"
        out.append(dict(name=name, category=cat, value=val, tolerance=tol, passed=bool(ok)))
    return out


def write_validation_record(report, path):
    from .provenance import environment_info
    with open(path, "w") as fh:
        json.dump({"checks": report, "environment": environment_info()}, fh, indent=2)
