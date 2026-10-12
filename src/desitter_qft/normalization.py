r"""Canonical normalisation (Wronskian) checks.

Derivation (docs/canonical_quantization.md):  with chi = a phi, pi_chi = chi',
    chi(eta,x) = \int d^3k/(2pi)^3 [a_k v_k(eta) e^{ik.x} + h.c.],   [a_k, a_k'^dag] = (2pi)^3 delta^3(k-k')
and [chi, pi_chi] = i delta^3 imply
    W[v] := v v*' - v* v' = i     (independent of eta for the real ODE).
The checker NEVER rescales a solution; it only reports residuals |W - i|.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np


def wronskian(v, dv):
    v = np.asarray(v)
    dv = np.asarray(dv)
    return v * np.conj(dv) - np.conj(v) * dv


def wronskian_residual(v, dv):
    """|W - i| pointwise."""
    return np.abs(wronskian(v, dv) - 1j)


@dataclass
class NormalizationReport:
    initial_residual: float
    max_residual: float
    final_residual: float
    max_residual_location: float
    rtol: float | None
    atol: float | None
    method: str | None
    growth_factor: float  # max_residual / max(initial_residual, eps)
    passes: bool
    tolerance_used: float
    note: str

    def as_dict(self):
        return asdict(self)


def check_normalization(eta, v, dv, rtol=None, atol=None, method=None, tolerance=None) -> NormalizationReport:
    """Assess whether the Wronskian condition holds along a trajectory.

    `tolerance` defaults to max(1e3*rtol, 1e-12) if rtol is given, else 1e-8.  The factor
    1e3 absorbs the accumulation of local error over the integration interval; the report
    includes the measured growth so a user can judge that choice.
    """
    res = wronskian_residual(v, dv)
    i = int(np.argmax(res))
    tol = tolerance if tolerance is not None else (max(1e3 * rtol, 1e-12) if rtol else 1e-8)
    r0 = float(res[0])
    return NormalizationReport(
        initial_residual=r0,
        max_residual=float(res[i]),
        final_residual=float(res[-1]),
        max_residual_location=float(np.asarray(eta)[i]),
        rtol=rtol,
        atol=atol,
        method=method,
        growth_factor=float(res[i] / max(r0, np.finfo(float).eps)),
        passes=bool(res[i] <= tol),
        tolerance_used=float(tol),
        note="no rescaling applied; residual is |v v*' - v* v' - i|",
    )
