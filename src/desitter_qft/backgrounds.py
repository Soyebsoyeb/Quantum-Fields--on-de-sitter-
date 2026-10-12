"""de Sitter background geometry (signature -+++, c = hbar = 1).

Conventions (see docs/conventions.md)
-------------------------------------
Christoffel:  Gamma^l_{mn} = 1/2 g^{ls}(d_m g_{sn} + d_n g_{sm} - d_s g_{mn})
Riemann:      R^r_{s m n} = d_m Gamma^r_{n s} - d_n Gamma^r_{m s}
                            + Gamma^r_{m l} Gamma^l_{n s} - Gamma^r_{n l} Gamma^l_{m s}
Ricci:        R_{s n} = R^r_{s r n},      R = g^{sn} R_{sn}
With these conventions de Sitter space has R = +12 H^2.

Units: H has dimension [mass]; t, eta, x have dimension [mass]^-1; a is dimensionless.
Only the spatially flat expanding chart (eta < 0) is implemented numerically.
The static patch metric is provided symbolically for the curvature check only.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import sympy as sp


@dataclass(frozen=True)
class DeSitter:
    """Spatially flat expanding patch of de Sitter, a(t) = exp(H t), eta in (-inf, 0)."""

    H: float = 1.0

    def __post_init__(self) -> None:
        if not (self.H > 0):
            raise ValueError("H must be positive")

    @property
    def ricci_scalar(self) -> float:
        return 12.0 * self.H**2

    @property
    def cosmological_constant(self) -> float:
        """Lambda = 3 H^2 so that G_{mn} + Lambda g_{mn} = 0 (vacuum)."""
        return 3.0 * self.H**2

    def scale_factor_cosmic(self, t):
        return np.exp(self.H * np.asarray(t, dtype=float))

    def scale_factor_conformal(self, eta):
        eta = np.asarray(eta, dtype=float)
        if np.any(eta >= 0):
            raise ValueError("conformal time must satisfy eta < 0 in the expanding patch")
        return -1.0 / (self.H * eta)

    def conformal_time(self, t):
        return -np.exp(-self.H * np.asarray(t, dtype=float)) / self.H

    def cosmic_time(self, eta):
        eta = np.asarray(eta, dtype=float)
        if np.any(eta >= 0):
            raise ValueError("eta must be negative")
        return -np.log(-self.H * eta) / self.H

    def efolds(self, eta, eta_ref):
        """N = ln[a(eta)/a(eta_ref)] = ln(eta_ref/eta)."""
        return np.log(np.asarray(eta_ref, dtype=float) / np.asarray(eta, dtype=float))

    def comoving_hubble_radius(self, eta):
        """1/(aH) = -eta."""
        return -np.asarray(eta, dtype=float)

    def horizon_exit_time(self, k):
        """Conformal time at which k = aH, i.e. k|eta| = 1."""
        return -1.0 / np.asarray(k, dtype=float)


# ----------------------------------------------------------------------
# Symbolic geometry
# ----------------------------------------------------------------------
def christoffel(g: sp.Matrix, X) -> list:
    n = len(X)
    ginv = g.inv()
    G = [[[0] * n for _ in range(n)] for _ in range(n)]
    for l in range(n):
        for m in range(n):
            for nn in range(n):
                G[l][m][nn] = sp.simplify(
                    sum(
                        ginv[l, s] * (sp.diff(g[s, nn], X[m]) + sp.diff(g[s, m], X[nn]) - sp.diff(g[m, nn], X[s]))
                        for s in range(n)
                    )
                    / 2
                )
    return G


def riemann(g: sp.Matrix, X, G=None):
    """Return R[r][s][m][n] = R^r_{s m n} in the documented convention."""
    n = len(X)
    G = G or christoffel(g, X)
    R = [[[[0] * n for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for r in range(n):
        for s in range(n):
            for m in range(n):
                for nn in range(n):
                    R[r][s][m][nn] = sp.simplify(
                        sp.diff(G[r][nn][s], X[m])
                        - sp.diff(G[r][m][s], X[nn])
                        + sum(G[r][m][l] * G[l][nn][s] - G[r][nn][l] * G[l][m][s] for l in range(n))
                    )
    return R


def ricci_tensor(g: sp.Matrix, X, G=None) -> sp.Matrix:
    n = len(X)
    Rm = riemann(g, X, G)
    return sp.Matrix(n, n, lambda s, nn: sp.simplify(sum(Rm[r][s][r][nn] for r in range(n))))


def ricci_scalar(g: sp.Matrix, X, G=None):
    Ric = ricci_tensor(g, X, G)
    ginv = g.inv()
    return sp.simplify(sum(ginv[i, j] * Ric[i, j] for i in range(len(X)) for j in range(len(X))))


def flat_chart_metric_cosmic(H=None):
    t, x, y, z = sp.symbols("t x y z", real=True)
    H = H if H is not None else sp.Symbol("H", positive=True)
    a = sp.exp(H * t)
    return sp.diag(-1, a**2, a**2, a**2), (t, x, y, z)


def flat_chart_metric_conformal(H=None):
    eta = sp.Symbol("eta", negative=True)
    x, y, z = sp.symbols("x y z", real=True)
    H = H if H is not None else sp.Symbol("H", positive=True)
    a = -1 / (H * eta)
    return sp.diag(-a**2, a**2, a**2, a**2), (eta, x, y, z)


def static_patch_metric(H=None):
    t, r, th, ph = sp.symbols("t r theta phi", positive=True)
    H = H if H is not None else sp.Symbol("H", positive=True)
    f = 1 - H**2 * r**2
    return sp.diag(-f, 1 / f, r**2, r**2 * sp.sin(th) ** 2), (t, r, th, ph)


def frw_conformal_metric():
    """Generic spatially flat FRW in conformal time with an unspecified a(eta)."""
    eta, x, y, z = sp.symbols("eta x y z", real=True)
    a = sp.Function("a", positive=True)(eta)
    return sp.diag(-a**2, a**2, a**2, a**2), (eta, x, y, z), a
