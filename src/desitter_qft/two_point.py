r"""Equal-time Wightman function of the Bunch-Davies state, from modes and in closed form.

Mode-sum (category 4/3):
    G^+(eta, r) = <phi(eta,x) phi(eta,x')> = \int d^3k/(2pi)^3 |phi_k|^2 e^{ik.(x-x')}
                = (1/2pi^2) \int_0^inf dk k^2 |phi_k|^2 sin(kr)/(kr),      phi_k = v_k/a.
The integral is not absolutely convergent (|phi_k|^2 ~ H^2 eta^2/(2k)): it defines a distribution.
For r>0 we use the *Abel* prescription, implemented by subtracting the conformal-vacuum piece
1/(2k), whose Abel sum is known exactly,
    (1/2pi^2) \int dk k^2 (1/2k) sin(kr)/(kr) = 1/(4 pi^2 r^2)   (Abel/distributional),
    G^+ = H^2 eta^2 [ 1/(4 pi^2 r^2) + (1/2pi^2) \int dk k (|v_k|^2 - 1/2k) sin(kr)/r ].
The remaining integrand decays like k^-2 and is handled by QUADPACK (QAWF for the oscillatory tail).

Closed form (Chernikov-Tagirov / Bunch-Davies), written in the flat chart with
    Z = (eta^2 + eta'^2 - r^2)/(2 eta eta'),
    G^+ = H^2/(16 pi^2) Gamma(3/2+nu) Gamma(3/2-nu) 2F1(3/2+nu, 3/2-nu; 2; (1+Z)/2).
It is checked here independently: (i) by substituting into the de Sitter-invariant ODE
    (1 - Z^2) G'' - 4 Z G' - (m^2/H^2 + 12 xi) G = 0,
(ii) against the mode sum, (iii) for nu=1/2 against H^2 eta eta'/(4 pi^2 (r^2 - (eta-eta')^2)).
"""
from __future__ import annotations

import mpmath as mp
import numpy as np
from scipy.integrate import quad
from scipy.special import jv, yv

from .analytic_modes import bunch_davies_mode


def mode_power_v(eta: float, k: float, nu) -> float:
    """|v_k(eta)|^2 for the BD mode."""
    nu = complex(nu)
    if abs(nu.imag) < 1e-15:
        x = -k * eta
        n = nu.real
        out = np.pi / 4 * (-eta) * (jv(n, x) ** 2 + yv(n, x) ** 2)
    else:
        out = np.abs(bunch_davies_mode(eta, k, nu)) ** 2
    return float(out) if np.ndim(out) == 0 else out


def mode_power_phi(eta, k, nu, H=1.0):
    return H**2 * eta**2 * mode_power_v(eta, k, nu)


def invariant_Z(eta, etap, r):
    return (eta**2 + etap**2 - r**2) / (2.0 * eta * etap)


def wightman_closed_form(Z, nu, H=1.0):
    """Closed-form BD Wightman function as function of the de Sitter invariant Z (spacelike: Z<1)."""
    nu = mp.mpmathify(complex(nu))
    a, b = mp.mpf(3) / 2 + nu, mp.mpf(3) / 2 - nu
    val = H**2 / (16 * mp.pi**2) * mp.gamma(a) * mp.gamma(b) * mp.hyp2f1(a, b, 2, (1 + mp.mpf(Z)) / 2)
    return float(val.real)


def wightman_closed_form_equal_time(eta, r, nu, H=1.0):
    return wightman_closed_form(invariant_Z(eta, eta, r), nu, H)


def invariant_ode_residual(Z, nu, H=1.0, dps=30):
    """Residual of (1-Z^2)G'' - 4ZG' - (9/4 - nu^2) G = 0 using mpmath numerical derivatives."""
    with mp.workdps(dps):
        nu_ = mp.mpmathify(complex(nu))
        M2 = mp.mpf(9) / 4 - nu_**2
        f = lambda z: H**2 / (16 * mp.pi**2) * mp.gamma(mp.mpf(3) / 2 + nu_) * mp.gamma(mp.mpf(3) / 2 - nu_) * mp.hyp2f1(
            mp.mpf(3) / 2 + nu_, mp.mpf(3) / 2 - nu_, 2, (1 + z) / 2)
        z = mp.mpf(Z)
        G, G1, G2 = f(z), mp.diff(f, z), mp.diff(f, z, 2)
        res = (1 - z**2) * G2 - 4 * z * G1 - M2 * G
        return float(abs(res)), float(abs(G))


def wightman_mode_sum(eta: float, r: float, nu, H=1.0, epsabs=1e-12, k_split=None):
    """Equal-time G^+(eta, r) from the Abel-regularised mode sum (see module docstring)."""
    k0 = k_split if k_split is not None else 10.0 / r
    if abs(complex(nu) - 0.5) < 1e-14:
        # |v_k|^2 = 1/(2k) exactly: the subtracted integrand vanishes identically (conformal vacuum).
        return H**2 * eta**2 / (4 * np.pi**2 * r**2)

    def delta_k(k):
        return k * (mode_power_v(eta, k, nu) - 1.0 / (2.0 * k)) / r

    head, _ = quad(lambda k: delta_k(k) * np.sin(k * r), 0.0, k0, limit=400, epsabs=epsabs)
    tail, _ = quad(delta_k, k0, np.inf, weight="sin", wvar=r, limlst=200, epsabs=epsabs)
    return H**2 * eta**2 * (1.0 / (4 * np.pi**2 * r**2) + (head + tail) / (2 * np.pi**2))


def short_distance_coefficients(nu, H=1.0):
    """Coefficients of G = A/u + B ln u + C + O(u ln u), with u = (1-Z)/2 -> 0 (DLMF 15.8.10, m = 1):
        A = H^2/(16 pi^2),  B = (9/4 - nu^2 - 2) A,  C = (9/4 - nu^2 - 2) A [psi(3/2+nu) + psi(3/2-nu) + 2 gamma_E - 1].
    The numerical values are verified against the closed form in tests/test_two_point.py."""
    nu = mp.mpmathify(complex(nu))
    A = H**2 / (16 * mp.pi**2)
    M2 = mp.mpf(9) / 4 - nu**2
    B = (M2 - 2) * A
    C = (M2 - 2) * A * (mp.digamma(mp.mpf(3) / 2 + nu) + mp.digamma(mp.mpf(3) / 2 - nu) + 2 * mp.euler - 1)
    return float(A), float(mp.re(B)), float(mp.re(C))


def coincident_phi2_pointsplit(nu, H=1.0, mu=None):
    """<phi^2> by point splitting with the *explicitly stated* subtraction
        G - 1/(4 pi^2 sigma^2) - B ln(mu^2 sigma^2)  as sigma -> 0,
    sigma = geodesic distance, u = sin^2(H sigma/2).  Result: A/3 + C + B ln(H^2/(4 mu^2)).
    This is a *scheme-dependent* finite part (mu, and the choice of subtracting 1/sigma^2 rather than
    the full Hadamard parametrix, shift it). It is NOT presented as a unique renormalised value."""
    A, B, C = short_distance_coefficients(nu, H)
    mu = H if mu is None else mu
    return A / 3 + C + B * np.log(H**2 / (4 * mu**2))


def phi2_cutoff(eta, k_uv, k_ir, nu, H=1.0, n_seg=None):
    r"""\int_{k_ir}^{k_uv} dk k^2 |phi_k|^2/(2 pi^2): a sharp comoving cutoff (NOT a renormalisation)."""
    f = lambda k: k**2 * mode_power_phi(eta, k, nu, H) / (2 * np.pi**2)
    n_seg = n_seg or max(60, int(3 * np.log(k_uv / k_ir)))
    pts = np.geomspace(k_ir, k_uv, n_seg)
    total = 0.0
    for a, b in zip(pts[:-1], pts[1:]):
        total += quad(f, a, b, epsabs=0, epsrel=1e-11)[0]
    return total


def phi2_cutoff_minimal_massless_analytic(eta, k_uv, k_ir, H=1.0):
    """Closed form for nu = 3/2:  (H^2/4pi^2) [ln(k_uv/k_ir) + eta^2 (k_uv^2 - k_ir^2)/2]."""
    return H**2 / (4 * np.pi**2) * (np.log(k_uv / k_ir) + eta**2 * (k_uv**2 - k_ir**2) / 2)


def wightman_closed_form_u(u, nu, H=1.0, dps=40):
    """Closed form at small u = (1 - Z)/2 evaluated in extended precision (avoids float cancellation in Z)."""
    with mp.workdps(dps):
        nu_ = mp.mpmathify(complex(nu))
        a, b = mp.mpf(3) / 2 + nu_, mp.mpf(3) / 2 - nu_
        return mp.re(H**2 / (16 * mp.pi**2) * mp.gamma(a) * mp.gamma(b) * mp.hyp2f1(a, b, 2, 1 - mp.mpf(u)))
