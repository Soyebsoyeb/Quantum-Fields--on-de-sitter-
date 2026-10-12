"""Analytic Bunch-Davies mode functions of v'' + [k^2 - (nu^2 - 1/4)/eta^2] v = 0.

    v_k(eta) = (sqrt(pi)/2) exp[i pi (nu/2 + 1/4)] sqrt(-eta) H^{(1)}_nu(-k eta)

Normalised so that  v v*' - v* v' = i  (see normalization.py) and
v -> exp(-i k eta)/sqrt(2k) as k|eta| -> infinity (positive frequency, DLMF 10.17.5).
Real nu uses scipy; imaginary nu (principal series) uses mpmath (complex order).
The physical field mode is phi_k = v_k / a = -H eta v_k.
"""
from __future__ import annotations

import numpy as np
import mpmath as mp
from scipy.special import hankel1


def _is_real_order(nu) -> bool:
    return abs(complex(nu).imag) < 1e-15


def _hankel1(order, x):
    x = np.asarray(x)
    if _is_real_order(order):
        return hankel1(complex(order).real, x)
    f = np.frompyfunc(lambda z: complex(mp.hankel1(complex(order), complex(z))), 1, 1)
    return np.array(f(x), dtype=complex)


def bunch_davies_mode(eta, k: float, nu):
    """v_k(eta) for eta<0."""
    eta = np.asarray(eta)
    x = -k * eta
    nu = complex(nu)
    pref = 0.5 * np.sqrt(np.pi) * np.exp(1j * np.pi * (nu / 2 + 0.25))
    return pref * np.sqrt(-eta) * _hankel1(nu, x)


def bunch_davies_mode_derivative(eta, k: float, nu):
    """dv/deta using H'_nu(x) = H_{nu-1}(x) - (nu/x) H_nu(x)  (DLMF 10.6.2)."""
    eta = np.asarray(eta)
    x = -k * eta
    nu = complex(nu)
    pref = 0.5 * np.sqrt(np.pi) * np.exp(1j * np.pi * (nu / 2 + 0.25))
    H = _hankel1(nu, x)
    Hm = _hankel1(nu - 1, x)
    dH = Hm - nu / x * H
    return pref * (-0.5 / np.sqrt(-eta) * H - k * np.sqrt(-eta) * dH)


# -- closed-form benchmarks ---------------------------------------------
def conformal_mode(eta, k):
    """m = 0, xi = 1/6  (nu = 1/2):  v = exp(-i k eta)/sqrt(2k)."""
    eta = np.asarray(eta, dtype=float)
    v = np.exp(-1j * k * eta) / np.sqrt(2 * k)
    return v, -1j * k * v


def minimal_massless_mode(eta, k):
    """m = 0, xi = 0  (nu = 3/2):  v = exp(-i k eta)/sqrt(2k) (1 - i/(k eta))."""
    eta = np.asarray(eta, dtype=float)
    e = np.exp(-1j * k * eta) / np.sqrt(2 * k)
    v = e * (1 - 1j / (k * eta))
    dv = e * (-1j * k * (1 - 1j / (k * eta)) + 1j / (k * eta**2))
    return v, dv


def field_mode(eta, k, nu, H=1.0):
    """phi_k = v_k/a = -H eta v_k."""
    return -H * np.asarray(eta) * bunch_davies_mode(eta, k, nu)


def late_time_amplitude_squared(k: float, nu: float, H: float = 1.0, eta: float = -1e-8):
    """Leading small-(-k eta) behaviour |phi_k|^2 ~ H^2 eta^2 (pi/4)(-eta)|H_nu|^2 using
    |H_nu(x)| ~ Gamma(nu)/pi (2/x)^nu for real nu>0 (DLMF 10.7.2)."""
    from math import gamma, pi

    x = -k * eta
    return H**2 * eta**2 * (pi / 4) * (-eta) * (gamma(nu) / pi) ** 2 * (2 / x) ** (2 * nu)
