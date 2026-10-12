r"""In-in (Schwinger-Keldysh) first-order perturbation theory for contact diagrams.

Setup (conformal time, interaction picture).  H_I(eta) = \int d^3x a^4 (lambda/4!) phi_I^4 = \int d^3x (lambda/4!) chi^4,
chi = a phi.  Initial state: the free Bunch-Davies vacuum prescribed in the infinite past via the i-epsilon
contour (eta' -> -infinity (1 + i eps) on the branch carrying e^{-i K eta'}).  Then for Q(eta) = prod_i phi_{k_i}(eta)
    <Q> = < U^dag Q U >,  U = T exp(-i \int_{-inf}^{eta} H_I)
    <Q>^{(1)} = i \int d eta' < [H_I(eta'), Q(eta)] >  = -2 Im X,   X = \int d eta' < H_I(eta') Q(eta) >.
Wick contraction (4! ways, all four legs attached) gives the connected, momentum-stripped correlator
    <phi_k1 ... phi_k4>'_conn = -2 lam Im[ prod_i phi*_{k_i}(eta) \int_{-inf}^{eta} d eta' a(eta')^4 prod_i phi_{k_i}(eta') ]
                              = -2 lam H^4 eta^4 Im[ prod_i v*_{k_i}(eta) I ],   I = \int d eta' prod_i v_{k_i}(eta').
(a^4 phi^4 = chi^4 for *any* mass/coupling, so I contains only the modes v.)
The Abel prescription is implemented by rotating the contour eta' = eta - (1+i) s, s in [0, inf).
This is a first-order, tree-level (contact) result: order lambda^1, no loops.
"""
from __future__ import annotations

import mpmath as mp
import numpy as np
from scipy.integrate import quad
from scipy.special import hankel1


# ---------------------------------------------------------------- closed forms
def conformal_contact_correlator(eta, ks, lam, H=1.0):
    """Analytic result for nu = 1/2 (m = 0, xi = 1/6):  -lam H^4 eta^4 / (8 K prod k_i), K = sum k_i."""
    ks = np.asarray(ks, dtype=float)
    return -lam * H**4 * eta**4 / (8.0 * ks.sum() * np.prod(ks))


def massless_minimal_log_coefficient(ks, lam, H=1.0):
    """Coefficient c in <phi^4>_conn ~ c ln(-eta) + const (eta -> 0^-), nu = 3/2:  c = lam H^4 sum k^3 / (24 prod k^3)."""
    ks = np.asarray(ks, dtype=float)
    return lam * H**4 * np.sum(ks**3) / (24.0 * np.prod(ks**3))


# ---------------------------------------------------------------- mode functions at complex time
def _v_complex_mp(eta_c, k, nu):
    k = mp.mpf(k)
    if nu == 0.5:
        return mp.exp(-1j * k * eta_c) / mp.sqrt(2 * k)
    if nu == 1.5:
        return mp.exp(-1j * k * eta_c) / mp.sqrt(2 * k) * (1 - 1j / (k * eta_c))
    nu_ = mp.mpf(nu)
    return mp.sqrt(mp.pi) / 2 * mp.exp(1j * mp.pi * (nu_ / 2 + mp.mpf(1) / 4)) * mp.sqrt(-eta_c) * mp.hankel1(nu_, -k * eta_c)


def _v_complex_np(eta_c, k, nu):
    if nu == 0.5:
        return np.exp(-1j * k * eta_c) / np.sqrt(2 * k)
    if nu == 1.5:
        return np.exp(-1j * k * eta_c) / np.sqrt(2 * k) * (1 - 1j / (k * eta_c))
    return 0.5 * np.sqrt(np.pi) * np.exp(1j * np.pi * (nu / 2 + 0.25)) * np.sqrt(-eta_c) * hankel1(nu, -k * eta_c)


def _s_edges(eta, K, smax_factor=45.0):
    a = abs(eta)
    smax = max(smax_factor / K, 50 * a)
    edges = [0.0]
    s = a
    while s < smax:
        edges.append(s)
        s *= 2.5
    edges.append(smax)
    return edges


def contact_integral_scipy(eta, ks, nu):
    """I = int_{-inf}^{eta} prod v_{k_i} d eta' via contour rotation (double precision)."""
    ks = list(map(float, ks))
    K = sum(ks)

    def f(s):
        z = eta - (1 + 1j) * s
        p = 1.0 + 0j
        for k in ks:
            p *= _v_complex_np(z, k, nu)
        return (1 + 1j) * p

    tot = 0j
    e = _s_edges(eta, K)
    for a, b in zip(e[:-1], e[1:]):
        re = quad(lambda s: f(s).real, a, b, limit=300, epsabs=0, epsrel=1e-12)[0]
        im = quad(lambda s: f(s).imag, a, b, limit=300, epsabs=0, epsrel=1e-12)[0]
        tot += re + 1j * im
    return tot


def contact_integral_mp(eta, ks, nu, dps=40):
    with mp.workdps(dps):
        eta_ = mp.mpf(eta)
        ks_ = [mp.mpf(k) for k in ks]
        K = sum(ks_)
        one_i = mp.mpc(1, 1)

        def f(s):
            z = eta_ - one_i * s
            p = mp.mpf(1)
            for k in ks_:
                p *= _v_complex_mp(z, k, nu)
            return one_i * p

        e = [mp.mpf(x) for x in _s_edges(float(eta), float(K))]
        return mp.quad(f, e)


def contact_correlator(eta, ks, nu, lam=1.0, H=1.0, backend="auto", dps=40):
    """First-order tree-level connected 4-point function (see module docstring).

    backend 'mp'    : extended precision (default for nu in {1/2, 3/2}: cancellation-safe at late times)
            'scipy' : double precision; loses ~ (k|eta|)^-3 relative digits in Im[...] as eta -> 0 for
                      nu > 1/2, so restrict to moderate k|eta| (documented in tests)."""
    if backend == "auto":
        backend = "mp" if nu in (0.5, 1.5) else "scipy"
    if backend == "mp":
        with mp.workdps(dps):
            I = contact_integral_mp(eta, ks, nu, dps)
            prod = mp.mpf(1)
            for k in ks:
                prod *= mp.conj(_v_complex_mp(mp.mpf(eta), k, nu))
            return float(-2 * lam * H**4 * mp.mpf(eta) ** 4 * mp.im(prod * I))
    I = contact_integral_scipy(eta, ks, nu)
    prod = np.prod([np.conj(_v_complex_np(eta + 0j, k, nu)) for k in ks])
    return float(-2 * lam * H**4 * eta**4 * np.imag(prod * I))


def contact_integral_real_axis_abel(eta, ks, nu, eps, n_max=None):
    """Independent check of the contour rotation: integrate along the REAL axis with an Abel damping factor
    exp(eps K eta') and return the integral for finite eps (nu in {1/2, 3/2}); extrapolate eps -> 0 externally."""
    ks = list(map(float, ks))
    K = sum(ks)
    T = 40.0 / (eps * K)

    def f(x):
        p = 1.0 + 0j
        for k in ks:
            p *= _v_complex_np(x, k, nu)
        return p * np.exp(eps * K * x)

    re = quad(lambda x: f(x).real, -T, eta, limit=4000, epsabs=0, epsrel=1e-11)[0]
    im = quad(lambda x: f(x).imag, -T, eta, limit=4000, epsabs=0, epsrel=1e-11)[0]
    return re + 1j * im


def schwinger_keldysh_branches(eta, ks, nu):
    """Plus (time-ordered) and minus (anti-time-ordered) branch contributions.
    X_plus = -i int <Q H_I>  (vertex on the time-ordered branch),  X_minus = +i int <H_I Q>.
    Their sum equals the commutator form -2 Im X; returned separately to display the equivalence."""
    I = contact_integral_mp(eta, ks, nu)
    prod = mp.mpf(1)
    for k in ks:
        prod *= mp.conj(_v_complex_mp(mp.mpf(eta), k, nu))
    X = prod * I  # int <H_I Q> / (lam H^4 eta^4 ...)-stripped
    return complex(-1j * mp.conj(X)), complex(1j * X)
