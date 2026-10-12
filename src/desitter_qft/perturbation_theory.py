r"""Perturbative-validity experiments and the independent 0+1-dimensional in-out check of the in-in formula.

QM analogue: H = p^2/2 + omega^2 q^2/2 + g q^4/4!.  The first-order in-in formula
    <q^4>_c^{(1)} = -2 g Im[ v*(t)^4 \int_{-inf}^t v(t')^4 dt' ] = -g/(32 omega^5)
is compared with exact diagonalisation of the ground state (adiabatic switch-on <-> true interacting ground state).
This validates the SIGN and NORMALISATION of the formula independently of any field-theory machinery.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import eigh


def qm_inin_prediction(g: float, omega: float = 1.0) -> float:
    return -g / (32.0 * omega**5)


def anharmonic_ground_state_cumulant(g: float, omega: float = 1.0, nmax: int = 120) -> float:
    """kappa_4 = <q^4> - 3<q^2>^2 in the exact ground state of H (truncated oscillator basis, size nmax)."""
    n = np.arange(nmax)
    a = np.diag(np.sqrt(n[1:]), 1)
    q = (a + a.T) / np.sqrt(2 * omega)
    q2 = q @ q
    q4 = q2 @ q2
    Hm = np.diag(omega * (n + 0.5)) + g / 24.0 * q4
    w, V = eigh(Hm)
    psi = V[:, 0]
    return float(psi @ q4 @ psi - 3 * (psi @ q2 @ psi) ** 2)


def qm_first_order_slope(omega: float = 1.0, h: float = 1e-3) -> float:
    """d kappa_4/dg at g = 0 from one-sided data (g > 0 only: for g < 0 the potential is unbounded below and the
    'ground state' of the truncated matrix is a basis artefact).  Richardson: s(h) = kappa_4(h)/h,
    slope = 2 s(h) - s(2h) removes the O(h) error."""
    s1 = anharmonic_ground_state_cumulant(h, omega) / h
    s2 = anharmonic_ground_state_cumulant(2 * h, omega) / (2 * h)
    return 2 * s1 - s2


def secular_ratio(eta_values, ks, nu, lam=1.0, H=1.0):
    """R(eta) = <phi^4>_conn^{(1)} / P_ref, P_ref = prod_i P_{k_i}^{1/2} (= P_k^2 for equal momenta), a measure of the relative size of the
    O(lambda) correction (Gaussian reference: product of two-point functions P_k).  For nu = 3/2 this grows
    ~ ln a; for light-but-massive fields it saturates or decays.  Numerical experiment, not a claim about the full theory."""
    from .in_in import contact_correlator
    from .analytic_modes import bunch_davies_mode

    out = []
    for eta in eta_values:
        # Gaussian reference scale: the disconnected 4-pt function is a sum of products of TWO spectra,
        # so P_ref = prod_i P_{k_i}^{1/2} (equal to P_k^2 for equal momenta).
        P = np.sqrt(np.prod([H**2 * eta**2 * abs(bunch_davies_mode(eta, k, nu)) ** 2 for k in ks]))
        out.append(contact_correlator(eta, ks, nu, lam, H) / P)
    return np.array(out)
