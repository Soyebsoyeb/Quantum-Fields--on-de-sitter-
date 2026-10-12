"""Reproducible numerical experiments. Each function returns plain arrays/dicts (no plotting) so that tests,
scripts and notebooks share one implementation. All randomness-free; results are deterministic.
"""
from __future__ import annotations

import numpy as np

from .analytic_modes import (bunch_davies_mode, bunch_davies_mode_derivative, field_mode,
                             minimal_massless_mode)
from .in_in import (conformal_contact_correlator, contact_correlator, massless_minimal_log_coefficient,
                    contact_integral_mp, contact_integral_real_axis_abel)
from .infrared import free_massive_stochastic_phi2, ir_variance
from .initial_conditions import adiabatic_data, exact_bd_data
from .mode_equations import ModeEquation, nu as nu_of
from .normalization import check_normalization, wronskian_residual
from .numerical_solver import solve_mode
from .perturbation_theory import (anharmonic_ground_state_cumulant, qm_first_order_slope, qm_inin_prediction,
                                  secular_ratio)
from .propagators import feynman_k, retarded_from_ode, retarded_k, wightman_k
from .tensor_modes import dimensionless_tensor_spectrum, tensor_spectrum_late_time
from .two_point import (coincident_phi2_pointsplit, invariant_ode_residual, mode_power_phi,
                        phi2_cutoff, phi2_cutoff_minimal_massless_analytic, short_distance_coefficients,
                        wightman_closed_form_equal_time, wightman_closed_form_u, wightman_mode_sum)

CASES = {
    "conformal  (m=0, xi=1/6)": (0.0, 1.0 / 6.0),
    "massless minimal (m=0, xi=0)": (0.0, 0.0),
    "light massive (m=0.8H, xi=0)": (0.8, 0.0),
    "heavy (m=1.7H, xi=0)": (1.7, 0.0),
}


def free_modes_experiment(k=1.0, eta0=-80.0, eta1=-0.02, rtol=1e-11, atol=1e-13, n=900):
    out = {}
    for name, (m, xi) in CASES.items():
        eq = ModeEquation(k, m, xi)
        v0, dv0 = exact_bd_data(eq, eta0)
        s = solve_mode(eq, eta0, eta1, v0, dv0, n_out=n, rtol=rtol, atol=atol)
        va = bunch_davies_mode(s.eta, k, eq.nu)
        rep = check_normalization(s.eta, s.v, s.dv, rtol=rtol, atol=atol, method="DOP853")
        out[name] = dict(eta=s.eta, v_num=s.v, v_ana=va, nu=eq.nu, nu2=eq.nu2, wronskian_report=rep.as_dict(),
                         max_abs_err=float(np.max(np.abs(s.v - va))), stats=s.stats)
    return out


def tolerance_sweep(rtols=(1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-13), k=1.0, m=0.0, xi=0.0, eta0=-100.0, eta1=-0.05,
                    methods=("RK45", "DOP853", "Radau")):
    eq = ModeEquation(k, m, xi)
    v0, dv0 = exact_bd_data(eq, eta0)
    eta = np.linspace(eta0, eta1, 600)
    ref = bunch_davies_mode(eta, k, eq.nu)
    res = {}
    for meth in methods:
        errs, wmax, nfev = [], [], []
        for rt in rtols:
            s = solve_mode(eq, eta0, eta1, v0, dv0, method=meth, rtol=rt, atol=rt * 1e-2, eta_eval=eta)
            errs.append(float(np.max(np.abs(s.v - ref) / np.abs(ref))))
            wmax.append(float(wronskian_residual(s.v, s.dv).max()))
            nfev.append(s.stats["nfev"])
        res[meth] = dict(err=np.array(errs), wronskian=np.array(wmax), nfev=np.array(nfev))
    return dict(rtols=np.array(rtols), results=res)


def adiabatic_convergence(eta0s=None, orders=(0, 2, 4), cases=(("massless minimal", 0.0, 0.0), ("light massive", 0.8, 0.0)),
                          eta_f=-0.5, k=1.0):
    eta0s = np.array([-3, -5, -8, -12, -20, -35, -60, -100, -200, -400.0]) if eta0s is None else eta0s
    out = {}
    for name, m, xi in cases:
        eq = ModeEquation(k, m, xi)
        ref = abs(bunch_davies_mode(eta_f, k, eq.nu)) ** 2
        for o in orders:
            e = []
            for eta0 in eta0s:
                v0, dv0 = adiabatic_data(eq, eta0, o)
                s = solve_mode(eq, eta0, eta_f, v0, dv0, rtol=1e-12, atol=1e-14, n_out=2)
                e.append(abs(abs(s.v[-1]) ** 2 / ref - 1))
            out[(name, o)] = np.array(e)
    return dict(eta0=eta0s, errors=out)


def power_spectrum_experiment(nus=(0.5, 1.0, 1.2689, 1.5, 0.8j), x=None, H=1.0):
    x = np.geomspace(1e-3, 1e2, 400) if x is None else x
    out = {}
    for nu in nus:
        eta = -x  # k = 1
        P = 1.0 / (2 * np.pi**2) * np.abs(field_mode(eta, 1.0, nu, H)) ** 2
        out[nu] = P
    return dict(x=x, spectra=out)


def two_point_experiment(nus=(0.5, 1.0, 1.3, 0.8j), eta=-1.0, H=1.0):
    data = {}
    for nu in nus:
        r_max = 5.0 if abs(complex(nu).imag) > 0 else 20.0  # mpmath-Hankel path: tail accuracy ~1e-11 absolute
        r = np.geomspace(0.05, r_max, 25)
        ms = np.array([wightman_mode_sum(eta, ri, nu, H) for ri in r])
        cf = np.array([wightman_closed_form_equal_time(eta, ri, nu, H) for ri in r])
        data[nu] = dict(r=r, mode_sum=ms, closed=cf, rel_err=np.abs(ms - cf) / np.abs(cf))
    # Hadamard leading singularity: 4 pi^2 d^2 G -> 1, d = a r proper distance
    a = -1.0 / (H * eta)
    rs = np.geomspace(1e-4, 1.0, 40)
    had = {}
    for nu in nus:
        had[nu] = np.array([4 * np.pi**2 * (a * ri) ** 2 * wightman_closed_form_equal_time(eta, ri, nu, H) for ri in rs])
    # next-order structure: coefficients A, B, C vs numerics at small u
    us = np.geomspace(1e-8, 1e-2, 30)
    coef = {}
    for nu in nus:
        A, B, C = short_distance_coefficients(nu, H)
        G = np.array([float(wightman_closed_form_u(u, nu, H)) for u in us])
        coef[nu] = dict(A=A, B=B, C=C, residual=G - (A / us + B * np.log(us) + C))
    return dict(eta=eta, data=data, r_had=rs, had=had, u=us, coef=coef)


def ultraviolet_cutoff_experiment(eta=-0.5, k_ir=1e-3, H=1.0):
    kuv = np.geomspace(2.0, 2e3, 24)
    num = np.array([phi2_cutoff(eta, ku, k_ir, 1.5, H) for ku in kuv])
    ana = np.array([phi2_cutoff_minimal_massless_analytic(eta, ku, k_ir, H) for ku in kuv])
    return dict(k_uv=kuv, numeric=num, analytic=ana, eta=eta, k_ir=k_ir)


def propagator_experiment(k=1.0, m=0.0, xi=0.0, etap=-6.0, H=1.0):
    eq = ModeEquation(k, m, xi, H)
    v = lambda e: bunch_davies_mode(np.asarray(e, float), k, eq.nu)
    eta = np.linspace(etap + 1e-9, -0.05, 400)
    GR = retarded_k(v, eta, etap)
    GR_ode = retarded_from_ode(eq, etap, eta)
    grid = np.linspace(-8, -0.2, 300)
    GF = feynman_k(v, grid, etap)
    Wp = wightman_k(v, grid, etap)
    # derivative jump of G_R at eta -> eta'^+
    d = 1e-6
    jump = float((retarded_k(v, etap + d, etap) - 0.0) / d)
    return dict(eta=eta, GR=GR, GR_ode=GR_ode, max_diff=float(np.max(np.abs(GR - GR_ode))), grid=grid, GF=GF, GW=Wp,
                jump_derivative=jump, etap=etap)


def tree_level_experiment(H=1.0, lam=1.0):
    ks = [1.0, 1.3, 0.7, 1.1]
    etas = -np.geomspace(0.02, 20, 18)
    conf_num = np.array([contact_correlator(e, ks, 0.5, lam, H) for e in etas])
    conf_ana = np.array([conformal_contact_correlator(e, ks, lam, H) for e in etas])
    # contour vs real-axis Abel extrapolation
    eta0 = -0.5
    Ic = complex(contact_integral_mp(eta0, ks, 1.5))
    eps = np.array([0.08, 0.04, 0.02, 0.01])
    Ia = np.array([contact_integral_real_axis_abel(eta0, ks, 1.5, e) for e in eps])
    # log growth for nu=3/2
    N = np.linspace(0.5, 8.0, 16)
    eta_log = -np.exp(-N)
    corr32 = np.array([contact_correlator(e, ks, 1.5, lam, H) for e in eta_log])
    c = massless_minimal_log_coefficient(ks, lam, H)
    # secular ratio for several nu (equal momenta) -- numerical experiment
    kk = [1.0] * 4
    Ns = np.linspace(0.3, 6.0, 12)
    et = -np.exp(-Ns)
    ratios = {}
    for nu in (0.5, 1.0, 1.25, 1.5):
        ratios[nu] = secular_ratio(et, kk, nu, lam, H)
    return dict(ks=ks, eta=etas, conf_num=conf_num, conf_ana=conf_ana, eps=eps, I_abel=Ia, I_contour=Ic,
                N=N, corr32=corr32, log_coeff=c, Ns=Ns, ratios=ratios)


def qm_inin_experiment():
    g = np.geomspace(1e-3, 3.0, 16)  # g > 0 only (bounded potential)
    exact = np.array([anharmonic_ground_state_cumulant(gi, 1.0, 140) for gi in g])
    return dict(g=g, exact=exact, first_order=np.array([qm_inin_prediction(gi, 1.0) for gi in g]),
                slope=qm_first_order_slope(1.0), predicted_slope=qm_inin_prediction(1.0, 1.0))


def tensor_experiment(H=1.0, Mp=1.0, k=1.0):
    eq = ModeEquation(k, 0.0, 0.0, H)
    eta0, eta1 = -150.0, -1e-3
    v0, dv0 = adiabatic_data(eq, eta0, 4)
    x = np.geomspace(1e-3, 150.0, 500)
    s = solve_mode(eq, eta0, eta1, v0, dv0, rtol=1e-12, atol=1e-14, eta_eval=np.sort(-x))
    eta = s.eta
    P = 2.0 / Mp**2 * H**2 * eta**2 * np.abs(s.v) ** 2
    Dt2 = k**3 * 4.0 * P / (2 * np.pi**2)
    return dict(x=-eta, Dt2_num=Dt2, Dt2_ana=dimensionless_tensor_spectrum(eta, k, H, Mp),
                late=tensor_spectrum_late_time(H, Mp), final_value=float(Dt2[-1]),
                wronskian=check_normalization(s.eta, s.v, s.dv, rtol=1e-12).as_dict())


def ir_experiment(k_ir=1.0, H=1.0):
    Ns = np.linspace(0.0, 40.0, 41)
    out = {}
    for label, (m, xi) in {"nu=3/2 (m=0)": (0.0, 0.0), "m=0.1H": (0.1, 0.0), "m=0.3H": (0.3, 0.0), "m=0.6H": (0.6, 0.0)}.items():
        nu = nu_of(m, xi, H).real
        vals = []
        for N in Ns:
            eta = -np.exp(-N) / k_ir
            vals.append(ir_variance(eta, k_ir, nu, H) if N > 0 else 0.0)
        out[label] = dict(var=np.array(vals), m=m)
    return dict(N=Ns, curves=out)
