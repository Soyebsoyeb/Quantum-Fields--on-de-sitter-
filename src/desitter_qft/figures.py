"""Figure generation. Every figure is produced from the functions in experiments.py (no hidden state)."""
from __future__ import annotations

import json
import os

import numpy as np

from . import experiments as ex
from .backgrounds import DeSitter
from .plotting import CYCLE, OI, apply_style, panel_label, plt, save
from .provenance import hash_arrays, write_record

LAB = {0.5: r"$\nu=1/2$", 1.0: r"$\nu=1$", 1.3: r"$\nu=1.3$", 1.25: r"$\nu=5/4$", 1.2689: r"$\nu\approx1.27$",
       1.5: r"$\nu=3/2$", 0.8j: r"$\nu=0.8\,i$"}


def _lab(nu):
    return LAB.get(nu, str(nu))


def fig01_geometry(outdir):
    apply_style()
    dS = DeSitter(1.0)
    N = np.linspace(-4, 10, 400)
    eta = -np.exp(-N)
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9))
    ax[0].plot(N, np.log10(dS.comoving_hubble_radius(eta)), color="k", label=r"$1/(aH)=-\eta$")
    for k, c in zip((1e-1, 1.0, 1e2), CYCLE):
        ax[0].axhline(np.log10(1 / k), color=c, ls="--", lw=1.0)
        ax[0].text(9.9, np.log10(1 / k) + 0.12, rf"$1/k$, $k={k:g}\,H$", color=c, ha="right", fontsize=8)
    ax[0].set(xlabel=r"$N=\ln a$", ylabel=r"$\log_{10}$ comoving length $[H^{-1}]$")
    ax[0].legend(loc="lower left")
    panel_label(ax[0], "(a)")
    ax[0].text(0.55, 0.62, "super-horizon\n$k|\\eta|<1$", transform=ax[0].transAxes, fontsize=8, color="0.3")
    t = np.linspace(-2, 8, 300)
    ax[1].plot(t, dS.conformal_time(t), color=CYCLE[0], label=r"$\eta(t)=-e^{-Ht}/H$")
    ax[1].plot(t, dS.scale_factor_cosmic(t) * 1e-3, color=CYCLE[1], label=r"$a(t)\times10^{-3}$")
    ax[1].set(xlabel=r"$Ht$", ylabel=r"$H\eta,\; a/10^{3}$", ylim=(-1.2, 3.0))
    ax[1].legend(loc="upper left")
    panel_label(ax[1], "(b)")
    ax[1].text(0.97, 0.12, r"$\eta\to0^-$ as $t\to\infty$" "\n" r"(future boundary $\mathscr{I}^+$)", transform=ax[1].transAxes,
               ha="right", fontsize=8)
    save(fig, outdir, "fig01_geometry_horizon")


def fig02_modes(outdir, data=None):
    apply_style()
    d = data or ex.free_modes_experiment()
    fig, axs = plt.subplots(2, 4, figsize=(11, 4.6), gridspec_kw={"height_ratios": [2.2, 1]}, sharex="col")
    for j, (name, r) in enumerate(d.items()):
        a, b = axs[0, j], axs[1, j]
        a.plot(r["eta"], r["v_ana"].real, color="k", lw=2.4, alpha=0.35, label="analytic")
        a.plot(r["eta"], r["v_num"].real, color=CYCLE[0], lw=1.0, ls="--", label="DOP853")
        a.set_title(name.replace("xi", r"$\xi$").replace("m=", r"$m=$"), fontsize=8.5, pad=14)
        a.set_xlim(-12, -0.02) if True else None
        mask = r["eta"] > -12
        lim = np.max(np.abs(r["v_ana"][mask].real)) * 1.1
        a.set_ylim(-min(lim, 6), min(lim, 6))
        if j == 0:
            a.set_ylabel(r"$\mathrm{Re}\,v_k(\eta)$")
            a.legend(loc="lower left")
        b.semilogy(r["eta"], np.abs(r["v_num"] - r["v_ana"]) + 1e-18, color=CYCLE[1], lw=1)
        b.semilogy(r["eta"], np.abs(r["v_num"] - r["v_ana"]) * 0 + r["wronskian_report"]["max_residual"], color="0.5", ls=":", lw=1)
        b.set_xlabel(r"$k\eta$")
        if j == 0:
            b.set_ylabel(r"$|v_{\rm num}-v_{\rm ana}|$")
        b.set_ylim(1e-14, 1e-5)
        panel_label(a, "(" + "abcd"[j] + ")")
    axs[1, 3].text(0.97, 0.9, "dotted: max Wronskian\nresidual", transform=axs[1, 3].transAxes, ha="right", va="top", fontsize=7.5, color="0.4")
    fig.tight_layout()
    save(fig, outdir, "fig02_mode_functions")


def fig03_solver(outdir, sweep=None, adi=None):
    apply_style()
    sw = sweep or ex.tolerance_sweep()
    ad = adi or ex.adiabatic_convergence()
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.3))
    for (meth, r), c in zip(sw["results"].items(), CYCLE):
        axs[0].loglog(sw["rtols"], r["err"], "o-", color=c, ms=4, label=meth)
        axs[1].loglog(sw["rtols"], r["wronskian"], "o-", color=c, ms=4, label=meth)
    axs[0].loglog(sw["rtols"], sw["rtols"] * 10, color="0.4", ls=":", label=r"$\propto$ rtol")
    axs[0].set(xlabel="rtol", ylabel=r"max rel. error vs analytic $v_k$")
    axs[1].set(xlabel="rtol", ylabel=r"max $|W-i|$")
    axs[0].legend()
    axs[0].invert_xaxis(); axs[1].invert_xaxis()
    panel_label(axs[0], "(a)"); panel_label(axs[1], "(b)")
    mk = {0: "o", 2: "s", 4: "^"}
    for (name, o), e in ad["errors"].items():
        if name != "massless minimal":
            continue
        axs[2].loglog(-ad["eta0"], np.maximum(e, 1e-17), mk[o] + "-", ms=4, label=f"adiabatic order {o}")
    x = -ad["eta0"]
    axs[2].loglog(x, 3 * x ** -4.0, color="0.4", ls=":", label=r"$\propto|\eta_0|^{-4}$")
    axs[2].loglog(x, 30 * x ** -6.0, color="0.6", ls="--", label=r"$\propto|\eta_0|^{-6}$")
    axs[2].set(xlabel=r"$k|\eta_0|$", ylabel=r"$\left|\,|v|^2/|v_{\rm BD}|^2-1\right|$ at $k\eta=-0.5$", ylim=(1e-16, 1))
    axs[2].legend(fontsize=7)
    panel_label(axs[2], "(c)")
    fig.tight_layout()
    save(fig, outdir, "fig03_solver_convergence_and_vacuum")
    return sw, ad


def fig04_power(outdir, data=None):
    from .analytic_modes import late_time_amplitude_squared
    apply_style()
    d = data or ex.power_spectrum_experiment()
    fig, ax = plt.subplots(figsize=(5.0, 3.7))
    x = d["x"]
    for (nu, P), c in zip(d["spectra"].items(), CYCLE):
        ax.loglog(x, P, color=c, label=_lab(nu))
        if abs(complex(nu).imag) < 1e-12 and nu.real > 0.5 if isinstance(nu, complex) else (nu > 0.5):
            xs = x[x < 0.03]
            asym = np.array([late_time_amplitude_squared(1.0, float(nu), 1.0, -xi) for xi in xs]) / (2 * np.pi**2)
            ax.loglog(xs, asym, color="k", ls="-", lw=3.2, alpha=0.22, zorder=0)
    ax.axhline(1 / (4 * np.pi**2), color="k", ls="--", lw=0.8)
    ax.text(1.3e-3, 1 / (4 * np.pi**2) * 1.35, r"$H^2/4\pi^2$", fontsize=8)
    ax.set(xlabel=r"$x=-k\eta$", ylabel=r"$\Delta^2_\phi=\frac{k^3}{2\pi^2}|\phi_k|^2\;[H^2]$", ylim=(1e-8, 1e2))
    ax.legend(loc="lower right", title="grey: small-$x$ Hankel asymptote", title_fontsize=7)
    save(fig, outdir, "fig04_power_spectra")


def fig05_two_point(outdir, tp=None):
    apply_style()
    d = tp or ex.two_point_experiment()
    fig, axs = plt.subplots(1, 2, figsize=(8.4, 3.4))
    for (nu, r), c in zip(d["data"].items(), CYCLE):
        axs[0].loglog(r["r"], r["closed"], color=c, label=_lab(nu))
        axs[0].loglog(r["r"], r["mode_sum"], "o", color=c, ms=3.2, mfc="none")
        axs[1].loglog(r["r"], np.maximum(r["rel_err"], 1e-17), "o-", color=c, ms=3, lw=0.8)
    axs[0].set(xlabel=r"$r\,[H^{-1}]$ (comoving)", ylabel=r"$G^+(\eta,r)\;[H^2]$, $H\eta=-1$")
    axs[0].legend(title="line: closed form\nmarkers: mode sum", fontsize=7.5)
    axs[1].set(xlabel=r"$r\,[H^{-1}]$", ylabel="relative difference")
    panel_label(axs[0], "(a)"); panel_label(axs[1], "(b)")
    fig.tight_layout()
    save(fig, outdir, "fig05_wightman_modesum_vs_closed_form")


def fig06_uv(outdir, tp=None, uv=None):
    apply_style()
    d = tp or ex.two_point_experiment()
    u = uv or ex.ultraviolet_cutoff_experiment()
    fig, axs = plt.subplots(1, 3, figsize=(11.5, 3.3))
    for (nu, v), c in zip(d["had"].items(), CYCLE):
        axs[0].semilogx(d["r_had"], v, color=c, label=_lab(nu))
    axs[0].axhline(1, color="k", ls=":")
    axs[0].set(xlabel=r"comoving $r$", ylabel=r"$4\pi^2 d^2\,G^+$, $d=a r$", ylim=(0.9, 1.1))
    axs[0].legend(fontsize=7.5)
    for (nu, v), c in zip(d["coef"].items(), CYCLE):
        if nu == 0.5:
            continue  # B = C = 0 identically: residual is round-off only
        axs[1].loglog(d["u"], np.abs(v["residual"]) + 1e-30, color=c, label=_lab(nu))
    axs[1].loglog(d["u"], 0.01 * d["u"] * np.abs(np.log(d["u"])), color="0.4", ls=":", label=r"$\propto u\ln u$")
    axs[1].set(xlabel=r"$u=(1-Z)/2$", ylabel=r"$|G-A/u-B\ln u-C|\;[H^2]$")
    axs[1].legend(fontsize=7.5, loc="upper left")
    axs[2].loglog(u["k_uv"], u["numeric"], "o", ms=3.5, mfc="none", label="numerical $k$-integral")
    axs[2].loglog(u["k_uv"], u["analytic"], color="k", lw=1, label="closed form")
    axs[2].loglog(u["k_uv"], u["eta"] ** 2 * u["k_uv"] ** 2 / (8 * np.pi**2), color=OI["vermillion"], ls="--", lw=1, label=r"$\eta^2\Lambda^2/8\pi^2$")
    axs[2].set(xlabel=r"comoving cutoff $\Lambda$", ylabel=r"$\langle\phi^2\rangle_\Lambda\;[H^2]$ ($\nu=3/2$)")
    axs[2].legend(fontsize=7.5)
    for a_, s in zip(axs, "abc"):
        panel_label(a_, f"({s})")
    fig.tight_layout()
    save(fig, outdir, "fig06_ultraviolet_structure")


def fig07_ir(outdir, ir=None):
    apply_style()
    d = ir or ex.ir_experiment()
    fig, ax = plt.subplots(figsize=(5.2, 3.7))
    for (lab, v), c in zip(d["curves"].items(), CYCLE):
        ax.plot(d["N"], v["var"] / 1.0, color=c, label=lab)
        if v["m"] > 0:
            ax.axhline(3 / (8 * np.pi**2 * v["m"] ** 2), color=c, ls=":", lw=1)
    ax.plot(d["N"], d["N"] / (4 * np.pi**2), color="k", ls="--", lw=0.9, label=r"$N/4\pi^2$")
    ax.set(xlabel=r"$N=\ln(1/k_{\rm IR}|\eta|)$ e-folds", ylabel=r"$\langle\phi^2\rangle_{k_{\rm IR}<k<1/|\eta|}\;[H^2]$", ylim=(0, 1.0))
    ax.legend(loc="upper left", fontsize=7.5, title="dotted: $3H^4/8\\pi^2m^2$", title_fontsize=7.5)
    save(fig, outdir, "fig07_infrared_variance")


def fig08_propagators(outdir, pr=None):
    apply_style()
    d = pr or ex.propagator_experiment()
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.2))
    axs[0].plot(d["eta"], d["GR"], color="k", lw=2.5, alpha=0.35, label=r"$-2\,\mathrm{Im}\,v(\eta)v^*(\eta')$")
    axs[0].plot(d["eta"], d["GR_ode"], color=CYCLE[0], ls="--", label="ODE, $y(\\eta')=0,\\,y'(\\eta')=1$")
    axs[0].set(xlabel=r"$\eta$", ylabel=r"$G^R_k(\eta,\eta')$")
    axs[0].legend(fontsize=7.5, loc="lower left")
    axs[0].set_ylim(-8, 2)
    axs[1].semilogy(d["eta"], np.abs(d["GR"] - d["GR_ode"]) + 1e-18, color=CYCLE[1])
    axs[1].set(xlabel=r"$\eta$", ylabel=r"$|G^R_{\rm Hankel}-G^R_{\rm ODE}|$")
    axs[2].plot(d["grid"], d["GF"].real, color=CYCLE[0], label=r"Re $G^F$")
    axs[2].plot(d["grid"], d["GF"].imag, color=CYCLE[1], label=r"Im $G^F$")
    axs[2].plot(d["grid"], d["GW"].real, color="k", ls=":", label=r"Re $G^+$")
    axs[2].axvline(d["etap"], color="0.5", lw=0.8)
    axs[2].set(xlabel=r"$\eta$", ylabel=r"$G_k(\eta,\eta')$")
    axs[2].legend(fontsize=7.5)
    axs[0].text(0.97, 0.95, rf"$\partial_\eta G^R|_{{\eta'^+}}={d['jump_derivative']:.6f}$", transform=axs[0].transAxes, ha="right", va="top", fontsize=8)
    for a_, s in zip(axs, "abc"):
        panel_label(a_, f"({s})")
    fig.tight_layout()
    save(fig, outdir, "fig08_mode_space_propagators")


def fig09_tree(outdir, tr=None):
    apply_style()
    d = tr or ex.tree_level_experiment()
    fig, axs = plt.subplots(1, 3, figsize=(11.5, 3.4))
    axs[0].loglog(-d["eta"], np.abs(d["conf_ana"]), color="k", lw=2.5, alpha=0.35, label=r"$\lambda H^4\eta^4/(8K\prod k_i)$")
    axs[0].loglog(-d["eta"], np.abs(d["conf_num"]), "o", color=CYCLE[0], ms=3.5, mfc="none", label="in-in contour integral")
    axs[0].set(xlabel=r"$-\eta\,[1/k]$", ylabel=r"$|\langle\phi^4\rangle'_{\rm conn}|$ ($\nu=1/2$)")
    axs[0].legend(fontsize=7.5)
    Ic = d["I_contour"]
    axs[1].plot(d["eps"], d["I_abel"].real, "o-", color=CYCLE[0], label=r"Re $I(\epsilon)$ (real axis)")
    axs[1].plot(d["eps"], d["I_abel"].imag, "s-", color=CYCLE[1], label=r"Im $I(\epsilon)$")
    axs[1].axhline(Ic.real, color=CYCLE[0], ls=":"); axs[1].axhline(Ic.imag, color=CYCLE[1], ls=":")
    axs[1].set(xlabel=r"Abel damping $\epsilon$", ylabel=r"$I=\int_{-\infty}^{\eta}\prod v_{k_i}\,d\eta'$", xlim=(0, 0.09))
    axs[1].legend(fontsize=7.5, title="dotted: contour value", title_fontsize=7.5, loc="lower left")
    axs[2].plot(d["N"], d["corr32"], "o", color=CYCLE[0], ms=3.5, label=r"numerical ($\nu=3/2$)")
    # analytic slope line anchored at the last point
    Np = d["N"]
    axs[2].plot(Np, d["corr32"][-1] + d["log_coeff"] * (-(Np - Np[-1])), color="k", lw=1, label=r"slope $-c\,N$, $c=\lambda H^4\sum k^3/24\prod k^3$")
    axs[2].set(xlabel=r"$N=\ln(1/k|\eta|)$", ylabel=r"$\langle\phi^4\rangle'_{\rm conn}\;[H^4]$")
    axs[2].legend(fontsize=7.5)
    for a_, s in zip(axs, "abc"):
        panel_label(a_, f"({s})")
    fig.tight_layout()
    save(fig, outdir, "fig09_tree_level_contact_correlator")


def fig10_secular(outdir, tr=None):
    apply_style()
    d = tr or ex.tree_level_experiment()
    fig, ax = plt.subplots(figsize=(5.0, 3.6))
    for (nu, r), c in zip(d["ratios"].items(), CYCLE):
        ax.plot(d["Ns"], r, "o-", ms=3.5, color=c, label=_lab(nu))
    ax.set(xlabel=r"$N=\ln(1/k|\eta|)$ e-folds after horizon exit",
           ylabel=r"$\langle\phi^4\rangle^{(1)}_{\rm conn}/P_k^{2}\;[\lambda]$ (all $k_i=k$)")
    ax.legend(loc="center right", fontsize=8)
    ax.text(0.03, 0.03, "linear growth in $N$ at $\\nu=3/2$:\nsecular term, fixed-order\nexpansion degrades when\n$\\lambda N\\sim1$", transform=ax.transAxes, ha="left", va="bottom", fontsize=7.5)
    ax.set_ylim(-5.4, 0.6)
    save(fig, outdir, "fig10_secular_growth")


def fig11_qm(outdir, q=None):
    apply_style()
    d = q or ex.qm_inin_experiment()
    fig, axs = plt.subplots(1, 2, figsize=(8.0, 3.2))
    axs[0].plot(d["g"], d["exact"], "o", color=CYCLE[0], ms=3.5, mfc="none", label="exact diagonalisation")
    axs[0].plot(d["g"], d["first_order"], color="k", lw=1, label=r"in-in, $O(g)$: $-g/32\omega^5$")
    axs[0].set(xlabel=r"$g/\omega^3$", ylabel=r"$\kappa_4=\langle q^4\rangle-3\langle q^2\rangle^2\;[\omega^{-2}]$")
    axs[0].legend(fontsize=7.5)
    m = d["g"] > 0
    rel = np.abs(d["exact"][m] - d["first_order"][m]) / np.abs(d["first_order"][m])
    axs[1].loglog(d["g"][m], rel, "o-", ms=3.5, color=CYCLE[1])
    axs[1].loglog(d["g"][m], 4 * d["g"][m], color="0.4", ls=":", label=r"$\propto g$ (next order)")
    axs[1].set(xlabel=r"$g/\omega^3$", ylabel="relative deviation from first order")
    axs[1].legend(fontsize=7.5)
    panel_label(axs[0], "(a)"); panel_label(axs[1], "(b)")
    fig.tight_layout()
    save(fig, outdir, "fig11_inin_vs_exact_diagonalisation")


def fig12_tensor(outdir, t=None):
    apply_style()
    d = t or ex.tensor_experiment()
    fig, axs = plt.subplots(1, 2, figsize=(8.0, 3.3))
    axs[0].loglog(d["x"], d["Dt2_ana"], color="k", lw=2.5, alpha=0.35, label="analytic $\\nu=3/2$")
    axs[0].loglog(d["x"], d["Dt2_num"], color=CYCLE[0], ls="--", label="numerical ODE")
    axs[0].axhline(d["late"], color=CYCLE[1], ls=":", label=r"$2H^2/\pi^2M_p^2$")
    axs[0].set(xlabel=r"$x=-k\eta$", ylabel=r"$\Delta_t^2(k,\eta)$ [$H=M_p=1$]")
    axs[0].legend(fontsize=7.5)
    axs[1].loglog(d["x"], np.abs(d["Dt2_num"] / d["Dt2_ana"] - 1) + 1e-18, color=CYCLE[2])
    axs[1].set(xlabel=r"$x=-k\eta$", ylabel="relative deviation (num/analytic)")
    panel_label(axs[0], "(a)"); panel_label(axs[1], "(b)")
    fig.tight_layout()
    save(fig, outdir, "fig12_tensor_spectrum")


def _record(outdir, name, payload: dict, data: dict):
    os.makedirs(os.path.join(os.path.dirname(os.path.normpath(outdir)), "results"), exist_ok=True)
    rd = os.path.join(os.path.dirname(os.path.normpath(outdir)), "results")
    write_record(os.path.join(rd, name + ".json"), {**payload, "data_sha256": hash_arrays(data)})


def make_free_mode_figures(outdir="figures"):
    fm = ex.free_modes_experiment()
    fig01_geometry(outdir)
    fig02_modes(outdir, fm)
    sw, ad = fig03_solver(outdir)
    fig04_power(outdir)
    _record(outdir, "free_modes", {
        "experiment": "free_modes", "category": "numerical experiment",
        "parameters": {"k": 1.0, "eta0": -80.0, "eta1": -0.02, "rtol": 1e-11, "atol": 1e-13, "method": "DOP853",
                       "cases": list(ex.CASES)},
        "diagnostics": {n: {"max_abs_err_vs_analytic": r["max_abs_err"], "wronskian": r["wronskian_report"], "nfev": r["stats"]["nfev"]} for n, r in fm.items()},
        "tolerance_sweep": {m: {"rtol": sw["rtols"].tolist(), "rel_err": r["err"].tolist(), "wronskian": r["wronskian"].tolist()} for m, r in sw["results"].items()},
    }, {n: {"eta": r["eta"], "v": r["v_num"]} for n, r in fm.items()})


def make_two_point_figures(outdir="figures"):
    tp = ex.two_point_experiment()
    fig05_two_point(outdir, tp)
    uv = ex.ultraviolet_cutoff_experiment()
    fig06_uv(outdir, tp, uv)
    fig07_ir(outdir)
    pr = ex.propagator_experiment()
    fig08_propagators(outdir, pr)
    _record(outdir, "two_point", {
        "experiment": "two_point", "category": "numerical experiment (mode sum) checked against analytic closed form",
        "parameters": {"eta": tp["eta"], "H": 1.0, "regularisation": "Abel (subtract conformal 1/2k, exact distributional integral)"},
        "diagnostics": {str(nu): {"max_rel_err_mode_sum_vs_closed_form": float(np.max(v["rel_err"]))} for nu, v in tp["data"].items()},
        "propagator": {"max_diff_G_R_hankel_vs_ode": pr["max_diff"], "jump_derivative": pr["jump_derivative"]},
    }, {str(nu): {"ms": v["mode_sum"]} for nu, v in tp["data"].items()})


def make_tree_figures(outdir="figures"):
    tr = ex.tree_level_experiment()
    fig09_tree(outdir, tr)
    fig10_secular(outdir, tr)
    q = ex.qm_inin_experiment()
    fig11_qm(outdir, q)
    t = ex.tensor_experiment()
    fig12_tensor(outdir, t)
    _record(outdir, "tree_level", {
        "experiment": "tree_level_contact_correlator", "category": "perturbative approximation (O(lambda), tree level) + numerical evaluation",
        "parameters": {"momenta": tr["ks"], "lambda": 1.0, "H": 1.0, "prescription": "Abel / contour rotation eta' = eta-(1+i)s, free Bunch-Davies vacuum"},
        "diagnostics": {
            "nu=1/2 max rel diff numeric vs analytic": float(np.max(np.abs(tr["conf_num"] / tr["conf_ana"] - 1))),
            "nu=3/2 analytic log coefficient": tr["log_coeff"],
            "qm_slope_exact": q["slope"], "qm_slope_inin": q["predicted_slope"],
            "tensor_Dt2_final": t["final_value"], "tensor_Dt2_analytic": t["late"],
        },
    }, {"conf": tr["conf_num"], "corr32": tr["corr32"], "ratios": {str(k): v for k, v in tr["ratios"].items()}})
