#!/usr/bin/env python
"""Generate and execute the notebooks. Notebooks only *call* the package (no essential logic lives here).
Usage: python scripts/build_notebooks.py [--no-execute]"""
import os, sys
import nbformat as nbf

HERE = os.path.dirname(os.path.abspath(__file__))
ND = os.path.join(HERE, "..", "notebooks")

SETUP = """import sys, os
sys.path.insert(0, os.path.abspath(os.path.join('..', 'src')))
import numpy as np, sympy as sp, matplotlib.pyplot as plt
from IPython.display import Image, display
from desitter_qft.plotting import apply_style; apply_style()
FIG = os.path.abspath(os.path.join('..', 'figures'))"""

NB = {}
NB["01_geometry_and_curvature"] = ("de Sitter geometry and curvature", "*Symbolic verification.* $R=+12H^2$ in three charts with the convention of `docs/conventions.md`.", [
"""from desitter_qft import backgrounds as bg
H = sp.Symbol('H', positive=True)
for name, f in [('flat (cosmic)', bg.flat_chart_metric_cosmic), ('flat (conformal)', bg.flat_chart_metric_conformal), ('static patch', bg.static_patch_metric)]:
    g, X = f(); print(f'{name:18s} R =', bg.ricci_scalar(g, X))
g, X = bg.flat_chart_metric_conformal()
print('R_mn - 3H^2 g_mn =', sp.simplify(bg.ricci_tensor(g, X) - 3*H**2*g) == sp.zeros(4, 4))""",
"""dS = bg.DeSitter(1.0)
t = np.array([0.0, 1.0, 5.0]); eta = dS.conformal_time(t)
print('eta(t) =', eta, ' a(eta) =', dS.scale_factor_conformal(eta), ' a(t) =', dS.scale_factor_cosmic(t))
display(Image(os.path.join(FIG, 'fig01_geometry_horizon.png')))"""])
NB["02_klein_gordon_derivation"] = ("Klein-Gordon equation from the action", "*Symbolic verification* of `docs/scalar_action.md` and `docs/klein_gordon_derivation.md`.", [
"""from desitter_qft import scalar_action as sa
EL, KG, d = sa.euler_lagrange_vs_klein_gordon(); print('E-L minus sqrt(-g)(Box - m^2 - xi R)phi =', d)
d, R, phi, a, _ = sa.conformal_time_equation(); print('conformal eq residual =', d, '; R =', R)
print('cosmic-time residual, R:', sa.cosmic_time_equation_difference())
print('chain-rule cross-check:', sa.cosmic_to_conformal_chain_rule())"""])
NB["03_free_mode_equations"] = ("Rescaled mode equation", "*Symbolic verification + analytic result.*", [
"""from desitter_qft import mode_equations as me
r = me.derive_rescaled_mode_equation(); print('v = a*phi substitution residual:', r['difference'])
print('de Sitter reduction residual:', me.de_sitter_reduction())""",
"""rows = [('conformal', 0, 1/6), ('massless minimal', 0, 0), ('m=0.8H', 0.8, 0), ('m=1.7H', 1.7, 0), ('xi=-0.1, m=0', 0, -0.1)]
for n, m, xi in rows:
    print(f'{n:18s} nu^2 = {me.nu_squared(m, xi):+.4f}  nu = {me.nu(m, xi):.4f}  ->  {me.regime(me.nu_squared(m, xi))}')
eta = np.linspace(-5, -0.5, 4); print('conformal omega^2 (should be k^2=1):', me.omega_squared(eta, 1.0, 0.0, 1/6))"""])
NB["04_analytic_mode_solutions"] = ("Analytic Bunch-Davies modes", "*Analytic results checked by substitution.*", [
"""import mpmath as mp
from desitter_qft.analytic_modes import *
from desitter_qft.normalization import wronskian_residual
eta = np.linspace(-40, -0.01, 400)
for nu in (0.5, 1.0, 1.2689, 1.5, 0.8j):
    v = bunch_davies_mode(eta, 1.0, nu); dv = bunch_davies_mode_derivative(eta, 1.0, nu)
    print(f'nu={nu}: max|W-i| = {wronskian_residual(v, dv).max():.2e}')
nu, k = 1.2689, 1.0
f = lambda e: 0.5*mp.sqrt(mp.pi)*mp.exp(1j*mp.pi*(nu/2+0.25))*mp.sqrt(-e)*mp.hankel1(nu, -k*e)
print('ODE residual by substitution:', [float(abs(mp.diff(f, e, 2) + (k**2 - (nu**2-0.25)/e**2)*f(e))) for e in (-3.0, -0.7, -0.1)])""",
"""display(Image(os.path.join(FIG, 'fig04_power_spectra.png')))"""])
NB["05_numerical_mode_evolution"] = ("Numerical mode evolution", "*Numerical experiment.* Parameters are printed with each result.", [
"""from desitter_qft import experiments as ex
fm = ex.free_modes_experiment(n=300)
for n, r in fm.items(): print(f"{n:32s} max|v_num - v_ana| = {r['max_abs_err']:.2e}   max|W-i| = {r['wronskian_report']['max_residual']:.2e}")""",
"""sw = ex.tolerance_sweep(rtols=(1e-4, 1e-6, 1e-8, 1e-10), methods=('DOP853',))
print('rtol  :', sw['rtols']); print('error :', sw['results']['DOP853']['err'])""",
"""ad = ex.adiabatic_convergence(eta0s=np.array([-5., -12., -35., -100.]), orders=(0, 2, 4), cases=(('massless minimal', 0.0, 0.0),))
for k, e in ad['errors'].items(): print(k, e)
display(Image(os.path.join(FIG, 'fig03_solver_convergence_and_vacuum.png')))"""])
NB["06_canonical_quantization"] = ("Canonical quantisation and the Wronskian", "*Analytic derivation + numerical normalisation check.*", [
"""from desitter_qft.analytic_modes import *
from desitter_qft.normalization import wronskian, wronskian_residual
eta = np.linspace(-20, -0.1, 5)
v = bunch_davies_mode(eta, 1.0, 1.5); dv = bunch_davies_mode_derivative(eta, 1.0, 1.5)
print('W (should be i, constant):', wronskian(v, dv))
print('2v has residual', wronskian_residual(2*v, 2*dv)[0], '(wrong normalisation is flagged, never silently rescaled)')""",
"""# Bogoliubov transform u = alpha v + beta v*: Wronskian preserved iff |alpha|^2 - |beta|^2 = 1
ch, sh = np.cosh(0.7), np.sinh(0.7)
u = ch*v + sh*np.conj(v); du = ch*dv + sh*np.conj(dv)
print('Bogoliubov (|a|^2-|b|^2=1):', wronskian(u, du)[0]); 
u2 = 1.3*v + 0.4*np.conj(v); du2 = 1.3*dv + 0.4*np.conj(dv)
print('violating the condition:', wronskian(u2, du2)[0])"""])
NB["07_vacuum_structure"] = ("Vacuum structure and infrared behaviour", "*Analytic + numerical.* IR integrability of $k^2|\\phi_k|^2$ as a function of $\\nu$.", [
"""from desitter_qft.two_point import mode_power_phi, phi2_cutoff, phi2_cutoff_minimal_massless_analytic
eta = -1.0
for nu in (0.5, 1.0, 1.3, 1.5, 1.7):
    k = np.array([1e-4, 1e-3]); p = k**2*mode_power_phi(eta, k, nu)
    slope = np.log(p[1]/p[0])/np.log(k[1]/k[0])
    print(f'nu={nu}: d ln(k^2|phi|^2)/d ln k = {slope:+.3f}   (integrable at k->0 iff slope > -1; expected 2-2nu = {2-2*nu:+.2f})')""",
"""print('cutoff <phi^2> (nu=3/2): numeric', phi2_cutoff(-0.5, 50, 1e-3, 1.5), ' closed form', phi2_cutoff_minimal_massless_analytic(-0.5, 50, 1e-3))
display(Image(os.path.join(FIG, 'fig07_infrared_variance.png')))"""])
NB["08_wightman_functions"] = ("Wightman functions", "*Numerical mode sum (Abel-regularised) vs analytic closed form; short-distance structure.*", [
"""from desitter_qft.two_point import *
for nu in (0.5, 1.0, 1.3, 0.8j):
    for r in (0.3, 2.0):
        a = wightman_mode_sum(-1.0, r, nu); b = wightman_closed_form_equal_time(-1.0, r, nu)
        print(f'nu={nu}, r={r}: mode sum {a:.12e}  closed form {b:.12e}  rel diff {abs(a/b-1):.1e}')""",
"""for nu in (1.0, 1.4):
    A, B, C = short_distance_coefficients(nu); u = 1e-6
    print(f'nu={nu}: A,B,C = {A:.5f}, {B:.5f}, {C:.5f}; residual at u=1e-6: {float(wightman_closed_form_u(u, nu)) - (A/u + B*np.log(u) + C):.2e}')
display(Image(os.path.join(FIG, 'fig06_ultraviolet_structure.png'))); display(Image(os.path.join(FIG, 'fig08_mode_space_propagators.png')))"""])
NB["09_interacting_scalar_fields"] = ("Interacting scalar field", "*Symbolic verification.*", [
"""from desitter_qft.interactions import *
Hint, expected = interaction_hamiltonian_density(); print('H_int =', Hint, '; equals -L_int:', sp.simplify(Hint - expected) == 0)
terms, ok = check_dimensions(); print(terms, '\\nall terms marginal, lambda dimensionless:', ok)"""])
NB["10_in_in_perturbation_theory"] = ("In-in perturbation theory: validation of the first-order formula", "*Perturbative approximation vs exact diagonalisation (0+1 d analogue).*", [
"""from desitter_qft.perturbation_theory import *
for w in (1.0, 2.0):
    print(f'omega={w}: exact slope {qm_first_order_slope(w):.8f}   in-in -g/(32 w^5) = {qm_inin_prediction(1.0, w):.8f}')""",
"""from desitter_qft.in_in import schwinger_keldysh_branches
ks = [1.0, 1.3, 0.7, 1.1]; plus, minus = schwinger_keldysh_branches(-0.5, ks, 0.5)
print('plus + minus =', plus + minus, ' expected', -1/(8*sum(ks)*np.prod(ks)))
display(Image(os.path.join(FIG, 'fig11_inin_vs_exact_diagonalisation.png')))"""])
NB["11_tree_level_correlators"] = ("Tree-level correlators", "*Perturbative approximation, $O(\\lambda)$.* Contact diagram.", [
"""from desitter_qft.in_in import *
ks = [1.0, 1.3, 0.7, 1.1]
for eta in (-2.0, -0.3):
    print(eta, 'numeric', contact_correlator(eta, ks, 0.5), ' analytic', conformal_contact_correlator(eta, ks, 1.0))""",
"""c = massless_minimal_log_coefficient(ks, 1.0)
e1, e2 = -1e-4, -1e-5
print('slope', (contact_correlator(e2, ks, 1.5) - contact_correlator(e1, ks, 1.5))/(np.log(-e2)-np.log(-e1)), ' analytic', c)
from desitter_qft.perturbation_theory import secular_ratio
N = np.array([1., 3., 5.]); 
for nu in (0.5, 1.0, 1.5): print(nu, secular_ratio(-np.exp(-N), [1.0]*4, nu))
display(Image(os.path.join(FIG, 'fig09_tree_level_contact_correlator.png'))); display(Image(os.path.join(FIG, 'fig10_secular_growth.png')))"""])
NB["12_tensor_perturbations"] = ("Tensor perturbations", "*Symbolic verification of the tensor equation; numerical spectrum.*", [
"""from desitter_qft.tensor_modes import *
r = derive_tensor_equation(); print(r)""",
"""from desitter_qft.experiments import tensor_experiment
d = tensor_experiment(); print('Delta_t^2 numerical', d['final_value'], ' analytic 2H^2/(pi^2 Mp^2) =', d['late'])
display(Image(os.path.join(FIG, 'fig12_tensor_spectrum.png')))"""])

def build(execute=True):
    os.makedirs(ND, exist_ok=True)
    for name, (title, intro, cells) in NB.items():
        nb = nbf.v4.new_notebook()
        nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
        nb.cells = [nbf.v4.new_markdown_cell(f"# {title}\n\n{intro}\n\nAll computations call `desitter_qft`; see the matching `docs/` file for derivations."),
                    nbf.v4.new_code_cell(SETUP)] + [nbf.v4.new_code_cell(c) for c in cells]
        path = os.path.join(ND, name + ".ipynb")
        if execute:
            from nbclient import NotebookClient
            NotebookClient(nb, timeout=900, kernel_name="python3", resources={"metadata": {"path": ND}}).execute()
        nbf.write(nb, path)
        print("wrote", path)

if __name__ == "__main__":
    build(execute="--no-execute" not in sys.argv)
