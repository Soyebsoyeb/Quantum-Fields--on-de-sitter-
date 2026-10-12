#!/usr/bin/env python
"""Assemble reports/technical_report.md. Numbers in section 17 are read from results/*.json so the report
cannot drift from the code. Run `desitter-qft all` first."""
import json, os

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
R = lambda *p: os.path.join(ROOT, *p)
val = json.load(open(R("results", "validation.json")))
fm = json.load(open(R("results", "free_modes.json")))
tp = json.load(open(R("results", "two_point.json")))
tl = json.load(open(R("results", "tree_level.json")))
env = val["environment"]

vt = "| # | Check | Category | Value | Tolerance | Result |\n|---|---|---|---|---|---|\n" + "\n".join(
    f"| {i+1} | {c['name']} | {c['category']} | {c['value']:.2e} | {c['tolerance']:.0e} | {'PASS' if c['passed'] else '**FAIL**'} |"
    for i, c in enumerate(val["checks"]))
sw = fm["tolerance_sweep"]["DOP853"]
swt = "| rtol | " + " | ".join(f"{r:.0e}" for r in sw["rtol"]) + " |\n|---|" + "---|" * len(sw["rtol"]) + "\n| max rel. error | " + " | ".join(f"{e:.1e}" for e in sw["rel_err"]) + " |\n| max \\|W−i\\| | " + " | ".join(f"{e:.1e}" for e in sw["wronskian"]) + " |"
fmt = "| Case | max \\|v_num − v_ana\\| | initial \\|W−i\\| | max \\|W−i\\| | n_fev |\n|---|---|---|---|---|\n" + "\n".join(
    f"| {n} | {d['max_abs_err_vs_analytic']:.1e} | {d['wronskian']['initial_residual']:.1e} | {d['wronskian']['max_residual']:.1e} | {d['nfev']} |"
    for n, d in fm["diagnostics"].items())
tpt = "| ν | max rel. diff, mode sum vs closed form |\n|---|---|\n" + "\n".join(f"| {k} | {v['max_rel_err_mode_sum_vs_closed_form']:.1e} |" for k, v in tp["diagnostics"].items())
dg = tl["diagnostics"]

TEXT = f"""# Quantum Fields on de Sitter Spacetime — Technical Report

*Scalar fields on a fixed four-dimensional de Sitter background: modes, vacuum structure, correlation functions, first-order in-in perturbation theory, tensor modes.*
Every statement is tagged: **[A]** analytic result, **[S]** symbolic verification, **[N]** numerical experiment, **[P]** perturbative approximation, **[C]** conceptual discussion, **[U]** unimplemented direction. Longer derivations are in `docs/`; each section names its file.

## 1. Scope and motivation

We study a quantum field on a *prescribed* de Sitter spacetime. This is quantum field theory in curved spacetime, not quantum gravity: the metric is a fixed classical background and no back-reaction is included (see §15, §18). De Sitter is the simplest maximally symmetric cosmological background and the leading-order geometry of inflation; it also exposes the basic conceptual issues of QFT in curved spacetime: the meaning of "vacuum", infrared sensitivity of light fields, and the breakdown of fixed-order perturbation theory by secular growth.

Implemented and validated: free scalar (all $\\nu$), Wightman function and its short-distance structure, mode-space Green functions, first-order tree-level contact correlator, tensor modes. Not implemented: loops, second order, renormalised $\\langle T_{{\\mu\\nu}}\\rangle$, graviton interactions (§19).

## 2. Mathematical conventions (`docs/conventions.md`)

Signature $(-+++)$, $c=\\hbar=1$; Riemann $R^\\rho{{}}_{{\\sigma\\mu\\nu}}=\\partial_\\mu\\Gamma^\\rho_{{\\nu\\sigma}}-\\dots$, $R_{{\\sigma\\nu}}=R^\\rho{{}}_{{\\sigma\\rho\\nu}}$, giving $R=+12H^2$. Conformal time $\\eta=-1/(aH)$. Fourier $f(\\mathbf x)=\\int\\frac{{d^3k}}{{(2\\pi)^3}}\\tilde f e^{{i\\mathbf k\\cdot\\mathbf x}}$. Mode normalisation $vv^{{*\\prime}}-v^*v'=i$. $\\xi R\\phi^2$ enters the action with $-\\tfrac12$; $\\xi=1/6$ is conformal. IEEE double precision; `mpmath` at 40 digits where stated; ODEs by DOP853 unless stated.

## 3. de Sitter geometry (`docs/de_sitter_geometry.md`) [A][S]

Hyperboloid $-X_0^2+\\sum X_i^2=H^{{-2}}$; flat chart $ds^2=-dt^2+e^{{2Ht}}d\\mathbf x^2$ covers half of it and is past-incomplete; static and global charts are different domains and are not interchangeable. SymPy computes $R=12H^2$ in the cosmic, conformal and static-patch charts and $R_{{\\mu\\nu}}=3H^2g_{{\\mu\\nu}}$ [S].

## 4. The scalar-field action (`docs/scalar_action.md`) [A][S]

$S=-\\frac12\\int\\sqrt{{-g}}[(\\nabla\\phi)^2+(m^2+\\xi R)\\phi^2]$. Variation with $\\delta\\phi|_\\partial=0$ gives $(\\Box-m^2-\\xi R)\\phi=0$; the Euler–Lagrange operator of a generic FRW metric equals $\\sqrt{{-g}}(\\Box-m^2-\\xi R)\\phi$ symbolically.

## 5. The Klein–Gordon equation (`docs/klein_gordon_derivation.md`) [A][S]

Cosmic time: $\\ddot\\phi+3H\\dot\\phi-\\nabla^2\\phi/a^2+(m^2+\\xi R)\\phi=0$; the friction term is $\\partial_t\\ln\\sqrt{{-g}}$. Conformal time: $\\phi''+2\\frac{{a'}}{{a}}\\phi'-\\nabla^2\\phi+a^2(m^2+\\xi R)\\phi=0$, derived from $\\Box$ and independently from the chain rule.

## 6. Fourier modes and canonical quantisation (`docs/mode_equation.md`, `docs/canonical_quantization.md`) [A]

With $v=a\\phi$: $v_k''+[k^2+a^2m^2+(6\\xi-1)a''/a]v_k=0$ → $v''+[k^2-(\\nu^2-\\frac14)/\\eta^2]v=0$, $\\nu^2=\\frac94-\\frac{{m^2}}{{H^2}}-12\\xi$. Canonical momentum $\\pi_\\chi=\\chi'$; $[\\chi,\\pi_\\chi]=i\\delta^{{(3)}}\\iff vv^{{*\\prime}}-v^*v'=i$ (constant in $\\eta$). A mode basis plus a Bogoliubov-invariant Wronskian defines an operator algebra; the Fock vacuum depends on the basis.

## 7. Analytic benchmark solutions [A]

| Case | $\\nu$ | $v_k$ | late-time $\\phi_k$ |
|---|---|---|---|
| $m=0,\\xi=\\frac16$ | $\\frac12$ | $e^{{-ik\\eta}}/\\sqrt{{2k}}$ | $\\propto\\eta$ (decays as $1/a$) |
| $m=0,\\xi=0$ | $\\frac32$ | $\\frac{{e^{{-ik\\eta}}}}{{\\sqrt{{2k}}}}(1-\\frac i{{k\\eta}})$ | $\\to H/\\sqrt{{2k^3}}$ (frozen) |
| massive | real/imag. | $\\frac{{\\sqrt\\pi}}2e^{{i\\pi(\\nu/2+1/4)}}\\sqrt{{-\\eta}}H^{{(1)}}_\\nu(-k\\eta)$ | $\\propto(-\\eta)^{{3/2\\mp\\dots}}$ |

Verified by direct substitution with `mpmath` differentiation (residual $<2\\times10^{{-13}}$) and by $W=i$ ($3\\times10^{{-13}}$) [S/N]. Regimes ($\\nu>3/2$, $=3/2$, real, $0$, imaginary) have different IR/asymptotic behaviour and are classified in `mode_equations.regime`.

## 8. Vacuum structure (`docs/vacuum_states.md`) [A][C]

Bunch–Davies: the Fock vacuum of the mode that is positive frequency at the chart's early boundary; equals the Euclidean vacuum and is $SO(1,4)$-invariant where the two-point function exists; one member of the $\\alpha$-vacuum family; for $\\nu\\ge3/2$ no invariant Hilbert-space vacuum (Allen–Folacci). A Hankel-free construction (adiabatic vacuum of order 0,2,4) converges to the same mode as $|\\eta_0|\\to\\infty$ (Fig. 3c) [N].

## 9. Wightman functions and propagators (`docs/wightman_functions.md`, `docs/propagators.md`) [A][N]

$G^+(\\eta,r)=\\int\\frac{{d^3k}}{{(2\\pi)^3}}\\phi_k\\phi_k^*e^{{i\\mathbf k\\cdot\\mathbf r}}$, defined as a distribution by Abel regularisation (subtract $1/2k$, exactly transformed). Closed form $\\frac{{H^2}}{{16\\pi^2}}\\Gamma(\\frac32+\\nu)\\Gamma(\\frac32-\\nu)\\,{{}}_2F_1(\\frac32+\\nu,\\frac32-\\nu;2;\\frac{{1+Z}}2)$, $Z=\\frac{{\\eta^2+\\eta'^2-r^2}}{{2\\eta\\eta'}}$, satisfying the invariant ODE. Mode-space $G^{{+,F,R,A}}$ with $DG^F=-i\\delta$, $DG^R=\\delta$; retarded function state-independent. $G^R$ vs an independent ODE solve: max diff ${tp['propagator']['max_diff_G_R_hankel_vs_ode']:.1e}$; derivative jump ${tp['propagator']['jump_derivative']:.10f}$.

## 10. Interacting scalar fields (`docs/interacting_scalar_fields.md`) [A][S]

$S_{{int}}=-\\int\\sqrt{{-g}}\\frac\\lambda{{4!}}\\phi^4$; $[\\lambda]=0$ (verified); $\\mathcal H_{{int}}=\\frac\\lambda{{4!}}\\chi^4$ for any $m,\\xi$ (all time dependence in the modes).

## 11. In-in perturbation theory (`docs/in_in_formalism.md`) [A][P]

$\\langle Q\\rangle=\\langle\\Omega|U^\\dagger QU|\\Omega\\rangle$, no final-state specification; first order $i\\int\\langle[H_I,Q]\\rangle$. Contour: Abel $i\\epsilon$, implemented as rotation $\\eta'=\\eta-(1+i)s$. Sign and normalisation validated by exact diagonalisation of the $0+1$-d analogue: slope ${dg['qm_slope_exact']:.8f}$ vs $-g/32\\omega^5=${dg['qm_slope_inin']:.8f}$.

## 12. A validated tree-level correlator (`docs/perturbative_correlators.md`) [P][N]

Vertex $-i\\frac\\lambda{{4!}}\\chi^4$; external modes $v_{{k_i}}$; $\\eta'\\in(-\\infty,\\eta]$; free BD vacuum; Fourier measure $(2\\pi)^3\\delta$ stripped; order $\\lambda^1$; no further regulator.
$$\\langle\\phi_{{k_1}}\\cdots\\phi_{{k_4}}\\rangle'_{{conn}}=-2\\lambda H^4\\eta^4\\,\\mathrm{{Im}}\\Big[\\prod v^*_{{k_i}}(\\eta)\\int_{{-\\infty}}^\\eta\\!d\\eta'\\prod v_{{k_i}}(\\eta')\\Big].$$
Analytic limit $\\nu=\\frac12$: $-\\lambda H^4\\eta^4/(8K\\prod k_i)$; numerical agreement ${dg['nu=1/2 max rel diff numeric vs analytic']:.1e}$ (relative). For $\\nu=\\frac32$: late-time $\\frac{{\\lambda H^4\\sum k_i^3}}{{24\\prod k_i^3}}\\ln(-\\eta)$, coefficient ${dg['nu=3/2 analytic log coefficient']:.6f}$ for the tested momenta, reproduced numerically (Fig. 9c).

## 13. Infrared and ultraviolet behaviour (`docs/infrared_structure.md`, `docs/ultraviolet_structure.md`) [A][N]

IR: $|\\phi_k|^2\\propto k^{{-2\\nu}}$; $\\int k^2|\\phi|^2dk$ converges iff $\\nu<\\frac32$; $\\nu=\\frac32$ gives $\\langle\\phi^2\\rangle=\\frac{{H^2}}{{4\\pi^2}}N$; $m\\ll H$ saturates near $3H^4/8\\pi^2m^2$ (matches the stochastic Fokker–Planck equilibrium as a benchmark). UV: $G=A/u+B\\ln u+C+O(u\\ln u)$ with $A=H^2/16\\pi^2$, $B=(\\frac94-\\nu^2-2)A$, $C=B[\\psi(\\frac32+\\nu)+\\psi(\\frac32-\\nu)+2\\gamma_E-1]$ (DLMF 15.8.10; $B$ cross-checked by Frobenius recursion; all three checked numerically). Leading singularity = flat space. Hadamard condition described [C]; sharp cutoffs shown to give $\\Lambda^2$ divergence and are not renormalisation.

## 14. Renormalisation issues (`docs/renormalization.md`) [C] + one explicit subtraction

One explicitly defined point-splitting subtraction ($1/4\\pi^2\\sigma^2$ and $B\\ln\\mu^2\\sigma^2$): $\\langle\\phi^2\\rangle_\\mu=A/3+C+B\\ln(H^2/4\\mu^2)$. It is **not** presented as a complete covariant renormalised quantity: the $\\mu$-dependence and the choice of parametrix ($\\Delta^{{1/2}}$ shifts the constant by $H^2/16\\pi^2$) are finite ambiguities; no $\\langle T_{{\\mu\\nu}}\\rangle$ and no loop renormalisation of $\\lambda$ is attempted.

## 15. Tensor perturbations (`docs/tensor_perturbations.md`) [S][A][C]

Linearised vacuum Einstein equations (SymPy) give $h_k''-\\frac2\\eta h_k'+k^2h_k=0$. With the standard quadratic action (quoted), each polarisation is a massless minimal scalar and $\\Delta_t^2=2H^2/\\pi^2M_p^2$; numerical integration gives {dg['tensor_Dt2_final']:.7f} vs {dg['tensor_Dt2_analytic']:.7f} ($H=M_p=1$, $x=10^{{-3}}$). Quantising a test scalar does not quantise gravity; gauge fixing, constraints, non-renormalisability and the EFT interpretation are conceptual only.

## 16. Numerical methodology [N]

Mode ODE: real 4-vector state, DOP853, `rtol` $10^{{-4}}$–$10^{{-13}}$, `atol`$\\approx10^{{-2}}$rtol. Integration $\\eta_0=-80\\to\\eta_1=-0.02$, $k=1$. Hankel functions: `scipy.special` (real order) / `mpmath` (imaginary order). In-in integrals: contour rotation with `scipy.quad` (double) or `mpmath.quad` (40 digits). Distributional $k$-integrals: QUADPACK QAWF. Environment: Python {env['python']}, NumPy {env['numpy']}, SciPy {env['scipy']}, SymPy {env['sympy']}, mpmath {env['mpmath']}, Matplotlib {env['matplotlib']}. Each experiment writes a JSON record (parameters, diagnostics, data hash, environment) to `results/`.

## 17. Validation results

Programmatic suite (`desitter-qft validate`; also 80+ `pytest` tests including negative controls):

{vt}

Free-mode runs ($k=1$, rtol $10^{{-11}}$, atol $10^{{-13}}$, DOP853; `results/free_modes.json`):

{fmt}

Solver-tolerance convergence (massless minimal, $\\eta\\in[-100,-0.05]$):

{swt}

The error scales linearly with rtol over nine decades before saturating near $10^{{-12}}$ (round-off). Wronskian residual is not rescaled; it grows with $|v||v'|$ (up to $2\\times10^{{-7}}$ for the tensor mode at $\\eta=-10^{{-3}}$).

Wightman mode sum vs closed form ($\\eta=-1$):

{tpt}

Figures: Fig. 1 horizon crossing; 2 modes vs analytic; 3 solver/vacuum convergence; 4 power spectra; 5 Wightman; 6 UV structure; 7 IR variance; 8 propagators; 9 contact correlator; 10 secular growth; 11 in-in vs exact diagonalisation; 12 tensor spectrum (`figures/`).

## 18. Limitations (`docs/limitations.md`)

Fixed background only; flat chart only; first-order tree-level interaction only; scheme-dependent point splitting; IR cutoff dependence for $\\nu\\ge3/2$; double-precision cancellation at late times for $\\nu\\ne\\frac12,\\frac32$ (mitigated by `mpmath`); tensor-action normalisation quoted; references not network-verified. Numerical agreement validates, but does not prove, the analytic statements.

## 19. Future research directions [U]

(i) Second-order and exchange diagrams, one-loop tadpole with covariant Hadamard renormalisation; (ii) $\\alpha$-vacua and their Hadamard violation; (iii) stochastic/Fokker–Planck resummation compared with fixed-order in-in; (iv) global-chart and static-patch (thermal) analysis; (v) $\\langle T_{{\\mu\\nu}}\\rangle$ and back-reaction; (vi) scalar and vector metric perturbations with constraints; (vii) graviton interactions and effective-field-theory loops; (viii) Euclidean/sphere methods and dS/CFT-like boundary correlators.

## 20. References

See `references.bib`. Key sources: Birrell–Davies 1982; Chernikov–Tagirov 1968; Bunch–Davies 1978; Allen 1985; Allen–Folacci 1987; Mottola 1985; Starobinsky 1986; Starobinsky–Yokoyama 1994; Schwinger 1961; Keldysh 1964; Weinberg 2005; Maldacena 2003; Radzikowski 1996; Hollands–Wald 2001; Wald 1994; Mukhanov–Feldman–Brandenberger 1992; Donoghue 1994; DLMF. The bibliography was written from bibliographic knowledge and has not been machine-verified against a database in this session.

### Source/assumption/check table

| Result | Source | Assumptions | Check |
|---|---|---|---|
| Hankel BD mode | Bunch–Davies 1978; DLMF 10.17 | $\\nu$ fixed by $m,\\xi,H$; chart boundary positive frequency | substitution, $W=i$, ODE integration, adiabatic vacuum |
| Closed-form $G^+$ | Chernikov–Tagirov 1968 | BD state, spacelike $Z<1$ | invariant ODE, mode sum, $\\nu=\\frac12$ limit |
| $(A,B,C)$ | DLMF 15.8.10 | $\\nu$ generic | numerics at $u\\to0$, Frobenius for $B$ |
| In-in first order | Schwinger 1961; Weinberg 2005 | adiabatic switch-on, free BD vacuum | exact diagonalisation of QM analogue |
| Contact correlator | derived here | $O(\\lambda)$, connected | analytic $\\nu=\\frac12$; contour vs Abel; log slope |
| Tensor equation | derived here (SymPy) | TT plus-polarisation, vacuum + $\\Lambda$ | symbolic residual 0 |
| $\\Delta_t^2=2H^2/\\pi^2M_p^2$ | Grishchuk 1974; Starobinsky 1979; MFB 1992 | quoted action normalisation | numerical ODE spectrum |
"""
open(R("reports", "technical_report.md"), "w").write(TEXT)
print("wrote report", len(TEXT))
