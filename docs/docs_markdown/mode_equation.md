# Rescaled mode equation and analytic benchmarks

*Status: derivations = symbolic verification (`mode_equations.derive_rescaled_mode_equation`); solutions = analytic results checked by substitution (`tests/test_analytic_modes.py`) and against numerical integration.*

## Derivation

Put $\phi=v/a$ and $\nabla^2\to-k^2$ in the conformal-time equation. Using $\phi'=v'/a-a'v/a^2$,
$$\phi''+2\tfrac{a'}{a}\phi'=\frac{v''}{a}-\frac{a''}{a^2}v.$$
Inserting $R=6a''/a^3$ so that $a^2\xi R=6\xi\,a''/a$:
$$v_k''+\Big[k^2+a^2m^2+(6\xi-1)\frac{a''}{a}\Big]v_k=0.$$
**Why rescale?** Writing $\chi=a\phi$ in the action and integrating by parts removes the friction term:
$$S=\frac12\int d\eta\,d^3x\Big[\chi'^2-(\nabla\chi)^2-\big(a^2m^2+(6\xi-1)\tfrac{a''}{a}\big)\chi^2\Big].$$
This is a flat-space scalar with a time-dependent mass: canonical momentum $\pi_\chi=\chi'$, so canonical commutators and the Wronskian are those of a unit-mass oscillator.

## de Sitter form

With $a''/a=2/\eta^2$ and $\nu^2=\frac94-\frac{m^2}{H^2}-12\xi$,
$$v_k''+\Big[k^2-\frac{\nu^2-\frac14}{\eta^2}\Big]v_k=0.$$
(The identity $a^2m^2+(6\xi-1)a''/a=-(\nu^2-\tfrac14)/\eta^2$ is verified symbolically.)

## Solution

The general solution is $\sqrt{-\eta}\,[c_1H^{(1)}_\nu(-k\eta)+c_2H^{(2)}_\nu(-k\eta)]$. Since $H^{(1)}_\nu(x)\sim\sqrt{2/\pi x}\,e^{i(x-\nu\pi/2-\pi/4)}$ (DLMF 10.17.5), positive frequency $e^{-ik\eta}$ at $k\vert\eta\vert\to\infty$ selects
$$\boxed{v_k=\frac{\sqrt\pi}{2}e^{i\pi(\nu/2+1/4)}\sqrt{-\eta}\,H^{(1)}_\nu(-k\eta)}$$
with $v v^{\ast\prime}-v^{\ast}v'=i$ (checked to $10^{-13}$).

## Regimes of $\nu$ (not interchangeable)

| $\nu^2$ | $\nu$ | Late-time behaviour of $\phi_k$ | IR integrability of $k^2\vert\phi_k\vert^2$ |
|---|---|---|---|
| $>9/4$ | real, $>3/2$ | growth $\propto(-\eta)^{3/2-\nu}$ | **not** integrable: $\vert\phi_k\vert^2\sim k^{-2\nu}$ |
| $=9/4$ | $3/2$ | frozen, $\vert\phi_k\vert^2\to H^2/2k^3$ | marginal: logarithmic IR divergence |
| $(0,9/4)$ | real | decay $(-\eta)^{3/2-\nu}$ (complementary series) | integrable |
| $=0$ | 0 | $(-\eta)^{3/2}\ln(-\eta)$ (limit) | integrable |
| $<0$ | $i\mu$ | $(-\eta)^{3/2}\cos(\mu\ln(-\eta)+\alpha)$ (principal series) | integrable |

Imaginary $\nu$ uses `mpmath` (complex-order Hankel); $\nu=0$ and the boundary $\nu=3/2$ are tested as limiting cases.

## Benchmarks

* **Conformal, $m=0,\xi=1/6$** ($\nu=1/2$). $6\xi-1=0$ and $m=0$: the time-dependent term vanishes, $v''+k^2v=0$, $v=e^{-ik\eta}/\sqrt{2k}$. $\phi=v/a=-H\eta\,v$. This is the statement that the massless conformally coupled field is Weyl-equivalent to a flat-space massless field; the in-in vertex $a^4\phi^4=\chi^4$ is then time-independent.
* **Massless minimal, $m=0,\xi=0$** ($\nu=3/2$). $v''+(k^2-2/\eta^2)v=0$, $v=\frac{e^{-ik\eta}}{\sqrt{2k}}\big(1-\frac{i}{k\eta}\big)$, $\phi_k=\frac{H}{\sqrt{2k^3}}(1+ik\eta)e^{-ik\eta}\to\frac{H}{\sqrt{2k^3}}$. Super-horizon modes freeze with a scale-invariant spectrum $\Delta^2_\phi=H^2/4\pi^2$ (Fig. 4). **IR subtlety:** $\int dk\,k^2\vert\phi_k\vert^2\sim\int dk/k$ diverges at small $k$. The formal mode expansion does not define a de Sitter-invariant state: there is no $SO(1,4)$-invariant Hilbert-space vacuum for this field (Allen–Folacci 1987); see `infrared_structure.md`.
* **Massive, general $\nu$.** As above; light fields ($\nu<3/2$) decay super-horizon, with $\phi_k\propto(-\eta)^{3/2-\nu}$, i.e. $\Delta^2_\phi\propto x^{3-2\nu}$ (grey asymptotes in Fig. 4).

## Numerical comparison

`solve_mode` integrates the real ODE from $\eta_0=-80$ to $-0.02$ using exact BD data; max abs error $\lesssim2\times10^{-9}$ (Fig. 2), error $\propto$ rtol over 8 decades (Fig. 3a). Using *adiabatic* (WKB) initial data of order 0/2/4 instead of Hankel data, the numerical solution converges to the analytic BD mode as $\vert\eta_0\vert\to\infty$ with higher order giving faster convergence (Fig. 3c): this is a Hankel-free confirmation that the positive-frequency boundary condition selects the same state.
