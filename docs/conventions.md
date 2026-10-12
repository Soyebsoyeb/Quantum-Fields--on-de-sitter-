# Conventions

All formulas, code and figures in this repository use the conventions below. Where a source uses another convention, the translation is stated where the source is used.

## Metric, units, curvature

| Item | Convention |
|---|---|
| Signature | $(-+++)$ |
| Units | $c=\hbar=1$; $[H]=[m]=$ mass, $[t]=[\eta]=[x]=$ mass$^{-1}$, $[a]=0$ |
| Christoffel | $\Gamma^\lambda_{\mu\nu}=\tfrac12 g^{\lambda\sigma}(\partial_\mu g_{\sigma\nu}+\partial_\nu g_{\sigma\mu}-\partial_\sigma g_{\mu\nu})$ |
| Riemann | $R^\rho{}_{\sigma\mu\nu}=\partial_\mu\Gamma^\rho_{\nu\sigma}-\partial_\nu\Gamma^\rho_{\mu\sigma}+\Gamma^\rho_{\mu\lambda}\Gamma^\lambda_{\nu\sigma}-\Gamma^\rho_{\nu\lambda}\Gamma^\lambda_{\mu\sigma}$ |
| Ricci | $R_{\sigma\nu}=R^\rho{}_{\sigma\rho\nu}$, $R=g^{\sigma\nu}R_{\sigma\nu}$ |
| Result | de Sitter has $R=+12H^2$, $R_{\mu\nu}=3H^2g_{\mu\nu}$, $\Lambda=3H^2$ (verified symbolically in three charts: `tests/test_curvature.py`) |
| d'Alembertian | $\Box\phi=\frac{1}{\sqrt{-g}}\partial_\mu\big(\sqrt{-g}\,g^{\mu\nu}\partial_\nu\phi\big)$ |
| Curvature coupling | action term $-\tfrac12\xi R\phi^2$; $\xi=0$ minimal, $\xi=1/6$ conformal in $d=4$. A positive $\xi R$ is a positive contribution to the effective mass squared. |

## Time variables

Cosmic time $t$, $a=e^{Ht}$. Conformal time $d\eta=dt/a$, i.e. $\eta=-1/(Ha)=-e^{-Ht}/H\in(-\infty,0)$ and $a(\eta)=-1/(H\eta)$.
Overdots are $d/dt$, primes $d/d\eta$. $a'/a=-1/\eta$, $a''/a=2/\eta^2$, $\mathcal H\equiv a'/a=aH$.

## Fourier transform and mode expansion

$$f(\mathbf x)=\int\frac{d^3k}{(2\pi)^3}\,\tilde f(\mathbf k)\,e^{i\mathbf k\cdot\mathbf x},\qquad \tilde f(\mathbf k)=\int d^3x\,f(\mathbf x)\,e^{-i\mathbf k\cdot\mathbf x}.$$

Operator expansion of the rescaled field $\chi=a\phi$:
$$\hat\chi(\eta,\mathbf x)=\int\frac{d^3k}{(2\pi)^3}\big[\hat a_{\mathbf k}\,v_k(\eta)\,e^{i\mathbf k\cdot\mathbf x}+\hat a^\dagger_{\mathbf k}\,v_k^*(\eta)\,e^{-i\mathbf k\cdot\mathbf x}\big],\quad [\hat a_{\mathbf k},\hat a^\dagger_{\mathbf k'}]=(2\pi)^3\delta^{(3)}(\mathbf k-\mathbf k').$$
Normalisation: $v_kv_k^{*\prime}-v_k^*v_k'=i$. Primed correlators $\langle\cdots\rangle'$ have $(2\pi)^3\delta^{(3)}(\sum\mathbf k)$ stripped.
The Bunch–Davies phase is fixed by $v_k\to e^{-ik\eta}/\sqrt{2k}$ as $k|\eta|\to\infty$; global phases never enter observables.

## Dimensions

| Quantity | Mass dimension |
|---|---|
| $\phi$ | 1 |
| $\chi=a\phi$, $v_k$ | 1 and $-1/2$ ($v_k$ per $\sqrt{(2\pi)^3\delta}$ normalisation) |
| $\lambda$ in $\lambda\phi^4/4!$ | 0 (checked in `interactions.check_dimensions`) |
| $\Delta^2_\phi=\frac{k^3}{2\pi^2}|\phi_k|^2$ | 2 |
| $M_p$ | 1 |

## Numerical conventions

* IEEE binary64 unless stated. `mpmath` (40 digits) is used for (i) Hankel functions of imaginary order, (ii) hypergeometric closed forms, (iii) the late-time in-in integrals where cancellations destroy double precision.
* ODE solver: `scipy.integrate.solve_ivp`, default `DOP853`, `rtol=1e-10…1e-12`, `atol≈1e-2·rtol` or `1e-13`. State vector $(\mathrm{Re}\,v,\mathrm{Im}\,v,\mathrm{Re}\,v',\mathrm{Im}\,v')$.
* Every solution is reported with the Wronskian residual $|v v^{*\prime}-v^*v'-i|$. **Nothing is rescaled after integration.**
* Mode-function numerics use $k=1$; the equation depends on $k$ only through $k\eta$ (scaling $v_k(\eta)=k^{-1/2}f(k\eta)$).
* Result classification (used in docs, docstrings and `results/*.json`): (1) analytic result, (2) symbolic verification, (3) numerical experiment, (4) perturbative approximation, (5) conceptual discussion, (6) unimplemented direction.
