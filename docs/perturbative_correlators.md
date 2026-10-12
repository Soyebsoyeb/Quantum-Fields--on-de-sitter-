# Tree-level correlators: the contact trispectrum

*Status: perturbative approximation (order $\lambda$, tree level, no loops), analytic + numerical.*

## Observable

The connected, momentum-stripped four-point function of the field at fixed time,
$$\langle\phi_{\mathbf k_1}\phi_{\mathbf k_2}\phi_{\mathbf k_3}\phi_{\mathbf k_4}\rangle'_{\rm conn}(\eta)=-2\lambda H^4\eta^4\,\mathrm{Im}\Big[\prod_{i=1}^4v^*_{k_i}(\eta)\,I(\eta)\Big],\quad I=\int_{-\infty}^{\eta}d\eta'\prod_iv_{k_i}(\eta').$$

**Derivation.** $H_I=\frac\lambda{4!}\int d^3x\,\chi^4$. Contract the four $\chi$'s of the vertex with the four external fields: $4!$ ways, cancelling the $1/4!$; each contraction gives $\langle\chi_{\mathbf k}(\eta')\chi_{\mathbf k}(\eta)\rangle'=v_k(\eta')v_k^*(\eta)$ (vertex to the left). Momentum conservation from $\int d^3x$. Convert $\chi=a\phi$ for the four external legs: $a^{-4}=H^4\eta^4$. Contractions of the vertex with itself are tadpole/disconnected and are not part of the connected correlator.

**Specification.** Vertex $-i\frac\lambda{4!}\chi^4$ (contact); external modes $v_{k_i}$; time range $(-\infty,\eta]$; vacuum free BD, $i\epsilon$ (Abel) via contour rotation $\eta'=\eta-(1+i)s$; Fourier normalisation $(2\pi)^3\delta^{(3)}$ stripped; order $\lambda^1$; no regulator needed (the integral converges on the rotated contour). Implemented in `in_in.contact_correlator`.

## Analytic limit (validation)

$\nu=1/2$: $\prod v=e^{-iK\eta'}/(4\sqrt{\prod k})$, $K=\sum k_i$, $I=ie^{-iK\eta}/(K\,4\sqrt{\prod k})$:
$$\boxed{\langle\phi^4\rangle'_{\rm conn}=-\frac{\lambda H^4\eta^4}{8K\prod_ik_i}}$$
The numerical contour integral reproduces this to $\sim10^{-16}$ at all tested $\eta$ (Fig. 9a); the $\eta^4$ scaling is $a^{-4}$ — the interaction does not evolve in conformal variables.

## Cross-checks of the machinery

1. **Contour rotation vs real-axis Abel damping.** $I_\epsilon=\int e^{\epsilon K\eta'}\prod v$ along the real axis tends linearly in $\epsilon$ to the rotated-contour value (Fig. 9b).
2. **Exact diagonalisation** of the QM analogue (above).
3. **Double vs extended precision:** the two backends agree to $10^{-8}$ at $\nu=1$.

## Massless minimal field: logarithm

For $\nu=3/2$, expand $\phi_k(\eta')=\frac{H}{\sqrt{2k^3}}\big(1+\tfrac{k^2\eta'^2}2-\tfrac{ik^3\eta'^3}3+\dots\big)$. The term $\propto\eta'^{3}$ in $a^4\prod\phi$ gives $\int d\eta'/\eta'$, so
$$\langle\phi^4\rangle'_{\rm conn}\ \xrightarrow{\eta\to0}\ \frac{\lambda H^4\sum k_i^3}{24\prod k_i^3}\ln(-\eta)+\text{const}.$$
Numerically, the slope with respect to $\ln(-\eta)$ equals this coefficient to better than $10^{-6}$ (relative) at $|\eta|\sim10^{-5}/k$ (Fig. 9c). This is a *secular* term: $\propto\lambda N$ with $N=\ln(1/k|\eta|)$.

## Perturbative validity (numerical experiment, Fig. 10)

Ratio $R=\langle\phi^4\rangle^{(1)}_{\rm conn}/P_k^2$ (equal momenta, $P_k=|\phi_k|^2$) vs $N$: $\nu=1/2$ constant $(-1/8)\lambda$ (conformal invariance); $\nu=1$ and $5/4$ saturate; $\nu=3/2$ grows linearly in $N$. Consequently for light fields the fixed-order result becomes unreliable once $\lambda N\sim1$, i.e. $N\sim1/\lambda$ e-folds; for the massive fields the saturation scale is set by $1/(3/2-\nu)\sim3H^2/2m^2$. **This is a limitation of the fixed-order expansion for these observables, not evidence that the theory is inconsistent**: resummation (Starobinsky stochastic dynamics, Starobinsky–Yokoyama 1994; Hollands 2013; Marolf–Morrison 2011) controls the IR sector. Those methods are described but not derived here.
