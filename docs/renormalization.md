# Renormalisation

*Status: conceptual discussion + one explicitly specified point-splitting subtraction. **No claim of a fully renormalised, covariant result.***

## Conceptual role

Divergences of $\langle\phi^2\rangle$, $\langle T_{\mu\nu}\rangle$ and of loop corrections are local: they are absorbed by counterterms built from local curvature invariants and the couplings ($m^2$, $\xi$, $\lambda$, $\Lambda_{\rm cc}$, $G$, higher-curvature terms). A *renormalisation condition* fixes the finite parts; different conditions differ by finite local terms, e.g. $\langle\phi^2\rangle\to\langle\phi^2\rangle+c_1m^2+c_2(\xi-\tfrac16)R+c_3\,\mu\text{-dependent}$. These ambiguities are physical inputs (measured couplings), not defects, but a number quoted without its scheme is meaningless. Covariant schemes: Hadamard point splitting (Wald 1977; Hollands–Wald 2001; Birrell–Davies 1982).

## The example implemented: $\langle\phi^2\rangle$ by point splitting

Definition (`two_point.coincident_phi2_pointsplit`):
$$\langle\phi^2\rangle_\mu\equiv\lim_{\sigma\to0}\Big[G^+(\sigma)-\frac1{4\pi^2\sigma^2}-B\ln(\mu^2\sigma^2)\Big].$$
Using $A/u=\frac{4A}{H^2\sigma^2}+\frac A3+O(\sigma^2)$ and $\ln u=\ln(H^2\sigma^2/4)+O(\sigma^2)$:
$$\langle\phi^2\rangle_\mu=\frac A3+C+B\ln\frac{H^2}{4\mu^2}.$$
Checks: for $\nu=1/2$ ($B=C=0$) $\langle\phi^2\rangle=A/3=H^2/48\pi^2$ for *any* $\mu$; for $m\to0$ with $\xi=0$, $\to3H^4/8\pi^2m^2$ (test), reproducing the leading IR physics.

**Finite renormalisation ambiguities (explicit):**
1. $\mu$: shifts the result by $B\ln(\mu^2/\mu_0^2)$ (tested).
2. Subtracting the full Hadamard parametrix with $\Delta^{1/2}$ rather than the bare $1/(4\pi^2\sigma^2)$ shifts the constant by $A=H^2/16\pi^2$ ($\Delta^{1/2}=(H\sigma/\sin H\sigma)^{3/2}\simeq1+H^2\sigma^2/4$ in $d=4$ gives a constant offset) — a different but equally legitimate scheme.
3. A term $\propto m^2$ or $(\xi-\frac16)R$ from the choice of renormalisation condition for $m^2,\xi$.

**What is NOT done:** a covariant (non-coordinate) Hadamard subtraction with the full local $V$ recursion; $\langle T_{\mu\nu}\rangle$; any loop renormalisation of $\lambda$ (beta function $\beta_\lambda=3\lambda^2/16\pi^2$ in flat space is quoted knowledge, not computed here; on dS curvature corrections and the IR sector alter the effective running and are not treated); no number from this file should be compared to literature values without reconciling the scheme. This is why `README` and the report label it as scheme-specific.
