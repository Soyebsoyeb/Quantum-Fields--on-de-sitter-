# Infrared structure

*Status: analytic + numerical experiment (Figs. 4, 7) + comparison with an effective stochastic benchmark.*

## Small-$k$ behaviour

For real $\nu>0$, $H^{(1)}_\nu(x)\simeq-\tfrac{i}{\pi}\Gamma(\nu)(2/x)^\nu$ as $x\to0$ (DLMF 10.7.2), hence
$$|\phi_k|^2\to\frac{H^2\eta^2(-\eta)}{4\pi}\Gamma(\nu)^2\,2^{2\nu}(-k\eta)^{-2\nu}\propto k^{-2\nu}.$$
The two-point integral $\int dk\,k^2|\phi_k|^2$ has IR behaviour $\int dk\,k^{2-2\nu}$: convergent iff $\nu<3/2$, logarithmic at $\nu=3/2$, power-divergent for $\nu>3/2$.

## Massless vs light massive

* $\nu=3/2$: $\Delta^2_\phi=H^2/4\pi^2$ frozen; variance of modes between a comoving IR scale $k_{\rm IR}$ and the horizon grows linearly, $\langle\phi^2\rangle=\frac{H^2}{4\pi^2}N$, i.e. $\frac{H^3t}{4\pi^2}$ (Fig. 7; closed form checked in `test_phi2_cutoff_numeric_matches_closed_form`).
* $\nu<3/2$ ($m\ll H$): $3/2-\nu\simeq m^2/3H^2$; the growth saturates after $N\sim3H^2/2m^2$ e-folds at $\sim3H^4/8\pi^2m^2$ (dotted lines in Fig. 7; the coarse-grained numerical plateau for $m=0.6H$ is $0.096\,H^2$ vs the leading-order estimate $0.105\,H^2$, as expected because that estimate is leading in $m/H$).

## Independent benchmark

The point-split $\langle\phi^2\rangle$ (see `renormalization.md`) for $m\to0$ tends to $3H^4/8\pi^2m^2$ (via $\psi(3/2-\nu)\to-1/(3/2-\nu)$), matching the stationary solution $P_{\rm eq}\propto e^{-8\pi^2V/3H^4}$ of the Starobinsky Fokker–Planck equation for $V=\frac12m^2\phi^2$ (tested: ratio $\to1+O(m^2/H^2)$). For $V=\lambda\phi^4/4!$ the same equilibrium gives $\langle\phi^2\rangle_{\rm eq}\approx0.323\,H^2/\sqrt{\lambda}$ (numerical, `stochastic_equilibrium_moments`; non-perturbative in $\lambda$, which shows the fixed-order series cannot reproduce it).

## What is and is not claimed

An IR-divergent *intermediate* expression ($\int dk/k$ with a formal cutoff) is not a demonstrably ill-defined *physical prediction*: observables insensitive to the zero mode (derivative correlators, $\langle\phi(x)-\phi(x')\rangle^2$) are IR finite at $\nu=3/2$; the dependence on $k_{\rm IR}$ in $\langle\phi^2\rangle$ reflects the dependence on the state/initial patch size. Conversely, the absence of an $SO(1,4)$-invariant state at $\nu=3/2$ is a theorem-level statement (Allen–Folacci 1987), not an artefact of the cutoff. Interacting IR behaviour: Starobinsky 1986; Starobinsky–Yokoyama 1994; Marolf–Morrison 2011; Hollands 2013; Polyakov 2012 (not endorsed here as settled).
