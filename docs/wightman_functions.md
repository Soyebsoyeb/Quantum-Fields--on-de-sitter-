# Wightman function

*Status: analytic derivation; numerical mode sum vs closed form (Fig. 5); short-distance expansion (Fig. 6b).*

## Definition

$G^+(x,x')=\langle\Omega|\hat\phi(x)\hat\phi(x')|\Omega\rangle$ with $|\Omega\rangle$ the BD state. From the expansion and $\langle\hat a_{\mathbf k}\hat a^\dagger_{\mathbf k'}\rangle=(2\pi)^3\delta(\mathbf k-\mathbf k')$:
$$G^+(\eta,\eta',\mathbf r)=\int\frac{d^3k}{(2\pi)^3}\phi_k(\eta)\phi_k^*(\eta')e^{i\mathbf k\cdot\mathbf r}=\frac1{2\pi^2}\int_0^\infty dk\,k^2\phi_k(\eta)\phi_k^*(\eta')\frac{\sin kr}{kr}.$$
Measure $d^3k/(2\pi)^3$; $\phi_k=v_k/a$. $G^+$ is a bi-**distribution**: the integrand decays only as $k^{-1}$ at equal times, so the integral is conditionally (Abel) summable for $r>0$ and singular at $r\to0$.

## Regularisation used

Subtract and add the conformal-vacuum piece $|v_k|^2\to1/2k$, whose Abel-regularised transform is exactly $1/(4\pi^2r^2)$ ($\int_0^\infty\sin(kr)\,dk=1/r$ in the Abel sense):
$$G^+(\eta,r)=H^2\eta^2\Big[\frac1{4\pi^2r^2}+\frac1{2\pi^2r}\int_0^\infty dk\,k\Big(|v_k|^2-\frac1{2k}\Big)\sin kr\Big].$$
The remainder decays as $k^{-2}$ and is computed by QUADPACK QAWF. This is an *explicitly stated distributional prescription*, not a momentum cutoff.

## Closed form and checks

$$G^+=\frac{H^2}{16\pi^2}\Gamma(\tfrac32+\nu)\Gamma(\tfrac32-\nu)\,{}_2F_1\!\big(\tfrac32+\nu,\tfrac32-\nu;2;\tfrac{1+Z}2\big),\qquad Z=\frac{\eta^2+\eta'^2-r^2}{2\eta\eta'}$$
(Chernikov–Tagirov; Bunch–Davies). Verified (i) as a solution of the invariant ODE $(1-Z^2)G''-4ZG'-(\tfrac94-\nu^2)G=0$ to $10^{-30}$ relative; (ii) $\nu=1/2$ reduces to $H^2\eta\eta'/[4\pi^2(r^2-\Delta\eta^2)]$; (iii) agrees with the mode sum to $\lesssim10^{-13}$ ($\nu$ real) and $\lesssim10^{-9}$ ($\nu=0.8i$, `mpmath` Hankel) — a test of the Hankel normalisation, the Fourier convention and $Z$ simultaneously.

## Dependence on parameters

Through $\nu=\nu(m^2/H^2+12\xi)$ only (and the state). $\nu$ real: large-separation power law $G\sim(r/|\eta|)^{-(3-2\nu)}$ (Fig. 5; slower decay as $\nu\to3/2$); imaginary $\nu=i\mu$: $G\sim(r/|\eta|)^{-3}$ with log-periodic oscillations $(r/|\eta|)^{\pm2i\mu}$. The state enters through $v_k$; changing to an $\alpha$-vacuum would add $e^{\pm}$ Bogoliubov terms and, away from BD, antipodal singularities.
