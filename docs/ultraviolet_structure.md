# Ultraviolet structure

*Status: analytic expansion verified numerically (Fig. 6).*

## Short-distance form

With $u=(1-Z)/2=\sin^2(H\sigma/2)$ ($\sigma$ geodesic distance), the closed-form $G^+$ has, from the ${}_2F_1$ connection formula with $c=a+b-1$ (DLMF 15.8.10),
$$G^+=\frac{A}{u}+B\ln u+C+O(u\ln u),\quad A=\frac{H^2}{16\pi^2},\ B=\big(\tfrac94-\nu^2-2\big)A,\ C=B\big[\psi(\tfrac32+\nu)+\psi(\tfrac32-\nu)+2\gamma_E-1\big].$$
The Frobenius recursion of the invariant ODE independently gives $B=(\mathrm{ab}-2)A$ — the same value. Numerically the residual $G-(A/u+B\ln u+C)$ scales like $u\ln u$ (Fig. 6b) for $\nu=1,1.3,0.8i$ (for $\nu=1/2$, $B=C=0$ and the residual is round-off).
Because $u\simeq H^2\sigma^2/4$, $A/u\to1/(4\pi^2\sigma^2)$: the **leading singularity is the flat-space one** (Fig. 6a: $4\pi^2d^2G^+\to1$, $d$ = proper distance), independent of $m,\xi$ and the (BD) state. The $\ln u$ coefficient $B$ is local geometry/mass data ($m^2+(\xi-\frac16)R$ times $H^2/16\pi^2$-type terms).

## Hadamard condition

A state is *Hadamard* if
$$G^+(x,x')=\frac{1}{8\pi^2}\Big[\frac{\Delta^{1/2}}{\sigma^2}+V\ln\sigma^2+W\Big],\ \ \sigma^2\to\sigma^2+i\epsilon\,(t-t')+\dots$$
with $\Delta^{1/2}$ the Van Vleck–Morette determinant and $V$ determined *locally* by the geometry and $m,\xi$ (via Hadamard recursion); only the smooth $W$ carries state information. Equivalent microlocal form: wavefront set condition (Radzikowski 1996). Hadamard states are those for which $\langle T_{\mu\nu}\rangle$ and Wick polynomials can be consistently defined (Wald 1977; Hollands–Wald 2001). BD satisfies it (the structure above, with $V\leftrightarrow B$), while generic $\alpha$-vacua do not (antipodal singularity). **Status:** the expansion above is checked; the full Hadamard recursion and the wavefront-set statement are *not* implemented, only described.

## Cutoffs are not renormalisation

Fig. 6c: with a sharp comoving cutoff $\Lambda$, $\langle\phi^2\rangle_\Lambda\simeq\frac{H^2}{4\pi^2}\big[\ln\frac\Lambda{k_{\rm IR}}+\frac{\eta^2\Lambda^2}{2}\big]$ ($\nu=3/2$, exact): a quadratic divergence $\propto(\Lambda/a)^2$ in *physical* units plus a log. A comoving cutoff breaks de Sitter invariance and time-translation covariance (physical cutoff $\Lambda/a$ changes with time), and the finite remainder depends on the cutoff shape. It is a diagnostic of the divergence structure only; no counterterm, renormalisation condition or scheme is attached to it.
