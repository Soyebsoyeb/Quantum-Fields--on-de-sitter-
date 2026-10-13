# de Sitter geometry

*Status: analytic result + symbolic verification (`backgrounds.py`, `tests/test_curvature.py`).*

## Embedding

$\mathrm{dS}_4$ is the hyperboloid $-X_0^2+X_1^2+X_2^2+X_3^2+X_4^2=H^{-2}$ in $\mathbb R^{1,4}$ with the induced metric. It is maximally symmetric ($SO(1,4)$), $R_{\mu\nu\rho\sigma}=H^2(g_{\mu\rho}g_{\nu\sigma}-g_{\mu\sigma}g_{\nu\rho})$, hence $R_{\mu\nu}=3H^2g_{\mu\nu}$ and $R=12H^2$.

The invariant $Z(x,x')=H^2X\!\cdot\!X'$ equals $\cos(H\sigma)$ for spacelike geodesic separation $\sigma$, $Z=1$ at coincidence, $Z=-1$ at the antipodal point (and light-cone of the antipode), and $Z>1$ for timelike separation.

## Three charts that must not be confused

| Chart | Coordinates | Covers | Used here |
|---|---|---|---|
| Global | $ds^2=-dT^2+H^{-2}\cosh^2(HT)d\Omega_3^2$ | whole manifold | not used; discussed in `vacuum_states.md` |
| Spatially flat (expanding) | $ds^2=-dt^2+e^{2Ht}d\mathbf x^2$, $\eta<0$ | the half $X_0+X_4>0$ | **all numerics** |
| Static | $ds^2=-(1-H^2r^2)dt^2+\frac{dr^2}{1-H^2r^2}+r^2d\Omega^2$ | one observer's causal diamond, $r<1/H$ | curvature check only |

The flat chart is obtained by $X_0=\frac1H\sinh Ht+\frac H2e^{Ht}x^2$, $X_4=\frac1H\cosh Ht-\frac H2e^{Ht}x^2$, $X_i=e^{Ht}x_i$ (one checks $X_4^2-X_0^2=H^{-2}-e^{2Ht}x^2$). It is **geodesically incomplete to the past**: $\eta\to-\infty$ is a coordinate horizon, not a past boundary of the manifold. A "state defined at $\eta\to-\infty$" is therefore a statement about the behaviour of modes near the boundary of the patch, not an initial condition on a Cauchy surface of the full spacetime.

In the flat chart $Z=\dfrac{\eta^2+\eta'^2-r^2}{2\eta\eta'}$, $r=\vert\mathbf x-\mathbf x'\vert$ (this expression is *tested*: the independent mode-sum and the closed form agree to $10^{-13}$ only if it is correct, `tests/test_two_point.py`).

## Horizon crossing

A mode with comoving wavenumber $k$ has physical wavenumber $k/a$; $k=aH\iff k\vert\eta\vert=1$. The comoving Hubble radius is $1/(aH)=-\eta$, shrinking in time, so every mode eventually leaves the horizon (Fig. 1).
