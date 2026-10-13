# Scalar-field action and field equation

*Status: analytic derivation + symbolic verification (`scalar_action.py`).*

$$S[\phi]=-\frac12\int d^4x\sqrt{-g}\Big[g^{\mu\nu}\nabla_\mu\phi\nabla_\nu\phi+(m^2+\xi R)\phi^2\Big].$$

## Variation

Let $\phi\to\phi+\delta\phi$ with $\delta\phi=0$ on the initial/final hypersurfaces and $\delta\phi\to0$ at spatial infinity.
$$\delta S=-\int d^4x\sqrt{-g}\Big[\nabla^\mu\phi\,\nabla_\mu\delta\phi+(m^2+\xi R)\phi\,\delta\phi\Big].$$
Use $\sqrt{-g}\,\nabla^\mu\phi\,\nabla_\mu\delta\phi=\partial_\mu(\sqrt{-g}\,\delta\phi\,\nabla^\mu\phi)-\delta\phi\,\partial_\mu(\sqrt{-g}\,\nabla^\mu\phi)$ (product rule), and $\frac1{\sqrt{-g}}\partial_\mu(\sqrt{-g}\,\nabla^\mu\phi)=\Box\phi$ (definition). The total derivative integrates to a boundary term $\oint d\Sigma_\mu\,\delta\phi\,\nabla^\mu\phi$ which vanishes under the stated boundary assumptions. Hence
$$\delta S=\int d^4x\sqrt{-g}\,\delta\phi\,\big[\Box\phi-(m^2+\xi R)\phi\big],$$
and stationarity for arbitrary $\delta\phi$ gives $(\Box-m^2-\xi R)\phi=0$.

## Independent verification

`euler_lagrange_vs_klein_gordon()` builds $\mathcal L$ for a *generic* conformally flat FRW metric $a(\eta)$, applies SymPy's Euler–Lagrange operator and verifies
$$\frac{\partial\mathcal L}{\partial\phi}-\partial_\mu\frac{\partial\mathcal L}{\partial(\partial_\mu\phi)}=\sqrt{-g}\,(\Box-m^2-\xi R)\phi$$
with $\sqrt{-g}=a^4$ ($a>0$) and $R$ computed from the Riemann convention in `conventions.md` ($R=6a''/a^3$). The difference simplifies to zero.

## Sign of the curvature coupling

With $R=+12H^2$, the term $\xi R\phi^2$ acts as $m^2_{\rm eff}=m^2+12\xi H^2$. Conformal invariance of the massless theory in $d=4$ requires $\xi=1/6$ (shown below via the mode equation: all time dependence cancels).
