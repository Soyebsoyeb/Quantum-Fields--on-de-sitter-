# Klein–Gordon equation on the flat chart

*Status: analytic derivation, symbolic verification (`cosmic_time_equation_difference`, `conformal_time_equation`, `cosmic_to_conformal_chain_rule`).*

## Cosmic time

$\sqrt{-g}=a^3$, $g^{tt}=-1$, $g^{ij}=a^{-2}\delta^{ij}$:
$$\Box\phi=\frac1{a^3}\partial_t\big(-a^3\dot\phi\big)+\frac{\nabla^2\phi}{a^2}=-\ddot\phi-3H\dot\phi+\frac{\nabla^2\phi}{a^2}.$$
Therefore
$$\ddot\phi+3H\dot\phi-\frac{\nabla^2}{a^2}\phi+(m^2+\xi R)\phi=0,\qquad R=12H^2.$$

**Origin of $3H\dot\phi$.** It is $\partial_t\ln\sqrt{-g}\,\dot\phi$: the physical volume of a comoving region grows as $a^3$, so the field momentum density $a^3\dot\phi$ is conserved only in the absence of the gradient and mass terms; equivalently, the friction term dilutes field gradients and, for a light field, freezes $\phi$ once $k/a\ll H$.

## Conformal time

With $\sqrt{-g}=a^4$, $g^{\eta\eta}=-a^{-2}$:
$$\Box\phi=\frac1{a^4}\Big[-\partial_\eta(a^2\phi')+a^2\nabla^2\phi\Big]=-\frac{\phi''}{a^2}-\frac{2a'}{a^3}\phi'+\frac{\nabla^2\phi}{a^2}.$$
Multiplying $(\Box-m^2-\xi R)\phi=0$ by $-a^2$:
$$\phi''+2\frac{a'}{a}\phi'-\nabla^2\phi+a^2(m^2+\xi R)\phi=0,\qquad R=\frac{6a''}{a^3}.$$
This was derived from the covariant $\Box$, and **independently** confirmed by the chain rule $\partial_t=a^{-1}\partial_\eta$ applied to the cosmic-time equation (both are asserted equal in `tests/test_mode_equations.py`).
