# Vacuum states

*Status: analytic derivation + conceptual discussion.*

## The distinct objects

| Object | What it is | Not to be confused with |
|---|---|---|
| Classical solution | a function $\phi(\eta,\mathbf x)$ obeying the KG equation | an operator |
| Normalised mode basis | solutions $v_k$ with $W=i$ | a state |
| Fock-space construction | $\hat a_k$, $\hat a_k^\dagger$ and $\vert0\rangle$ built from a basis | a unique physical vacuum |
| Quantum state | positive normalised linear functional $\omega$ on the field algebra | a coordinate chart |
| Wightman function | $\omega(\hat\phi(x)\hat\phi(x'))$, a bi-distribution | a function with a coincident limit |
| Time-ordered propagator | $\omega(T\hat\phi\hat\phi)$; also depends on the state | the Green function of $\Box$ alone |
| Physical observable | Hermitian operator, smeared (e.g. $\langle\phi^2\rangle$ after renormalisation, $\langle T_{\mu\nu}\rangle$) | a formal mode-sum |

## Bunch–Davies state

**Definition used here.** For $\nu$ such that the modes exist as above, the BD state is the Fock vacuum of $v_k=\frac{\sqrt\pi}2e^{i\pi(\nu/2+1/4)}\sqrt{-\eta}H^{(1)}_\nu(-k\eta)$, which satisfies $v_k\to e^{-ik\eta}/\sqrt{2k}$ for $k\vert\eta\vert\to\infty$: at early times each mode has physical wavelength $a/k\ll H^{-1}$ and the background is locally flat, so the state looks like the Minkowski vacuum.

**Conditions and limitations.**

1. The chart is incomplete to the past; "early time" means near the chart's horizon $\eta\to-\infty$. The BD state coincides with the Euclidean (Hartle–Hawking-like) vacuum obtained by analytic continuation from the sphere, and is invariant under the full $SO(1,4)$ (Chernikov–Tagirov 1968; Bunch–Davies 1978; Allen 1985) — for fields whose two-point function exists (below).
2. It is part of the one-parameter family of "$\alpha$-vacua" $\vert\alpha\rangle$ (Mottola 1985; Allen 1985). Only $\alpha=$BD is Hadamard at short distance with no antipodal singularity; the others generally are not Hadamard (singular at antipodal points) — see `ultraviolet_structure.md`.
3. For $\nu>3/2$ (tachyonic $m^2+12\xi H^2<0$) and for the massless minimally coupled field ($\nu=3/2$) the BD two-point function is infrared divergent: $\int dk\,k^2\vert\phi_k\vert^2$ diverges at $k\to0$. For $\nu=3/2$ no de Sitter-invariant Fock vacuum exists (Allen–Folacci 1987); one may keep the mode functions but must accept broken dS invariance (e.g., via an IR cutoff or by restricting to shift-invariant observables such as derivatives).
4. A statement "the vacuum in conformal coordinates" is meaningless without the mode basis. The code constructs states only from explicit modes, and `initial_conditions.adiabatic_data` offers a Hankel-free, local construction (adiabatic vacuum of order 0/2/4) whose convergence to BD is shown numerically.

## Static-patch interpretation (conceptual)

Restricted to a static patch the BD state is thermal with temperature $T=H/2\pi$ with respect to static time (Gibbons–Hawking 1977). That is a statement about the *restriction* of one global state, not a separate vacuum. Not implemented.
