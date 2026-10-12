# Interacting scalar field

*Status: analytic derivations; symbolic check of $H_{\rm int}$ and dimensions (`interactions.py`).*

$$S_{\rm int}=-\int d^4x\sqrt{-g}\,\frac{\lambda}{4!}\phi^4.$$

**Dimension.** $[d^4x]=-4$, $[\phi]=1$ $\Rightarrow[\lambda]=0$: $\lambda\phi^4$ is classically marginal in $d=4$. (Quantum mechanically, $\lambda$ runs and $\phi^4$ needs counterterms; see `renormalization.md`.)

**Interaction Hamiltonian.** The interaction contains no derivatives so $\pi_\phi=a^2\phi'$ is unchanged by it, and the Legendre transform gives $\mathcal H_{\rm int}=-\mathcal L_{\rm int}=a^4\frac\lambda{4!}\phi^4$ (conformal time, density per unit comoving volume). Verified symbolically (`interaction_hamiltonian_density`). In terms of $\chi=a\phi$:
$$\mathcal H_{\rm int}=\frac{\lambda}{4!}\chi^4,$$
explicitly time-independent for *any* $m,\xi$: all time dependence of the interaction sits in the mode functions $v_k$.

**Free vs interacting evolution.** Write $\hat\chi_H(\eta)=U^\dagger(\eta,\eta_0)\hat\chi_I(\eta)U(\eta,\eta_0)$, $U=T\exp[-i\int_{\eta_0}^\eta H_I]$, with $\hat\chi_I$ the free field built from $v_k$. Perturbation theory is an expansion in $\lambda$ **around a specified initial state** $|\Omega\rangle$ (here the free BD vacuum, switched on adiabatically in the infinite past). Different initial states or $\eta_0$ give different perturbative series; the "interacting vacuum" is not unique and may not exist on the whole patch (IR/secular issues).

**Conformal special case.** $m=0,\xi=1/6$: $\nu=1/2$, $v=e^{-ik\eta}/\sqrt{2k}$, the problem is flat-space $\chi^4$ theory on $\eta\in(-\infty,0)$; one expects (and we find) no secular growth.
