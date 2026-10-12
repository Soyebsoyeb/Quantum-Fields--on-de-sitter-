# Propagators and Green functions

*Status: analytic definitions; mode-space identities verified numerically (`propagators.py`, Fig. 8). Position-space Feynman function: defined, not implemented.*

For the operator $D_k=\partial_\eta^2+w_k^2(\eta)$ acting on $\chi_k$ and a normalised BD mode $v_k$:

| Function | Definition | Equation | Boundary condition |
|---|---|---|---|
| Wightman | $G^+=v(\eta)v^*(\eta')$ | $D G^+=0$ | state-dependent (BD) |
| Feynman | $G^F=\theta(\eta-\eta')v(\eta)v^*(\eta')+\theta(\eta'-\eta)v^*(\eta)v(\eta')$ | $DG^F=-i\delta(\eta-\eta')$ | BD positive-frequency (the $i\epsilon$ at $\eta\to-\infty$) |
| Retarded | $G^R=i\theta(\eta-\eta')[v(\eta)v^*(\eta')-v^*(\eta)v(\eta')]$ | $DG^R=\delta(\eta-\eta')$ | $G^R=0$ for $\eta<\eta'$ — **state independent** |
| Advanced | $G^A(\eta,\eta')=G^R(\eta',\eta)$ | $DG^A=\delta$ | $G^A=0$ for $\eta>\eta'$ |

**Delta strengths.** The derivative jump across $\eta=\eta'$ is, for $G^R$: $i(v'v^*-v^{*\prime}v)=i(-i)=1$; for $G^F$: $v'v^*-v^{*\prime}v=-i$ (uses $W=i$). Hence the signs above.

**State independence.** $G^R=-2\theta\,\mathrm{Im}[v(\eta)v^*(\eta')]$ is built from the commutator $\langle[\hat\chi(\eta),\hat\chi(\eta')]\rangle=v v^{*\prime}$-antisymmetrisation, a c-number fixed by the CCR; every normalised mode basis (every Fock vacuum) gives the same $G^R$. $G^F$ and $G^+$ do depend on the state.

## Numerical verification (mode space, $\nu=3/2$, $\eta'=-6$)

1. $G^R$ from Hankel/closed-form modes vs. an *independent* numerical solution of $D y=0$ with $y(\eta')=0,y'(\eta')=1$ (no Hankel functions): max difference $4\times10^{-12}$.
2. $\partial_\eta G^R|_{\eta'^+}=1.0000000001$.
3. $G^F$ is continuous at $\eta=\eta'$ with a kink of slope $-i$; its real part equals $\mathrm{Re}\,G^+$.

## Position space

The position-space Feynman function is the distribution $G^F=\theta G^++\theta' G^-$; its singularity at $x=x'$ is the same as in flat space: $\frac{1}{4\pi^2}\frac1{\sigma^2+i\epsilon}$-type, so **coincident limits are not finite numbers** (see `ultraviolet_structure.md`). We therefore never evaluate $G(x,x)$ without a stated subtraction.
