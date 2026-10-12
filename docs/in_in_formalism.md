# In-in (Schwinger–Keldysh) formalism

*Status: derivation; implemented at first order for contact diagrams (`in_in.py`). Higher orders: unimplemented.*

## Why not in–out

An in–out amplitude $\langle\text{out}|S|\text{in}\rangle$ requires asymptotic regions with well-defined particle states and a vacuum-persistence structure. In an expanding patch there is no late-time asymptotic region with a free-particle interpretation ($\eta\to0^-$ is a spacelike boundary), and what is measured are expectation values $\langle\Omega|Q(\eta)|\Omega\rangle$ at a *finite* time in a *given* initial state, not transition amplitudes. The in-in expectation value uses $U^\dagger QU$ — both the forward and backward evolution — so the final state is never specified.

## Formula

$$\langle Q(\eta)\rangle=\Big\langle\Omega\Big|\,\bar T e^{\,i\int_{-\infty}^\eta H_I}\;Q_I(\eta)\;Te^{-i\int_{-\infty}^\eta H_I}\Big|\Omega\Big\rangle.$$
Closed-time path: the contour runs $-\infty\to\eta$ (time-ordered branch "+", vertices $-iH_I$) then back $\eta\to-\infty$ (anti-time-ordered "−", vertices $+iH_I$); the vertical segment closing the contour at $-\infty$ is controlled by the $i\epsilon$ (Abel) prescription $\eta'\to-\infty(1\mp i\epsilon)$ on the $\pm$ branches so that the interacting vacuum is projected out of the free one.

**Assumptions.** (i) $|\Omega\rangle$ is the free BD vacuum; (ii) coupling switched on adiabatically; (iii) $H_I$ has no derivative couplings (so $H_I=-L_I$); (iv) first order in $\lambda$ — higher orders require time-ordered nested integrals and renormalisation.

## First order

$$\langle Q\rangle^{(1)}=i\int_{-\infty}^\eta d\eta'\langle[H_I(\eta'),Q(\eta)]\rangle=-2\,\mathrm{Im}\int^\eta d\eta'\langle H_I(\eta')Q(\eta)\rangle.$$
The commutator ("retarded") form manifestly has support only on $\eta'<\eta$, and each term is finite when $\mathrm{Im}\to$ cancels the real divergent pieces; this is why the $\mathrm{Im}$ is computed *before* integrating in the extended-precision backend (cancellation $\sim(k|\eta|)^{-7}$ at late times in double precision).

## Equivalence of the two forms

`schwinger_keldysh_branches` returns the "+" branch contribution $-i\int\langle QH\rangle=-iX^*$ and the "−" branch $+i\int\langle HQ\rangle=iX$; their sum is $-2\,\mathrm{Im}X$, identical to the commutator form (tested).

## Independent validation

The 0+1-dimensional analogue $H=p^2/2+\omega^2q^2/2+gq^4/4!$ has the same formula, $\kappa_4^{(1)}=-g/32\omega^5$, and the exact ground state obtained by diagonalising $H$ in a truncated oscillator basis reproduces the slope to $1.6\times10^{-6}$ relative (Fig. 11). Sign and normalisation of the in-in expression are therefore not merely derived but cross-checked against a non-perturbative computation.
