# Limitations

Honest scope statement; every item is either a stated assumption or an unimplemented extension.

**Scope.** Fixed four-dimensional de Sitter, spatially flat chart only. Free scalar: fully implemented and validated for $\nu$ real and imaginary. Interactions: first-order, tree-level contact diagram (the conformal case analytically, the massless minimal case with its log). No loops, no second order, no exchange diagrams.

**Not implemented / only described**: $\alpha$-vacua; unequal-time position-space Wightman via mode sum (the closed form is available); position-space Feynman propagator; the global or static chart; Hadamard recursion and a covariant renormalised $\langle\phi^2\rangle$ or $\langle T_{\mu\nu}\rangle$; loop corrections, running of $\lambda$, stochastic resummation beyond the equilibrium benchmark; scalar and vector metric perturbations, graviton interactions, ghosts.

**Numerical.** Double precision Hankel functions at late times lose $\sim(k\vert\eta\vert)^{-3}$ digits in the commutator combination; `mpmath` is used where this matters and `scipy` results are only trusted away from $\eta\to0$ (a dedicated test cross-checks both). The ODE Wronskian residual grows with $\vert v\vert\vert v'\vert$. Figures show numerics, which validate (they do not prove) the stated analytic results.

**Physics caveats.** (i) The BD state is a statement about the chart boundary, not an initial condition on a complete Cauchy surface. (ii) For $\nu\ge3/2$ the BD two-point function is IR divergent; results with a cutoff $k_{\rm IR}$ depend on it. (iii) Point-splitting constants are scheme-dependent. (iv) The first-order secular growth signals breakdown of the *fixed-order expansion*, not an inconsistency of the theory. (v) The tensor action normalisation is quoted from the literature (standard) and validated by reproducing $\Delta_t^2=2H^2/\pi^2M_p^2$ numerically, not rederived symbolically from the Einstein–Hilbert action.

**References.** `references.bib` was typed from standard bibliographic knowledge; there was no network access to verify DOIs/pages in this session. Check each entry (especially volume/page numbers) before citing in a publication.
