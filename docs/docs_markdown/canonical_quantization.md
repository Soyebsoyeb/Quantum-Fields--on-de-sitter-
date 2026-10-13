# Canonical quantisation

*Status: analytic derivation; numerical Wronskian checker in `normalization.py`.*

## Phase space

From $S=\frac12\int[\chi'^2-(\nabla\chi)^2-w^2\!\text{-term}]$ the canonical momentum is $\pi_\chi=\partial\mathcal L/\partial\chi'=\chi'$. For the original field, $\pi_\phi=a^2\phi'$ and $[\hat\phi(\eta,\mathbf x),\hat\pi_\phi(\eta,\mathbf y)]=i\delta^{(3)}(\mathbf x-\mathbf y)$. Since $\chi=a\phi$, $\pi_\chi=\chi'=\pi_\phi/a+\tfrac{a'}{a}\chi$; the shift term commutes with $\chi$, so $[\hat\chi,\hat\pi_\chi]=[\hat\phi,\hat\pi_\phi]=i\delta^{(3)}$ and the two descriptions are canonically equivalent.
Equal-time $[\hat\phi,\hat\phi]=[\hat\pi,\hat\pi]=0$.

## Expansion and Wronskian

Expand $\hat\chi=\int\frac{d^3k}{(2\pi)^3}[\hat a_{\mathbf k}v_ke^{i\mathbf k\cdot\mathbf x}+\text{h.c.}]$ with $[\hat a,\hat a^\dagger]=(2\pi)^3\delta$. Then
$$[\hat\chi(\mathbf x),\hat\pi_\chi(\mathbf y)]=\int\frac{d^3k}{(2\pi)^3}\big(v_kv_k^{\ast\prime}-v_k^{\ast}v_k'\big)e^{i\mathbf k\cdot(\mathbf x-\mathbf y)}=i\delta^{(3)}(\mathbf x-\mathbf y)$$
iff $\boxed{v_kv_k^{\ast\prime}-v_k^{\ast}v_k'=i}$. For a real second-order ODE $v''+w^2v=0$ the Wronskian is $\eta$-independent ($W'=vv^{\ast\prime\prime}-v^{\ast}v''=-w^2(vv^{\ast}-v^{\ast}v)=0$), so checking it at one time suffices analytically; numerically it measures integration error.

## Mode basis vs state

Another basis $u_k=\alpha_kv_k+\beta_kv_k^{\ast}$ with $\vert\alpha\vert^2-\vert\beta\vert^2=1$ gives the same operator algebra and the same Wronskian; the corresponding annihilation operators are Bogoliubov transforms $\hat b=\alpha^{\ast}\hat a-\beta^{\ast}\hat a^\dagger$. Different bases define **different Fock vacua** $\vert0_v\rangle\neq\vert0_u\rangle$ (inequivalent in general), hence different states. A mode basis is data used to build a Fock space; the state is the one annihilated by the $\hat a_k$ attached to the chosen basis.

## Numerical normalisation checker

`check_normalization` returns the pointwise residual, initial/final/max residual and location, growth factor, solver settings, and a pass/fail against $10^3\,$rtol. No rescaling. Observed behaviour: with exact initial data the residual stays at the $10^{-11}$ level for moderate $\vert v\vert\vert v'\vert$; for modes that grow like $1/\vert\eta\vert$ (massless minimal, $\eta\to-10^{-3}$) the absolute residual rises to $\sim10^{-7}$ at rtol $=10^{-12}$ because the residual scales with $\vert v\vert\vert v'\vert\sim10^9$ — accuracy loss that is reported, not hidden (`tests/test_tensor.py`).
