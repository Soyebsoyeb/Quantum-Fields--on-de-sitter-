# Tensor perturbations and the boundary with quantum gravity

*Status: tensor equation = symbolic verification; action normalisation and spectrum = quoted standard result + numerical check; everything on quantisation of gravity = conceptual discussion.*

## Decomposition

$g_{\mu\nu}=\bar g_{\mu\nu}+h_{\mu\nu}$, $\bar g$ = de Sitter. Gauge freedom $h_{\mu\nu}\to h_{\mu\nu}+\bar\nabla_\mu\xi_\nu+\bar\nabla_\nu\xi_\mu$ (linearised diffeomorphisms). On the flat chart, decompose by $SO(3)$ helicity: scalars (2 metric functions in a gauge-invariant combination plus constraints), vectors (decaying, no propagating dof in pure gravity), tensors $h_{ij}^{TT}$ (2 propagating polarisations). Tensors are gauge invariant at linear order and decouple from scalars/vectors at linear order.

## Tensor equation (derived)

With $ds^2=a^2[-d\eta^2+(\delta_{ij}+h_{ij})dx^idx^j]$, $R_{\mu\nu}=\Lambda g_{\mu\nu}$, $\Lambda=3H^2$, and $h_{xx}=-h_{yy}=\epsilon h(\eta,z)$, the SymPy computation shows the background residual is zero and
$$\delta(R_{xx}-\Lambda g_{xx})=\tfrac12\epsilon\big[h''+2\tfrac{a'}{a}h'-\partial_z^2h\big]\ \Rightarrow\ h_k''-\tfrac2\eta h_k'+k^2h_k=0.$$
(`tensor_modes.derive_tensor_equation`; the first-order part is *derived*, not recalled.)

## Canonical normalisation (quoted)

$S_2=\frac{M_p^2}8\int d\eta\,d^3x\,a^2[(h'_{ij})^2-(\partial_kh_{ij})^2]$, $M_p^2=1/8\pi G$. With $h_{ij}=\sum_\lambda e^\lambda_{ij}h_\lambda$, $e^\lambda_{ij}e^{\lambda'}_{ij}=2\delta^{\lambda\lambda'}$, $v_\lambda=M_pa\,h_\lambda/\sqrt2$ obeys $v''+(k^2-a''/a)v=0$: each polarisation is a **massless minimally coupled scalar** ($\nu=3/2$, same $v_k$ as in `mode_equation.md`). Thus
$$\Delta_t^2=\frac{k^3}{2\pi^2}\langle h_{ij}h_{ij}\rangle'=\frac{2H^2}{\pi^2M_p^2}.$$
Numerically integrating the ODE from adiabatic data reproduces this plateau and the full $x$-dependence to $<10^{-6}$ (Fig. 12). The scale invariance at late times is the graviton counterpart of the frozen scalar spectrum. (Grishchuk 1974; Starobinsky 1979; Mukhanov–Feldman–Brandenberger 1992.)

## What changes when the metric is quantised

* A test scalar on fixed $\bar g$ has no back-reaction: $\delta g_{\mu\nu}=0$. Quantising $h_{\mu\nu}$ promotes the geometry itself to operators; the scalar's $T_{\mu\nu}$ sources $h$ ($\Box h\sim8\pi G\,T$), giving $O(G)$ corrections that are cubic/quartic gravitational and matter-gravity vertices.
* Gauge: quantisation requires gauge fixing (e.g., de Donder) with Faddeev–Popov ghosts, or reduction to gauge-invariant variables; the Hamiltonian and momentum constraints remove the unphysical components. In dS there is no spatial infinity/time-translation Killing energy: "energy" and the S-matrix have no standard meaning (Abbott–Deser 1982 define charges perturbatively).
* Renormalisation: GR is perturbatively non-renormalisable ('t Hooft–Veltman 1974; Goroff–Sagnotti 1986 two-loop divergence). It is a valid *effective field theory* below $M_p$ (Donoghue 1994) with controlled corrections $\sim H^2/M_p^2$, but with dS-specific IR issues (secular growth of graviton loops, Tsamis–Woodard line of work; contested).

## Not claimed

Quantising $\phi$ does not quantise gravity. We solve free tensor modes on a fixed background only; no graviton loops, no non-linear constraint solution, no nonperturbative de Sitter quantum gravity. The boundary of the fixed-background calculation is exactly $O(H^2/M_p^2)$ corrections, back-reaction and global issues (finiteness of the dS entropy, observables in a closed universe).
