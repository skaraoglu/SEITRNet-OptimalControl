<div align="center">

# SEITR \& Optimal Control

### Optimal control of SEITR epidemic dynamics on contact networks

*Optimal control extension on [SEITRNet](https://github.com/skaraoglu/SEITRNet). A companion implementing Pontryagin-based optimal control of the deterministic SEITR system and a numerically-optimized treatment policy for its stochastic network analogue, with cross-formulation comparison.*

[![R](https://img.shields.io/badge/R-%E2%89%A54.2-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-Academic_Use-lightgrey)]()
[![Status](https://img.shields.io/badge/Status-Pilot_Complete-brightgreen)]()
[![Topologies](https://img.shields.io/badge/Topologies-ER%20%C2%B7%20BA%20%C2%B7%20WS%20%C2%B7%20LN%20%C2%B7%20RR-orange)]()
[![DOI](https://img.shields.io/badge/DOI-EPJ-blue)](https://doi.org/10.1140/epjp/s13360-025-06481-z)

[Overview](#overview) ·
[Formulation](#mathematical-formulation) ·
[Methodology](#methodology) ·
[Repository](#repository-contents) ·
[Run](#reproducing-the-experiments) ·
[Results](#experimental-results) ·
[Cite](#citation)

</div>

*[SEITRNet++](https://github.com/skaraoglu/SEITRNet-Cpp) presents the C++ refactored library with up to x200 performance improvements.

---

## Overview

This repository accompanies a study of **optimal treatment policy** for an SEITR-class epidemic, formulated and solved at two levels:

- **Continuous level.** A deterministic SEITR ODE system with a time-varying treatment-uptake control $u(t) \in [0, \zeta]$, solved via a forward-backward sweep that combines RK4 integration of state and adjoint equations with the Pontryagin characterization of the optimal control.
- **Discrete level.** The same epidemic process simulated stochastically on an explicit contact graph with a *neighbor-resolved* force of infection. The control profile is parameterized as a piecewise-constant function over $K$ intervals and optimized numerically against the Monte Carlo objective.

The two solutions are compared head-to-head: the ODE-derived $u^\star(t)$ is evaluated on the network, and a network-fitted $u$ is evaluated against the ODE benchmark. The intent is to quantify the gap between mean-field and network-resolved optimal interventions and to expose the effect of contact-graph structure on policy performance.

This repository is the experimental complement to a baseline SEITR-on-networks pipeline; it is a script-and-notebook research artifact rather than a packaged R library.

---

## Mathematical Formulation

### Controlled deterministic system

Let $S(t), E(t), I(t), T(t), R(t)$ denote compartment sizes with $N = S + E + I + T + R$, and let $u: [0, t_f] \to [0, \zeta]$ be a measurable control. The controlled SEITR system is

$$
\begin{aligned}
\dot S &= \Lambda - \frac{\beta_1 S I}{N} - \mu S, \\
\dot E &= \frac{\beta_1 S I}{N} - (\beta_2 + \mu)E, \\
\dot I &= \beta_2 E - \bigl(\beta_3 + \mu + \delta_I + u(t)\bigr) I, \\
\dot T &= u(t)I - (\mu + \delta_T + \alpha_2)T, \\
\dot R &= \beta_3 I + \alpha_2 T - \mu R, \\
\dot N &= \Lambda - \mu N - \delta_I I - \delta_T T.
\end{aligned}
$$

In this controlled formulation $u(t)$ takes the place of the constant treatment rate $\alpha_1$ in the I-equation; the uncontrolled baseline is recovered by fixing $u \equiv \alpha_1$.

### Objective functional

$$
\mathcal{J}(u) = \int_{0}^{t_f}\Bigl[E(t) + I(t) + \tfrac{1}{2} w_1 u(t)^2\Bigr]dt,
$$

with weight $w_1 > 0$ penalizing the squared control effort. The problem is to find $u^\star \in \mathcal{U}_\text{ad} = \{u \text{ measurable} : 0 \le u(t) \le \zeta\}$ that minimizes $\mathcal{J}$ subject to the controlled dynamics above.

### Pontryagin characterization

Introducing adjoints $\boldsymbol{\lambda} = (\lambda_1,\ldots,\lambda_5)$ for $(S, E, I, T, R)$, the Hamiltonian is

$$
\mathcal{H} = E + I + \tfrac{1}{2}w_1 u^2 + \sum_{k=1}^{5} \lambda_k f_k(\mathbf{x}, u),
$$

where $f_k$ are the right-hand sides of the state equations. The adjoint system $\dot\lambda_k = -\partial\mathcal{H}/\partial x_k$ with terminal condition $\lambda_k(t_f) = 0$ for $k = 1,\ldots,5$ is integrated backward in time. Stationarity of $\mathcal{H}$ in $u$ together with the box constraint yields the closed-form characterization

$$
u^\star(t) = \min\!\Bigl(\zeta,\ \max\!\Bigl(0,\ \tfrac{I(t)\bigl(\lambda_3(t) - \lambda_4(t)\bigr)}{w_1}\Bigr)\Bigr).
$$

The forward-backward sweep iterates state forward, adjoint backward, control update via the expression above, and a convex-combination smoothing $u \leftarrow \tfrac{1}{2}(u^\star + u^\text{old})$ until a relative-change criterion is met across the eleven trajectories.

### Stochastic network analogue

The deterministic system is shadowed by a discrete-time stochastic process on a contact graph $G_t = (V_t, E_t)$. At each step:

- $S \to E$ for node $i$ with probability $\beta_1 \cdot |\mathcal{N}_i \cap I_t| / |V_t|$, where $\mathcal{N}_i$ is the neighborhood of $i$. *Force of infection is local, not mean-field.*
- $E \to I$ with probability $\beta_2$.
- $I \to R$ with probability $\beta_3$, else $I \to T$ with probability $\alpha_1 + u(t)$.
- $T \to R$ with probability $\alpha_2$.
- Disease deaths from $I$ and $T$ are sampled as $\mathrm{Binomial}(|I_t|, \delta_I)$ and $\mathrm{Binomial}(|T_t|, \delta_T)$; natural deaths as $\mathrm{Binomial}(|V_t|, \mu)$; recruitment as $\mathrm{Poisson}(\Lambda)$, with newborn nodes attached according to the host topology's generative rule.

### Basic reproduction number

For the *uncontrolled* system ($u \equiv 0$), the next-generation matrix construction yields

$$
\mathcal{R}_0 = \frac{\beta_1\beta_2}{(\beta_2 + \mu)(\beta_3 + \mu + \delta_I + \alpha_1)},
$$

reported alongside each parameter set in the experiments.

---

## Methodology

### ODE optimal control

| Step | Method |
|------|--------|
| State integration | Classical fourth-order Runge–Kutta, forward in time |
| Adjoint integration | RK4, backward in time, with terminal $\lambda_k(t_f) = 0$ |
| Control update | Pontryagin characterization, projected onto $[0, \zeta]$ |
| Smoothing | Convex combination $u \leftarrow \tfrac{1}{2}(u^\star + u^\text{old})$ |
| Cost integration | Composite Simpson's rule |

### Network optimal control

| Step | Method |
|------|--------|
| Control parameterization | Piecewise-constant on $K$ equal intervals over $[0, t_f]$ |
| Inner objective | Mean total cost over `num_exp` Monte Carlo replicates |
| Outer optimizer | L-BFGS-B with box constraints $[0, \zeta]$ via `optimParallel` |
| Parallel RNG | `L'Ecuyer-CMRG` streams seeded for reproducibility |

The ODE-derived control $u^\star(t)$ is rasterized to the network's time grid and evaluated against the Monte Carlo network process. Symmetrically, the network-fitted control is interpolated onto the fine ODE grid and evaluated. Min-max bands across replicates are reported alongside ensemble means.

---

## Repository Contents

```
.
├── README.md
├── output                           # Folder: contains experiment results.
├── SEITR_OC.R                       # Core: network simulator + ODE OC solver + utilities
├── SEITR_function.R                 # Legacy baseline SEITR simulator (no control)
├── OptControl_SEITR.ipynb           # Experiment notebook (R kernel)
└── supplement.ipynb                 # Supplement experiment and plots
```

---

## Installation

This repository ships R source and a Jupyter notebook with an R kernel. There is no installable package; functions are sourced directly.

### Required R packages

```r
install.packages(c(
  "igraph", "deSolve", "ggplot2", "dplyr", "tidyr", "gridExtra",
  "GA", "optimx", "optimParallel", "doParallel"
))
```

### R kernel for Jupyter

```r
install.packages("IRkernel")
IRkernel::installspec(user = TRUE)
```

R $\geq$ 4.2 is recommended.

---

## Reproducing the Experiments

### Quick start

```r
source("SEITR_OC.R")

# Pick a parameter set
Lambda  <- 0.4
beta1   <- 0.9;  beta2 <- 0.059; beta3 <- 0.2
alpha1  <- 0.03; alpha2 <- 0.055
delta_I <- 0.03; delta_T <- 0.03; mu <- 0.02

# ODE optimal control
oc <- solve_seitr_optimal_control(
  Lambda, beta1, beta2, beta3, alpha1, alpha2, delta_I, delta_T, mu,
  h = 0.01, t = seq(0, 100, by = 0.01),
  w1 = 0.2, zeta = 1, delta = 1e-3,
  S0 = 50, E0 = 25, I0 = 15, T0 = 5, R0 = 5
)

# Apply ODE-derived u1* to a stochastic ER network
res <- run_seitr_network_simulation(
  network_type = "ER", n = 100, n_par1 = 0.2, n_par2 = NA,
  Lambda, beta1, beta2, beta3, alpha1, alpha2, delta_I, delta_T, mu,
  init_S = 50, init_E = 25, init_I = 15, init_T = 5, init_R = 5,
  t_max = 100, u1_profile = oc$u1, num_exp = 20
)
J <- calculate_objective_functional(res$avg)$J_total
```

The notebook is structured as:

1. **Library and parameter setup.**
2. **ODE optimal control at $n = 100$** — full forward-backward sweep, objective trace, and before/after compartment trajectories.
3. **Network optimal control on Erdős–Rényi graphs** — separate sections for $p = 0.2$ and $p = 0.5$, each comprising:
   - Network optimization via L-BFGS-B in parallel,
   - Network simulation with the optimized $u$,
   - Network simulation without control,
   - Network–ODE overlays with min/max bands.

### Reproducibility

Random number generation in parallel scenarios uses `RNGkind("L'Ecuyer-CMRG")` with `clusterSetRNGStream(cl, 12345)`.

---

## Experimental Results

For each experiment the notebook reports: the optimized $u$ profile, the objective-functional trace across iterations, ensemble-mean compartment trajectories with min/max bands across the Monte Carlo replicates, and an overlay against the matched ODE solution with and without control. Numerical values of $\mathcal{J}$, $\mathcal{J}_E$, $\mathcal{J}_I$, $\mathcal{J}_W$ are tabulated for each scenario.

> **Reproducibility caveat.** Because the notebook runs Monte Carlo network simulations, exact figure replication requires the same seed *and* the same parallel-stream configuration. Numerical values may shift at the second-decimal level across hardware.

---


## Citation

If this software contributes to a publication, please cite: Karaoglu S, Imran M, McKinney BA. Network-based SEITR epidemiological model with contact heterogeneity: comparison with homogeneous models for random, scale-free and small-world networks. The European Physical Journal Plus. 2025 Jun 18;140(6):551, https://doi.org/10.1140/epjp/s13360-025-06481-z.

---

## License

Released for academic use. See `LICENSE` for full terms. For commercial use, contact the maintainer.

---

<div align="center">

<sub>Built with <a href="https://igraph.org/r/">igraph</a>, <a href="https://desolve.r-forge.r-project.org/">deSolve</a>, and <a href="https://cran.r-project.org/package=optimParallel">optimParallel</a>

</div>
