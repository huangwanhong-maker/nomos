# 0017. Population Games and Evolutionary Dynamics

- **Status:** Accepted, provisionally (2026-11-23).
- **Date:** 2026-11-23
- **Decided by:** the infrastructure's survey *Evolutionary economics models and methods* (its §7.2, item 5, and §7.5, step 1.2), for the economics of evolution the user asked for (the infrastructure's ADR 0165).

## 1. Context

The economics of evolution needs selection as a dynamic, not only as an accounting (tyche's `census`): market shares that grow with fitness, strategies revised by imitation or comparison, conventions that settle. Evolutionary game theory gives these dynamics a common form, a revision protocol and the mean dynamic it defines, with incentive properties and Lyapunov functions that are theorems (Sandholm 2015). They name agents and choices, so they are nomos's, beside its games and agents.

## 2. Decision

`games/include/nomos/games/dynamics.hpp`:

- population games, linear (F(x) = Ax) or any callable, and Rock–Paper–Scissors in its three forms;
- the protocols of Sandholm's Table 13.1 (replicator, best response, logit, BNN, Smith) as matrices, the mean dynamic of eq. 13.3, each dynamic's closed-form field, and the Maynard Smith replicator;
- trajectories by the classical Runge–Kutta method, the incentive properties (positive correlation, Nash stationarity) at a point, the Lyapunov functions of Table 13.4 for contractive games;
- the finite population's Markov jump process, simulated exactly.

## 3. Alternatives Considered

- **In `agents`.** The finite population is an agent-based model, but its protocols and mean dynamics are game theory and are checked as such; `agents`' schedules and logs remain for models with more structure.
- **calx's Runge–Kutta.** It takes states of fixed size; the populations here have any number of strategies, and a step of four evaluations is written inline.
- **Best response as an exact inclusion.** Its trajectory takes the uniform law on ties, one selection of the inclusion; near the equilibrium the step chatters within its size, which the checks allow.

## 4. Consequences

- `check_nomos_dynamics`, 22 checks: each protocol's mean dynamic equal to its closed form on 200 random games; the replicator's closed orbits in standard Rock–Paper–Scissors (x₁x₂x₃ kept to 10⁻¹³ over 30 time units); convergence and the Lyapunov functions in good Rock–Paper–Scissors; the approach to the heteroclinic cycle in bad Rock–Paper–Scissors; positive correlation on 300 random games; Nash stationarity failing for the replicator at a vertex; the potential never falling in a coordination game, with BNN's algebraic approach to a strict equilibrium (x₁ t → 1/2); Hawk–Dove's V/C; a strictly dominated twin surviving under the Smith dynamic above a tenth of the population (Hofbauer and Sandholm 2011); the logit's rest point; the finite population of 20,000 within 0.005 of the mean dynamic (Kurtz).
- Not built: several populations, the projection and imitative logit dynamics, exact stationary distributions and stochastic stability (next), learning rules (fictitious play, reinforcement, experience-weighted attraction).
