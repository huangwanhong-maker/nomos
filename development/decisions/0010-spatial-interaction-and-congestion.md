# 0010. Spatial Interaction and Congestion

- **Status:** Accepted, provisionally (2026-11-12).
- **Decided by:** the infrastructure's ADR 0094 (the urban priorities' stage 6, land use and transport interacting), at the user's word ("continue the next feature").

## 1. Context

Where people travel and the way they take are two questions that the social sciences answer with two families of models, used wherever people, goods or money move between places: commuting, shopping, migration, trade and routing. The spatial interaction models (Wilson's family) say how flows between origins and destinations fall with the cost between them; the routing games say how travellers who share a congested network choose their paths, and what their choices cost them all. Each names no object in nature: a flow is a table of origins by destinations, a network a list of arcs with times. The infrastructure's land-use and transport model needs both, and so does any work that asks where flows go on a network that slows as it fills.

## 2. Decision

- **A new module, `interaction`:**
  - `gravity.hpp`, Wilson's family:
    - the production-, attraction- and doubly constrained models;
    - Furness's balancing of any seed table to its margins;
    - Hyman's calibration of β to a mean cost, or to an observed table, which is the Poisson maximum likelihood estimate;
    - the common part of commuters and the SRMSE;
    - the logsum accessibility.
  - `lowry.hpp`: Lowry's model of a metropolis by iteration, and Garin's closed form.
  - `combined.hpp`: Evans's combined equilibrium of distribution and assignment, by partial linearisation with its exact lower bound.
- **`games` gains `congestion.hpp`:** routing games with arcs whose times rise with flow, the Bureau of Public Roads' among them. It finds Wardrop's equilibrium and the system optimum by Frank–Wolfe on Beckmann's potential, with the relative gap and its lower bound, skims of least times, and the price of anarchy. Shortest paths are morphe's Dijkstra, which `games` now uses.

## 3. Alternatives Considered

- **In urban science.** Rejected: trade, migration, retail and freight use the same models, and nothing in them is a street.
- **In tyche.** Rejected: these are models of choice and equilibrium among people, nomos's by decision 0002; tyche's Poisson regression is reached through `econometrics`, whose PPML calibrates the same gravity model by maximum likelihood.
- **Conjugate Frank–Wolfe as the default** (Mitradjieva and Lindberg 2013). Tried, and set aside: on the checks' grid of streets its gap stalled near 4 × 10⁻⁴ for thousands of iterations, while plain Frank–Wolfe kept falling. Path-based methods are the way to faster equilibria, and are not yet written.

## 4. Consequences

- `check_nomos_interaction` (27 checks):
  - Furness's table against its margins and its log-linear form;
  - Wilson's maximum entropy against perturbations that keep the margins and the total cost;
  - Hyman's β against statsmodels' Poisson GLM (fixtures by `make_interaction_fixtures.py`) and against nomos's own PPML;
  - the logsum's derivatives;
  - Lowry against Garin;
  - Evans's equilibrium against its two conditions, and against distributing at free flow and assigning once;
  - Evans's equilibrium on Braess's network, whose road of no time leaves a node and its parent in a shortest-path tree at one distance.
- `check_nomos_congestion` (23 checks):
  - Pigou's and Braess's networks in closed form, and the price of anarchy 4/3;
  - two parallel roads by bisection;
  - a grid of streets against scipy's path-based solution (`make_congestion_fixtures.py`), held within the distance Frank–Wolfe's own bound allows;
  - Roughgarden and Tardos's bound on forty random affine networks.
- Frank–Wolfe and Evans's iteration converge sublinearly: a gap of 10⁻⁴, which forecasting practice asks for, comes quickly, while 10⁻⁸ comes only slowly. The open questions of the two READMEs name the faster methods.
