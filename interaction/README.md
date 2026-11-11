# interaction: Spatial Interaction

*Status: draft (2026-11-12), [decision 0010](../development/decisions/0010-spatial-interaction-and-congestion.md).*

## 1. Concepts

- **Flows that fall with cost.** Origins send flows to destinations, fewer the dearer the journey: commuters to jobs, shoppers to shops, migrants to cities, goods to markets. Wilson (1967, 1970) derived the models from entropy. Of all tables with the observed margins and the observed total cost, the most probable is T_ij = A_i O_i B_j D_j e^{−βc_ij}: the doubly constrained model. Relaxing one margin gives the production-constrained model (shares of each origin's outflow) or the attraction-constrained model.
- **Balancing.** The factors A and B are found by Furness's iteration, which scales rows and columns in turn. It is iterative proportional fitting, Sinkhorn's scaling, and the cyclic Bregman projection on the entropy, and converges whenever a table with the seed's support meets the margins.
- **Calibration.** The mean cost falls as β rises, so one β reproduces an observed mean cost; Hyman (1969) finds it by a secant. Matching the margins and the total cost of an observed table is exactly the score equations of its Poisson likelihood with origin and destination effects (Flowerdew and Aitkin 1982). Hyman's β is therefore the maximum likelihood estimate, which `econometrics`' PPML also gives.
- **Accessibility as a logsum.** A_i = (1/β) ln Σ_j W_j e^{−βc_ij} is the expected greatest utility of a traveller choosing among destinations by the logit, in units of cost (Williams 1977). Its derivative in a cost is minus the share of travellers going there.
- **Lowry's metropolis.** Basic employment brings workers, who bring residents, who need services, whose workers bring residents in turn. Located by interaction shares, the series sums to Garin's closed form e = e^B (I − ab R S)^{−1}.
- **Where travellers go and the way they take, together.** Evans (1976) wrote the doubly constrained table at congested costs and Wardrop's equilibrium for that table as one convex programme: Beckmann's potential plus the entropy over β. Its partial linearisation solves it with an exact lower bound at every iteration.

## 2. Content

| Header | Content |
|---|---|
| `gravity.hpp` | `Table` (origins by destinations); `deterrence`; `production_constrained`, `attraction_constrained`; `balance` (Furness, any seed), `doubly_constrained` → `Balanced{flows, a, b, iterations, residual, converged}`; `mean_cost`; `calibrate(O, D, cost, target)` and `calibrate(observed, cost)` → `Calibrated{beta, fit, mean_cost, …}`; `common_part`, `srmse`; `logsum` |
| `lowry.hpp` | `lowry` (iteration) and `garin` (closed form) → `Metropolis{employment, service, population, …}` |
| `combined.hpp` | `combined(network, origins, destinations, O, D, β, options)` → `Combined{flows, costs, arc_flow, arc_time, objective, lower_bound, gap, …}`, on `games::Network` |

Costs are a dense table; +∞ marks a pair that cannot interact. Units are the caller's: minutes, metres or money.

## 3. Checks

`tests/check_interaction.cpp` (27 checks):

- **Balancing.** Furness's table meets both margins, and ln T + βc has no interaction of row and column.
- **Wilson's theorem.** Fifty perturbations in the null space of the constraints: none nearby has more entropy.
- **The singly constrained tables** meet their margins; with the balancing factors as attractions, the production-constrained table is the doubly constrained one.
- **Calibration.**
  - The mean cost falls as β rises.
  - Hyman's β and fitted table are statsmodels' Poisson GLM with origin and destination indicators, to 10⁻¹⁰ (fixtures by `tests/fixtures/make_interaction_fixtures.py`).
  - The same β comes from `econometrics`' PPML.
- **The logsum's derivatives,** by central differences, and its bounds.
- **The fit measures** at their limits.
- **Lowry.** The iteration agrees with Garin to 10⁻¹³, and the economic base multiplier holds.
- **Evans's equilibrium** on a grid of sixteen crossings: the table is the gravity table at the final least times, the arc flows are Wardrop's equilibrium for it, the bound lies below the objective, and it beats distributing at free flow and assigning once; on Braess's network, whose added road takes no time, its arc flows carry the table and the journey takes 80 minutes.

## 4. Dependencies

`games` (routing games, for `combined.hpp`) and calx `linalg` (Garin's solve). In the checks, `econometrics` (PPML), calx `random` and `verification`.

## 5. References

- Ben-Akiva, M. and Lerman, S. R. (1985), *Discrete Choice Analysis*, MIT Press.
- Bregman, L. M. (1967), The relaxation method of finding the common point of convex sets, *USSR Comput. Math. Math. Phys.* 7, 200–217.
- Evans, S. P. (1976), Derivation and analysis of some models for combining trip distribution and assignment, *Transp. Res.* 10, 37–57.
- Flowerdew, R. and Aitkin, M. (1982), A method of fitting the gravity model based on the Poisson distribution, *J. Reg. Sci.* 22, 191–202.
- Furness, K. P. (1965), Time function iteration, *Traffic Eng. Control* 7, 458–460.
- Garin, R. A. (1966), A matrix formulation of the Lowry model for intrametropolitan activity allocation, *J. Am. Inst. Plann.* 32, 361–364.
- Hyman, G. M. (1969), The calibration of trip distribution models, *Environ. Plan.* 1, 105–112.
- Lenormand, M., Bassolas, A. and Ramasco, J. J. (2016), Systematic comparison of trip distribution laws and models, *J. Transp. Geogr.* 51, 158–169.
- Lowry, I. S. (1964), *A Model of Metropolis*, RM-4035-RC, RAND Corporation.
- Sinkhorn, R. and Knopp, P. (1967), Concerning nonnegative matrices and doubly stochastic matrices, *Pacific J. Math.* 21, 343–348.
- Williams, H. C. W. L. (1977), On the formation of travel demand models and economic evaluation measures of user benefit, *Environ. Plan. A* 9, 285–344.
- Wilson, A. G. (1967), A statistical theory of spatial distribution models, *Transp. Res.* 1, 253–269.
- Wilson, A. G. (1970), *Entropy in Urban and Regional Modelling*, Pion.

## 6. Open Questions

1. **Faster combined equilibria.** Evans's iteration reaches a gap of 10⁻⁴ quickly and 10⁻⁸ slowly. Path-based assignment inside it, or Bar-Gera and Boyce's origin-based algorithm, would converge faster.
2. **Other deterrences.** Power and combined (Tanner) functions, and the radiation and intervening-opportunities laws, when a work needs them.
3. **Harris and Wilson's dynamics** of retail centres, whose stationary states this module's tables already give.
