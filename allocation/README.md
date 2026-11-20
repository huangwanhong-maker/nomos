# allocation

Status: **draft** (2026-11-20). Module `nomos::allocation` ([decision 0015](../development/decisions/0015-allocation-and-welfare-orders.md)).

## 1. Concepts

When there is not enough, someone decides who is served first, what each receives and how much. These are not one optimisation: each is decided by a rule, and a rule is known by the properties it has and lacks. This module holds the rules that decide each part, each with its source and with its properties checked:

- **How much** of a divisible sum, when the claims on it exceed it. These are the rationing rules beyond `games`'s claims problems: Piniles's rule, adjusted proportional, the TAL family between equal losses and equal awards, duality, weighted proportional and priority classes. They come with the formulas by which public bodies divide a sum among places:
  - sufficiency first, closing shortfalls to a declared threshold and switching from shares in plenty to survival shares in scarcity;
  - Porto Alegre's banded points;
  - fiscal equalisation of needs against capacities, with its envelope reported apart.
- **Who first.** Reserve systems give each category of units its own priority order and process the categories in a declared precedence (Pathak, Sönmez, Ünver and Yenmez). Their three axioms are checked, and the beneficiaries under every precedence order are reported, since the order is a choice.
- **Who gets which indivisible good.** The options are serial dictatorship, the random serial dictatorship, probabilistic serial and picking sequences. A random assignment is carried out as a lottery over deterministic ones (Birkhoff).
- **Whether a division of goods is fair.** It is tested for envy-freeness, envy-freeness up to one good, EFX, proportionality and the maximin share. The procedures here (round-robin, envy-cycle elimination, the maximum Nash welfare) have known guarantees.
- **What to fund from a budget by votes.** The rules are utilitarian greedy, the Method of Equal Shares and sequential Phragmén, with a check of extended justified representation.

Every rule returns its parts and its reasons, so that an award can be contested on its own ground. The infrastructure's survey *Allocative justice under scarcity* (2026, Part II) collects the literature.

## 2. Content

- `rationing.hpp`: `piniles`, `adjusted_proportional`, `tal(θ)`, `dual(R)`, `weighted_proportional`, `priority_classes`; `Parts`, `sufficiency_first`; `Indicator`, `band_grade`, `Points`, `banded_points`; `Unit`, `Envelope`, `equalisation`.
- `reserves.hpp`: `Category`, `ReserveMatching`, `sequential_reserve`, `AxiomReport`, `check_axioms` (eligibility, non-wastefulness, respect of priorities), `OrderOutcome`, `every_order`.
- `assignment.hpp`: `serial_dictatorship`, `RandomAssignment`, `random_serial` (exact to eight agents, else sampled), `probabilistic_serial`, `picking_sequence`, `round_robin_sequence`; `Draw`, `birkhoff`, `Squared`, `square`, `draw_index`.
- `division.hpp`: `Values`, `Division`, `bundles_of`, `value_of`, `FairnessReport`, `maximin_shares`, `check_fairness`, `round_robin`, `envy_cycle`, `max_nash_welfare` (an optimisation, by enumeration, said so).
- `budgeting.hpp`: `Election`, `Satisfaction`, `Selection`, `utilitarian_greedy`, `equal_shares` (with a declared completion), `phragmen`, `EjrReport`, `check_ejr`.

*Caveat.*
- Enumerations (maximin shares, the maximum Nash welfare, the random serial dictatorship beyond eight agents, EJR) are for small instances, and refuse larger ones or say they sample.
- With unequal costs the Method of Equal Shares guarantees EJR up to one project, not EJR.
- The three axioms of a reserve system admit waste of a reserved unit (hard reserves processed open first). The smart reserve matching that avoids it is not yet built.

## 3. Checks

37 (`tests/check_allocation.cpp`):

- **Rationing.** The TAL family at θ = ½ is the Talmud rule on the Talmud's three estates, and at 0 and 1 the equal losses and awards. Piniles's rule is checked by hand. Adjusted proportional equals the contested garment for two claimants. On 200 random problems, duality holds and every rule is feasible. Weighted proportional is capped at the claims, and classes are served in order.
- **Formulas.** Sufficiency first behaves correctly in plenty and in scarcity. Banded points are checked by hand. Equalisation is checked by hand, and its grants sum to zero. An envelope's floor is kept.
- **Reserves.** On random instances the sequential matching satisfies the three axioms, and it is among all matchings that do, found by enumerating them. Hard reserves processed open first leave the reserved unit idle, though the axioms hold. The beneficiaries under every order are reported.
- **Assignment.** Probabilistic serial gives two agents of one preference half of each object, and is stochastic-dominance envy-free on random instances. Random assignments keep within their bounds, and their lotteries sum back. The random serial dictatorship is exact for eight agents and sampled for nine, against counting.
- **Division.** Round-robin, envy-cycle elimination and the maximum Nash welfare are EF1 on 80 random instances. An unfair division is named, with its witnesses.
- **Budgeting.** A majority of six and a minority of four: utilitarian greedy gives the minority nothing and breaks EJR, while equal shares and Phragmén give it four, and equal shares satisfies EJR.

## 4. Dependencies

`nomos::games` (the claims rules); calx's `random`.

## 5. References

- **Claims.** O'Neill (1982); Aumann and Maschler (1985); Curiel, Maschler and Tijs (1987); Moulin (2000); Thomson (2003, 2015); Moreno-Ternero and Villar (2006).
- **Formulas.** Santos (1998).
- **Reserves.** Pathak, Sönmez, Ünver and Yenmez (2024; arXiv:2008.00374).
- **Assignment.** Bogomolnaia and Moulin (2001); Birkhoff (1946); Bouveret and Lang (2011).
- **Division.** Lipton, Markakis, Mossel and Saberi (2004); Budish (2011); Caragiannis et al. (2019); Kurokawa, Procaccia and Wang (2018).
- **Budgeting.** Peters and Skowron (2020); Peters, Pierczyński and Skowron (2021); Brill, Freeman, Janson and Lackner (2017).
