# 0018. The Structure of Regional Economies

- **Status:** Accepted, provisionally (2026-11-23).
- **Date:** 2026-11-23
- **Decided by:** the user's question of what economics remains for urban science, and the infrastructure's survey *Evolutionary economics models and methods* (its Part IV and §7.2, item 12), for the economics domain's planned `diversification` (the infrastructure's ADR 0165).

## 1. Context

A city's or region's economy is read by how its employment changed against its reference, how varied and how related its activities are, and which activities lie near those it holds. The measures are general: they name units and categories and no particular place, and they serve geography, economics and urban science alike. tyche already holds the ordering that the economic complexity index is (its decision 0033) and nomos the location quotients; shift-share, related and unrelated variety, relatedness and fitness were missing.

## 2. Decision

`measurement/include/nomos/measurement/structure.hpp`:

- `shift_share(e0, e1, E0, E1)`: Dunn's national share, industrial mix and competitive terms, and Esteban-Marquillas's competitive effect and allocation effect where every base is positive;
- `variety(amounts, group, base)`: entropy as unrelated variety plus related variety (Frenken, Van Oort and Verburg 2007), the chain rule exact;
- `incidence(table, threshold)` and `from_binary(M)`, `proximity`, `density` (Hidalgo et al. 2007; the activity itself left out of its density, a stated convention);
- `complexity(incidence)`: the indices of regions and activities on tyche's `bipartite_ordering`, the activities' as the mean index of the regions that hold them; an incidence that falls apart into blocks refused;
- `fitness_complexity(incidence, steps)` (Tacchella et al. 2012), with the number of steps declared and reported.

`nomos::measurement` now links tyche's `networks`.

## 3. Alternatives Considered

- **The economic complexity index written again.** tyche holds it as the bipartite ordering; nomos composes it.
- **The activities' index by its own eigenvector.** Its sign is free; taking it as the mean regional index of the regions that hold an activity (Hidalgo and Hausmann's relation) fixes the sign with the regions', and the relation is exact.
- **Servedio et al.'s stable form of fitness.** Its exact equations were not read at the source; the original iteration is given with its caveat.

## 4. Consequences

- `check_nomos_structure`, 19 checks: shift-share and its split worked by hand and exact on random regions; variety's chain rule and its values on uniform shares; proximity and density by hand; the raw complexity index orthogonal to diversity and the regions' index affine in the mean index of their activities, on random tables; fitness respecting nestedness; refusals, among them an incidence fallen apart into blocks, found while checking: such an incidence had given an index that only named its blocks.
