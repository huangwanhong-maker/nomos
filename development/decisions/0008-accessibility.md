# 0008. Accessibility

- **Status:** Accepted, provisionally (2026-11-08).
- **Decided by:** the infrastructure's ADR 0086 (the urban priorities' stage 2), at the user's word ("please implement the next priority one").

## 1. Context

How much of what people need lies within their reach, and who is left without it, is asked of clinics, schools, shops, jobs and parks, in cities and in the countryside, by planners, geographers, health services and development banks alike. The measures name no object in nature. They take places of demand, places of supply and the travel costs between them, however those costs were found. They are measures of the social sciences, beside inequality and poverty, and so are nomos's.

## 2. Decision

`measurement` gains `accessibility.hpp`, on a sparse list of the pairs within reach (a demand, a supply, a cost) that the caller computes on a network:

- **The decays of a cost:** binary, exponential and Gaussian within a limit (Dai 2010).
- **Cumulative opportunities and gravity** (Wachs and Kumagai 1973; Hansen 1959).
- **The two-step floating catchment** with a decay in both steps (Luo and Wang 2003; Luo and Qi 2009). Each supply's ratio to the demand that reaches it, then each place's share of the ratios within its reach. It conserves supply: Σ_i P_i A_i is the supply some demand reaches.
- **Unmet need:** the people with nothing within reach, those below a declared share of the mean, and the population-weighted Gini coefficient of access.
- **The Rural Access Index** (Roberts, KC and Rastogi 2006; SDG 9.1.1): the share of rural people within 2 km of an all-season road.

## 3. Alternatives Considered

- **In urban science.** Rejected: rural access, hospitals and labour markets need the same measures, and nothing in them is a street.
- **In tyche.** Rejected: these are measures of a population, as inequality is, not inference.
- **Travel costs inside nomos.** Rejected: the network is the caller's (urban science's `access` on morphe's shortest paths), so the measures serve any cost: minutes, metres or money.

## 4. Consequences

- `check_accessibility.cpp` holds the measures to a case by hand, to the conservation law under every decay, to gravity's limit as cumulative opportunities, and to unmet need and the index by hand.
- Later, if asked: three-step and variable catchments, competition for supply (the inverted 2SFCA), and access over time of day.
