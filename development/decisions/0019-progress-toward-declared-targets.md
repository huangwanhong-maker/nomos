# 0019. Progress Toward Declared Targets

- **Status:** Accepted, provisionally (2026-11-26).
- **Date:** 2026-11-26
- **Decided by:** the infrastructure's case of a book written by communities (its design record of 2026-11-26), whose analysts assess progress toward the targets of the Sustainable Development Goals with the shared bench.

## 1. Context

How far a quantity has come toward a target it should reach by a date is assessed in one standard way across the social sciences' monitoring of goals: from a baseline and a latest observation, a compound annual rate, the value it would reach at the target year if the rate held, and the share of the distance to the target that this would cover, classed by declared floors. The method names no goal and no country; nomos lacked it.

## 2. Decision

`measurement/include/nomos/measurement/progress.hpp`:

- `compound_rate(baseline, latest, years)`;
- `trend(b, x_b, t, x_t, desired, scheme, has_target, target, target_year)` → `Progress{rate, expected, ratio, met, level, trend, basis}`: the ratio ρ = (x̂_T − x_b)/(x* − x_b) classed by the scheme's floors, the target met when the latest value is at or beyond it, the best or the worst class when the baseline equals the target or already lay beyond it (a baseline beyond the target that the latest value no longer meets is regression, where ρ's sign would read the reversal as progress; found by the infrastructure's reading for error of ADR 0172); without a target the rate is classed, its sign turned for a quantity to be lowered;
- `baseline_index`, `latest_index`, `trend_of_series(years, values, …)`: the baseline the declared base year if observed, else the earliest year after it;
- `ProgressScheme`: the floors and the names of the classes, and the base year; its default the United Nations Statistics Division's for the global assessment (UN DESA Statistics Division, *Technical Note for Progress Assessment 2025*, Section 3.1): floors 0.95, 0.50, 0.10, 0 on ρ and 2, 1.25, 0.5, 0 per cent a year on r, five classes, base year 2015.

The note's scores from 0 to 5 and their means by target are left to a caller.

## 3. Alternatives Considered

- **Linear rates.** Some regional commissions extrapolate linearly; the global method is compound, and a caller who wants another form composes it. The caveat is stated in the header.
- **Targets of zero replaced by a small number.** The ratio is defined for a target of zero (the distance is the baseline); no substitute is needed.

## 4. Consequences

- The API's `science/progress.hpp` composes it with the spread among units, convergence, breaks and links, for any work.
- Checks: `check_nomos_progress` (20): the rate by hand; a quantity halved toward zero (ρ = 0.875); one raised toward 100; regression; a baseline beyond the target that fell back from it (regression, where the ratio's sign would read progress); the classes at their floors; the rate's classes both ways; a baseline at the target; a target met; the baseline year chosen; invariance under a change of unit; a compound path reaching its target giving ρ = 1; refusals.
