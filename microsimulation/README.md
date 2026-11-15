# Microsimulation

Status: **draft** (2026-11-15). Module `nomos::microsimulation`, [decision 0013](../development/decisions/0013-microsimulation-and-psychometrics.md).

## 1. Concepts

A static tax-benefit microsimulation applies declared rules to every household of a weighted sample: taxes and contributions, universal amounts, a means-tested benefit, and from them disposable and equivalised income (Sutherland and Figari 2013; Bourguignon and Spadaro 2006).

- **Comparing a reform.** A reform is compared with the baseline in revenue, in inequality and poverty over persons (nomos's `measurement`, the poverty line held fixed), by decile, and in gainers and losers.
- **Incentives.** Marginal effective and participation tax rates are found by applying every rule again, so that interactions are included (Immervoll et al. 2007).
- **Static.** Incomes and behaviour are held as they are, so only the rules change. No rule is any jurisdiction's.

## 2. Content

- `microsimulation.hpp`:
  - `Person`, `Household`, `Bracket`, `Schedule`, `tax`, `Benefit`, `Policy`, `Scale`, `Outcome`, `apply`;
  - `PovertyLine`, `Result`, `simulate`, `Comparison`, `compare`;
  - `marginal_effective_rate`, `participation_rate`, `uprate`.

## 3. Checks

`check_nomos_microsimulation` (16):
- a schedule and a household's outcome under every rule, by hand;
- marginal effective rates equal to what the rules combine to;
- a participation tax rate;
- the accounting identity of a population;
- inequality and poverty equal to nomos's measures;
- a reform's gainers, losers, deciles and revenue consistent;
- uprating as a change of scale that leaves inequality and poverty as they were.

## 4. Not Yet

Behavioural responses (labour supply); dynamic microsimulation (ageing, transitions); reweighting a sample to new totals (peira's raking, from a work); non-take-up.

## 5. References

- Sutherland, H. and Figari, F. (2013), EUROMOD, *International Journal of Microsimulation* 6.
- Bourguignon, F. and Spadaro, A. (2006), Microsimulation as a tool for evaluating redistribution policies, *Journal of Economic Inequality* 4.
- Immervoll, H., Kleven, H. J., Kreiner, C. T. and Saez, E. (2007), Welfare reform in European countries, *Economic Journal* 117.
