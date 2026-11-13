# 0012. Appraisal, Futures and the Agreement of Coders

- **Status:** Accepted, provisionally (2026-11-13).
- **Date:** 2026-11-13
- **Decided by:** the infrastructure's programme of the handbook *Material and Ecological Conditions of Life*, at the user's request to implement its methods. The handbook's §9 asks for multi-criteria decision analysis and cost–benefit analysis. It warns against "treating everything as a compensable weighted score". Its §10 asks for Delphi and morphological analysis, and its §8 rests on coded interviews. These are methods of the social sciences, nomos's by decision 0002.

## 1. Context

nomos held composite indicators (`measurement`, decision 0006), social choice and the Shapley value. It had no way to compare alternatives on several criteria with declared weights and test how far the weights may move. It had no present values, no record of a panel's rounds or of the space a set of exclusions leaves, and no measure of whether coders agree.

## 2. Decision

- **A new module `appraisal`:**
  - `mcda.hpp`:
    - `screen`, non-compensatory limits applied first;
    - `weighted_sum`, `topsis`, `promethee` (PROMETHEE II with five preference functions);
    - `ahp`, weights by calx's Perron vector, with λ_max, CI and CR against Saaty's random indices;
    - `weight_stability`, the interval of each weight over which the best stays best, exactly;
    - `acceptability`, SMAA-2's rank acceptability under weights uniform on the simplex.
  - `cba.hpp`:
    - `Schedule`, a constant or declining rate year by year;
    - `appraise`: present values, NPV, benefit–cost ratio, discounted payback, equivalent annual value, IRR where the net stream changes sign once, switching values;
    - `distributional`, weights (ȳ/y)^η.
  - `futures.hpp`:
    - `delphi`, rounds with their medians, interquartile ranges, trimmed means, consensus and stability by declared criteria;
    - `morphological`, Zwicky's box reduced by Ritchey's cross-consistency to the configurations with no excluded pair, counted.
- **`measurement/agreement.hpp`:**
  - `cohen`, κ unweighted and weighted linearly or quadratically;
  - `fleiss`, κ for many coders;
  - `krippendorff`, α at the nominal, ordinal, interval and ratio levels, with missing values.

## 3. Alternatives Considered

- **MCDA in `measurement`, beside composite indicators.** Rejected: a composite indicator measures units; an appraisal chooses among alternatives, with screening and stability that a measure does not have.
- **Cost–benefit analysis in the economics domain.** Rejected: discounting and distributional weights are general methods of welfare analysis, used by every domain that appraises; the domain declares rates and values.
- **Agreement in tyche.** The statistics are general, but their established use is the reliability of coding in the social sciences: nomos's by decision 0002.

## 4. Consequences

- **Checks.**
  - Agreement (10): held to statsmodels and the krippendorff package, Krippendorff's own example among them.
  - Appraisal (24): by hand, and the AHP against numpy's eigensolver.
- **Not built:**
  - ELECTRE's outranking with vetoes;
  - scenario discovery (PRIM), whose place is tyche's `uncertainty`;
  - real options;
  - the robustness of decisions across scenarios (minimax regret).
