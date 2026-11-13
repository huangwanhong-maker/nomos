# appraisal: Choosing Among Alternatives, Over Time, and Under Futures

*Status: draft (2026-11-13; [decision 0012](../development/decisions/0012-appraisal-and-agreement.md)).*

## 1. Concepts

An appraisal compares alternatives, and its methods differ in what they presume. Three presumptions matter most:

- **Compensation.** A weighted sum lets much of one criterion make up for little of another. TOPSIS does too, by distance to an ideal. PROMETHEE allows it only above declared thresholds. A limit applied first allows none: an alternative that breaks it is set aside, whatever else it offers.
- **Time.** Cost–benefit analysis brings benefits and costs to the present by a declared rate, and can weight them by whom they fall on.
- **Futures.** A Delphi panel's rounds record what experts come to agree on. A morphological field records what a declared set of exclusions leaves open. Neither is a forecast.

Weights are the analyst's to declare. The stability of the answer to them is computed, not assumed.

## 2. Content

| Header | What it holds |
|---|---|
| `mcda.hpp` | `screen`, `weighted_sum`, `topsis`, `ahp` (Perron weights, λ_max, CI, CR), `promethee` (PROMETHEE II), `weight_stability`, `acceptability` (SMAA-2) |
| `cba.hpp` | `Schedule` (constant or declining rates), `appraise` (NPV, ratio, payback, annual value, IRR, switching values), `distributional` |
| `futures.hpp` | `delphi` (rounds, consensus, stability), `morphological` (Zwicky's box under Ritchey's cross-consistency) |

## 3. Checks

`check_appraisal.cpp` (24):

- screening;
- the weighted sum, TOPSIS and PROMETHEE II by hand;
- the AHP against numpy's eigensolver, and a consistent matrix exactly;
- the weight stability's crossing by hand;
- SMAA's acceptabilities;
- present values, IRR, payback and switching values by hand, and a declining schedule;
- distributional weights;
- Delphi's rounds;
- a morphological field counted by hand.

## 4. Dependencies and References

- **Depends on:** calx (`linalg`: the Perron vector; `random`; `solvers`: Brent's method for the IRR).
- **References:**
  - Hwang and Yoon (1981), *Multiple Attribute Decision Making*;
  - Saaty (1980), *The Analytic Hierarchy Process*;
  - Brans and Vincke (1985), *Manage. Sci.* 31, 647–656;
  - Lahdelma and Salminen (2001), *Oper. Res.* 49, 444–454;
  - Boardman et al. (2018), *Cost–Benefit Analysis*, 5th ed.;
  - HM Treasury (2022), *The Green Book*;
  - Ritchey (2011), *Wicked Problems – Social Messes*;
  - von der Gracht (2012), *Technol. Forecast. Soc. Change* 79, 1525–1536.
