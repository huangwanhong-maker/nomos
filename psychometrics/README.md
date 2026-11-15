# Psychometrics

Status: **draft** (2026-11-15). Module `nomos::psychometrics`, [decision 0013](../development/decisions/0013-microsimulation-and-psychometrics.md).

## 1. Concepts

- **Classical test theory.** How reliably a set of items measures one trait: Cronbach's α, each item's correlation with the rest, α without it, Spearman and Brown's prophecy, split-half reliability, the standard error of measurement (Lord and Novick 1968).
- **Item response theory.** Rasch's model and the two-parameter logistic model, fitted by marginal maximum likelihood (Bock and Aitkin's EM on a grid of 61 points over [−6, 6] weighted by the normal density, as mirt does). Its log-likelihood never decreases.
- **Measuring persons.** Each person's trait by its expected a posteriori value. Item and test information say where the scale measures well.

## 2. Content

- `psychometrics.hpp`:
  - `cronbach_alpha`, `ItemAnalysis`, `item_analysis`, `spearman_brown`, `split_half`, `standard_error_of_measurement`;
  - `Item`, `probability`, `information`, `Grid`, `grid`, `Model`, `Fit`, `marginal_log_likelihood`, `fit`, `Score`, `score`.

## 3. Checks

`check_nomos_psychometrics` (17):
- α = 10/11 by hand, with rest correlations and α without each item;
- the prophecy and the standard error of measurement;
- information greatest at an item's difficulty;
- ten items recovered from 2,000 simulated persons, the EM monotone and its estimate at least as likely as the truth;
- Rasch's model recovering a common discrimination of one, and persons with the same number right given the same trait;
- symmetric items, and an empty record taking the prior.

## 4. Not Yet

The three-parameter model; polytomous items (graded response, partial credit); McDonald's ω (on tyche's `latent`); differential item functioning; adaptive testing.

## 5. References

- Cronbach, L. J. (1951), Coefficient alpha and the internal structure of tests, *Psychometrika* 16.
- Lord, F. M. and Novick, M. R. (1968), *Statistical Theories of Mental Test Scores*.
- Bock, R. D. and Aitkin, M. (1981), Marginal maximum likelihood estimation of item parameters, *Psychometrika* 46.
- Chalmers, R. P. (2012), mirt, *Journal of Statistical Software* 48(6).
