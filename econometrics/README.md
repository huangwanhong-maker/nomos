# econometrics: Discrete Choice and Fixed-Effects Poisson

*Status: draft (2026-10-25; reduced 2026-11-08 by [decision 0007](../development/decisions/0007-panels-instruments-and-differences-move.md)).* Built on tyche's `regression`: every generalised linear model is tyche's `glm`, and the absorption of fixed effects and the sandwich of robust covariances are tyche's (`tyche::panel`, `tyche::regression`), named through `common.hpp`. Panels moved to tyche's `panel`, instruments and differences in differences to peira's `quasi`, with their checks.

## 1. Concepts

- **Discrete choice.** Average marginal effects of binary logit and probit models, derivatives for continuous regressors and differences for binary ones, with delta-method standard errors from the fit's covariance (Wooldridge 2010, ch. 15).
- **The multinomial logit.** The probabilities of J categories as a softmax of linear predictors, fitted by quasi-Newton maximum likelihood, with model-based, robust and clustered covariances and the marginal effects of each regressor on each probability.
- **Poisson pseudo-maximum likelihood with two-way fixed effects.** For a two-way table of nonnegative outcomes with E(y | x) = exp(xβ + α_r + γ_c), the Poisson quasi-likelihood is consistent whatever the distribution of y and keeps the zeros (Gourieroux, Monfort and Trognon 1984; Santos Silva and Tenreyro 2006). It is fitted by iteratively reweighted least squares, each step absorbing the row and column effects by weighted alternating projections (Correia, Guimarães and Zylkin 2020), so a table of 200 × 200 cells never forms its 400 dummies. Rows or columns with no positive outcome are separated and dropped. Its robust (HC0) or clustered covariance is exactly the β block of the full dummy-variable sandwich.

## 2. Content

| Header | Content |
|---|---|
| `common.hpp` | The names read from tyche: `Factor`, `factor`, `AbsorbOptions`, `demean`, `recover_effects`; the sandwich's pieces |
| `discrete_choice.hpp` | `average_marginal_effects(glm_fit, X, discrete_columns, confidence)` → `MarginalEffects{effect, se, z, p, ci, discrete, covariance}` |
| `multinomial.hpp` | `MultinomialCovariance {model, hc0, cluster}`, `MultinomialOptions`; `multinomial_logit(X, y, categories, options)` → `MultinomialFit`; `multinomial_marginal_effects(fit, X)` |
| `poisson_fixed_effects.hpp` | `PoissonCovariance {hc0, cluster}`, `PoissonTwoWayOptions{covariance, clusters, offset, tolerance, max_iterations, confidence, absorb}`; `poisson_two_way(X, y, row, column, options)` → `PoissonTwoWayFit{beta, se, z, p, ci, covariance, aliased, row_effects, column_effects, fitted, dropped, used, dropped_rows, dropped_columns, clusters, deviance, iterations, converged, sweeps}` |

## 3. Checks

- `tests/check_discrete_choice.cpp`: for logit and probit, the continuous effect as the derivative of the average probability and as β times the average density, the binary effect as an average difference, the delta-method errors against a gradient by central differences, the intercept without effect, the truth within four standard errors.
- `tests/check_multinomial.cpp`: the multinomial logit's coefficients, its covariances and marginal effects.
- `tests/check_poisson_fixed_effects.cpp`: the independence model μ_rc = y_r· y_·c / y_·· with an all-zero row and column dropped; against tyche's GLM with dense dummies on a 10 × 12 table with zeros (coefficients, deviance, fitted values, effects, HC0 errors, clustered errors in the ratio of the small-sample factors), a regressor absorbed by the row effects; a gravity-like table of 39,800 cells, converged with the truth within four standard errors.

## 4. Dependencies

tyche `panel` (absorbed effects), `regression` (`ols`, `glm`, the sandwich) and `probability`; calx `linalg` and `optimisation`. In the checks, calx `random` and `verification`.

## 5. References

- Cameron, A. C. and Trivedi, P. K. (2005), *Microeconometrics*, Cambridge University Press.
- Correia, S., Guimarães, P. and Zylkin, T. (2019), Verifying the existence of maximum likelihood estimates for generalized linear models, arXiv:1903.01633; (2020), Fast Poisson estimation with high-dimensional fixed effects, *Stata J.* 20, 95–115.
- Fernández-Val, I. and Weidner, M. (2016), Individual and time effects in nonlinear panel models with large N, T, *J. Econometrics* 192, 291–312.
- Gourieroux, C., Monfort, A. and Trognon, A. (1984), Pseudo maximum likelihood methods: applications to Poisson models, *Econometrica* 52, 701–720.
- McCullagh, P. and Nelder, J. A. (1989), *Generalized Linear Models*, 2nd ed., Chapman and Hall.
- Santos Silva, J. M. C. and Tenreyro, S. (2006), The log of gravity, *Rev. Econ. Stat.* 88, 641–658; (2010), On the existence of the maximum likelihood estimates in Poisson regression, *Econ. Lett.* 107, 310–312.
- Weidner, M. and Zylkin, T. (2021), Bias and consistency in three-way gravity models, *J. Int. Econ.* 132, 103513.
- White, H. (1982), Maximum likelihood estimation of misspecified models, *Econometrica* 50, 1–25.
- Wooldridge, J. M. (2010), *Econometric Analysis of Cross Section and Panel Data*, 2nd ed., MIT Press.

## 6. Open Questions

1. **Acceleration** of the projections (Aitken or conjugate gradients, as reghdfe and ppmlhdfe do) and singleton dropping, when tables grow beyond a few hundred thousand cells; it belongs to tyche's `panel`.
2. **Conditional and ordered logit** and their marginal effects.
