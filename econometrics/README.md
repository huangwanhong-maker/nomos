# econometrics: Panels, Instruments, Differences in Differences, Discrete Choice and Fixed-Effects Poisson

*Status: draft (2026-10-25), [decision 0003](../development/decisions/0003-the-first-modules.md).* Built on tyche's `regression`: every least-squares fit is tyche's `ols` (pivoted QR, never the normal equations) and every generalised linear model tyche's `glm`; this module adds the designs of econometrics around them. Multinomial logit, regression discontinuity, durations and the estimators of staggered adoption are planned.

## 1. Concepts

- **Absorbed fixed effects.** The effects of one or more factors (units, periods, origins, destinations) are removed by the theorem of Frisch, Waugh and Lovell: the coefficients on X of a regression on X and the factors' indicators are those of the regression of the residualised outcome on the residualised X, with the same residuals and the same sandwich. The residualising is done without forming the indicators, by alternating weighted projections (Guimarães and Portugal 2010; Gaure 2013), exact in one sweep for one factor. The parameters absorbed are L for one factor and L₁ + L₂ − (connected components) for two (Abowd, Creecy and Kramarz 2002).
- **Panel data.** The within estimator treats unit (and period) effects as parameters; its classical standard errors count the absorbed parameters in the degrees of freedom, and its cluster-robust ones count only those not nested in the clusters (Cameron and Miller 2015), as Stata's `xtreg` and `reghdfe` do. The between estimator regresses unit means. Random effects is feasible GLS on quasi-demeaned data with Swamy and Arora's (1972) variance components; in a balanced panel it is a matrix-weighted average of within and between (Maddala 1971). Hausman's (1978) test compares fixed and random effects.
- **Instrumental variables.** Two-stage least squares replaces the endogenous regressors by their projections on the instruments (Theil 1953); just identified it is (ZᵀX)⁻¹Zᵀy. The first stage reports the F statistic of the excluded instruments and their partial R² (weak instruments: Staiger and Stock 1997; Stock and Yogo 2005). Overidentifying restrictions are tested by Sargan's (1958) statistic and by Hansen's (1982) J at the two-step efficient GMM estimate.
- **Differences in differences.** The 2 × 2 estimate is the difference of the treated and control groups' changes, equal to the interaction coefficient of the saturated regression, with standard errors clustered by unit (Bertrand, Duflo and Mullainathan 2004); with many periods, the two-way fixed-effects regression, which with staggered adoption and heterogeneous effects need not estimate an average effect (Goodman-Bacon 2021).
- **Discrete choice.** Average marginal effects of binary logit and probit models, derivatives for continuous regressors and differences for binary ones, with delta-method standard errors from the fit's covariance (Wooldridge 2010, ch. 15).
- **Poisson pseudo-maximum likelihood with two-way fixed effects.** For a two-way table of nonnegative outcomes with E(y | x) = exp(xβ + α_r + γ_c), the Poisson quasi-likelihood is consistent whatever the distribution of y and keeps the zeros (Gourieroux, Monfort and Trognon 1984; Santos Silva and Tenreyro 2006). It is fitted by iteratively reweighted least squares, each step absorbing the row and column effects by weighted alternating projections (Correia, Guimarães and Zylkin 2020), so a table of 200 × 200 cells never forms its 400 dummies. Rows or columns with no positive outcome are separated and dropped. Its robust (HC0) or clustered covariance is exactly the β block of the full dummy-variable sandwich.

## 2. Content

| Header | Content |
|---|---|
| `absorb.hpp` | `Factor{level, levels}`, `factor(labels)`; `AbsorbOptions{tolerance, max_sweeps}`; `demean(columns, factors, w, options)`; `components(a, b)`, `occurring(f)`, `absorbed_rank(factors)`, `nested(f, clusters)`; `recover_effects(d, factors, w, options)` |
| `panel.hpp` | `PanelEffects {unit, two_way}`, `PanelCovariance {classical, cluster}`, `PanelOptions`; `within(X, y, unit, period, options)` → `WithinFit{fit, unit_effects, period_effects, unit_labels, period_labels, absorbed, absorbed_not_nested, sweeps}`; `between(X, y, unit, ols_options)` → `tyche::regression::LinearFit`; `random_effects(X, y, unit, ols_options)` → `RandomEffectsFit{fit, sigma_u, sigma_e, rho, theta}`; `hausman(fe, re)` → `HausmanTest{statistic, df, p_value, positive_definite}` |
| `instrumental.hpp` | `IvCovariance {classical, hc0, hc1, cluster}`, `IvOptions`; `two_stage_least_squares(y, exogenous, endogenous, instruments, options)` → `IvFit{beta, se, t, p, ci, covariance, residuals, first_stage[FirstStage{f_statistic, df1, df2, p_value, partial_r_squared}], sargan, hansen_j, gmm_beta, …}` |
| `difference_in_differences.hpp` | `difference_in_differences(y, treated, post, units, confidence)` → `DifferenceInDifferences{estimate, se, t, p, ci, mean, count, regression}`; `two_way_fixed_effects(y, treatment, unit, period, controls, confidence)` → `WithinFit` |
| `discrete_choice.hpp` | `average_marginal_effects(glm_fit, X, discrete_columns, confidence)` → `MarginalEffects{effect, se, z, p, ci, discrete, covariance}` |
| `poisson_fixed_effects.hpp` | `PoissonCovariance {hc0, cluster}`, `PoissonTwoWayOptions{covariance, clusters, offset, tolerance, max_iterations, confidence, absorb}`; `poisson_two_way(X, y, row, column, options)` → `PoissonTwoWayFit{beta, se, z, p, ci, covariance, aliased, row_effects, column_effects, fitted, dropped, used, dropped_rows, dropped_columns, clusters, deviance, iterations, converged, sweeps}` |

The within estimator's `fit` is a `tyche::regression::LinearFit` whose covariance, degrees of freedom and inference have been corrected for the absorbed parameters, so tyche's `wald` applies to it; its `r_squared` is the within R².

## 3. Checks

- `tests/check_panel.cpp`: one-way within against the dummy-variable regression (coefficients, classical errors, residuals, the effects; clustered errors in the exact ratio √((n − k − N)/(n − k))); two-way on an unbalanced panel (coefficients, classical errors, fitted values, the effects rebuilding them, N + T − 1 absorbed, T − 1 of them not nested in unit clusters, the truth recovered, a time-invariant regressor absorbed); random effects as Maddala's weighted average, exactly; Swamy–Arora components near the truth; Hausman accepting and rejecting; differences in differences, 2 × 2 against the four means and the truth, and two-way fixed effects giving the same estimate.
- `tests/check_instrumental.cpp`: just identified, (ZᵀX)⁻¹Zᵀy and its classical covariance; overidentified, OLS biased and 2SLS recovering the truth, a strong first stage, Sargan and Hansen accepting and nearly agreeing, clustered errors; an invalid instrument rejected by both; a weak instrument's first-stage F below 10.
- `tests/check_discrete_choice.cpp`: for logit and probit, the continuous effect as the derivative of the average probability and as β times the average density, the binary effect as an average difference, the delta-method errors against a gradient by central differences, the intercept without effect, the truth within four standard errors.
- `tests/check_poisson_fixed_effects.cpp`: the independence model μ_rc = y_r· y_·c / y_·· with an all-zero row and column dropped; against tyche's GLM with dense dummies on a 10 × 12 table with zeros (coefficients, deviance, fitted values, effects, HC0 errors, clustered errors in the ratio of the small-sample factors), a regressor absorbed by the row effects; a gravity-like table of 39,800 cells, converged with the truth within four standard errors.

The counts are in the release notes of the commit that adds them.

## 4. Dependencies

tyche `regression` (`ols`, `wald`, `glm`) and `probability` (Student's t, χ², the normal); calx `linalg` (`Rectangular`, `Matrix`, `solve`, `cholesky`). In the checks, calx `random` and `verification`.

## 5. References

- Abowd, J. M., Creecy, R. H. and Kramarz, F. (2002), Computing person and firm effects using linked longitudinal employer-employee data, US Census Bureau, Technical Paper 2002-06.
- Angrist, J. D. and Pischke, J.-S. (2009), *Mostly Harmless Econometrics*, Princeton University Press.
- Baltagi, B. H. (2013), *Econometric Analysis of Panel Data*, 5th ed., Wiley.
- Bertrand, M., Duflo, E. and Mullainathan, S. (2004), How much should we trust differences-in-differences estimates?, *Q. J. Econ.* 119, 249–275.
- Bound, J., Jaeger, D. A. and Baker, R. M. (1995), Problems with instrumental variables estimation when the correlation between the instruments and the endogenous explanatory variable is weak, *J. Am. Stat. Assoc.* 90, 443–450.
- Cameron, A. C. and Miller, D. L. (2015), A practitioner's guide to cluster-robust inference, *J. Hum. Resour.* 50, 317–372.
- Cameron, A. C. and Trivedi, P. K. (2005), *Microeconometrics*, Cambridge University Press.
- Card, D. and Krueger, A. B. (1994), Minimum wages and employment, *Am. Econ. Rev.* 84, 772–793.
- Correia, S. (2016), Linear models with high-dimensional fixed effects: an efficient and feasible estimator, working paper (reghdfe).
- Correia, S., Guimarães, P. and Zylkin, T. (2019), Verifying the existence of maximum likelihood estimates for generalized linear models, arXiv:1903.01633; (2020), Fast Poisson estimation with high-dimensional fixed effects, *Stata J.* 20, 95–115.
- de Chaisemartin, C. and D'Haultfœuille, X. (2020), Two-way fixed effects estimators with heterogeneous treatment effects, *Am. Econ. Rev.* 110, 2964–2996.
- Fernández-Val, I. and Weidner, M. (2016), Individual and time effects in nonlinear panel models with large N, T, *J. Econometrics* 192, 291–312.
- Frisch, R. and Waugh, F. V. (1933), Partial time regressions as compared with individual trends, *Econometrica* 1, 387–401; Lovell, M. C. (1963), Seasonal adjustment of economic time series and multiple regression analysis, *J. Am. Stat. Assoc.* 58, 993–1010.
- Gaure, S. (2013), OLS with multiple high dimensional category variables, *Comput. Stat. Data Anal.* 66, 8–18.
- Goodman-Bacon, A. (2021), Difference-in-differences with variation in treatment timing, *J. Econometrics* 225, 254–277.
- Gourieroux, C., Monfort, A. and Trognon, A. (1984), Pseudo maximum likelihood methods: applications to Poisson models, *Econometrica* 52, 701–720.
- Guimarães, P. and Portugal, P. (2010), A simple feasible procedure to fit models with high-dimensional fixed effects, *Stata J.* 10, 628–649.
- Hansen, L. P. (1982), Large sample properties of generalized method of moments estimators, *Econometrica* 50, 1029–1054.
- Hausman, J. A. (1978), Specification tests in econometrics, *Econometrica* 46, 1251–1271.
- Hayashi, F. (2000), *Econometrics*, Princeton University Press.
- Maddala, G. S. (1971), The use of variance components models in pooling cross section and time series data, *Econometrica* 39, 341–358.
- McCullagh, P. and Nelder, J. A. (1989), *Generalized Linear Models*, 2nd ed., Chapman and Hall.
- Santos Silva, J. M. C. and Tenreyro, S. (2006), The log of gravity, *Rev. Econ. Stat.* 88, 641–658; (2010), On the existence of the maximum likelihood estimates in Poisson regression, *Econ. Lett.* 107, 310–312.
- Sargan, J. D. (1958), The estimation of economic relationships using instrumental variables, *Econometrica* 26, 393–415.
- Staiger, D. and Stock, J. H. (1997), Instrumental variables regression with weak instruments, *Econometrica* 65, 557–586.
- Stock, J. H. and Yogo, M. (2005), Testing for weak instruments in linear IV regression, in D. W. K. Andrews and J. H. Stock (eds), *Identification and Inference for Econometric Models*, Cambridge University Press, 80–108.
- Swamy, P. A. V. B. and Arora, S. S. (1972), The exact finite sample properties of the estimators of coefficients in the error components regression models, *Econometrica* 40, 261–275.
- Theil, H. (1953), Repeated least squares applied to complete equation systems, Central Planning Bureau, The Hague (mimeo).
- Weidner, M. and Zylkin, T. (2021), Bias and consistency in three-way gravity models, *J. Int. Econ.* 132, 103513.
- White, H. (1982), Maximum likelihood estimation of misspecified models, *Econometrica* 50, 1–25.
- Wooldridge, J. M. (2010), *Econometric Analysis of Cross Section and Panel Data*, 2nd ed., MIT Press.

## 6. Open Questions

1. **Where absorption belongs.** The alternating projections name no social object and could serve any science's regression with many categorical effects; if another domain needs them, they move to tyche's `regression` and this module uses them there.
2. **Acceleration** of the projections (Aitken or conjugate gradients, as reghdfe and ppmlhdfe do) and singleton dropping, when tables grow beyond a few hundred thousand cells.
3. **Weak-instrument-robust inference** (Anderson–Rubin, the Cragg–Donald and Kleibergen–Paap statistics with Stock and Yogo's tables) and LIML.
4. **Multinomial and conditional logit**, ordered models and their marginal effects; staggered differences in differences (Callaway and Sant'Anna 2021).
