# 0002. Scope, Boundaries and Dependencies

- **Status:** Accepted
- **Date:** 2026-10-25
- **Decided by:** the user, who approved a foundation for the methods of the social sciences under the name nomos; Claude drew the boundaries, following the infrastructure's map of the mathematics (`docs/mathematics.md`, §2 and §7) and its ADR 0056.

## 1. Context

The infrastructure's map of the mathematics planned, beneath the domains and above the inference foundation, a library for the methods of the quantitative social sciences: measurement, surveys, econometrics, games and social choice, social networks and agent-based models. Whether they should form a library of their own or modules of the inference foundation was open. The economics domain needs them now: it measures concentration and specialisation of flows, and will estimate gravity equations and measure inequality. The inference foundation is being written at the same time as tyche.

## 2. Decision

1. **A library of its own.** The methods of the social sciences are the library nomos (Greek *νόμος*, law, custom), namespace `nomos`, one target `nomos::<module>` for each module, in `foundations/nomos`.
2. **The boundary with tyche.** tyche holds general probability, descriptive statistics, estimation and testing, regression (linear and generalised linear models with their standard errors), time series, and (planned) causal graphs and potential outcomes in general. nomos holds the methods of the social sciences built on them. The rule: a general statistical model or procedure, used across the sciences, is tyche's; a method whose formulation is that of agents, choices, institutions, welfare and distribution, or whose established use is in the social sciences, is nomos's. Examples:
   - ordinary least squares with heteroskedasticity-robust standard errors, logit and probit fits, quantiles: tyche;
   - fixed-effects and random-effects panel estimators, two-stage least squares, differences in differences, average marginal effects, pseudo-maximum likelihood with absorbed fixed effects: nomos `econometrics`, on tyche's regressions;
   - the Gini coefficient, the Lorenz curve, poverty and concentration indices, index numbers: nomos;
   - games, matching and social choice: nomos, on calx's optimisation and, later, morphe's graphs and orders.
3. **Dependencies.** nomos depends on tyche, calx, morphe and physis, and never on the infrastructure's domains, edges or works. It does not nest them as submodules: inside the infrastructure it uses the siblings under `foundations/`; built alone, it would fetch them.
4. **No particular objects.** nomos names no country, currency, commodity, election or census. A distribution is a list of amounts with weights; a two-way table has rows and columns; a market has applicants and institutions. The domains supply the meaning.
5. **Language.** C++20, header-only while the modules are templates and small functions, following the infrastructure's ADR 0003 and ADR 0035.

## 3. Alternatives Considered

- **Modules of tyche.** Rejected. Games, matching and social choice are not statistics, and measurement and index numbers are axiomatic before they are statistical; putting them in tyche would make the inference foundation a library of everything quantitative.
- **Methods in the economics domain.** Rejected. Inequality, concentration, games and choice serve sociology, political science, demography and history as much as economics, and none names an object in nature; the infrastructure's principle is that such methods are written once, in a foundation.
- **A library per discipline** (econometrics, psychometrics, game theory). Rejected as premature; modules can be extracted later under the graduation rule of the infrastructure's ADR 0026.

## 4. Consequences

- A change that spans libraries is committed in dependency order: calx and morphe, then physis, then tyche, then nomos, then the infrastructure.
- The domains of the social sciences apply nomos's methods to particular societies; a method one of them needs that names no social object goes to tyche or calx instead.
