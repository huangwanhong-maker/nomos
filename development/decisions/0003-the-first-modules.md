# 0003. The First Modules

- **Status:** Accepted, provisionally
- **Date:** 2026-10-25
- **Decided by:** the user, who asked for the first modules of nomos and, through the infrastructure's coordinating session, for Poisson pseudo-maximum likelihood with two-way fixed effects; Claude chose their content, layout and checks.

## 1. Context

The economics domain measures the concentration and specialisation of flows between places and will estimate gravity equations, measure inequality and compare prices over time; the history domain and the planned domains of political science and sociology need games, choice and matching. None of these names an object in nature, so by the infrastructure's principle of reuse they belong in a foundation, nomos ([decision 0002](0002-scope-boundaries-and-dependencies.md)). tyche's `regression` (ordinary and weighted least squares with robust and clustered covariances, Wald tests, generalised linear models) and calx's simplex method were written at the same time, and the econometrics here is built on them.

## 2. Decision

### 2.1 The Modules

Six modules, each a header-only target `nomos::<module>` with its checks:

| Module | Content |
|---|---|
| `measurement` | inequality (Gini, Lorenz, generalised entropy, Atkinson, Palma, quantile shares, CV), poverty (FGT, Watts), concentration (HHI, CR_k, entropy, Hannah–Kay), specialisation (location quotients and Balassa's index, dissimilarity, Grubel–Lloyd), decompositions by groups |
| `indices` | bilateral price and quantity indices (Laspeyres, Paasche, Fisher, Törnqvist, Walsh, Marshall–Edgeworth, Lowe, Young), elementary indices, direct and chained series |
| `games` | bimatrix and n-player normal forms, best responses, pure equilibria, iterated dominance (pure and mixed, by linear programming), support enumeration, zero-sum games by linear programming |
| `matching` | deferred acceptance (one-to-one, many-to-one, either side proposing), stability, top trading cycles |
| `choice` | profiles and the pairwise majority matrix, plurality, Borda, approval, instant runoff, Condorcet winners and losers, Copeland, Schulze, Kemeny–Young |
| `econometrics` | absorbed fixed effects by alternating projections; within (one- and two-way), between, random effects, Hausman; two-stage least squares with the first stage, Sargan and Hansen; differences in differences (2 × 2 and two-way fixed effects); average marginal effects of logit and probit; Poisson pseudo-maximum likelihood with two-way fixed effects |

Planned: `psychometrics`, `surveys`, `networks`, `mechanisms`, `agents`.

### 2.2 Their Order

| Module | Depends on |
|---|---|
| `measurement` | calx `linalg` |
| `indices` | — |
| `games` | calx `linalg`, `optimisation` |
| `matching` | — |
| `choice` | calx `linalg` |
| `econometrics` | tyche `regression`, `probability`; calx `linalg` |

No module of nomos uses another yet; the root `CMakeLists.txt` adds them in the order of the table.

### 2.3 Conventions

- Weights are frequency or sampling weights, and indices of distributions use population divisors.
- Preferences and ballots are lists, best first; unlisted partners are unacceptable and unranked candidates tied below the ranked.
- Payoffs and tables are `calx::linalg::Rectangular`; square matrices (the pairwise majority matrix, covariances) `calx::linalg::Matrix`.
- Econometric fits return tyche's `LinearFit` where the estimator is a least-squares fit (within, between, random effects), corrected where the degrees of freedom differ, so that tyche's Wald tests apply to them.
- Equation, theorem and section numbers of sources are cited only where they are certain; otherwise the work, or its chapter.

### 2.4 Checks

Each module is checked against published examples where their numbers are certain (the CPI Manual's artificial data set, Gale and Shapley's example, Schulze's example, Gibbons' and von Stengel's games), closed forms (the exponential distribution's inequality indices, the independence model of a two-way table, the simple IV estimator), brute force (every matching, every coalition, every ranking, every false list) and the axioms each method satisfies or is known to fail; the econometric estimators against the equivalent dense regressions of tyche (dummy variables for absorbed effects) and against simulated truths within their standard errors.

## 3. Alternatives Considered

- **Absorbing fixed effects in tyche.** The alternating projections name no social object; they are here because econometrics is their established home and tyche's `regression` was being written at the same time. If another science needs them, they move down to tyche (the module's open question 1).
- **Dense dummies for the two-way Poisson model.** Rejected: a table of 200 × 200 cells would need 400 dummy columns and a QR of a 40,000 × 400 matrix in every iteration, where the projections cost a few passes over the cells.
- **Exact rational arithmetic for games.** Deferred: morphe's rationals would make equilibria exact, at a cost in speed and code; floating point with stated tolerances serves the first uses.

## 4. Consequences

- The economics domain can measure inequality, concentration and specialisation, compare prices, and estimate gravity equations with nomos, writing none of these methods itself.
- tyche's interfaces (`ols`, `glm`, `wald`, `LinearFit`, `GlmFit`) are now used by a sibling; changing them means changing nomos too.
