# nomos

*The methods of the social sciences.* The name comes from the Greek *νόμος*, "law, custom", the *-nomy* of economy.

## 1. Purpose

nomos implements the general methods by which the social sciences measure, model and infer: how a distribution is unequal and a total concentrated, how prices are compared across time, how strategic agents reach equilibrium, how two sides of a market are matched, how individual preferences are aggregated into a collective choice, and how the effects of choices and institutions are estimated from observational data. It names no object in nature: a distribution is a list of amounts with weights, a market a list of preferences, a panel a table of units and periods. The domains of The Science Infrastructures (economics first) apply its methods to particular economies, polities and societies.

It stands on the project's other foundations:

- [tyche](https://github.com/huangwanhong-maker/tyche), for probability, estimation, regression and time series;
- [calx](https://github.com/huangwanhong-maker/calx), for numerical linear algebra, optimisation (the simplex method) and random streams;
- [morphe](https://github.com/huangwanhong-maker/morphe), for the exact and combinatorial structures (graphs, orders) the planned modules will need;
- [physis](https://github.com/huangwanhong-maker/physis), through tyche.

nomos is one of the foundation libraries of [The Science Infrastructures](https://github.com/huangwanhong-maker/The-Science-Infrastructures). It stands on its own.

## 2. Modules

| Module | Content | Status |
|---|---|---|
| [`measurement`](measurement/) | Inequality: the Gini coefficient (weighted), the Lorenz curve and its generalisation, the generalised entropy family (Theil's T and L), Atkinson's indices, quantile shares, the Palma ratio, the coefficient of variation. Poverty: the Foster–Greer–Thorbecke family, the Watts index. Concentration: Herfindahl–Hirschman, concentration ratios, the entropy index, Hannah and Kay's numbers equivalents. Specialisation: location quotients (Balassa's index), the index of dissimilarity, the Grubel–Lloyd index and its adjusted aggregate. Decompositions by groups | draft; composite indicators ([decision 0006](development/decisions/0006-composite-indicators.md)) |
| [`indices`](indices/) | Index numbers: Laspeyres, Paasche, Fisher, Törnqvist, Walsh, Marshall–Edgeworth, Lowe and Young; the elementary indices of Jevons, Dutot and Carli, the harmonic and CSWD; quantity indices; direct and chained series | draft |
| [`games`](games/) | Finite games in normal form: payoffs and best responses, pure equilibria (two and n players), iterated elimination of strictly dominated strategies (by pure strategies and by mixtures), all equilibria of nondegenerate bimatrix games by support enumeration, zero-sum games by linear programming | draft |
| [`matching`](matching/) | Two-sided matching: deferred acceptance, one-to-one and many-to-one with capacities, either side proposing; individual rationality and blocking pairs. Housing markets: top trading cycles | draft |
| [`choice`](choice/) | Social choice: profiles of ranked ballots with weights, the pairwise majority matrix, plurality, Borda, approval, instant runoff; Condorcet winners and losers, Copeland, Schulze, Kemeny–Young | draft |
| [`econometrics`](econometrics/) | Average marginal effects of logit and probit, the multinomial logit, Poisson pseudo-maximum likelihood with two sets of absorbed fixed effects. Panels moved to tyche, instruments and differences in differences to peira ([decision 0007](development/decisions/0007-panels-instruments-and-differences-move.md)) | draft |
| `psychometrics` | Scales and their reliability (Cronbach's α), factor models of items, item response theory (Rasch, two- and three-parameter models) | planned |
| `surveys` | Sampling designs, design weights, calibration and raking, variance estimation (linearisation, replication) | planned |
| `networks` | Social networks on morphe's graphs: centrality, communities, homophily, diffusion | planned |
| `mechanisms` | Mechanism design: auctions, Vickrey–Clarke–Groves, implementation | planned |
| [`agents`](agents/) | Agent-based models: agents, schedules and interactions on calx's random streams | draft: synchronous (on a pool, the same on any number of threads), random sequential and Gillespie's continuous-time schedules, draws keyed by agent and tick, the event log, a rate tree, ensembles; lattices, graphs and the plane; logit choice; Schelling's model; the ODD+D description ([decision 0004](development/decisions/0004-agent-based-models.md)); Kirman's herding, with its exact law ([decision 0005](development/decisions/0005-kirmans-herding.md)) |

## 3. Principles

- **General methods, not particular societies.** nomos knows distributions, preferences, games and panels. It does not know a country, a currency, a census or an election; those are inputs from the domains.
- **Reuse before writing.** What tyche, calx or morphe hold is used, never written again: regressions are tyche's, the simplex method and linear algebra calx's.
- **Cited sources, checked results.** Every implemented formula cites its source; every result that can be checked numerically is checked, against published examples, closed forms, brute force and the axioms the method satisfies.

## 4. Relation to Sibling Libraries

nomos depends on tyche, calx, morphe and physis, the last two through tyche so far.

- **Inside the infrastructure**, it uses the copies already present under `foundations/`.
- **Built on its own**, it would fetch them through the build system (not yet written).

No library of the project depends on the infrastructure's domains or applications.

## 5. Status

**Draft (2026-10-25).** The first six modules are written and checked ([decision 0003](development/decisions/0003-the-first-modules.md)); the rest are planned. Inside The Science Infrastructures, build and check with CMake 3.20 or later (calx, morphe, physis and tyche are taken from the sibling directories):

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

or, from the infrastructure's root, `node development/tools/run-cpp.mjs foundations/nomos foundations/nomos`.

## 6. License

nomos is released under the MIT License. See [`LICENSE`](LICENSE).
