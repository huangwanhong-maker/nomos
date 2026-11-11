# games: Finite Games in Normal Form

*Status: draft (2026-10-25), [decision 0003](../development/decisions/0003-the-first-modules.md); routing games with congestion (2026-11-12, [decision 0010](../development/decisions/0010-spatial-interaction-and-congestion.md)).* Extensive forms, Lemke–Howson, evolutionary dynamics and mechanism design are planned (the last in `mechanisms`).

## 1. Concepts

- **Normal form.** Each player chooses a strategy, and every profile of choices gives every player a payoff (von Neumann and Morgenstern 1944). A two-player game is a bimatrix (A, B): the row player receives A_ij and the column player B_ij. Mixed strategies are probability vectors x and y, with expected payoffs xᵀAy and xᵀBy.
- **Best responses and equilibrium.** A Nash equilibrium is a profile in which each strategy is a best response to the others (Nash 1951). In a bimatrix game, (x, y) is one exactly when every pure strategy in the support of each is a best response to the other: the best response condition (von Stengel 2007, Proposition 3.1).
- **Dominance.** A strategy strictly dominated by another, pure or mixed, is never a best response and is used in no equilibrium; eliminating such strategies, for both players in turn, leaves the equilibria unchanged, and for strict dominance the order of elimination does not matter (Gilboa, Kalai and Zemel 1990). Dominance by a mixture is a linear programme.
- **Support enumeration.** In a nondegenerate game every equilibrium has supports of equal size k, on which the strategies solve square linear systems of indifference; enumerating the pairs of supports finds them all (von Stengel 2007, §3.2), and their number is odd (Lemke and Howson 1964; Wilson 1971).
- **Zero-sum games.** When B = −A, the maximin and minimax values meet (von Neumann 1928); the value and the optimal strategies are a pair of dual linear programmes (Dantzig 1951).
- **Many players.** Pure equilibria of an n-player game are found by checking every profile against every unilateral deviation; one need not exist.
- **Routing games with congestion.** A continuum of travellers between origins and destinations each choose a path, every arc's time rising with its flow. At Wardrop's equilibrium (1952) every path in use takes the least time of its pair. The equilibrium arc flows minimise Beckmann's potential Σ_a ∫ t_a (Beckmann, McGuire and Winsten 1956), Rosenthal's potential made continuous, so they exist and, for strictly increasing times, are unique. The flows least in total time are the equilibrium of the same game at marginal costs. Selfish routing costs at most 4/3 of the optimum when times are affine (Roughgarden and Tardos 2002), and a road added can lengthen every journey (Braess 1968). Frank–Wolfe solves both: shortest paths at the current times, a step by exact line search, and a lower bound from convexity.

## 2. Content

| Header | Content |
|---|---|
| `bimatrix.hpp` | `Bimatrix{a, b}` with `zero_sum(A)` and `symmetric(A)`; `row_payoff(g, x, y)`, `column_payoff`; `best_responses_row(g, y)`, `best_responses_column(g, x)`; `Cell{row, col}`; `pure_equilibria(g)` |
| `dominance.hpp` | `enum class Dominance {pure, mixed}`; `eliminate_dominated(g, kind)` → `Reduction{rows, cols, removed}` |
| `equilibria.hpp` | `Equilibrium{x, y, row_value, column_value}`; `is_equilibrium(g, x, y)`; `support_enumeration(g)` |
| `zero_sum.hpp` | `solve_zero_sum(A)` → `ZeroSumSolution{lower_value, upper_value, value, x, y}` |
| `normal_form.hpp` | `NormalForm{strategies, payoffs}` with `index`, `profile`, `payoff`; `normal_form(Bimatrix)`; `pure_equilibria(NormalForm)` |
| `congestion.hpp` | `Arc{tail, head, a, b, power}` (t = a + b x^p) and `bpr(…)`; `Network{nodes, arcs, demands}`; `assign(network, {objective, gap, …})` → `Assignment{flow, time, total_time, objective, lower_bound, gap, …}`; `relative_gap`, `skims`, `times`, `beckmann`, `total_time`, `price_of_anarchy` |

Payoff matrices are `calx::linalg::Rectangular`. Ties within `tolerance` (default 10⁻⁹) count as ties.

## 3. Checks

`tests/check_games.cpp` (30 checks): the prisoner's dilemma (Gibbons 1992, §1.1.A), with one pure equilibrium reached by dominance; Gibbons' dominance-solvable game (§1.1.B) and the order of its eliminations; matching pennies (no pure equilibrium; (½, ½); value 0 by linear programming); the battle of the sexes (three equilibria, the mixed one paying ⅔); rock, paper, scissors (value 0, uniform); a 2 × 2 zero-sum game of value ⅕; von Stengel's 3 × 2 game with three equilibria; a row dominated only by a mixture; forty random nondegenerate games (the best response condition, an odd number of equilibria, the pure ones among them); ten random zero-sum games (the minimax theorem, the guarantees, agreement with support enumeration); a three-player coordination game, Jordan's (1993) game without a pure equilibrium, and the numbering of profiles.

`tests/check_congestion.cpp` (23 checks):

- **Closed forms.**
  - Pigou's two roads: the equilibrium, the optimum, and the price of anarchy 4/3; with x⁴, 1/(1 − 4 · 5^(−5/4)).
  - Braess's network: 65 minutes, then 80 with the road added.
  - Two parallel roads of the Bureau of Public Roads, by bisection.
- **A grid of three by three crossings** against scipy's path-based solution, enumerated by networkx and polished by gradient projection (fixtures by `tests/fixtures/make_congestion_fixtures.py`):
  - scipy's flows have a gap of zero by the code here;
  - Frank–Wolfe at a gap of 10⁻⁵ lies within the distance its own bound allows;
  - the bound lies below the least potential;
  - the skim is the least time of every simple path.
- **Forty random affine networks:** the optimum at most the equilibrium, and their ratio at most 4/3.

## 4. Dependencies

calx `linalg` (`Rectangular`, `solve`) and `optimisation` (the simplex method); morphe `combinatorics` (Dijkstra's shortest paths, for routing games). In the checks, calx `random` and `verification`.

## 5. References

- Beckmann, M., McGuire, C. B. and Winsten, C. B. (1956), *Studies in the Economics of Transportation*, Yale University Press.
- Braess, D. (1968), Über ein Paradoxon aus der Verkehrsplanung, *Unternehmensforschung* 12, 258–268.
- Dantzig, G. B. (1951), A proof of the equivalence of the programming problem and the game problem, in T. C. Koopmans (ed.), *Activity Analysis of Production and Allocation*, Wiley.
- Chvátal, V. (1983), *Linear Programming*, Freeman.
- Dickhaut, J. and Kaplan, T. (1991), A program for finding Nash equilibria, *Mathematica Journal* 1, 87–93.
- Frank, M. and Wolfe, P. (1956), An algorithm for quadratic programming, *Naval Res. Logist. Q.* 3, 95–110.
- Gibbons, R. (1992), *Game Theory for Applied Economists*, Princeton University Press.
- Gilboa, I., Kalai, E. and Zemel, E. (1990), On the order of eliminating dominated strategies, *Oper. Res. Lett.* 9, 85–89.
- Jordan, J. S. (1993), Three problems in learning mixed-strategy Nash equilibria, *Games Econ. Behav.* 5, 368–386.
- LeBlanc, L. J., Morlok, E. K. and Pierskalla, W. P. (1975), An efficient approach to solving the road network equilibrium traffic assignment problem, *Transp. Res.* 9, 309–318.
- Lemke, C. E. and Howson, J. T. (1964), Equilibrium points of bimatrix games, *J. SIAM* 12, 413–423.
- Monderer, D. and Shapley, L. S. (1996), Potential games, *Games Econ. Behav.* 14, 124–143.
- Nash, J. (1951), Non-cooperative games, *Ann. Math.* 54, 286–295.
- Osborne, M. J. and Rubinstein, A. (1994), *A Course in Game Theory*, MIT Press.
- Pearce, D. G. (1984), Rationalizable strategic behavior and the problem of perfection, *Econometrica* 52, 1029–1050.
- Rosenthal, R. W. (1973), A class of games possessing pure-strategy Nash equilibria, *Int. J. Game Theory* 2, 65–67.
- Roughgarden, T. and Tardos, É. (2002), How bad is selfish routing?, *J. ACM* 49, 236–259.
- Shoham, Y. and Leyton-Brown, K. (2009), *Multiagent Systems*, Cambridge University Press.
- von Neumann, J. (1928), Zur Theorie der Gesellschaftsspiele, *Math. Ann.* 100, 295–320.
- von Neumann, J. and Morgenstern, O. (1944), *Theory of Games and Economic Behavior*, Princeton University Press.
- von Stengel, B. (2007), Equilibrium computation for two-player games in strategic and extensive form, in N. Nisan et al. (eds), *Algorithmic Game Theory*, Cambridge University Press, ch. 3.
- Wardrop, J. G. (1952), Some theoretical aspects of road traffic research, *Proc. Inst. Civ. Eng.* II 1, 325–378.
- Wilson, R. (1971), Computing equilibria of n-person games, *SIAM J. Appl. Math.* 21, 80–87.

## 6. Open Questions

1. **Exact equilibria** on morphe's rationals (support enumeration and Lemke–Howson by integer pivoting), for degenerate games and for certainty.
2. **Extensive forms**, the sequence form and subgame perfection, when a domain models sequential interaction.
3. **Faster routing equilibria.** Frank–Wolfe reaches a gap of 10⁻⁴ quickly and 10⁻⁸ slowly. Path-based gradient projection (Jayakrishnan et al. 1994) or Bar-Gera's origin-based algorithm would reach tight gaps. Conjugate Frank–Wolfe (Mitradjieva and Lindberg 2013) was tried and stalled on the checks' grid.
