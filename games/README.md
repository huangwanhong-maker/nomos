# games: Finite Games in Normal Form

*Status: draft (2026-10-25), [decision 0003](../development/decisions/0003-the-first-modules.md).* Extensive forms, Lemke–Howson, evolutionary dynamics and mechanism design are planned (the last in `mechanisms`).

## 1. Concepts

- **Normal form.** Each player chooses a strategy, and every profile of choices gives every player a payoff (von Neumann and Morgenstern 1944). A two-player game is a bimatrix (A, B): the row player receives A_ij and the column player B_ij. Mixed strategies are probability vectors x and y, with expected payoffs xᵀAy and xᵀBy.
- **Best responses and equilibrium.** A Nash equilibrium is a profile in which each strategy is a best response to the others (Nash 1951). In a bimatrix game, (x, y) is one exactly when every pure strategy in the support of each is a best response to the other: the best response condition (von Stengel 2007, Proposition 3.1).
- **Dominance.** A strategy strictly dominated by another, pure or mixed, is never a best response and is used in no equilibrium; eliminating such strategies, for both players in turn, leaves the equilibria unchanged, and for strict dominance the order of elimination does not matter (Gilboa, Kalai and Zemel 1990). Dominance by a mixture is a linear programme.
- **Support enumeration.** In a nondegenerate game every equilibrium has supports of equal size k, on which the strategies solve square linear systems of indifference; enumerating the pairs of supports finds them all (von Stengel 2007, §3.2), and their number is odd (Lemke and Howson 1964; Wilson 1971).
- **Zero-sum games.** When B = −A, the maximin and minimax values meet (von Neumann 1928); the value and the optimal strategies are a pair of dual linear programmes (Dantzig 1951).
- **Many players.** Pure equilibria of an n-player game are found by checking every profile against every unilateral deviation; one need not exist.

## 2. Content

| Header | Content |
|---|---|
| `bimatrix.hpp` | `Bimatrix{a, b}` with `zero_sum(A)` and `symmetric(A)`; `row_payoff(g, x, y)`, `column_payoff`; `best_responses_row(g, y)`, `best_responses_column(g, x)`; `Cell{row, col}`; `pure_equilibria(g)` |
| `dominance.hpp` | `enum class Dominance {pure, mixed}`; `eliminate_dominated(g, kind)` → `Reduction{rows, cols, removed}` |
| `equilibria.hpp` | `Equilibrium{x, y, row_value, column_value}`; `is_equilibrium(g, x, y)`; `support_enumeration(g)` |
| `zero_sum.hpp` | `solve_zero_sum(A)` → `ZeroSumSolution{lower_value, upper_value, value, x, y}` |
| `normal_form.hpp` | `NormalForm{strategies, payoffs}` with `index`, `profile`, `payoff`; `normal_form(Bimatrix)`; `pure_equilibria(NormalForm)` |

Payoff matrices are `calx::linalg::Rectangular`. Ties within `tolerance` (default 10⁻⁹) count as ties.

## 3. Checks

`tests/check_games.cpp` (30 checks): the prisoner's dilemma (Gibbons 1992, §1.1.A), with one pure equilibrium reached by dominance; Gibbons' dominance-solvable game (§1.1.B) and the order of its eliminations; matching pennies (no pure equilibrium; (½, ½); value 0 by linear programming); the battle of the sexes (three equilibria, the mixed one paying ⅔); rock, paper, scissors (value 0, uniform); a 2 × 2 zero-sum game of value ⅕; von Stengel's 3 × 2 game with three equilibria; a row dominated only by a mixture; forty random nondegenerate games (the best response condition, an odd number of equilibria, the pure ones among them); ten random zero-sum games (the minimax theorem, the guarantees, agreement with support enumeration); a three-player coordination game, Jordan's (1993) game without a pure equilibrium, and the numbering of profiles.

## 4. Dependencies

calx `linalg` (`Rectangular`, `solve`) and `optimisation` (the simplex method). In the checks, calx `random` and `verification`.

## 5. References

- Dantzig, G. B. (1951), A proof of the equivalence of the programming problem and the game problem, in T. C. Koopmans (ed.), *Activity Analysis of Production and Allocation*, Wiley.
- Chvátal, V. (1983), *Linear Programming*, Freeman.
- Dickhaut, J. and Kaplan, T. (1991), A program for finding Nash equilibria, *Mathematica Journal* 1, 87–93.
- Gibbons, R. (1992), *Game Theory for Applied Economists*, Princeton University Press.
- Gilboa, I., Kalai, E. and Zemel, E. (1990), On the order of eliminating dominated strategies, *Oper. Res. Lett.* 9, 85–89.
- Jordan, J. S. (1993), Three problems in learning mixed-strategy Nash equilibria, *Games Econ. Behav.* 5, 368–386.
- Lemke, C. E. and Howson, J. T. (1964), Equilibrium points of bimatrix games, *J. SIAM* 12, 413–423.
- Monderer, D. and Shapley, L. S. (1996), Potential games, *Games Econ. Behav.* 14, 124–143.
- Nash, J. (1951), Non-cooperative games, *Ann. Math.* 54, 286–295.
- Osborne, M. J. and Rubinstein, A. (1994), *A Course in Game Theory*, MIT Press.
- Pearce, D. G. (1984), Rationalizable strategic behavior and the problem of perfection, *Econometrica* 52, 1029–1050.
- Shoham, Y. and Leyton-Brown, K. (2009), *Multiagent Systems*, Cambridge University Press.
- von Neumann, J. (1928), Zur Theorie der Gesellschaftsspiele, *Math. Ann.* 100, 295–320.
- von Neumann, J. and Morgenstern, O. (1944), *Theory of Games and Economic Behavior*, Princeton University Press.
- von Stengel, B. (2007), Equilibrium computation for two-player games in strategic and extensive form, in N. Nisan et al. (eds), *Algorithmic Game Theory*, Cambridge University Press, ch. 3.
- Wilson, R. (1971), Computing equilibria of n-person games, *SIAM J. Appl. Math.* 21, 80–87.

## 6. Open Questions

1. **Exact equilibria** on morphe's rationals (support enumeration and Lemke–Howson by integer pivoting), for degenerate games and for certainty.
2. **Extensive forms**, the sequence form and subgame perfection, when a domain models sequential interaction.
