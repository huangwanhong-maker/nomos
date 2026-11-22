# games: Finite Games in Normal Form

*Status: draft (2026-10-25), [decision 0003](../development/decisions/0003-the-first-modules.md); routing games with congestion (2026-11-12, [decision 0010](../development/decisions/0010-spatial-interaction-and-congestion.md)); voting power, bargaining and the nucleolus (2026-11-19, [decision 0014](../development/decisions/0014-voting-power-bargaining-and-the-nucleolus.md)).* Extensive forms, Lemke–Howson, evolutionary dynamics and mechanism design are planned (the last in `mechanisms`).

## 1. Concepts

- **Normal form.** Each player chooses a strategy, and every profile of choices gives every player a payoff (von Neumann and Morgenstern 1944). A two-player game is a bimatrix (A, B): the row player receives A_ij and the column player B_ij. Mixed strategies are probability vectors x and y, with expected payoffs xᵀAy and xᵀBy.
- **Best responses and equilibrium.** A Nash equilibrium is a profile in which each strategy is a best response to the others (Nash 1951). In a bimatrix game, (x, y) is one exactly when every pure strategy in the support of each is a best response to the other: the best response condition (von Stengel 2007, Proposition 3.1).
- **Dominance.** A strategy strictly dominated by another, pure or mixed, is never a best response and is used in no equilibrium; eliminating such strategies, for both players in turn, leaves the equilibria unchanged, and for strict dominance the order of elimination does not matter (Gilboa, Kalai and Zemel 1990). Dominance by a mixture is a linear programme.
- **Support enumeration.** In a nondegenerate game every equilibrium has supports of equal size k, on which the strategies solve square linear systems of indifference; enumerating the pairs of supports finds them all (von Stengel 2007, §3.2), and their number is odd (Lemke and Howson 1964; Wilson 1971).
- **Zero-sum games.** When B = −A, the maximin and minimax values meet (von Neumann 1928); the value and the optimal strategies are a pair of dual linear programmes (Dantzig 1951).
- **Many players.** Pure equilibria of an n-player game are found by checking every profile against every unilateral deviation; one need not exist.
- **Coalitional games and claims.** A coalitional game gives a worth to every set of players; the Shapley value (1953) gives each the mean of what it adds over all orders of arrival. A claims problem is an estate short of the claims on it (O'Neill 1982): paid in an order, the order decides; the mean over orders is the Shapley value of what each set secures after the others are paid; paid together, the rule depends on the rates (in proportion to the claims, or equal); the Talmud's rule is the one consistent with the contested-garment division of every pair (Aumann and Maschler 1985).
- **Routing games with congestion.** A continuum of travellers between origins and destinations each choose a path, every arc's time rising with its flow. At Wardrop's equilibrium (1952) every path in use takes the least time of its pair. The equilibrium arc flows minimise Beckmann's potential Σ_a ∫ t_a (Beckmann, McGuire and Winsten 1956), Rosenthal's potential made continuous, so they exist and, for strictly increasing times, are unique. The flows least in total time are the equilibrium of the same game at marginal costs. Selfish routing costs at most 4/3 of the optimum when times are affine (Roughgarden and Tardos 2002), and a road added can lengthen every journey (Braess 1968). Frank–Wolfe solves both: shortest paths at the current times, a step by exact line search, and a lower bound from convexity.
- **Simple games and a priori voting power.** A decision rule is a simple game: which coalitions win. What the rule alone gives a member, before any preference is known, is its a priori power, a measure of the rule and never of influence (Felsenthal and Machover 1998). A swing of i is a losing coalition that i makes winning. The four measures:
  - Penrose's measure (1946): the probability that i is decisive when the others vote yes or no with equal chance, its swings over 2^(n−1);
  - Banzhaf's index (1965): the swings normalised;
  - Shapley and Shubik's index (1954): the Shapley value of the game;
  - Coleman's measures (1971): the powers to prevent and to initiate action.

  In a two-tier system, council weights proportional to the square root of each unit's size make citizens' indirect power nearly equal (Penrose 1946).
- **Bargaining.** Two parties choose a point of a convex set of utilities, with what each has if they do not agree. The solutions:
  - Nash's (1950) maximises the product of the gains. It is the one solution that is Pareto optimal, symmetric, invariant under changes of scale and independent of irrelevant alternatives.
  - Kalai and Smorodinsky's (1975) gives up the last for individual monotonicity.
  - The egalitarian solution (Kalai 1977) equalises gains, and needs utilities comparable between the parties.

  Brams and Taylor's Adjusted Winner (1996) divides issues between two parties. Their points come out equal, neither envies the other, and no other division is better for one and no worse for the other.
- **The core and the nucleolus.** The core of a coalitional game holds the divisions no coalition can improve on by leaving (Gillies 1959); it is non-empty exactly when the game is balanced (Bondareva 1963; Shapley 1967). The nucleolus makes the largest dissatisfactions, lexicographically, as small as possible (Schmeidler 1969); it lies in the core when the core is not empty, and for a claims problem it is the Talmud's rule (Aumann and Maschler 1985).

## 2. Content

| Header | Content |
|---|---|
| `bimatrix.hpp` | `Bimatrix{a, b}` with `zero_sum(A)` and `symmetric(A)`; `row_payoff(g, x, y)`, `column_payoff`; `best_responses_row(g, y)`, `best_responses_column(g, x)`; `Cell{row, col}`; `pure_equilibria(g)` |
| `dominance.hpp` | `enum class Dominance {pure, mixed}`; `eliminate_dominated(g, kind)` → `Reduction{rows, cols, removed}` |
| `equilibria.hpp` | `Equilibrium{x, y, row_value, column_value}`; `is_equilibrium(g, x, y)`; `support_enumeration(g)` |
| `zero_sum.hpp` | `solve_zero_sum(A)` → `ZeroSumSolution{lower_value, upper_value, value, x, y}` |
| `normal_form.hpp` | `NormalForm{strategies, payoffs}` with `index`, `profile`, `payoff`; `normal_form(Bimatrix)`; `pure_equilibria(NormalForm)` |
| `claims.hpp` | `shapley(n, v)`, the Shapley value of a coalitional game given on sets as masks; `Claims{estate, claims}`; `priority(P, order)` and `priority_range(P)`; `random_arrival(P)`; `proportional`, `constrained_equal_awards`, `constrained_equal_losses`, `talmud`; `contested_garment(E, c₁, c₂)`; `interleaved(P, order, rates, step)` ([decision 0011](../development/decisions/0011-claims-and-the-shapley-value.md)) |
| `congestion.hpp` | `Arc{tail, head, a, b, power}` (t = a + b x^p) and `bpr(…)`; `Network{nodes, arcs, demands}`; `assign(network, {objective, gap, …})` → `Assignment{flow, time, total_time, objective, lower_bound, gap, …}`; `relative_gap`, `skims`, `times`, `beckmann`, `total_time`, `price_of_anarchy` |
| `power.hpp` | `SimpleGame{players, wins}`, `WeightedGame{weights, quota}`, `simple(g)`, `all_of`, `any_of`; `properties(g)` → `Properties{monotone, proper, strong, decisive(), vetoes, dummies, minimal_winning}`; `power(g)` → `Power{swings, penrose, banzhaf, shapley_shubik, coleman_prevent, coleman_initiate, power_to_act}`, by enumeration (n ≤ 24) or, for a weighted game, by generating functions (any n); `decisive_in_majority(m)`, `indirect_power(sizes, council)` ([decision 0014](../development/decisions/0014-voting-power-bargaining-and-the-nucleolus.md)) |
| `bargaining.hpp` | `Bargaining{outcomes, disagreement}`; `nash`, `kalai_smorodinsky`, `egalitarian` → `BargainingSolution{point, ideal, frontier}`; `adjusted_winner(a, b)` → `Division{share, points_a, points_b, split}` |
| `dynamics.hpp` | Population games and evolutionary dynamics ([decision 0017](../development/decisions/0017-population-games-and-evolutionary-dynamics.md)): `linear(A)`, `rock_paper_scissors(w, l)`; `Protocol {replicator, best_response, logit, bnn, smith, maynard_smith}`; `protocol(p, π, x)` (the matrix ρ), `mean_dynamic(ρ, x)` (eq. 13.3), `field(p, F, x)` (Table 13.1's closed forms); `trajectory(p, F, x0, T, h)`; `incentive(p, F, x)` (positive correlation, Nash stationarity), `is_nash`; `replicator_entropy`, `bnn_lyapunov`, `smith_lyapunov`, `best_response_lyapunov` (Table 13.4); `finite_population(p, F, counts, times, stream)`, exact |
| `nucleolus.hpp` | `least_core(n, v)` → `LeastCore{epsilon, x, core_nonempty()}`; `in_core(x, v)`; `excesses(x, v)`; `prenucleolus(n, v)` → `Nucleolus{x, levels, rounds, imputation}` (n ≤ 10) |

Payoff matrices are `calx::linalg::Rectangular`. Ties within `tolerance` (default 10⁻⁹) count as ties.

## 3. Checks

`tests/check_games.cpp` (30 checks): the prisoner's dilemma (Gibbons 1992, §1.1.A), with one pure equilibrium reached by dominance; Gibbons' dominance-solvable game (§1.1.B) and the order of its eliminations; matching pennies (no pure equilibrium; (½, ½); value 0 by linear programming); the battle of the sexes (three equilibria, the mixed one paying ⅔); rock, paper, scissors (value 0, uniform); a 2 × 2 zero-sum game of value ⅕; von Stengel's 3 × 2 game with three equilibria; a row dominated only by a mixture; forty random nondegenerate games (the best response condition, an odd number of equilibria, the pure ones among them); ten random zero-sum games (the minimax theorem, the guarantees, agreement with support enumeration); a three-player coordination game, Jordan's (1993) game without a pure equilibrium, and the numbering of profiles.

`tests/check_claims.cpp` (28 checks): the glove game; efficiency, symmetry and additivity of the value; the Talmud's three estates and the contested garment; consistency on every pair; random arrival against every order by brute force; every rule efficient, bounded and order-preserving; payment by turns, which tends to the proportional rule at rates in proportion to the claims and to the equal awards at equal rates.

`tests/check_power.cpp` (26 checks):

- the UN Security Council's rule against swings counted in closed form (848 and 84; Shapley–Shubik's 0.0018648 for a non-permanent member), and the same rule as a weighted game, by generating functions;
- Nassau County's board of 1964 and its three dummies (Banzhaf 1965);
- the two routes agreeing on 300 random weighted games;
- 101 and 301 equal members;
- properties of games, a double majority with a veto among them;
- the square-root rule on ten units of 1 to 100 million: indirect power within 10% under weights ∝ √m, and the smallest unit a dummy under weights ∝ m.

`tests/check_bargaining.cpp` (20 checks):

- closed forms;
- affine invariance (Nash, Kalai–Smorodinsky) and its failure (egalitarian);
- independence of irrelevant alternatives (Nash) and its failure (Kalai–Smorodinsky);
- individual monotonicity (Kalai–Smorodinsky) and its failure (Nash);
- Nash's product maximal on 200 random problems;
- Adjusted Winner worked by hand, and on 300 random problems equitable, envy-free, splitting at most one issue, and efficient against a linear programme.

`tests/check_nucleolus.cpp` (11 checks):

- the glove game;
- the majority game of three: an empty core, the least core's value 1/3;
- the nucleolus of a claims game equal to the Talmud rule, on Aumann and Maschler's estate of 200 and on 150 random problems;
- convex games, whose core holds the nucleolus and the Shapley value;
- no nearby pre-imputation lexicographically better;
- symmetric players alike.

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

- **Evolutionary dynamics** (`check_nomos_dynamics`, 22): each protocol's mean dynamic its closed form; the replicator's closed orbits in standard Rock–Paper–Scissors; convergence with the Lyapunov functions never rising in good Rock–Paper–Scissors; the heteroclinic approach in bad; positive correlation on random games; Nash stationarity failing for the replicator at a vertex; a potential that never falls, and BNN's algebraic approach to a strict equilibrium; Hawk–Dove's V/C; a dominated twin surviving under Smith (Hofbauer and Sandholm 2011); the logit's rest point; the finite population within 0.005 of the mean dynamic at N = 20,000.

## 4. Dependencies

calx `linalg` (`Rectangular`, `solve`), `optimisation` (the simplex method) and `special_functions` (log-binomials, for the square-root rule); morphe `combinatorics` (Dijkstra's shortest paths, for routing games). In the checks, calx `random` and `verification`.

## 5. References

- Aumann, R. J. and Maschler, M. (1985), Game theoretic analysis of a bankruptcy problem from the Talmud, *J. Econ. Theory* 36, 195–213.
- Banzhaf, J. F. (1965), Weighted voting doesn't work: a mathematical analysis, *Rutgers Law Review* 19, 317–343.
- Bondareva, O. N. (1963), Some applications of linear programming methods to the theory of cooperative games, *Problemy Kibernetiki* 10, 119–139.
- Brams, S. J. and Taylor, A. D. (1996), *Fair Division: From Cake-Cutting to Dispute Resolution*, Cambridge University Press.
- Coleman, J. S. (1971), Control of collectivities and the power of a collectivity to act, in B. Lieberman (ed.), *Social Choice*, Gordon and Breach.
- Felsenthal, D. S. and Machover, M. (1998), *The Measurement of Voting Power*, Edward Elgar.
- Gillies, D. B. (1959), Solutions to general non-zero-sum games, *Contributions to the Theory of Games* IV, 47–85.
- Kalai, E. (1977), Proportional solutions to bargaining situations: interpersonal utility comparisons, *Econometrica* 45, 1623–1630.
- Kalai, E. and Smorodinsky, M. (1975), Other solutions to Nash's bargaining problem, *Econometrica* 43, 513–518.
- Kopelowitz, A. (1967), Computation of the kernels of simple games and the nucleolus of n-person games, Research Memorandum 31, Hebrew University.
- Maschler, M., Peleg, B. and Shapley, L. S. (1979), Geometric properties of the kernel, nucleolus and related solution concepts, *Math. Oper. Res.* 4, 303–338.
- Nash, J. F. (1950), The bargaining problem, *Econometrica* 18, 155–162.
- Penrose, L. S. (1946), The elementary statistics of majority voting, *J. Royal Statistical Soc.* 109, 53–57.
- Schmeidler, D. (1969), The nucleolus of a characteristic function game, *SIAM J. Appl. Math.* 17, 1163–1170.
- Shapley, L. S. (1967), On balanced sets and cores, *Naval Res. Logist. Q.* 14, 453–460.
- Shapley, L. S. and Shubik, M. (1954), A method for evaluating the distribution of power in a committee system, *Amer. Polit. Sci. Rev.* 48, 787–792.
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
- Sandholm, W. H. (2015), Population games and deterministic evolutionary dynamics, in *Handbook of Game Theory with Economic Applications* 4, ch. 13, 703–778 (doi 10.1016/B978-0-444-53766-9.00013-6).
- Hofbauer, J. and Sandholm, W. H. (2011), Survival of dominated strategies under evolutionary dynamics, *Theoretical Economics* 6, 341–377.
- Taylor, P. D. and Jonker, L. B. (1978), Evolutionarily stable strategies and game dynamics, *Mathematical Biosciences* 40, 145–156.
- Smith, M. J. (1984), The stability of a dynamic model of traffic assignment, *Transportation Science* 18, 245–252.

## 6. Open Questions

1. **Exact equilibria** on morphe's rationals (support enumeration and Lemke–Howson by integer pivoting), for degenerate games and for certainty.
2. **Extensive forms**, the sequence form and subgame perfection, when a domain models sequential interaction.
3. **Faster routing equilibria.** Frank–Wolfe reaches a gap of 10⁻⁴ quickly and 10⁻⁸ slowly. Path-based gradient projection (Jayakrishnan et al. 1994) or Bar-Gera's origin-based algorithm would reach tight gaps. Conjugate Frank–Wolfe (Mitradjieva and Lindberg 2013) was tried and stalled on the checks' grid.
