# 0014. Voting Power, Bargaining and the Nucleolus

- **Status:** Accepted, provisionally (2026-11-19).
- **Date:** 2026-11-19
- **Decided by:** the infrastructure's ADR 0132, which founds its diplomacy domain and asks nomos for the a priori power of decision rules, axiomatic bargaining and the nucleolus (its survey's §A.11.2); the nucleolus and the core were named as next in decision 0011.

## 1. Context

A multilateral body decides by a rule: a weighted vote with a quota, a double majority, vetoes, a rule deemed met when too few would block it. What the rule alone gives each member, before any preference is known, is its a priori power, a property of the rule and not of influence. Negotiations between parties are modelled by bargaining solutions over the utilities they can reach together, and the division of a joint gain by the nucleolus or the core. All of them name no society and are established in the social sciences, so they are nomos's (decision 0002). `games` held the Shapley value and claims problems (decision 0011), but no simple games, no bargaining and no nucleolus.

## 2. Decision

`games` gains three headers.

1. **`power.hpp`: simple games and a priori voting power.**
   - Simple games given as predicates on coalitions (n ≤ 24), weighted games with a quota, and compound games (`all_of`, `any_of`: a double majority, a rule deemed met).
   - Their properties: monotone, proper, strong, decisive; vetoes, dummies and minimal winning coalitions.
   - The measures: swings; the Penrose measure (Banzhaf's absolute index), the probability of being decisive under independent equiprobable votes; Banzhaf's normalised index; Shapley–Shubik's index, through `shapley`; Coleman's powers to prevent and to initiate action; the collectivity's power to act.
   - Two routes: any simple game on n ≤ 24 players exactly over its coalitions; a weighted game on any number of players by generating functions, pseudo-polynomial in the quota, its counts kept as probabilities so that nothing overflows.
   - The square-root rule: `decisive_in_majority(m)` = C(m − 1, ⌊m/2⌋)/2^(m−1), and `indirect_power`, a citizen's power through their unit's majority and their delegate's council (Penrose 1946).
2. **`bargaining.hpp`: two-party bargaining and fair division.**
   - Over the convex hull of finitely many pairs of utilities with a disagreement point: the Pareto frontier at or above it, the ideal point, and the Nash, Kalai–Smorodinsky and egalitarian solutions, each in closed form on the frontier's segments.
   - Brams and Taylor's Adjusted Winner for two parties who each spread 100 points over the issues.
3. **`nucleolus.hpp`: the core, the least core and the nucleolus**, for n ≤ 10, on calx's simplex.
   - `least_core` (the core is non-empty exactly when its value is at most 0), `in_core`, `excesses`.
   - `prenucleolus`, by Kopelowitz's sequence of programmes: the coalitions with positive duals fixed at each level, those whose excess the fixed ones determine dropped, until the fixed coalitions span ℝⁿ. It says whether the point is an imputation, in which case it is the nucleolus.

## 3. Verification

- **`check_nomos_power` (26 checks).**
  - The UN Security Council's rule, against swings counted in closed form: 848 and 84, Banzhaf's 848/5080 and 84/5080, and Shapley–Shubik's C(9, 3) 8! 6!/15! = 0.0018648 for a non-permanent member. The same rule as the weighted game [39; 7 × 5, 1 × 10] gives the same by generating functions.
  - Nassau County's board of 1964 with its three dummies (Banzhaf 1965).
  - The two routes on 300 random weighted games, agreeing to 10⁻¹⁴.
  - 101 and 301 equal members, against C(n − 1, ⌊n/2⌋)/2^(n−1) and 1/n.
  - Properties of games, among them a double majority with a veto.
  - The square-root rule: ten units of 1 to 100 million. With council weights ∝ √m, citizens' indirect power lies within 10% (max/min 1.10). With weights ∝ m the smallest unit is a dummy, and its citizens have no power at all.
- **`check_nomos_bargaining` (20 checks).**
  - Closed forms on a triangle, a linear frontier and a disagreement point inside.
  - Invariance under positive affine changes for Nash and Kalai–Smorodinsky, not for the egalitarian solution.
  - Independence of irrelevant alternatives for Nash, not for Kalai–Smorodinsky.
  - Individual monotonicity for Kalai–Smorodinsky, not for Nash, on a problem of the kind Kalai and Smorodinsky gave.
  - On 200 random problems, Nash's product maximal over sampled points of the frontier.
  - Adjusted Winner worked by hand. On 300 random problems it is equitable, envy-free and splits at most one issue, and it is efficient, held to a linear programme.
- **`check_nomos_nucleolus` (11 checks).**
  - The glove game and the majority game of three.
  - The nucleolus of a claims game equals the Talmud rule (Aumann and Maschler 1985): on their estate of 200, and on 150 random problems to 10⁻¹⁵.
  - Convex games, whose core holds the nucleolus and the Shapley value.
  - No nearby pre-imputation lexicographically better, in 12,000 perturbations of 40 random games.
  - Symmetric players treated alike.

nomos's checks rise from 408 to 465.

## 4. Alternatives Considered

- **Power indices in the domain.** Rejected: a decision rule's power names no body; the diplomacy domain declares the rules of a body and asks nomos.
- **The nucleolus by enumerating the kernel or by Kohlberg's criterion.** Kopelowitz's sequence is the standard and reuses calx's simplex; the checks hold it to the Talmud rule, which Aumann and Maschler proved to be the nucleolus.
- **Bargaining among n parties.** Deferred: two parties carry the axioms and the cases the survey names. More parties need a convex programme over a general set, or the transferable case, which a caller can state as a coalitional game.

## 5. Consequences

- Large weighted games with large weights (votes in the hundreds of thousands) are the caller's to round, the generating functions being pseudo-polynomial; sampling over orders, and the power indices of games too large for both routes, can follow.
- The kernel, the bargaining set, the egalitarian solution for n parties and Rubinstein's alternating offers are left until a module needs them.
