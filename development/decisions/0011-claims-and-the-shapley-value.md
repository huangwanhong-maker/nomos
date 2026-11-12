# 0011. Claims on a Short Estate, and the Shapley Value

- **Status:** Accepted, provisionally (2026-11-12).
- **Decided by:** the infrastructure's ADR 0099 (the user's generative–relational market formalism, whose case on justice reads the order of procedural steps through the division of an estate short of its claims).

## 1. Decision

`games` gains `claims.hpp`:

- **the Shapley value** of a coalitional game on up to 24 players, exactly over the $2^n$ sets;
- **claims problems** (an estate $E$ and claims with $\sum c \ge E$) and their rules:
  - `priority`, the claimants paid in an order, with `priority_range`, the least and the greatest each can receive by order alone, in closed form;
  - `random_arrival`, the mean of priority over all orders: the Shapley value of $w(S) = \max(0, E - \sum_{j \notin S} c_j)$ (O'Neill 1982);
  - `proportional`, `constrained_equal_awards`, `constrained_equal_losses`;
  - `talmud`, Aumann and Maschler's rule, and `contested_garment` for two;
  - `interleaved`, the claimants paid by turns in small steps at given rates.

## 2. Rationale

Dividing what is short among claims is a problem of every institution that settles: bankruptcy, rationing, taxation, the sharing of a cost. The rules and the value name no society, so they are nomos's, beside its games and its social choice.

## 3. Verification

`check_nomos_claims` (28): the glove game's value $(2/3, 1/6, 1/6)$; efficiency, symmetry and additivity of the value on random games; the Talmud's three estates against claims of 100, 200 and 300 (Aumann and Maschler's Table 1) and the contested garment; the Talmud rule's consistency with the contested garment on every pair, on 200 problems; random arrival against the mean over every order by brute force; the three coinciding for two claimants; every rule efficient, bounded by the claims and order-preserving; and payment by turns.

*Found while checking.* Fine interleaving removes the order of the turns but not their rates: at rates in proportion to the claims it tends to the proportional rule, at equal rates to the constrained equal awards. "Paid together" is a family of rules.

nomos's checks rise from 311 to 339.

## 4. Consequences

The value is exact and exponential in the players. Sampling over orders, the nucleolus and the core can follow.
