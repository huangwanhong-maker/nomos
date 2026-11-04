# matching: Two-Sided Matching and Housing Markets

*Status: draft (2026-10-25), [decision 0003](../development/decisions/0003-the-first-modules.md).*

## 1. Concepts

- **Two-sided markets.** Applicants are to be matched with institutions, each applicant to at most one, each institution to at most its capacity; every agent lists the partners it finds acceptable, best first, and an institution's preferences over sets are responsive to its list. With every capacity one this is the marriage market (Gale and Shapley 1962; Roth and Sotomayor 1990, chs 2 and 5).
- **Stability.** An assignment is stable when it is individually rational and no applicant and institution would both rather be together: she prefers it to her assignment, and it has a free place or an applicant it likes less. With responsive preferences, pairwise stability is stability.
- **Deferred acceptance.** One side proposes in order of preference, the other holds the best offers it has received and rejects the rest, until no proposal is left to make (Gale and Shapley 1962). The result is stable; it is the best stable assignment for every member of the proposing side and the worst for every member of the other (Gale and Shapley 1962, Theorems 1 and 2). The same agents are matched, and every institution fills as many places, in every stable assignment (McVitie and Wilson 1970; Roth 1986). Truth is a dominant strategy for the proposing side (Dubins and Freedman 1981; Roth 1982a), but not for the other.
- **Housing markets.** Each agent owns one house and ranks them all. Gale's top trading cycles, in which each agent points to the owner of its favourite remaining house and every cycle trades, gives the unique allocation in the core (Shapley and Scarf 1974; Roth and Postlewaite 1977), and truth is a dominant strategy (Roth 1982b).

## 2. Content

| Header | Content |
|---|---|
| `deferred_acceptance.hpp` | `Preferences`; `unmatched`; `Market{applicants, institutions, capacity}`; `Assignment{institution_of, applicants_of}`; `enum class Proposing {applicants, institutions}`; `deferred_acceptance(market, side)`; the one-to-one `stable_matching(proposers, receivers)` → `Matching{partner_of_proposer, partner_of_receiver}` |
| `stability.hpp` | `individually_rational(market, x)`, `blocking_pairs(market, x)` → `BlockingPair{applicant, institution}`, `is_stable(market, x)` |
| `top_trading_cycles.hpp` | `top_trading_cycles(preferences)` → the house each agent receives |

## 3. Checks

`tests/check_matching.cpp` (19 checks): Gale and Shapley's three men and three women (the men's and the women's first choices from the two algorithms, exactly three stable matchings, a woman gaining by truncating her list); sixty random one-to-one markets with incomplete lists against every stable matching found by enumeration (stability, the proposers' optimum and the receivers' pessimum and conversely, the same agents matched); every false list of every proposer in thirty markets (no gain); forty random many-to-one markets against every stable assignment (stability, optimality and pessimality, the same number of places filled, and one copy of each institution per place giving the same assignment); top trading cycles by hand, in the core and alone in it on twenty-five random markets, and strategy-proof against all twenty-four lists of each agent.

## 4. Dependencies

None. In the checks, calx `random` and `verification`.

## 5. References

- Dubins, L. E. and Freedman, D. A. (1981), Machiavelli and the Gale–Shapley algorithm, *Am. Math. Mon.* 88, 485–494.
- Gale, D. and Shapley, L. S. (1962), College admissions and the stability of marriage, *Am. Math. Mon.* 69, 9–15.
- McVitie, D. G. and Wilson, L. B. (1970), Stable marriage assignment for unequal sets, *BIT* 10, 295–309.
- Roth, A. E. (1982a), The economics of matching: stability and incentives, *Math. Oper. Res.* 7, 617–628.
- Roth, A. E. (1982b), Incentive compatibility in a market with indivisible goods, *Econ. Lett.* 9, 127–132.
- Roth, A. E. (1986), On the allocation of residents to rural hospitals, *Econometrica* 54, 425–427.
- Roth, A. E. and Postlewaite, A. (1977), Weak versus strong domination in a market with indivisible goods, *J. Math. Econ.* 4, 131–137.
- Roth, A. E. and Sotomayor, M. A. O. (1990), *Two-Sided Matching*, Cambridge University Press.
- Shapley, L. and Scarf, H. (1974), On cores and indivisibility, *J. Math. Econ.* 1, 23–37.

## 6. Open Questions

1. **The lattice of stable matchings** (its enumeration by rotations; Gusfield and Irving 1989) and couples, when a domain needs them.
2. **School choice** with priorities and ties (Abdulkadiroğlu and Sönmez 2003), and the roommates problem (Irving 1985).
