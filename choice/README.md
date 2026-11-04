# choice: Social Choice

*Status: draft (2026-10-25), [decision 0003](../development/decisions/0003-the-first-modules.md).* Proportional representation (single transferable vote, divisor methods of apportionment) and the axioms of impossibility as executable checks are planned.

## 1. Concepts

- **Profiles.** Each ballot ranks candidates, best first, and carries a weight (its number of voters); unranked candidates are tied below the ranked ones. The pairwise majority matrix d(a, b), the weight of ballots ranking a above b, is all a Condorcet method reads (Condorcet 1785).
- **Positional and sequential rules.** Plurality counts first preferences; Borda gives each candidate a point for every candidate ranked below it and half a point for every tie (Borda 1781; Young 1974), its scores the row sums of d; approval counts approvals (Brams and Fishburn 1978); instant runoff eliminates the candidate with fewest first preferences until one holds a majority of the ballots still counting (Hare 1859).
- **Condorcet methods.** A Condorcet winner beats every other candidate head to head; the majority relation may cycle (the paradox of Condorcet 1785). Copeland counts pairwise victories (Copeland 1951); Schulze compares the strengths of the strongest paths of pairwise victories, a transitive relation that elects the Condorcet winner when there is one and is monotonic and independent of clones (Schulze 2011); Kemeny–Young chooses the rankings that agree most with the ballots, the maximum likelihood rankings of Condorcet's model (Kemeny 1959; Young and Levenglick 1978; Young 1988), NP-hard in general (Bartholdi, Tovey and Trick 1989) and found here by trying every ranking.
- **Rules disagree.** On one profile plurality, Borda, instant runoff and the Condorcet methods may elect four different candidates, and the plurality winner may be the Condorcet loser.

## 2. Content

| Header | Content |
|---|---|
| `profile.hpp` | `Ballot{ranking, weight}`, `Profile{candidates, ballots}`; `total_weight(p)`, `pairwise(p)` → `calx::linalg::Matrix`; `winners(scores)` |
| `scoring.hpp` | `plurality(p)`, `borda(p)`, `approval(candidates, approved, weights)`, `instant_runoff(p)` → `Runoff{winner, eliminated, rounds, tie_broken}` |
| `condorcet.hpp` | `condorcet_winner(d)`, `condorcet_loser(d)` → `std::optional`; `copeland(d)`; `schulze(d)` → `Schulze{strength, winners, ranking}`; `kemeny_support(d, ranking)`, `kemeny_young(d)` → `Kemeny{rankings, support}` |

Ties for elimination in instant runoff are broken against the candidate of larger index, and the result says when that happened.

## 3. Checks

`tests/check_choice.cpp` (20 checks): the Condorcet paradox (the cycle, no winner or loser, all three Schulze winners, three Kemeny rankings of support 5); the election of 37 voters among four candidates from COMAP's *For All Practical Purposes* (plurality A, Borda B with scores 42, 69, 67, 44, instant runoff D after B and C, the Condorcet winner C elected by Copeland, Schulze and Kemeny–Young, and A the Condorcet loser); Schulze's (2011) example of 45 voters (its pairwise matrix, strongest paths and ranking E > A > C > B > D); two hundred random profiles (Condorcet consistency of Copeland, Schulze and Kemeny–Young, transitivity of Schulze's relation, Kemeny's rankings minimising the total Kendall distance, Borda as row sums adding to W k(k − 1)/2 with truncated ballots, weights as repeated ballots); approval with weights.

## 4. Dependencies

calx `linalg` (`Matrix`). In the checks, calx `random` and `verification`.

## 5. References

- Bartholdi, J., Tovey, C. A. and Trick, M. A. (1989), Voting schemes for which it can be difficult to tell who won the election, *Soc. Choice Welf.* 6, 157–165.
- Borda, J.-C. de (1781), Mémoire sur les élections au scrutin, *Histoire de l'Académie royale des sciences*, Paris.
- Brams, S. J. and Fishburn, P. C. (1978), Approval voting, *Am. Polit. Sci. Rev.* 72, 831–847.
- COMAP (Consortium for Mathematics and Its Applications), *For All Practical Purposes: Mathematical Literacy in Today's World*, W. H. Freeman, the chapter on social choice.
- Condorcet, Marquis de (1785), *Essai sur l'application de l'analyse à la probabilité des décisions rendues à la pluralité des voix*, Paris.
- Copeland, A. H. (1951), A "reasonable" social welfare function, Seminar on Applications of Mathematics to the Social Sciences, University of Michigan (mimeo).
- Hare, T. (1859), *A Treatise on the Election of Representatives*, London.
- Kemeny, J. G. (1959), Mathematics without numbers, *Daedalus* 88, 577–591.
- Saari, D. G. (1995), *Basic Geometry of Voting*, Springer; Saari, D. G. and Merlin, V. R. (1996), The Copeland method I, *Econ. Theory* 8, 51–76.
- Schulze, M. (2011), A new monotonic, clone-independent, reversal symmetric, and Condorcet-consistent single-winner election method, *Soc. Choice Welf.* 36, 267–303.
- Young, H. P. (1974), An axiomatization of Borda's rule, *J. Econ. Theory* 9, 43–52; (1988), Condorcet's theory of voting, *Am. Polit. Sci. Rev.* 82, 1231–1244.
- Young, H. P. and Levenglick, A. (1978), A consistent extension of Condorcet's election principle, *SIAM J. Appl. Math.* 35, 285–300.

## 6. Open Questions

1. **Kemeny rankings of many candidates** by integer programming on calx's simplex method or by branch and bound.
2. **Ranked pairs, minimax and the Smith set**, and proportional rules, when a domain counts real elections.
