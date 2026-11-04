// Condorcet methods: rules read from the pairwise majority matrix d (profile.hpp).
//
// - A Condorcet winner beats every other candidate, d(a, b) > d(b, a) for all b ≠ a; a Condorcet
//   loser is beaten by every other (Condorcet 1785). Neither need exist: the majority relation
//   may cycle (the Condorcet paradox). A rule is Condorcet consistent when it elects the
//   Condorcet winner whenever there is one.
// - Copeland's rule (Copeland 1951; Saari and Merlin 1996) scores each candidate by its pairwise
//   victories, a pairwise tie counting half.
// - Schulze's rule (Schulze 2011): the strength of the direct link from a to b is d(a, b) when a
//   beats b and 0 otherwise (winning votes); the strength of a path is that of its weakest link,
//   and p(a, b) is the strength of the strongest path from a to b, computed by the widest-path
//   variant of the Floyd–Warshall algorithm in O(k³) (Schulze 2011). a is preferred to b
//   when p(a, b) > p(b, a). This relation is transitive (Schulze 2011), so it has maximal
//   elements, the winners (p(a, b) ≥ p(b, a) for every b), and the candidates ordered by how many
//   others they are preferred to give a ranking that extends it. The rule is Condorcet
//   consistent, monotonic and independent of clones (Schulze 2011).
// - The Kemeny–Young rule (Kemeny 1959; Young and Levenglick 1978) chooses the rankings that
//   maximise the support Σ_{a above b} d(a, b), equivalently minimise the total weighted Kendall
//   distance to the ballots. It is the maximum likelihood ranking of Condorcet's model (Young
//   1988) and is Condorcet consistent, but finding it is NP-hard (Bartholdi, Tovey and Trick
//   1989); here every one of the k! rankings is tried, for k ≤ 9 by default.

#pragma once

#include <algorithm>
#include <numeric>
#include <optional>
#include <vector>

#include "profile.hpp"

namespace nomos::choice {

// The candidate that beats every other, if any.
inline std::optional<std::size_t> condorcet_winner(const Matrix& d) {
  for (std::size_t a = 0; a < d.n; a++) {
    bool all = true;
    for (std::size_t b = 0; b < d.n && all; b++) all = a == b || d(a, b) > d(b, a);
    if (all) return a;
  }
  return std::nullopt;
}

// The candidate beaten by every other, if any.
inline std::optional<std::size_t> condorcet_loser(const Matrix& d) {
  for (std::size_t a = 0; a < d.n; a++) {
    bool all = true;
    for (std::size_t b = 0; b < d.n && all; b++) all = a == b || d(a, b) < d(b, a);
    if (all) return a;
  }
  return std::nullopt;
}

// Copeland scores: victories plus half the ties.
inline std::vector<double> copeland(const Matrix& d) {
  std::vector<double> s(d.n, 0.0);
  for (std::size_t a = 0; a < d.n; a++)
    for (std::size_t b = 0; b < d.n; b++)
      if (a != b) s[a] += d(a, b) > d(b, a) ? 1.0 : d(a, b) == d(b, a) ? 0.5 : 0.0;
  return s;
}

struct Schulze {
  Matrix strength;                   // p(a, b), the strongest path from a to b (0 on the diagonal)
  std::vector<std::size_t> winners;  // p(a, b) ≥ p(b, a) for every b
  std::vector<std::size_t> ranking;  // all candidates, a linear extension of the Schulze relation
};

inline Schulze schulze(const Matrix& d) {
  const std::size_t k = d.n;
  Schulze r;
  r.strength = Matrix(k);
  Matrix& p = r.strength;
  for (std::size_t a = 0; a < k; a++)
    for (std::size_t b = 0; b < k; b++)
      if (a != b && d(a, b) > d(b, a)) p(a, b) = d(a, b);
  for (std::size_t i = 0; i < k; i++)
    for (std::size_t a = 0; a < k; a++) {
      if (a == i) continue;
      for (std::size_t b = 0; b < k; b++)
        if (b != a && b != i) p(a, b) = std::max(p(a, b), std::min(p(a, i), p(i, b)));
    }
  std::vector<std::size_t> beaten(k, 0);
  for (std::size_t a = 0; a < k; a++) {
    bool best = true;
    for (std::size_t b = 0; b < k; b++) {
      if (a == b) continue;
      if (p(a, b) > p(b, a)) beaten[a]++;
      if (p(b, a) > p(a, b)) best = false;
    }
    if (best) r.winners.push_back(a);
  }
  r.ranking.resize(k);
  std::iota(r.ranking.begin(), r.ranking.end(), std::size_t{0});
  std::stable_sort(r.ranking.begin(), r.ranking.end(), [&](std::size_t a, std::size_t b) { return beaten[a] > beaten[b]; });
  return r;
}

struct Kemeny {
  std::vector<std::vector<std::size_t>> rankings;  // every ranking of greatest support, in lexicographic order
  double support = 0;                              // Σ_{a above b} d(a, b)
};

// The support Σ_{a above b} d(a, b) of a ranking of all candidates.
inline double kemeny_support(const Matrix& d, std::span<const std::size_t> ranking) {
  double s = 0;
  for (std::size_t i = 0; i < ranking.size(); i++)
    for (std::size_t j = i + 1; j < ranking.size(); j++) s += d(ranking[i], ranking[j]);
  return s;
}

inline Kemeny kemeny_young(const Matrix& d, std::size_t max_candidates = 9, double tolerance = 1e-9) {
  detail::require(d.n > 0 && d.n <= max_candidates, "few enough candidates to try every ranking");
  std::vector<std::size_t> r(d.n);
  std::iota(r.begin(), r.end(), std::size_t{0});
  Kemeny best;
  best.support = -1;
  do {
    const double s = kemeny_support(d, r);
    if (s > best.support + tolerance) {
      best.support = s;
      best.rankings = {r};
    } else if (s >= best.support - tolerance) {
      best.rankings.push_back(r);
    }
  } while (std::next_permutation(r.begin(), r.end()));
  return best;
}

}  // namespace nomos::choice
