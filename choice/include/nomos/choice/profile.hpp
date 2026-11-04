// Profiles of ballots and their pairwise majority matrix.
//
// A ballot ranks candidates 0, …, k − 1, best first, and carries a weight: the number of voters
// who cast it (or a survey weight). It may rank only some of the candidates; those left unranked
// are tied with one another below every ranked candidate (the common reading of truncated
// ballots, as in Schulze 2011). A profile is the list of ballots.
//
// The pairwise majority matrix d(a, b) is the total weight of ballots ranking a above b (Condorcet
// 1785; Schulze 2011). With complete ballots d(a, b) + d(b, a) = W, the total weight; with
// truncated ones the remainder is the weight of ballots tying a and b. a beats b when
// d(a, b) > d(b, a). Every Condorcet method reads the profile through d alone.
//
// Cost: O(Σ_ballots k²).

#pragma once

#include <calx/linalg/dense.hpp>

#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace nomos::choice {

using calx::linalg::Matrix;

struct Ballot {
  std::vector<std::size_t> ranking;  // candidates, best first; the unranked are tied below
  double weight = 1;
};

struct Profile {
  std::size_t candidates = 0;
  std::vector<Ballot> ballots;
};

namespace detail {

inline void require(bool ok, const char* what) {
  if (!ok) throw std::invalid_argument(std::string("nomos::choice: ") + what);
}

// position[c]: c's place on the ballot, or k (below every ranked candidate) when unranked.
inline std::vector<std::size_t> positions(const Ballot& b, std::size_t k) {
  std::vector<std::size_t> pos(k, k);
  for (std::size_t r = 0; r < b.ranking.size(); r++) {
    require(b.ranking[r] < k, "a ranked candidate within range");
    require(pos[b.ranking[r]] == k, "no candidate ranked twice");
    pos[b.ranking[r]] = r;
  }
  return pos;
}

inline void check(const Profile& p) {
  require(p.candidates > 0, "at least one candidate");
  for (const auto& b : p.ballots) require(b.weight >= 0, "nonnegative weights");
}

}  // namespace detail

// The total weight W of the ballots.
inline double total_weight(const Profile& p) {
  double w = 0;
  for (const auto& b : p.ballots) w += b.weight;
  return w;
}

// d(a, b): the weight of ballots ranking a above b.
inline Matrix pairwise(const Profile& p) {
  detail::check(p);
  const std::size_t k = p.candidates;
  Matrix d(k);
  for (const auto& b : p.ballots) {
    const auto pos = detail::positions(b, k);
    for (std::size_t x = 0; x < k; x++)
      for (std::size_t y = 0; y < k; y++)
        if (pos[x] < pos[y]) d(x, y) += b.weight;
  }
  return d;
}

// The candidates with the largest score (within `tolerance`), in increasing order.
inline std::vector<std::size_t> winners(std::span<const double> score, double tolerance = 1e-9) {
  std::vector<std::size_t> w;
  if (score.empty()) return w;
  double best = score[0];
  for (double s : score) best = std::max(best, s);
  for (std::size_t c = 0; c < score.size(); c++)
    if (score[c] >= best - tolerance) w.push_back(c);
  return w;
}

}  // namespace nomos::choice
