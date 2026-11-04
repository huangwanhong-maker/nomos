// Positional and sequential rules: plurality, Borda, approval and instant runoff.
//
// - Plurality counts first preferences.
// - The Borda count (Borda 1781) gives a candidate, on each ballot, one point for every candidate
//   ranked below it and half a point for every candidate tied with it; on a complete ballot of k
//   candidates the one in place r (from 0) receives k − 1 − r. Equivalently, the Borda score of a
//   is Σ_{b≠a} (d(a, b) + ½ t(a, b)), t(a, b) the weight of ballots tying them: a row sum of the
//   pairwise matrix, which with complete ballots is Σ_b d(a, b) (Saari 1995; Young 1974
//   characterises the rule). The scores sum to W k(k − 1)/2.
// - Approval voting (Brams and Fishburn 1978): each voter approves a set of candidates, and each
//   approved candidate receives the voter's weight.
// - Instant runoff, the alternative vote (Hare 1859, for one seat): the candidate with the fewest
//   first preferences among those remaining is eliminated, each ballot counting for its highest
//   remaining candidate, until one candidate holds more than half the weight of the ballots still
//   counting (a ballot whose ranked candidates are all eliminated is exhausted and counts no
//   more), or one remains. Ties for elimination are broken against the candidate of larger index
//   (a stated, arbitrary rule; the result records that a tie was broken).
//
// Plurality, Borda and instant runoff may each elect a different candidate, and none of them is
// Condorcet consistent; the plurality winner may even be the Condorcet loser (condorcet.hpp; the
// checks give such a profile).

#pragma once

#include <vector>

#include "profile.hpp"

namespace nomos::choice {

// First preferences.
inline std::vector<double> plurality(const Profile& p) {
  detail::check(p);
  std::vector<double> s(p.candidates, 0.0);
  for (const auto& b : p.ballots)
    if (!b.ranking.empty()) {
      detail::require(b.ranking[0] < p.candidates, "a ranked candidate within range");
      s[b.ranking[0]] += b.weight;
    }
  return s;
}

// Borda scores, ties counting half.
inline std::vector<double> borda(const Profile& p) {
  detail::check(p);
  const std::size_t k = p.candidates;
  std::vector<double> s(k, 0.0);
  for (const auto& b : p.ballots) {
    const auto pos = detail::positions(b, k);
    const std::size_t ranked = b.ranking.size(), unranked = k - ranked;
    for (std::size_t c = 0; c < k; c++)
      s[c] += b.weight * (pos[c] < k ? static_cast<double>(k - 1 - pos[c]) : 0.5 * static_cast<double>(unranked - 1));
  }
  return s;
}

// Approval scores: approved[v] is the set approved by voter v, of weight weights[v] (1 if none given).
inline std::vector<double> approval(std::size_t candidates, const std::vector<std::vector<std::size_t>>& approved, std::span<const double> weights = {}) {
  detail::require(weights.empty() || weights.size() == approved.size(), "a weight for each voter");
  std::vector<double> s(candidates, 0.0);
  for (std::size_t v = 0; v < approved.size(); v++) {
    std::vector<bool> seen(candidates, false);
    for (std::size_t c : approved[v]) {
      detail::require(c < candidates && !seen[c], "approved candidates within range, each once");
      seen[c] = true;
      s[c] += weights.empty() ? 1.0 : weights[v];
    }
  }
  return s;
}

struct Runoff {
  std::size_t winner = 0;
  std::vector<std::size_t> eliminated;    // in order of elimination
  std::vector<std::vector<double>> rounds;  // the count of each round (eliminated candidates at 0)
  bool tie_broken = false;                // whether a tie for elimination was broken by index
};

inline Runoff instant_runoff(const Profile& p) {
  detail::check(p);
  const std::size_t k = p.candidates;
  std::vector<bool> out(k, false);
  Runoff r;
  for (std::size_t left = k;; left--) {
    std::vector<double> count(k, 0.0);
    double active = 0;
    for (const auto& b : p.ballots)
      for (std::size_t c : b.ranking) {
        detail::require(c < k, "a ranked candidate within range");
        if (!out[c]) {
          count[c] += b.weight;
          active += b.weight;
          break;
        }
      }
    r.rounds.push_back(count);
    std::size_t lead = k, last = k, tied = 0;
    for (std::size_t c = 0; c < k; c++) {
      if (out[c]) continue;
      if (lead == k || count[c] > count[lead]) lead = c;
      if (last == k || count[c] < count[last]) last = c, tied = 1;
      else if (count[c] == count[last]) last = c, tied++;  // the larger index goes first
    }
    if (left == 1 || count[lead] > active / 2) {
      r.winner = lead;
      return r;
    }
    if (tied > 1) r.tie_broken = true;
    out[last] = true;
    r.eliminated.push_back(last);
  }
}

}  // namespace nomos::choice
