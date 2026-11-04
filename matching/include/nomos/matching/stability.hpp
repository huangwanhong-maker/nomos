// Stability of an assignment: individual rationality and the absence of blocking pairs.
//
// An assignment of applicants to institutions (deferred_acceptance.hpp) is individually rational
// when every applicant is assigned only to an institution she lists, every institution only to
// applicants it lists, and no institution beyond its capacity. An applicant a and an institution
// c block it when they find each other acceptable, a is not assigned to c, a prefers c to her
// assignment (or is unassigned), and c either has a free place or prefers a to one of its assigned
// applicants, equivalently to the worst of them, its preferences being responsive (Gale and
// Shapley 1962; Roth and Sotomayor 1990, ch. 5, where with responsive preferences this pairwise
// stability is stability). It is stable when individually rational and unblocked.
//
// Cost: O(Σ list lengths) after O(n_a n_c) tables of ranks.

#pragma once

#include <vector>

#include "deferred_acceptance.hpp"

namespace nomos::matching {

struct BlockingPair {
  std::size_t applicant = 0, institution = 0;
  bool operator==(const BlockingPair&) const = default;
};

inline bool individually_rational(const Market& m, const Assignment& x) {
  const std::size_t na = m.applicants.size(), nc = m.institutions.size();
  if (x.institution_of.size() != na || x.applicants_of.size() != nc) return false;
  const auto ra = detail::ranks(m.applicants, nc), rc = detail::ranks(m.institutions, na);
  const auto cap = detail::capacities(m);
  for (std::size_t c = 0; c < nc; c++) {
    if (x.applicants_of[c].size() > cap[c]) return false;
    for (std::size_t a : x.applicants_of[c])
      if (a >= na || x.institution_of[a] != c || rc[c][a] == unmatched) return false;
  }
  for (std::size_t a = 0; a < na; a++) {
    const std::size_t c = x.institution_of[a];
    if (c == unmatched) continue;
    if (c >= nc || ra[a][c] == unmatched) return false;
    if (std::find(x.applicants_of[c].begin(), x.applicants_of[c].end(), a) == x.applicants_of[c].end()) return false;
  }
  return true;
}

// Every blocking pair, in order of applicant and then of institution.
inline std::vector<BlockingPair> blocking_pairs(const Market& m, const Assignment& x) {
  const std::size_t na = m.applicants.size(), nc = m.institutions.size();
  const auto ra = detail::ranks(m.applicants, nc), rc = detail::ranks(m.institutions, na);
  const auto cap = detail::capacities(m);
  std::vector<BlockingPair> out;
  for (std::size_t a = 0; a < na; a++) {
    const std::size_t mine = x.institution_of[a];
    const std::size_t my_rank = mine == unmatched ? unmatched : ra[a][mine];
    for (std::size_t c : m.applicants[a]) {
      if (c == mine || ra[a][c] >= my_rank || rc[c][a] == unmatched) continue;
      bool wants = x.applicants_of[c].size() < cap[c];
      for (std::size_t b : x.applicants_of[c]) wants = wants || rc[c][a] < rc[c][b];
      if (wants) out.push_back({a, c});
    }
  }
  return out;
}

inline bool is_stable(const Market& m, const Assignment& x) { return individually_rational(m, x) && blocking_pairs(m, x).empty(); }

}  // namespace nomos::matching
