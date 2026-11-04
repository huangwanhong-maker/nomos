// Two-sided matching: Gale and Shapley's deferred acceptance, one-to-one and many-to-one.
//
// A market has applicants a (each to be matched with at most one institution) and institutions c
// (each with a capacity q_c ≥ 1). Each lists the partners it finds acceptable, best first; an
// unlisted partner is unacceptable. An institution's preferences over sets of applicants are
// responsive to its list (it prefers a set with a better applicant in place of a worse one, and
// adding an acceptable applicant to a set below capacity), the model of college admissions (Gale
// and Shapley 1962; Roth and Sotomayor 1990, ch. 5). With every capacity one it is the marriage
// market of the one-to-one case (Gale and Shapley 1962; Roth and Sotomayor 1990, ch. 2).
//
// An assignment is stable when it is individually rational (everyone matched only to acceptable
// partners, within capacity) and no applicant and institution block it: prefer each other to
// what they have, the institution having a free place or an assigned applicant it likes less
// (stability.hpp checks this).
//
// Deferred acceptance, with applicants proposing: every unassigned applicant proposes to the best
// institution on their list that has not yet rejected them; each institution holds the best of
// its proposers up to its capacity, among those it finds acceptable, and rejects the rest;
// repeated until no unassigned applicant has an institution left to propose to (Gale and Shapley
// 1962). With institutions proposing, each institution offers places to the best applicants
// on its list that have not rejected it, keeping as many offers outstanding as it has places, and
// each applicant holds the best offer received, rejecting the rest.
//
// The result is stable whatever the order of proposals; when applicants propose it is the stable
// assignment every applicant likes best among all stable assignments (applicant-optimal), and the
// one every institution likes least; when institutions propose, the reverse (Gale and Shapley
// 1962, Theorems 1 and 2; Roth and Sotomayor 1990, ch. 2 and 5). The same agents are matched,
// and every institution fills the same number of places, in every stable assignment (McVitie and
// Wilson 1970; Roth 1986, the "rural hospitals" theorem). Proposing truthfully is a dominant
// strategy for the proposing side (Dubins and Freedman 1981; Roth 1982a), though not for the other.
//
// Cost: O(L) proposals, L the total length of the proposing side's lists, each handled in
// O(q log q) for a capacity q (the held applicants are kept in order of rank by sorting).

#pragma once

#include <algorithm>
#include <cstddef>
#include <deque>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace nomos::matching {

// For each agent of one side, the acceptable agents of the other side, best first.
using Preferences = std::vector<std::vector<std::size_t>>;

inline constexpr std::size_t unmatched = std::numeric_limits<std::size_t>::max();

enum class Proposing { applicants, institutions };

struct Market {
  Preferences applicants;              // each applicant's acceptable institutions, best first
  Preferences institutions;            // each institution's acceptable applicants, best first
  std::vector<std::size_t> capacity;   // q_c; empty means one place each (the one-to-one case)
};

struct Assignment {
  std::vector<std::size_t> institution_of;          // each applicant's institution, or unmatched
  std::vector<std::vector<std::size_t>> applicants_of;  // each institution's applicants, in its order of preference
};

namespace detail {

inline void require(bool ok, const char* what) {
  if (!ok) throw std::invalid_argument(std::string("nomos::matching: ") + what);
}

// rank[x][y]: the position of y on x's list, or `unmatched` when y is unacceptable to x.
inline std::vector<std::vector<std::size_t>> ranks(const Preferences& p, std::size_t others) {
  std::vector<std::vector<std::size_t>> r(p.size(), std::vector<std::size_t>(others, unmatched));
  for (std::size_t x = 0; x < p.size(); x++)
    for (std::size_t k = 0; k < p[x].size(); k++) {
      require(p[x][k] < others, "a listed partner within range");
      require(r[x][p[x][k]] == unmatched, "no partner listed twice");
      r[x][p[x][k]] = k;
    }
  return r;
}

inline std::vector<std::size_t> capacities(const Market& m) {
  if (m.capacity.empty()) return std::vector<std::size_t>(m.institutions.size(), 1);
  require(m.capacity.size() == m.institutions.size(), "a capacity for each institution");
  return m.capacity;
}

}  // namespace detail

// The stable assignment that the proposing side likes best.
inline Assignment deferred_acceptance(const Market& m, Proposing side = Proposing::applicants) {
  const std::size_t na = m.applicants.size(), nc = m.institutions.size();
  const auto ra = detail::ranks(m.applicants, nc), rc = detail::ranks(m.institutions, na);
  const auto cap = detail::capacities(m);
  Assignment out;
  out.institution_of.assign(na, unmatched);
  out.applicants_of.assign(nc, {});
  auto by_rank = [&](std::size_t c) {
    std::sort(out.applicants_of[c].begin(), out.applicants_of[c].end(), [&](std::size_t x, std::size_t y) { return rc[c][x] < rc[c][y]; });
  };
  if (side == Proposing::applicants) {
    std::vector<std::size_t> next(na, 0);
    std::deque<std::size_t> free;
    for (std::size_t a = 0; a < na; a++) free.push_back(a);
    while (!free.empty()) {
      const std::size_t a = free.front();
      free.pop_front();
      while (next[a] < m.applicants[a].size()) {
        const std::size_t c = m.applicants[a][next[a]++];
        if (rc[c][a] == unmatched || cap[c] == 0) continue;  // unacceptable to c: rejected at once
        auto& held = out.applicants_of[c];
        if (held.size() < cap[c]) {
          held.push_back(a);
          out.institution_of[a] = c;
          by_rank(c);
          break;
        }
        const std::size_t worst = held.back();
        if (rc[c][a] < rc[c][worst]) {
          held.back() = a;
          out.institution_of[a] = c;
          out.institution_of[worst] = unmatched;
          by_rank(c);
          free.push_back(worst);
          break;
        }
      }
    }
  } else {
    // Institutions propose; each applicant holds her best offer.
    std::vector<std::size_t> next(nc, 0), held(nc, 0);
    for (bool progress = true; progress;) {
      progress = false;
      for (std::size_t c = 0; c < nc; c++)
        while (held[c] < cap[c] && next[c] < m.institutions[c].size()) {
          const std::size_t a = m.institutions[c][next[c]++];
          progress = true;
          if (ra[a][c] == unmatched) continue;  // c unacceptable to a
          const std::size_t d = out.institution_of[a];
          if (d == unmatched) {
            out.institution_of[a] = c;
            held[c]++;
          } else if (ra[a][c] < ra[a][d]) {
            out.institution_of[a] = c;
            held[c]++;
            held[d]--;
          }
        }
    }
    for (std::size_t a = 0; a < na; a++)
      if (out.institution_of[a] != unmatched) out.applicants_of[out.institution_of[a]].push_back(a);
    for (std::size_t c = 0; c < nc; c++) by_rank(c);
  }
  return out;
}

// The one-to-one case: proposers and receivers, each matched to at most one; the stable matching
// the proposers like best. partner_of_proposer[i] is i's receiver (or unmatched), and conversely.
struct Matching {
  std::vector<std::size_t> partner_of_proposer, partner_of_receiver;
};

inline Matching stable_matching(const Preferences& proposers, const Preferences& receivers) {
  const Market m{proposers, receivers, {}};
  const auto a = deferred_acceptance(m, Proposing::applicants);
  Matching out;
  out.partner_of_proposer = a.institution_of;
  out.partner_of_receiver.assign(receivers.size(), unmatched);
  for (std::size_t r = 0; r < receivers.size(); r++)
    if (!a.applicants_of[r].empty()) out.partner_of_receiver[r] = a.applicants_of[r][0];
  return out;
}

}  // namespace nomos::matching
