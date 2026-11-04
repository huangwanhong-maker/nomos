// Housing markets: Gale's top trading cycles.
//
// In a housing market (Shapley and Scarf 1974) each of n agents owns one indivisible good, its
// house (agent i owns house i), and ranks the houses strictly. An allocation gives each agent one
// house, a permutation. Top trading cycles, which Shapley and Scarf attribute to David Gale: every
// remaining agent points to the owner of the best remaining house on its list; following the
// pointers from any agent leads into a cycle, since each agent points to exactly one; every agent
// on a cycle receives the house it points to, and leaves with it; repeated until none remain.
//
// The result is in the core: no coalition can reallocate its own houses so that each member is at
// least as well off and one better off (Shapley and Scarf 1974); with strict
// preferences it is the unique allocation in the core, indeed in the strict core (Roth and
// Postlewaite 1977), and truthful reporting is a dominant strategy (Roth 1982b). It is therefore
// individually rational and Pareto efficient.
//
// Lists: an agent may list only the houses it prefers to its own; its own house is appended at
// the end when not listed, and houses after it never matter (no agent ends below its own).
// Cost: O(n²) by following pointers in each round, n rounds at most.

#pragma once

#include <vector>

#include "deferred_acceptance.hpp"

namespace nomos::matching {

// The house each agent receives (house i being agent i's to begin with).
inline std::vector<std::size_t> top_trading_cycles(const Preferences& preferences) {
  const std::size_t n = preferences.size();
  std::vector<std::vector<std::size_t>> list(n);
  for (std::size_t i = 0; i < n; i++) {
    for (std::size_t h : preferences[i]) {
      detail::require(h < n, "a house within range");
      list[i].push_back(h);
      if (h == i) break;  // houses below its own never matter
    }
    if (list[i].empty() || list[i].back() != i) list[i].push_back(i);
  }
  std::vector<std::size_t> gets(n, unmatched), next(n, 0), points(n, unmatched);
  std::vector<bool> gone(n, false);
  std::size_t remaining = n;
  while (remaining > 0) {
    // Each remaining agent points to the owner of its best remaining house (owners keep their houses until they leave).
    for (std::size_t i = 0; i < n; i++) {
      if (gone[i]) continue;
      while (gone[list[i][next[i]]]) next[i]++;
      points[i] = list[i][next[i]];
    }
    // Find every cycle: walk from each remaining agent, marking the walk.
    std::vector<std::size_t> mark(n, unmatched);
    std::vector<std::size_t> leaving;
    for (std::size_t s = 0; s < n; s++) {
      if (gone[s] || mark[s] != unmatched) continue;
      std::size_t v = s;
      while (mark[v] == unmatched) {
        mark[v] = s;
        v = points[v];
      }
      if (mark[v] == s) {  // a new cycle through v
        std::size_t u = v;
        do {
          leaving.push_back(u);
          u = points[u];
        } while (u != v);
      }
    }
    for (std::size_t u : leaving) gets[u] = points[u];
    for (std::size_t u : leaving) gone[u] = true;
    remaining -= leaving.size();
  }
  return gets;
}

}  // namespace nomos::matching
