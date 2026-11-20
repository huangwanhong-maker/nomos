// What each receives of indivisible goods: the fairness of a division by additive valuations, and the procedures that give
// divisions with known guarantees.
//
// n agents, m goods, agent i valuing good g at v_ig ≥ 0, a bundle at the sum of its goods (additive valuations). A division
// gives each good to one agent. Its tests (the survey, Part II §4):
//
//   envy-free (EF)      every agent values its own bundle at least as much as anyone else's
//   EF1                 any envy disappears on removing some good from the envied bundle (Budish 2011; Lipton et al. 2004)
//   EFX                 any envy disappears on removing any good from the envied bundle (Caragiannis et al. 2019)
//   proportional (PROP) every agent receives at least 1/n of its value of all the goods
//   maximin share (MMS) every agent receives at least the most it could guarantee itself by dividing the goods into n bundles
//                       and taking the least (Budish 2011); it can fail to exist (Kurokawa, Procaccia and Wang 2018)
//
// And the procedures: round-robin picking (EF1); envy-cycle elimination (Lipton, Markakis, Mossel and Saberi 2004: EF1 for
// any monotone valuation); the division that maximises the product of values, the maximum Nash welfare, which is EF1 and
// Pareto optimal (Caragiannis, Kurokawa, Moulin, Procaccia, Shah and Wang 2019, ACM TEAC 7(3)) — an optimisation, and said
// to be one, here by enumeration for small instances.
#pragma once

#include <nomos/allocation/assignment.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace nomos::allocation {

using Values = std::vector<std::vector<double>>;   // v[agent][good]
using Division = std::vector<std::size_t>;         // for each good, its agent

inline std::vector<std::vector<std::size_t>> bundles_of(const Division& d, std::size_t agents) {
  std::vector<std::vector<std::size_t>> b(agents);
  for (std::size_t g = 0; g < d.size(); g++) b.at(d[g]).push_back(g);
  return b;
}

inline double value_of(const Values& v, std::size_t i, const std::vector<std::size_t>& bundle) {
  double s = 0;
  for (std::size_t g : bundle) s += v[i][g];
  return s;
}

struct FairnessReport {
  bool envy_free = true, ef1 = true, efx = true, proportional = true, mms = true;
  std::vector<double> mms_share;      // each agent's maximin share
  std::vector<std::string> witnesses; // the first pair that breaks each property
};

// Each agent's maximin share by enumerating the n^m divisions (m ≤ 12, n^m ≤ 2·10^7).
inline std::vector<double> maximin_shares(const Values& v) {
  const std::size_t n = v.size(), m = n ? v[0].size() : 0;
  if (std::pow(static_cast<double>(n), static_cast<double>(m)) > 2e7) throw std::invalid_argument("nomos::allocation::maximin_shares: too many divisions to enumerate");
  std::vector<double> out(n, 0);
  std::vector<std::size_t> d(m, 0);
  std::vector<double> sums(n);
  for (;;) {
    for (std::size_t i = 0; i < n; i++) {
      std::fill(sums.begin(), sums.end(), 0.0);
      for (std::size_t g = 0; g < m; g++) sums[d[g]] += v[i][g];
      out[i] = std::max(out[i], *std::min_element(sums.begin(), sums.end()));
    }
    std::size_t k = 0;
    while (k < m && ++d[k] == n) d[k++] = 0;
    if (k == m) break;
  }
  return out;
}

inline FairnessReport check_fairness(const Values& v, const Division& d, bool with_mms = true) {
  const std::size_t n = v.size();
  const auto B = bundles_of(d, n);
  FairnessReport r;
  for (std::size_t i = 0; i < n; i++) {
    const double own = value_of(v, i, B[i]);
    double all = 0;
    for (double x : v[i]) all += x;
    if (own < all / static_cast<double>(n) - 1e-12) {
      if (r.proportional) r.witnesses.push_back("agent " + std::to_string(i) + " receives less than 1/n of its value of all the goods");
      r.proportional = false;
    }
    for (std::size_t j = 0; j < n; j++) {
      if (j == i) continue;
      const double other = value_of(v, i, B[j]);
      if (own >= other - 1e-12) continue;
      if (r.envy_free) r.witnesses.push_back("agent " + std::to_string(i) + " envies agent " + std::to_string(j));
      r.envy_free = false;
      double most = 0, least = 1e300;
      for (std::size_t g : B[j]) most = std::max(most, v[i][g]), least = std::min(least, v[i][g]);
      if (own < other - most - 1e-12) {
        if (r.ef1) r.witnesses.push_back("agent " + std::to_string(i) + " envies agent " + std::to_string(j) + " even without its most valued good");
        r.ef1 = false;
      }
      if (own < other - least - 1e-12) {
        if (r.efx) r.witnesses.push_back("agent " + std::to_string(i) + " envies agent " + std::to_string(j) + " even without its least valued good");
        r.efx = false;
      }
    }
  }
  if (with_mms) {
    r.mms_share = maximin_shares(v);
    for (std::size_t i = 0; i < n; i++)
      if (value_of(v, i, B[i]) < r.mms_share[i] - 1e-12) {
        if (r.mms) r.witnesses.push_back("agent " + std::to_string(i) + " receives less than its maximin share");
        r.mms = false;
      }
  }
  return r;
}

inline Division round_robin(const Values& v) {
  const std::size_t n = v.size(), m = n ? v[0].size() : 0;
  const auto b = picking_sequence(v, round_robin_sequence(n, m));
  Division d(m, 0);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t g : b[i]) d[g] = i;
  return d;
}

// Lipton et al.'s envy-cycle elimination: goods in a declared order, each to an agent no one envies; when everyone is envied,
// a cycle of envy is found and its bundles rotated, which removes it.
inline Division envy_cycle(const Values& v, std::vector<std::size_t> goods_order = {}) {
  const std::size_t n = v.size(), m = n ? v[0].size() : 0;
  if (goods_order.empty())
    for (std::size_t g = 0; g < m; g++) goods_order.push_back(g);
  std::vector<std::vector<std::size_t>> B(n);
  auto envies = [&](std::size_t i, std::size_t j) { return value_of(v, i, B[i]) < value_of(v, i, B[j]) - 1e-12; };
  for (std::size_t g : goods_order) {
    for (int guard = 0; guard < static_cast<int>(n * n + 1); guard++) {
      std::size_t free = none;
      for (std::size_t j = 0; j < n && free == none; j++) {
        bool envied = false;
        for (std::size_t i = 0; i < n && !envied; i++) envied = i != j && envies(i, j);
        if (!envied) free = j;
      }
      if (free != none) {
        B[free].push_back(g);
        break;
      }
      // Everyone is envied, so every agent has someone who envies it: walk back from agent 0 to an agent that envies it, and
      // on, until an agent repeats. On the cycle so found each agent envies the one before it; each takes that one's bundle.
      std::vector<std::size_t> path{0}, at(n, none);
      at[0] = 0;
      for (;;) {
        const std::size_t u = path.back();
        std::size_t envier = none;
        for (std::size_t w = 0; w < n && envier == none; w++)
          if (w != u && envies(w, u)) envier = w;
        if (at[envier] != none) {
          // path[at[envier]] … path.back(): path[k + 1] envies path[k], and path[at[envier]] (= envier) envies path.back().
          const std::size_t from = at[envier];
          std::vector<std::size_t> cycle(path.begin() + static_cast<std::ptrdiff_t>(from), path.end());
          const std::size_t L = cycle.size();
          std::vector<std::vector<std::size_t>> taken(L);
          for (std::size_t k = 0; k < L; k++) taken[k] = B[k == 0 ? cycle[L - 1] : cycle[k - 1]];
          for (std::size_t k = 0; k < L; k++) B[cycle[k]] = taken[k];
          break;
        }
        at[envier] = path.size();
        path.push_back(envier);
      }
    }
  }
  Division d(m, 0);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t g : B[i]) d[g] = i;
  return d;
}

// The maximum Nash welfare by enumeration (n^m ≤ 2·10^7): first the most agents with a positive value, then the largest product
// among them, as Caragiannis et al. define it when some agent must receive nothing of value.
inline Division max_nash_welfare(const Values& v) {
  const std::size_t n = v.size(), m = n ? v[0].size() : 0;
  if (std::pow(static_cast<double>(n), static_cast<double>(m)) > 2e7) throw std::invalid_argument("nomos::allocation::max_nash_welfare: too many divisions to enumerate");
  Division d(m, 0), best = d;
  std::size_t best_positive = 0;
  double best_log = -1e300;
  std::vector<double> u(n);
  for (;;) {
    std::fill(u.begin(), u.end(), 0.0);
    for (std::size_t g = 0; g < m; g++) u[d[g]] += v[d[g]][g];
    std::size_t positive = 0;
    double lg = 0;
    for (double x : u)
      if (x > 0) positive++, lg += std::log(x);
    if (positive > best_positive || (positive == best_positive && lg > best_log + 1e-12)) best_positive = positive, best_log = lg, best = d;
    std::size_t k = 0;
    while (k < m && ++d[k] == n) d[k++] = 0;
    if (k == m) break;
  }
  return best;
}

}  // namespace nomos::allocation
