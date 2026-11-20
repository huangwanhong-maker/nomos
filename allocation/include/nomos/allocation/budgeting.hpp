// What to fund from a budget when voters approve projects: participatory budgeting, and whether the funded projects
// represent groups of voters in proportion to their size.
//
// n voters, projects j with costs c_j, a budget B; voter i approves the set A_i. A rule selects projects W with cost(W) ≤ B.
//
//   utilitarian_greedy   projects by approvals (or approvals per unit of cost), each funded if it still fits: the rule most
//                        cities use, which lets a majority take everything
//   equal_shares         the Method of Equal Shares (Peters, Pierczyński and Skowron 2021, NeurIPS; equalshares.net): each
//                        voter holds B/n; a project is bought by its supporters, each paying min(their money, ρ u_ij), with the
//                        least price per unit of satisfaction ρ among the projects its supporters can afford; under cost
//                        satisfaction u_ij = c_j (each pays at most ρ c_j), under cardinality u_ij = 1. With unit costs
//                        it satisfies extended justified representation (EJR; Peters and Skowron 2020); with unequal costs
//                        EJR can fail to be satisfiable and its guarantee is EJR up to one project (Peters, Pierczyński and
//                        Skowron 2021). It can leave money unspent, so a declared completion
//                        (here, utilitarian greedy with what is left) follows; it is neither exhaustive nor committee-monotone
//   phragmen             sequential Phragmén with costs: voters earn money at one rate, and a project is bought as soon as its
//                        supporters together hold its cost; its load spread so that the supporters' loads are equal (Brill,
//                        Freeman, Janson and Lackner 2017 for committees; Los, Christoff and Grossi 2022 for costs)
//   check_ejr            whether a selection gives every cohesive group its due (approval ballots, cardinality satisfaction):
//                        a group S of voters that could afford a set T of projects all of them approve, |S| ≥ n cost(T)/B,
//                        must contain a voter with at least |T| approved projects funded (Peters et al. 2021's definition for
//                        costs); checked by enumerating groups, for small instances
#pragma once

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace nomos::allocation {

struct Election {
  double budget = 0;
  std::vector<double> costs;                       // per project
  std::vector<std::vector<std::size_t>> approvals; // per voter, the projects approved
};

enum class Satisfaction { cardinality, cost };

struct Selection {
  std::vector<std::size_t> funded;   // in the order they were bought
  double spent = 0;
  std::vector<double> paid;          // equal shares: what each voter paid
  std::vector<std::size_t> completed;  // projects added by the completion step
  std::string said;
};

namespace detail {

inline void check_election(const Election& e) {
  if (!(e.budget >= 0)) throw std::invalid_argument("nomos::allocation: a budget that is not negative");
  for (double c : e.costs)
    if (!(c > 0)) throw std::invalid_argument("nomos::allocation: project costs above zero");
  for (const auto& a : e.approvals)
    for (std::size_t j : a)
      if (j >= e.costs.size()) throw std::invalid_argument("nomos::allocation: a ballot approves a project that is not one");
}

inline std::vector<std::vector<std::size_t>> supporters(const Election& e) {
  std::vector<std::vector<std::size_t>> s(e.costs.size());
  for (std::size_t i = 0; i < e.approvals.size(); i++)
    for (std::size_t j : e.approvals[i]) s[j].push_back(i);
  return s;
}

}  // namespace detail

inline Selection utilitarian_greedy(const Election& e, bool per_cost = false, std::vector<std::size_t> already = {}, double spent = 0) {
  detail::check_election(e);
  const auto S = detail::supporters(e);
  std::vector<std::size_t> order(e.costs.size());
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
    const double ka = static_cast<double>(S[a].size()) / (per_cost ? e.costs[a] : 1.0), kb = static_cast<double>(S[b].size()) / (per_cost ? e.costs[b] : 1.0);
    return ka > kb;
  });
  Selection s;
  s.funded = already, s.spent = spent;
  for (std::size_t j : order) {
    if (std::find(s.funded.begin(), s.funded.end(), j) != s.funded.end()) continue;
    if (S[j].empty()) continue;
    if (s.spent + e.costs[j] <= e.budget + 1e-9) s.funded.push_back(j), s.spent += e.costs[j];
  }
  s.said = per_cost ? "projects by approvals per unit of cost, each funded if it still fits" : "projects by approvals, each funded if it still fits";
  return s;
}

inline Selection equal_shares(const Election& e, Satisfaction sat = Satisfaction::cost, bool complete = true) {
  detail::check_election(e);
  const std::size_t n = e.approvals.size(), m = e.costs.size();
  if (n == 0) throw std::invalid_argument("nomos::allocation::equal_shares: at least one voter");
  const auto S = detail::supporters(e);
  std::vector<double> money(n, e.budget / static_cast<double>(n));
  std::vector<bool> funded(m, false);
  Selection s;
  s.paid.assign(n, 0);
  for (;;) {
    std::size_t best = static_cast<std::size_t>(-1);
    double best_rho = 1e300;
    for (std::size_t j = 0; j < m; j++) {
      if (funded[j] || S[j].empty()) continue;
      double held = 0;
      for (std::size_t i : S[j]) held += money[i];
      if (held < e.costs[j] - 1e-9) continue;   // its supporters cannot afford it
      const double u = sat == Satisfaction::cost ? e.costs[j] : 1.0;
      // The least ρ with Σ min(money_i, ρ u) = c_j: supporters sorted by money, the poorest paying all they have.
      std::vector<double> ms;
      for (std::size_t i : S[j]) ms.push_back(money[i]);
      std::sort(ms.begin(), ms.end());
      double paid = 0, rho = 1e300;
      for (std::size_t k = 0; k < ms.size(); k++) {
        const double rest = static_cast<double>(ms.size() - k);
        const double r = (e.costs[j] - paid) / (rest * u);
        if (r * u <= ms[k] + 1e-12) {
          rho = r;
          break;
        }
        paid += ms[k];
      }
      // Ties broken by the smaller index, then by more approvals (stable over projects in their order).
      if (rho < best_rho - 1e-15) best_rho = rho, best = j;
    }
    if (best == static_cast<std::size_t>(-1)) break;
    funded[best] = true;
    s.funded.push_back(best);
    s.spent += e.costs[best];
    const double u = sat == Satisfaction::cost ? e.costs[best] : 1.0;
    for (std::size_t i : S[best]) {
      const double pay = std::min(money[i], best_rho * u);
      money[i] -= pay, s.paid[i] += pay;
    }
  }
  s.said = std::string("the Method of Equal Shares under ") + (sat == Satisfaction::cost ? "cost" : "cardinality") + " satisfaction";
  if (complete) {
    const auto before = s.funded;
    const auto c = utilitarian_greedy(e, false, s.funded, s.spent);
    for (std::size_t j : c.funded)
      if (std::find(before.begin(), before.end(), j) == before.end()) s.completed.push_back(j);
    s.funded = c.funded, s.spent = c.spent;
    s.said += ", completed by utilitarian greedy with the money left (" + std::to_string(s.completed.size()) + " projects added)";
  }
  return s;
}

inline Selection phragmen(const Election& e) {
  detail::check_election(e);
  const std::size_t n = e.approvals.size(), m = e.costs.size();
  const auto S = detail::supporters(e);
  std::vector<double> load(n, 0);
  std::vector<bool> out(m, false);
  Selection s;
  for (;;) {
    std::size_t best = static_cast<std::size_t>(-1);
    double best_load = 1e300;
    for (std::size_t j = 0; j < m; j++) {
      if (out[j] || S[j].empty() || s.spent + e.costs[j] > e.budget + 1e-9) continue;
      double sum = e.costs[j];
      for (std::size_t i : S[j]) sum += load[i];
      const double t = sum / static_cast<double>(S[j].size());
      if (t < best_load - 1e-15) best_load = t, best = j;
    }
    if (best == static_cast<std::size_t>(-1)) break;
    out[best] = true;
    s.funded.push_back(best);
    s.spent += e.costs[best];
    for (std::size_t i : S[best]) load[i] = best_load;
  }
  s.said = "sequential Phragmén with costs: each project bought when its supporters' equalised load is least";
  return s;
}

struct EjrReport {
  bool satisfied = true;
  std::vector<std::size_t> group;      // the cohesive group left furthest short of its due (the largest such group on ties)
  std::vector<std::size_t> projects;   // the set it could afford
  std::size_t shortfall = 0;           // |T| less the most any of its members has funded
};

// EJR for approval ballots with cardinality satisfaction, by enumerating groups of voters (n ≤ 16) and the projects they all
// approve (at most 12 in common).
inline EjrReport check_ejr(const Election& e, const std::vector<std::size_t>& funded) {
  detail::check_election(e);
  const std::size_t n = e.approvals.size(), m = e.costs.size();
  if (n > 16) throw std::invalid_argument("nomos::allocation::check_ejr: at most 16 voters to enumerate");
  std::vector<bool> in(m, false);
  for (std::size_t j : funded) in[j] = true;
  std::vector<std::size_t> sat(n, 0);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j : e.approvals[i]) sat[i] += in[j];
  EjrReport r;
  for (std::uint32_t mask = 1; mask < (1u << n); mask++) {
    std::vector<std::size_t> group;
    for (std::size_t i = 0; i < n; i++)
      if (mask >> i & 1u) group.push_back(i);
    // Projects every member approves.
    std::vector<std::size_t> common;
    for (std::size_t j = 0; j < m; j++) {
      bool all = true;
      for (std::size_t i : group) all = all && std::find(e.approvals[i].begin(), e.approvals[i].end(), j) != e.approvals[i].end();
      if (all) common.push_back(j);
    }
    if (common.empty() || common.size() > 12) continue;
    std::size_t best_sat = 0;
    for (std::size_t i : group) best_sat = std::max(best_sat, sat[i]);
    const double share = static_cast<double>(group.size()) / static_cast<double>(n) * e.budget;
    // The largest T ⊆ common affordable with the group's share: the cheapest projects first.
    std::vector<double> cs;
    for (std::size_t j : common) cs.push_back(e.costs[j]);
    std::sort(cs.begin(), cs.end());
    double c = 0;
    std::size_t t = 0;
    while (t < cs.size() && c + cs[t] <= share + 1e-9) c += cs[t++];
    if (t > best_sat) {
      r.satisfied = false;
      const std::size_t gap = t - best_sat;
      if (gap > r.shortfall || (gap == r.shortfall && group.size() > r.group.size())) {
        r.shortfall = gap;
        r.group = group;
        r.projects.clear();
        std::vector<std::pair<double, std::size_t>> by;
        for (std::size_t j : common) by.push_back({e.costs[j], j});
        std::sort(by.begin(), by.end());
        for (std::size_t k = 0; k < t; k++) r.projects.push_back(by[k].second);
      }
    }
  }
  return r;
}

}  // namespace nomos::allocation
