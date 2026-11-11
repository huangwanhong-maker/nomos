// Where travellers go and the way they take, in one equilibrium: the doubly constrained flows at the network's congested
// costs, assigned to the network by Wardrop's principle at the times those flows make.
//
// THE PROBLEM (Evans 1976, "Derivation and analysis of some models for combining trip distribution and assignment", Transp.
// Res. 10, 37–57; Florian, Nguyen and Ferland 1975). Minimise
//   Z(T, x) = Σ_a ∫_0^{x_a} t_a(s) ds + (1/β) Σ_ij T_ij (ln T_ij − 1)
// over the tables with margins O and D and the arc flows x that some assignment of T to paths makes. Its optimality
// conditions are the two models at once: T is the doubly constrained table at the least path times c_ij(x), and x is
// Wardrop's equilibrium for T. Z is strictly convex in T and convex in x, so the table is unique and so are the arc flows
// when the times are strictly increasing. Distributing at free-flow times and then assigning, or alternating the two, is
// not this equilibrium unless it is iterated to it.
//
// THE ALGORITHM is Evans's partial linearisation: at (T, x), with times t = t(x) and least times c, the subproblem keeps the
// entropy and linearises the potential, min Σ_ij c_ij T'_ij + (1/β) Σ T' (ln T' − 1) over the margins, whose solution T* is
// the doubly constrained table at c; y is T* loaded all-or-nothing on the shortest paths. The step along (T* − T, y − x)
// minimises Z (its derivative is monotone; bisection). Because the subproblem minimises a function below Z (convexity of the
// potential), its value
//   L = Σ_a ∫_0^{x_a} t_a + Σ_ij c_ij T*_ij − Σ_a t_a x_a + (1/β) Σ T* (ln T* − 1)
// is a lower bound on the least Z, and the gap (Z − L) / Σ_a t_a x_a is zero exactly at the equilibrium. Evans proved the
// iteration converges; like Frank–Wolfe it does so sublinearly near the solution.
//
// Conventions. Origin i leaves from node `origins[i]` and destination j is node `destinations[j]`; an optional table is added
// to every least time (terminal times, the walk to and from the network); a pair that cannot be joined has an infinite cost
// and no flow, and margins that only such pairs could meet leave the balancing unconverged.

#pragma once

#include <nomos/games/congestion.hpp>
#include <nomos/interaction/gravity.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <span>
#include <vector>

namespace nomos::interaction {

struct CombinedOptions {
  double gap = 1e-6;
  int max_iterations = 10000;
  const Table* added = nullptr;   // added to every least time, origins by destinations
};

struct Combined {
  Table flows;                    // T
  Table costs;                    // c at the end: least times plus what is added
  std::vector<double> arc_flow, arc_time;
  double objective = 0, lower_bound = 0;
  double gap = 0;
  int iterations = 0;
  bool converged = false;
};

namespace detail {
// The shortest-path trees from every origin's node at the arcs' times, once an iteration: the least costs and the all-or-nothing
// loading are both read from them (between parallel arcs the cheapest, ties to the lower index, as games' loading takes it).
struct Trees {
  std::vector<std::size_t> origin_node;                    // distinct origin nodes
  std::vector<std::size_t> tree_of;                        // origin → its tree
  std::vector<std::vector<double>> distance;               // tree → node
  std::vector<std::vector<std::size_t>> arc_into;          // tree → node → the arc its shortest path enters by (none at the root)
  std::vector<std::vector<std::size_t>> order;             // tree → its reached nodes, deepest first: every node before its parent
};

inline Trees trees(const games::Network& g, const std::vector<double>& t, const std::vector<std::size_t>& origins) {
  namespace mc = morphe::combinatorics;
  constexpr std::size_t none = static_cast<std::size_t>(-1);
  std::map<std::pair<std::size_t, std::size_t>, std::size_t> best;
  for (std::size_t a = 0; a < g.arcs.size(); a++) {
    const auto key = std::make_pair(g.arcs[a].tail, g.arcs[a].head);
    const auto it = best.find(key);
    if (it == best.end() || t[a] < t[it->second]) best[key] = a;
  }
  mc::WeightedAdjacency adj(g.nodes);
  for (const auto& [key, a] : best) adj[key.first].push_back({key.second, t[a]});
  Trees out;
  std::map<std::size_t, std::size_t> seen;
  for (const std::size_t o : origins) {
    const auto it = seen.find(o);
    if (it != seen.end()) { out.tree_of.push_back(it->second); continue; }
    const std::size_t k = out.origin_node.size();
    seen.emplace(o, k);
    out.tree_of.push_back(k);
    out.origin_node.push_back(o);
    const auto sp = mc::dijkstra(adj, o);
    std::vector<std::size_t> into(g.nodes, none);
    for (std::size_t v = 0; v < g.nodes; v++)
      if (sp.parent[v] != mc::ShortestPaths::none) into[v] = best.at({sp.parent[v], v});
    // Depths in the tree, not distances, order the loading: an arc of no time leaves a node and its parent at one distance.
    std::vector<std::size_t> depth(g.nodes, none), chain, order;
    for (std::size_t v = 0; v < g.nodes; v++) {
      if (!std::isfinite(sp.distance[v])) continue;
      order.push_back(v);
      chain.clear();
      std::size_t q = v;
      while (depth[q] == none && into[q] != none) chain.push_back(q), q = g.arcs[into[q]].tail;
      if (depth[q] == none) depth[q] = 0;
      std::size_t d = depth[q];
      for (auto it = chain.rbegin(); it != chain.rend(); ++it) depth[*it] = ++d;
    }
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return depth[a] > depth[b]; });
    out.distance.push_back(sp.distance);
    out.arc_into.push_back(std::move(into));
    out.order.push_back(std::move(order));
  }
  return out;
}

inline Table least_costs(const Trees& tr, const std::vector<std::size_t>& destinations, const Table* added) {
  Table c(tr.tree_of.size(), destinations.size());
  for (std::size_t i = 0; i < c.rows; i++)
    for (std::size_t j = 0; j < c.cols; j++) c(i, j) = tr.distance[tr.tree_of[i]][destinations[j]] + (added ? (*added)(i, j) : 0.0);
  return c;
}

// T loaded all or nothing on the trees' shortest paths: each origin's flows summed down its tree from the destinations.
inline std::vector<double> load(const games::Network& g, const Trees& tr, const Table& T, const std::vector<std::size_t>& destinations) {
  constexpr std::size_t none = static_cast<std::size_t>(-1);
  std::vector<double> x(g.arcs.size(), 0.0);
  for (std::size_t k = 0; k < tr.origin_node.size(); k++) {
    std::vector<double> at(g.nodes, 0.0);
    bool any = false;
    for (std::size_t i = 0; i < T.rows; i++) {
      if (tr.tree_of[i] != k) continue;
      for (std::size_t j = 0; j < T.cols; j++)
        if (T(i, j) > 0) at[destinations[j]] += T(i, j), any = true;
    }
    if (!any) continue;
    // The nodes deepest first: each passes what reaches it to the arc it is entered by, before its parent is visited.
    for (const std::size_t v : tr.order[k]) {
      if (at[v] == 0) continue;
      const std::size_t a = tr.arc_into[k][v];
      if (a == none) continue;
      x[a] += at[v];
      at[g.arcs[a].tail] += at[v];
    }
  }
  return x;
}

inline double entropy_term(const Table& T, double beta) {
  double s = 0;
  for (const double x : T.v)
    if (x > 0) s += x * (std::log(x) - 1);
  return s / beta;
}
}  // namespace detail

inline Combined combined(const games::Network& network, const std::vector<std::size_t>& origins, const std::vector<std::size_t>& destinations,
                         std::span<const double> outflow, std::span<const double> inflow, double beta, const CombinedOptions& options = {}) {
  detail::require(origins.size() == outflow.size() && destinations.size() == inflow.size(), "a node and an outflow per origin, a node and an inflow per destination");
  detail::require(beta > 0, "a positive β");
  const std::size_t A = network.arcs.size();
  if (options.added) detail::require(options.added->rows == origins.size() && options.added->cols == destinations.size(), "the added costs have the table's shape");
  Combined out;
  // From the doubly constrained table at free flow, loaded all or nothing.
  std::vector<double> x(A, 0.0), t = games::times(network, x);
  auto tr = detail::trees(network, t, origins);
  Table c = detail::least_costs(tr, destinations, options.added);
  Table T = doubly_constrained(outflow, inflow, c, beta).flows;
  x = detail::load(network, tr, T, destinations);
  for (out.iterations = 0;; out.iterations++) {
    t = games::times(network, x);
    tr = detail::trees(network, t, origins);
    c = detail::least_costs(tr, destinations, options.added);
    const Table S = doubly_constrained(outflow, inflow, c, beta).flows;
    // The added costs are linear in T and enter Z through Σ T·added; the potential's part is over the arcs.
    double added_T = 0;
    if (options.added)
      for (std::size_t k = 0; k < T.v.size(); k++)
        if (T.v[k] > 0) added_T += T.v[k] * options.added->v[k];
    struct { std::vector<double> flow; } y{detail::load(network, tr, S, destinations)};
    double used = 0, potential = 0;
    for (std::size_t a = 0; a < A; a++) used += t[a] * x[a], potential += games::arc_integral(network.arcs[a], x[a]);
    const double z = potential + added_T + detail::entropy_term(T, beta);
    double cs = 0;
    for (std::size_t k = 0; k < S.v.size(); k++)
      if (S.v[k] > 0) cs += c.v[k] * S.v[k];   // Σ c T*, the added costs included
    const double lower = potential + cs - used + detail::entropy_term(S, beta);
    out.objective = z;
    out.lower_bound = lower;
    const double scale = used + added_T;
    out.gap = scale > 0 ? (z - lower) / scale : 0.0;
    if (out.gap <= options.gap) { out.converged = true; break; }
    if (out.iterations >= options.max_iterations) break;
    // The step minimising Z along (S − T, y − x).
    const auto derivative = [&](double l) {
      double s = 0;
      for (std::size_t a = 0; a < A; a++) {
        const double d = y.flow[a] - x[a];
        if (d != 0) s += d * games::arc_time(network.arcs[a], x[a] + l * d);
      }
      for (std::size_t k = 0; k < T.v.size(); k++) {
        const double d = S.v[k] - T.v[k];
        if (d == 0) continue;
        const double at = T.v[k] + l * d;
        s += d * ((at > 0 ? std::log(at) / beta : -1e300) + (options.added ? options.added->v[k] : 0.0));
      }
      return s;
    };
    double lo = 0, hi = 1, l = 1;
    if (derivative(1) > 0) {
      for (int k = 0; k < 64; k++) {
        const double mid = 0.5 * (lo + hi);
        (derivative(mid) > 0 ? hi : lo) = mid;
      }
      l = 0.5 * (lo + hi);
    }
    for (std::size_t a = 0; a < A; a++) x[a] += l * (y.flow[a] - x[a]);
    for (std::size_t k = 0; k < T.v.size(); k++) T.v[k] += l * (S.v[k] - T.v[k]);
  }
  out.flows = T;
  out.costs = c;
  out.arc_flow = x;
  out.arc_time = t;
  return out;
}

}  // namespace nomos::interaction
