// Routing games with congestion: travellers between origins and destinations of a network each choose a path, every
// arc's time rising with the flow on it.
//
// THE GAME. A directed network with arcs a, each with a time t_a(x_a) that rises with its flow x_a, and demands T_w of
// travellers between pairs w = (o, d). Travellers are a continuum (a nonatomic game): one traveller's choice changes no
// arc's time. Path flows h_p ≥ 0 sum over the paths of each pair to its demand; an arc's flow is the sum of the flows of
// the paths through it.
//
// EQUILIBRIUM. Wardrop's first principle (Wardrop 1952, "Some theoretical aspects of road traffic research", Proc. Inst.
// Civ. Eng. II 1, 325–378): every path in use takes the least time of its pair, so that no traveller can shorten a
// journey by changing path. It is the Nash equilibrium of the nonatomic game. When every t_a is continuous and
// nondecreasing, the equilibrium arc flows are exactly the minimisers of Beckmann's potential
//   Z(x) = Σ_a ∫_0^{x_a} t_a(s) ds
// over the feasible flows (Beckmann, McGuire and Winsten 1956, Studies in the Economics of Transportation, Yale, ch. 3):
// its optimality conditions are Wardrop's. Z is convex, so equilibria exist; with strictly increasing times the arc flows
// are unique (path flows need not be). This is Rosenthal's (1973) potential, made continuous (Monderer and Shapley 1996).
//
// THE SYSTEM OPTIMUM. The flows least in the total time S(x) = Σ_a x_a t_a(x_a) are Wardrop's second principle. They are
// an equilibrium of the same game with each arc's time replaced by its marginal cost t_a + x_a t_a'(x_a), so that the same
// algorithm finds them. The price of anarchy S(equilibrium)/S(optimum) is at most 4/3 for affine times, attained by
// Pigou's two roads (Roughgarden and Tardos 2002, "How bad is selfish routing?", J. ACM 49, 236–259); Braess's network
// (Braess 1968, "Über ein Paradoxon aus der Verkehrsplanung", Unternehmensforschung 12, 258–268) shows that a road added
// can lengthen every journey.
//
// TIMES. t_a(x) = a + b x^p, b ≥ 0, p > 0: the Bureau of Public Roads' function t₀ (1 + α (x/C)^p) is a = t₀ and
// b = t₀ α / C^p (Bureau of Public Roads 1964, Traffic Assignment Manual, US Department of Commerce; α = 0.15 and p = 4 by
// convention), and affine times are p = 1.
//
// THE ALGORITHM. Frank and Wolfe's (1956, "An algorithm for quadratic programming", Naval Res. Logist. Q. 3, 95–110), as
// LeBlanc, Morlok and Pierskalla (1975, "An efficient approach to solving the road network equilibrium traffic assignment
// problem", Transp. Res. 9, 309–318) applied it: from the all-or-nothing loading at free flow, each iteration loads every
// demand on a shortest path at the current times (y, the minimiser of the potential's linearisation over the feasible
// flows), and moves x toward y by the step that minimises the potential along the segment (its derivative there is
// monotone, found by bisection). By convexity Z(x) + ∇Z(x)·(y − x) is a lower bound on the least potential; the
// relative gap (Σ_a x_a c_a − Σ_w T_w κ_w) / Σ_a x_a c_a, κ_w the least time of pair w at the current costs c, is zero
// exactly at equilibrium (Boyce, Ralevic-Dekic and Bar-Gera 2004, "Convergence of traffic assignment: how much is enough?",
// J. Transp. Eng. 130, 49–55). Frank–Wolfe converges sublinearly near the solution: a gap of 10⁻⁴, which practice asks of
// a forecast, is reached quickly, 10⁻⁸ only slowly. Path-based methods converge faster and are not yet written.
//
// Conventions. Shortest paths by morphe's Dijkstra; between parallel arcs the least cost, ties to the lower index. A
// demand that cannot reach its destination is an error.

#pragma once

#include <morphe/combinatorics/shortest_paths.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace nomos::games {

struct Arc {
  std::size_t tail = 0, head = 0;
  double a = 0;      // the time at no flow
  double b = 0;      // t(x) = a + b x^p
  double power = 1;
};

// The Bureau of Public Roads' arc: free-flow time t₀, capacity C, t = t₀ (1 + α (x/C)^p).
inline Arc bpr(std::size_t tail, std::size_t head, double free_time, double capacity, double alpha = 0.15, double power = 4) {
  if (!(capacity > 0) || !(free_time >= 0)) throw std::invalid_argument("nomos::games (congestion): a road has a capacity and a time");
  return {tail, head, free_time, free_time * alpha / std::pow(capacity, power), power};
}

inline double arc_time(const Arc& e, double x) { return e.a + (e.b != 0 ? e.b * std::pow(std::max(x, 0.0), e.power) : 0.0); }
// ∫_0^x t(s) ds.
inline double arc_integral(const Arc& e, double x) {
  x = std::max(x, 0.0);
  return e.a * x + (e.b != 0 ? e.b * std::pow(x, e.power + 1) / (e.power + 1) : 0.0);
}
// d/dx [x t(x)] = a + (p + 1) b x^p.
inline double arc_marginal(const Arc& e, double x) { return e.a + (e.b != 0 ? (e.power + 1) * e.b * std::pow(std::max(x, 0.0), e.power) : 0.0); }
struct Demand {
  std::size_t origin = 0, destination = 0;
  double trips = 0;
};

struct Network {
  std::size_t nodes = 0;
  std::vector<Arc> arcs;
  std::vector<Demand> demands;
};

enum class Objective { user, system };

namespace detail {
inline void congestion_require(bool ok, const std::string& what) {
  if (!ok) throw std::invalid_argument("nomos::games (congestion): " + what);
}

// Shortest-path trees from each origin at the given arc costs, and the all-or-nothing loading of the demands on them.
struct Loading {
  std::vector<double> flow;    // per arc
  double least = 0;            // Σ_w T_w κ_w
};

inline Loading all_or_nothing(const Network& g, const std::vector<double>& cost) {
  namespace mc = morphe::combinatorics;
  mc::WeightedAdjacency adj(g.nodes);
  // The cheapest arc between each ordered pair of nodes (ties to the lower index).
  std::map<std::pair<std::size_t, std::size_t>, std::size_t> best;
  for (std::size_t i = 0; i < g.arcs.size(); i++) {
    const auto key = std::make_pair(g.arcs[i].tail, g.arcs[i].head);
    const auto it = best.find(key);
    if (it == best.end() || cost[i] < cost[it->second]) best[key] = i;
  }
  for (const auto& [key, i] : best) adj[key.first].push_back({key.second, cost[i]});
  Loading out{std::vector<double>(g.arcs.size(), 0.0), 0.0};
  std::map<std::size_t, std::vector<const Demand*>> by_origin;
  for (const auto& d : g.demands)
    if (d.trips > 0) by_origin[d.origin].push_back(&d);
  for (const auto& [o, ds] : by_origin) {
    const auto sp = mc::dijkstra(adj, o);
    for (const Demand* d : ds) {
      congestion_require(std::isfinite(sp.distance[d->destination]), "a demand cannot reach its destination");
      out.least += d->trips * sp.distance[d->destination];
      for (std::size_t v = d->destination; sp.parent[v] != mc::ShortestPaths::none; v = sp.parent[v])
        out.flow[best.at({sp.parent[v], v})] += d->trips;
    }
  }
  return out;
}

inline std::vector<double> costs(const Network& g, const std::vector<double>& x, Objective o) {
  std::vector<double> c(g.arcs.size());
  for (std::size_t i = 0; i < c.size(); i++) c[i] = o == Objective::user ? arc_time(g.arcs[i], x[i]) : arc_marginal(g.arcs[i], x[i]);
  return c;
}

inline double objective(const Network& g, const std::vector<double>& x, Objective o) {
  double z = 0;
  for (std::size_t i = 0; i < x.size(); i++) z += o == Objective::user ? arc_integral(g.arcs[i], x[i]) : x[i] * arc_time(g.arcs[i], x[i]);
  return z;
}
}  // namespace detail

// The least time between every origin and destination node given, at the arcs' times (+inf where unreached): a skim.
inline std::vector<std::vector<double>> skims(const Network& g, const std::vector<double>& times, const std::vector<std::size_t>& origins,
                                              const std::vector<std::size_t>& destinations) {
  namespace mc = morphe::combinatorics;
  mc::WeightedAdjacency adj(g.nodes);
  for (std::size_t i = 0; i < g.arcs.size(); i++) adj[g.arcs[i].tail].push_back({g.arcs[i].head, times[i]});
  std::vector<std::vector<double>> out(origins.size(), std::vector<double>(destinations.size()));
  for (std::size_t r = 0; r < origins.size(); r++) {
    const auto sp = mc::dijkstra(adj, origins[r]);
    for (std::size_t c = 0; c < destinations.size(); c++) out[r][c] = sp.distance[destinations[c]];
  }
  return out;
}

inline std::vector<double> times(const Network& g, const std::vector<double>& flow) { return detail::costs(g, flow, Objective::user); }

inline double total_time(const Network& g, const std::vector<double>& flow) { return detail::objective(g, flow, Objective::system); }

inline double beckmann(const Network& g, const std::vector<double>& flow) { return detail::objective(g, flow, Objective::user); }

// The relative gap of arc flows at the given objective: (Σ x c − Σ_w T_w κ_w) / Σ x c.
inline double relative_gap(const Network& g, const std::vector<double>& flow, Objective o = Objective::user) {
  const auto c = detail::costs(g, flow, o);
  const auto y = detail::all_or_nothing(g, c);
  double used = 0;
  for (std::size_t i = 0; i < c.size(); i++) used += flow[i] * c[i];
  return used > 0 ? (used - y.least) / used : 0.0;
}

struct AssignmentOptions {
  Objective objective = Objective::user;
  double gap = 1e-6;              // stop when the relative gap falls below it
  int max_iterations = 100000;
  std::vector<double> start;      // feasible arc flows to start from; empty: all or nothing at free flow
};

struct Assignment {
  std::vector<double> flow, time;  // per arc
  double total_time = 0;           // Σ x t(x)
  double objective = 0;            // Beckmann's potential (user) or the total time (system)
  double lower_bound = -std::numeric_limits<double>::infinity();   // the best Frank–Wolfe bound on the least objective
  double gap = 0;                  // the relative gap at the end
  int iterations = 0;
  bool converged = false;
};

// Wardrop's equilibrium (Objective::user) or the system optimum (Objective::system), by Frank–Wolfe.
inline Assignment assign(const Network& g, const AssignmentOptions& options = {}) {
  for (const auto& e : g.arcs) {
    detail::congestion_require(e.tail < g.nodes && e.head < g.nodes, "an arc's ends are nodes");
    detail::congestion_require(e.a >= 0 && e.b >= 0 && e.power > 0, "an arc's time is nondecreasing and not negative");
  }
  for (const auto& d : g.demands) detail::congestion_require(d.origin < g.nodes && d.destination < g.nodes && d.trips >= 0, "a demand joins two nodes");
  const Objective o = options.objective;
  const std::size_t A = g.arcs.size();
  Assignment out;
  std::vector<double> x = options.start;
  if (x.empty()) x = detail::all_or_nothing(g, detail::costs(g, std::vector<double>(A, 0.0), o)).flow;
  detail::congestion_require(x.size() == A, "one starting flow per arc");
  for (out.iterations = 0;; out.iterations++) {
    const auto c = detail::costs(g, x, o);
    const auto y = detail::all_or_nothing(g, c);
    double used = 0, slope = 0;
    for (std::size_t i = 0; i < A; i++) used += x[i] * c[i], slope += c[i] * (y.flow[i] - x[i]);
    const double z = detail::objective(g, x, o);
    out.lower_bound = std::max(out.lower_bound, z + slope);
    out.gap = used > 0 ? (used - y.least) / used : 0.0;
    if (out.gap <= options.gap) { out.converged = true; break; }
    if (out.iterations >= options.max_iterations) break;
    // The step λ ∈ [0, 1] minimising the objective along x + λ (y − x): its derivative Σ c_a(x + λd) d_a is nondecreasing.
    const auto derivative = [&](double l) {
      double s = 0;
      for (std::size_t i = 0; i < A; i++) {
        const double d = y.flow[i] - x[i];
        if (d == 0) continue;
        const double at = x[i] + l * d;
        s += d * (o == Objective::user ? arc_time(g.arcs[i], at) : arc_marginal(g.arcs[i], at));
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
    for (std::size_t i = 0; i < A; i++) x[i] += l * (y.flow[i] - x[i]);
  }
  out.flow = x;
  out.time = times(g, x);
  out.total_time = total_time(g, x);
  out.objective = detail::objective(g, x, o);
  return out;
}

// S(equilibrium) / S(optimum).
inline double price_of_anarchy(const Network& g, double gap = 1e-9, int max_iterations = 1000000) {
  const auto ue = assign(g, {Objective::user, gap, max_iterations, {}});
  const auto so = assign(g, {Objective::system, gap, max_iterations, {}});
  return ue.total_time / so.total_time;
}

}  // namespace nomos::games
