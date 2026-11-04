// How agents choose: the multinomial logit, and Schelling's model of segregation.
//
// Logit choice (McFadden 1974, "Conditional logit analysis of qualitative choice behavior", in Zarembka (ed.),
// Frontiers in Econometrics, 105–142): among options of utilities u_j, option j with probability
// exp(β u_j) / Σ_k exp(β u_k), computed after subtracting max_k u_k so that no exponential overflows. β → ∞ chooses
// the best; β = 0 chooses uniformly. nomos's econometrics estimates the utilities; here they are given.
//
// Schelling's model (Schelling 1971, "Dynamic models of segregation", J. Math. Sociol. 1, 143–186): agents of two
// kinds on a lattice with vacant cells; an agent is content when at least the fraction f of its occupied neighbours
// are of its kind (an agent with no occupied neighbour is content), and a discontented agent moves to a vacant cell
// drawn uniformly among those where it would be content, or stays if there is none. Under random sequential
// updating, mild preferences (f near 1/3) give strongly sorted neighbourhoods: the canonical check that micro-motives
// need not resemble macro-patterns. Segregation is measured by the mean fraction of like neighbours.
//
// Kirman's model of herding (Kirman 1993, "Ants, rationality, and recruitment", Q. J. Econ. 108, 137–156): N agents
// each hold one of two views (two food sources, a market's optimism or pessimism, two suppliers); an agent changes its
// mind of itself at rate ε, or is recruited by another agent of the other view at rate δ times the share of the others
// who hold it: agent i switches at the rate ε + δ (others of the other view) / (N − 1). The number k holding the first
// view is a birth–death chain, with rates up (N − k)(ε + δ k / (N − 1)) and down k (ε + δ (N − k)/(N − 1)), whose
// stationary law is π_k ∝ Π_{j<k} up_j / down_{j+1} exactly. When ε < δ / (N − 1) it is bimodal: the market herds to
// one view, holds it, and switches, though no agent prefers either.

#pragma once

#include <nomos/agents/engine.hpp>
#include <nomos/agents/space.hpp>

#include <calx/random/categorical.hpp>
#include <calx/random/stream.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace nomos::agents {

inline std::vector<double> logit(std::span<const double> utility, double beta) {
  if (utility.empty()) throw std::invalid_argument("nomos::agents::logit: at least one option");
  const double top = *std::max_element(utility.begin(), utility.end());
  std::vector<double> p(utility.size());
  double total = 0;
  for (std::size_t j = 0; j < utility.size(); j++) total += (p[j] = std::exp(beta * (utility[j] - top)));
  for (double& x : p) x /= total;
  return p;
}

inline std::size_t choose(calx::random::Stream& rng, std::span<const double> utility, double beta) {
  const auto p = logit(utility, beta);
  return calx::random::categorical(rng, p);
}

// ---- Schelling --------------------------------------------------------------------------------------------------

struct Schelling {
  struct State {
    std::vector<int> cell;                 // −1 vacant, else the kind (0 or 1) of the agent there
    std::vector<std::size_t> where;        // agent → its cell
  };
  Lattice space;
  double threshold = 1.0 / 3;              // f

  std::size_t agents(const State& x) const { return x.where.size(); }

  // The fraction of an agent of kind k's occupied neighbours that are of its kind, were it at cell c (and its own
  // cell `self` counted as vacant); 1 when it would have no occupied neighbour.
  double like(const State& x, std::size_t c, int k, std::size_t self) const {
    int same = 0, other = 0;
    for (const std::size_t m : space.neighbours(c)) {
      if (m == self || x.cell[m] < 0) continue;
      (x.cell[m] == k ? same : other)++;
    }
    return same + other == 0 ? 1.0 : static_cast<double>(same) / (same + other);
  }
  bool content(const State& x, std::size_t agent) const {
    const std::size_t c = x.where[agent];
    return like(x, c, x.cell[c], c) >= threshold;
  }

  void update(std::size_t i, const State&, State& x, calx::random::Stream& rng, Log* log) const {
    const std::size_t c = x.where[i];
    const int k = x.cell[c];
    if (like(x, c, k, c) >= threshold) return;
    std::vector<std::size_t> good;
    for (std::size_t m = 0; m < x.cell.size(); m++)
      if (x.cell[m] < 0 && like(x, m, k, c) >= threshold) good.push_back(m);
    if (good.empty()) return;
    const std::size_t to = good[static_cast<std::size_t>(rng.below(good.size()))];
    x.cell[c] = -1, x.cell[to] = k, x.where[i] = to;
    if (log) log->add(0, i, 0, {c, to});
  }

  // A random placement of n agents, half of each kind, on the lattice.
  State place(std::size_t n, std::uint64_t seed) const {
    if (n > space.cells()) throw std::invalid_argument("nomos::agents::Schelling: more agents than cells");
    State x;
    x.cell.assign(space.cells(), -1);
    std::vector<std::size_t> cells(space.cells());
    for (std::size_t c = 0; c < cells.size(); c++) cells[c] = c;
    calx::random::Stream rng(seed, 0x706c6163u, 0, 0);  // "plac"
    for (std::size_t k = cells.size(); k > 1; k--) std::swap(cells[k - 1], cells[static_cast<std::size_t>(rng.below(k))]);
    for (std::size_t i = 0; i < n; i++) x.cell[cells[i]] = static_cast<int>(i % 2), x.where.push_back(cells[i]);
    return x;
  }

  // The mean fraction of like neighbours over agents with at least one occupied neighbour.
  double segregation(const State& x) const {
    double sum = 0;
    int n = 0;
    for (std::size_t i = 0; i < x.where.size(); i++) {
      const std::size_t c = x.where[i];
      int same = 0, other = 0;
      for (const std::size_t m : space.neighbours(c))
        if (x.cell[m] >= 0) (x.cell[m] == x.cell[c] ? same : other)++;
      if (same + other > 0) sum += static_cast<double>(same) / (same + other), n++;
    }
    return n ? sum / n : 0;
  }
  double share_content(const State& x) const {
    int n = 0;
    for (std::size_t i = 0; i < x.where.size(); i++) n += content(x, i);
    return static_cast<double>(n) / static_cast<double>(x.where.size());
  }

  Description describe() const {
    return {"How mild preferences about the kind of one's neighbours sort a population (Schelling 1971).",
            "Agents of two kinds; the cells of a lattice, vacant or occupied.",
            "One tick: every agent visited once; the lattice's cells, Moore's neighbourhood.",
            "A discontented agent moves to a vacant cell where it would be content, drawn uniformly; else it stays.",
            "Emergence: sorted neighbourhoods from mild preferences. Sensing: the kinds of the occupied neighbours.",
            "A threshold rule on the fraction of like neighbours; no planning.",
            "Agents placed uniformly at random, half of each kind.",
            "None.",
            "like(): the fraction of like occupied neighbours; segregation(): its mean over agents.",
            "Random sequential, as Schelling moved one agent at a time."};
  }
};

// ---- Kirman ------------------------------------------------------------------------------------------------------

struct Kirman {
  using State = std::vector<int>;  // each agent's view, 0 or 1
  double epsilon = 0.002, delta = 0.3;

  std::size_t agents(const State& x) const { return x.size(); }
  std::size_t rules() const { return 1; }
  static int ones(const State& x) {
    int k = 0;
    for (const int v : x) k += v;
    return k;
  }
  double rate(const State& x, std::size_t i, std::size_t) const {
    const double n = static_cast<double>(x.size());
    const int k = ones(x);
    const double other = x[i] == 1 ? n - k : k - 0.0;  // others of the other view (i itself is not of it)
    return epsilon + delta * other / (n - 1);
  }
  void fire(State& x, std::size_t i, std::size_t, calx::random::Stream&, std::vector<std::size_t>& affected) const {
    x[i] = 1 - x[i];
    for (std::size_t j = 0; j < x.size(); j++) affected.push_back(j);  // every rate depends on the count
  }

  // The exact stationary law of the number holding view 1, for N agents.
  std::vector<double> stationary(std::size_t N) const {
    const double n = static_cast<double>(N);
    std::vector<double> p(N + 1, 1.0);
    for (std::size_t k = 1; k <= N; k++) {
      const double j = static_cast<double>(k - 1);
      const double up = (n - j) * (epsilon + delta * j / (n - 1)), down = (j + 1) * (epsilon + delta * (n - j - 1) / (n - 1));
      p[k] = p[k - 1] * up / down;
    }
    double t = 0;
    for (const double v : p) t += v;
    for (double& v : p) v /= t;
    return p;
  }

  Description describe() const {
    return {"How a population herds between two views that no one prefers (Kirman 1993).",
            "Agents, each holding one of two views.",
            "Continuous time; no space: every agent may meet every other.",
            "An agent changes its mind of itself at rate ε, or is recruited at rate δ times the share of the others who think otherwise.",
            "Emergence: herding and sudden switches of the whole population. Interaction: recruitment by random meetings.",
            "A rate rule; no planning.",
            "Each agent's view drawn at random.",
            "None.",
            "stationary(): the exact law of the number holding view 1, a birth–death chain.",
            "Gillespie's direct method: the model is a jump process."};
  }
};

}  // namespace nomos::agents
