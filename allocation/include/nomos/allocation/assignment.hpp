// Who first, when goods are indivisible and there is no money: serial dictatorship, the random serial dictatorship, the
// probabilistic serial rule, picking sequences, and a random assignment carried out as a lottery over deterministic ones.
//
// Agents have strict preferences over objects (best first), each object a number of copies; an agent receives at most one
// object (unit demand) unless a picking sequence says otherwise.
//
//   serial_dictatorship    the agents in a declared order, each taking its most preferred object left (strategy-proof and
//                          Pareto efficient; the order is the whole of its justice)
//   random_serial          the probabilities of the serial dictatorship over a uniformly random order: exact over all n!
//                          orders for n ≤ 8, else by sampling with its error; computing them exactly is #P-hard in general
//                          (Aziz, Brandt and Brill 2013, Theor. Comput. Sci. 507, as the survey reports it)
//   probabilistic_serial   every agent eats its most preferred object still available at unit speed for one unit of time;
//                          an agent's share of an object is what it ate of it (Bogomolnaia and Moulin 2001, J. Econ. Theory
//                          100, 295–328): ordinally efficient and envy-free in the sense of stochastic dominance, but not
//                          strategy-proof, and no rule has all three for n ≥ 4 (their impossibility)
//   picking_sequence       agents pick in a declared sequence, each its most valued item left (Brams and Taylor; Bouveret
//                          and Lang 2011): round-robin is the sequence 0, 1, …, n − 1 repeated, which is envy-free up to one
//                          good for additive valuations
//   birkhoff               a bistochastic matrix (a random assignment of n agents to n objects) written as a lottery over
//                          permutation matrices (Birkhoff 1946; von Neumann 1953): the chances made into draws. Copies of an
//                          object, and agents left unserved, are made square by expanding the copies and adding a dummy.
#pragma once

#include <calx/random/stream.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace nomos::allocation {

using Preferences = std::vector<std::vector<std::size_t>>;   // for each agent, the objects it accepts, best first
inline constexpr std::size_t none = static_cast<std::size_t>(-1);

inline std::vector<std::size_t> serial_dictatorship(const Preferences& prefs, const std::vector<std::size_t>& copies, const std::vector<std::size_t>& order) {
  std::vector<std::size_t> left = copies, got(prefs.size(), none);
  for (std::size_t a : order) {
    if (a >= prefs.size()) throw std::invalid_argument("nomos::allocation::serial_dictatorship: an order of the agents");
    for (std::size_t o : prefs[a])
      if (o < left.size() && left[o] > 0) {
        got[a] = o, left[o]--;
        break;
      }
  }
  return got;
}

struct RandomAssignment {
  std::vector<std::vector<double>> p;   // p[agent][object]
  bool exact = true;
  std::size_t samples = 0;              // when not exact
};

inline RandomAssignment random_serial(const Preferences& prefs, const std::vector<std::size_t>& copies, calx::random::Stream* rng = nullptr, std::size_t samples = 20000) {
  const std::size_t n = prefs.size(), m = copies.size();
  RandomAssignment r;
  r.p.assign(n, std::vector<double>(m, 0));
  std::vector<std::size_t> order(n);
  std::iota(order.begin(), order.end(), 0);
  if (n <= 8) {
    double count = 0;
    do {
      const auto got = serial_dictatorship(prefs, copies, order);
      for (std::size_t a = 0; a < n; a++)
        if (got[a] != none) r.p[a][got[a]] += 1;
      count += 1;
    } while (std::next_permutation(order.begin(), order.end()));
    for (auto& row : r.p)
      for (double& v : row) v /= count;
    return r;
  }
  if (!rng) throw std::invalid_argument("nomos::allocation::random_serial: more than 8 agents need a random stream to sample orders");
  r.exact = false, r.samples = samples;
  for (std::size_t s = 0; s < samples; s++) {
    for (std::size_t i = n; i > 1; i--) std::swap(order[i - 1], order[rng->below(i)]);
    const auto got = serial_dictatorship(prefs, copies, order);
    for (std::size_t a = 0; a < n; a++)
      if (got[a] != none) r.p[a][got[a]] += 1;
  }
  for (auto& row : r.p)
    for (double& v : row) v /= static_cast<double>(samples);
  return r;
}

inline RandomAssignment probabilistic_serial(const Preferences& prefs, const std::vector<std::size_t>& copies) {
  const std::size_t n = prefs.size(), m = copies.size();
  RandomAssignment r;
  r.p.assign(n, std::vector<double>(m, 0));
  std::vector<double> supply(m), eaten(n, 0);
  for (std::size_t o = 0; o < m; o++) supply[o] = static_cast<double>(copies[o]);
  double t = 0;
  for (int phase = 0; phase < static_cast<int>(n + m + 2) && t < 1 - 1e-12; phase++) {
    // Each agent eats its best object with supply left; agents with none left stop.
    std::vector<std::size_t> eating(n, none);
    std::vector<double> eaters(m, 0);
    for (std::size_t a = 0; a < n; a++)
      for (std::size_t o : prefs[a])
        if (o < m && supply[o] > 1e-12) {
          eating[a] = o, eaters[o] += 1;
          break;
        }
    double step = 1 - t;
    bool anyone = false;
    for (std::size_t o = 0; o < m; o++)
      if (eaters[o] > 0) step = std::min(step, supply[o] / eaters[o]), anyone = true;
    if (!anyone) break;
    for (std::size_t a = 0; a < n; a++)
      if (eating[a] != none) r.p[a][eating[a]] += step, supply[eating[a]] -= step, eaten[a] += step;
    for (double& s : supply)
      if (s < 1e-12) s = 0;
    t += step;
  }
  return r;
}

// Agents pick in the sequence given, each its most valued item left (values[agent][item], additive; ties to the lower index).
inline std::vector<std::vector<std::size_t>> picking_sequence(const std::vector<std::vector<double>>& values, const std::vector<std::size_t>& sequence) {
  const std::size_t n = values.size(), g = n ? values[0].size() : 0;
  std::vector<bool> taken(g, false);
  std::vector<std::vector<std::size_t>> bundles(n);
  for (std::size_t a : sequence) {
    if (a >= n) throw std::invalid_argument("nomos::allocation::picking_sequence: a sequence of agents");
    std::size_t best = none;
    for (std::size_t j = 0; j < g; j++)
      if (!taken[j] && (best == none || values[a][j] > values[a][best])) best = j;
    if (best == none) break;
    taken[best] = true;
    bundles[a].push_back(best);
  }
  return bundles;
}

inline std::vector<std::size_t> round_robin_sequence(std::size_t agents, std::size_t items) {
  std::vector<std::size_t> s(items);
  for (std::size_t k = 0; k < items; k++) s[k] = k % agents;
  return s;
}

// ---- a random assignment as a lottery ---------------------------------------------------------------------------

struct Draw {
  double weight = 0;
  std::vector<std::size_t> object_of;   // for each row (agent), its column (object)
};

namespace detail {

// A perfect matching on the positive entries of a square matrix, by augmenting paths (Kuhn 1955).
inline bool perfect_matching(const std::vector<std::vector<double>>& M, std::vector<std::size_t>& col_of, double eps) {
  const std::size_t n = M.size();
  std::vector<std::size_t> row_of(n, none);
  col_of.assign(n, none);
  for (std::size_t r = 0; r < n; r++) {
    std::vector<bool> seen(n, false);
    std::function<bool(std::size_t)> augment = [&](std::size_t u) {
      for (std::size_t c = 0; c < n; c++)
        if (M[u][c] > eps && !seen[c]) {
          seen[c] = true;
          if (row_of[c] == none || augment(row_of[c])) {
            row_of[c] = u, col_of[u] = c;
            return true;
          }
        }
      return false;
    };
    if (!augment(r)) return false;
  }
  return true;
}

}  // namespace detail

// A bistochastic n × n matrix as a lottery over permutations whose weights sum to one (within `eps`).
inline std::vector<Draw> birkhoff(std::vector<std::vector<double>> M, double eps = 1e-9) {
  const std::size_t n = M.size();
  for (const auto& row : M) {
    if (row.size() != n) throw std::invalid_argument("nomos::allocation::birkhoff: a square matrix");
    double s = 0;
    for (double v : row) {
      if (v < -eps) throw std::invalid_argument("nomos::allocation::birkhoff: entries are not negative");
      s += v;
    }
    if (std::abs(s - 1) > 1e-6) throw std::invalid_argument("nomos::allocation::birkhoff: each row sums to one");
  }
  std::vector<Draw> out;
  double left = 1;
  for (std::size_t k = 0; k < n * n + 1 && left > eps; k++) {
    std::vector<std::size_t> col;
    if (!detail::perfect_matching(M, col, eps)) break;
    double w = 1;
    for (std::size_t r = 0; r < n; r++) w = std::min(w, M[r][col[r]]);
    for (std::size_t r = 0; r < n; r++) M[r][col[r]] -= w;
    out.push_back({w, col});
    left -= w;
  }
  return out;
}

// A random assignment of agents to objects with copies (rows sum to at most one, columns to at most the copies) made square:
// each copy a column, a dummy object for what an agent does not receive, dummy agents for copies no one receives. Returns
// the square matrix and, for each of its columns, the object it is a copy of (`none` for the dummy).
struct Squared {
  std::vector<std::vector<double>> matrix;
  std::vector<std::size_t> object_of_column;
  std::size_t agents = 0;
};

inline Squared square(const RandomAssignment& r, const std::vector<std::size_t>& copies) {
  const std::size_t n = r.p.size(), m = copies.size();
  std::size_t total = 0;
  for (std::size_t c : copies) total += c;
  const std::size_t N = std::max(n, total) + n;   // copies, then dummy columns for the unserved
  Squared s;
  s.agents = n;
  s.matrix.assign(N, std::vector<double>(N, 0));
  std::vector<std::size_t> first(m);
  std::size_t col = 0;
  for (std::size_t o = 0; o < m; o++) {
    first[o] = col;
    for (std::size_t k = 0; k < copies[o]; k++) s.object_of_column.push_back(o), col++;
  }
  while (s.object_of_column.size() < N) s.object_of_column.push_back(none);
  const std::size_t dummies = N - col;
  // Real agents: their probability of an object spread evenly over its copies; what they lack goes to the dummy columns.
  std::vector<double> col_load(N, 0);
  for (std::size_t a = 0; a < n; a++) {
    double got = 0;
    for (std::size_t o = 0; o < m; o++) {
      if (copies[o] == 0) continue;
      for (std::size_t k = 0; k < copies[o]; k++) s.matrix[a][first[o] + k] = r.p[a][o] / static_cast<double>(copies[o]);
      got += r.p[a][o];
    }
    for (std::size_t d = col; d < N; d++) s.matrix[a][d] = (1 - got) / static_cast<double>(dummies);
  }
  // Dummy agents take up what is left in every column, so that rows and columns all sum to one.
  for (std::size_t c = 0; c < N; c++) {
    double used = 0;
    for (std::size_t a = 0; a < n; a++) used += s.matrix[a][c];
    col_load[c] = 1 - used;
  }
  const std::size_t fill = N - n;
  for (std::size_t a = n; a < N; a++)
    for (std::size_t c = 0; c < N; c++) s.matrix[a][c] = col_load[c] / static_cast<double>(fill);
  return s;
}

// A draw from a lottery, by a random stream: the deterministic assignment chosen and its index.
inline std::size_t draw_index(const std::vector<Draw>& lottery, calx::random::Stream& rng) {
  double u = rng.uniform(), acc = 0, total = 0;
  for (const auto& d : lottery) total += d.weight;
  u *= total;
  for (std::size_t k = 0; k < lottery.size(); k++) {
    acc += lottery[k].weight;
    if (u <= acc) return k;
  }
  return lottery.empty() ? none : lottery.size() - 1;
}

}  // namespace nomos::allocation
