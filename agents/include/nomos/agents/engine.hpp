// Agent-based models: schedules, reproducible draws, the event log, and Gillespie's exact simulation of models in
// continuous time.
//
// An agent-based model is a Markov process on configurations ξ = (the agents' states, the environment's) generated
// by local rules, each reading an agent's neighbourhood and changing that agent (and what it touches). How the rules
// compose is the schedule, and it is part of the model (Huberman and Glance 1993, "Evolutionary games and computer
// simulations", PNAS 90, 7716–7718, show that synchronous and asynchronous updating of one model can behave
// differently):
// - synchronous: every agent reads ξ_t and writes its own part of ξ_{t+1} (a product of kernels; order-free, and so
//   run in parallel);
// - random sequential: in each tick the agents are visited in an order drawn afresh, each reading what the earlier
//   ones wrote (a composition of kernels);
// - asynchronous in continuous time: rule r fires at agent i at rate λ_r(ξ), the generator
//   L = Σ_i Σ_r λ_r(ξ) (ξ^{i,r} − ξ); simulated exactly by Gillespie's direct method (Gillespie 1977, "Exact
//   stochastic simulation of coupled chemical reactions", J. Phys. Chem. 81, 2340–2361): the waiting time is
//   exponential with the total rate, the event drawn with probability proportional to its rate. Its law is the
//   master equation ṗ = pL (physis's statistical_mechanics), against which small models are checked exactly.
//
// Draws. Every draw comes from a calx stream keyed by what it is for: agent i's draws at tick t from
// (seed, "agnt", i, t), the visiting order of tick t from (seed, "ordr", t, 0), Gillespie's from (seed, "gill",
// 0, 0). A synchronous run is therefore the same on any number of threads and in any order of evaluation; a random
// sequential or Gillespie run depends on the seed alone.
//
// The event log. A rule that does something records an Event: when, which agent, which rule, and what it touched
// (agents, cells, couplings: the caller's numbering). A domain reads the log into its own record of events.
//
// The model's contract (discrete time):
//   struct Model {
//     using State = ...;
//     std::size_t agents(const State&) const;
//     // Agent i's rule, reading `read` and writing (its own part of) `write`; in random sequential order the two
//     // are one state.
//     void update(std::size_t i, const State& read, State& write, calx::random::Stream&, Log*) const;
//   };
// and in continuous time:
//   struct JumpModel {
//     using State = ...;
//     std::size_t agents(const State&) const;
//     std::size_t rules() const;
//     double rate(const State&, std::size_t i, std::size_t r) const;
//     // Fires rule r at agent i; appends to `affected` the agents whose rates may have changed (i itself included).
//     void fire(State&, std::size_t i, std::size_t r, calx::random::Stream&, std::vector<std::size_t>& affected) const;
//   };

#pragma once

#include <calx/random/stream.hpp>
#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace nomos::agents {

enum class Schedule { synchronous, random_sequential };

struct Event {
  double time = 0;
  std::size_t agent = 0;
  std::uint32_t rule = 0;
  std::vector<std::size_t> touched;
};

struct Log {
  std::vector<Event> events;
  void add(double time, std::size_t agent, std::uint32_t rule, std::vector<std::size_t> touched = {}) { events.push_back({time, agent, rule, std::move(touched)}); }
};

constexpr std::uint32_t agent_stream = 0x61676e74u;  // "agnt"
constexpr std::uint32_t order_stream = 0x6f726472u;  // "ordr"
constexpr std::uint32_t jump_stream = 0x67696c6cu;   // "gill"

// One tick of a discrete-time model. Synchronous ticks may run on a pool (each agent writing only its own part);
// their log is then kept per agent and joined in the agents' order, so that it too does not depend on the threads.
template <class Model>
void tick(const Model& M, typename Model::State& state, Schedule schedule, std::uint64_t seed, std::uint32_t t, Log* log = nullptr,
          calx::scheduler::ThreadPool* pool = nullptr) {
  const std::size_t n = M.agents(state);
  if (schedule == Schedule::synchronous) {
    const typename Model::State read = state;
    if (pool && n > 1) {
      std::vector<Log> part(log ? n : 0);
      calx::scheduler::parallel_for(*pool, n, 64, [&](std::size_t begin, std::size_t end, std::size_t) {
        for (std::size_t i = begin; i < end; i++) {
          calx::random::Stream rng(seed, agent_stream, static_cast<std::uint32_t>(i), t);
          M.update(i, read, state, rng, log ? &part[i] : nullptr);
        }
      });
      if (log)
        for (auto& p : part)
          for (auto& e : p.events) log->events.push_back(std::move(e));
    } else {
      for (std::size_t i = 0; i < n; i++) {
        calx::random::Stream rng(seed, agent_stream, static_cast<std::uint32_t>(i), t);
        M.update(i, read, state, rng, log);
      }
    }
    return;
  }
  // Random sequential: a uniformly random order (Fisher and Yates, by the tick's own stream).
  std::vector<std::size_t> order(n);
  std::iota(order.begin(), order.end(), 0);
  calx::random::Stream ord(seed, order_stream, t, 0);
  for (std::size_t k = n; k > 1; k--) std::swap(order[k - 1], order[static_cast<std::size_t>(ord.below(k))]);
  for (const std::size_t i : order) {
    calx::random::Stream rng(seed, agent_stream, static_cast<std::uint32_t>(i), t);
    M.update(i, state, state, rng, log);
  }
}

// ---- continuous time ----------------------------------------------------------------------------------------------

// A Fenwick tree of non-negative rates: their total, a point's update and the index at a cumulative mass, each
// O(log n) (Fenwick 1994, "A new data structure for cumulative frequency tables", Softw. Pract. Exp. 24, 327–336).
class RateTree {
 public:
  explicit RateTree(std::size_t n = 0) : n_(n), tree_(n + 1, 0.0), value_(n, 0.0) {}
  std::size_t size() const { return n_; }
  double operator[](std::size_t i) const { return value_[i]; }
  void set(std::size_t i, double rate) {
    if (!(rate >= 0) || !std::isfinite(rate)) throw std::invalid_argument("nomos::agents: rates are non-negative and finite");
    const double d = rate - value_[i];
    value_[i] = rate;
    for (std::size_t k = i + 1; k <= n_; k += k & (~k + 1)) tree_[k] += d;
  }
  double total() const {
    double t = 0;
    for (std::size_t k = n_; k > 0; k -= k & (~k + 1)) t += tree_[k];
    return t;
  }
  // The least i with Σ_{j ≤ i} rate_j > u (0 ≤ u < total), skipping rates of zero.
  std::size_t find(double u) const {
    std::size_t pos = 0, step = 1;
    while (step * 2 <= n_) step *= 2;
    for (; step > 0; step /= 2)
      if (pos + step <= n_ && tree_[pos + step] <= u) pos += step, u -= tree_[pos];
    std::size_t i = std::min(pos, n_ - 1);
    while (value_[i] == 0 && i + 1 < n_) i++;  // rounding at a boundary
    return i;
  }

 private:
  std::size_t n_;
  std::vector<double> tree_, value_;
};

struct JumpRun {
  double time = 0;
  long events = 0;
  bool absorbed = false;  // every rate zero before the end
};

// Gillespie's direct method from `state` at time `t0` until `t_end`, with an optional event log; the rates of the
// agents a firing affects are recomputed, the rest kept.
template <class Model>
JumpRun gillespie(const Model& M, typename Model::State& state, double t0, double t_end, std::uint64_t seed, Log* log = nullptr, long max_events = -1) {
  const std::size_t n = M.agents(state), R = M.rules();
  RateTree rates(n * R);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t r = 0; r < R; r++) rates.set(i * R + r, M.rate(state, i, r));
  calx::random::Stream rng(seed, jump_stream, 0, 0);
  JumpRun out;
  out.time = t0;
  std::vector<std::size_t> affected;
  while (max_events < 0 || out.events < max_events) {
    const double total = rates.total();
    if (!(total > 0)) { out.absorbed = true; break; }
    const double dt = rng.exponential() / total;
    if (out.time + dt > t_end) break;
    out.time += dt;
    const std::size_t k = rates.find(rng.uniform() * total);
    const std::size_t i = k / R, r = k % R;
    affected.clear();
    M.fire(state, i, r, rng, affected);
    if (log) log->add(out.time, i, static_cast<std::uint32_t>(r));
    out.events++;
    for (const std::size_t j : affected)
      for (std::size_t q = 0; q < R; q++) rates.set(j * R + q, M.rate(state, j, q));
  }
  // The process stands in its state until t_end (the next event would come after it, or none ever will).
  if (max_events < 0 || out.events < max_events) out.time = t_end;
  return out;
}

// ---- ensembles ----------------------------------------------------------------------------------------------------

// The mean and standard error of a scalar over independent replicates (replicate k run with seed + k).
struct Summary {
  double mean = 0, error = 0;
  int runs = 0;
};

template <class Run>
Summary replicate(int runs, std::uint64_t seed, Run&& run) {
  double sum = 0, sq = 0;
  for (int k = 0; k < runs; k++) {
    const double x = run(seed + static_cast<std::uint64_t>(k));
    sum += x, sq += x * x;
  }
  Summary s;
  s.runs = runs;
  s.mean = sum / runs;
  s.error = runs > 1 ? std::sqrt(std::max(0.0, (sq - sum * sum / runs) / (runs - 1)) / runs) : 0;
  return s;
}

// ---- the description ----------------------------------------------------------------------------------------------

// A model's description by the ODD protocol (Grimm et al. 2006, 2010, 2020), with the decision section of ODD+D
// (Müller et al. 2013): kept beside the model, not as prose apart from it.
struct Description {
  std::string purpose, entities, scales, process, design_concepts, decisions, initialisation, input, submodels;
  std::string schedule;  // which of the schedules, and why
};

}  // namespace nomos::agents
