// Checks of agent-based models (agents/engine.hpp, space.hpp, choice.hpp).
//
// - The rate tree: totals and the index at a cumulative mass, against a linear search, after random updates.
// - Gillespie's simulation against the master equation, exactly: an SIS model of 5 agents on the complete graph
//   (32 configurations), the law of the number infected at t = 1.5 from physis's evolve, against 40 000 runs
//   (Pearson's χ²).
// - Kurtz's mean field: the same model on 1000 agents, its fraction infected at t = 4 against the logistic solution
//   of di/dt = β i (1 − i) − δ i.
// - The voter model, random sequential on the complete graph of 20: the chance of consensus on an opinion is its
//   initial share (a martingale), 0.25, over 4000 runs.
// - Synchronous ticks: the same states on one thread and on four, and the same event log.
// - Schelling's model on a torus of 50 × 50, 2000 agents: from about ½ like neighbours, a threshold of ⅓ sorts them
//   to more than 0.7, everyone content; a threshold of 0 moves no one; ½ sorts more than ⅓.
// - Logit choice: the law, its ratios, the limits β → ∞ and β = 0, no overflow; draws against the law.
// - The plane's radius queries against brute force.

#include <nomos/agents/choice.hpp>
#include <nomos/agents/engine.hpp>
#include <nomos/agents/space.hpp>

#include <physis/statistical_mechanics/master_equation.hpp>

#include <calx/linalg/dense.hpp>
#include <calx/random/stream.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <string>
#include <vector>

namespace na = nomos::agents;

namespace {

// SIS on the complete graph: a susceptible agent is infected at rate β · (infected others) / (n − 1); an infected
// one recovers at rate δ.
struct Sis {
  using State = std::vector<int>;
  double beta = 2, delta = 1;
  std::size_t agents(const State& x) const { return x.size(); }
  std::size_t rules() const { return 2; }
  int infected(const State& x) const {
    int k = 0;
    for (const int v : x) k += v;
    return k;
  }
  double rate(const State& x, std::size_t i, std::size_t r) const {
    const double n = static_cast<double>(x.size());
    if (r == 0) return x[i] == 0 ? beta * (infected(x) - x[i]) / (n - 1) : 0.0;
    return x[i] == 1 ? delta : 0.0;
  }
  void fire(State& x, std::size_t i, std::size_t r, calx::random::Stream&, std::vector<std::size_t>& affected) const {
    x[i] = r == 0 ? 1 : 0;
    for (std::size_t j = 0; j < x.size(); j++) affected.push_back(j);  // every infection rate depends on the count
  }
};

// The same SIS on counts, its rates kept in O(1) per event: infections at β k (n − k)/(n − 1), recoveries at δ k.
// (For the mean-field check on 1000 agents; it is the same process as Sis by exchangeability.)
struct SisCount {
  using State = std::vector<int>;  // one entry: the number infected
  double beta = 2, delta = 1;
  std::size_t n = 1000;
  std::size_t agents(const State&) const { return 1; }
  std::size_t rules() const { return 2; }
  double rate(const State& x, std::size_t, std::size_t r) const {
    const double k = x[0], N = static_cast<double>(n);
    return r == 0 ? beta * k * (N - k) / (N - 1) : delta * k;
  }
  void fire(State& x, std::size_t, std::size_t r, calx::random::Stream&, std::vector<std::size_t>& affected) const {
    x[0] += r == 0 ? 1 : -1;
    affected.push_back(0);
  }
};

// The voter model on the complete graph: an agent adopts the opinion of a neighbour drawn uniformly.
struct Voter {
  using State = std::vector<int>;
  na::Graph graph;
  std::size_t agents(const State& x) const { return x.size(); }
  void update(std::size_t i, const State& read, State& write, calx::random::Stream& rng, na::Log*) const {
    const auto& nb = graph.neighbours(i);
    write[i] = read[nb[static_cast<std::size_t>(rng.below(nb.size()))]];
  }
};

// A noisy majority rule on a lattice, synchronous: each cell takes its neighbours' majority, or with probability ε a
// fair coin; an event logged when it changes.
struct Majority {
  using State = std::vector<int>;
  na::Lattice space;
  double noise = 0.1;
  std::size_t agents(const State& x) const { return x.size(); }
  void update(std::size_t i, const State& read, State& write, calx::random::Stream& rng, na::Log* log) const {
    int ones = 0;
    const auto nb = space.neighbours(i);
    for (const std::size_t m : nb) ones += read[m];
    int v = 2 * ones > static_cast<int>(nb.size()) ? 1 : 2 * ones < static_cast<int>(nb.size()) ? 0 : read[i];
    if (rng.uniform() < noise) v = rng.bernoulli(0.5) ? 1 : 0;
    write[i] = v;
    if (log && v != read[i]) log->add(0, i, 0, {i});
  }
};

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- the rate tree --------------------------------------------------------------------------------------------
  {
    calx::random::Stream rng(20261104, 0x74726565u, 0, 0);  // "tree"
    na::RateTree tree(37);
    std::vector<double> r(37, 0.0);
    bool ok = true;
    for (int round = 0; round < 2000; round++) {
      const std::size_t i = static_cast<std::size_t>(rng.below(37));
      r[i] = rng.bernoulli(0.2) ? 0.0 : rng.exponential();
      tree.set(i, r[i]);
      double total = 0;
      for (const double x : r) total += x;
      ok = ok && std::fabs(tree.total() - total) <= 1e-12 * (1 + total);
      if (total > 0) {
        const double u = rng.uniform() * total;
        double c = 0;
        std::size_t want = 0;
        for (std::size_t j = 0; j < r.size(); j++)
          if ((c += r[j]) > u) { want = j; break; }
        ok = ok && tree.find(u) == want;
      }
    }
    s.require("the rate tree: totals and the index at a cumulative mass agree with a linear search over 2000 updates", ok);
  }

  // ---- Gillespie against the master equation --------------------------------------------------------------------
  {
    const Sis M;
    const std::size_t n = 5, C = std::size_t{1} << n;
    calx::linalg::Matrix k(C);
    for (std::size_t c = 0; c < C; c++) {
      std::vector<int> x(n);
      for (std::size_t i = 0; i < n; i++) x[i] = (c >> i) & 1;
      for (std::size_t i = 0; i < n; i++)
        for (std::size_t r = 0; r < 2; r++) {
          const double rate = M.rate(x, i, r);
          if (rate > 0) k(c, c ^ (std::size_t{1} << i)) += rate;
        }
    }
    std::vector<double> p0(C, 0.0);
    p0[1] = 1;  // agent 0 infected
    const auto p = physis::statistical_mechanics::evolve(k, p0, 1.5);
    std::vector<double> exact(n + 1, 0.0), seen(n + 1, 0.0);
    for (std::size_t c = 0; c < C; c++) exact[static_cast<std::size_t>(__builtin_popcountll(c))] += p[c];
    const int runs = 40000;
    for (int run = 0; run < runs; run++) {
      std::vector<int> x(n, 0);
      x[0] = 1;
      na::gillespie(M, x, 0, 1.5, 20261104 + static_cast<std::uint64_t>(run));
      seen[static_cast<std::size_t>(M.infected(x))] += 1;
    }
    double chi2 = 0;
    int df = -1;
    for (std::size_t m = 0; m <= n; m++)
      if (exact[m] * runs > 5) chi2 += (seen[m] - runs * exact[m]) * (seen[m] - runs * exact[m]) / (runs * exact[m]), df++;
    // χ² on 5 degrees of freedom: its 99.9 % quantile is 20.5.
    s.require("Gillespie's SIS of 5 agents against the master equation's law of the number infected at t = 1.5 (Pearson's χ², " + std::to_string(df) + " degrees of freedom)",
              chi2 < 20.5, std::to_string(chi2));
  }

  // ---- Kurtz's mean field ---------------------------------------------------------------------------------------
  {
    const SisCount M;
    const double i0 = 0.1, istar = 1 - M.delta / M.beta, t = 4;
    const double ode = istar / (1 + (istar / i0 - 1) * std::exp(-(M.beta - M.delta) * t));
    const auto sum = na::replicate(40, 77, [&](std::uint64_t seed) {
      std::vector<int> x{static_cast<int>(i0 * M.n)};
      na::gillespie(M, x, 0, t, seed);
      return x[0] / static_cast<double>(M.n);
    });
    s.require("SIS on 1000 agents follows Kurtz's mean field: the fraction infected at t = 4 against the logistic solution " + std::to_string(ode),
              std::fabs(sum.mean - ode) <= 4 * sum.error + 3.0 / M.n, std::to_string(sum.mean) + " ± " + std::to_string(sum.error));
  }

  // ---- the voter model's martingale -----------------------------------------------------------------------------
  {
    const Voter M{na::Graph::complete(20)};
    const auto sum = na::replicate(4000, 1000, [&](std::uint64_t seed) {
      std::vector<int> x(20, 0);
      for (int i = 0; i < 5; i++) x[static_cast<std::size_t>(i)] = 1;
      for (std::uint32_t t = 0; t < 100000; t++) {
        int ones = 0;
        for (const int v : x) ones += v;
        if (ones == 0 || ones == 20) return ones == 20 ? 1.0 : 0.0;
        na::tick(M, x, na::Schedule::random_sequential, seed, t);
      }
      return -1.0;
    });
    s.require("the voter model reaches consensus on an opinion with the chance of its initial share, 0.25", std::fabs(sum.mean - 0.25) <= 4 * sum.error,
              std::to_string(sum.mean) + " ± " + std::to_string(sum.error));
  }

  // ---- synchronous ticks on any number of threads -----------------------------------------------------------------
  {
    const Majority M{{40, 40, true, na::Neighbourhood::moore}};
    std::vector<int> a(M.space.cells()), b;
    calx::random::Stream rng(20261104, 0x6d616a6fu, 0, 0);  // "majo"
    for (auto& v : a) v = rng.bernoulli(0.5) ? 1 : 0;
    b = a;
    na::Log la, lb;
    calx::scheduler::ThreadPool pool(4);
    for (std::uint32_t t = 0; t < 20; t++) na::tick(M, a, na::Schedule::synchronous, 9, t, &la), na::tick(M, b, na::Schedule::synchronous, 9, t, &lb, &pool);
    bool same_log = la.events.size() == lb.events.size();
    for (std::size_t e = 0; same_log && e < la.events.size(); e++) same_log = la.events[e].agent == lb.events[e].agent;
    s.require("synchronous ticks give the same states on one thread and on four", a == b);
    s.require("... and the same event log (" + std::to_string(la.events.size()) + " changes)", same_log);
  }

  // ---- Schelling --------------------------------------------------------------------------------------------------
  {
    const auto run = [&](double f, std::uint64_t seed, double* content) {
      na::Schelling M{{50, 50, true, na::Neighbourhood::moore}, f};
      auto x = M.place(2000, seed);
      for (std::uint32_t t = 0; t < 40; t++) na::tick(M, x, na::Schedule::random_sequential, seed, t);
      if (content) *content = M.share_content(x);
      return M.segregation(x);
    };
    na::Schelling M0{{50, 50, true, na::Neighbourhood::moore}, 1.0 / 3};
    const double before = M0.segregation(M0.place(2000, 5));
    double content = 0;
    const double third = run(1.0 / 3, 5, &content), half = run(0.5, 5, nullptr), none = run(0.0, 5, nullptr);
    s.require("Schelling, threshold ⅓: from " + std::to_string(before) + " like neighbours to more than 0.7, all content",
              before < 0.55 && third > 0.7 && content == 1.0, std::to_string(third) + ", content " + std::to_string(content));
    s.check("... a threshold of 0 moves no one", std::fabs(none - before), 1e-15);
    s.require("... and ½ sorts more than ⅓", half > third, std::to_string(half));
    s.require("its description fills ODD's sections", !M0.describe().purpose.empty() && !M0.describe().schedule.empty());
  }

  // ---- logit choice -----------------------------------------------------------------------------------------------
  {
    const std::vector<double> u{1.0, 2.0, 0.5};
    const auto p = na::logit(u, 1.5);
    double total = 0;
    for (const double x : p) total += x;
    s.check("logit: a law whose ratios are exp(β Δu)", std::max(std::fabs(total - 1), std::fabs(p[1] / p[0] - std::exp(1.5))), 1e-14);
    const auto big = na::logit(std::vector<double>{1000, 1001, 999}, 50), flat = na::logit(u, 0);
    s.require("... β → ∞ takes the best, β = 0 is uniform, and utilities of 1000 do not overflow",
              big[1] > 1 - 1e-15 && std::fabs(flat[2] - 1.0 / 3) < 1e-15);
    calx::random::Stream rng(20261104, 0x6c6f6769u, 0, 0);  // "logi"
    std::vector<double> seen(3, 0.0);
    const int n = 200000;
    for (int k = 0; k < n; k++) seen[na::choose(rng, u, 1.5)] += 1;
    double chi2 = 0;
    for (std::size_t j = 0; j < 3; j++) chi2 += (seen[j] - n * p[j]) * (seen[j] - n * p[j]) / (n * p[j]);
    s.require("... its draws follow its law (χ² on 2 degrees of freedom below 13.8)", chi2 < 13.8, std::to_string(chi2));
  }

  // ---- the plane ----------------------------------------------------------------------------------------------------
  {
    calx::random::Stream rng(20261104, 0x706c616eu, 0, 0);  // "plan"
    na::Plane plane(0, 0, 100, 60, 7);
    for (int k = 0; k < 3000; k++) plane.insert(rng.uniform(0, 100), rng.uniform(0, 60));
    bool ok = true;
    for (int q = 0; q < 200; q++) {
      const double x = rng.uniform(-5, 105), y = rng.uniform(-5, 65), r = rng.uniform(0, 15);
      std::vector<std::size_t> brute;
      for (std::size_t k = 0; k < plane.points.size(); k++) {
        const double dx = plane.points[k].first - x, dy = plane.points[k].second - y;
        if (dx * dx + dy * dy <= r * r) brute.push_back(k);
      }
      ok = ok && plane.within(x, y, r) == brute;
    }
    s.require("the plane's radius queries equal brute force at 200 queries over 3000 points", ok);
  }
  return s.done();
}
