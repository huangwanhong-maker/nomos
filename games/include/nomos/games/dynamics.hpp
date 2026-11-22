// Population games and their evolutionary dynamics: revision protocols, the mean dynamics they define, their incentive properties, their
// Lyapunov functions, and the finite population's Markov process (decision 0017).
//
// A population of mass one plays n strategies; its state x lies on the simplex X, and a population game F gives each strategy's payoff
// F_i(x). A revision protocol ρ(π, x) is a nonnegative n × n matrix: an agent playing i who receives a revision opportunity switches to j
// at rate ρ_ij. The mean dynamic is ẋ_i = Σ_j x_j ρ_ji − x_i Σ_j ρ_ij (Sandholm 2015, eq. 13.3), and in a population of N agents each
// revising at the rate of an exponential clock, the state is a Markov jump process whose paths approach the mean dynamic's solutions over
// finite horizons as N grows (Kurtz's theorem; Sandholm 2015, Theorem 13.1). The protocols of Sandholm's Table 13.1, with F̂_i = F_i − F̄,
// F̄ = Σ x_i F_i:
//
//   replicator      ρ_ij = x_j [π_j − π_i]_+ (pairwise proportional imitation)        ẋ_i = x_i F̂_i
//   best response   ρ_i· = M(π), the uniform law on the best strategies                ẋ ∈ M(F(x)) − x
//   logit (η)       ρ_ij = exp(π_j/η) / Σ_k exp(π_k/η)                                 ẋ_i = exp(F_i/η)/Σ_k exp(F_k/η) − x_i
//   BNN             ρ_ij = [π_j − Σ_k x_k π_k]_+ (excess payoff)                       ẋ_i = [F̂_i]_+ − x_i Σ_j [F̂_j]_+
//   Smith           ρ_ij = [π_j − π_i]_+ (pairwise comparison)                         ẋ_i = Σ_j x_j [F_i − F_j]_+ − x_i Σ_j [F_j − F_i]_+
//
// and the Maynard Smith replicator ẋ_i = x_i F̂_i / F̄ (positive payoffs; one population: the replicator's orbits at another speed).
// Incentive properties (§13.4.2): positive correlation, V(x) ≠ 0 ⇒ V(x)′F(x) > 0; Nash stationarity, V(x) = 0 ⇔ x is a Nash equilibrium.
// The replicator satisfies the first and not the second (every vertex is a rest point); best response, BNN and Smith both; logit both
// only approximately (Table 13.3). In contractive games, (y − x)′(F(y) − F(x)) ≤ 0, each dynamic has a Lyapunov function (Table 13.4):
// the replicator's H(x) = Σ_{i ∈ supp x*} x*_i log(x*_i / x_i), BNN's Γ(x) = ½ Σ_i [F̂_i]_+², Smith's Ψ(x) = ½ Σ_i Σ_j x_i [F_j − F_i]_+².
// In a potential game F = ∇f, f is a strict Lyapunov function of every dynamic with positive correlation (Theorem 13.6).
//
// Accuracy: trajectories by the classical Runge–Kutta method with a declared step, clipped to the simplex and renormalised after each
// step (the dynamics keep the simplex; clipping removes rounding only). Best response is a differential inclusion: its trajectory takes
// the uniform law on ties, which is one selection of the inclusion. The finite population is simulated exactly (Gillespie).
//
// Source: W. H. Sandholm, "Population games and deterministic evolutionary dynamics", Handbook of Game Theory with Economic Applications
// 4 (2015), ch. 13, 703–778 (doi 10.1016/B978-0-444-53766-9.00013-6); W. H. Sandholm, Population Games and Evolutionary Dynamics, MIT
// Press 2010; P. D. Taylor and L. B. Jonker, "Evolutionarily stable strategies and game dynamics", Math. Biosci. 40 (1978) 145–156;
// M. J. Smith, "The stability of a dynamic model of traffic assignment", Transp. Sci. 18 (1984) 245–252; J. Hofbauer and W. H. Sandholm,
// "Survival of dominated strategies under evolutionary dynamics", Theor. Econ. 6 (2011) 341–377.

#pragma once

#include <calx/random/stream.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace nomos::games::dynamics {

namespace detail {
inline void require(bool ok, const std::string& what) {
  if (!ok) throw std::invalid_argument("nomos::games::dynamics: " + what);
}
}  // namespace detail

using State = std::vector<double>;
using Game = std::function<std::vector<double>(const State&)>;
using Matrix = std::vector<std::vector<double>>;

// A linear population game F(x) = A x (random matching in a symmetric two-player game).
inline Game linear(const Matrix& A) {
  const std::size_t n = A.size();
  for (const auto& row : A) detail::require(row.size() == n, "a square payoff matrix");
  return [A, n](const State& x) {
    std::vector<double> f(n, 0.0);
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j < n; j++) f[i] += A[i][j] * x[j];
    return f;
  };
}

// Rock–Paper–Scissors with winning payoff w and losing payoff l: standard (w = l), good (w > l), bad (w < l) (Sandholm 2015, Example 13.20).
inline Matrix rock_paper_scissors(double w, double l) { return {{0, -l, w}, {w, 0, -l}, {-l, w, 0}}; }

inline void check_state(const State& x) {
  detail::require(!x.empty(), "a state has at least one strategy");
  double s = 0;
  for (double v : x) detail::require(v >= -1e-12, "a state's shares are not negative"), s += v;
  detail::require(std::fabs(s - 1) <= 1e-9, "a state's shares sum to one");
}

inline double mean_payoff(const State& x, const std::vector<double>& f) {
  double m = 0;
  for (std::size_t i = 0; i < x.size(); i++) m += x[i] * f[i];
  return m;
}

enum class Protocol { replicator, best_response, logit, bnn, smith, maynard_smith };

inline const char* name(Protocol p) {
  switch (p) {
    case Protocol::replicator: return "replicator";
    case Protocol::best_response: return "best response";
    case Protocol::logit: return "logit";
    case Protocol::bnn: return "Brown–von Neumann–Nash";
    case Protocol::smith: return "Smith";
    case Protocol::maynard_smith: return "Maynard Smith replicator";
  }
  return "";
}

// The uniform law on the strategies whose payoff is within `tie` of the best.
inline std::vector<double> best_responses(const std::vector<double>& f, double tie = 1e-12) {
  const double m = *std::max_element(f.begin(), f.end());
  std::vector<double> b(f.size(), 0.0);
  double k = 0;
  for (std::size_t i = 0; i < f.size(); i++)
    if (f[i] >= m - tie) b[i] = 1, k++;
  for (auto& v : b) v /= k;
  return b;
}

inline std::vector<double> logit_choice(const std::vector<double>& f, double eta) {
  detail::require(eta > 0, "the logit's noise level is positive");
  const double m = *std::max_element(f.begin(), f.end());
  std::vector<double> p(f.size());
  double z = 0;
  for (std::size_t i = 0; i < f.size(); i++) p[i] = std::exp((f[i] - m) / eta), z += p[i];
  for (auto& v : p) v /= z;
  return p;
}

// The protocol's matrix ρ(π, x): conditional switch rates from i (row) to j (column). The Maynard Smith replicator has no protocol here
// (it is imitation of success with repeated sampling); its mean dynamic is given directly.
inline Matrix protocol(Protocol p, const std::vector<double>& pi, const State& x, double eta = 0.1) {
  const std::size_t n = x.size();
  detail::require(pi.size() == n, "a payoff for every strategy");
  Matrix r(n, std::vector<double>(n, 0.0));
  switch (p) {
    case Protocol::replicator:
      for (std::size_t i = 0; i < n; i++)
        for (std::size_t j = 0; j < n; j++) r[i][j] = x[j] * std::max(0.0, pi[j] - pi[i]);
      break;
    case Protocol::best_response: {
      const auto b = best_responses(pi);
      for (std::size_t i = 0; i < n; i++) r[i] = b;
      break;
    }
    case Protocol::logit: {
      const auto c = logit_choice(pi, eta);
      for (std::size_t i = 0; i < n; i++) r[i] = c;
      break;
    }
    case Protocol::bnn: {
      const double m = mean_payoff(x, pi);
      for (std::size_t i = 0; i < n; i++)
        for (std::size_t j = 0; j < n; j++) r[i][j] = std::max(0.0, pi[j] - m);
      break;
    }
    case Protocol::smith:
      for (std::size_t i = 0; i < n; i++)
        for (std::size_t j = 0; j < n; j++) r[i][j] = std::max(0.0, pi[j] - pi[i]);
      break;
    case Protocol::maynard_smith:
      detail::require(false, "the Maynard Smith replicator is given by its mean dynamic, not a protocol");
  }
  return r;
}

// The mean dynamic of a protocol's matrix at x: ẋ_i = Σ_j x_j ρ_ji − x_i Σ_j ρ_ij (eq. 13.3).
inline std::vector<double> mean_dynamic(const Matrix& r, const State& x) {
  const std::size_t n = x.size();
  std::vector<double> v(n, 0.0);
  for (std::size_t i = 0; i < n; i++) {
    double in = 0, out = 0;
    for (std::size_t j = 0; j < n; j++) in += x[j] * r[j][i], out += r[i][j];
    v[i] = in - x[i] * out;
  }
  return v;
}

// The vector field V(x) of a dynamic, by its closed form (Table 13.1), for the game F.
inline std::vector<double> field(Protocol p, const Game& F, const State& x, double eta = 0.1) {
  const auto f = F(x);
  const std::size_t n = x.size();
  detail::require(f.size() == n, "the game gives a payoff for every strategy");
  const double m = mean_payoff(x, f);
  std::vector<double> v(n, 0.0);
  switch (p) {
    case Protocol::replicator:
      for (std::size_t i = 0; i < n; i++) v[i] = x[i] * (f[i] - m);
      break;
    case Protocol::maynard_smith:
      detail::require(m > 0, "the Maynard Smith replicator needs a positive mean payoff");
      for (std::size_t i = 0; i < n; i++) v[i] = x[i] * (f[i] - m) / m;
      break;
    case Protocol::best_response: {
      const auto b = best_responses(f);
      for (std::size_t i = 0; i < n; i++) v[i] = b[i] - x[i];
      break;
    }
    case Protocol::logit: {
      const auto c = logit_choice(f, eta);
      for (std::size_t i = 0; i < n; i++) v[i] = c[i] - x[i];
      break;
    }
    case Protocol::bnn: {
      double s = 0;
      for (std::size_t j = 0; j < n; j++) s += std::max(0.0, f[j] - m);
      for (std::size_t i = 0; i < n; i++) v[i] = std::max(0.0, f[i] - m) - x[i] * s;
      break;
    }
    case Protocol::smith:
      for (std::size_t i = 0; i < n; i++) {
        double in = 0, out = 0;
        for (std::size_t j = 0; j < n; j++) in += x[j] * std::max(0.0, f[i] - f[j]), out += std::max(0.0, f[j] - f[i]);
        v[i] = in - x[i] * out;
      }
      break;
  }
  return v;
}

// A trajectory of the dynamic from x0 over [0, T] by the classical Runge–Kutta method with step h; the states every `every` steps.
struct Trajectory {
  std::vector<double> t;
  std::vector<State> x;
};

inline Trajectory trajectory(Protocol p, const Game& F, State x, double T, double h, std::size_t every = 1, double eta = 0.1) {
  check_state(x);
  detail::require(T >= 0 && h > 0 && every > 0, "a horizon not negative, a positive step and a positive stride");
  const std::size_t n = x.size(), steps = static_cast<std::size_t>(std::ceil(T / h - 1e-9));
  Trajectory out;
  out.t.push_back(0), out.x.push_back(x);
  auto add = [n](const State& a, const std::vector<double>& b, double c) {
    State y(n);
    for (std::size_t i = 0; i < n; i++) y[i] = a[i] + c * b[i];
    return y;
  };
  for (std::size_t s = 1; s <= steps; s++) {
    const auto k1 = field(p, F, x, eta), k2 = field(p, F, add(x, k1, h / 2), eta), k3 = field(p, F, add(x, k2, h / 2), eta), k4 = field(p, F, add(x, k3, h), eta);
    double tot = 0;
    for (std::size_t i = 0; i < n; i++) x[i] = std::max(0.0, x[i] + h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])), tot += x[i];
    for (auto& v : x) v /= tot;
    if (s % every == 0 || s == steps) out.t.push_back(s * h), out.x.push_back(x);
  }
  return out;
}

// Whether x is a Nash equilibrium of F: every strategy in use earns the best payoff, within `tol`.
inline bool is_nash(const Game& F, const State& x, double tol = 1e-9) {
  const auto f = F(x);
  const double m = *std::max_element(f.begin(), f.end());
  for (std::size_t i = 0; i < x.size(); i++)
    if (x[i] > tol && f[i] < m - tol) return false;
  return true;
}

// The incentive properties at a point: V(x)′F(x), and whether V(x) vanishes.
struct Incentive {
  double inner = 0;    // V(x)′F(x)
  double speed = 0;    // |V(x)|
  bool rest = false;   // V(x) = 0 within tol
  bool nash = false;
};
inline Incentive incentive(Protocol p, const Game& F, const State& x, double eta = 0.1, double tol = 1e-9) {
  const auto v = field(p, F, x, eta), f = F(x);
  Incentive I;
  for (std::size_t i = 0; i < x.size(); i++) I.inner += v[i] * f[i], I.speed += v[i] * v[i];
  I.speed = std::sqrt(I.speed), I.rest = I.speed <= tol, I.nash = is_nash(F, x, tol);
  return I;
}

// The Lyapunov functions of Table 13.4, for contractive games with equilibrium x*.
inline double replicator_entropy(const State& xstar, const State& x) {
  double h = 0;
  for (std::size_t i = 0; i < x.size(); i++)
    if (xstar[i] > 0) h += xstar[i] * std::log(xstar[i] / x[i]);
  return h;
}
inline double bnn_lyapunov(const Game& F, const State& x) {
  const auto f = F(x);
  const double m = mean_payoff(x, f);
  double g = 0;
  for (double v : f) g += 0.5 * std::pow(std::max(0.0, v - m), 2);
  return g;
}
inline double smith_lyapunov(const Game& F, const State& x) {
  const auto f = F(x);
  double s = 0;
  for (std::size_t i = 0; i < x.size(); i++)
    for (std::size_t j = 0; j < x.size(); j++) s += 0.5 * x[i] * std::pow(std::max(0.0, f[j] - f[i]), 2);
  return s;
}
inline double best_response_lyapunov(const Game& F, const State& x) {
  const auto f = F(x);
  return *std::max_element(f.begin(), f.end()) - mean_payoff(x, f);
}

// The finite population: N agents, each with an exponential clock of rate one; an agent playing i switches to j with probability
// ρ_ij / R, R a bound on the rows' sums (taken as the largest row sum over the visited states, raised when exceeded: uniformisation that
// leaves the process's law unchanged while R bounds the rates). Simulated exactly by Gillespie's direct method on the counts; the shares
// at the times asked.
inline std::vector<State> finite_population(Protocol p, const Game& F, const std::vector<long long>& counts, const std::vector<double>& times,
                                            calx::random::Stream& g, double eta = 0.1) {
  detail::require(p != Protocol::maynard_smith, "the finite population needs a protocol");
  detail::require(!counts.empty(), "counts of every strategy");
  long long N = 0;
  for (auto c : counts) detail::require(c >= 0, "counts are not negative"), N += c;
  detail::require(N > 0, "a population with agents");
  for (std::size_t k = 1; k < times.size(); k++) detail::require(times[k] >= times[k - 1], "times in order");
  const std::size_t n = counts.size();
  std::vector<long long> c = counts;
  std::vector<State> out;
  double t = 0;
  std::size_t next = 0;
  auto shares = [&] {
    State x(n);
    for (std::size_t i = 0; i < n; i++) x[i] = static_cast<double>(c[i]) / static_cast<double>(N);
    return x;
  };
  while (next < times.size()) {
    const State x = shares();
    const auto r = protocol(p, F(x), x, eta);
    // Rates of each move i → j: c_i ρ_ij.
    double total = 0;
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j < n; j++)
        if (i != j) total += static_cast<double>(c[i]) * r[i][j];
    const double dt = total > 0 ? g.exponential() / total : std::numeric_limits<double>::infinity();
    while (next < times.size() && times[next] <= t + dt) out.push_back(x), next++;
    if (next >= times.size() || !(total > 0)) break;
    t += dt;
    double u = g.uniform() * total;
    for (std::size_t i = 0; i < n && u >= 0; i++)
      for (std::size_t j = 0; j < n; j++) {
        if (i == j) continue;
        const double rate = static_cast<double>(c[i]) * r[i][j];
        if (u < rate) {
          c[i]--, c[j]++;
          u = -1;
          break;
        }
        u -= rate;
      }
  }
  while (out.size() < times.size()) out.push_back(shares());
  return out;
}

}  // namespace nomos::games::dynamics
