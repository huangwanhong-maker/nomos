// Lowry's model of a metropolis: basic employment, the residents its workers bring, the services those residents need and
// the workers those services employ in turn, located by spatial interaction (Lowry 1964, A Model of Metropolis, RM-4035-RC,
// RAND Corporation).
//
// THE MODEL. Zones 1..n. Basic employment e^B (row vector) is given: what the region sells beyond itself. A worker at
// workplace j lives at i with probability R_ji, and a resident of i buys services at k with probability S_ik (both rows
// stochastic, typically production-constrained gravity shares). With a residents per worker and b service jobs per
// resident,
//   p = a e R,    e^S = b p S,    e = e^B + e^S.
// Lowry solved it by iteration, from e = e^B; Garin (1966, "A matrix formulation of the Lowry model for intrametropolitan
// activity allocation", J. Am. Inst. Plann. 32, 361–364) wrote its solution in closed form,
//   e = e^B (I − ab R S)^{−1} = e^B Σ_{k≥0} (ab R S)^k,
// the series converging whenever ab < 1, since R S is row-stochastic and its spectral radius is 1. Totals follow without the
// geography: Σ e = Σ e^B / (1 − ab) (the economic base multiplier) and Σ p = a Σ e.
//
// Caveat. Lowry's own model also capped residential density and dropped service centres below a least size, which make it
// nonlinear; those constraints are not held here.

#pragma once

#include <nomos/interaction/gravity.hpp>

#include <calx/linalg/dense.hpp>

#include <cmath>
#include <cstddef>
#include <span>
#include <vector>

namespace nomos::interaction {

struct Metropolis {
  std::vector<double> employment, service, population;   // e, e^S, p by zone
  int iterations = 0;
  double residual = 0;   // the last change of e, relative to Σ e^B
  bool converged = false;
};

namespace detail {
inline void stochastic(const Table& m, std::size_t n, const char* what) {
  require(m.rows == n && m.cols == n, std::string(what) + " is square, one row and column per zone");
  for (std::size_t i = 0; i < n; i++) {
    double s = 0;
    for (std::size_t j = 0; j < n; j++) require(m(i, j) >= 0, std::string(what) + " has no negative entry"), s += m(i, j);
    require(std::fabs(s - 1) <= 1e-9, std::string(what) + "'s rows each sum to one");
  }
}
inline std::vector<double> times(std::span<const double> x, const Table& m) {
  std::vector<double> y(m.cols, 0.0);
  for (std::size_t i = 0; i < m.rows; i++)
    if (x[i] != 0)
      for (std::size_t j = 0; j < m.cols; j++) y[j] += x[i] * m(i, j);
  return y;
}
}  // namespace detail

// Lowry's iteration from e = e^B.
inline Metropolis lowry(std::span<const double> basic, const Table& residence, const Table& service, double residents_per_worker,
                        double service_per_resident, double tolerance = 1e-13, int max_iterations = 100000) {
  const std::size_t n = basic.size();
  detail::stochastic(residence, n, "the residence matrix");
  detail::stochastic(service, n, "the service matrix");
  const double a = residents_per_worker, b = service_per_resident;
  detail::require(a >= 0 && b >= 0 && a * b < 1, "residents per worker times service jobs per resident is below one");
  double scale = 0;
  for (const double x : basic) detail::require(x >= 0, "basic employment is not negative"), scale += x;
  Metropolis out;
  out.employment.assign(basic.begin(), basic.end());
  for (out.iterations = 0; out.iterations < max_iterations; out.iterations++) {
    out.population = detail::times(out.employment, residence);
    for (double& x : out.population) x *= a;
    out.service = detail::times(out.population, service);
    for (double& x : out.service) x *= b;
    double change = 0;
    for (std::size_t j = 0; j < n; j++) {
      const double e = basic[j] + out.service[j];
      change = std::max(change, std::fabs(e - out.employment[j]));
      out.employment[j] = e;
    }
    out.residual = scale > 0 ? change / scale : 0.0;
    if (out.residual <= tolerance) { out.converged = true; break; }
  }
  out.population = detail::times(out.employment, residence);
  for (double& x : out.population) x *= a;
  return out;
}

// Garin's closed form: e solves (I − ab R S)ᵀ eᵀ = e^Bᵀ.
inline Metropolis garin(std::span<const double> basic, const Table& residence, const Table& service, double residents_per_worker,
                        double service_per_resident) {
  const std::size_t n = basic.size();
  detail::stochastic(residence, n, "the residence matrix");
  detail::stochastic(service, n, "the service matrix");
  const double ab = residents_per_worker * service_per_resident;
  detail::require(residents_per_worker >= 0 && service_per_resident >= 0 && ab < 1, "residents per worker times service jobs per resident is below one");
  calx::linalg::Matrix M(n), B(n);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++) {
      double rs = 0;
      for (std::size_t k = 0; k < n; k++) rs += residence(i, k) * service(k, j);
      M(j, i) = (i == j ? 1.0 : 0.0) - ab * rs;   // the transpose
    }
  for (std::size_t i = 0; i < n; i++) B(i, 0) = basic[i];
  const auto X = calx::linalg::solve(M, B);
  Metropolis out;
  out.employment.resize(n);
  for (std::size_t i = 0; i < n; i++) out.employment[i] = X(i, 0);
  out.population = detail::times(out.employment, residence);
  for (double& x : out.population) x *= residents_per_worker;
  out.service = detail::times(out.population, service);
  for (double& x : out.service) x *= service_per_resident;
  out.converged = true;
  return out;
}

}  // namespace nomos::interaction
