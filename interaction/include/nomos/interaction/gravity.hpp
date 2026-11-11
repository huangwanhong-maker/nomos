// Spatial interaction: flows between origins and destinations that fall with the cost between them, and the
// accessibility that the same costs give.
//
// WILSON'S FAMILY. Origins i with outflows O_i, destinations j with inflows D_j (or attractions W_j), costs c_ij ≥ 0 and a
// deterrence f(c) = e^{−βc}, β > 0 (Wilson 1970, Entropy in Urban and Regional Modelling, Pion, ch. 2; Wilson 1971, "A
// family of spatial interaction models, and associated developments", Environ. Plan. 3, 1–32):
//   production-constrained   T_ij = O_i W_j f_ij / Σ_k W_k f_ik,
//   attraction-constrained   T_ij = D_j W_i f_ij / Σ_k W_k f_kj,
//   doubly constrained       T_ij = A_i O_i B_j D_j f_ij, with A and B such that Σ_j T_ij = O_i and Σ_i T_ij = D_j.
// The doubly constrained flows are the unique maximiser of the entropy −Σ T_ij (ln T_ij − 1) subject to the two margins
// and a total cost Σ T_ij c_ij = C, β its multiplier (Wilson 1967, "A statistical theory of spatial distribution models",
// Transp. Res. 1, 253–269); equivalently, of −Σ T ln(T/f) subject to the margins alone, the I-projection of f onto the
// tables with those margins.
//
// BALANCING. A and B are found by Furness's iteration (Furness 1965, "Time function iteration", Traffic Eng. Control 7,
// 458–460), which alternately scales the rows and the columns of f to their margins: it is Deming and Stephan's (1940)
// iterative proportional fitting, Sinkhorn's (1964) scaling, and cyclic Bregman projection on the entropy (Bregman 1967).
// It converges to the balanced table whenever one exists with the support of f (Sinkhorn and Knopp 1967; Csiszár 1975),
// and needs Σ O = Σ D. Its residual is the largest relative error of a margin.
//
// CALIBRATION. The mean cost of the doubly constrained flows, C(β)/ΣO, falls as β rises (its derivative is minus a
// variance of the costs under the flows, the dual being convex), so that the β which reproduces an observed mean cost is
// unique; Hyman's method (Hyman 1969, "The calibration of trip distribution models", Environ. Plan. 1, 105–112) finds it
// by a secant on β, safeguarded here by a bracket. Fitted to an observed table, the flows then match its two margins and
// its total cost, which are exactly the score equations of the Poisson likelihood of the table with origin and destination
// effects and the cost as regressor: Hyman's β is the maximum likelihood estimate (Flowerdew and Aitkin 1982, "A method of
// fitting the gravity model based on the Poisson distribution", J. Reg. Sci. 22, 191–202), nomos::econometrics'
// `poisson_two_way` with −β the cost's coefficient.
//
// FIT. The common part of commuters, CPC = 2 Σ min(T_obs, T_mod) / (Σ T_obs + Σ T_mod), the share of flows the model
// places where they are observed (Gargiulo et al. 2012, "Commuting network models: getting the essentials", J. Artif. Soc.
// Soc. Simul. 15(2), 6; Lenormand, Bassolas and Ramasco 2016, "Systematic comparison of trip distribution laws and models",
// J. Transp. Geogr. 51, 158–169), and the standardised root mean square error (Knudsen and Fotheringham 1986).
//
// ACCESSIBILITY. The logsum A_i = (1/β) ln Σ_j W_j e^{−β c_ij} is the expected greatest utility of a traveller from i
// choosing among the destinations by the logit with utilities ln W_j − β c_ij, divided by β (Williams 1977, "On the
// formation of travel demand models and economic evaluation measures of user benefit", Environ. Plan. A 9, 285–344;
// Ben-Akiva and Lerman 1985, Discrete Choice Analysis, MIT Press, §10.4): in units of cost, ∂A_i/∂c_ij = −P_ij and
// ∂A_i/∂ln W_j = P_ij/β, P_i the production-constrained shares. Hansen's (1959) gravity accessibility is β A_i's exponent.
//
// Costs are a dense table, origins by destinations; +∞ is a pair that cannot interact (f = 0). Units are the caller's.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace nomos::interaction {

namespace detail {
inline void require(bool ok, const std::string& what) {
  if (!ok) throw std::invalid_argument("nomos::interaction: " + what);
}
}  // namespace detail

// A dense table, origins (rows) by destinations (columns).
struct Table {
  std::size_t rows = 0, cols = 0;
  std::vector<double> v;

  Table() = default;
  Table(std::size_t r, std::size_t c, double fill = 0.0) : rows(r), cols(c), v(r * c, fill) {}
  double& operator()(std::size_t i, std::size_t j) { return v[i * cols + j]; }
  double operator()(std::size_t i, std::size_t j) const { return v[i * cols + j]; }
  double total() const {
    double s = 0;
    for (const double x : v) s += x;
    return s;
  }
  std::vector<double> row_sums() const {
    std::vector<double> s(rows, 0.0);
    for (std::size_t i = 0; i < rows; i++)
      for (std::size_t j = 0; j < cols; j++) s[i] += (*this)(i, j);
    return s;
  }
  std::vector<double> col_sums() const {
    std::vector<double> s(cols, 0.0);
    for (std::size_t i = 0; i < rows; i++)
      for (std::size_t j = 0; j < cols; j++) s[j] += (*this)(i, j);
    return s;
  }
};

// f_ij = e^{−β c_ij}, zero where the cost is infinite.
inline Table deterrence(const Table& cost, double beta) {
  detail::require(beta >= 0, "a deterrence's rate is not negative");
  Table f(cost.rows, cost.cols);
  for (std::size_t k = 0; k < f.v.size(); k++) {
    detail::require(cost.v[k] >= 0, "a cost is not negative");
    f.v[k] = std::isfinite(cost.v[k]) ? std::exp(-beta * cost.v[k]) : 0.0;
  }
  return f;
}

// T_ij = O_i W_j f_ij / Σ_k W_k f_ik: each origin's outflow shared among the destinations in proportion to W_j f_ij.
inline Table production_constrained(std::span<const double> outflow, std::span<const double> attraction, const Table& cost, double beta) {
  detail::require(outflow.size() == cost.rows && attraction.size() == cost.cols, "an outflow per origin and an attraction per destination");
  const Table f = deterrence(cost, beta);
  Table t(cost.rows, cost.cols);
  for (std::size_t i = 0; i < cost.rows; i++) {
    double s = 0;
    for (std::size_t j = 0; j < cost.cols; j++) s += attraction[j] * f(i, j);
    if (!(outflow[i] > 0)) continue;
    detail::require(s > 0, "an origin with an outflow reaches no attraction");
    for (std::size_t j = 0; j < cost.cols; j++) t(i, j) = outflow[i] * attraction[j] * f(i, j) / s;
  }
  return t;
}

// T_ij = D_j W_i f_ij / Σ_k W_k f_kj: each destination's inflow drawn from the origins in proportion to W_i f_ij.
inline Table attraction_constrained(std::span<const double> inflow, std::span<const double> weight, const Table& cost, double beta) {
  detail::require(inflow.size() == cost.cols && weight.size() == cost.rows, "an inflow per destination and a weight per origin");
  const Table f = deterrence(cost, beta);
  Table t(cost.rows, cost.cols);
  for (std::size_t j = 0; j < cost.cols; j++) {
    double s = 0;
    for (std::size_t i = 0; i < cost.rows; i++) s += weight[i] * f(i, j);
    if (!(inflow[j] > 0)) continue;
    detail::require(s > 0, "a destination with an inflow is reached by no origin");
    for (std::size_t i = 0; i < cost.rows; i++) t(i, j) = inflow[j] * weight[i] * f(i, j) / s;
  }
  return t;
}

struct Balanced {
  Table flows;
  std::vector<double> a, b;   // the balancing factors: T_ij = a_i b_j seed_ij (a_i = A_i O_i, b_j = B_j D_j)
  int iterations = 0;
  double residual = 0;        // the largest relative error of a margin
  bool converged = false;
};

// Furness's balancing of a seed table to the margins (Σ outflow = Σ inflow), from factors `start` if given.
inline Balanced balance(const Table& seed, std::span<const double> outflow, std::span<const double> inflow, double tolerance = 1e-12,
                        int max_iterations = 100000, std::vector<double> start_b = {}) {
  const std::size_t R = seed.rows, C = seed.cols;
  detail::require(outflow.size() == R && inflow.size() == C, "an outflow per origin and an inflow per destination");
  double so = 0, si = 0;
  for (const double x : outflow) detail::require(x >= 0, "an outflow is not negative"), so += x;
  for (const double x : inflow) detail::require(x >= 0, "an inflow is not negative"), si += x;
  detail::require(std::fabs(so - si) <= 1e-9 * std::max(so, si), "the outflows and the inflows have one total");
  Balanced out;
  out.a.assign(R, 0.0);
  out.b = start_b.size() == C ? std::move(start_b) : std::vector<double>(C, 1.0);
  for (out.iterations = 0; out.iterations < max_iterations; out.iterations++) {
    for (std::size_t i = 0; i < R; i++) {
      double s = 0;
      for (std::size_t j = 0; j < C; j++) s += seed(i, j) * out.b[j];
      out.a[i] = s > 0 ? outflow[i] / s : 0.0;
    }
    std::vector<double> col(C, 0.0);
    for (std::size_t i = 0; i < R; i++)
      for (std::size_t j = 0; j < C; j++) col[j] += out.a[i] * seed(i, j);
    double worst = 0;
    for (std::size_t j = 0; j < C; j++) {
      // After the rows are exact, the columns' error is the margin's.
      const double have = col[j] * out.b[j];
      if (inflow[j] > 0) worst = std::max(worst, std::fabs(have - inflow[j]) / inflow[j]);
      else if (have > 0) worst = std::max(worst, 1.0);
      out.b[j] = col[j] > 0 ? inflow[j] / col[j] : 0.0;
    }
    out.residual = worst;
    if (worst <= tolerance) { out.converged = true; break; }
  }
  out.flows = Table(R, C);
  for (std::size_t i = 0; i < R; i++)
    for (std::size_t j = 0; j < C; j++) out.flows(i, j) = out.a[i] * seed(i, j) * out.b[j];
  return out;
}

// The doubly constrained flows. The seed is e^{−β (c_ij − min_k c_ik)}: scaling a row of the seed changes no balanced
// table, and the shift keeps every row from underflowing at a large β; the factors returned are the seed's.
inline Balanced doubly_constrained(std::span<const double> outflow, std::span<const double> inflow, const Table& cost, double beta,
                                   double tolerance = 1e-12, int max_iterations = 100000) {
  detail::require(outflow.size() == cost.rows && inflow.size() == cost.cols, "an outflow per origin and an inflow per destination");
  Table shifted = cost;
  for (std::size_t i = 0; i < cost.rows; i++) {
    double least = std::numeric_limits<double>::infinity();
    for (std::size_t j = 0; j < cost.cols; j++) least = std::min(least, cost(i, j));
    if (std::isfinite(least))
      for (std::size_t j = 0; j < cost.cols; j++) shifted(i, j) = cost(i, j) - least;
  }
  return balance(deterrence(shifted, beta), outflow, inflow, tolerance, max_iterations);
}

// Σ T c / Σ T over the pairs that interact.
inline double mean_cost(const Table& flows, const Table& cost) {
  double s = 0, n = 0;
  for (std::size_t k = 0; k < flows.v.size(); k++)
    if (flows.v[k] > 0) s += flows.v[k] * cost.v[k], n += flows.v[k];
  return n > 0 ? s / n : 0.0;
}

struct Calibrated {
  double beta = 0;
  Balanced fit;
  double mean_cost = 0;   // the model's, at β
  int iterations = 0;
  bool converged = false;
};

// Hyman's calibration: the β > 0 at which the doubly constrained flows' mean cost is `target`.
inline Calibrated calibrate(std::span<const double> outflow, std::span<const double> inflow, const Table& cost, double target,
                            double tolerance = 1e-10, int max_iterations = 200) {
  detail::require(target > 0, "a mean cost to reproduce");
  const auto at = [&](double beta, Calibrated& c) {
    c.fit = doubly_constrained(outflow, inflow, cost, beta, 1e-13);
    c.beta = beta;
    c.mean_cost = mean_cost(c.fit.flows, cost);
    return c.mean_cost - target;
  };
  Calibrated out;
  // β₀ = 1/target (Hyman), then a bracket [lo, hi] on which the mean cost falls through the target.
  double b0 = 1 / target, g0 = at(b0, out);
  double lo = b0, hi = b0, glo = g0, ghi = g0;
  for (int k = 0; k < 200 && glo < 0; k++) glo = at(lo /= 2, out);   // mean too low: β too high
  for (int k = 0; k < 200 && ghi > 0; k++) ghi = at(hi *= 2, out);   // mean too high: β too low
  detail::require(glo >= 0 && ghi <= 0, "no β reproduces the mean cost (it lies outside what the margins allow)");
  // A secant on ln β, kept within the bracket (bisection when it would leave it).
  double x0 = std::log(lo), x1 = std::log(hi), f0 = glo, f1 = ghi;
  for (out.iterations = 0; out.iterations < max_iterations; out.iterations++) {
    double x = f1 != f0 ? x1 - f1 * (x1 - x0) / (f1 - f0) : 0.5 * (x0 + x1);
    if (!(x > std::min(std::log(lo), std::log(hi)) && x < std::max(std::log(lo), std::log(hi)))) x = 0.5 * (std::log(lo) + std::log(hi));
    const double g = at(std::exp(x), out);
    if (std::fabs(g) <= tolerance * target) { out.converged = true; break; }
    (g > 0 ? lo : hi) = std::exp(x);
    x0 = x1, f0 = f1, x1 = x, f1 = g;
  }
  return out;
}

// Calibrated to an observed table: its margins and its mean cost (the Poisson maximum likelihood estimate).
inline Calibrated calibrate(const Table& observed, const Table& cost, double tolerance = 1e-10, int max_iterations = 200) {
  detail::require(observed.rows == cost.rows && observed.cols == cost.cols, "the observed table and the costs have one shape");
  const auto o = observed.row_sums(), d = observed.col_sums();
  return calibrate(o, d, cost, mean_cost(observed, cost), tolerance, max_iterations);
}

// 2 Σ min(T_obs, T_mod) / (Σ T_obs + Σ T_mod).
inline double common_part(const Table& observed, const Table& modelled) {
  double common = 0, a = 0, b = 0;
  for (std::size_t k = 0; k < observed.v.size(); k++) common += std::min(observed.v[k], modelled.v[k]), a += observed.v[k], b += modelled.v[k];
  return a + b > 0 ? 2 * common / (a + b) : 1.0;
}

// The standardised root mean square error: √(Σ (T_obs − T_mod)² / n) / (Σ T_obs / n).
inline double srmse(const Table& observed, const Table& modelled) {
  const double n = static_cast<double>(observed.v.size());
  double se = 0, s = 0;
  for (std::size_t k = 0; k < observed.v.size(); k++) se += (observed.v[k] - modelled.v[k]) * (observed.v[k] - modelled.v[k]), s += observed.v[k];
  return s > 0 ? std::sqrt(se / n) / (s / n) : 0.0;
}

// A_i = (1/β) ln Σ_j W_j e^{−β c_ij}, by the log-sum-exp; −∞ where nothing is reached.
inline std::vector<double> logsum(std::span<const double> attraction, const Table& cost, double beta) {
  detail::require(beta > 0, "a logsum needs a positive β");
  detail::require(attraction.size() == cost.cols, "an attraction per destination");
  std::vector<double> out(cost.rows, -std::numeric_limits<double>::infinity());
  for (std::size_t i = 0; i < cost.rows; i++) {
    double top = -std::numeric_limits<double>::infinity();
    for (std::size_t j = 0; j < cost.cols; j++)
      if (attraction[j] > 0 && std::isfinite(cost(i, j))) top = std::max(top, std::log(attraction[j]) - beta * cost(i, j));
    if (!std::isfinite(top)) continue;
    double s = 0;
    for (std::size_t j = 0; j < cost.cols; j++)
      if (attraction[j] > 0 && std::isfinite(cost(i, j))) s += std::exp(std::log(attraction[j]) - beta * cost(i, j) - top);
    out[i] = (top + std::log(s)) / beta;
  }
  return out;
}

}  // namespace nomos::interaction
