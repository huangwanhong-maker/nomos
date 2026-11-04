// Poisson pseudo-maximum likelihood with two sets of absorbed fixed effects.
//
// THE MODEL. Observations are cells (r, c) of a two-way table (origins and destinations, rows and
// columns), possibly incomplete, with an outcome y ≥ 0, regressors x, and
//   E(y | x) = exp(o + xβ + α_r + γ_c),
// o an optional offset. The Poisson quasi-likelihood estimator (PPML) is consistent whenever this
// conditional mean is right, whatever the distribution of y, which need not be integer, and keeps
// the cells with y = 0 (Gourieroux, Monfort and Trognon 1984; Santos Silva and Tenreyro 2006, the
// estimator of gravity equations). With the two sets of effects the Poisson estimator is also free
// of the incidental parameter bias in β (Fernández-Val and Weidner 2016 for the bias of others;
// Weidner and Zylkin 2021 for the Poisson case's consistency and its small bias in the standard
// errors).
//
// THE ALGORITHM. Iteratively reweighted least squares (McCullagh and Nelder 1989, ch. 2): the
// working response z = η − o + (y − μ)/μ and weights w = μ; each step is the weighted least squares
// of z on [X, row indicators, column indicators], solved by absorbing the indicators (absorb.hpp):
// z and X are residualised by alternating weighted projections, β is the weighted OLS of the
// residualised z on the residualised X (by tyche), and the fitted working response is z minus the
// common residual (Frisch–Waugh–Lovell), so the dense dummies (hundreds of columns for a table of
// a few hundred rows and columns) are never formed (Guimarães and Portugal 2010; Gaure 2013;
// Correia, Guimarães and Zylkin 2020, ppmlhdfe, whose algorithm this is, without its acceleration).
// It starts from μ = (y + ȳ)/2 and stops when the deviance D = 2 Σ [y ln(y/μ) − (y − μ)] changes
// by less than tolerance · (|D| + 0.1), halving the step if the deviance is not finite.
//
// SEPARATION. A row (or column) whose outcomes are all zero has its effect at −∞ in the limit: its
// cells are dropped before fitting (their fitted value is 0, their effect −∞), as ppmlhdfe drops
// them; dropping them changes nothing else, since they contribute nothing to the other score
// equations (Santos Silva and Tenreyro 2010; Correia, Guimarães and Zylkin 2019). Separation by
// the regressors is not detected here.
//
// THE COVARIANCE of β̂ is the sandwich B (Σ_g s_g s_gᵀ) B with B = (X̃ᵀWX̃)⁻¹, X̃ the W-residualised
// regressors at convergence and s_g = Σ_{i∈g} x̃_i (y_i − μ_i), each cell its own group (HC0) or
// the caller's clusters (by pair, by origin, …) with the factor G/(G − 1). By Frisch–Waugh–Lovell
// this is exactly the β block of the sandwich of the full dummy-variable Poisson regression
// (checked against tyche's GLM). Inference is normal.
//
// THE EFFECTS α̂_r, γ̂_c are recovered from η − o − Xβ by backfitting (absorb.hpp), with the first
// column of each connected component at zero.

#pragma once

#include <tyche/probability/distributions.hpp>
#include <tyche/probability/p_values.hpp>
#include <tyche/regression/linear.hpp>

#include <cmath>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

#include "absorb.hpp"
#include "instrumental.hpp"

namespace nomos::econometrics {

enum class PoissonCovariance { hc0, cluster };

struct PoissonTwoWayOptions {
  PoissonCovariance covariance = PoissonCovariance::hc0;
  std::vector<std::int64_t> clusters;  // one label per observation, for PoissonCovariance::cluster
  std::vector<double> offset;          // empty, or one oᵢ per observation
  double tolerance = 1e-10;            // on |ΔD| / (|D| + 0.1)
  int max_iterations = 200;
  double confidence = 0.95;
  AbsorbOptions absorb;
};

struct PoissonTwoWayFit {
  std::vector<double> beta, se, z, p, ci_lower, ci_upper;  // one per column of X; NaN where absorbed by the effects
  calx::linalg::Matrix covariance;
  std::vector<bool> aliased;               // columns of X absorbed by the effects (constant within rows or columns)
  std::vector<double> row_effects;         // α̂_r for r = 0 … R − 1; −∞ for dropped rows, NaN for rows without cells
  std::vector<double> column_effects;      // γ̂_c likewise
  std::vector<double> fitted;              // μ̂ for every observation (0 for dropped cells)
  std::vector<bool> dropped;               // cells dropped because their row or column has no positive outcome
  std::size_t used = 0, dropped_rows = 0, dropped_columns = 0;
  std::size_t clusters = 0;
  double deviance = 0;                     // over the cells used
  int iterations = 0;
  bool converged = false;
  int sweeps = 0;                          // the most sweeps of the projections in any iteration
};

namespace detail {
inline double poisson_deviance(std::span<const double> y, std::span<const double> mu) {
  double d = 0;
  for (std::size_t i = 0; i < y.size(); i++) d += 2 * ((y[i] > 0 ? y[i] * std::log(y[i] / mu[i]) : 0.0) - (y[i] - mu[i]));
  return d;
}
}  // namespace detail

// PPML of y on X with row and column fixed effects: row[i] and column[i] are the cell's levels.
inline PoissonTwoWayFit poisson_two_way(const Rectangular& X, std::span<const double> y, std::span<const std::size_t> row, std::span<const std::size_t> column,
                                        const PoissonTwoWayOptions& options = {}) {
  const std::size_t n = y.size(), k = X.cols;
  detail::require(row.size() == n && column.size() == n && (k == 0 || X.rows == n), "a row, a column and regressors for every observation");
  detail::require(options.offset.empty() || options.offset.size() == n, "one offset per observation");
  detail::require(options.covariance != PoissonCovariance::cluster || options.clusters.size() == n, "a cluster label per observation");
  std::size_t R = 0, C = 0;
  for (std::size_t i = 0; i < n; i++) {
    detail::require(y[i] >= 0 && std::isfinite(y[i]), "outcomes y ≥ 0");
    R = std::max(R, row[i] + 1);
    C = std::max(C, column[i] + 1);
  }
  // Separation: rows and columns without a positive outcome.
  std::vector<double> rt(R, 0.0), ct(C, 0.0);
  std::vector<bool> rseen(R, false), cseen(C, false);
  for (std::size_t i = 0; i < n; i++) {
    rt[row[i]] += y[i];
    ct[column[i]] += y[i];
    rseen[row[i]] = cseen[column[i]] = true;
  }
  PoissonTwoWayFit fit;
  for (std::size_t r = 0; r < R; r++) fit.dropped_rows += rseen[r] && rt[r] == 0;
  for (std::size_t c = 0; c < C; c++) fit.dropped_columns += cseen[c] && ct[c] == 0;
  fit.dropped.assign(n, false);
  std::vector<std::size_t> kept;
  for (std::size_t i = 0; i < n; i++) {
    fit.dropped[i] = rt[row[i]] == 0 || ct[column[i]] == 0;
    if (!fit.dropped[i]) kept.push_back(i);
  }
  const std::size_t m = kept.size();
  fit.used = m;
  detail::require(m > k, "more cells than regressors after dropping the separated ones");
  Factor rf, cf;
  rf.levels = R;
  cf.levels = C;
  std::vector<double> yk(m), ok(m, 0.0);
  Rectangular Xk(m, k);
  for (std::size_t t = 0; t < m; t++) {
    const std::size_t i = kept[t];
    rf.level.push_back(row[i]);
    cf.level.push_back(column[i]);
    yk[t] = y[i];
    if (!options.offset.empty()) ok[t] = options.offset[i];
    for (std::size_t j = 0; j < k; j++) Xk(t, j) = X(i, j);
  }
  const std::vector<Factor> factors{rf, cf};
  std::vector<double> xnorm(k, 0.0);
  for (std::size_t j = 0; j < k; j++)
    for (std::size_t t = 0; t < m; t++) xnorm[j] = std::max(xnorm[j], std::fabs(Xk(t, j)));

  // Start.
  double ybar = 0;
  for (double v : yk) ybar += v;
  ybar /= static_cast<double>(m);
  std::vector<double> mu(m), eta(m), z(m);
  for (std::size_t t = 0; t < m; t++) {
    mu[t] = 0.5 * (yk[t] + ybar);
    eta[t] = std::log(mu[t]);
  }
  double dev_old = detail::poisson_deviance(yk, mu);
  std::vector<double> beta(k, 0.0);
  fit.aliased.assign(k, false);
  // One weighted least-squares step: residualise z and X under weights mu, regress; returns the residuals.
  auto step = [&](std::vector<double>& b, Rectangular& Xt, tr::LinearFit* out) {
    std::vector<std::vector<double>> cols(k + 1, std::vector<double>(m));
    for (std::size_t t = 0; t < m; t++) {
      z[t] = eta[t] - ok[t] + (yk[t] - mu[t]) / mu[t];
      cols[0][t] = z[t];
      for (std::size_t j = 0; j < k; j++) cols[j + 1][t] = Xk(t, j);
    }
    fit.sweeps = std::max(fit.sweeps, demean(cols, factors, mu, options.absorb));
    Xt = Rectangular(m, k);
    for (std::size_t j = 0; j < k; j++) {
      double left = 0;
      for (double v : cols[j + 1]) left = std::max(left, std::fabs(v));
      fit.aliased[j] = left <= 1e-9 * xnorm[j];
      for (std::size_t t = 0; t < m; t++) Xt(t, j) = fit.aliased[j] ? 0.0 : cols[j + 1][t];
    }
    std::vector<double> e = cols[0];
    if (k > 0) {
      tr::OlsOptions o;
      o.weights = mu;
      o.drop_collinear = true;
      auto f = tr::ols(Xt, cols[0], o);
      b = f.beta;
      e = f.residuals;
      if (out) *out = std::move(f);
    }
    return e;
  };
  Rectangular Xt;
  for (int it = 1; it <= options.max_iterations; it++) {
    const auto e = step(beta, Xt, nullptr);
    std::vector<double> eta_new(m), mu_new(m);
    double dev = 0;
    for (int half = 0; half < 40; half++) {
      bool ok_step = true;
      for (std::size_t t = 0; t < m && ok_step; t++) {
        const double target = ok[t] + z[t] - e[t];
        eta_new[t] = half == 0 ? target : eta[t] + (target - eta[t]) * std::ldexp(1.0, -half);
        mu_new[t] = std::exp(eta_new[t]);
        ok_step = std::isfinite(mu_new[t]) && mu_new[t] > 0;
      }
      if (ok_step) {
        dev = detail::poisson_deviance(yk, mu_new);
        if (std::isfinite(dev)) break;
      }
    }
    eta = eta_new;
    mu = mu_new;
    fit.iterations = it;
    if (std::fabs(dev - dev_old) / (std::fabs(dev) + 0.1) < options.tolerance) {
      fit.converged = true;
      dev_old = dev;
      break;
    }
    dev_old = dev;
  }
  fit.deviance = dev_old;

  // The covariance at convergence, from the residualised regressors under the final weights.
  const double nan = std::numeric_limits<double>::quiet_NaN();
  fit.beta.assign(k, nan);
  fit.covariance = calx::linalg::Matrix(k);
  if (k > 0) {
    tr::LinearFit last;
    std::vector<double> b2;
    (void)step(b2, Xt, &last);
    for (std::size_t j = 0; j < k; j++)
      if (!fit.aliased[j] && !last.aliased[j]) fit.beta[j] = beta[j];
    std::vector<double> u(m);
    for (std::size_t t = 0; t < m; t++) u[t] = yk[t] - mu[t];
    std::vector<std::int64_t> cl;
    if (options.covariance == PoissonCovariance::cluster)
      for (std::size_t i : kept) cl.push_back(options.clusters[i]);
    std::size_t groups = 0;
    const auto M = detail::sandwich_meat(Xt, u, cl, groups);
    fit.covariance = detail::sandwich(last.xtx_inverse, M);
    if (options.covariance == PoissonCovariance::cluster) {
      detail::require(groups >= 2, "at least two clusters");
      const double G = static_cast<double>(groups);
      for (double& v : fit.covariance.a) v *= G / (G - 1);
      fit.clusters = groups;
    }
  }
  const double crit = tyche::probability::Normal{}.upper_quantile((1 - options.confidence) / 2);
  fit.se = fit.z = fit.p = fit.ci_lower = fit.ci_upper = std::vector<double>(k, nan);
  for (std::size_t j = 0; j < k; j++) {
    if (std::isnan(fit.beta[j])) continue;
    fit.se[j] = std::sqrt(fit.covariance(j, j));
    fit.z[j] = fit.beta[j] / fit.se[j];
    fit.p[j] = tyche::probability::two_sided_normal(fit.z[j]);
    fit.ci_lower[j] = fit.beta[j] - crit * fit.se[j];
    fit.ci_upper[j] = fit.beta[j] + crit * fit.se[j];
  }
  // Fitted values and the effects.
  fit.fitted.assign(n, 0.0);
  std::vector<double> d(m);
  for (std::size_t t = 0; t < m; t++) {
    fit.fitted[kept[t]] = mu[t];
    double xb = 0;
    for (std::size_t j = 0; j < k; j++)
      if (!std::isnan(fit.beta[j])) xb += Xk(t, j) * fit.beta[j];
    d[t] = eta[t] - ok[t] - xb;
  }
  const auto effects = recover_effects(d, factors, {}, options.absorb);
  fit.row_effects.assign(R, nan);
  fit.column_effects.assign(C, nan);
  for (std::size_t t = 0; t < m; t++) {
    fit.row_effects[rf.level[t]] = effects[0][rf.level[t]];
    fit.column_effects[cf.level[t]] = effects[1][cf.level[t]];
  }
  for (std::size_t r = 0; r < R; r++)
    if (rseen[r] && rt[r] == 0) fit.row_effects[r] = -std::numeric_limits<double>::infinity();
  for (std::size_t c = 0; c < C; c++)
    if (cseen[c] && ct[c] == 0) fit.column_effects[c] = -std::numeric_limits<double>::infinity();
  return fit;
}

}  // namespace nomos::econometrics
