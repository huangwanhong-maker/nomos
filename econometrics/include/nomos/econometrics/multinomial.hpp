// Multinomial logit: unordered choice among J alternatives, and its average marginal effects.
//
// THE MODEL (McFadden 1974; Greene 2012, ch. 18; Cameron and Trivedi 2005, ch. 15). Observation i
// chooses category y_i ∈ {0, …, J − 1} with probability
//   p_ij = exp(x_i β_j) / Σ_m exp(x_i β_m),   β_0 = 0 (category 0 the base),
// so that β_j measures the log odds of j against the base. The log-likelihood
//   ℓ(β) = Σ_i [x_i β_{y_i} − ln Σ_m exp(x_i β_m)]
// is concave, with gradient ∂ℓ/∂β_j = Σ_i x_i (1(y_i = j) − p_ij) and Hessian blocks
// ∂²ℓ/∂β_j∂β_lᵀ = −Σ_i p_ij (1(j = l) − p_il) x_i x_iᵀ. It is maximised by calx's BFGS on −ℓ/n
// with this gradient, from β = 0. At the maximum, with an intercept, the mean predicted
// probability of each category equals its observed share (the first-order condition of the
// intercepts; checked).
//
// THE COVARIANCE of β̂ (the J − 1 blocks of k, in order): the inverse of the information −∇²ℓ
// (model), or the sandwich H⁻¹ (Σ_g s_g s_gᵀ) H⁻¹ with the score contributions s_i, HC0 or
// clustered with the factor G/(G − 1). Inference is normal.
//
// AVERAGE MARGINAL EFFECTS. For a continuous regressor l, ∂p_ij/∂x_il = p_ij (β_jl − Σ_m p_im β_ml),
// averaged over the observations; the effects on the J categories sum to zero. Their standard
// errors are by the delta method, with the Jacobian of the effects in β taken by central
// differences (calx's gradient_by_differences; as Stata's margins does by default).
// Constant columns have no effect (NaN).
//
// With J = 2 the model is the binary logit, and agrees with tyche's GLM (checked).

#pragma once

#include <tyche/probability/distributions.hpp>
#include <tyche/probability/p_values.hpp>

#include <calx/linalg/cholesky.hpp>
#include <calx/linalg/rectangular.hpp>
#include <calx/optimisation/quasi_newton.hpp>

#include <cmath>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

#include "absorb.hpp"
#include "instrumental.hpp"

namespace nomos::econometrics {

enum class MultinomialCovariance { model, hc0, cluster };

struct MultinomialOptions {
  MultinomialCovariance covariance = MultinomialCovariance::model;
  std::vector<std::int64_t> clusters;  // one label per observation, for MultinomialCovariance::cluster
  double gradient_tolerance = 1e-10;   // on the gradient of the mean log-likelihood
  double confidence = 0.95;
};

struct MultinomialFit {
  std::size_t n = 0, k = 0, categories = 0;
  std::vector<double> beta;            // (J − 1) blocks of k: beta[(j − 1) k + l] for categories j = 1 … J − 1
  std::vector<double> se, z, p, ci_lower, ci_upper;
  calx::linalg::Matrix covariance;     // of beta, (J − 1) k square
  calx::linalg::Rectangular probabilities;  // n × J, the fitted p_ij
  double log_likelihood = 0, null_log_likelihood = 0, pseudo_r_squared = 0, aic = 0, bic = 0;
  int iterations = 0;
  bool converged = false;
  std::size_t clusters = 0;
};

namespace detail {

// p_ij for every observation, from the (J − 1) k coefficients.
inline calx::linalg::Rectangular choice_probabilities(const Rectangular& X, std::span<const double> b, std::size_t J) {
  const std::size_t n = X.rows, k = X.cols;
  calx::linalg::Rectangular P(n, J);
  std::vector<double> u(J);
  for (std::size_t i = 0; i < n; i++) {
    u[0] = 0;
    double top = 0;
    for (std::size_t j = 1; j < J; j++) {
      double s = 0;
      for (std::size_t l = 0; l < k; l++) s += X(i, l) * b[(j - 1) * k + l];
      u[j] = s;
      top = std::max(top, s);
    }
    double z = 0;
    for (std::size_t j = 0; j < J; j++) z += std::exp(u[j] - top);
    for (std::size_t j = 0; j < J; j++) P(i, j) = std::exp(u[j] - top) / z;
  }
  return P;
}

}  // namespace detail

inline MultinomialFit multinomial_logit(const Rectangular& X, std::span<const std::size_t> y, std::size_t categories, const MultinomialOptions& options = {}) {
  const std::size_t n = X.rows, k = X.cols, J = categories;
  detail::require(y.size() == n && J >= 2 && k > 0, "a choice for every row of X, two categories or more");
  detail::require(options.covariance != MultinomialCovariance::cluster || options.clusters.size() == n, "a cluster label per observation");
  std::vector<double> count(J, 0.0);
  for (std::size_t c : y) {
    detail::require(c < J, "choices among the categories");
    count[c] += 1;
  }
  const std::size_t q = (J - 1) * k;
  const double dn = static_cast<double>(n);
  auto negative_mean_ll = [&](std::span<const double> b) {
    const auto P = detail::choice_probabilities(X, b, J);
    double s = 0;
    for (std::size_t i = 0; i < n; i++) s -= std::log(std::max(P(i, y[i]), std::numeric_limits<double>::min()));
    return s / dn;
  };
  auto gradient = [&](std::span<const double> b) {
    const auto P = detail::choice_probabilities(X, b, J);
    std::vector<double> g(q, 0.0);
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 1; j < J; j++) {
        const double r = (y[i] == j ? 1.0 : 0.0) - P(i, j);
        for (std::size_t l = 0; l < k; l++) g[(j - 1) * k + l] -= X(i, l) * r / dn;
      }
    return g;
  };
  calx::optimisation::BfgsOptions bo;
  bo.gradient_tolerance = options.gradient_tolerance;
  const auto res = calx::optimisation::bfgs(negative_mean_ll, gradient, std::vector<double>(q, 0.0), bo);
  MultinomialFit fit;
  fit.n = n;
  fit.k = k;
  fit.categories = J;
  fit.beta = res.x;
  fit.iterations = res.iterations;
  fit.converged = res.converged;
  fit.probabilities = detail::choice_probabilities(X, fit.beta, J);
  const auto& P = fit.probabilities;
  // The information −∇²ℓ, and the scores.
  calx::linalg::Matrix I(q);
  Rectangular S(n, q);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 1; j < J; j++) {
      for (std::size_t l = 0; l < k; l++) S(i, (j - 1) * k + l) = X(i, l) * ((y[i] == j ? 1.0 : 0.0) - P(i, j));
      for (std::size_t m = 1; m < J; m++) {
        const double w = P(i, j) * ((j == m ? 1.0 : 0.0) - P(i, m));
        for (std::size_t a = 0; a < k; a++)
          for (std::size_t c = 0; c < k; c++) I((j - 1) * k + a, (m - 1) * k + c) += w * X(i, a) * X(i, c);
      }
    }
  const auto Hinv = calx::linalg::solve(I, calx::linalg::Matrix::identity(q));
  if (options.covariance == MultinomialCovariance::model) {
    fit.covariance = Hinv;
  } else {
    std::size_t groups = 0;
    const std::vector<double> ones(n, 1.0);
    const auto M = detail::sandwich_meat(S, ones, options.covariance == MultinomialCovariance::cluster ? std::span<const std::int64_t>(options.clusters) : std::span<const std::int64_t>{}, groups);
    fit.covariance = detail::sandwich(Hinv, M);
    if (options.covariance == MultinomialCovariance::cluster) {
      detail::require(groups >= 2, "at least two clusters");
      const double G = static_cast<double>(groups);
      for (double& v : fit.covariance.a) v *= G / (G - 1);
      fit.clusters = groups;
    }
  }
  const double crit = tyche::probability::Normal{}.upper_quantile((1 - options.confidence) / 2);
  for (std::size_t a = 0; a < q; a++) {
    const double se = std::sqrt(fit.covariance(a, a));
    fit.se.push_back(se);
    fit.z.push_back(fit.beta[a] / se);
    fit.p.push_back(tyche::probability::two_sided_normal(fit.beta[a] / se));
    fit.ci_lower.push_back(fit.beta[a] - crit * se);
    fit.ci_upper.push_back(fit.beta[a] + crit * se);
  }
  fit.log_likelihood = -dn * negative_mean_ll(fit.beta);
  for (std::size_t j = 0; j < J; j++)
    if (count[j] > 0) fit.null_log_likelihood += count[j] * std::log(count[j] / dn);
  fit.pseudo_r_squared = 1 - fit.log_likelihood / fit.null_log_likelihood;
  fit.aic = -2 * fit.log_likelihood + 2 * static_cast<double>(q);
  fit.bic = -2 * fit.log_likelihood + static_cast<double>(q) * std::log(dn);
  return fit;
}

struct MultinomialEffects {
  calx::linalg::Rectangular effect, se;  // J × k: the average effect of regressor l on the probability of category j
};

inline MultinomialEffects multinomial_marginal_effects(const MultinomialFit& fit, const Rectangular& X) {
  const std::size_t n = X.rows, k = X.cols, J = fit.categories, q = (J - 1) * k;
  detail::require(k == fit.k && n == fit.n, "the design the model was fitted on");
  std::vector<bool> constant(k, true);
  for (std::size_t l = 0; l < k; l++)
    for (std::size_t i = 1; i < n && constant[l]; i++) constant[l] = X(i, l) == X(0, l);
  // The average effect (j, l) at coefficients b.
  auto ame = [&](std::span<const double> b, std::size_t j, std::size_t l) {
    const auto P = detail::choice_probabilities(X, b, J);
    double s = 0;
    for (std::size_t i = 0; i < n; i++) {
      double mean_beta = 0;
      for (std::size_t m = 1; m < J; m++) mean_beta += P(i, m) * b[(m - 1) * k + l];
      const double bj = j == 0 ? 0.0 : b[(j - 1) * k + l];
      s += P(i, j) * (bj - mean_beta);
    }
    return s / static_cast<double>(n);
  };
  MultinomialEffects e;
  const double nan = std::numeric_limits<double>::quiet_NaN();
  e.effect = calx::linalg::Rectangular(J, k);
  e.se = calx::linalg::Rectangular(J, k);
  for (std::size_t j = 0; j < J; j++)
    for (std::size_t l = 0; l < k; l++) {
      if (constant[l]) {
        e.effect(j, l) = e.se(j, l) = nan;
        continue;
      }
      e.effect(j, l) = ame(fit.beta, j, l);
      const auto g = calx::optimisation::gradient_by_differences([&](std::span<const double> b) { return ame(b, j, l); }, fit.beta);
      double v = 0;
      for (std::size_t a = 0; a < q; a++)
        for (std::size_t c = 0; c < q; c++) v += g[a] * fit.covariance(a, c) * g[c];
      e.se(j, l) = std::sqrt(v);
    }
  return e;
}

}  // namespace nomos::econometrics
