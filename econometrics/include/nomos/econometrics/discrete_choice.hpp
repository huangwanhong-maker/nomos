// Discrete choice: average marginal effects of binary logit and probit models.
//
// A binary choice model P(y = 1 | x) = F(xβ) is fitted by tyche's GLM (binomial family, logit
// link F = Λ, Λ(η) = 1/(1 + e^{−η}), or probit link F = Φ). Its coefficients are not effects on the
// probability; the average marginal effect of regressor j is (Wooldridge 2010, ch. 15; Cameron and
// Trivedi 2005, ch. 14):
//
//   continuous j:  AME_j = (1/n) Σ_i f(η_i) β_j,           f = F′;
//   discrete j:    AME_j = (1/n) Σ_i [F(η_i⁽¹⁾) − F(η_i⁽⁰⁾)],  η_i⁽ᵛ⁾ the index with x_ij set to v ∈ {0, 1},
//
// with η_i the fitted linear predictor (offset included). Their covariance is by the delta method,
// G V Gᵀ with V the fit's covariance of β̂ (model-based or robust, as fitted) and the gradient
//
//   continuous:  ∂AME_j/∂β_l = (1/n) Σ_i [f(η_i) 1(j = l) + β_j f′(η_i) x_il],
//   discrete:    ∂AME_j/∂β_l = (1/n) Σ_i [f(η_i⁽¹⁾) x_il⁽¹⁾ − f(η_i⁽⁰⁾) x_il⁽⁰⁾],
//
// where for the logit f = Λ(1 − Λ) and f′ = Λ(1 − Λ)(1 − 2Λ), and for the probit f = φ and
// f′(η) = −η φ(η). Inference is normal. A constant column (the intercept) has no marginal effect
// (NaN). This is what Stata's margins, dydx(*) reports at its default (average) settings.

#pragma once

#include <tyche/probability/distributions.hpp>
#include <tyche/probability/p_values.hpp>
#include <tyche/regression/glm.hpp>

#include <cmath>
#include <limits>
#include <span>
#include <vector>

#include "common.hpp"

namespace nomos::econometrics {

namespace tr = tyche::regression;
using calx::linalg::Rectangular;

struct MarginalEffects {
  std::vector<double> effect, se, z, p, ci_lower, ci_upper;  // one per column of X; NaN for constant columns
  std::vector<bool> discrete;
  calx::linalg::Matrix covariance;                             // G V Gᵀ
};

namespace detail {
struct Choice {
  tr::Link link;
  double F(double eta) const {
    if (link == tr::Link::logit) return 1 / (1 + std::exp(-eta));
    return tyche::probability::Normal{}.cdf(eta);
  }
  double f(double eta) const {
    if (link == tr::Link::logit) {
      const double L = F(eta);
      return L * (1 - L);
    }
    return tyche::probability::Normal{}.pdf(eta);
  }
  double fprime(double eta) const {
    if (link == tr::Link::logit) {
      const double L = F(eta);
      return L * (1 - L) * (1 - 2 * L);
    }
    return -eta * tyche::probability::Normal{}.pdf(eta);
  }
};
}  // namespace detail

// The average marginal effects of a binary logit or probit fit on design X; `discrete` lists the
// columns of 0/1 regressors to be treated by differences.
inline MarginalEffects average_marginal_effects(const tr::GlmFit& fit, const Rectangular& X, std::span<const std::size_t> discrete = {}, double confidence = 0.95) {
  detail::require(fit.family == tr::Family::binomial && (fit.link == tr::Link::logit || fit.link == tr::Link::probit), "a binomial fit with the logit or probit link");
  const std::size_t n = X.rows, k = X.cols;
  detail::require(fit.k == k && fit.linear_predictor.size() == n, "the design the model was fitted on");
  const detail::Choice c{fit.link};
  MarginalEffects m;
  m.discrete.assign(k, false);
  for (std::size_t j : discrete) {
    detail::require(j < k, "a discrete column within range");
    m.discrete[j] = true;
  }
  const double nan = std::numeric_limits<double>::quiet_NaN();
  std::vector<bool> constant(k, true);
  for (std::size_t j = 0; j < k; j++)
    for (std::size_t i = 1; i < n && constant[j]; i++) constant[j] = X(i, j) == X(0, j);
  Rectangular G(k, k);
  m.effect.assign(k, nan);
  const double dn = static_cast<double>(n);
  for (std::size_t j = 0; j < k; j++) {
    if (constant[j] || fit.aliased[j]) continue;
    double ame = 0;
    for (std::size_t i = 0; i < n; i++) {
      const double eta = fit.linear_predictor[i];
      if (!m.discrete[j]) {
        ame += c.f(eta) * fit.beta[j];
        for (std::size_t l = 0; l < k; l++)
          if (!fit.aliased[l]) G(j, l) += (fit.beta[j] * c.fprime(eta) * X(i, l) + (j == l ? c.f(eta) : 0.0)) / dn;
      } else {
        const double e1 = eta + (1 - X(i, j)) * fit.beta[j], e0 = eta - X(i, j) * fit.beta[j];
        ame += c.F(e1) - c.F(e0);
        for (std::size_t l = 0; l < k; l++) {
          if (fit.aliased[l]) continue;
          const double x1 = l == j ? 1.0 : X(i, l), x0 = l == j ? 0.0 : X(i, l);
          G(j, l) += (c.f(e1) * x1 - c.f(e0) * x0) / dn;
        }
      }
    }
    m.effect[j] = ame / dn;
  }
  // G V Gᵀ over the estimated coefficients.
  m.covariance = calx::linalg::Matrix(k);
  for (std::size_t a = 0; a < k; a++)
    for (std::size_t b = 0; b < k; b++) {
      double s = 0;
      for (std::size_t l = 0; l < k; l++) {
        if (fit.aliased[l] || G(a, l) == 0) continue;
        for (std::size_t r = 0; r < k; r++)
          if (!fit.aliased[r]) s += G(a, l) * fit.covariance(l, r) * G(b, r);
      }
      m.covariance(a, b) = s;
    }
  const double crit = tyche::probability::Normal{}.upper_quantile((1 - confidence) / 2);
  m.se = m.z = m.p = m.ci_lower = m.ci_upper = std::vector<double>(k, nan);
  for (std::size_t j = 0; j < k; j++) {
    if (std::isnan(m.effect[j])) continue;
    m.se[j] = std::sqrt(m.covariance(j, j));
    m.z[j] = m.effect[j] / m.se[j];
    m.p[j] = tyche::probability::two_sided_normal(m.z[j]);
    m.ci_lower[j] = m.effect[j] - crit * m.se[j];
    m.ci_upper[j] = m.effect[j] + crit * m.se[j];
  }
  return m;
}

}  // namespace nomos::econometrics
