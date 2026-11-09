// The nested logit: alternatives in nests share unobserved utility (nomos's decision 0009).
//
// THE MODEL (McFadden 1978, "Modelling the choice of residential location"; Train 2009, ch. 4). Each alternative belongs to a nest m;
// within a nest the errors are correlated, which the dissimilarity λ_m ∈ (0, 1] measures (1 − λ_m is roughly their correlation). With
// V_sa = x_sa β and I_sm = Σ_{b ∈ A_s ∩ m} exp(V_sb/λ_m),
//   P_sa = exp(V_sa/λ_m) I_sm^{λ_m − 1} / Σ_l I_sl^{λ_l},
// the probability of the nest, I_sm^{λ_m}/Σ_l I_sl^{λ_l}, times the conditional logit within it. With every λ_m = 1 it is the
// conditional logit (checked). A nest with one alternative has no λ to estimate; it is held at 1. λ_m outside (0, 1] is reported, not
// forbidden: it says the nesting is not consistent with random utility for all values of the attributes (Train 2009, §4.2.1).
//
// THE ESTIMATE. The log-likelihood is not concave in (β, λ). It is maximised by calx's BFGS from the conditional logit's β and λ = 1,
// with gradients by central differences. The covariance is the inverse of the Hessian by differences. Inference is normal.

#pragma once

#include <calx/optimisation/quasi_newton.hpp>

#include "conditional_logit.hpp"

namespace nomos::econometrics {

struct NestedFit {
  std::vector<double> beta, se, z, p;         // the attributes' coefficients
  std::vector<double> lambda, lambda_se;      // per nest (1 and NaN for a nest of one alternative)
  calx::linalg::Matrix covariance;            // of (β, the λ estimated)
  double log_likelihood = 0, aic = 0;
  int iterations = 0;
  bool converged = false;
};

namespace detail_nested {

// The log-likelihood at (β, λ); `nest` per row, `M` nests.
inline double log_likelihood(const ChoiceData& d, std::span<const std::size_t> nest, std::size_t M, std::span<const double> beta, std::span<const double> lambda) {
  const std::size_t k = d.X.cols;
  double ll = 0;
  std::vector<double> V, top(M), I(M);
  for (std::size_t s = 0; s < d.situations(); s++) {
    const std::size_t a0 = d.start[s], a1 = d.start[s + 1];
    V.assign(a1 - a0, 0.0);
    std::fill(top.begin(), top.end(), -std::numeric_limits<double>::infinity());
    std::fill(I.begin(), I.end(), 0.0);
    for (std::size_t a = a0; a < a1; a++) {
      double v = 0;
      for (std::size_t l = 0; l < k; l++) v += d.X(a, l) * beta[l];
      V[a - a0] = v;
      top[nest[a]] = std::max(top[nest[a]], v / lambda[nest[a]]);
    }
    for (std::size_t a = a0; a < a1; a++) I[nest[a]] += std::exp(V[a - a0] / lambda[nest[a]] - top[nest[a]]);
    // ln Σ_l I_l^{λ_l} with I_l = e^{top_l} Ĩ_l.
    double big = -std::numeric_limits<double>::infinity();
    std::vector<double> inc(M, -std::numeric_limits<double>::infinity());
    for (std::size_t m = 0; m < M; m++)
      if (I[m] > 0) inc[m] = lambda[m] * (top[m] + std::log(I[m])), big = std::max(big, inc[m]);
    double z = 0;
    for (std::size_t m = 0; m < M; m++)
      if (I[m] > 0) z += std::exp(inc[m] - big);
    const std::size_t c = d.chosen[s], m = nest[c];
    const double lnP = V[c - a0] / lambda[m] + (lambda[m] - 1) * (top[m] + std::log(I[m])) - (big + std::log(z));
    ll += d.w(s) * lnP;
  }
  return ll;
}

}  // namespace detail_nested

inline NestedFit nested_logit(const ChoiceData& d, std::span<const std::size_t> nest, double confidence = 0.95) {
  detail_choice::check(d);
  detail::require(nest.size() == d.X.rows, "a nest for every row");
  std::size_t M = 0;
  for (std::size_t m : nest) M = std::max(M, m + 1);
  // Nests with more than one alternative in some situation have a λ to estimate.
  std::vector<char> free(M, 0);
  for (std::size_t s = 0; s < d.situations(); s++) {
    std::vector<int> count(M, 0);
    for (std::size_t a = d.start[s]; a < d.start[s + 1]; a++)
      if (++count[nest[a]] > 1) free[nest[a]] = 1;
  }
  std::vector<std::size_t> which;
  for (std::size_t m = 0; m < M; m++)
    if (free[m]) which.push_back(m);
  const std::size_t k = d.X.cols, q = k + which.size();
  const auto unpack = [&](std::span<const double> th, std::vector<double>& lambda) {
    lambda.assign(M, 1.0);
    for (std::size_t j = 0; j < which.size(); j++) lambda[which[j]] = th[k + j];
  };
  double wsum = 0;
  for (std::size_t s = 0; s < d.situations(); s++) wsum += d.w(s);
  const auto f = [&](std::span<const double> th) {
    std::vector<double> lambda;
    unpack(th, lambda);
    for (double l : lambda)
      if (!(l > 0.02)) return 1e30;
    return -detail_nested::log_likelihood(d, nest, M, th.subspan(0, k), lambda) / wsum;
  };
  const auto start = conditional_logit(d);
  std::vector<double> th(start.beta);
  th.resize(q, 1.0);
  calx::optimisation::BfgsOptions bo;
  bo.gradient_tolerance = 1e-9;
  const auto res = calx::optimisation::bfgs(f, [&](std::span<const double> x) { return calx::optimisation::gradient_by_differences(f, x); }, th, bo);
  NestedFit fit;
  fit.iterations = res.iterations, fit.converged = res.converged;
  th = res.x;
  auto H = calx::optimisation::hessian_by_differences(f, th);
  for (double& v : H.a) v *= wsum;  // the information of the sum
  fit.covariance = calx::linalg::solve(H, calx::linalg::Matrix::identity(q));
  fit.beta.assign(th.begin(), th.begin() + static_cast<std::ptrdiff_t>(k));
  for (std::size_t l = 0; l < k; l++) {
    const double se = std::sqrt(std::max(0.0, fit.covariance(l, l)));
    fit.se.push_back(se), fit.z.push_back(fit.beta[l] / se), fit.p.push_back(tyche::probability::two_sided_normal(fit.beta[l] / se));
  }
  unpack(th, fit.lambda);
  fit.lambda_se.assign(M, std::numeric_limits<double>::quiet_NaN());
  for (std::size_t j = 0; j < which.size(); j++) fit.lambda_se[which[j]] = std::sqrt(std::max(0.0, fit.covariance(k + j, k + j)));
  fit.log_likelihood = -f(th) * wsum;
  fit.aic = -2 * fit.log_likelihood + 2 * static_cast<double>(q);
  (void)confidence;
  return fit;
}

}  // namespace nomos::econometrics
