// Choice among alternatives described by their attributes: the conditional logit, and willingness to pay (nomos's decision 0009).
//
// THE DATA, in long form. A choice situation s offers alternatives a ∈ A_s, each a row of attributes x_sa (its price, its travel time,
// whether it has a protected lane, …); exactly one is chosen. Situations may carry frequency weights w_s (a situation seen w_s times)
// and the decision maker who faced them (for clustering here; for the mixed logit's panel in mixed_logit.hpp). The rows of a
// situation are contiguous: situation s holds rows start[s] … start[s + 1] − 1.
//
// THE MODEL (McFadden 1974, "Conditional logit analysis of qualitative choice behavior", in Zarembka (ed.), Frontiers in Econometrics;
// Train 2009, Discrete Choice Methods with Simulation, 2nd ed., ch. 3). With utility U_sa = x_sa β + ε_sa and ε i.i.d. extreme value,
//   P_sa = exp(x_sa β) / Σ_{b ∈ A_s} exp(x_sb β).
// Only differences of utility matter: an attribute constant within every situation is not identified, and alternative-specific
// constants are columns of indicators with one alternative left out. The log-likelihood ℓ = Σ_s w_s ln P_s,chosen is concave, with
// gradient Σ_s w_s (x_s,chosen − Σ_a P_sa x_sa) and Hessian −Σ_s w_s Σ_a P_sa (x_sa − x̄_s)(x_sa − x̄_s)ᵀ, x̄_s = Σ_a P_sa x_sa. It is
// maximised by Newton's method from β = 0, with step halving.
//
// THE COVARIANCE: the inverse information (model), or the sandwich with the situations' scores, clustered by decision maker with the
// factor G/(G − 1). Inference is normal.
//
// WILLINGNESS TO PAY (Train 2009, §2.3; Hensher, Rose and Greene 2015, Applied Choice Analysis, ch. 12). With a price coefficient β_p,
// the money that leaves utility unchanged for a unit of attribute k is WTP_k = −β_k / β_p. Its standard error is by the delta method,
// with gradient (−1/β_p, β_k/β_p²) in (β_k, β_p); the ratio of two normals has no moments, so its interval is reported only where β_p
// is many standard errors from zero, and stated to be.

#pragma once

#include <tyche/probability/distributions.hpp>
#include <tyche/probability/p_values.hpp>

#include <calx/linalg/dense.hpp>
#include <calx/linalg/rectangular.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "common.hpp"

namespace nomos::econometrics {

struct ChoiceData {
  Rectangular X;                      // one row per alternative of each situation
  std::vector<std::size_t> start;     // situation s holds rows start[s] … start[s + 1] − 1; start.back() = rows
  std::vector<std::size_t> chosen;    // per situation: the row chosen (an index of X)
  std::vector<double> weight;         // per situation, frequency weights; empty for ones
  std::vector<std::int64_t> person;   // per situation, the decision maker; empty if each situation is its own
  std::size_t situations() const { return chosen.size(); }
  double w(std::size_t s) const { return weight.empty() ? 1.0 : weight[s]; }
};

enum class ChoiceCovariance { model, robust, cluster };

struct ChoiceFit {
  std::vector<double> beta, se, z, p, ci_lower, ci_upper;
  calx::linalg::Matrix covariance;
  double log_likelihood = 0, null_log_likelihood = 0, rho_squared = 0, aic = 0, bic = 0;
  std::vector<double> probabilities;  // per row, P_sa at β̂
  int iterations = 0;
  bool converged = false;
  std::size_t clusters = 0;
};

struct WillingnessToPay {
  double value = 0, se = 0, lower = 0, upper = 0;
  double price_z = 0;                 // β_p / se(β_p): the interval means little where it is small
};

namespace detail_choice {

inline void check(const ChoiceData& d) {
  detail::require(!d.start.empty() && d.start.size() == d.chosen.size() + 1 && d.start.back() == d.X.rows, "situations that cover the rows of X");
  for (std::size_t s = 0; s < d.situations(); s++) {
    detail::require(d.start[s + 1] > d.start[s] + 1, "two alternatives or more in every situation");
    detail::require(d.chosen[s] >= d.start[s] && d.chosen[s] < d.start[s + 1], "the chosen row within its situation");
  }
  detail::require(d.weight.empty() || d.weight.size() == d.situations(), "a weight per situation, or none");
  detail::require(d.person.empty() || d.person.size() == d.situations(), "a decision maker per situation, or none");
}

// P_sa for the rows of one situation at β (stable in the largest utility).
inline void probabilities(const ChoiceData& d, std::size_t s, std::span<const double> b, std::vector<double>& out) {
  const std::size_t a0 = d.start[s], a1 = d.start[s + 1], k = d.X.cols;
  out.assign(a1 - a0, 0.0);
  double top = -std::numeric_limits<double>::infinity();
  for (std::size_t a = a0; a < a1; a++) {
    double v = 0;
    for (std::size_t l = 0; l < k; l++) v += d.X(a, l) * b[l];
    out[a - a0] = v;
    top = std::max(top, v);
  }
  double z = 0;
  for (double& v : out) v = std::exp(v - top), z += v;
  for (double& v : out) v /= z;
}

inline double log_likelihood(const ChoiceData& d, std::span<const double> b) {
  std::vector<double> P;
  double ll = 0;
  for (std::size_t s = 0; s < d.situations(); s++) {
    probabilities(d, s, b, P);
    ll += d.w(s) * std::log(std::max(P[d.chosen[s] - d.start[s]], std::numeric_limits<double>::min()));
  }
  return ll;
}

}  // namespace detail_choice

inline ChoiceFit conditional_logit(const ChoiceData& d, ChoiceCovariance covariance = ChoiceCovariance::model, double confidence = 0.95, int max_iterations = 200) {
  namespace dc = detail_choice;
  dc::check(d);
  const std::size_t k = d.X.cols, S = d.situations();
  std::vector<double> b(k, 0.0), P;
  ChoiceFit fit;
  double ll = dc::log_likelihood(d, b);
  fit.null_log_likelihood = 0;
  for (std::size_t s = 0; s < S; s++) fit.null_log_likelihood -= d.w(s) * std::log(static_cast<double>(d.start[s + 1] - d.start[s]));
  calx::linalg::Matrix I(k);
  std::vector<double> g(k);
  const auto derivatives = [&](std::span<const double> beta) {
    std::fill(g.begin(), g.end(), 0.0);
    I = calx::linalg::Matrix(k);
    std::vector<double> mean(k);
    for (std::size_t s = 0; s < S; s++) {
      dc::probabilities(d, s, beta, P);
      std::fill(mean.begin(), mean.end(), 0.0);
      for (std::size_t a = d.start[s]; a < d.start[s + 1]; a++)
        for (std::size_t l = 0; l < k; l++) mean[l] += P[a - d.start[s]] * d.X(a, l);
      const double w = d.w(s);
      for (std::size_t l = 0; l < k; l++) g[l] += w * (d.X(d.chosen[s], l) - mean[l]);
      for (std::size_t a = d.start[s]; a < d.start[s + 1]; a++) {
        const double pa = P[a - d.start[s]];
        for (std::size_t l = 0; l < k; l++)
          for (std::size_t m = 0; m < k; m++) I(l, m) += w * pa * (d.X(a, l) - mean[l]) * (d.X(a, m) - mean[m]);
      }
    }
  };
  for (fit.iterations = 0; fit.iterations < max_iterations; fit.iterations++) {
    derivatives(b);
    calx::linalg::Matrix G(k);
    for (std::size_t l = 0; l < k; l++) G(l, 0) = g[l];
    const auto step = calx::linalg::solve(I, G);
    double t = 1, gnorm = 0;
    for (double v : g) gnorm = std::max(gnorm, std::fabs(v));
    if (gnorm < 1e-9 * std::max(1.0, std::fabs(ll))) { fit.converged = true; break; }
    std::vector<double> next(k);
    for (int h = 0; h < 40; h++, t /= 2) {
      for (std::size_t l = 0; l < k; l++) next[l] = b[l] + t * step(l, 0);
      const double llt = dc::log_likelihood(d, next);
      if (llt >= ll - 1e-12 * std::fabs(ll)) { ll = llt; break; }
    }
    b = next;
  }
  derivatives(b);
  fit.beta = b;
  fit.log_likelihood = ll;
  const auto Hinv = calx::linalg::solve(I, calx::linalg::Matrix::identity(k));
  if (covariance == ChoiceCovariance::model) {
    fit.covariance = Hinv;
  } else {
    // Scores by situation, summed within decision makers for clustering.
    std::map<std::int64_t, std::vector<double>> by;
    for (std::size_t s = 0; s < S; s++) {
      dc::probabilities(d, s, b, P);
      const std::int64_t key = (covariance == ChoiceCovariance::cluster && !d.person.empty()) ? d.person[s] : static_cast<std::int64_t>(s);
      auto& sc = by[key];
      sc.resize(k, 0.0);
      for (std::size_t l = 0; l < k; l++) {
        double mean = 0;
        for (std::size_t a = d.start[s]; a < d.start[s + 1]; a++) mean += P[a - d.start[s]] * d.X(a, l);
        sc[l] += d.w(s) * (d.X(d.chosen[s], l) - mean);
      }
    }
    calx::linalg::Matrix M(k);
    for (const auto& [key, sc] : by)
      for (std::size_t l = 0; l < k; l++)
        for (std::size_t m = 0; m < k; m++) M(l, m) += sc[l] * sc[m];
    fit.covariance = Hinv * M * Hinv;
    if (covariance == ChoiceCovariance::cluster) {
      const double G = static_cast<double>(by.size());
      detail::require(G >= 2, "at least two decision makers to cluster by");
      for (double& v : fit.covariance.a) v *= G / (G - 1);
      fit.clusters = by.size();
    }
  }
  const double crit = tyche::probability::Normal{}.upper_quantile((1 - confidence) / 2);
  for (std::size_t l = 0; l < k; l++) {
    const double se = std::sqrt(fit.covariance(l, l));
    fit.se.push_back(se), fit.z.push_back(b[l] / se), fit.p.push_back(tyche::probability::two_sided_normal(b[l] / se));
    fit.ci_lower.push_back(b[l] - crit * se), fit.ci_upper.push_back(b[l] + crit * se);
  }
  double n = 0;
  for (std::size_t s = 0; s < S; s++) n += d.w(s);
  fit.rho_squared = 1 - fit.log_likelihood / fit.null_log_likelihood;
  fit.aic = -2 * fit.log_likelihood + 2 * static_cast<double>(k);
  fit.bic = -2 * fit.log_likelihood + static_cast<double>(k) * std::log(n);
  fit.probabilities.assign(d.X.rows, 0.0);
  for (std::size_t s = 0; s < S; s++) {
    dc::probabilities(d, s, b, P);
    for (std::size_t a = d.start[s]; a < d.start[s + 1]; a++) fit.probabilities[a] = P[a - d.start[s]];
  }
  return fit;
}

// −β_k/β_p and its delta-method interval, from any fit's coefficients and covariance.
inline WillingnessToPay willingness_to_pay(std::span<const double> beta, const calx::linalg::Matrix& covariance, std::size_t k, std::size_t price, double confidence = 0.95) {
  detail::require(k < beta.size() && price < beta.size() && k != price, "an attribute and a price among the coefficients");
  const double bk = beta[k], bp = beta[price];
  WillingnessToPay w;
  w.value = -bk / bp;
  const double gk = -1 / bp, gp = bk / (bp * bp);
  const double var = gk * gk * covariance(k, k) + 2 * gk * gp * covariance(k, price) + gp * gp * covariance(price, price);
  w.se = std::sqrt(std::max(0.0, var));
  const double crit = tyche::probability::Normal{}.upper_quantile((1 - confidence) / 2);
  w.lower = w.value - crit * w.se, w.upper = w.value + crit * w.se;
  w.price_z = bp / std::sqrt(covariance(price, price));
  return w;
}

}  // namespace nomos::econometrics
