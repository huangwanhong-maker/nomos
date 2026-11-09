// The mixed logit: tastes that vary between decision makers (nomos's decision 0009).
//
// THE MODEL (Revelt and Train 1998, Rev. Econ. Stat. 80, 647–657; McFadden and Train 2000, J. Appl. Econometrics 15, 447–470; Train
// 2009, ch. 6). Decision maker n has coefficients β_n: for the random columns, β_nk = b_k + s_k η_nk with η_nk ~ N(0, 1) independent;
// the other columns are fixed, β_nk = b_k. Given β_n, its choices in its situations are conditional logits, so its probability is
//   P_n = ∫ Π_{s ∈ n} P_s(β) φ(η) dη,
// simulated by R draws: P̌_n = (1/R) Σ_r Π_s P_s(b + s ⊙ η_nr). The simulated log-likelihood Σ_n w_n ln P̌_n is maximised by calx's BFGS
// with its analytic gradient, ∂ ln P̌_n/∂θ = Σ_r Π_s P_s(β_nr) Σ_s ∂ ln P_s(β_nr)/∂θ / Σ_r Π_s P_s(β_nr), where ∂ ln P_s/∂β =
// x_s,chosen − Σ_a P_sa x_sa, ∂β/∂b = 1 and ∂β_k/∂s_k = η_k. The sign of s_k is not identified; |s_k| is reported.
//
// THE DRAWS. Halton sequences (Train 2009, §9.3.3): one sequence per random column, in prime bases 2, 3, 5, …, its first ten points
// dropped, consecutive blocks of R given to consecutive decision makers, taken to normals by the inverse of Φ. With R draws the
// simulated likelihood is biased for finite R (ln is concave); the bias falls as R grows (Train 2009, §10.5), and R is the caller's.
//
// THE COVARIANCE: the inverse of the Hessian of the simulated log-likelihood, by differences of its analytic gradient.

#pragma once

#include <calx/optimisation/quasi_newton.hpp>

#include "conditional_logit.hpp"

namespace nomos::econometrics {

struct MixedFit {
  std::vector<double> mean, mean_se;          // b, per column
  std::vector<double> spread, spread_se;      // |s|, per random column (in the order given)
  std::vector<std::size_t> random;            // the random columns
  calx::linalg::Matrix covariance;            // of (b, s)
  double log_likelihood = 0, aic = 0;
  std::size_t draws = 0, people = 0;
  int iterations = 0;
  bool converged = false;
};

namespace detail_mixed {

// The radical inverse of i in a base.
inline double halton(std::size_t i, std::size_t base) {
  double f = 1, r = 0;
  while (i > 0) {
    f /= static_cast<double>(base);
    r += f * static_cast<double>(i % base);
    i /= base;
  }
  return r;
}

inline std::size_t prime(std::size_t j) {
  static const std::size_t P[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53};
  detail::require(j < 16, "at most sixteen random columns");
  return P[j];
}

}  // namespace detail_mixed

// The mixed logit with normal coefficients on the columns `random`, R Halton draws per decision maker.
inline MixedFit mixed_logit(const ChoiceData& d, std::span<const std::size_t> random, std::size_t R = 500) {
  detail_choice::check(d);
  const std::size_t k = d.X.cols, K = random.size(), q = k + K, S = d.situations();
  for (std::size_t r : random) detail::require(r < k, "random columns among X's");
  // Decision makers and their situations.
  std::vector<std::vector<std::size_t>> people;
  {
    std::map<std::int64_t, std::size_t> index;
    for (std::size_t s = 0; s < S; s++) {
      const std::int64_t key = d.person.empty() ? static_cast<std::int64_t>(s) : d.person[s];
      auto [it, fresh] = index.emplace(key, people.size());
      if (fresh) people.emplace_back();
      people[it->second].push_back(s);
    }
  }
  const std::size_t N = people.size();
  // η_nr, dimension K, from Halton sequences.
  std::vector<double> eta(N * R * K);
  const tyche::probability::Normal phi{};
  for (std::size_t j = 0; j < K; j++)
    for (std::size_t t = 0; t < N * R; t++) eta[t * K + j] = phi.quantile(detail_mixed::halton(t + 11, detail_mixed::prime(j)));
  double wsum = 0;
  for (std::size_t s = 0; s < S; s++) wsum += d.w(s);
  // The negative mean simulated log-likelihood and its gradient, together.
  const auto evaluate = [&](std::span<const double> th, std::vector<double>* grad) {
    double total = 0;
    if (grad) grad->assign(q, 0.0);
    std::vector<double> beta(k), P, g(q), acc(q);
    for (std::size_t n = 0; n < N; n++) {
      double sum = 0;
      std::fill(acc.begin(), acc.end(), 0.0);
      // Products over a person's situations can underflow: carried in logarithms, shifted by the largest.
      std::vector<double> logprod(R), dl(R * q);
      for (std::size_t r = 0; r < R; r++) {
        const double* e = &eta[(n * R + r) * K];
        for (std::size_t l = 0; l < k; l++) beta[l] = th[l];
        for (std::size_t j = 0; j < K; j++) beta[random[j]] += th[k + j] * e[j];
        double lp = 0;
        std::fill(g.begin(), g.end(), 0.0);
        for (const std::size_t s : people[n]) {
          detail_choice::probabilities(d, s, beta, P);
          lp += d.w(s) * std::log(std::max(P[d.chosen[s] - d.start[s]], std::numeric_limits<double>::min()));
          if (grad)
            for (std::size_t l = 0; l < k; l++) {
              double mean = 0;
              for (std::size_t a = d.start[s]; a < d.start[s + 1]; a++) mean += P[a - d.start[s]] * d.X(a, l);
              g[l] += d.w(s) * (d.X(d.chosen[s], l) - mean);
            }
        }
        logprod[r] = lp;
        if (grad) {
          for (std::size_t j = 0; j < K; j++) g[k + j] = g[random[j]] * e[j];
          for (std::size_t l = 0; l < q; l++) dl[r * q + l] = g[l];
        }
      }
      double top = -std::numeric_limits<double>::infinity();
      for (double v : logprod) top = std::max(top, v);
      for (std::size_t r = 0; r < R; r++) {
        const double w = std::exp(logprod[r] - top);
        sum += w;
        if (grad)
          for (std::size_t l = 0; l < q; l++) acc[l] += w * dl[r * q + l];
      }
      total += top + std::log(sum / static_cast<double>(R));
      if (grad)
        for (std::size_t l = 0; l < q; l++) (*grad)[l] -= acc[l] / sum / wsum;
    }
    return -total / wsum;
  };
  const auto f = [&](std::span<const double> th) { return evaluate(th, nullptr); };
  const auto gradient = [&](std::span<const double> th) {
    std::vector<double> g;
    evaluate(th, &g);
    return g;
  };
  const auto start = conditional_logit(d);
  std::vector<double> th(start.beta);
  for (std::size_t j = 0; j < K; j++) th.push_back(0.5 * std::max(0.1, std::fabs(start.beta[random[j]])));
  calx::optimisation::BfgsOptions bo;
  bo.gradient_tolerance = 1e-8;
  const auto res = calx::optimisation::bfgs(f, gradient, th, bo);
  th = res.x;
  // The Hessian by central differences of the gradient.
  calx::linalg::Matrix H(q);
  for (std::size_t c = 0; c < q; c++) {
    const double h = 1e-5 * std::max(1.0, std::fabs(th[c]));
    auto up = th, down = th;
    up[c] += h, down[c] -= h;
    const auto gu = gradient(up), gd = gradient(down);
    for (std::size_t r = 0; r < q; r++) H(r, c) = (gu[r] - gd[r]) / (2 * h) * wsum;
  }
  for (std::size_t r = 0; r < q; r++)
    for (std::size_t c = r + 1; c < q; c++) H(r, c) = H(c, r) = 0.5 * (H(r, c) + H(c, r));
  MixedFit fit;
  fit.covariance = calx::linalg::solve(H, calx::linalg::Matrix::identity(q));
  fit.random.assign(random.begin(), random.end());
  fit.draws = R, fit.people = N, fit.iterations = res.iterations, fit.converged = res.converged;
  for (std::size_t l = 0; l < k; l++) fit.mean.push_back(th[l]), fit.mean_se.push_back(std::sqrt(std::max(0.0, fit.covariance(l, l))));
  for (std::size_t j = 0; j < K; j++) fit.spread.push_back(std::fabs(th[k + j])), fit.spread_se.push_back(std::sqrt(std::max(0.0, fit.covariance(k + j, k + j))));
  fit.log_likelihood = -f(th) * wsum;
  fit.aic = -2 * fit.log_likelihood + 2 * static_cast<double>(q);
  return fit;
}

}  // namespace nomos::econometrics
