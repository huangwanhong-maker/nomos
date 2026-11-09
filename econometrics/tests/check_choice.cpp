// Checks of the models of choice among alternatives (conditional_logit.hpp, nested_logit.hpp, mixed_logit.hpp).
//
// - The conditional logit against statsmodels' ConditionalLogit on 300 situations of three alternatives (fixtures): β̂, its standard
//   errors and the log-likelihood; willingness to pay for a lane and its delta-method standard error against numpy.
// - A situation weighted 2 is the same as the situation twice; clustered by a decision maker per situation, the covariance is the
//   robust one times G/(G − 1).
// - The nested logit's likelihood with every λ = 1 is the conditional logit's, exactly; on 4,000 situations drawn from a nested logit
//   with λ = 0.5, β̂ and λ̂ lie within four standard errors.
// - The mixed logit on 300 people with 8 situations each, the price coefficient normal (−1, 0.5²): b̂ and |ŝ| within four standard
//   errors; with no random column it is the conditional logit.

#include <nomos/econometrics/conditional_logit.hpp>
#include <nomos/econometrics/mixed_logit.hpp>
#include <nomos/econometrics/nested_logit.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "fixtures/choice_fixtures.hpp"

namespace ne = nomos::econometrics;

namespace {
// Situations of J alternatives from per-alternative attribute rows.
ne::ChoiceData situations(const std::vector<std::vector<double>>& rows, std::size_t J, const std::vector<std::size_t>& chosen) {
  ne::ChoiceData d;
  const std::size_t k = rows.empty() ? 0 : rows[0].size();
  d.X = calx::linalg::Rectangular(rows.size(), k);
  for (std::size_t r = 0; r < rows.size(); r++)
    for (std::size_t l = 0; l < k; l++) d.X(r, l) = rows[r][l];
  for (std::size_t s = 0; s * J < rows.size(); s++) d.start.push_back(s * J), d.chosen.push_back(s * J + chosen[s]);
  d.start.push_back(rows.size());
  return d;
}

std::size_t draw(const std::vector<double>& p, calx::random::Stream& st) {
  double u = st.uniform(), c = 0;
  for (std::size_t j = 0; j < p.size(); j++) {
    c += p[j];
    if (u <= c) return j;
  }
  return p.size() - 1;
}
}  // namespace

int main() {
  calx::verification::Suite s;
  using namespace fixtures;
  {
    std::vector<std::vector<double>> rows;
    std::vector<std::size_t> chosen;
    for (std::size_t i = 0; i < choice_situations * choice_alternatives; i++) rows.push_back({choice_price[i], choice_time[i], choice_lane[i]});
    for (double c : choice_chosen) chosen.push_back(static_cast<std::size_t>(c));
    const auto d = situations(rows, choice_alternatives, chosen);
    const auto f = ne::conditional_logit(d);
    double eb = 0, es = 0;
    for (std::size_t l = 0; l < 3; l++) eb = std::max(eb, std::fabs(f.beta[l] - choice_beta[l]) / std::fabs(choice_beta[l])), es = std::max(es, std::fabs(f.se[l] - choice_se[l]) / choice_se[l]);
    s.check("conditional logit: β̂ against statsmodels' ConditionalLogit (relative)", eb, 1e-6);
    s.check("... its standard errors (relative)", es, 1e-5);
    s.check("... and the log-likelihood", std::fabs(f.log_likelihood - choice_loglik), 1e-8);
    const auto w = ne::willingness_to_pay(f.beta, f.covariance, 2, 0);
    s.check("willingness to pay for a lane, −β_lane/β_price, and its delta-method standard error against numpy",
            std::fabs(w.value - choice_wtp_lane) + std::fabs(w.se - choice_wtp_lane_se), 1e-6);
    // Weights as replication.
    auto twice = d, weighted = d;
    {
      calx::linalg::Rectangular X(d.X.rows + 3, 3);
      for (std::size_t i = 0; i < d.X.rows; i++)
        for (std::size_t l = 0; l < 3; l++) X(i, l) = d.X(i, l);
      for (std::size_t r = 0; r < 3; r++)
        for (std::size_t l = 0; l < 3; l++) X(d.X.rows + r, l) = d.X(r, l);
      twice.X = X;
      twice.start.back() = d.X.rows;
      twice.start.push_back(d.X.rows + 3);
      twice.chosen.push_back(d.X.rows + (d.chosen[0] - d.start[0]));
      weighted.weight.assign(d.situations(), 1.0);
      weighted.weight[0] = 2;
    }
    const auto ft = ne::conditional_logit(twice), fw = ne::conditional_logit(weighted);
    double gap = 0;
    for (std::size_t l = 0; l < 3; l++) gap = std::max(gap, std::fabs(ft.beta[l] - fw.beta[l]) + std::fabs(ft.se[l] - fw.se[l]));
    s.check("a situation weighted 2 is the situation twice: β̂ and standard errors", gap, 1e-9);
    auto people = d;
    for (std::size_t i = 0; i < d.situations(); i++) people.person.push_back(static_cast<std::int64_t>(i));
    const auto fr = ne::conditional_logit(people, ne::ChoiceCovariance::robust), fc = ne::conditional_logit(people, ne::ChoiceCovariance::cluster);
    const double G = static_cast<double>(d.situations());
    s.check("clustered by a decision maker per situation, the covariance is the robust one times G/(G − 1)", std::fabs(fc.covariance(0, 0) - fr.covariance(0, 0) * G / (G - 1)) / fr.covariance(0, 0), 1e-12);
  }
  {
    // Nested logit: four alternatives in two nests of two, λ = 0.5 in both, β = (−1, 0.8).
    calx::random::Stream st(20261109, 0x6e657374u, 0, 0);  // "nest"
    const std::size_t S = 4000, J = 4;
    const std::vector<std::size_t> nests{0, 0, 1, 1};
    const double b0 = -1, b1 = 0.8, lam = 0.5;
    std::vector<std::vector<double>> rows;
    std::vector<std::size_t> chosen;
    for (std::size_t sit = 0; sit < S; sit++) {
      std::vector<double> V(J);
      for (std::size_t j = 0; j < J; j++) {
        const double x0 = st.uniform(0, 3), x1 = static_cast<double>(st.below(2));
        rows.push_back({x0, x1});
        V[j] = b0 * x0 + b1 * x1;
      }
      double I[2] = {0, 0};
      for (std::size_t j = 0; j < J; j++) I[nests[j]] += std::exp(V[j] / lam);
      const double den = std::pow(I[0], lam) + std::pow(I[1], lam);
      std::vector<double> P(J);
      for (std::size_t j = 0; j < J; j++) P[j] = std::exp(V[j] / lam) * std::pow(I[nests[j]], lam - 1) / den;
      chosen.push_back(draw(P, st));
    }
    const auto d = situations(rows, J, chosen);
    std::vector<std::size_t> nest(d.X.rows);
    for (std::size_t r = 0; r < d.X.rows; r++) nest[r] = nests[r % J];
    const std::vector<double> ones{1.0, 1.0};
    const std::vector<double> bb{b0, b1};
    s.check("nested logit: with every λ = 1 its likelihood is the conditional logit's",
            std::fabs(ne::detail_nested::log_likelihood(d, nest, 2, bb, ones) - ne::detail_choice::log_likelihood(d, bb)), 1e-8);
    const auto f = ne::nested_logit(d, nest);
    const double worst = std::max({std::fabs(f.beta[0] - b0) / f.se[0], std::fabs(f.beta[1] - b1) / f.se[1], std::fabs(f.lambda[0] - lam) / f.lambda_se[0], std::fabs(f.lambda[1] - lam) / f.lambda_se[1]});
    s.check("... on 4,000 situations drawn with λ = 0.5: β̂ and λ̂ within four standard errors (largest |z|)", worst, 4);
  }
  {
    // Mixed logit: 300 people, 8 situations each of three alternatives; price ~ N(−1, 0.5²), time −0.05, lane 0.6.
    calx::random::Stream st(20261109, 0x6d697864u, 0, 0);  // "mixd"
    const std::size_t N = 300, T = 8, J = 3;
    std::vector<std::vector<double>> rows;
    std::vector<std::size_t> chosen;
    std::vector<std::int64_t> person;
    for (std::size_t n = 0; n < N; n++) {
      const double bp = -1 + 0.5 * st.normal();
      for (std::size_t t = 0; t < T; t++) {
        std::vector<double> V(J);
        for (std::size_t j = 0; j < J; j++) {
          const double p = st.uniform(1, 5), tm = st.uniform(5, 40), ln = static_cast<double>(st.below(2));
          rows.push_back({p, tm, ln});
          V[j] = bp * p - 0.05 * tm + 0.6 * ln;
        }
        double top = *std::max_element(V.begin(), V.end()), z = 0;
        std::vector<double> P(J);
        for (std::size_t j = 0; j < J; j++) P[j] = std::exp(V[j] - top), z += P[j];
        for (double& v : P) v /= z;
        chosen.push_back(draw(P, st));
        person.push_back(static_cast<std::int64_t>(n));
      }
    }
    auto d = situations(rows, J, chosen);
    d.person = person;
    const std::vector<std::size_t> random{0};
    const auto f = ne::mixed_logit(d, random, 200);
    const double worst = std::max({std::fabs(f.mean[0] + 1) / f.mean_se[0], std::fabs(f.mean[1] + 0.05) / f.mean_se[1], std::fabs(f.mean[2] - 0.6) / f.mean_se[2],
                                   std::fabs(f.spread[0] - 0.5) / f.spread_se[0]});
    s.check("mixed logit, 300 people of 8 situations, the price coefficient N(−1, 0.5²): b̂ and |ŝ| within four standard errors (largest |z|)", worst, 4);
    s.require("... and it fits better than the conditional logit (simulated log-likelihood)", f.log_likelihood > ne::conditional_logit(d).log_likelihood);
    const auto plain = ne::mixed_logit(d, std::vector<std::size_t>{}, 10);
    const auto cl = ne::conditional_logit(d);
    double gap = 0;
    for (std::size_t l = 0; l < 3; l++) gap = std::max(gap, std::fabs(plain.mean[l] - cl.beta[l]));
    s.check("... with no random column it is the conditional logit", gap, 1e-5);
  }
  return s.done();
}
