// Checks of average marginal effects (discrete_choice.hpp) on simulated logit and probit data.
//
// For a continuous regressor, the effect is the derivative of the average predicted probability
// (by central differences); for a binary one, the average difference of the probabilities at 1
// and at 0. The delta-method covariance against the gradient taken by central differences in β.
// The intercept has no effect (NaN). The effects at the true coefficients lie within four standard
// errors of the estimates. A probit alike.

#include <nomos/econometrics/discrete_choice.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <string>
#include <vector>

namespace ne = nomos::econometrics;
namespace tr = tyche::regression;
using calx::linalg::Rectangular;
using calx::verification::relative;

namespace {

calx::random::Stream z(20261025, 0x6c6f6769u, 0, 0);  // "logi"

double F(tr::Link l, double e) { return l == tr::Link::logit ? 1 / (1 + std::exp(-e)) : 0.5 * std::erfc(-e / std::sqrt(2.0)); }

double f(tr::Link l, double e) {
  if (l == tr::Link::logit) return F(l, e) * (1 - F(l, e));
  return std::exp(-0.5 * e * e) / std::sqrt(2 * 3.14159265358979323846);
}

// The average marginal effects at β, from their definitions (column 2 binary).
std::vector<double> ames(tr::Link l, const Rectangular& X, const std::vector<double>& b) {
  const std::size_t n = X.rows;
  std::vector<double> out(3, 0.0);
  for (std::size_t i = 0; i < n; i++) {
    out[1] += b[1] * f(l, b[0] + b[1] * X(i, 1) + b[2] * X(i, 2));
    out[2] += F(l, b[0] + b[1] * X(i, 1) + b[2]) - F(l, b[0] + b[1] * X(i, 1));
  }
  for (double& v : out) v /= static_cast<double>(n);
  return out;
}

// The derivative of the average probability in x₁, by central differences.
double slope(tr::Link l, const Rectangular& X, const std::vector<double>& b, double h = 1e-5) {
  double s = 0;
  for (std::size_t i = 0; i < X.rows; i++) {
    const double e = b[0] + b[1] * X(i, 1) + b[2] * X(i, 2);
    s += (F(l, e + b[1] * h) - F(l, e - b[1] * h)) / (2 * h);
  }
  return s / static_cast<double>(X.rows);
}

}  // namespace

int main() {
  calx::verification::Suite s;
  for (tr::Link link : {tr::Link::logit, tr::Link::probit}) {
    const std::string name = link == tr::Link::logit ? "logit" : "probit";
    const std::size_t n = 3000;
    const std::vector<double> truth{-0.5, 0.8, -0.6};
    Rectangular X(n, 3);
    std::vector<double> y(n);
    for (std::size_t i = 0; i < n; i++) {
      X(i, 0) = 1;
      X(i, 1) = z.normal();
      X(i, 2) = z.uniform() < 0.4 ? 1 : 0;
      y[i] = z.uniform() < F(link, truth[0] + truth[1] * X(i, 1) + truth[2] * X(i, 2)) ? 1 : 0;
    }
    tr::GlmOptions o;
    o.family = tr::Family::binomial;
    o.link = link;
    o.covariance = tr::GlmCovariance::hc0;
    const auto fit = tr::glm(X, y, o);
    const std::vector<std::size_t> binary{2};
    const auto m = ne::average_marginal_effects(fit, X, binary);
    const auto direct = ames(link, X, fit.beta);
    s.check(name + " (n = 3000): the continuous effect is the derivative of the average probability", std::fabs(m.effect[1] - slope(link, X, fit.beta)), 1e-9);
    s.check("... and β₁ times the average density", std::fabs(m.effect[1] - direct[1]), 1e-14);
    s.check("... the binary effect is the average difference of probabilities at 1 and 0", std::fabs(m.effect[2] - direct[2]), 1e-12);
    // The delta method against central differences in β.
    const double h = 1e-5;
    double worst = 0;
    std::vector<std::vector<double>> G(3, std::vector<double>(3));
    for (std::size_t l = 0; l < 3; l++) {
      auto bp = fit.beta, bm = fit.beta;
      bp[l] += h;
      bm[l] -= h;
      const auto ap = ames(link, X, bp), am = ames(link, X, bm);
      for (std::size_t j = 1; j < 3; j++) G[j][l] = (ap[j] - am[j]) / (2 * h);
    }
    for (std::size_t j = 1; j < 3; j++) {
      double v = 0;
      for (std::size_t a = 0; a < 3; a++)
        for (std::size_t b = 0; b < 3; b++) v += G[j][a] * fit.covariance(a, b) * G[j][b];
      worst = std::max(worst, relative(m.se[j], std::sqrt(v)));
    }
    s.check("... the delta-method standard errors against a gradient by central differences (relative)", worst, 1e-5);
    s.require("... the intercept has no marginal effect", std::isnan(m.effect[0]) && std::isnan(m.se[0]));
    const auto at_truth = ames(link, X, truth);
    s.require("... the effects at the true coefficients lie within four standard errors",
              std::fabs(m.effect[1] - at_truth[1]) < 4 * m.se[1] && std::fabs(m.effect[2] - at_truth[2]) < 4 * m.se[2],
              calx::verification::exponential(m.effect[1]) + " ± " + calx::verification::exponential(m.se[1]));
  }
  return s.done();
}
