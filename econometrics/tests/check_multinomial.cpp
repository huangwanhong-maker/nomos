// Checks of the multinomial logit (multinomial.hpp) on simulated choices.
//
// With two categories it is the binary logit: the coefficients and model-based standard errors of
// tyche's GLM, and the average marginal effects and their standard errors of
// discrete_choice.hpp. With three categories: the true coefficients within four standard errors;
// the mean predicted probability of each category equal to its share (the intercepts' first-order
// condition); the average marginal effects summing to zero over the categories and equal to the
// derivatives of the mean predicted probabilities by central differences.

#include <nomos/econometrics/discrete_choice.hpp>
#include <nomos/econometrics/multinomial.hpp>

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
calx::random::Stream z(20261025, 0x6d6e6c67u, 0, 0);  // "mnlg"

// Choices from the multinomial logit with coefficients b ((J − 1) blocks of 3) on X = [1, x1, x2].
std::vector<std::size_t> choose(const Rectangular& X, const std::vector<double>& b, std::size_t J) {
  const auto P = ne::detail::choice_probabilities(X, b, J);
  std::vector<std::size_t> y(X.rows);
  for (std::size_t i = 0; i < X.rows; i++) {
    double u = z.uniform(), c = 0;
    std::size_t j = 0;
    for (; j + 1 < J; j++) {
      c += P(i, j);
      if (u < c) break;
    }
    y[i] = j;
  }
  return y;
}

Rectangular design(std::size_t n) {
  Rectangular X(n, 3);
  for (std::size_t i = 0; i < n; i++) {
    X(i, 0) = 1;
    X(i, 1) = z.normal();
    X(i, 2) = z.uniform(-1, 1);
  }
  return X;
}
}  // namespace

int main() {
  calx::verification::Suite s;

  // Two categories: the binary logit.
  {
    const auto X = design(2000);
    const auto y = choose(X, {-0.3, 0.9, -0.7}, 2);
    const auto m = ne::multinomial_logit(X, y, 2);
    std::vector<double> yd(y.begin(), y.end());
    tr::GlmOptions o;
    o.family = tr::Family::binomial;
    o.tolerance = 1e-14;
    const auto g = tr::glm(X, yd, o);
    double db = 0, dse = 0;
    for (std::size_t l = 0; l < 3; l++) {
      db = std::max(db, std::fabs(m.beta[l] - g.beta[l]));
      dse = std::max(dse, relative(m.se[l], g.se[l]));
    }
    s.check("two categories: the coefficients of the binary logit (tyche's GLM)", db, 1e-6);
    s.check("... its model-based standard errors (relative)", dse, 1e-5);
    s.require("... converged", m.converged, std::to_string(m.iterations) + " iterations");
    const auto e = ne::multinomial_marginal_effects(m, X);
    const auto b = ne::average_marginal_effects(g, X);
    double de = 0, des = 0;
    for (std::size_t l = 1; l < 3; l++) {
      de = std::max(de, std::fabs(e.effect(1, l) - b.effect[l]));
      des = std::max(des, relative(e.se(1, l), b.se[l]));
    }
    s.check("... and the binary average marginal effects", de, 1e-6);
    s.check("... with their delta-method standard errors (relative)", des, 1e-4);
  }
  // Three categories.
  {
    const std::vector<double> truth{0.5, 1.0, -0.5, -0.3, -0.8, 0.7};
    const auto X = design(4000);
    const auto y = choose(X, truth, 3);
    const auto m = ne::multinomial_logit(X, y, 3);
    bool within = true;
    for (std::size_t a = 0; a < 6; a++) within = within && std::fabs(m.beta[a] - truth[a]) < 4 * m.se[a];
    s.require("three categories (n = 4000): the true coefficients within four standard errors", within && m.converged,
              std::to_string(m.iterations) + " iterations");
    std::vector<double> share(3, 0.0), mean(3, 0.0);
    for (std::size_t i = 0; i < X.rows; i++) {
      share[y[i]] += 1.0 / static_cast<double>(X.rows);
      for (std::size_t j = 0; j < 3; j++) mean[j] += m.probabilities(i, j) / static_cast<double>(X.rows);
    }
    double dshare = 0;
    for (std::size_t j = 0; j < 3; j++) dshare = std::max(dshare, std::fabs(mean[j] - share[j]));
    s.check("... the mean predicted probabilities are the observed shares", dshare, 1e-8);
    const auto e = ne::multinomial_marginal_effects(m, X);
    double sum = 0;
    for (std::size_t l = 1; l < 3; l++) sum = std::max(sum, std::fabs(e.effect(0, l) + e.effect(1, l) + e.effect(2, l)));
    s.check("... the average marginal effects sum to zero over the categories", sum, 1e-14);
    // Central differences of the mean predicted probability in x₁.
    double worst = 0;
    const double h = 1e-5;
    auto Xp = X, Xm = X;
    for (std::size_t i = 0; i < X.rows; i++) Xp(i, 1) += h, Xm(i, 1) -= h;
    const auto Pp = ne::detail::choice_probabilities(Xp, m.beta, 3), Pm = ne::detail::choice_probabilities(Xm, m.beta, 3);
    for (std::size_t j = 0; j < 3; j++) {
      double d = 0;
      for (std::size_t i = 0; i < X.rows; i++) d += (Pp(i, j) - Pm(i, j)) / (2 * h) / static_cast<double>(X.rows);
      worst = std::max(worst, std::fabs(d - e.effect(j, 1)));
    }
    s.check("... and are the derivatives of the mean predicted probabilities", worst, 1e-8);
    s.require("... the intercept has none", std::isnan(e.effect(1, 0)));
    ne::MultinomialOptions oc;
    oc.covariance = ne::MultinomialCovariance::hc0;
    const auto r = ne::multinomial_logit(X, y, 3, oc);
    s.require("... HC0 standard errors: the same estimate, positive and different variances",
              std::fabs(r.beta[1] - m.beta[1]) < 1e-12 && r.se[1] > 0 && r.se[1] != m.se[1]);
  }
  return s.done();
}
