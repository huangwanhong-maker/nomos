// Checks of Poisson pseudo-maximum likelihood with two absorbed sets of fixed effects
// (poisson_fixed_effects.hpp).
//
// The independence model: with no regressors, the fitted table is the product of its margins,
// μ_rc = y_r· y_·c / y_·· (Bishop, Fienberg and Holland 1975, ch. 3), here with an all-zero row
// and column dropped (fitted 0, effects −∞). Against the dense dummy-variable Poisson regression of
// tyche's GLM on a 10 × 12 table with zeros: the same coefficients, deviance and fitted values, the
// same HC0 standard errors (the β block of the full sandwich), and clustered standard errors in
// the ratio of the two small-sample factors; the recovered effects rebuild the linear predictor.
// A regressor constant within rows is absorbed. A gravity-like table of 200 × 200 cells less the
// diagonal (39,800 cells, continuous outcomes with zeros): converged, and the true coefficients
// within four standard errors.

#include <nomos/econometrics/poisson_fixed_effects.hpp>

#include <tyche/regression/glm.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <chrono>
#include <cmath>
#include <string>
#include <vector>

namespace ne = nomos::econometrics;
namespace tr = tyche::regression;
using calx::linalg::Rectangular;
using calx::verification::exponential;
using calx::verification::relative;

namespace {
calx::random::Stream z(20261025, 0x70706d6cu, 0, 0);  // "ppml"
}

int main() {
  calx::verification::Suite s;

  // The independence model, with an all-zero row (3) and column (5).
  {
    const std::size_t R = 6, C = 8;
    std::vector<double> y;
    std::vector<std::size_t> row, col;
    for (std::size_t r = 0; r < R; r++)
      for (std::size_t c = 0; c < C; c++) {
        row.push_back(r);
        col.push_back(c);
        y.push_back(r == 3 || c == 5 ? 0.0 : (z.uniform() < 0.2 ? 0.0 : 10 * z.uniform()));
      }
    const auto fit = ne::poisson_two_way(Rectangular(), y, row, col);
    std::vector<double> rt(R, 0), ct(C, 0);
    double total = 0;
    for (std::size_t i = 0; i < y.size(); i++) rt[row[i]] += y[i], ct[col[i]] += y[i], total += y[i];
    double worst = 0;
    for (std::size_t i = 0; i < y.size(); i++)
      if (!fit.dropped[i]) worst = std::max(worst, relative(fit.fitted[i], rt[row[i]] * ct[col[i]] / total));
    s.check("the independence model (6 × 8, zeros): μ_rc = y_r· y_·c / y_·· (relative)", worst, 1e-9);
    s.require("... the all-zero row and column dropped, fitted 0, effects −∞",
              fit.dropped_rows == 1 && fit.dropped_columns == 1 && fit.used == 5 * 7 && fit.fitted[3 * C] == 0 && std::isinf(fit.row_effects[3]) && std::isinf(fit.column_effects[5]));
    s.require("... converged", fit.converged, std::to_string(fit.iterations) + " iterations");
  }
  // Against the dense dummy-variable fit.
  {
    const std::size_t R = 10, C = 12;
    std::vector<double> alpha(R), gamma(C);
    for (double& a : alpha) a = z.normal();
    for (double& g : gamma) g = z.normal();
    std::vector<double> y;
    std::vector<std::size_t> row, col;
    std::vector<std::int64_t> cluster;
    std::vector<std::vector<double>> xs;
    for (std::size_t r = 0; r < R; r++)
      for (std::size_t c = 0; c < C; c++) {
        const double x1 = z.normal(), x2 = z.uniform() < 0.3 ? 1 : 0;
        const double mu = std::exp(0.5 + alpha[r] + gamma[c] - 0.8 * x1 + 0.4 * x2);
        y.push_back(static_cast<double>(z.poisson(mu)));
        xs.push_back({x1, x2, static_cast<double>(r % 3)});
        row.push_back(r);
        col.push_back(c);
        cluster.push_back(static_cast<std::int64_t>(r));
      }
    std::size_t zeros = 0;
    for (double v : y) zeros += v == 0;
    Rectangular X(y.size(), 2), Xw(y.size(), 3), D(y.size(), 2 + R + C - 1);
    for (std::size_t i = 0; i < y.size(); i++) {
      for (std::size_t j = 0; j < 2; j++) X(i, j) = D(i, j) = xs[i][j];
      for (std::size_t j = 0; j < 3; j++) Xw(i, j) = xs[i][j];
      D(i, 2 + row[i]) = 1;
      if (col[i] > 0) D(i, 2 + R + col[i] - 1) = 1;
    }
    const auto fit = ne::poisson_two_way(X, y, row, col);
    tr::GlmOptions go;
    go.family = tr::Family::poisson;
    go.covariance = tr::GlmCovariance::hc0;
    go.tolerance = 1e-14;
    const auto dense = tr::glm(D, y, go);
    s.check("against dense dummies (10 × 12, " + std::to_string(zeros) + " zeros): coefficients", std::max(std::fabs(fit.beta[0] - dense.beta[0]), std::fabs(fit.beta[1] - dense.beta[1])), 1e-8);
    s.check("... deviance (relative)", relative(fit.deviance, dense.deviance), 1e-10);
    double df = 0, deta = 0;
    for (std::size_t i = 0; i < y.size(); i++) {
      df = std::max(df, relative(fit.fitted[i], dense.fitted[i]));
      const double eta = X(i, 0) * fit.beta[0] + X(i, 1) * fit.beta[1] + fit.row_effects[row[i]] + fit.column_effects[col[i]];
      deta = std::max(deta, std::fabs(eta - dense.linear_predictor[i]));
    }
    s.check("... fitted values (relative)", df, 1e-8);
    s.check("... and the recovered effects rebuild the linear predictor", deta, 1e-7);
    s.check("... HC0 standard errors: the β block of the full sandwich (relative)", std::max(relative(fit.se[0], dense.se[0]), relative(fit.se[1], dense.se[1])), 1e-6);
    ne::PoissonTwoWayOptions oc;
    oc.covariance = ne::PoissonCovariance::cluster;
    oc.clusters = cluster;
    const auto fc = ne::poisson_two_way(X, y, row, col, oc);
    go.covariance = tr::GlmCovariance::cluster;
    go.clusters = cluster;
    const auto densec = tr::glm(D, y, go);
    const double n = static_cast<double>(y.size()), r = static_cast<double>(densec.rank);
    s.check("... clustered by row: in the ratio √((n − r)/(n − 1)) of the small-sample factors (relative)",
            std::max(relative(fc.se[0] / densec.se[0], std::sqrt((n - r) / (n - 1))), relative(fc.se[1] / densec.se[1], std::sqrt((n - r) / (n - 1)))), 1e-6);
    const auto fw = ne::poisson_two_way(Xw, y, row, col);
    s.require("... a regressor constant within rows is absorbed by the row effects", fw.aliased[2] && std::isnan(fw.beta[2]) && std::fabs(fw.beta[0] - fit.beta[0]) < 1e-8);
  }
  // A gravity-like table: 200 × 200 less the diagonal.
  {
    const std::size_t N = 200;
    std::vector<double> lat(N), lon(N), a(N), g(N);
    for (std::size_t i = 0; i < N; i++) {
      lat[i] = z.uniform(-1, 1);
      lon[i] = z.uniform(-3, 3);
      a[i] = 2 + z.normal();
      g[i] = z.normal();
    }
    std::vector<double> y;
    std::vector<std::size_t> row, col;
    std::vector<std::vector<double>> xs;
    for (std::size_t i = 0; i < N; i++)
      for (std::size_t j = 0; j < N; j++) {
        if (i == j) continue;
        const double dist = std::sqrt((lat[i] - lat[j]) * (lat[i] - lat[j]) + (lon[i] - lon[j]) * (lon[i] - lon[j])) + 0.01;
        const double near = dist < 0.4 ? 1 : 0;
        const double mu = std::exp(a[i] + g[j] - 1.0 * std::log(dist) + 0.5 * near);
        // A continuous outcome with mean μ and zeros: μ times a Bernoulli(½) times a gamma of mean 2.
        y.push_back(z.uniform() < 0.5 ? 0.0 : mu * z.gamma(2.0));
        xs.push_back({std::log(dist), near});
        row.push_back(i);
        col.push_back(j);
      }
    const auto X = Rectangular::from_rows(xs);
    const auto start = std::chrono::steady_clock::now();
    const auto fit = ne::poisson_two_way(X, y, row, col);
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    s.require("gravity-like, 200 × 200 less the diagonal (39,800 cells, half zero): converged", fit.converged,
              std::to_string(fit.iterations) + " iterations, at most " + std::to_string(fit.sweeps) + " sweeps, " + exponential(seconds) + " s");
    s.require("... the true coefficients −1 and 0.5 within four standard errors",
              std::fabs(fit.beta[0] + 1) < 4 * fit.se[0] && std::fabs(fit.beta[1] - 0.5) < 4 * fit.se[1],
              exponential(fit.beta[0]) + " ± " + exponential(fit.se[0]) + ", " + exponential(fit.beta[1]) + " ± " + exponential(fit.se[1]));
  }
  return s.done();
}
