// Checks of the indices of inequality (inequality.hpp).
//
// Closed forms: equal amounts give zero; one unit holding everything among n gives a Gini
// coefficient of (n − 1)/n, Theil's T of ln n and GE(2) of (n − 1)/2; the amounts 1, …, 5 have
// mean difference 2 and a Gini coefficient of 4/15; ten amounts 1, …, 10 have a Palma ratio of 1.
// The exponential distribution, sampled at its midpoint quantiles, against its indices (Gini ½,
// Theil's T 1 − γ, Theil's L γ, GE(2) ½, Atkinson's A(½) 1 − π/4 and A(1) 1 − e^{−γ}, the Lorenz
// curve p + (1 − p) ln(1 − p), from E[x ln x] = Γ′(2) = 1 − γ and E[√x] = Γ(3/2)). Identities: the
// O(n log n) Gini against the double sum of the mean difference, and against 1 − 2∫L; Atkinson
// against GE(1 − ε); CV² = 2 GE(2); GE continuous at α = 0 and 1. Invariances: scale, replication,
// integer weights as replication, order of tied amounts. The principle of transfers.

#include <nomos/measurement/inequality.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <numbers>
#include <string>
#include <vector>

namespace nm = nomos::measurement;
using calx::verification::relative;

namespace {

double mean_difference_gini(const std::vector<double>& x, const std::vector<double>& w) {
  double W = 0, S = 0, D = 0;
  for (std::size_t i = 0; i < x.size(); i++) {
    W += w[i];
    S += w[i] * x[i];
    for (std::size_t j = 0; j < x.size(); j++) D += w[i] * w[j] * std::fabs(x[i] - x[j]);
  }
  return D / (2 * W * W * (S / W));
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream z(20261024, 0x6e6f6d6fu, 1, 0);  // "nomo"
  const std::vector<double> alphas{-1, 0, 0.5, 1, 2, 3};

  // Equal amounts: every index vanishes.
  {
    const std::vector<double> x(17, 3.25);
    double worst = std::fabs(nm::gini(x));
    for (double a : alphas) worst = std::max(worst, std::fabs(nm::generalised_entropy(x, a)));
    for (double e : {0.5, 1.0, 2.0}) worst = std::max(worst, std::fabs(nm::atkinson(x, e)));
    worst = std::max(worst, nm::coefficient_of_variation(x));
    s.check("equal amounts: Gini, GE(α) for α ∈ {−1, 0, ½, 1, 2, 3}, A(ε) and CV all vanish", worst, 1e-15);
    const auto c = nm::lorenz(x);
    s.check("... and the Lorenz curve is the diagonal: L(0.37) = 0.37", c.at(0.37) - 0.37, 1e-15);
  }
  // One unit holds everything among n = 10.
  {
    std::vector<double> x(10, 0.0);
    x[6] = 42;
    s.check("one holds all among 10: Gini = (n − 1)/n = 0.9", nm::gini(x) - 0.9, 1e-15);
    s.check("... the Lorenz area gives the same", nm::gini(nm::lorenz(x)) - 0.9, 1e-15);
    s.check("... Theil's T = ln n", nm::theil_t(x) - std::log(10.0), 1e-15);
    s.check("... GE(2) = (n − 1)/2", nm::generalised_entropy(x, 2) - 4.5, 1e-13);
    s.check("... the richest tenth holds everything", nm::top_share(nm::lorenz(x), 0.1) - 1, 1e-15);
    bool threw = false;
    try {
      (void)nm::theil_l(x);
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    s.require("... and Theil's L, which needs positive amounts, refuses zeros", threw);
  }
  // 1, …, 5: Σ_ij |x_i − x_j| = 40, μ = 3, so G = 40/(2·25·3) = 4/15.
  {
    const std::vector<double> x{3, 1, 5, 2, 4};
    s.check("1, …, 5: Gini = 4/15", nm::gini(x) - 4.0 / 15, 1e-15);
    const auto c = nm::lorenz(x);
    s.check("... L(0.4) = 3/15", c.at(0.4) - 0.2, 1e-15);
    s.check("... L(0.5) = 4.5/15, interpolated", c.at(0.5) - 0.3, 1e-15);
    s.check("... generalised Lorenz at 1 is the mean", c.generalised_at(1) - 3, 1e-15);
  }
  // 1, …, 10: the poorest 40 per cent hold 10/55, the richest 10 per cent hold 10/55.
  {
    std::vector<double> x;
    for (int i = 10; i >= 1; i--) x.push_back(i);
    s.check("1, …, 10: Palma ratio = (10/55)/(10/55) = 1", nm::palma_ratio(x) - 1, 1e-14);
    s.check("... the middle half (0.25 to 0.75) holds 27.5/55, by interpolation",
            nm::quantile_share(nm::lorenz(x), 0.25, 0.75) - (0.5 * 3 + 4 + 5 + 6 + 7 + 0.5 * 8) / 55, 1e-15);
  }
  // The Gini coefficient against the double sum, with random weights and ties.
  {
    std::vector<double> x(300), w(300);
    for (std::size_t i = 0; i < x.size(); i++) {
      x[i] = std::floor(20 * z.uniform());  // many ties
      w[i] = 0.1 + 3 * z.uniform();
    }
    const double g = nm::gini(x, w), ref = mean_difference_gini(x, w);
    s.check("300 weighted amounts with ties: Gini by covariance against the double sum (relative)", relative(g, ref), 1e-13);
    s.check("... and against 1 − 2∫L", relative(nm::gini(nm::lorenz(x, w)), ref), 1e-13);
    // Reversing the order (so that tied amounts come in the other order) changes nothing.
    std::vector<double> xr(x.rbegin(), x.rend()), wr(w.rbegin(), w.rend());
    s.check("... and the order of tied amounts does not matter", relative(nm::gini(xr, wr), g), 1e-14);
  }
  // Integer weights are replication; replication and scale leave every index unchanged.
  {
    std::vector<double> x, w, rep, big, scaled;
    for (int i = 0; i < 40; i++) {
      x.push_back(0.5 + 10 * z.uniform() * z.uniform());
      w.push_back(static_cast<double>(1 + z.below(4)));
      for (int k = 0; k < static_cast<int>(w.back()); k++) rep.push_back(x.back());
    }
    for (int k = 0; k < 3; k++) big.insert(big.end(), x.begin(), x.end());
    for (double v : x) scaled.push_back(7.3 * v);
    double worst_weight = relative(nm::gini(x, w), nm::gini(rep)), worst_rep = relative(nm::gini(big), nm::gini(x)),
           worst_scale = relative(nm::gini(scaled), nm::gini(x));
    for (double a : alphas) {
      worst_weight = std::max(worst_weight, relative(nm::generalised_entropy(x, a, w), nm::generalised_entropy(rep, a)));
      worst_rep = std::max(worst_rep, relative(nm::generalised_entropy(big, a), nm::generalised_entropy(x, a)));
      worst_scale = std::max(worst_scale, relative(nm::generalised_entropy(scaled, a), nm::generalised_entropy(x, a)));
    }
    for (double e : {0.5, 1.0, 2.0}) {
      worst_weight = std::max(worst_weight, relative(nm::atkinson(x, e, w), nm::atkinson(rep, e)));
      worst_rep = std::max(worst_rep, relative(nm::atkinson(big, e), nm::atkinson(x, e)));
      worst_scale = std::max(worst_scale, relative(nm::atkinson(scaled, e), nm::atkinson(x, e)));
    }
    worst_weight = std::max(worst_weight, relative(nm::coefficient_of_variation(x, w), nm::coefficient_of_variation(rep)));
    worst_weight = std::max(worst_weight, relative(nm::palma_ratio(x, w), nm::palma_ratio(rep)));
    worst_rep = std::max(worst_rep, relative(nm::palma_ratio(big), nm::palma_ratio(x)));
    s.check("integer weights are replication: Gini, GE(α), A(ε), CV, Palma (relative)", worst_weight, 1e-12);
    s.check("the principle of population: three copies of the population", worst_rep, 1e-12);
    s.check("scale invariance: every amount times 7.3", worst_scale, 1e-12);
    std::vector<double> w0 = w, x0 = x;
    x0.push_back(1e6);
    w0.push_back(0);
    s.check("a unit of zero weight counts for nothing", relative(nm::gini(x0, w0), nm::gini(x, w)), 1e-14);
  }
  // Identities among the indices, on lognormal-like amounts.
  {
    std::vector<double> x(500), w(500);
    for (std::size_t i = 0; i < x.size(); i++) {
      x[i] = std::exp(0.8 * z.normal());
      w[i] = 0.5 + z.uniform();
    }
    double worst = 0;
    for (double e : {0.25, 0.5, 2.0, 3.0}) {
      const double ge = nm::generalised_entropy(x, 1 - e, w);
      const double a = 1 - std::pow(1 + e * (e - 1) * ge, 1 / (1 - e));
      double m = 0, W = 0, mu = nm::mean(x, w);
      for (std::size_t i = 0; i < x.size(); i++) {
        m += w[i] * std::pow(x[i] / mu, 1 - e);
        W += w[i];
      }
      const double direct = 1 - std::pow(m / W, 1 / (1 - e));
      worst = std::max({worst, relative(nm::atkinson(x, e, w), a), relative(nm::atkinson(x, e, w), direct)});
    }
    s.check("Atkinson's A(ε) as the power mean and through GE(1 − ε), ε ∈ {¼, ½, 2, 3} (relative)", worst, 1e-12);
    double gm = 0, W = 0;
    for (std::size_t i = 0; i < x.size(); i++) {
      gm += w[i] * std::log(x[i]);
      W += w[i];
    }
    s.check("A(1) = 1 − geometric mean / mean", relative(nm::atkinson(x, 1, w), 1 - std::exp(gm / W) / nm::mean(x, w)), 1e-12);
    const double cv = nm::coefficient_of_variation(x, w);
    s.check("CV² = 2 GE(2)", relative(cv * cv, 2 * nm::generalised_entropy(x, 2, w)), 1e-12);
    s.check("GE is continuous at α = 0: GE(10⁻⁹) against Theil's L", relative(nm::generalised_entropy(x, 1e-9, w), nm::theil_l(x, w)), 1e-8);
    s.check("GE is continuous at α = 1: GE(1 + 10⁻⁹) against Theil's T", relative(nm::generalised_entropy(x, 1 + 1e-9, w), nm::theil_t(x, w)), 1e-8);
    s.check("... and from below, GE(1 − 10⁻⁹)", relative(nm::generalised_entropy(x, 1 - 1e-9, w), nm::theil_t(x, w)), 1e-8);
    // The principle of transfers: move 0.1 from the richest to the poorest.
    std::vector<double> y = x;
    std::size_t lo = 0, hi = 0;
    for (std::size_t i = 0; i < y.size(); i++) {
      if (y[i] < y[lo]) lo = i;
      if (y[i] > y[hi]) hi = i;
    }
    y[lo] += 0.1 / w[lo];
    y[hi] -= 0.1 / w[hi];
    bool lower = nm::gini(y, w) < nm::gini(x, w) && nm::coefficient_of_variation(y, w) < cv;
    for (double a : alphas) lower = lower && nm::generalised_entropy(y, a, w) < nm::generalised_entropy(x, a, w);
    for (double e : {0.5, 1.0, 2.0}) lower = lower && nm::atkinson(y, e, w) < nm::atkinson(x, e, w);
    s.require("the principle of transfers: a transfer from the richest to the poorest lowers every index", lower);
  }
  // The exponential distribution at its midpoint quantiles x_i = −ln(1 − (i − ½)/n).
  {
    const std::size_t n = 1000000;
    std::vector<double> x(n);
    for (std::size_t i = 0; i < n; i++) x[i] = -std::log1p(-(static_cast<double>(i) + 0.5) / static_cast<double>(n));
    const double gamma = 0.57721566490153286;
    s.check("the exponential distribution (10⁶ quantiles): Gini ½", nm::gini(x) - 0.5, 1e-6);
    s.check("... Theil's T = 1 − γ", nm::theil_t(x) - (1 - gamma), 5e-6);
    s.check("... Theil's L = γ", nm::theil_l(x) - gamma, 5e-6);
    s.check("... GE(2) = ½", nm::generalised_entropy(x, 2) - 0.5, 2e-5);
    s.check("... A(½) = 1 − Γ(3/2)² = 1 − π/4", nm::atkinson(x, 0.5) - (1 - std::numbers::pi / 4), 1e-6);
    s.check("... A(1) = 1 − e^{−γ}", nm::atkinson(x, 1) - (1 - std::exp(-gamma)), 2e-6);
    const auto c = nm::lorenz(x);
    double worst = 0;
    for (double p : {0.1, 0.4, 0.5, 0.9, 0.99}) worst = std::max(worst, std::fabs(c.at(p) - (p + (1 - p) * std::log1p(-p))));
    s.check("... the Lorenz curve p + (1 − p) ln(1 − p) at 0.1, 0.4, 0.5, 0.9, 0.99", worst, 2e-6);
  }
  return s.done();
}
