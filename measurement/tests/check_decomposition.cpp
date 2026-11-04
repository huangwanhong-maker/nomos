// Checks of the decompositions by groups (decomposition.hpp).
//
// On a random weighted population in four groups (one label unused): within + between = total
// for GE(α), α ∈ {−1, 0, ½, 1, 2, 3}; the weights of Theil's L are the population shares and
// those of Theil's T the shares of the total (Bourguignon 1979; Shorrocks 1980); the between part
// vanishes when every group has the same mean, and the within part when every group is equal
// inside; the contributions of the groups to FGT(α) add up to it.

#include <nomos/measurement/decomposition.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <vector>

namespace nm = nomos::measurement;
using calx::verification::relative;

int main() {
  calx::verification::Suite s;
  calx::random::Stream z(20261024, 0x6e6f6d6fu, 3, 0);

  std::vector<double> x(400), w(400);
  std::vector<std::size_t> g(400);
  for (std::size_t i = 0; i < x.size(); i++) {
    g[i] = std::vector<std::size_t>{0, 1, 3, 4}[z.below(4)];  // label 2 unused
    x[i] = (1 + static_cast<double>(g[i])) * std::exp(0.6 * z.normal());
    w[i] = 0.2 + z.uniform();
  }
  double worst = 0;
  for (double a : {-1.0, 0.0, 0.5, 1.0, 2.0, 3.0}) {
    const auto d = nm::decompose_generalised_entropy(x, g, a, w);
    worst = std::max(worst, relative(d.within + d.between, d.total));
    worst = std::max(worst, relative(d.total, nm::generalised_entropy(x, a, w)));
  }
  s.check("four groups (one empty), α ∈ {−1, 0, ½, 1, 2, 3}: within + between = total (relative)", worst, 1e-12);
  const auto L = nm::decompose_generalised_entropy(x, g, 0, w);
  const auto T = nm::decompose_generalised_entropy(x, g, 1, w);
  double dl = 0, dt = 0, sum_v = 0, sum_s = 0;
  for (std::size_t k = 0; k < L.weight.size(); k++) {
    dl = std::max(dl, std::fabs(L.weight[k] - L.population_share[k]));
    dt = std::max(dt, std::fabs(T.weight[k] - T.amount_share[k]));
    sum_v += L.population_share[k];
    sum_s += T.amount_share[k];
  }
  s.check("Theil's L weighs the groups by their population shares", dl, 1e-15);
  s.check("Theil's T weighs them by their shares of the total", dt, 1e-15);
  s.check("... both sets of shares sum to one", std::fabs(sum_v - 1) + std::fabs(sum_s - 1), 1e-14);
  s.require("the unused label is an empty group of weight zero", L.population_share[2] == 0 && L.weight[2] == 0 && L.weight.size() == 5);

  // Every group the same distribution: no inequality between them.
  {
    std::vector<double> y, wy;
    std::vector<std::size_t> gy;
    for (std::size_t k = 0; k < 3; k++)
      for (double v : {1.0, 2.0, 4.0, 8.0}) {
        y.push_back(v);
        wy.push_back(1);
        gy.push_back(k);
      }
    const auto d = nm::decompose_generalised_entropy(y, gy, 2, wy);
    s.check("groups with one distribution: the between part vanishes", d.between, 1e-15);
    s.check("... and the within part is the whole", relative(d.within, d.total), 1e-14);
  }
  // Every group equal inside: all inequality is between.
  {
    const std::vector<double> y{1, 1, 1, 5, 5, 9};
    const std::vector<std::size_t> gy{0, 0, 0, 1, 1, 2};
    const auto d = nm::decompose_generalised_entropy(y, gy, 0.5);
    s.check("groups equal inside: the within part vanishes", d.within, 1e-15);
    s.check("... and the between part is the whole", relative(d.between, d.total), 1e-14);
  }
  // Poverty: contributions add up.
  {
    double worst_p = 0;
    for (double a : {0.0, 1.0, 2.0}) {
      const auto d = nm::decompose_fgt(x, g, 2.0, a, w);
      double sum = 0;
      for (double c : d.contribution) sum += c;
      worst_p = std::max({worst_p, std::fabs(sum - d.total), std::fabs(d.total - nm::fgt(x, 2.0, a, w))});
    }
    s.check("FGT(0, 1, 2) at the line 2: the contributions of the groups add up", worst_p, 1e-15);
  }
  return s.done();
}
