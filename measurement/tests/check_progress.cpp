// Checks of nomos/measurement/progress: the compound rate by hand; a quantity halved in five years toward zero (ρ = 0.875, moderate);
// one raised from 80 to 90 toward 100 (ρ ≈ 0.773); movement away (regression); the classes at their floors; the rate's classes for a
// quantity raised and one lowered; a baseline at the target; a target met; the baseline year chosen from a series; invariance under a
// change of unit; a path that reaches its target exactly giving ρ = 1; refusals.

#include <nomos/measurement/progress.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <limits>
#include <vector>

using namespace nomos::measurement;

namespace {
template <class F>
bool refused(F f) {
  try {
    f();
  } catch (const std::invalid_argument&) {
    return true;
  }
  return false;
}
}  // namespace

int main() {
  calx::verification::Suite s;
  const ProgressScheme un;
  s.check("the compound rate of 10 to 12.1 over two years is 0.1", std::fabs(compound_rate(10, 12.1, 2) - 0.1), 1e-12);

  {  // halved in five years, toward zero by 2030: x̂ = 50 · 0.5² = 12.5, ρ = (12.5 − 100) / (0 − 100) = 0.875
    const auto p = trend(2015, 100, 2020, 50, Desired::decrease, un, true, 0, 2030);
    s.check("halved in five years toward zero: x̂_2030 = 12.5", std::fabs(p.expected - 12.5), 1e-9);
    s.check("... ρ = 0.875", std::fabs(p.ratio - 0.875), 1e-12);
    s.require("... moderate progress, but acceleration needed", p.level == 1 && p.basis == "ratio" && p.trend == un.names[1]);
  }
  {  // 80 to 90 from 2015 to 2025 toward 100: x̂ = 90 (9/8)^{1/2}, ρ = (x̂ − 80) / 20
    const auto p = trend(2015, 80, 2025, 90, Desired::increase, un, true, 100, 2030);
    const double x = 90 * std::sqrt(9.0 / 8.0);
    s.check("raised from 80 to 90 toward 100: ρ = (90·√(9/8) − 80)/20", std::fabs(p.ratio - (x - 80) / 20), 1e-12);
    s.require("... moderate progress (ρ ≈ 0.773)", p.level == 1);
  }
  {
    const auto p = trend(2015, 30, 2022, 36, Desired::decrease, un, true, 10, 2030);
    s.require("a quantity to be lowered that rose is regression (ρ < 0)", p.ratio < 0 && p.level == 4 && p.trend == "regression");
  }
  s.require("the classes at their floors: 0.95 on track, below it moderate; 0.5 moderate; 0.1 marginal; 0 stagnation; below 0 regression",
               level_of_ratio(0.95) == 0 && level_of_ratio(0.9499) == 1 && level_of_ratio(0.5) == 1 && level_of_ratio(0.4999) == 2 &&
                   level_of_ratio(0.1) == 2 && level_of_ratio(0.0999) == 3 && level_of_ratio(0) == 3 && level_of_ratio(-1e-12) == 4);
  s.require("without a target, a quantity raised: 2% on track, 1.25% moderate, 0.5% marginal, 0 stagnation, a fall regression",
               level_of_rate(0.02, Desired::increase) == 0 && level_of_rate(0.0125, Desired::increase) == 1 &&
                   level_of_rate(0.005, Desired::increase) == 2 && level_of_rate(0, Desired::increase) == 3 && level_of_rate(-0.001, Desired::increase) == 4);
  s.require("... a quantity lowered: −2% on track, −1.25% moderate, a rise regression",
               level_of_rate(-0.02, Desired::decrease) == 0 && level_of_rate(-0.0125, Desired::decrease) == 1 && level_of_rate(0.001, Desired::decrease) == 4);
  {
    const auto p = trend(2015, 20, 2023, 20 * std::pow(1.03, 8), Desired::increase, un);
    s.require("a series without a target is classed by its rate (3% a year raised: on track)", !p.has_target && p.level == 0 && p.basis == "rate");
  }
  {
    const auto good = trend(2015, 70, 2022, 65, Desired::decrease, un, true, 70, 2030);
    const auto bad = trend(2015, 70, 2022, 75, Desired::decrease, un, true, 70, 2030);
    s.require("a baseline at the target: the best class if the latest is at or beyond it, else regression",
                 good.level == 0 && good.basis == "baseline at target" && bad.level == 4 && std::isnan(good.ratio));
  }
  {
    const auto p = trend(2015, 98, 2023, 100, Desired::increase, un, true, 99, 2030);
    s.require("a latest value beyond the target is the target met", p.met && p.level == 0 && p.basis == "target met");
    const auto q = trend(2015, 101, 2023, 99.5, Desired::increase, un, true, 99, 2030);
    s.require("... met even while it falls toward the target", q.met && q.level == 0);
    const auto back = trend(2015, 10, 2023, 80, Desired::decrease, un, true, 70, 2030);
    const auto up = trend(2015, 95, 2023, 85, Desired::increase, un, true, 90, 2030);
    s.require("a baseline beyond the target that fell back from it is regression, not progress", back.level == 4 && up.level == 4 && back.basis == "fallen back from the target");
  }
  {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const std::vector<double> years{2010, 2012, 2015, 2018, 2022}, v{5, 6, 7, nan, 9};
    const auto p = trend_of_series(years, v, Desired::increase, un);
    s.require("the baseline is the base year when observed, the latest the last observed", p.baseline_year == 2015 && p.latest_year == 2022 && p.baseline == 7);
    const std::vector<double> later{2017, 2019, 2023}, w{4, 5, 6};
    s.require("... else the earliest year after it", trend_of_series(later, w, Desired::increase, un).baseline_year == 2017);
    s.require("a series with nothing after the base year is refused",
                 refused([] { trend_of_series({2005, 2010, 2014}, {1, 2, 3}, Desired::increase); }));
    s.require("... and one with a single year from it", refused([] { trend_of_series({2010, 2016}, {1, 2}, Desired::increase); }));
  }
  {
    calx::random::Stream g(20261126, 19, 1, 2);
    double worst_scale = 0, worst_exact = 0;
    for (int k = 0; k < 300; k++) {
      const double b = 2015, t = 2016 + std::floor(9 * g.uniform()), xb = 1 + 99 * g.uniform(), xt = 1 + 99 * g.uniform(), target = 1 + 99 * g.uniform();
      if (std::fabs(target - xb) < 1e-6) continue;
      const Desired d = target > xb ? Desired::increase : Desired::decrease;
      const auto p = trend(b, xb, t, xt, d, un, true, target, 2030);
      const double c = 0.001 + 1000 * g.uniform();
      const auto q = trend(b, c * xb, t, c * xt, d, un, true, c * target, 2030);
      if (!p.met && std::fabs(target - xb) > 1) worst_scale = std::max(worst_scale, std::fabs(p.ratio - q.ratio) / std::max(1.0, std::fabs(p.ratio)) + std::fabs(p.rate - q.rate) + (p.level != q.level));
      // a path that reaches the target exactly in 2030 covers the whole distance: ρ = 1
      const double r = std::pow(target / xb, 1.0 / 15) - 1, on = xb * std::pow(1 + r, t - b);
      const auto e = trend(b, xb, t, on, d, un, true, target, 2030);
      if (!e.met) worst_exact = std::max(worst_exact, std::fabs(e.ratio - 1));
    }
    s.check("the rate, ρ (relatively) and the class do not change with the unit (300 random series)", worst_scale, 1e-9);
    s.check("a compound path that reaches its target in 2030 gives ρ = 1 (300 random series)", worst_exact, 1e-9);
  }
  s.require("a baseline at zero, or a negative value, is refused", refused([] { compound_rate(0, 1, 3); }) && refused([] { compound_rate(1, -1, 3); }));
  s.require("a target year before the latest is refused", refused([] { trend(2015, 1, 2024, 2, Desired::increase, {}, true, 3, 2020); }));
  ProgressScheme bad = un;
  bad.names.pop_back();
  s.require("a scheme with as many names as floors is refused", refused([&] { level_of_ratio(0.5, bad); }));
  return s.done();
}
