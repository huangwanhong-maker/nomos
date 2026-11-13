// Cost–benefit analysis: streams of benefits and costs over years brought to the present, with rates declared, weights for whom
// they fall on, and the values at which a decision would change (decision 0012; Boardman et al. 2018, Cost–Benefit Analysis:
// Concepts and Practice, 5th ed., ch. 9–10; HM Treasury 2022, The Green Book).
//
//   Schedule                     a discount rate year by year: constant, or declining by stated bands (as the Green Book's 3.5% for
//                                years 0–30, 3.0% for 31–75, … ; HM Treasury 2022), each year's factor the product of
//                                1/(1 + r_t) over the years before it, so that a later band never revalues the earlier years
//   appraise(b, c, schedule)     present values of benefits and costs, the net present value, the benefit–cost ratio, the discounted
//                                payback year, the equivalent annual value; the internal rate of return where the net stream changes
//                                sign once (Brent's method on NPV(r) = 0); the switching values: the factor on benefits, and on costs,
//                                at which the net present value would be zero
//   distributional(groups, η)    weights for whom benefits fall on: w_g = (ȳ / y_g)^η, the elasticity of the marginal utility of
//                                income η declared (the Green Book takes about 1.3), so that a pound to the poorer
//                                counts for more; the weighted net present value
//
// Conventions: year 0 is undiscounted; flows fall at the ends of their years; amounts in one currency at constant prices.

#pragma once

#include <calx/solvers/roots.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace nomos::appraisal {

struct Schedule {
  std::vector<std::size_t> from;  // the first year of each band (from[0] = 0), increasing
  std::vector<double> rate;       // its annual rate
  static Schedule constant(double r) { return {{0}, {r}}; }
  double rate_in(std::size_t year) const {
    double r = rate.at(0);
    for (std::size_t b = 0; b < from.size(); b++)
      if (year >= from[b]) r = rate[b];
    return r;
  }
  // Factors for years 0 .. T − 1: D(0) = 1, D(t) = D(t − 1) / (1 + r_t) with r_t the rate of year t.
  std::vector<double> factors(std::size_t T) const {
    std::vector<double> d(T, 1.0);
    for (std::size_t t = 1; t < T; t++) d[t] = d[t - 1] / (1 + rate_in(t));
    return d;
  }
};

struct Appraisal {
  double benefits = 0, costs = 0, npv = 0;  // present values
  double ratio = NAN;                       // benefits / costs
  double irr = NAN;                         // where the net stream changes sign once
  double payback = NAN;                     // the first year the cumulative discounted net value is ≥ 0
  double annual = NAN;                      // the equivalent annual value at the first band's rate
  double benefit_switch = NAN, cost_switch = NAN;  // the factors on benefits and on costs that make NPV zero
};

inline Appraisal appraise(const std::vector<double>& benefit, const std::vector<double>& cost, const Schedule& s) {
  if (benefit.size() != cost.size() || benefit.empty()) throw std::invalid_argument("nomos::appraisal::appraise: one benefit and one cost a year");
  const std::size_t T = benefit.size();
  const auto d = s.factors(T);
  Appraisal a;
  double cum = 0;
  for (std::size_t t = 0; t < T; t++) {
    a.benefits += benefit[t] * d[t], a.costs += cost[t] * d[t];
    cum += (benefit[t] - cost[t]) * d[t];
    if (std::isnan(a.payback) && cum >= 0 && t > 0) a.payback = static_cast<double>(t);
  }
  a.npv = a.benefits - a.costs;
  if (a.costs > 0) a.ratio = a.benefits / a.costs;
  if (a.benefits > 0) a.benefit_switch = a.costs / a.benefits;
  if (a.costs > 0) a.cost_switch = a.benefits / a.costs;
  const double r0 = s.rate.at(0);
  if (r0 > 0) {
    double annuity = 0;  // Σ_{t=1}^{T−1} (1 + r)^{−t}
    for (std::size_t t = 1; t < T; t++) annuity += std::pow(1 + r0, -static_cast<double>(t));
    if (annuity > 0) a.annual = a.npv / annuity;
  }
  // IRR: the net stream must change sign exactly once (Descartes's rule makes the root unique).
  std::size_t changes = 0;
  double last = 0;
  for (std::size_t t = 0; t < T; t++) {
    const double n = benefit[t] - cost[t];
    if (n == 0) continue;
    if (last != 0 && (n > 0) != (last > 0)) changes++;
    last = n;
  }
  if (changes == 1) {
    auto npv = [&](double r) {
      double v = 0, f = 1;
      for (std::size_t t = 0; t < T; t++, f /= 1 + r) v += (benefit[t] - cost[t]) * f;
      return v;
    };
    double lo = -0.99, hi = 1.0;
    for (int i = 0; i < 60 && npv(lo) * npv(hi) > 0; i++) hi *= 2;
    if (npv(lo) * npv(hi) <= 0) a.irr = calx::solvers::brent(npv, lo, hi, 1e-14, 300).x;
  }
  return a;
}

struct Distributional {
  std::vector<double> weight;  // per group
  double weighted_npv = 0;
};

// Each group's mean income y_g and the present value of its net benefit; the mean income ȳ of the population (by its shares).
inline Distributional distributional(const std::vector<double>& income, const std::vector<double>& net_present_value, double eta, double mean_income) {
  if (income.size() != net_present_value.size()) throw std::invalid_argument("nomos::appraisal::distributional: one income and one value per group");
  Distributional d;
  for (std::size_t g = 0; g < income.size(); g++) {
    if (!(income[g] > 0) || !(mean_income > 0)) throw std::invalid_argument("nomos::appraisal::distributional: incomes > 0");
    d.weight.push_back(std::pow(mean_income / income[g], eta));
    d.weighted_npv += d.weight.back() * net_present_value[g];
  }
  return d;
}

}  // namespace nomos::appraisal
