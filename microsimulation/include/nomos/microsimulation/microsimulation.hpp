// Static microsimulation of taxes and benefits: a declared policy applied to every household of a weighted sample, and a
// reform compared with the baseline in revenue, in incomes, in inequality and poverty, and in the incentive to earn.
//
// A tax-benefit microsimulation model computes, for each unit of a representative sample, what a set of rules makes of its
// gross incomes: taxes and contributions paid, benefits received, disposable income (Sutherland and Figari 2013, EUROMOD:
// the European Union tax-benefit microsimulation model, International Journal of Microsimulation 6, 4–26; Bourguignon and
// Spadaro 2006, Microsimulation as a tool for evaluating redistribution policies, Journal of Economic Inequality 4, 77–106).
// It is static: incomes and behaviour are held as they are, so that what changes is the rules alone (the "morning after"
// effect of a reform). Behavioural responses, and the ageing of a population, are dynamic models' and are not here.
//
// The rules (Policy) are the caller's, and no jurisdiction's:
//   - an income tax on each person's taxable income (earnings and other income, less a personal allowance), by marginal
//     rates above thresholds;
//   - a social contribution on earnings at a flat rate up to a ceiling;
//   - universal amounts: a basic income per adult and an allowance per child;
//   - a means-tested benefit per household: a maximum by its adults and children, withdrawn at a taper rate on the
//     household's income after tax and contributions above a disregard, never below zero.
// Equivalised income divides disposable income by an equivalence scale, 1 for the first adult and stated weights for the
// other adults and for children (the modified OECD scale is 1, 0.5, 0.3: Hagenaars, de Vos and Zaidi 1994). Inequality and
// poverty are nomos's measures, over persons: each household's equivalised income counts once for each of its members,
// weighted by its survey weight. The poverty line is a share of the median, or a stated amount, or the baseline's line held
// fixed when a reform is compared (anchored).
//
// Incentives: the marginal effective tax rate of a person is 1 − Δdisposable/Δearnings for a small rise in their earnings,
// with every rule of the household applied; the participation tax rate is 1 − (disposable working − disposable not
// working)/earnings (Immervoll, Kleven, Kreiner and Saez 2007, Welfare reform in European countries, Economic Journal 117,
// 1–44). Both are computed by applying the rules again, not by formula, so that interactions among rules are included.
//
// Units: amounts per year in one currency; rates as fractions.

#pragma once

#include <nomos/measurement/inequality.hpp>
#include <nomos/measurement/poverty.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nomos::microsimulation {

struct Person {
  double earnings = 0, other = 0;  // gross, per year
  bool child = false;
};

struct Household {
  std::string id;
  double weight = 1;  // the survey weight: how many households of the population it stands for
  std::vector<Person> persons;
};

struct Bracket {
  double from = 0;  // the threshold of taxable income above which the rate applies
  double rate = 0;
};

struct Schedule {
  double allowance = 0;            // taxable income = max(0, income − allowance)
  std::vector<Bracket> brackets;   // thresholds strictly increasing, the first at 0
};

struct Benefit {
  double base = 0, per_adult = 0, per_child = 0;  // the maximum: base + per adult + per child
  double disregard = 0, taper = 0;                 // withdrawn at `taper` on income after tax above the disregard
};

struct Policy {
  std::string name;
  Schedule income_tax;
  double contribution_rate = 0, contribution_ceiling = std::numeric_limits<double>::infinity();
  double basic_income = 0, child_allowance = 0;  // per adult, per child
  Benefit benefit;
};

struct Scale {
  double first = 1, adult = 0.5, child = 0.3;  // the modified OECD scale by default
};

struct Outcome {
  double gross = 0, tax = 0, contributions = 0, universal = 0, benefit = 0, disposable = 0;
  double equivalised = 0;
  std::size_t persons = 0;
};

namespace detail {
inline void validate(const Schedule& s) {
  if (s.brackets.empty() || s.brackets[0].from != 0) throw std::invalid_argument("nomos::microsimulation: a schedule's first bracket starts at 0");
  for (std::size_t k = 1; k < s.brackets.size(); k++)
    if (!(s.brackets[k].from > s.brackets[k - 1].from)) throw std::invalid_argument("nomos::microsimulation: a schedule's thresholds strictly increase");
  for (const auto& b : s.brackets)
    if (!(b.rate >= 0 && b.rate <= 1)) throw std::invalid_argument("nomos::microsimulation: a marginal rate lies in [0, 1]");
  if (!(s.allowance >= 0)) throw std::invalid_argument("nomos::microsimulation: a negative allowance");
}
}  // namespace detail

// The tax on an income under a schedule.
inline double tax(const Schedule& s, double income) {
  detail::validate(s);
  const double taxable = std::max(0.0, income - s.allowance);
  double t = 0;
  for (std::size_t k = 0; k < s.brackets.size(); k++) {
    const double lo = s.brackets[k].from, hi = k + 1 < s.brackets.size() ? s.brackets[k + 1].from : std::numeric_limits<double>::infinity();
    if (taxable > lo) t += (std::min(taxable, hi) - lo) * s.brackets[k].rate;
  }
  return t;
}

// The household's outcome under a policy.
inline Outcome apply(const Policy& p, const Household& h, const Scale& scale = {}) {
  if (h.persons.empty()) throw std::invalid_argument("nomos::microsimulation: the household " + h.id + " has no members");
  if (!(p.contribution_rate >= 0 && p.contribution_rate <= 1) || !(p.benefit.taper >= 0 && p.benefit.taper <= 1))
    throw std::invalid_argument("nomos::microsimulation: rates of contribution and withdrawal lie in [0, 1]");
  Outcome o;
  std::size_t adults = 0, children = 0;
  for (const auto& m : h.persons) {
    if (m.earnings < 0 || m.other < 0) throw std::invalid_argument("nomos::microsimulation: a negative income in " + h.id);
    (m.child ? children : adults)++;
    o.gross += m.earnings + m.other;
    o.tax += tax(p.income_tax, m.earnings + m.other);
    o.contributions += p.contribution_rate * std::min(m.earnings, p.contribution_ceiling);
  }
  o.universal = p.basic_income * static_cast<double>(adults) + p.child_allowance * static_cast<double>(children);
  const double net = o.gross - o.tax - o.contributions;
  const double maximum = p.benefit.base + p.benefit.per_adult * static_cast<double>(adults) + p.benefit.per_child * static_cast<double>(children);
  o.benefit = std::max(0.0, maximum - p.benefit.taper * std::max(0.0, net - p.benefit.disregard));
  o.disposable = net + o.universal + o.benefit;
  o.persons = h.persons.size();
  const double size = adults == 0 ? scale.first + scale.child * (static_cast<double>(children) - 1)
                                  : scale.first + scale.adult * (static_cast<double>(adults) - 1) + scale.child * static_cast<double>(children);
  o.equivalised = o.disposable / size;
  return o;
}

struct PovertyLine {
  enum class Kind { share_of_median, amount } kind = Kind::share_of_median;
  double value = 0.6;
};

struct Result {
  std::string policy;
  std::vector<Outcome> households;
  double gross = 0, taxes = 0, contributions = 0, universal = 0, benefits = 0, disposable = 0;  // weighted totals
  double net_revenue = 0;  // taxes + contributions − universal − benefits
  double gini = 0, poverty_line = 0, poverty_rate = 0, median = 0;  // over persons
};

namespace detail {
// Each household's equivalised income with the weight of its members.
inline std::pair<std::vector<double>, std::vector<double>> persons(const std::vector<Household>& pop, const std::vector<Outcome>& out) {
  std::vector<double> x, w;
  for (std::size_t i = 0; i < pop.size(); i++) {
    x.push_back(out[i].equivalised);
    w.push_back(pop[i].weight * static_cast<double>(out[i].persons));
  }
  return {x, w};
}
inline double weighted_median(std::vector<double> x, std::vector<double> w) {
  std::vector<std::size_t> idx(x.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) { return x[a] < x[b]; });
  const double total = std::accumulate(w.begin(), w.end(), 0.0);
  double run = 0;
  for (const auto i : idx) {
    run += w[i];
    if (run >= total / 2) return x[i];
  }
  return x.empty() ? 0 : x[idx.back()];
}
}  // namespace detail

inline Result simulate(const Policy& p, const std::vector<Household>& pop, const PovertyLine& line = {}, const Scale& scale = {}) {
  if (pop.empty()) throw std::invalid_argument("nomos::microsimulation::simulate: an empty population");
  Result r;
  r.policy = p.name;
  for (const auto& h : pop) {
    if (!(h.weight > 0)) throw std::invalid_argument("nomos::microsimulation::simulate: the household " + h.id + " has no positive weight");
    const auto o = apply(p, h, scale);
    r.households.push_back(o);
    r.gross += h.weight * o.gross;
    r.taxes += h.weight * o.tax;
    r.contributions += h.weight * o.contributions;
    r.universal += h.weight * o.universal;
    r.benefits += h.weight * o.benefit;
    r.disposable += h.weight * o.disposable;
  }
  r.net_revenue = r.taxes + r.contributions - r.universal - r.benefits;
  const auto [x, w] = detail::persons(pop, r.households);
  r.gini = nomos::measurement::gini(x, w);
  r.median = detail::weighted_median(x, w);
  r.poverty_line = line.kind == PovertyLine::Kind::amount ? line.value : line.value * r.median;
  r.poverty_rate = nomos::measurement::headcount_ratio(x, r.poverty_line, w);
  return r;
}

struct Comparison {
  double revenue_change = 0;               // the reform's net revenue less the baseline's
  double gini_before = 0, gini_after = 0;
  double poverty_before = 0, poverty_after = 0;  // the reform measured against the baseline's line (anchored)
  std::vector<double> decile_change;       // mean change of equivalised income by decile of the baseline's, over persons
  double gainers = 0, losers = 0, unchanged = 0;  // shares of persons, by a change beyond `tolerance`
};

inline Comparison compare(const Policy& baseline, const Policy& reform, const std::vector<Household>& pop, const PovertyLine& line = {}, const Scale& scale = {},
                          double tolerance = 1) {
  const auto b = simulate(baseline, pop, line, scale);
  const auto r = simulate(reform, pop, PovertyLine{PovertyLine::Kind::amount, b.poverty_line}, scale);
  Comparison c;
  c.revenue_change = r.net_revenue - b.net_revenue;
  c.gini_before = b.gini;
  c.gini_after = r.gini;
  c.poverty_before = b.poverty_rate;
  c.poverty_after = r.poverty_rate;
  const auto [x, w] = detail::persons(pop, b.households);
  std::vector<std::size_t> idx(x.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t bb) { return x[a] < x[bb]; });
  const double total = std::accumulate(w.begin(), w.end(), 0.0);
  std::vector<double> sum(10, 0), weight(10, 0);
  double run = 0;
  for (const auto i : idx) {
    // A household belongs to the decile of the midpoint of its persons' weight in the order of the baseline's incomes.
    const double mid = (run + w[i] / 2) / total;
    run += w[i];
    const std::size_t d = std::min<std::size_t>(9, static_cast<std::size_t>(mid * 10));
    const double change = r.households[i].equivalised - b.households[i].equivalised;
    sum[d] += w[i] * change;
    weight[d] += w[i];
    if (change > tolerance) c.gainers += w[i];
    else if (change < -tolerance) c.losers += w[i];
    else c.unchanged += w[i];
  }
  for (std::size_t d = 0; d < 10; d++) c.decile_change.push_back(weight[d] > 0 ? sum[d] / weight[d] : std::numeric_limits<double>::quiet_NaN());
  c.gainers /= total, c.losers /= total, c.unchanged /= total;
  return c;
}

// The marginal effective tax rate of a member of a household, for a rise of `delta` in their earnings.
inline double marginal_effective_rate(const Policy& p, const Household& h, std::size_t member, double delta = 1) {
  if (member >= h.persons.size()) throw std::invalid_argument("nomos::microsimulation: no such member");
  if (!(delta > 0)) throw std::invalid_argument("nomos::microsimulation: a positive rise in earnings");
  Household up = h;
  up.persons[member].earnings += delta;
  return 1 - (apply(p, up).disposable - apply(p, h).disposable) / delta;
}

// The participation tax rate of a member: the share of their earnings lost to taxes and withdrawn benefits when they work.
inline double participation_rate(const Policy& p, const Household& h, std::size_t member) {
  if (member >= h.persons.size()) throw std::invalid_argument("nomos::microsimulation: no such member");
  if (!(h.persons[member].earnings > 0)) throw std::invalid_argument("nomos::microsimulation: the member earns nothing");
  Household idle = h;
  idle.persons[member].earnings = 0;
  return 1 - (apply(p, h).disposable - apply(p, idle).disposable) / h.persons[member].earnings;
}

// Incomes uprated: earnings and other income multiplied by their factors.
inline std::vector<Household> uprate(std::vector<Household> pop, double earnings, double other) {
  if (!(earnings > 0) || !(other > 0)) throw std::invalid_argument("nomos::microsimulation::uprate: positive factors");
  for (auto& h : pop)
    for (auto& m : h.persons) m.earnings *= earnings, m.other *= other;
  return pop;
}

}  // namespace nomos::microsimulation
