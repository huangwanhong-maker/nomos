// Checks of tax-benefit microsimulation (microsimulation.hpp): a schedule by hand; a household's outcome under every rule;
// marginal effective tax rates equal to the rates the rules combine to (the tax rate and contribution in work, the taper
// applied to what is left inside the benefit's withdrawal), found by applying the rules again; a participation tax rate by
// hand; the accounting identity of a population's totals; inequality and poverty over persons equal to nomos's measures;
// a reform's gainers, losers and deciles consistent with its mean change; uprating as a change of scale.

#include <nomos/microsimulation/microsimulation.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <string>
#include <vector>

using namespace nomos::microsimulation;

int main() {
  calx::verification::Suite s;

  const Schedule sch{5000, {{0, 0.1}, {10000, 0.2}, {40000, 0.4}}};
  // Taxable 45,000: 10,000 at 10%, 30,000 at 20%, 5,000 at 40% = 1,000 + 6,000 + 2,000.
  s.check("an income of 50,000 with an allowance of 5,000 pays 9,000", std::abs(tax(sch, 50000) - 9000), 1e-9);
  s.require("an income below the allowance pays nothing", tax(sch, 4000) == 0);

  Policy p;
  p.name = "an example baseline";
  p.income_tax = sch;
  p.contribution_rate = 0.08;
  p.contribution_ceiling = 60000;
  p.child_allowance = 1200;
  p.benefit = {3000, 2000, 1500, 8000, 0.6};  // maximum 3,000 + 2,000 per adult + 1,500 per child, withdrawn at 60% above 8,000

  // A couple, one earning 14,000, with two children.
  const Household h{"h1", 1, {{14000, 0, false}, {0, 0, false}, {0, 0, true}, {0, 0, true}}};
  const auto o = apply(p, h);
  // Tax: taxable 9,000 at 10% = 900. Contributions 8% of 14,000 = 1,120. Net 11,980. Maximum 3,000 + 4,000 + 3,000 = 10,000;
  // withdrawn 0.6 × (11,980 − 8,000) = 2,388; benefit 7,612. Child allowance 2,400. Disposable 21,992; scale 1 + 0.5 + 0.6 = 2.1.
  s.require("the couple's outcome by hand: tax 900, contributions 1,120, benefit 7,612, disposable 21,992, equivalised 21,992/2.1",
            std::abs(o.tax - 900) < 1e-9 && std::abs(o.contributions - 1120) < 1e-9 && std::abs(o.benefit - 7612) < 1e-9 && std::abs(o.disposable - 21992) < 1e-9 &&
                std::abs(o.equivalised - 21992 / 2.1) < 1e-9);
  // Inside the withdrawal: 1 − (1 − 0.1 − 0.08)(1 − 0.6) = 0.672.
  s.check("the marginal effective rate inside the benefit's withdrawal is 1 − (1 − 0.1 − 0.08)(1 − 0.6) = 0.672",
          std::abs(marginal_effective_rate(p, h, 0, 1) - 0.672), 1e-9);
  const Household rich{"h2", 1, {{70000, 0, false}}};
  s.check("above the contribution ceiling and the benefit, a single earner of 70,000 faces the top rate alone, 0.4", std::abs(marginal_effective_rate(p, rich, 0) - 0.4), 1e-9);
  const Household middle{"h3", 1, {{30000, 0, false}}};
  s.check("a single earner of 30,000 faces 0.2 + 0.08", std::abs(marginal_effective_rate(p, middle, 0) - 0.28), 1e-9);
  // Participation of the couple's earner: not working, net 0, benefit 10,000, allowance 2,400 → 12,400; working 21,992.
  s.check("the couple's earner keeps 9,592 of 14,000 by working: a participation tax rate of 1 − 9,592/14,000",
          std::abs(participation_rate(p, h, 0) - (1 - 9592.0 / 14000)), 1e-12);

  // ---- A population ------------------------------------------------------------------------------------------------------------
  calx::random::Stream rng(20261115, 21, 0, 0);
  std::vector<Household> pop;
  for (int i = 0; i < 300; i++) {
    Household x{"h" + std::to_string(i), 50 + 100 * rng.uniform(), {}};
    const int adults = 1 + static_cast<int>(rng.below(2)), children = static_cast<int>(rng.below(4));
    for (int a = 0; a < adults; a++) x.persons.push_back({rng.bernoulli(0.7) ? rng.lognormal(10, 0.7) : 0.0, rng.bernoulli(0.2) ? rng.lognormal(8, 1) : 0.0, false});
    for (int c = 0; c < children; c++) x.persons.push_back({0, 0, true});
    pop.push_back(x);
  }
  const auto b = simulate(p, pop);
  s.check("gross less disposable is what the state keeps net: taxes and contributions less what it pays out", std::abs((b.gross - b.disposable) - b.net_revenue) / b.gross, 1e-12);
  std::vector<double> x, w;
  for (std::size_t i = 0; i < pop.size(); i++) x.push_back(b.households[i].equivalised), w.push_back(pop[i].weight * b.households[i].persons);
  s.check("the Gini over persons is nomos's, each household's equivalised income weighted by its members", std::abs(b.gini - nomos::measurement::gini(x, w)), 1e-14);
  s.check("the poverty rate is the headcount below 60% of the median over persons", std::abs(b.poverty_rate - nomos::measurement::headcount_ratio(x, 0.6 * b.median, w)), 1e-14);

  // A reform: the child allowance doubled, paid for by a higher top rate.
  Policy r = p;
  r.name = "an example reform";
  r.child_allowance = 2400;
  r.income_tax.brackets[2].rate = 0.45;
  const auto c = compare(p, r, pop);
  double mean_change = 0, total = 0;
  const auto after = simulate(r, pop);
  for (std::size_t i = 0; i < pop.size(); i++) {
    mean_change += w[i] * (after.households[i].equivalised - b.households[i].equivalised);
    total += w[i];
  }
  mean_change /= total;
  s.check("gainers, losers and the unchanged make up every person", std::abs(c.gainers + c.losers + c.unchanged - 1), 1e-12);
  s.require("the reform raises the poorest decile's income and lowers inequality", c.decile_change[0] > 0 && c.gini_after < c.gini_before, std::to_string(c.decile_change[0]));
  s.check("the reform's revenue change is the fall in what households keep", std::abs(c.revenue_change + (after.disposable - b.disposable)), 1e-6);
  s.require("its mean change over persons lies between the lowest and highest decile's", mean_change >= *std::min_element(c.decile_change.begin(), c.decile_change.end()) - 1e-9 &&
                                                                                          mean_change <= *std::max_element(c.decile_change.begin(), c.decile_change.end()) + 1e-9);
  // Uprating with a policy scaled alike changes every amount by the same factor and leaves the Gini as it was.
  Policy scaled = p;
  scaled.income_tax.allowance *= 1.1;
  for (auto& br : scaled.income_tax.brackets) br.from *= 1.1;
  scaled.contribution_ceiling *= 1.1;
  scaled.child_allowance *= 1.1;
  scaled.benefit.base *= 1.1, scaled.benefit.per_adult *= 1.1, scaled.benefit.per_child *= 1.1, scaled.benefit.disregard *= 1.1;
  const auto up = simulate(scaled, uprate(pop, 1.1, 1.1));
  s.check("incomes and every threshold uprated by 10%: the Gini and the poverty rate stay as they were", std::abs(up.gini - b.gini) + std::abs(up.poverty_rate - b.poverty_rate), 1e-12);
  bool refused = false;
  try { (void)tax(Schedule{0, {{100, 0.1}}}, 50); } catch (const std::invalid_argument&) { refused = true; }
  s.require("a schedule whose first bracket does not start at 0 is refused", refused);
  return s.done();
}
