// Checks of psychometrics (psychometrics.hpp): Cronbach's α worked by hand (10/11), the rest correlations and α without
// each item; Spearman and Brown's prophecy and the standard error of measurement; an item's information greatest, a²/4, at
// its difficulty; marginal maximum likelihood recovering the items that simulated 2,000 persons' responses, its
// log-likelihood never falling; Rasch's model recovering a common discrimination of one, and persons with the same number
// right given the same trait (the sum score sufficient); symmetric items giving symmetric traits; an empty record the prior.

#include <nomos/psychometrics/psychometrics.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <map>
#include <string>
#include <vector>

using namespace nomos::psychometrics;

namespace {

Responses simulate(const std::vector<Item>& items, std::size_t n, calx::random::Stream& rng) {
  Responses r(n, std::vector<int>(items.size()));
  for (auto& row : r) {
    const double theta = rng.normal();
    for (std::size_t j = 0; j < items.size(); j++) row[j] = rng.bernoulli(probability(items[j], theta));
  }
  return r;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- Classical test theory ------------------------------------------------------------------------------------------------
  const std::vector<std::vector<double>> sc = {{1, 2, 3}, {2, 3, 4}, {3, 3, 3}, {4, 4, 5}};
  // Item variances 5/3, 2/3, 2.75/3 (Σ 3.25); the totals 6, 9, 9, 13 have variance 8.25; α = 1.5 (1 − 3.25/8.25) = 10/11.
  s.check("Cronbach's α of four persons on three items is 10/11, by hand", std::abs(cronbach_alpha(sc) - 10.0 / 11), 1e-14);
  const auto ia = item_analysis(sc);
  s.require("each item correlates positively with the rest, and α without an item is that of the two left",
            ia.rest_correlation.size() == 3 && ia.rest_correlation[0] > 0 && std::abs(ia.alpha_without[2] - cronbach_alpha({{1, 2}, {2, 3}, {3, 3}, {4, 4}})) < 1e-14);
  s.check("a test of reliability 0.6 doubled in length has 0.75", std::abs(spearman_brown(0.6, 2) - 0.75), 1e-15);
  s.check("with a standard deviation of 10 and reliability 0.91 the standard error of measurement is 3", std::abs(standard_error_of_measurement(10, 0.91) - 3), 1e-12);
  const Item it{1.7, 0.4};
  s.check("an item's information is greatest at its difficulty, a²/4", std::abs(information(it, 0.4) - 1.7 * 1.7 / 4), 1e-15);
  s.require("and smaller on either side", information(it, 0.0) < information(it, 0.4) && information(it, 0.9) < information(it, 0.4));

  // ---- Marginal maximum likelihood ------------------------------------------------------------------------------------------------
  calx::random::Stream rng(20261115, 31, 0, 0);
  const std::vector<Item> truth = {{0.8, -1.5}, {1.0, -1.0}, {1.2, -0.6}, {1.5, -0.2}, {2.0, 0.0}, {1.0, 0.3}, {1.3, 0.6}, {0.9, 1.0}, {1.6, 1.3}, {1.1, 1.6}};
  const auto r = simulate(truth, 2000, rng);
  const auto f = fit(r, Model::two_parameter);
  double worst_a = 0, worst_b = 0;
  for (std::size_t j = 0; j < truth.size(); j++) worst_a = std::max(worst_a, std::abs(f.items[j].a - truth[j].a)), worst_b = std::max(worst_b, std::abs(f.items[j].b - truth[j].b));
  s.require("the EM algorithm converges", f.converged, std::to_string(f.iterations) + " iterations");
  bool monotone = true;
  for (std::size_t t = 1; t < f.trace.size(); t++) monotone = monotone && f.trace[t] >= f.trace[t - 1] - 1e-9;
  s.require("its marginal log-likelihood never falls", monotone);
  s.check("the ten discriminations are recovered from 2,000 simulated persons within 0.3 (a few standard errors)", worst_a, 0.3);
  s.check("and the ten difficulties within 0.25", worst_b, 0.25);
  s.require("the fitted items give the data a higher marginal likelihood than the true ones (they are its maximum)",
            f.log_likelihood >= marginal_log_likelihood(r, truth) - 1e-6);

  // ---- Rasch's model ----------------------------------------------------------------------------------------------------------------------
  std::vector<Item> rasch_truth;
  for (double b : {-1.5, -0.8, -0.2, 0.3, 0.9, 1.4}) rasch_truth.push_back({1, b});
  const auto rr = simulate(rasch_truth, 2000, rng);
  const auto fr = fit(rr, Model::rasch);
  s.check("Rasch's model recovers a common discrimination of one", std::abs(fr.items[0].a - 1), 0.12);
  bool common = true;
  for (const auto& x : fr.items) common = common && x.a == fr.items[0].a;
  s.require("and every item shares it", common);
  const auto scores = score(rr, fr.items);
  std::map<int, double> by_sum;
  double spread = 0;
  for (std::size_t i = 0; i < rr.size(); i++) {
    int sum = 0;
    for (const int x : rr[i]) sum += x;
    if (!by_sum.count(sum)) by_sum[sum] = scores[i].eap;
    spread = std::max(spread, std::abs(by_sum[sum] - scores[i].eap));
  }
  s.check("under Rasch's model persons with the same number right have the same trait, whichever items they answered", spread, 1e-12);
  const std::vector<Item> sym = {{1, -1}, {1, 0}, {1, 1}};
  const auto ends = score({{1, 1, 1}, {0, 0, 0}, {-1, -1, -1}}, sym);
  s.check("symmetric items: everything right and everything wrong give traits of opposite sign", std::abs(ends[0].eap + ends[1].eap), 1e-12);
  s.require("and a record with no response gives the prior: mean 0, standard deviation near 1", std::abs(ends[2].eap) < 1e-12 && std::abs(ends[2].sd - 1) < 1e-3);
  bool refused = false;
  try { (void)fit({{1, 2}}, Model::rasch); } catch (const std::invalid_argument&) { refused = true; }
  s.require("a response that is not 1, 0 or −1 is refused", refused);
  return s.done();
}
