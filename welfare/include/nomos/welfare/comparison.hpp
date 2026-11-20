// Comparing candidate allocations under several principles at once, without adding the principles together.
//
// Each principle p orders the options; their intersection x ≽ y ⟺ x ≽_p y for every admissible p is a partial order, and its
// maximal set, the options no other option beats under every principle, is what can be defended without choosing among the
// principles (Sen 1997's maximal sets; the survey, Part I §2 and §7.1). For a pair of undominated options the comparison turns
// on which principles favour which: the pivot. For a principle with a parameter (inequality aversion ε, a threshold τ) the
// value at which two options exchange places is found, so that the choice of the parameter is seen as the choice it is.
//
// Only when a person declares a weighted aggregate are scores added; then the stochastic multicriteria acceptability analysis
// (SMAA-2, nomos's appraisal module) reports how often each option ranks first, second, … as the weights range over the
// simplex, beside the undominated set.
#pragma once

#include <nomos/appraisal/mcda.hpp>
#include <nomos/welfare/orders.hpp>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace nomos::welfare {

struct Pairwise {
  std::size_t a = 0, b = 0;
  std::vector<int> by_principle;   // +1 a better, −1 b better, 0 equal, for each principle
  Relation relation = Relation::unordered;  // under the intersection of the principles
};

struct Comparison {
  std::vector<std::string> options, principles;
  std::vector<std::vector<double>> values;     // [option][principle]
  std::vector<Pairwise> pairs;                 // every unordered pair
  std::vector<std::size_t> maximal;            // the undominated options
  std::vector<std::vector<std::size_t>> best_under;  // for each principle, the options best under it
};

inline Comparison compare_all(const std::vector<Outcome>& options, const std::vector<Principle>& principles) {
  Comparison c;
  for (const auto& o : options) c.options.push_back(o.name);
  for (const auto& p : principles) c.principles.push_back(p.name);
  c.values.assign(options.size(), std::vector<double>(principles.size(), 0));
  for (std::size_t i = 0; i < options.size(); i++)
    for (std::size_t k = 0; k < principles.size(); k++) c.values[i][k] = evaluate(principles[k], options[i]);
  std::vector<bool> dominated(options.size(), false);
  for (std::size_t i = 0; i < options.size(); i++)
    for (std::size_t j = i + 1; j < options.size(); j++) {
      Pairwise pw{i, j, {}, Relation::equal};
      bool better = false, worse = false;
      for (const auto& p : principles) {
        const int r = compare(p, options[i], options[j]);
        pw.by_principle.push_back(r);
        better |= r > 0, worse |= r < 0;
      }
      pw.relation = better && worse ? Relation::unordered : better ? Relation::dominates : worse ? Relation::dominated : Relation::equal;
      if (pw.relation == Relation::dominates) dominated[j] = true;
      if (pw.relation == Relation::dominated) dominated[i] = true;
      c.pairs.push_back(pw);
    }
  for (std::size_t i = 0; i < options.size(); i++)
    if (!dominated[i]) c.maximal.push_back(i);
  for (std::size_t k = 0; k < principles.size(); k++) {
    std::vector<std::size_t> best{0};
    for (std::size_t i = 1; i < options.size(); i++) {
      const int r = compare(principles[k], options[i], options[best.front()]);
      if (r > 0) best = {i};
      else if (r == 0) best.push_back(i);
    }
    c.best_under.push_back(best);
  }
  return c;
}

// The pivot of a pair: the principles favouring each option, in words.
inline std::string pivot(const Comparison& c, std::size_t a, std::size_t b) {
  for (const auto& pw : c.pairs) {
    if (!((pw.a == a && pw.b == b) || (pw.a == b && pw.b == a))) continue;
    std::string fa, fb;
    for (std::size_t k = 0; k < pw.by_principle.size(); k++) {
      const int r = pw.a == a ? pw.by_principle[k] : -pw.by_principle[k];
      if (r > 0) fa += (fa.empty() ? "" : ", ") + c.principles[k];
      if (r < 0) fb += (fb.empty() ? "" : ", ") + c.principles[k];
    }
    if (fa.empty() && fb.empty()) return c.options[a] + " and " + c.options[b] + " are equal under every principle";
    if (fb.empty()) return c.options[a] + " is better than " + c.options[b] + " under every principle that tells them apart";
    if (fa.empty()) return c.options[b] + " is better than " + c.options[a] + " under every principle that tells them apart";
    return c.options[a] + " is better under " + fa + "; " + c.options[b] + " under " + fb + ": the choice turns on which of these principles is given precedence";
  }
  return "";
}

// The parameter values in [lo, hi] at which two options exchange places under a principle with a parameter (ε or τ), found
// on a grid of `steps` and refined by bisection.
inline std::vector<double> flips(Principle p, const Outcome& a, const Outcome& b, double lo, double hi, int steps = 200) {
  std::vector<double> out;
  auto sign = [&](double x) {
    p.parameter = x;
    return compare(p, a, b, 1e-12);
  };
  double x0 = lo;
  int s0 = sign(x0);
  for (int k = 1; k <= steps; k++) {
    const double x1 = lo + (hi - lo) * k / steps;
    const int s1 = sign(x1);
    if (s0 != 0 && s1 != 0 && s0 != s1) {
      double l = x0, r = x1;
      for (int it = 0; it < 60; it++) {
        const double m = 0.5 * (l + r);
        (sign(m) == s0 ? l : r) = m;
      }
      out.push_back(0.5 * (l + r));
    }
    if (s1 != 0) x0 = x1, s0 = s1;
  }
  return out;
}

// Rank acceptabilities of a weighted sum of the principles' values as the weights range uniformly over the simplex: nomos's
// appraisal module's SMAA-2, given the options as alternatives and the principles as criteria.
inline std::vector<std::vector<double>> acceptability(const Comparison& c, const std::vector<bool>& lower_better, std::size_t draws = 10000, std::uint64_t seed = 1) {
  appraisal::Matrix m;
  m.alternatives = c.options, m.criteria = c.principles, m.score = c.values;
  m.benefit.resize(c.principles.size());
  for (std::size_t k = 0; k < c.principles.size(); k++) m.benefit[k] = !(k < lower_better.size() && lower_better[k]);
  m.weight.assign(c.principles.size(), 1.0 / static_cast<double>(c.principles.size()));
  return appraisal::acceptability(m, draws, seed);
}

}  // namespace nomos::welfare
