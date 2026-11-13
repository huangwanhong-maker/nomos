// Checks of nomos::appraisal: screening sets aside what breaks a limit whatever its other scores; the weighted sum and TOPSIS by
// hand on three alternatives; the AHP against numpy's eigensolver (Saaty's 3 × 3 and an inconsistent 4 × 4; fixtures/
// make_appraisal_references.py) and a consistent matrix exactly; PROMETHEE II's net flows summing to zero, and by hand; the weight
// stability interval at the crossing by hand; SMAA's acceptabilities summing to one, a dominated alternative never first; present
// values, IRR, payback and switching values by hand, a declining schedule's factors; distributional weights; Delphi's rounds; a
// morphological field counted by hand.

#include <nomos/appraisal/appraisal.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <vector>

using namespace nomos::appraisal;

int main() {
  calx::verification::Suite s;

  // Three alternatives on cost (lower is better), benefit and risk (lower is better).
  Matrix m;
  m.alternatives = {"A", "B", "C"};
  m.criteria = {"cost", "benefit", "risk"};
  m.score = {{100, 7, 2}, {150, 9, 1}, {120, 5, 3}};
  m.benefit = {false, true, false};
  m.weight = {0.5, 0.3, 0.2};
  {
    const auto sc = screen(m, {{2, 2.5, false}});  // risk at most 2.5
    s.require("screen: C, whose risk breaks the limit, is set aside whatever its cost", sc.set_aside.size() == 1 && sc.set_aside[0] == "C" && sc.kept.alternatives.size() == 2);
    // Min–max: cost A 1, B 0, C 0.6; benefit A 0.5, B 1, C 0; risk A 0.5, B 1, C 0.
    const auto ws = weighted_sum(m);
    s.check("weighted sum by hand: A 0.75, B 0.5, C 0.3", std::fabs(ws.value[0] - 0.75) + std::fabs(ws.value[1] - 0.5) + std::fabs(ws.value[2] - 0.3), 1e-15);
    const auto tp = topsis(m);
    // By hand: vector norms, weighted, ideal (min cost, max benefit, min risk) and anti-ideal.
    const double nc = std::sqrt(100.0 * 100 + 150 * 150 + 120 * 120), nb = std::sqrt(49.0 + 81 + 25), nr = std::sqrt(4.0 + 1 + 9);
    double c[3];
    for (int i = 0; i < 3; i++) {
      const double v[3] = {0.5 * m.score[i][0] / nc, 0.3 * m.score[i][1] / nb, 0.2 * m.score[i][2] / nr};
      const double best[3] = {0.5 * 100 / nc, 0.3 * 9 / nb, 0.2 * 1 / nr}, worst[3] = {0.5 * 150 / nc, 0.3 * 5 / nb, 0.2 * 3 / nr};
      double dp = 0, dm = 0;
      for (int j = 0; j < 3; j++) dp += (v[j] - best[j]) * (v[j] - best[j]), dm += (v[j] - worst[j]) * (v[j] - worst[j]);
      c[i] = std::sqrt(dm) / (std::sqrt(dp) + std::sqrt(dm));
    }
    s.check("TOPSIS by hand: the closeness of each to the ideal", std::fabs(tp.value[0] - c[0]) + std::fabs(tp.value[1] - c[1]) + std::fabs(tp.value[2] - c[2]), 1e-14);
  }
  {
    const auto r3 = ahp({{1, 3, 5}, {1.0 / 3, 1, 3}, {1.0 / 5, 1.0 / 3, 1}});
    s.check("AHP, Saaty's 3 × 3, against numpy: λ = 3.03851", std::fabs(r3.lambda - 3.0385110905581745), 1e-10);
    s.check("AHP: the weights 0.63699, 0.25828, 0.10473", std::fabs(r3.weight[0] - 0.6369855717447569) + std::fabs(r3.weight[1] - 0.2582849943744952) + std::fabs(r3.weight[2] - 0.10472943388074793), 1e-10);
    s.check("AHP: CR = CI / 0.58, and it is consistent (CR ≤ 0.1)", std::fabs(r3.cr - (3.0385110905581745 - 3) / 2 / 0.58) + (r3.consistent ? 0.0 : 1.0), 1e-10);
    const auto r4 = ahp({{1, 2, 4, 3}, {0.5, 1, 3, 0.5}, {0.25, 1.0 / 3, 1, 0.2}, {1.0 / 3, 2, 5, 1}});
    s.check("AHP, a 4 × 4, against numpy: λ = 4.20264, CR = CI/0.90 = 0.075, consistent enough", std::fabs(r4.lambda - 4.202640902929386) + std::fabs(r4.cr - (r4.lambda - 4) / 3 / 0.90) + (r4.consistent ? 0.0 : 1.0), 1e-10);
    // A cycle of preferences, 1 ≻ 2 ≻ 3 ≻ 1 at 9 each: λ = 1 + (a b / c)^{1/3} + (c / a b)^{1/3} = 1 + 9 + 1/9 for a = b = 9, c = 1/9.
    const auto rc3 = ahp({{1, 9, 1.0 / 9}, {1.0 / 9, 1, 9}, {9, 1.0 / 9, 1}});
    s.check("AHP: a cycle of judgements has λ = 1 + 9 + 1/9 in closed form", std::fabs(rc3.lambda - (10 + 1.0 / 9)), 1e-10);
    s.require("AHP: its CR far exceeds 0.1, and it says so", !rc3.consistent && rc3.cr > 6);
    const std::vector<double> w{0.5, 0.3, 0.2};
    std::vector<std::vector<double>> P(3, std::vector<double>(3));
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) P[i][j] = w[i] / w[j];
    const auto rc = ahp(P);
    s.check("AHP: a consistent matrix gives back its weights, λ = n, CI = 0", std::fabs(rc.weight[0] - 0.5) + std::fabs(rc.weight[2] - 0.2) + std::fabs(rc.lambda - 3) + std::fabs(rc.ci), 1e-12);
  }
  {
    const auto f = promethee(m, {Preference{}, Preference{}, Preference{}});
    s.check("PROMETHEE II: the net flows sum to zero", std::fabs(f.net[0] + f.net[1] + f.net[2]), 1e-15);
    // Usual criteria: A beats B on cost (0.5), loses on benefit and risk; π(A, B) = 0.5, π(B, A) = 0.5; A beats C on all: π(A, C) = 1;
    // B beats C on benefit and risk (0.5), loses on cost: π(B, C) = 0.5, π(C, B) = 0.5; π(C, A) = 0.
    // φ(A) = (0.5 + 1 − 0.5 − 0)/2 = 0.5; φ(B) = (0.5 + 0.5 − 0.5 − 0.5)/2 = 0; φ(C) = (0 + 0.5 − 1 − 0.5)/2 = −0.5.
    s.check("PROMETHEE II by hand, usual criteria: A 0.5, B 0, C −0.5", std::fabs(f.net[0] - 0.5) + std::fabs(f.net[1]) + std::fabs(f.net[2] + 0.5), 1e-15);
    const auto st = weight_stability(m);
    // A stays first while its weighted sum beats B's. In cost's weight t (benefit, risk scaled by (1 − t)/0.5):
    // A = t + (1 − t)(0.3·0.5 + 0.2·0.5)/0.5 = t + 0.5(1 − t); B = (1 − t)(0.3 + 0.2)/0.5 = 1 − t. Equal at t = 1/3.
    s.check("weight stability: A stays first for cost's weight above 1/3, by hand", std::fabs(st.lower[0] - 1.0 / 3) + (st.best == 0 ? 0.0 : 1.0), 1e-12);
    Matrix dm = m;
    dm.score.push_back({160, 4, 3.5});  // D is worse than C on every criterion
    dm.alternatives.push_back("D");
    const auto acc = acceptability(dm, 20000, 3);
    double sum = 0;
    for (const auto& row : acc) sum += row[0];
    s.check("SMAA: the first-rank acceptabilities sum to one", std::fabs(sum - 1), 1e-12);
    s.require("SMAA: an alternative dominated on every criterion is never first", acc[3][0] == 0);
  }
  {
    // Cost 100 now, benefits 30 a year for five years, at 5%.
    const std::vector<double> b{0, 30, 30, 30, 30, 30}, c{100, 0, 0, 0, 0, 0};
    const auto a = appraise(b, c, Schedule::constant(0.05));
    double pv = 0;
    for (int t = 1; t <= 5; t++) pv += 30 / std::pow(1.05, t);
    s.check("cost–benefit by hand: NPV = Σ 30/1.05ᵗ − 100", std::fabs(a.npv - (pv - 100)), 1e-12);
    s.check("cost–benefit: the IRR solves Σ 30/(1 + r)ᵗ = 100 (15.2382%)", std::fabs(30 * (1 - std::pow(1 + a.irr, -5)) / a.irr - 100), 1e-9);
    s.check("cost–benefit: the switching value of benefits is 100/PV, of costs PV/100", std::fabs(a.benefit_switch - 100 / pv) + std::fabs(a.cost_switch - pv / 100), 1e-14);
    s.check("cost–benefit: discounted payback in year 4", std::fabs(a.payback - 4), 0.0);
    Schedule decl{{0, 31, 76}, {0.035, 0.030, 0.025}};
    const auto d = decl.factors(80);
    s.check("a declining schedule: year 31's factor is year 30's over 1.03, year 76's year 75's over 1.025",
            std::fabs(d[31] / d[30] - 1 / 1.03) + std::fabs(d[76] / d[75] - 1 / 1.025) + std::fabs(d[30] - std::pow(1.035, -30)), 1e-15);
    const auto w = distributional({10000, 40000}, {100, 100}, 1.3, 20000);
    s.check("distributional weights: (ȳ/y)^1.3, so 100 to the poorer group counts 2^1.3 times", std::fabs(w.weight[0] - std::pow(2.0, 1.3)) + std::fabs(w.weighted_npv - 100 * (std::pow(2.0, 1.3) + std::pow(0.5, 1.3))), 1e-12);
  }
  {
    const auto r = delphi({{10, 20, 30, 40, 100}, {25, 28, 30, 32, 40}, {29, 30, 30, 31, 32}});
    s.check("Delphi: the first round's median 30 and IQR 20", std::fabs(r[0].median - 30) + std::fabs(r[0].iqr - 20), 1e-15);
    s.require("Delphi: no consensus at first; consensus and a stable median by the third round", !r[0].consensus && r[2].consensus && r[2].stable);
    Field f;
    f.parameters = {"energy", "ownership", "provision"};
    f.values = {{"grid", "local", "none"}, {"public", "private"}, {"universal", "means-tested"}};
    f.excluded = {{{0, 2}, {2, 0}}, {{1, 1}, {2, 0}}};  // no energy with universal provision; private ownership not universal
    const auto mm = morphological(f);
    // 12 in all; universal needs public ownership and energy: 2 × 1 × 1 = 2; means-tested: 3 × 2 = 6. Consistent: 8.
    s.require("morphological field by hand: 12 configurations, 8 consistent", mm.size == 12 && mm.count == 8);
  }
  return s.done();
}
