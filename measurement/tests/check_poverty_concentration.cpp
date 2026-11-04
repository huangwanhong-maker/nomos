// Checks of poverty (poverty.hpp), concentration (concentration.hpp) and specialisation
// (specialisation.hpp).
//
// Poverty: the amounts 1, …, 10 at the line 5.5, by hand: headcount ½, gap 12.5/55, severity
// 41.25/302.5, Watts (1/10) Σ ln(5.5/x); a unit exactly at the line is not poor; the axioms of
// monotonicity and of transfers (Sen 1976; Foster, Greer and Thorbecke 1984).
// Concentration: shares 50, 30, 15, 5 per cent: H = 0.365 (3650 in per cent), CR₂ = 0.8; equal
// shares give H = 1/n, E = ln n and every Hannah–Kay number n; Hannah and Kay's numbers decrease
// in α, equal 1/H at α = 2 and are continuous at α = 1, where they are exp(E).
// Specialisation: location quotients of a 2 × 3 table by hand, and their two weighted means of
// one on a random table; the index of dissimilarity at its ends; the Grubel–Lloyd indices of a
// three-category example by hand, the adjusted aggregate one when every category is as balanced
// as the totals allow.

#include <nomos/measurement/concentration.hpp>
#include <nomos/measurement/poverty.hpp>
#include <nomos/measurement/specialisation.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <vector>

namespace nm = nomos::measurement;
using calx::verification::relative;

int main() {
  calx::verification::Suite s;
  calx::random::Stream z(20261024, 0x6e6f6d6fu, 2, 0);

  // Poverty.
  {
    const std::vector<double> x{7, 1, 9, 2, 10, 3, 8, 4, 6, 5};
    const double line = 5.5;
    s.check("1, …, 10 at the line 5.5: headcount ratio ½", nm::headcount_ratio(x, line) - 0.5, 1e-15);
    s.check("... poverty gap (4.5 + 3.5 + 2.5 + 1.5 + 0.5)/(5.5·10)", nm::poverty_gap_index(x, line) - 12.5 / 55, 1e-15);
    s.check("... severity (4.5² + … + 0.5²)/(5.5²·10)", nm::poverty_severity(x, line) - 41.25 / 302.5, 1e-15);
    double watts = 0;
    for (double v : {1.0, 2.0, 3.0, 4.0, 5.0}) watts += std::log(line / v);
    s.check("... Watts (1/10) Σ ln(5.5/x)", nm::watts_index(x, line) - watts / 10, 1e-15);
    s.check("... and FGT(α) at α = 3 by its definition", nm::fgt(x, line, 3) - (std::pow(4.5, 3) + std::pow(3.5, 3) + std::pow(2.5, 3) + std::pow(1.5, 3) + std::pow(0.5, 3)) / std::pow(5.5, 3) / 10, 1e-15);
    s.check("a unit exactly at the line 5 is not poor: headcount 4/10", nm::headcount_ratio(x, 5.0) - 0.4, 1e-15);
    // Weights as replication.
    const std::vector<double> w{1, 2, 1, 3, 1, 1, 1, 2, 1, 1};
    std::vector<double> rep;
    for (std::size_t i = 0; i < x.size(); i++)
      for (int k = 0; k < static_cast<int>(w[i]); k++) rep.push_back(x[i]);
    double worst = 0;
    for (double a : {0.0, 1.0, 2.0, 2.5}) worst = std::max(worst, relative(nm::fgt(x, line, a, w), nm::fgt(rep, line, a)));
    worst = std::max(worst, relative(nm::watts_index(x, line, w), nm::watts_index(rep, line)));
    s.check("integer weights are replication: FGT(0, 1, 2, 2.5) and Watts", worst, 1e-14);
    // Monotonicity: a poor unit loses 0.5. Transfers: 0.5 from the poor unit at 4 to the one at 1.
    std::vector<double> poorer = x, moved = x;
    poorer[7] -= 0.5;  // 4 → 3.5
    moved[7] -= 0.5;   // 4 → 3.5
    moved[1] += 0.5;   // 1 → 1.5
    s.require("monotonicity: a poorer poor unit raises the gap, the severity and Watts",
              nm::poverty_gap_index(poorer, line) > nm::poverty_gap_index(x, line) && nm::poverty_severity(poorer, line) > nm::poverty_severity(x, line) &&
                  nm::watts_index(poorer, line) > nm::watts_index(x, line));
    s.check("transfers: a transfer among the poor leaves the gap unchanged", nm::poverty_gap_index(moved, line) - nm::poverty_gap_index(x, line), 1e-15);
    s.require("... and lowers the severity and Watts",
              nm::poverty_severity(moved, line) < nm::poverty_severity(x, line) && nm::watts_index(moved, line) < nm::watts_index(x, line));
  }
  // Concentration.
  {
    const std::vector<double> x{30, 5, 50, 15};
    s.check("shares 50, 30, 15, 5 per cent: H = 0.365 (3650 on per cent)", nm::herfindahl_hirschman(x) - 0.365, 1e-15);
    s.check("... normalised (0.365 − ¼)/(1 − ¼)", nm::herfindahl_hirschman_normalised(x) - (0.365 - 0.25) / 0.75, 1e-15);
    s.check("... CR₂ = 0.8", nm::concentration_ratio(x, 2) - 0.8, 1e-15);
    s.check("... CR₁₀ = 1 (k beyond n)", nm::concentration_ratio(x, 10) - 1, 1e-15);
    s.check("... Hannah–Kay at α = 2 is 1/H", relative(nm::hannah_kay(x, 2), 1 / 0.365), 1e-14);
    s.check("... Hannah–Kay is continuous at α = 1, where it is exp(E)", relative(nm::hannah_kay(x, 1 + 1e-8), nm::hannah_kay(x, 1)), 1e-7);
    bool decreasing = true;
    double prev = nm::hannah_kay(x, 0.25);
    for (double a : {0.5, 1.0, 1.5, 2.0, 3.0, 8.0}) {
      const double v = nm::hannah_kay(x, a);
      decreasing = decreasing && v < prev && v >= 1 && v <= 4;
      prev = v;
    }
    s.require("... Hannah–Kay decreases in α within [1, n]", decreasing);
    const std::vector<double> eq(7, 2.5);
    s.check("seven equal holders: H = 1/7", nm::herfindahl_hirschman(eq) - 1.0 / 7, 1e-15);
    s.check("... normalised H = 0", nm::herfindahl_hirschman_normalised(eq), 1e-15);
    s.check("... E = ln 7", nm::entropy_index(eq) - std::log(7.0), 1e-15);
    s.check("... Hannah–Kay at α = ½ is 7", nm::hannah_kay(eq, 0.5) - 7, 1e-13);
    const std::vector<double> one{0, 0, 9, 0};
    s.require("a single holder: H = 1, normalised 1, E = 0",
              nm::herfindahl_hirschman(one) == 1 && nm::herfindahl_hirschman_normalised(one) == 1 && nm::entropy_index(one) == 0);
  }
  // Specialisation: location quotients.
  {
    // Rows 100 and 100; columns 40, 60, 100; total 200.
    const auto t = calx::linalg::Rectangular::from_rows({{10, 20, 70}, {30, 40, 30}});
    const auto q = nm::location_quotients(t);
    s.check("a 2 × 3 table: LQ₀₀ = (10/100)/(40/200) = ½", q(0, 0) - 0.5, 1e-15);
    s.check("... LQ₀₂ = 0.7/0.5 = 1.4", q(0, 2) - 1.4, 1e-15);
    s.check("... LQ₁₀ = 0.3/0.2 = 1.5", q(1, 0) - 1.5, 1e-15);
    const auto sym = nm::symmetric_location_quotients(q);
    s.check("... symmetric (½ − 1)/(½ + 1) = −⅓", sym(0, 0) + 1.0 / 3, 1e-15);
    // Random 6 × 9 table: the two weighted means of one.
    const std::size_t R = 6, C = 9;
    calx::linalg::Rectangular x(R, C);
    std::vector<double> row(R, 0), col(C, 0);
    double total = 0;
    for (std::size_t r = 0; r < R; r++)
      for (std::size_t c = 0; c < C; c++) {
        x(r, c) = 100 * z.uniform();
        row[r] += x(r, c);
        col[c] += x(r, c);
        total += x(r, c);
      }
    const auto lq = nm::balassa_index(x);
    double worst = 0;
    for (std::size_t c = 0; c < C; c++) {
      double m = 0;
      for (std::size_t r = 0; r < R; r++) m += row[r] / total * lq(r, c);
      worst = std::max(worst, std::fabs(m - 1));
    }
    for (std::size_t r = 0; r < R; r++) {
      double m = 0;
      for (std::size_t c = 0; c < C; c++) m += col[c] / total * lq(r, c);
      worst = std::max(worst, std::fabs(m - 1));
    }
    s.check("a random 6 × 9 table: Σ_r (x_r·/x··) LQ_rc = 1 and Σ_c (x_·c/x··) LQ_rc = 1", worst, 1e-13);
    auto holes = t;
    holes(0, 0) = holes(1, 0) = 0;  // column 0 empty
    s.require("... an empty column gives an undefined quotient (NaN)", std::isnan(nm::location_quotients(holes)(1, 0)));
  }
  // Dissimilarity and Grubel–Lloyd.
  {
    const std::vector<double> a{1, 0}, b{0, 3}, c{2, 6, 4}, d{1, 3, 2};
    s.check("dissimilarity of disjoint distributions is 1", nm::dissimilarity(a, b) - 1, 1e-15);
    s.check("... of proportional ones 0", nm::dissimilarity(c, d), 1e-15);
    s.check("... of (1, 1) and (1, 3): ½|½ − ¼| · 2 = ¼", nm::dissimilarity(std::vector<double>{1, 1}, std::vector<double>{1, 3}) - 0.25, 1e-15);
    const std::vector<double> out{10, 0, 5}, in{10, 8, 3};
    s.require("Grubel–Lloyd by category: 1, 0 and 1 − 2/8",
              nm::grubel_lloyd(10, 10) == 1 && nm::grubel_lloyd(0, 8) == 0 && nm::grubel_lloyd(5, 3) == 0.75);
    s.check("... aggregate (36 − 10)/36", nm::grubel_lloyd(out, in) - 26.0 / 36, 1e-15);
    s.check("... adjusted for the imbalance 21 − 15: 26/(36 − 6)", nm::grubel_lloyd_adjusted(out, in) - 26.0 / 30, 1e-15);
    const std::vector<double> o2{10, 5}, i2{4, 5};
    s.check("... as balanced as the totals allow: the adjusted aggregate is 1", nm::grubel_lloyd_adjusted(o2, i2) - 1, 1e-15);
  }
  return s.done();
}
