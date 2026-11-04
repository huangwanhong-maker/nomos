// Checks of index numbers (bilateral.hpp, series.hpp).
//
// The artificial data set of the Consumer Price Index Manual (ILO et al. 2004, ch. 19: six
// commodities over five periods, its Tables 19.1 and 19.2): the fixed-base Laspeyres, Paasche,
// Fisher and Törnqvist indices and the chained Laspeyres index against the values the chapter
// prints to four decimals; the chained superlative indices agree with one another closely
// (Diewert 1978). The tests of the axiomatic approach (Diewert 1992; ILO et al. 2004, ch. 16) on
// random data: identity, proportionality in current prices, invariance to proportional changes
// in current quantities, commensurability, the mean value test, time reversal (Fisher, Törnqvist,
// Walsh and Marshall–Edgeworth pass; Laspeyres fails, its product with the reversed index being
// P_L / P_P) and factor reversal (Fisher passes, Törnqvist fails). Lowe and Young with the basket
// and shares of period 0 are Laspeyres. Elementary indices: Jevons is transitive and passes time
// reversal, Carli is biased upward, CSWD passes time reversal, Dutot fails commensurability.
// Prices that bounce back: direct indices and chained Jevons return to one, chained Laspeyres
// drifts.

#include <nomos/indices/bilateral.hpp>
#include <nomos/indices/series.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <vector>

namespace ni = nomos::indices;
using ni::Formula;
using calx::verification::relative;

int main() {
  calx::verification::Suite s;
  calx::random::Stream z(20261024, 0x696e6478u, 0, 0);  // "indx"
  const std::vector<Formula> all{Formula::laspeyres, Formula::paasche, Formula::fisher, Formula::tornqvist, Formula::walsh, Formula::marshall_edgeworth};

  // The Manual's artificial data set: prices (Table 19.1) and quantities (Table 19.2).
  const ni::Periods P{{1.0, 1.0, 1.0, 1.0, 1.0, 1.0}, {1.2, 3.0, 1.3, 0.7, 1.4, 0.8}, {1.0, 1.0, 1.5, 0.5, 1.7, 0.6}, {0.8, 0.5, 1.6, 0.3, 1.9, 0.4}, {1.0, 1.0, 1.6, 0.1, 2.0, 0.2}};
  const ni::Periods Q{{1.0, 1.0, 2.0, 1.0, 4.5, 0.5}, {0.8, 0.9, 1.9, 1.3, 4.7, 0.6}, {1.0, 1.1, 1.8, 3.0, 5.0, 0.8}, {1.2, 1.2, 1.9, 6.0, 5.6, 1.3}, {0.9, 1.2, 2.0, 12.0, 6.5, 2.5}};
  {
    struct Row {
      Formula f;
      const char* name;
      std::vector<double> v;
      bool chained;
    };
    const std::vector<Row> printed{
        {Formula::laspeyres, "fixed-base Laspeyres", {1.0000, 1.4200, 1.3450, 1.3550, 1.4400}, false},
        {Formula::paasche, "fixed-base Paasche", {1.0000, 1.3824, 1.2031, 1.0209, 0.7968}, false},
        {Formula::fisher, "fixed-base Fisher", {1.0000, 1.4011, 1.2721, 1.1762, 1.0712}, false},
        {Formula::tornqvist, "fixed-base Törnqvist", {1.0000, 1.4052, 1.2890, 1.2268, 1.2477}, false},
        {Formula::laspeyres, "chained Laspeyres", {1.0000, 1.4200, 1.3646, 1.3351, 1.3306}, true},
    };
    for (const auto& r : printed) {
      const auto I = r.chained ? ni::chained_series(r.f, P, Q) : ni::direct_series(r.f, P, Q);
      double worst = 0;
      for (std::size_t t = 0; t < I.size(); t++) worst = std::max(worst, std::fabs(I[t] - r.v[t]));
      s.check(std::string("the CPI Manual's artificial data: ") + r.name + ", periods 1–5, to the printed four decimals", worst, 5.01e-5);
    }
    const auto F = ni::chained_series(Formula::fisher, P, Q), T = ni::chained_series(Formula::tornqvist, P, Q);
    double worst = 0;
    for (std::size_t t = 0; t < F.size(); t++) worst = std::max(worst, std::fabs(F[t] - T[t]));
    s.check("... chained Fisher and chained Törnqvist stay within 0.01 of each other (Diewert 1978)", worst, 1e-2);
    s.check("... and end within 3·10⁻⁴ (1.2226 and 1.2224 in period 5)", F[4] - T[4], 3e-4);
    const auto Ld = ni::direct_series(Formula::laspeyres, P, Q), Pd = ni::direct_series(Formula::paasche, P, Q);
    s.require("... the spread of fixed-base Laspeyres and Paasche (1.4400 against 0.7968) is far wider", Ld[4] - Pd[4] > 0.6);
    const auto rebased = ni::direct_series(Formula::fisher, P, Q, 2);
    s.check("... rebased to period 3, the direct Fisher series at period 1 is 1/P_F(1 → 3), by time reversal",
            relative(rebased[0], 1 / ni::fisher(P[0], P[2], Q[0], Q[2])), 1e-15);
  }

  // The axiomatic tests on random data, N = 12 items.
  const std::size_t N = 12;
  auto vec = [&](double lo, double hi) {
    std::vector<double> v(N);
    for (double& x : v) x = z.uniform(lo, hi);
    return v;
  };
  const auto p0 = vec(0.5, 2), p1 = vec(0.5, 3), q0 = vec(0.1, 5), q1 = vec(0.1, 5);
  {
    double identity = 0, proportional = 0, quantities = 0, units = 0;
    bool mean_value = true;
    std::vector<double> lam = vec(0.01, 100), p0u(N), p1u(N), q0u(N), q1u(N), lp0(N), lq1(N);
    for (std::size_t i = 0; i < N; i++) {
      p0u[i] = lam[i] * p0[i];
      p1u[i] = lam[i] * p1[i];
      q0u[i] = q0[i] / lam[i];
      q1u[i] = q1[i] / lam[i];
      lp0[i] = 1.7 * p0[i];
      lq1[i] = 3.1 * q1[i];
    }
    double lo = 1e300, hi = 0;
    for (std::size_t i = 0; i < N; i++) lo = std::min(lo, p1[i] / p0[i]), hi = std::max(hi, p1[i] / p0[i]);
    for (Formula f : all) {
      const double v = ni::price_index(f, p0, p1, q0, q1);
      identity = std::max(identity, std::fabs(ni::price_index(f, p0, p0, q0, q1) - 1));
      proportional = std::max(proportional, relative(ni::price_index(f, p0, lp0, q0, q1), 1.7));
      if (f != Formula::marshall_edgeworth) quantities = std::max(quantities, relative(ni::price_index(f, p0, p1, q0, lq1), v));
      units = std::max(units, relative(ni::price_index(f, p0u, p1u, q0u, q1u), v));
      mean_value = mean_value && lo <= v && v <= hi;
    }
    s.check("identity: P(p, p, q⁰, q¹) = 1 for all six formulas", identity, 1e-15);
    s.check("proportionality: P(p⁰, 1.7 p⁰, q⁰, q¹) = 1.7 for all six", proportional, 1e-15);
    s.check("invariance to proportional changes in current quantities (all but Marshall–Edgeworth)", quantities, 1e-14);
    s.require("... which Marshall–Edgeworth fails", relative(ni::marshall_edgeworth(p0, p1, q0, lq1), ni::marshall_edgeworth(p0, p1, q0, q1)) > 1e-4);
    s.check("commensurability: a change of the units of each item changes none of the six", units, 1e-13);
    s.require("the mean value test: every index lies between the smallest and largest price relative", mean_value);
    s.require("... while Dutot's elementary index changes with the units", relative(ni::dutot(p0u, p1u), ni::dutot(p0, p1)) > 1e-3);
  }
  {
    double reversal = 0;
    for (Formula f : {Formula::fisher, Formula::tornqvist, Formula::walsh, Formula::marshall_edgeworth})
      reversal = std::max(reversal, std::fabs(ni::price_index(f, p0, p1, q0, q1) * ni::price_index(f, p1, p0, q1, q0) - 1));
    s.check("time reversal: Fisher, Törnqvist, Walsh and Marshall–Edgeworth", reversal, 1e-14);
    const double LL = ni::laspeyres(p0, p1, q0, q1) * ni::laspeyres(p1, p0, q1, q0);
    s.check("... Laspeyres fails it: P_L(0, 1) P_L(1, 0) = P_L / P_P", relative(LL, ni::laspeyres(p0, p1, q0, q1) / ni::paasche(p0, p1, q0, q1)), 1e-14);
    s.require("... and is not 1 here", std::fabs(LL - 1) > 1e-4);
    const double V = ni::value_ratio(p0, p1, q0, q1);
    s.check("factor reversal: Fisher's price and quantity indices multiply to the value ratio",
            relative(ni::fisher(p0, p1, q0, q1) * ni::quantity_index(Formula::fisher, p0, p1, q0, q1), V), 1e-14);
    s.require("... Törnqvist's do not", relative(ni::tornqvist(p0, p1, q0, q1) * ni::quantity_index(Formula::tornqvist, p0, p1, q0, q1), V) > 1e-5);
    s.check("Laspeyres' quantity index is Σ p⁰q¹ / Σ p⁰q⁰",
            relative(ni::quantity_index(Formula::laspeyres, p0, p1, q0, q1), [&] {
              double a = 0, b = 0;
              for (std::size_t i = 0; i < N; i++) a += p0[i] * q1[i], b += p0[i] * q0[i];
              return a / b;
            }()),
            1e-15);
    std::vector<double> s0(N);
    for (std::size_t i = 0; i < N; i++) s0[i] = p0[i] * q0[i];
    s.check("Lowe with the basket of period 0 is Laspeyres", relative(ni::lowe(p0, p1, q0), ni::laspeyres(p0, p1, q0, q1)), 1e-15);
    s.check("Young with the shares of period 0 is Laspeyres", relative(ni::young(p0, p1, s0), ni::laspeyres(p0, p1, q0, q1)), 1e-14);
  }
  // Elementary indices.
  {
    const auto e0 = vec(1, 2), e1 = vec(1, 3), e2 = vec(0.5, 2);
    s.check("Jevons is transitive: J(0, 2) = J(0, 1) J(1, 2)", relative(ni::jevons(e0, e2), ni::jevons(e0, e1) * ni::jevons(e1, e2)), 1e-14);
    s.check("... and passes time reversal", std::fabs(ni::jevons(e0, e1) * ni::jevons(e1, e0) - 1), 1e-15);
    s.require("Carli fails it upward: C(0, 1) C(1, 0) > 1", ni::carli(e0, e1) * ni::carli(e1, e0) > 1 + 1e-6);
    s.check("the harmonic index is the reciprocal of the reversed Carli", relative(ni::harmonic(e0, e1), 1 / ni::carli(e1, e0)), 1e-15);
    s.check("CSWD passes time reversal", std::fabs(ni::cswd(e0, e1) * ni::cswd(e1, e0) - 1), 1e-15);
    s.require("harmonic ≤ Jevons ≤ Carli (the means)", ni::harmonic(e0, e1) <= ni::jevons(e0, e1) && ni::jevons(e0, e1) <= ni::carli(e0, e1));
  }
  // Bouncing prices: period 3 repeats period 0.
  {
    const ni::Periods p{{1, 1, 1}, {1.5, 0.8, 1.2}, {0.9, 1.4, 1.1}, {1, 1, 1}};
    const ni::Periods q{{2, 3, 1}, {1.6, 3.9, 1.0}, {2.4, 2.1, 1.3}, {2, 3, 1}};
    double direct = 0;
    for (Formula f : all) direct = std::max(direct, std::fabs(ni::direct_series(f, p, q)[3] - 1));
    s.check("prices and quantities that bounce back: every direct index returns to 1", direct, 1e-15);
    const double L = ni::chained_series(Formula::laspeyres, p, q)[3], F = ni::chained_series(Formula::fisher, p, q)[3];
    s.require("... chained Laspeyres drifts (Szulc 1983)", std::fabs(L - 1) > 1e-3, std::to_string(L));
    s.require("... and chained Fisher drifts less", std::fabs(F - 1) < std::fabs(L - 1), std::to_string(F));
    double J = 1;
    for (std::size_t t = 1; t < p.size(); t++) J *= ni::jevons(p[t - 1], p[t]);
    s.check("... while chained Jevons, transitive, returns to 1", std::fabs(J - 1), 1e-15);
  }
  return s.done();
}
