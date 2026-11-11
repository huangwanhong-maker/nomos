// Checks of spatial interaction (gravity.hpp, lowry.hpp, combined.hpp).
//
// - Furness's balancing meets both margins, and the balanced table is a_i b_j e^{−βc_ij}: ln T + βc is a sum of a row's term
//   and a column's (its double-centred residual vanishes).
// - Wilson's theorem: the doubly constrained table has more entropy than every nearby table with the same margins and the
//   same total cost (perturbations drawn in the null space of the constraints).
// - The production- and attraction-constrained tables meet their one margin; the production-constrained table with the
//   balancing factors as attractions is the doubly constrained table.
// - The mean cost falls as β rises.
// - Hyman's calibration to an observed table gives statsmodels' Poisson GLM with origin and destination effects (β and the
//   fitted table; fixtures/make_interaction_fixtures.py), and nomos::econometrics' PPML with absorbed effects.
// - The logsum: ∂A_i/∂c_ij = −P_ij and ∂A_i/∂ln W_j = P_ij/β by central differences; it lies between the best destination's
//   utility and that plus ln J / β.
// - The common part of commuters and the SRMSE at their bounds.
// - Lowry's iteration is Garin's closed form, Σ e = Σ e^B / (1 − ab) and Σ p = a Σ e.
// - Evans's combined model on a grid of streets: at the end the table is the doubly constrained table at the network's
//   least times and the arc flows are Wardrop's equilibrium for it; the lower bound lies below the objective; distributing
//   at free flow and assigning once leaves a larger objective.
// - On Braess's network, whose added road takes no time, the combined model's arc flows carry the table (a node and its
//   parent in a shortest-path tree at one distance) and the journey takes Braess's 80 minutes.

#include <nomos/econometrics/poisson_fixed_effects.hpp>
#include <nomos/interaction/combined.hpp>
#include <nomos/interaction/gravity.hpp>
#include <nomos/interaction/lowry.hpp>

#include <calx/linalg/dense.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "fixtures/interaction_fixtures.hpp"

namespace ni = nomos::interaction;
namespace ng = nomos::games;

namespace {
double entropy(const ni::Table& t) {
  double h = 0;
  for (const double x : t.v)
    if (x > 0) h -= x * (std::log(x) - 1);
  return h;
}
}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream st(20261111, 0x696e7472u, 0, 0);  // "intr"
  const std::size_t R = 7, C = 6;
  std::vector<double> ox(R), oy(R), dx(C), dy(C), O(R), D(C);
  for (std::size_t i = 0; i < R; i++) ox[i] = st.uniform(0, 10), oy[i] = st.uniform(0, 10), O[i] = st.uniform(100, 900);
  for (std::size_t j = 0; j < C; j++) dx[j] = st.uniform(0, 10), dy[j] = st.uniform(0, 10), D[j] = st.uniform(100, 900);
  double so = 0, sd = 0;
  for (const double x : O) so += x;
  for (const double x : D) sd += x;
  for (double& x : D) x *= so / sd;
  ni::Table cost(R, C);
  for (std::size_t i = 0; i < R; i++)
    for (std::size_t j = 0; j < C; j++) cost(i, j) = std::hypot(ox[i] - dx[j], oy[i] - dy[j]) + 0.3;
  const double beta = 0.4;
  const auto dc = ni::doubly_constrained(O, D, cost, beta);
  {
    const auto r = dc.flows.row_sums(), c = dc.flows.col_sums();
    double e = 0;
    for (std::size_t i = 0; i < R; i++) e = std::max(e, std::fabs(r[i] - O[i]) / O[i]);
    for (std::size_t j = 0; j < C; j++) e = std::max(e, std::fabs(c[j] - D[j]) / D[j]);
    s.require("Furness's balancing converges", dc.converged);
    s.check("Furness's balancing meets both margins", e, 1e-11);
    // ln T_ij + β c_ij = α_i + γ_j: the double-centred residual.
    ni::Table L(R, C);
    for (std::size_t k = 0; k < L.v.size(); k++) L.v[k] = std::log(dc.flows.v[k]) + beta * cost.v[k];
    std::vector<double> rm(R, 0), cm(C, 0);
    double all = 0;
    for (std::size_t i = 0; i < R; i++)
      for (std::size_t j = 0; j < C; j++) rm[i] += L(i, j) / C, cm[j] += L(i, j) / R, all += L(i, j) / (R * C);
    double worst = 0;
    for (std::size_t i = 0; i < R; i++)
      for (std::size_t j = 0; j < C; j++) worst = std::max(worst, std::fabs(L(i, j) - rm[i] - cm[j] + all));
    s.check("The balanced table is a_i b_j e^{−βc}: ln T + βc has no interaction of row and column", worst, 1e-11);
  }
  {
    // The null space of the constraints (R row sums, C column sums, the total cost), by projection.
    const std::size_t n = R * C, m = R + C + 1;
    std::vector<std::vector<double>> K(m, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < R; i++)
      for (std::size_t j = 0; j < C; j++) K[i][i * C + j] = 1, K[R + j][i * C + j] = 1, K[R + C][i * C + j] = cost(i, j);
    calx::linalg::Matrix G(m);
    for (std::size_t a = 0; a < m; a++)
      for (std::size_t b = 0; b < m; b++)
        for (std::size_t k = 0; k < n; k++) G(a, b) += K[a][k] * K[b][k];
    // One row sum is implied by the others and the columns: G is singular; regularise the redundant direction.
    for (std::size_t a = 0; a < m; a++) G(a, a) += 1e-9;
    const double h0 = entropy(dc.flows);
    bool beaten = false;
    double drift = 0;
    for (int trial = 0; trial < 50; trial++) {
      std::vector<double> z(n);
      for (double& x : z) x = st.uniform(-1, 1);
      calx::linalg::Matrix B(m);
      for (std::size_t a = 0; a < m; a++)
        for (std::size_t k = 0; k < n; k++) B(a, 0) += K[a][k] * z[k];
      const auto y = calx::linalg::solve(G, B);
      for (std::size_t k = 0; k < n; k++)
        for (std::size_t a = 0; a < m; a++) z[k] -= K[a][k] * y(a, 0);
      ni::Table t = dc.flows;
      for (std::size_t k = 0; k < n; k++) t.v[k] += 0.5 * z[k];
      double total = 0, base = 0;
      for (std::size_t k = 0; k < n; k++) total += t.v[k] * cost.v[k], base += dc.flows.v[k] * cost.v[k];
      drift = std::max(drift, std::fabs(total - base) / base);
      if (entropy(t) >= h0) beaten = true;
    }
    s.check("Wilson: the perturbations keep the total cost", drift, 1e-9);
    s.require("Wilson: no table nearby with the same margins and total cost has more entropy than the doubly constrained one", !beaten);
  }
  {
    const auto pc = ni::production_constrained(O, D, cost, beta);
    const auto ac = ni::attraction_constrained(D, O, cost, beta);
    const auto r = pc.row_sums(), c = ac.col_sums();
    double e = 0;
    for (std::size_t i = 0; i < R; i++) e = std::max(e, std::fabs(r[i] - O[i]) / O[i]);
    for (std::size_t j = 0; j < C; j++) e = std::max(e, std::fabs(c[j] - D[j]) / D[j]);
    s.check("The production- and attraction-constrained tables meet their margins", e, 1e-12);
    // The balanced seed differs from e^{−βc} by a factor in each row, which the production constraint absorbs: W = b.
    const auto same = ni::production_constrained(O, dc.b, cost, beta);
    double d = 0;
    for (std::size_t k = 0; k < same.v.size(); k++) d = std::max(d, std::fabs(same.v[k] - dc.flows.v[k]) / dc.flows.v[k]);
    s.check("The production-constrained table with the balancing factors as attractions is the doubly constrained table", d, 1e-10);
  }
  {
    bool falls = true;
    double last = 1e300;
    for (double b = 0.02; b < 3; b *= 1.3) {
      const double m = ni::mean_cost(ni::doubly_constrained(O, D, cost, b).flows, cost);
      falls = falls && m < last;
      last = m;
    }
    s.require("The doubly constrained mean cost falls as β rises", falls);
  }
  {
    namespace f = interaction_fixtures;
    ni::Table cost_f(f::rows, f::cols), obs(f::rows, f::cols);
    for (std::size_t k = 0; k < cost_f.v.size(); k++) cost_f.v[k] = f::cost[k], obs.v[k] = f::observed[k];
    const auto cal = ni::calibrate(obs, cost_f);
    s.require("Hyman's calibration converges", cal.converged);
    s.check("Hyman's β is statsmodels' Poisson maximum likelihood estimate", (cal.beta + f::coefficient) / f::coefficient, 1e-8);
    double e = 0;
    for (std::size_t k = 0; k < cost_f.v.size(); k++) e = std::max(e, std::fabs(cal.fit.flows.v[k] - f::fitted[k]) / f::fitted[k]);
    s.check("Hyman's table is statsmodels' fitted table", e, 1e-7);
    calx::linalg::Rectangular X(cost_f.v.size(), 1);
    std::vector<std::size_t> row(cost_f.v.size()), col(cost_f.v.size());
    for (std::size_t i = 0; i < cost_f.rows; i++)
      for (std::size_t j = 0; j < cost_f.cols; j++) X(i * cost_f.cols + j, 0) = cost_f(i, j), row[i * cost_f.cols + j] = i, col[i * cost_f.cols + j] = j;
    nomos::econometrics::PoissonTwoWayOptions po;
    po.tolerance = 1e-14;
    const auto ppml = nomos::econometrics::poisson_two_way(X, obs.v, row, col, po);
    s.check("Hyman's β is nomos's PPML with origin and destination effects", (cal.beta + ppml.beta[0]) / ppml.beta[0], 1e-7);
  }
  {
    std::vector<double> W(C);
    for (double& w : W) w = st.uniform(10, 200);
    const auto A = ni::logsum(W, cost, beta);
    std::vector<double> ones(R, 1.0);
    const auto P = ni::production_constrained(ones, W, cost, beta);
    double e1 = 0, e2 = 0, bound = 0;
    const double h = 1e-5;
    for (std::size_t i = 0; i < R; i++)
      for (std::size_t j = 0; j < C; j++) {
        ni::Table up = cost, down = cost;
        up(i, j) += h, down(i, j) -= h;
        e1 = std::max(e1, std::fabs((ni::logsum(W, up, beta)[i] - ni::logsum(W, down, beta)[i]) / (2 * h) + P(i, j)));
        auto Wu = W, Wd = W;
        Wu[j] *= std::exp(h), Wd[j] *= std::exp(-h);
        e2 = std::max(e2, std::fabs((ni::logsum(Wu, cost, beta)[i] - ni::logsum(Wd, cost, beta)[i]) / (2 * h) - P(i, j) / beta));
      }
    for (std::size_t i = 0; i < R; i++) {
      double best = -1e300;
      for (std::size_t j = 0; j < C; j++) best = std::max(best, std::log(W[j]) / beta - cost(i, j));
      if (A[i] < best - 1e-12 || A[i] > best + std::log(static_cast<double>(C)) / beta + 1e-12) bound = 1;
    }
    s.check("The logsum: ∂A_i/∂c_ij = −P_ij", e1, 1e-8);
    s.check("The logsum: ∂A_i/∂ln W_j = P_ij/β", e2, 1e-8);
    s.check("The logsum lies between the best destination's utility and that plus ln J / β", bound, 0);
  }
  {
    ni::Table other = dc.flows;
    for (double& x : other.v) x *= 0;
    other(0, 0) = 1;
    ni::Table lone(R, C);
    lone(1, 1) = 1;
    s.check("The common part of a table with itself is 1, and its SRMSE 0", std::fabs(ni::common_part(dc.flows, dc.flows) - 1) + ni::srmse(dc.flows, dc.flows), 1e-15);
    s.check("The common part of two tables with no pair in common is 0", ni::common_part(other, lone), 0);
  }
  {
    const std::size_t n = 8;
    ni::Table c8(n, n);
    std::vector<double> x(n), y(n), basic(n), weight(n), ones(n, 1.0);
    for (std::size_t i = 0; i < n; i++) x[i] = st.uniform(0, 10), y[i] = st.uniform(0, 10), basic[i] = st.uniform(0, 1000), weight[i] = st.uniform(1, 10);
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j < n; j++) c8(i, j) = std::hypot(x[i] - x[j], y[i] - y[j]) + 0.5;
    const auto Rm = ni::production_constrained(ones, weight, c8, 0.3);
    const auto Sm = ni::production_constrained(ones, basic, c8, 0.6);
    const double a = 2.2, b = 0.18;
    const auto it = ni::lowry(basic, Rm, Sm, a, b);
    const auto closed = ni::garin(basic, Rm, Sm, a, b);
    double e = 0, tb = 0, te = 0, tp = 0;
    for (std::size_t i = 0; i < n; i++) e = std::max(e, std::fabs(it.employment[i] - closed.employment[i]) / closed.employment[i]), tb += basic[i], te += closed.employment[i], tp += closed.population[i];
    s.check("Lowry's iteration is Garin's closed form", e, 1e-11);
    s.check("Lowry: total employment is basic employment over 1 − ab", (te - tb / (1 - a * b)) / te, 1e-12);
    s.check("Lowry: total population is a times total employment", (tp - a * te) / tp, 1e-12);
  }
  {
    // A grid of four by four crossings, streets both ways; four zones of homes at the corners, three of work inside.
    ng::Network g{16, {}, {}};
    for (std::size_t r = 0; r < 4; r++)
      for (std::size_t c = 0; c < 4; c++) {
        const std::size_t v = r * 4 + c;
        if (c + 1 < 4) g.arcs.push_back(ng::bpr(v, v + 1, st.uniform(1, 3), st.uniform(300, 900))), g.arcs.push_back(ng::bpr(v + 1, v, st.uniform(1, 3), st.uniform(300, 900)));
        if (r + 1 < 4) g.arcs.push_back(ng::bpr(v, v + 4, st.uniform(1, 3), st.uniform(300, 900))), g.arcs.push_back(ng::bpr(v + 4, v, st.uniform(1, 3), st.uniform(300, 900)));
      }
    const std::vector<std::size_t> homes{0, 3, 12, 15}, work{5, 6, 10};
    const std::vector<double> Oh{1200, 800, 900, 1100}, Dw{1500, 1000, 1500};
    const auto eq = ni::combined(g, homes, work, Oh, Dw, 0.25, {1e-6, 1000000});
    s.require("Evans: the combined model converges to a gap of 10⁻⁶", eq.converged, std::to_string(eq.iterations) + " iterations, gap " + std::to_string(eq.gap));
    const auto again = ni::doubly_constrained(Oh, Dw, eq.costs, 0.25).flows;
    double e = 0;
    for (std::size_t k = 0; k < again.v.size(); k++) e = std::max(e, std::fabs(again.v[k] - eq.flows.v[k]) / eq.flows.v[k]);
    s.check("Evans: the table is the doubly constrained table at the network's least times", e, 5e-3);
    ng::Network loaded = g;
    for (std::size_t i = 0; i < homes.size(); i++)
      for (std::size_t j = 0; j < work.size(); j++) loaded.demands.push_back({homes[i], work[j], eq.flows(i, j)});
    s.check("Evans: the arc flows are Wardrop's equilibrium for that table (relative gap)", ng::relative_gap(loaded, eq.arc_flow), 1e-5);
    s.require("Evans: the lower bound lies below the objective", eq.lower_bound <= eq.objective);
    // Distribute at free flow, assign once.
    const auto free = ni::doubly_constrained(Oh, Dw, [&] {
      ni::Table c(homes.size(), work.size());
      const auto sk = ng::skims(g, ng::times(g, std::vector<double>(g.arcs.size(), 0.0)), homes, work);
      for (std::size_t i = 0; i < c.rows; i++)
        for (std::size_t j = 0; j < c.cols; j++) c(i, j) = sk[i][j];
      return c;
    }(), 0.25).flows;
    ng::Network once = g;
    for (std::size_t i = 0; i < homes.size(); i++)
      for (std::size_t j = 0; j < work.size(); j++) once.demands.push_back({homes[i], work[j], free(i, j)});
    const auto a1 = ng::assign(once, {ng::Objective::user, 1e-9});
    double z1 = a1.objective;
    for (const double x : free.v) z1 += x * (std::log(x) - 1) / 0.25;
    s.require("Evans: distributing at free flow and assigning once leaves a larger objective", z1 > eq.objective);
  }
  {
    // Braess's network with its road of no time, 4,000 travellers from S to E: one pair, so the combined model is Wardrop's
    // equilibrium, everyone on S→A→B→E at 40 + 0 + 40. A node and its parent at one distance must not lose the flow.
    const ng::Network g{4, {{0, 1, 0, 0.01, 1}, {1, 3, 45, 0, 1}, {0, 2, 45, 0, 1}, {2, 3, 0, 0.01, 1}, {1, 2, 0, 0, 1}}, {}};
    const std::vector<double> O{4000}, D{4000};
    const auto eq = ni::combined(g, {0}, {3}, O, D, 0.05, {1e-10, 1000000});
    std::vector<double> net(4, 0.0);
    for (std::size_t a = 0; a < g.arcs.size(); a++) net[g.arcs[a].tail] -= eq.arc_flow[a], net[g.arcs[a].head] += eq.arc_flow[a];
    s.check("Evans on Braess's network, with an arc of no time: the arc flows carry the 4,000 from S to E", std::fabs(net[0] + 4000) + std::fabs(net[1]) + std::fabs(net[2]) + std::fabs(net[3] - 4000), 1e-9);
    s.check("... at 80 minutes, Braess's equilibrium", eq.costs(0, 0) - 80, 1e-6);
  }
  return s.done();
}
