// Checks of nomos/games/dynamics: population games and their evolutionary dynamics (Sandholm 2015).
//
// Each protocol of Table 13.1 gives, through the mean dynamic of eq. 13.3, the closed-form field of its dynamic, on random games and
// states; the Maynard Smith replicator is the replicator over the mean payoff. Standard Rock–Paper–Scissors: the replicator keeps
// x₁x₂x₃ (closed orbits). Good Rock–Paper–Scissors (strictly contractive): every dynamic reaches (1/3, 1/3, 1/3), the replicator's,
// BNN's and Smith's Lyapunov functions of Table 13.4 never rising. Bad Rock–Paper–Scissors: the replicator's x₁x₂x₃ falls toward the
// boundary (the heteroclinic cycle). Positive correlation at random interior states of random games for the replicator, best response,
// BNN and Smith; Nash stationarity failing for the replicator at a vertex and holding for BNN and Smith. A coordination game, a potential
// game: the potential never falls. Hawk–Dove: the mixed state V/C reached. The feeble twin: under the Smith dynamic in bad
// Rock–Paper–Scissors a strictly dominated twin of scissors is recurrently played by more than a tenth of the population at d = 1/10
// (Hofbauer and Sandholm 2011). The logit's rest point a fixed point of the logit choice. The finite population's paths within 0.02 of
// the mean dynamic at N = 20,000 (Kurtz). Refusals.

#include <nomos/games/dynamics.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

using namespace nomos::games::dynamics;

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
double distance(const State& a, const State& b) {
  double d = 0;
  for (std::size_t i = 0; i < a.size(); i++) d = std::max(d, std::fabs(a[i] - b[i]));
  return d;
}
State random_state(calx::random::Stream& g, std::size_t n) {
  State x(n);
  double s = 0;
  for (auto& v : x) v = 0.05 + g.uniform(), s += v;
  for (auto& v : x) v /= s;
  return x;
}
}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream g(20261123, 7, 8, 9);
  const State centre{1.0 / 3, 1.0 / 3, 1.0 / 3};

  {
    double worst = 0;
    for (int t = 0; t < 200; t++) {
      const std::size_t n = 2 + g.below(4);
      Matrix A(n, std::vector<double>(n));
      for (auto& row : A)
        for (auto& a : row) a = g.uniform(-2, 2);
      const auto F = linear(A);
      const auto x = random_state(g, n);
      for (Protocol p : {Protocol::replicator, Protocol::best_response, Protocol::logit, Protocol::bnn, Protocol::smith}) {
        const auto v = mean_dynamic(protocol(p, F(x), x, 0.3), x), w = field(p, F, x, 0.3);
        worst = std::max(worst, distance(v, w));
      }
    }
    s.check("each protocol of Table 13.1 gives its closed-form dynamic through the mean dynamic, on 200 random games and states", worst, 1e-13);
    Matrix A = rock_paper_scissors(1, 1);
    for (auto& row : A)
      for (auto& a : row) a += 3;
    const auto F = linear(A);
    double gap = 0;
    for (int t = 0; t < 50; t++) {
      const auto x = random_state(g, 3);
      const auto r = field(Protocol::replicator, F, x), m = field(Protocol::maynard_smith, F, x);
      const double fbar = mean_payoff(x, F(x));
      for (std::size_t i = 0; i < 3; i++) gap = std::max(gap, std::fabs(m[i] - r[i] / fbar));
    }
    s.check("the Maynard Smith replicator is the replicator over the mean payoff (positive payoffs)", gap, 1e-15);
  }
  {
    const auto F = linear(rock_paper_scissors(1, 1));
    const auto T = trajectory(Protocol::replicator, F, {0.5, 0.3, 0.2}, 30, 0.005, 100);
    const double p0 = 0.5 * 0.3 * 0.2;
    double drift = 0;
    for (const auto& x : T.x) drift = std::max(drift, std::fabs(x[0] * x[1] * x[2] / p0 - 1));
    s.check("standard Rock–Paper–Scissors: the replicator keeps x₁x₂x₃ over 30 time units (closed orbits)", drift, 1e-7);
    double away = 0;
    for (const auto& x : T.x) away = std::max(away, distance(x, centre));
    s.require("... and does not approach the centre", away > 0.1);
  }
  {
    const auto F = linear(rock_paper_scissors(2, 1));
    const State x0{0.6, 0.3, 0.1};
    bool reached = true;
    for (Protocol p : {Protocol::replicator, Protocol::logit, Protocol::bnn, Protocol::smith}) {
      const auto T = trajectory(p, F, x0, p == Protocol::replicator ? 150 : 60, 0.01, 1000, 0.2);
      reached = reached && distance(T.x.back(), centre) < 1e-4;
    }
    const auto B = trajectory(Protocol::best_response, F, x0, 60, 0.01, 1000);
    s.require("good Rock–Paper–Scissors: the replicator, logit, BNN and Smith dynamics reach the centre", reached);
    s.check("... and best response, a differential inclusion stepped at 0.01, comes within 0.05", distance(B.x.back(), centre), 0.05);
    double rise_h = 0, rise_g = 0, rise_p = 0;
    const auto R = trajectory(Protocol::replicator, F, x0, 40, 0.01), N = trajectory(Protocol::bnn, F, x0, 40, 0.01), S = trajectory(Protocol::smith, F, x0, 40, 0.01);
    for (std::size_t k = 1; k < R.x.size(); k++) rise_h = std::max(rise_h, replicator_entropy(centre, R.x[k]) - replicator_entropy(centre, R.x[k - 1]));
    for (std::size_t k = 1; k < N.x.size(); k++) rise_g = std::max(rise_g, bnn_lyapunov(F, N.x[k]) - bnn_lyapunov(F, N.x[k - 1]));
    for (std::size_t k = 1; k < S.x.size(); k++) rise_p = std::max(rise_p, smith_lyapunov(F, S.x[k]) - smith_lyapunov(F, S.x[k - 1]));
    s.check("... the replicator's H, BNN's Γ and Smith's Ψ (Table 13.4) never rise along their paths", std::max({rise_h, rise_g, rise_p}), 1e-12);
  }
  {
    const auto F = linear(rock_paper_scissors(1, 2));
    const auto T = trajectory(Protocol::replicator, F, {0.4, 0.35, 0.25}, 80, 0.005, 200);
    double lowest = 1, rise = 0;
    for (std::size_t k = 1; k < T.x.size(); k++) {
      const auto& x = T.x[k];
      if (T.t[k] > 40) lowest = std::min(lowest, *std::min_element(x.begin(), x.end()));
      rise = std::max(rise, x[0] * x[1] * x[2] - T.x[k - 1][0] * T.x[k - 1][1] * T.x[k - 1][2]);
    }
    s.require("bad Rock–Paper–Scissors: the replicator's x₁x₂x₃ never rises and a share falls below 10⁻³ (toward the heteroclinic cycle)", rise <= 1e-15 && lowest < 1e-3,
              "lowest " + std::to_string(lowest));
  }
  {
    bool pc = true;
    for (int t = 0; t < 300; t++) {
      const std::size_t n = 2 + g.below(4);
      Matrix A(n, std::vector<double>(n));
      for (auto& row : A)
        for (auto& a : row) a = g.uniform(-2, 2);
      const auto F = linear(A);
      const auto x = random_state(g, n);
      for (Protocol p : {Protocol::replicator, Protocol::best_response, Protocol::bnn, Protocol::smith}) {
        const auto I = incentive(p, F, x);
        pc = pc && (I.rest || I.inner > 0);
      }
    }
    s.require("positive correlation: V(x)′F(x) > 0 wherever V(x) ≠ 0, for the replicator, best response, BNN and Smith, on 300 random games and states", pc);
    const auto F = linear(rock_paper_scissors(1, 1));
    const State vertex{1, 0, 0};
    const auto Ir = incentive(Protocol::replicator, F, vertex), Ib = incentive(Protocol::bnn, F, vertex), Is = incentive(Protocol::smith, F, vertex);
    s.require("Nash stationarity: a vertex of Rock–Paper–Scissors is a rest point of the replicator and no equilibrium; BNN and Smith move from it",
              Ir.rest && !Ir.nash && !Ib.rest && !Is.rest);
    s.require("... the centre is a rest point of all three and the equilibrium",
              incentive(Protocol::replicator, F, centre).rest && incentive(Protocol::bnn, F, centre).rest && incentive(Protocol::smith, F, centre).rest && is_nash(F, centre));
  }
  {
    const Matrix A{{1, 0}, {0, 2}};
    const auto F = linear(A);
    auto f = [&](const State& x) { return 0.5 * (A[0][0] * x[0] * x[0] + A[1][1] * x[1] * x[1]); };
    double fall = 0, bnn_x1 = 0;
    bool to_e2 = true;
    for (Protocol p : {Protocol::replicator, Protocol::bnn, Protocol::smith}) {
      const auto T = trajectory(p, F, {0.6, 0.4}, 40, 0.01);
      for (std::size_t k = 1; k < T.x.size(); k++) fall = std::max(fall, f(T.x[k - 1]) - f(T.x[k]));
      if (p == Protocol::bnn) bnn_x1 = T.x.back()[0];
      else to_e2 = to_e2 && T.x.back()[1] > 0.999;
    }
    s.check("a coordination game, a potential game: the potential never falls under the replicator, BNN and Smith", fall, 1e-13);
    s.require("... from (0.6, 0.4), below the mixed equilibrium's 2/3, the replicator and Smith reach the second strategy", to_e2);
    s.require("... and BNN approaches it algebraically, ẋ₁ ≈ −2x₁² near the vertex: x₁ t between 1/2 and 0.6 at t = 40", bnn_x1 * 40 > 0.5 && bnn_x1 * 40 < 0.6,
              "x1 t = " + std::to_string(bnn_x1 * 40));
    const auto L = trajectory(Protocol::logit, F, {0.6, 0.4}, 60, 0.01, 1000, 0.05);
    const auto c = logit_choice(F(L.x.back()), 0.05);
    s.check("... the logit's rest point at η = 0.05 is a fixed point of the logit choice", distance(L.x.back(), c), 1e-8);
  }
  {
    const auto F = linear({{-1, 2}, {0, 1}});  // Hawk–Dove, V = 2, C = 4
    bool ok = true;
    for (Protocol p : {Protocol::replicator, Protocol::bnn, Protocol::smith}) ok = ok && std::fabs(trajectory(p, F, {0.9, 0.1}, 80, 0.01, 1000).x.back()[0] - 0.5) < 1e-6;
    s.require("Hawk–Dove (V = 2, C = 4): the replicator, BNN and Smith reach hawks at V/C = 1/2", ok);
  }
  {
    const double d = 0.1;
    const auto F = linear({{0, -2, 1, 1}, {1, 0, -2, -2}, {-2, 1, 0, 0}, {-2 - d, 1 - d, -d, -d}});
    const auto T = trajectory(Protocol::smith, F, {0.3, 0.3, 0.2, 0.2}, 300, 0.01, 10);
    double most = 0;
    for (std::size_t k = 0; k < T.x.size(); k++)
      if (T.t[k] > 150) most = std::max(most, T.x[k][3]);
    s.require("the feeble twin: under Smith in bad Rock–Paper–Scissors, a twin of scissors dominated by 1/10 recurrently holds more than a tenth (Hofbauer and Sandholm 2011)",
              most > 0.10, "most " + std::to_string(most));
  }
  {
    const auto F = linear(rock_paper_scissors(2, 1));
    const std::vector<double> times{0.5, 1, 2};
    const auto T = trajectory(Protocol::replicator, F, {0.6, 0.3, 0.1}, 2, 0.001, 500);
    const auto Ts = trajectory(Protocol::smith, F, {0.6, 0.3, 0.1}, 2, 0.001, 500);
    const auto X = finite_population(Protocol::replicator, F, {12000, 6000, 2000}, times, g);
    const auto Y = finite_population(Protocol::smith, F, {12000, 6000, 2000}, times, g);
    double worst = 0;
    const std::size_t at[] = {1, 2, 4};  // the strides of 0.5 at t = 0.5, 1 and 2
    for (std::size_t k = 0; k < 3; k++) worst = std::max({worst, distance(X[k], T.x[at[k]]), distance(Y[k], Ts.x[at[k]])});
    s.check("the finite population of 20,000 under the replicator's and Smith's protocols stays within 0.02 of the mean dynamic (Kurtz)", worst, 0.02);
  }
  {
    const auto F = linear(rock_paper_scissors(1, 1));
    s.require("a state whose shares do not sum to one is refused", refused([&] { trajectory(Protocol::smith, F, {0.5, 0.5, 0.5}, 1, 0.1); }));
    s.require("a logit without noise is refused", refused([&] { logit_choice({1, 2}, 0); }));
    s.require("the Maynard Smith replicator has no protocol", refused([&] { protocol(Protocol::maynard_smith, {1, 1, 1}, centre); }));
    s.require("a finite population with a negative count is refused", refused([&] { finite_population(Protocol::smith, F, {1, -1, 2}, {1}, g); }));
  }
  return s.done();
}
