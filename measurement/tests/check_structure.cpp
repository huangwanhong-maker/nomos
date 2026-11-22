// Checks of nomos/measurement/structure: shift-share worked by hand (a region of 30 growing to 33 in a nation of 400 growing to 430:
// national share 2.25, mix 0.75, competitive 0, Esteban-Marquillas's competitive effect −0.75 and allocation 0.75), its identities on
// random regions, an industry new to the region; related and unrelated variety summing to the entropy, their values on uniform shares
// and in base 2; proximity and density worked by hand; the complexity index's raw form orthogonal to diversity and its affine relation
// with the mean index of a region's activities; fitness and complexity respecting nestedness; refusals.

#include <nomos/measurement/structure.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <string>
#include <vector>

using namespace nomos::measurement;

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
double correlation(const std::vector<double>& a, const std::vector<double>& b) {
  double ma = 0, mb = 0, sab = 0, saa = 0, sbb = 0;
  for (std::size_t i = 0; i < a.size(); i++) ma += a[i] / a.size(), mb += b[i] / b.size();
  for (std::size_t i = 0; i < a.size(); i++) sab += (a[i] - ma) * (b[i] - mb), saa += (a[i] - ma) * (a[i] - ma), sbb += (b[i] - mb) * (b[i] - mb);
  return sab / std::sqrt(saa * sbb);
}
}  // namespace

int main() {
  calx::verification::Suite s;
  {
    const auto S = shift_share({10, 20}, {15, 18}, {100, 300}, {130, 300});
    s.check("shift-share by hand: national share 0.75 + 1.5, mix 2.25 − 1.5, competitive 2 − 2, the change 3",
            std::fabs(S.national_total - 2.25) + std::fabs(S.mix_total - 0.75) + std::fabs(S.competitive_total) + std::fabs(S.change - 3) +
                std::fabs(S.competitive[0] - 2) + std::fabs(S.mix[1] + 1.5),
            1e-12);
    s.check("... Esteban-Marquillas: homothetic 7.5 and 22.5; competitive effect 1.5 and −2.25, allocation 0.5 and 0.25",
            std::fabs(S.competitive_effect[0] - 1.5) + std::fabs(S.competitive_effect[1] + 2.25) + std::fabs(S.allocation[0] - 0.5) + std::fabs(S.allocation[1] - 0.25), 1e-12);
    calx::random::Stream g(20261123, 11, 12, 13);
    double worst = 0;
    for (int t = 0; t < 200; t++) {
      const std::size_t n = 2 + g.below(6);
      std::vector<double> e0(n), e1(n), E0(n), E1(n);
      for (std::size_t i = 0; i < n; i++) e0[i] = g.uniform(1, 50), e1[i] = g.uniform(0, 60), E0[i] = e0[i] + g.uniform(10, 500), E1[i] = g.uniform(5, 700);
      const auto R = shift_share(e0, e1, E0, E1);
      worst = std::max(worst, std::fabs(R.national_total + R.mix_total + R.competitive_total - R.change));
      for (std::size_t i = 0; i < n; i++) worst = std::max(worst, std::fabs(R.competitive_effect[i] + R.allocation[i] - R.competitive[i]));
    }
    s.check("... on 200 random regions: the three terms sum to the change, Esteban-Marquillas's two to the competitive", worst, 1e-9);
    const auto N = shift_share({10, 0}, {12, 5}, {100, 50}, {110, 60});
    s.require("... an industry new to the region: its competitive term all its employment, Esteban-Marquillas's split undefined",
              std::fabs(N.competitive[1] - 5) < 1e-12 && !N.homothetic && std::isnan(N.allocation[0]));
  }
  {
    const auto U = variety({1, 1, 1, 1, 1, 1}, {"a", "a", "a", "b", "b", "b"});
    s.check("variety: six equal classes in two groups: entropy log 6, unrelated log 2, related log 3",
            std::fabs(U.entropy - std::log(6.0)) + std::fabs(U.unrelated - std::log(2.0)) + std::fabs(U.related - std::log(3.0)), 1e-12);
    const auto B = variety({2, 2, 2, 2}, {"x", "x", "y", "y"}, 2);
    s.check("... in base 2, four equal classes in two groups: 2 = 1 + 1", std::fabs(B.entropy - 2) + std::fabs(B.unrelated - 1) + std::fabs(B.related - 1), 1e-12);
    const auto O = variety({3, 1, 4}, {"g", "g", "g"});
    s.require("... one group: no unrelated variety, all of it related", std::fabs(O.unrelated) < 1e-15 && std::fabs(O.related - O.entropy) < 1e-12);
    calx::random::Stream g(20261123, 14, 15, 16);
    double worst = 0;
    for (int t = 0; t < 200; t++) {
      const std::size_t n = 2 + g.below(10);
      std::vector<double> a(n);
      std::vector<std::string> gr(n);
      for (std::size_t i = 0; i < n; i++) a[i] = g.below(4) ? g.uniform(0, 10) : 0, gr[i] = std::string(1, static_cast<char>('a' + g.below(3)));
      a[0] += 1;
      const auto V = variety(a, gr);
      worst = std::max(worst, std::fabs(V.unrelated + V.related - V.entropy));
    }
    s.check("... on 200 random structures, entropy = unrelated + related (the chain rule)", worst, 1e-12);
  }
  {
    calx::linalg::Rectangular M(4, 3);
    const double m[4][3] = {{1, 1, 0}, {1, 0, 1}, {1, 1, 1}, {0, 1, 0}};
    for (int r = 0; r < 4; r++)
      for (int p = 0; p < 3; p++) M(r, p) = m[r][p];
    const auto I = from_binary(M);
    const auto phi = proximity(I);
    s.check("proximity by hand: φ₀₁ = 2/3, φ₀₂ = 2/3, φ₁₂ = 1/3, symmetric, φ_pp = 1",
            std::fabs(phi(0, 1) - 2.0 / 3) + std::fabs(phi(0, 2) - 2.0 / 3) + std::fabs(phi(1, 2) - 1.0 / 3) + std::fabs(phi(2, 1) - phi(1, 2)) + std::fabs(phi(1, 1) - 1), 1e-15);
    const auto w = density(I, phi);
    s.check("... the density of the fourth region around the first activity: (2/3)/(4/3) = 1/2", std::fabs(w(3, 0) - 0.5), 1e-15);
  }
  {
    calx::random::Stream g(20261123, 17, 18, 19);
    double orth = 0, affine = 0;
    int done = 0;
    for (int t = 0; t < 60; t++) {
      const std::size_t R = 5 + g.below(10), P = 5 + g.below(12);
      calx::linalg::Rectangular x(R, P);
      for (auto& v : x.a) v = g.below(3) ? g.uniform(0, 100) : 0;
      auto I = incidence(x);
      bool ok = true;
      for (double d : I.diversity) ok = ok && d > 0;
      for (double u : I.ubiquity) ok = ok && u > 0;
      if (!ok) continue;
      Complexity C;
      try {
        C = complexity(I);
      } catch (const std::invalid_argument&) {
        continue;  // fallen apart into blocks
      }
      double dot = 0, norm = 0;
      for (std::size_t r = 0; r < R; r++) dot += I.diversity[r] * C.regions_raw[r], norm += std::fabs(I.diversity[r] * C.regions_raw[r]);
      orth = std::max(orth, std::fabs(dot) / norm);
      std::vector<double> mean_act(R, 0.0);
      for (std::size_t r = 0; r < R; r++) {
        for (std::size_t p = 0; p < P; p++) mean_act[r] += I.M(r, p) * C.activities[p];
        mean_act[r] /= I.diversity[r];
      }
      affine = std::max(affine, 1 - correlation(C.regions, mean_act));
      done++;
    }
    s.check("economic complexity on random tables: the raw index of regions has zero dot product with diversity", orth, 1e-9);
    s.check("... and the regions' index is an increasing affine function of the mean index of their activities (Hidalgo and Hausmann; Mealy et al.)", affine, 1e-9);
    s.require("... on enough tables", done >= 30, std::to_string(done));
  }
  {
    // A nested incidence: region r holds the first r + 1 activities, so that diversity and ubiquity are nested.
    const std::size_t R = 7, P = 7;
    calx::linalg::Rectangular M(R, P);
    for (std::size_t r = 0; r < R; r++)
      for (std::size_t p = 0; p <= r; p++) M(r, p) = 1;
    const auto F = fitness_complexity(from_binary(M), 30);
    bool nested = true;
    double mf = 0, mq = 0;
    for (std::size_t r = 1; r < R; r++) nested = nested && F.fitness[r] >= F.fitness[r - 1];
    for (std::size_t p = 1; p < P; p++) nested = nested && F.complexity[p] >= F.complexity[p - 1];
    for (double v : F.fitness) mf += v / R;
    for (double v : F.complexity) mq += v / P;
    s.require("fitness and complexity on a nested incidence: fitness rises with diversity, complexity with rarity, each with mean one",
              nested && std::fabs(mf - 1) < 1e-12 && std::fabs(mq - 1) < 1e-12);
  }
  {
    s.require("shift-share with an industry the reference does not hold is refused", refused([] { shift_share({1, 2}, {1, 2}, {0, 5}, {1, 5}); }));
    s.require("variety with no total is refused", refused([] { variety({0, 0}, {"a", "b"}); }));
    calx::linalg::Rectangular M(3, 3);
    M(0, 0) = 1, M(1, 1) = 1;
    s.require("complexity with a region holding nothing is refused", refused([&] { complexity(from_binary(M)); }));
    calx::linalg::Rectangular B(4, 4);
    B(0, 0) = B(1, 0) = B(1, 1) = 1, B(2, 2) = B(3, 2) = B(3, 3) = 1;
    s.require("complexity on an incidence fallen apart into two blocks is refused", refused([&] { complexity(from_binary(B)); }));
    M(0, 1) = 0.5;
    s.require("an incidence that is not of zeros and ones is refused", refused([&] { from_binary(M); }));
  }
  return s.done();
}
