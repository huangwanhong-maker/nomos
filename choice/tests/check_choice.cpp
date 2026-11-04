// Checks of social choice (profile.hpp, scoring.hpp, condorcet.hpp).
//
// The Condorcet paradox: three voters, a majority cycle, no Condorcet winner, every candidate a
// Schulze winner, three Kemeny rankings (the rotations) of support 5. A profile of 37 voters and
// four candidates (COMAP, For All Practical Purposes, the election of the Math Appreciation
// Society) in which four rules elect four candidates: plurality A, Borda B (scores 42, 69, 67,
// 44), instant runoff D (B, then C eliminated, D beating A 23 to 14) and C, the Condorcet winner,
// elected by Copeland, Schulze and Kemeny–Young, while the plurality winner A is the Condorcet
// loser. Schulze's (2011) example of 45 voters and five candidates: its pairwise matrix, its
// strongest paths and the ranking E > A > C > B > D. On random profiles: Copeland, Schulze and
// Kemeny–Young are Condorcet consistent; the Schulze relation is transitive; Kemeny's rankings
// minimise the total Kendall distance to the ballots; Borda scores are the row sums of the
// pairwise matrix and add up to W k(k − 1)/2, truncated ballots included; weights are repeated
// ballots. Approval by hand.

#include <nomos/choice/condorcet.hpp>
#include <nomos/choice/profile.hpp>
#include <nomos/choice/scoring.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <numeric>
#include <string>
#include <vector>

namespace nc = nomos::choice;

namespace {

calx::random::Stream z(20261024, 0x766f7465u, 0, 0);  // "vote"

// Ballots from strings of letters ("ACBED") with weights.
nc::Profile letters(std::size_t k, const std::vector<std::pair<double, std::string>>& rows) {
  nc::Profile p;
  p.candidates = k;
  for (const auto& [w, s] : rows) {
    nc::Ballot b;
    b.weight = w;
    for (char c : s) b.ranking.push_back(static_cast<std::size_t>(c - 'A'));
    p.ballots.push_back(b);
  }
  return p;
}

nc::Profile random_profile(std::size_t k, std::size_t voters, bool truncate) {
  nc::Profile p;
  p.candidates = k;
  for (std::size_t v = 0; v < voters; v++) {
    nc::Ballot b;
    b.ranking.resize(k);
    std::iota(b.ranking.begin(), b.ranking.end(), std::size_t{0});
    for (std::size_t i = k; i > 1; i--) std::swap(b.ranking[i - 1], b.ranking[z.below(i)]);
    if (truncate) b.ranking.resize(z.below(k + 1));
    b.weight = static_cast<double>(1 + z.below(3));
    p.ballots.push_back(b);
  }
  return p;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // The Condorcet paradox.
  {
    const auto p = letters(3, {{1, "ABC"}, {1, "BCA"}, {1, "CAB"}});
    const auto d = nc::pairwise(p);
    s.require("the Condorcet paradox: A beats B, B beats C and C beats A, each 2 to 1",
              d(0, 1) == 2 && d(1, 2) == 2 && d(2, 0) == 2 && d(1, 0) == 1 && d(2, 1) == 1 && d(0, 2) == 1);
    s.require("... no Condorcet winner, no Condorcet loser", !nc::condorcet_winner(d) && !nc::condorcet_loser(d));
    s.require("... Copeland and Borda tie all three, and all three are Schulze winners",
              nc::winners(nc::copeland(d)).size() == 3 && nc::winners(nc::borda(p)).size() == 3 && nc::schulze(d).winners.size() == 3);
    const auto k = nc::kemeny_young(d);
    s.require("... three Kemeny rankings, the rotations ABC, BCA, CAB, of support 5",
              k.rankings.size() == 3 && k.support == 5 && k.rankings[0] == std::vector<std::size_t>{0, 1, 2} && k.rankings[1] == std::vector<std::size_t>{1, 2, 0} &&
                  k.rankings[2] == std::vector<std::size_t>{2, 0, 1});
  }
  // Four rules, four winners (COMAP's Math Appreciation Society).
  {
    const auto p = letters(4, {{14, "ABCD"}, {10, "CBDA"}, {8, "DCBA"}, {4, "BDCA"}, {1, "CDBA"}});
    const auto d = nc::pairwise(p);
    const auto plu = nc::plurality(p), bor = nc::borda(p);
    s.require("37 voters, four candidates: plurality elects A (14 first preferences)", nc::winners(plu) == std::vector<std::size_t>{0});
    s.require("... Borda elects B, with scores 42, 69, 67, 44", bor == std::vector<double>{42, 69, 67, 44});
    const auto irv = nc::instant_runoff(p);
    s.require("... instant runoff elects D: B eliminated, then C, and D beats A 23 to 14",
              irv.winner == 3 && irv.eliminated == std::vector<std::size_t>{1, 2} && irv.rounds[2][3] == 23 && irv.rounds[2][0] == 14 && !irv.tie_broken);
    const auto cw = nc::condorcet_winner(d);
    s.require("... C is the Condorcet winner (beating A 23–14, B 19–18, D 25–12)",
              cw && *cw == 2 && d(2, 0) == 23 && d(2, 1) == 19 && d(1, 2) == 18 && d(2, 3) == 25);
    s.require("... and Copeland, Schulze and Kemeny–Young elect C",
              nc::winners(nc::copeland(d)) == std::vector<std::size_t>{2} && nc::schulze(d).winners == std::vector<std::size_t>{2} &&
                  nc::kemeny_young(d).rankings.size() == 1 && nc::kemeny_young(d).rankings[0][0] == 2);
    const auto cl = nc::condorcet_loser(d);
    s.require("... while the plurality winner A is the Condorcet loser", cl && *cl == 0);
  }
  // Schulze's example: 45 voters, five candidates A–E.
  {
    const auto p = letters(5, {{5, "ACBED"}, {5, "ADECB"}, {8, "BEDAC"}, {3, "CABED"}, {7, "CAEBD"}, {2, "CBADE"}, {7, "DCEBA"}, {8, "EBADC"}});
    const auto d = nc::pairwise(p);
    const std::vector<std::vector<double>> dref{{0, 20, 26, 30, 22}, {25, 0, 16, 33, 18}, {19, 29, 0, 17, 24}, {15, 12, 28, 0, 14}, {23, 27, 21, 31, 0}};
    const std::vector<std::vector<double>> pref{{0, 28, 28, 30, 24}, {25, 0, 28, 33, 24}, {25, 29, 0, 29, 24}, {25, 28, 28, 0, 24}, {25, 28, 28, 31, 0}};
    const auto r = nc::schulze(d);
    bool dm = true, pm = true;
    for (std::size_t a = 0; a < 5; a++)
      for (std::size_t b = 0; b < 5; b++) {
        dm = dm && d(a, b) == dref[a][b];
        pm = pm && r.strength(a, b) == pref[a][b];
      }
    s.require("Schulze's example of 45 voters: the pairwise matrix", dm);
    s.require("... the strongest paths", pm);
    s.require("... the winner E and the ranking E > A > C > B > D", r.winners == std::vector<std::size_t>{4} && r.ranking == std::vector<std::size_t>{4, 0, 2, 1, 3});
  }
  // Random profiles.
  {
    bool consistent = true, transitive = true, kendall = true, borda_rows = true, borda_sum = true, weights = true;
    std::size_t with_winner = 0;
    for (int t = 0; t < 200; t++) {
      const std::size_t k = 3 + z.below(3);
      const bool truncate = t % 2 == 1;
      const auto p = random_profile(k, 3 + z.below(9), truncate);
      const auto d = nc::pairwise(p);
      const auto sch = nc::schulze(d);
      const auto kem = nc::kemeny_young(d);
      if (const auto cw = nc::condorcet_winner(d)) {
        with_winner++;
        consistent = consistent && nc::winners(nc::copeland(d)) == std::vector<std::size_t>{*cw} && sch.winners == std::vector<std::size_t>{*cw};
        for (const auto& rk : kem.rankings) consistent = consistent && rk[0] == *cw;
      }
      // Transitivity of a ≻ b ⟺ p(a, b) > p(b, a).
      const auto& q = sch.strength;
      for (std::size_t a = 0; a < k; a++)
        for (std::size_t b = 0; b < k; b++)
          for (std::size_t c = 0; c < k; c++)
            if (q(a, b) > q(b, a) && q(b, c) > q(c, b)) transitive = transitive && q(a, c) > q(c, a);
      // Kemeny: the least total Kendall distance, computed from the ballots, over every ranking.
      if (!truncate) {
        auto distance = [&](const std::vector<std::size_t>& rk) {
          std::vector<std::size_t> at(k);
          for (std::size_t i = 0; i < k; i++) at[rk[i]] = i;
          double total = 0;
          for (const auto& b : p.ballots) {
            const auto pos = nc::detail::positions(b, k);
            for (std::size_t x = 0; x < k; x++)
              for (std::size_t y = 0; y < k; y++)
                if (pos[x] < pos[y] && at[x] > at[y]) total += b.weight;
          }
          return total;
        };
        std::vector<std::size_t> rk(k);
        std::iota(rk.begin(), rk.end(), std::size_t{0});
        double least = 1e300;
        do least = std::min(least, distance(rk));
        while (std::next_permutation(rk.begin(), rk.end()));
        for (const auto& r : kem.rankings) kendall = kendall && distance(r) == least;
        const auto b = nc::borda(p);
        for (std::size_t a = 0; a < k; a++) {
          double row = 0;
          for (std::size_t c = 0; c < k; c++) row += d(a, c);
          borda_rows = borda_rows && b[a] == row;
        }
      }
      const auto b = nc::borda(p);
      borda_sum = borda_sum && std::fabs(std::accumulate(b.begin(), b.end(), 0.0) - nc::total_weight(p) * static_cast<double>(k * (k - 1)) / 2) < 1e-9;
      // Weights as repeated ballots.
      nc::Profile rep;
      rep.candidates = k;
      for (const auto& bl : p.ballots)
        for (int c = 0; c < static_cast<int>(bl.weight); c++) rep.ballots.push_back({bl.ranking, 1});
      weights = weights && nc::borda(rep) == nc::borda(p) && nc::plurality(rep) == nc::plurality(p) && nc::instant_runoff(rep).winner == nc::instant_runoff(p).winner &&
                nc::pairwise(rep).a == d.a;
    }
    s.require("200 random profiles: Copeland, Schulze and Kemeny–Young elect the Condorcet winner when there is one", consistent,
              std::to_string(with_winner) + " with a Condorcet winner");
    s.require("... the Schulze relation is transitive", transitive);
    s.require("... Kemeny's rankings minimise the total Kendall distance to the ballots", kendall);
    s.require("... Borda scores are the row sums of the pairwise matrix (complete ballots)", borda_rows);
    s.require("... and sum to W k(k − 1)/2, truncated ballots included", borda_sum);
    s.require("... weights count as repeated ballots", weights);
  }
  // Approval.
  {
    const auto a = nc::approval(4, {{0, 2}, {1}, {0, 1, 3}, {2, 0}}, std::vector<double>{1, 2, 1, 3});
    s.require("approval with weights 1, 2, 1, 3: scores 5, 3, 4, 1", a == std::vector<double>{5, 3, 4, 1});
  }
  return s.done();
}
