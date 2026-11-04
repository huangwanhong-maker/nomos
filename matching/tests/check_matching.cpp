// Checks of two-sided matching and housing markets (deferred_acceptance.hpp, stability.hpp,
// top_trading_cycles.hpp), against Gale and Shapley's example and brute force.
//
// Gale and Shapley's (1962) three men and three women with cyclic preferences: the men-proposing
// algorithm gives every man his first choice, the women-proposing every woman hers, and exactly
// three matchings are stable, the third giving everyone their second choice; a woman gains by
// truncating her list (the receiving side may manipulate; Dubins and Freedman 1981). On random
// one-to-one markets with incomplete lists, against every stable matching found by enumerating
// all matchings: the proposers' algorithm is stable, best for every proposer and worst for every
// receiver; the same agents are matched in every stable matching (McVitie and Wilson 1970); no
// proposer gains by any false list (Dubins and Freedman 1981; Roth 1982a). On random many-to-one
// markets, against every stable assignment: both algorithms stable, applicant-optimal and
// applicant-pessimal; every institution fills the same number of places in each (Roth 1986); and
// replacing each institution by copies, one per place, gives the same assignment. Top trading
// cycles: a four-agent example by hand; on random markets, the allocation in the core and the
// only one there (Roth and Postlewaite 1977), and no agent gains by a false list (Roth 1982b).

#include <nomos/matching/deferred_acceptance.hpp>
#include <nomos/matching/stability.hpp>
#include <nomos/matching/top_trading_cycles.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <functional>
#include <numeric>
#include <string>
#include <vector>

namespace nm = nomos::matching;
using nm::unmatched;

namespace {

calx::random::Stream z(20261024, 0x6d617463u, 0, 0);  // "matc"

// A random list: each of `n` partners acceptable with probability `p`, in random order.
std::vector<std::size_t> random_list(std::size_t n, double p) {
  std::vector<std::size_t> l;
  for (std::size_t k = 0; k < n; k++)
    if (z.uniform() < p) l.push_back(k);
  for (std::size_t k = l.size(); k > 1; k--) std::swap(l[k - 1], l[z.below(k)]);
  return l;
}

nm::Market random_market(std::size_t na, std::size_t nc, double p, std::vector<std::size_t> capacity) {
  nm::Market m;
  for (std::size_t a = 0; a < na; a++) m.applicants.push_back(random_list(nc, p));
  for (std::size_t c = 0; c < nc; c++) m.institutions.push_back(random_list(na, p));
  m.capacity = std::move(capacity);
  return m;
}

// Every stable assignment, by enumerating each applicant's institution (or none).
std::vector<nm::Assignment> all_stable(const nm::Market& m) {
  const std::size_t na = m.applicants.size(), nc = m.institutions.size();
  std::vector<nm::Assignment> out;
  std::vector<std::size_t> choice(na, 0);  // 0 … nc, nc meaning none
  while (true) {
    nm::Assignment x;
    x.institution_of.assign(na, unmatched);
    x.applicants_of.assign(nc, {});
    for (std::size_t a = 0; a < na; a++)
      if (choice[a] < nc) x.institution_of[a] = choice[a], x.applicants_of[choice[a]].push_back(a);
    if (nm::is_stable(m, x)) out.push_back(x);
    std::size_t a = 0;
    while (a < na && ++choice[a] > nc) choice[a++] = 0;
    if (a == na) break;
  }
  return out;
}

// The position of `c` on `list` (unmatched, or not listed, ranks below everything listed).
std::size_t rank(const std::vector<std::size_t>& list, std::size_t c) {
  const auto it = std::find(list.begin(), list.end(), c);
  return c == unmatched || it == list.end() ? list.size() + 1 : static_cast<std::size_t>(it - list.begin());
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // Gale and Shapley's example: men α, β, γ (0, 1, 2) and women A, B, C (0, 1, 2).
  {
    const nm::Preferences men{{0, 1, 2}, {1, 2, 0}, {2, 0, 1}};
    const nm::Preferences women{{1, 2, 0}, {2, 0, 1}, {0, 1, 2}};
    const auto by_men = nm::stable_matching(men, women);
    const auto by_women = nm::stable_matching(women, men);
    s.require("Gale and Shapley's example: men proposing, every man gets his first choice (αA, βB, γC)",
              by_men.partner_of_proposer == std::vector<std::size_t>{0, 1, 2});
    s.require("... women proposing, every woman gets hers (βA, γB, αC)", by_women.partner_of_proposer == std::vector<std::size_t>{1, 2, 0});
    const auto stable = all_stable({men, women, {}});
    bool second = false;
    for (const auto& x : stable) second = second || x.institution_of == std::vector<std::size_t>{1, 2, 0};
    s.require("... exactly three matchings are stable, among them everyone's second choice (αB, βC, γA)", stable.size() == 3 && second,
              std::to_string(stable.size()) + " stable");
    // Woman A truncates her list [β, γ, α] to [β, γ] when the men propose.
    nm::Preferences lie = women;
    lie[0] = {1, 2};
    const auto after = nm::stable_matching(men, lie);
    s.require("... and woman A, matched to her last choice α, gets γ by declaring α unacceptable", after.partner_of_receiver[0] == 2);
  }
  // Random one-to-one markets.
  {
    bool stable_ok = true, proposer_best = true, receiver_worst = true, same_matched = true, receivers_best = true;
    std::size_t count = 0;
    for (int t = 0; t < 60; t++) {
      const std::size_t np = 3 + z.below(3), nr = 3 + z.below(3);
      const auto m = random_market(np, nr, 0.8, {});
      const auto da = nm::deferred_acceptance(m, nm::Proposing::applicants);
      const auto rda = nm::deferred_acceptance(m, nm::Proposing::institutions);
      stable_ok = stable_ok && nm::is_stable(m, da) && nm::is_stable(m, rda);
      const auto stable = all_stable(m);
      count += stable.size();
      for (const auto& x : stable) {
        for (std::size_t a = 0; a < np; a++) {
          proposer_best = proposer_best && rank(m.applicants[a], da.institution_of[a]) <= rank(m.applicants[a], x.institution_of[a]);
          same_matched = same_matched && (da.institution_of[a] == unmatched) == (x.institution_of[a] == unmatched);
        }
        for (std::size_t c = 0; c < nr; c++) {
          const std::size_t mine = da.applicants_of[c].empty() ? unmatched : da.applicants_of[c][0];
          const std::size_t theirs = x.applicants_of[c].empty() ? unmatched : x.applicants_of[c][0];
          const std::size_t best = rda.applicants_of[c].empty() ? unmatched : rda.applicants_of[c][0];
          receiver_worst = receiver_worst && rank(m.institutions[c], mine) >= rank(m.institutions[c], theirs);
          receivers_best = receivers_best && rank(m.institutions[c], best) <= rank(m.institutions[c], theirs);
          same_matched = same_matched && (mine == unmatched) == (theirs == unmatched);
        }
      }
    }
    s.require("sixty random one-to-one markets (3–5 a side, incomplete lists): both algorithms give stable matchings", stable_ok,
              std::to_string(count) + " stable matchings enumerated");
    s.require("... the proposers' is the best stable matching for every proposer", proposer_best);
    s.require("... and the worst for every receiver", receiver_worst);
    s.require("... the receivers' is the best for every receiver", receivers_best);
    s.require("... the same agents are matched in every stable matching", same_matched);
  }
  // Truthful proposing is dominant: every false list of one proposer, on random 3 × 3 markets.
  {
    bool truthful = true;
    std::size_t tried = 0;
    for (int t = 0; t < 30; t++) {
      const auto m = random_market(3, 3, 1.0, {});
      const auto truth = nm::deferred_acceptance(m);
      for (std::size_t a = 0; a < 3; a++) {
        // Every ordered list of a subset of {0, 1, 2}.
        for (std::size_t mask = 0; mask < 8; mask++) {
          std::vector<std::size_t> sub;
          for (std::size_t c = 0; c < 3; c++)
            if (mask >> c & 1) sub.push_back(c);
          do {
            auto lie = m;
            lie.applicants[a] = sub;
            const auto r = nm::deferred_acceptance(lie);
            tried++;
            truthful = truthful && rank(m.applicants[a], r.institution_of[a]) >= rank(m.applicants[a], truth.institution_of[a]);
          } while (std::next_permutation(sub.begin(), sub.end()));
        }
      }
    }
    s.require("no proposer gains by a false list (30 markets, every list of every proposer)", truthful, std::to_string(tried) + " lists");
  }
  // Random many-to-one markets.
  {
    bool stable_ok = true, optimal = true, pessimal = true, places = true, copies = true;
    std::size_t count = 0;
    for (int t = 0; t < 40; t++) {
      const std::vector<std::size_t> cap{2, 1, 2};
      const auto m = random_market(6, 3, 0.75, cap);
      const auto ap = nm::deferred_acceptance(m, nm::Proposing::applicants);
      const auto ip = nm::deferred_acceptance(m, nm::Proposing::institutions);
      stable_ok = stable_ok && nm::is_stable(m, ap) && nm::is_stable(m, ip);
      const auto stable = all_stable(m);
      count += stable.size();
      for (const auto& x : stable) {
        for (std::size_t a = 0; a < 6; a++) {
          optimal = optimal && rank(m.applicants[a], ap.institution_of[a]) <= rank(m.applicants[a], x.institution_of[a]);
          pessimal = pessimal && rank(m.applicants[a], ip.institution_of[a]) >= rank(m.applicants[a], x.institution_of[a]);
        }
        for (std::size_t c = 0; c < 3; c++) places = places && x.applicants_of[c].size() == ap.applicants_of[c].size();
      }
      // Copies: institution c becomes cap[c] institutions with its list; applicants list copies in place of c.
      nm::Market one;
      std::vector<std::size_t> owner;
      std::vector<std::vector<std::size_t>> copy_of(3);
      for (std::size_t c = 0; c < 3; c++)
        for (std::size_t k = 0; k < cap[c]; k++) {
          copy_of[c].push_back(owner.size());
          owner.push_back(c);
          one.institutions.push_back(m.institutions[c]);
        }
      for (std::size_t a = 0; a < 6; a++) {
        std::vector<std::size_t> l;
        for (std::size_t c : m.applicants[a]) l.insert(l.end(), copy_of[c].begin(), copy_of[c].end());
        one.applicants.push_back(l);
      }
      const auto r = nm::deferred_acceptance(one);
      for (std::size_t a = 0; a < 6; a++) copies = copies && (r.institution_of[a] == unmatched ? unmatched : owner[r.institution_of[a]]) == ap.institution_of[a];
    }
    s.require("forty random many-to-one markets (6 applicants, capacities 2, 1, 2): both algorithms stable", stable_ok,
              std::to_string(count) + " stable assignments enumerated");
    s.require("... applicants proposing is best for every applicant among stable assignments", optimal);
    s.require("... institutions proposing is worst for every applicant", pessimal);
    s.require("... every institution fills the same number of places in each (rural hospitals)", places);
    s.require("... and one copy of each institution per place gives the same assignment", copies);
  }
  // Top trading cycles.
  {
    const nm::Preferences p{{2, 1, 0, 3}, {0, 2, 1, 3}, {0, 1, 2, 3}, {3, 0, 1, 2}};
    s.require("top trading cycles by hand: 0 and 2 swap, 3 keeps its own, then 1 keeps its own",
              nm::top_trading_cycles(p) == std::vector<std::size_t>{2, 1, 0, 3});
    // The core, by every coalition and every reallocation of its houses.
    auto in_core = [](const nm::Preferences& pr, const std::vector<std::size_t>& alloc) {
      const std::size_t n = pr.size();
      for (std::size_t mask = 1; mask < (std::size_t{1} << n); mask++) {
        std::vector<std::size_t> S;
        for (std::size_t i = 0; i < n; i++)
          if (mask >> i & 1) S.push_back(i);
        std::vector<std::size_t> houses = S;
        do {
          bool weak = true, strict = false;
          for (std::size_t k = 0; k < S.size(); k++) {
            const std::size_t now = rank(pr[S[k]], alloc[S[k]]), then = rank(pr[S[k]], houses[k]);
            weak = weak && then <= now;
            strict = strict || then < now;
          }
          if (weak && strict) return false;
        } while (std::next_permutation(houses.begin(), houses.end()));
      }
      return true;
    };
    auto random_preferences = [](std::size_t n) {
      nm::Preferences pr(n);
      for (auto& l : pr) {
        l.resize(n);
        std::iota(l.begin(), l.end(), std::size_t{0});
        for (std::size_t k = n; k > 1; k--) std::swap(l[k - 1], l[z.below(k)]);
      }
      return pr;
    };
    bool core = true, unique = true;
    for (int t = 0; t < 25; t++) {
      const auto pr = random_preferences(4);
      const auto ttc = nm::top_trading_cycles(pr);
      core = core && in_core(pr, ttc);
      std::vector<std::size_t> alloc(4);
      std::iota(alloc.begin(), alloc.end(), std::size_t{0});
      do
        if (alloc != ttc) unique = unique && !in_core(pr, alloc);
      while (std::next_permutation(alloc.begin(), alloc.end()));
    }
    s.require("twenty-five random four-agent markets: the allocation is in the core", core);
    s.require("... and is the only allocation in it (Roth and Postlewaite 1977)", unique);
    bool truthful = true;
    for (int t = 0; t < 25; t++) {
      const auto pr = random_preferences(4);
      const auto truth = nm::top_trading_cycles(pr);
      for (std::size_t i = 0; i < 4; i++) {
        std::vector<std::size_t> lie(4);
        std::iota(lie.begin(), lie.end(), std::size_t{0});
        do {
          auto q = pr;
          q[i] = lie;
          truthful = truthful && rank(pr[i], nm::top_trading_cycles(q)[i]) >= rank(pr[i], truth[i]);
        } while (std::next_permutation(lie.begin(), lie.end()));
      }
    }
    s.require("... and no agent gains by any of the 24 lists it could declare (Roth 1982b)", truthful);
  }
  return s.done();
}
