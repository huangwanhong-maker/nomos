// Who is served first: priorities, reserve categories and the order in which they are processed.
//
// A reserve system (Pathak, Sönmez, Ünver and Yenmez, "Fair allocation of vaccines, ventilators and antiviral treatments:
// leaving no ethical value behind in health care rationing", arXiv:2008.00374, read in full by the infrastructure's survey
// *Allocative justice under scarcity*, Part II §2) divides the units of a scarce good into categories, each with its own
// capacity, the individuals eligible for it and its own priority order over them: an open category ordered by one principle
// (say, need), reserves for groups ordered by another (essential workers, the worse off), each reserve a value the rationing
// is meant to honour. A matching gives each individual at most one category's unit. Its three axioms:
//
//   eligibility          an individual holds a unit only of a category it is eligible for;
//   non-wastefulness     a unit is left unused only if no unmatched individual is eligible for its category;
//   respect of priorities an individual holds a unit of a category while another, eligible for it and higher in its
//                         priority, holds none.
//
// The paper's Theorem 1: a matching satisfies the three if and only if it is a cutoff equilibrium; its Proposition 1: these
// are the outcomes of deferred acceptance for some artificial ordering of the categories. The sequential reserve matching
// here processes the categories in a declared precedence order, each taking its highest-priority eligible individuals not
// yet served, and satisfies the three axioms; which matching it gives depends on the order, and the paper's Example 1 (six
// units, seven patients) shows that processing a reserve earlier can help its beneficiaries or not. So the order is a
// declared parameter whose effect is computed: `every_order` gives the beneficiaries under every precedence order.
#pragma once

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace nomos::allocation {

struct Category {
  std::string name;
  std::size_t capacity = 0;
  std::vector<std::size_t> priority;   // the eligible individuals, highest priority first; the others are ineligible
};

inline constexpr std::size_t unserved = static_cast<std::size_t>(-1);

struct ReserveMatching {
  std::vector<std::size_t> category_of;   // for each individual, the category whose unit it holds, or `unserved`
  std::vector<std::vector<std::size_t>> holders;  // for each category, its holders in the order they were served
  std::vector<std::size_t> served_order;  // the individuals in the order they were served
};

// The categories processed in the order given (indices into `categories`), each taking its highest-priority eligible
// individuals not yet served, up to its capacity.
inline ReserveMatching sequential_reserve(std::size_t individuals, const std::vector<Category>& categories, const std::vector<std::size_t>& precedence) {
  for (std::size_t c : precedence)
    if (c >= categories.size()) throw std::invalid_argument("nomos::allocation::sequential_reserve: a precedence order of the categories");
  ReserveMatching m;
  m.category_of.assign(individuals, unserved);
  m.holders.assign(categories.size(), {});
  for (std::size_t c : precedence) {
    for (std::size_t i : categories[c].priority) {
      if (m.holders[c].size() >= categories[c].capacity) break;
      if (i >= individuals) throw std::invalid_argument("nomos::allocation::sequential_reserve: a priority names an individual that is not one");
      if (m.category_of[i] != unserved) continue;
      m.category_of[i] = c;
      m.holders[c].push_back(i);
      m.served_order.push_back(i);
    }
  }
  return m;
}

struct AxiomReport {
  bool eligible = true, non_wasteful = true, respects_priorities = true;
  std::vector<std::string> violations;
};

inline AxiomReport check_axioms(std::size_t individuals, const std::vector<Category>& categories, const std::vector<std::size_t>& category_of) {
  AxiomReport r;
  auto rank = [&](std::size_t c, std::size_t i) {
    const auto& p = categories[c].priority;
    const auto it = std::find(p.begin(), p.end(), i);
    return it == p.end() ? unserved : static_cast<std::size_t>(it - p.begin());
  };
  std::vector<std::size_t> used(categories.size(), 0);
  for (std::size_t i = 0; i < individuals; i++) {
    const std::size_t c = category_of[i];
    if (c == unserved) continue;
    used[c]++;
    if (rank(c, i) == unserved) r.eligible = false, r.violations.push_back("individual " + std::to_string(i) + " holds a unit of " + categories[c].name + " without being eligible");
  }
  for (std::size_t c = 0; c < categories.size(); c++) {
    if (used[c] > categories[c].capacity) r.eligible = false, r.violations.push_back(categories[c].name + " holds more than its capacity");
    if (used[c] < categories[c].capacity)
      for (std::size_t i : categories[c].priority)
        if (category_of[i] == unserved) {
          r.non_wasteful = false;
          r.violations.push_back("a unit of " + categories[c].name + " is unused while individual " + std::to_string(i) + ", eligible for it, is unserved");
          break;
        }
    for (std::size_t i = 0; i < individuals; i++) {
      if (category_of[i] != c) continue;
      for (std::size_t j : categories[c].priority) {
        if (j == i) break;  // only those ranked above i
        if (category_of[j] == unserved) {
          r.respects_priorities = false;
          r.violations.push_back("individual " + std::to_string(i) + " holds a unit of " + categories[c].name + " while " + std::to_string(j) +
                                 ", higher in its priority, holds none");
        }
      }
    }
  }
  return r;
}

struct OrderOutcome {
  std::vector<std::size_t> precedence;
  std::vector<std::size_t> served;   // the individuals served, ascending
};

// The individuals served under every precedence order of the categories (at most 8 categories).
inline std::vector<OrderOutcome> every_order(std::size_t individuals, const std::vector<Category>& categories) {
  if (categories.size() > 8) throw std::invalid_argument("nomos::allocation::every_order: at most 8 categories");
  std::vector<std::size_t> order(categories.size());
  std::iota(order.begin(), order.end(), 0);
  std::vector<OrderOutcome> out;
  do {
    const auto m = sequential_reserve(individuals, categories, order);
    OrderOutcome o{order, {}};
    for (std::size_t i = 0; i < individuals; i++)
      if (m.category_of[i] != unserved) o.served.push_back(i);
    out.push_back(o);
  } while (std::next_permutation(order.begin(), order.end()));
  return out;
}

}  // namespace nomos::allocation
