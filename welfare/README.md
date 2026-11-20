# welfare

Status: **draft** (2026-11-20). Module `nomos::welfare` ([decision 0015](../development/decisions/0015-allocation-and-welfare-orders.md)).

## 1. Concepts

Principles of distribution order whole allocations, and they disagree. No complete aggregation of them escapes Fleurbaey and Tungodden's two tyrannies: tiny gains to the worst off outweighing everything, or small gains to many outweighing a large loss to the worst off. So this module does not add principles into one score. It does four things:

- **evaluates** each principle on each option;
- **reports the undominated set**: the options no other option beats under every principle (Sen's maximal sets);
- **reports the pivot of a pair**: which principles favour which option;
- **finds where a principle's parameter flips a comparison**, so that choosing the parameter is seen as the choice it is.

Dominance relations (Lorenz, generalised Lorenz, sequential for groups ordered by need) give partial orders that hold whatever the parameter. Only when a person declares a weighted aggregate are scores added. Even then the rank acceptabilities (the appraisal module's SMAA-2) are reported beside the undominated set.

The principles:
- **total**: the utilitarian sum;
- **prioritarian(ε)**: Atkinson's equally distributed equivalent;
- **leximin**;
- **shortfall(τ, α)**: the Foster–Greer–Thorbecke gaps below a threshold, sufficientarianism's measure;
- **opportunity**: the least mean over declared types of circumstance, a simplification of Roemer's equality of opportunity, said so.

## 2. Content

- `orders.hpp`: `Outcome`, `Kind`, `Principle`, `evaluate`, `compare`; `generalised_lorenz_at`, `Relation`, `generalised_lorenz`, `lorenz`, `sequential`.
- `comparison.hpp`: `Pairwise`, `Comparison`, `compare_all` (values, pairwise relations, the maximal set, the best under each principle), `pivot`, `flips`, `acceptability`.

*Caveat.*
- The opportunity principle compares type means, not tranches of effort.
- The prioritarian value with aversion above zero needs outcomes above zero.
- Flips are found on a grid, so two flips closer than its step are missed.

## 3. Checks

13 (`tests/check_welfare.cpp`):
- **Principles by hand:** the prioritarian value from the mean towards the minimum, with the geometric mean at ε = 1; the shortfall's three measures; leximin; the opportunity principle.
- **Dominance:** generalised Lorenz dominance implying a higher prioritarian value at every ε on random pairs; a progressive transfer giving Lorenz dominance; crossing curves unordered; sequential dominance by hand.
- **Comparison:** a dominated option left out of the maximal set; the pivot in words; the aversion at which (1, 9) and (4, 5) exchange places; rank acceptabilities summing to one, with a dominated option never first.

## 4. Dependencies

`nomos::measurement` (the mean, Atkinson's index); `nomos::appraisal` (SMAA-2); calx's `random`.

## 5. References

Atkinson (1970); Shorrocks (1983); Foster, Greer and Thorbecke (1984); Atkinson and Bourguignon (1987); Sen (1997); Roemer (1998); Fleurbaey and Tungodden (2010); Adler (2012); Lahdelma and Salminen (2001).
