# 0015. Allocation, and Principles as Orders over Allocations

- **Status:** Accepted, provisionally (2026-11-20).
- **Date:** 2026-11-20
- **Decided by:** the infrastructure's ADR 0136, which founds its allocation domain at the user's request on allocative justice. The user asked, when resources are limited, how to allocate, to which entity first and then which, what and how much, holding that "it is not an formal optimization".

## 1. Context

`games` held claims problems (proportional, equal awards and losses, the Talmud rule, random arrival) and the Shapley value. `matching` held deferred acceptance and top trading cycles; `measurement` inequality and poverty; `appraisal` multi-criteria analysis with SMAA. The infrastructure's survey *Allocative justice under scarcity* (2026, five parts) found three points:

- **Rules, not an optimum.** Just allocation is decided by rules whose properties are known: rationing rules, priority and reserve systems, sequences and lotteries, proportional budgeting. Principles of distribution are orders over whole allocations that conflict, and a complete aggregation of them falls into one of two tyrannies.
- **What replaces maximisation:** constraints, sufficiency first, priorities with declared precedence, partial orders and their undominated sets, the pivots of a comparison, and lotteries among options on a par.
- **Formulas from practice:** sufficiency first, Porto Alegre's banded points, fiscal equalisation of needs against capacities.

## 2. Decision

- **A module `allocation`** holds the rules:
  - more rationing rules, among them TAL, duality and priority classes;
  - the formulas of practice, each reporting its parts apart;
  - reserve systems with their three axioms checked and every precedence order's beneficiaries;
  - sequential and random assignment, with a random assignment carried out as a lottery;
  - tests of fair division and procedures with guarantees;
  - participatory budgeting with a check of EJR.
- **A module `welfare`** holds principles as orders over allocations:
  - the principles: total, prioritarian, leximin, shortfall, opportunity;
  - the dominance relations: Lorenz, generalised Lorenz, sequential;
  - the comparison of options under all principles at once: the maximal set, the pivots, where a parameter flips a comparison, and SMAA-2 from `appraisal` only when weights are declared.
- **Reuse.** The claims rules of `games`, the inequality measures of `measurement` and the SMAA of `appraisal` are reused, not written again.

## 3. Alternatives Considered

- **Extending `games`.** Rejected: allocation is a subject of its own, and its rules reach beyond cooperative games.
- **A single social welfare function as the module's answer.** Rejected for the reasons the survey gives; a weighted aggregate is available only when declared, and shown beside the undominated set.
- **The smart reserve matching, EFX and MMS algorithms, the bihierarchy decomposition, and Roemer's tranches.** Deferred.

## 4. Consequences

- **Checks:** 50 (37 and 13), held to closed forms, published cases and axioms on random instances. Among them: the Talmud's three estates; the hard reserve that open-first processing leaves idle; a majority that utilitarian greedy gives everything and equal shares does not.
- **The infrastructure's allocation domain** names what is allocated, to whom and by what procedure, and uses these rules.
- **Provisional:** the modules' names; the tolerance of comparisons; the grid of flips.
