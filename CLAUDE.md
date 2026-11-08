# CLAUDE.md

This file guides Claude Code in the nomos repository.

nomos is one of the foundation libraries of [The Science Infrastructures](https://github.com/huangwanhong-maker/The-Science-Infrastructures), which includes it as the git submodule `foundations/nomos`. It is a project in its own right and must make sense without the infrastructure.

## 1. Purpose and Scope

nomos implements the general methods of the social sciences:

- measurement: inequality, poverty, concentration, specialisation, and later scales and item response theory;
- index numbers;
- games: equilibria, dominance, zero-sum games, and later mechanism design;
- matching markets;
- social choice;
- econometrics: discrete choice (marginal effects, the multinomial logit) and fixed-effects pseudo-maximum likelihood; panels moved to tyche, instruments and differences in differences to peira (decision 0007);
- agent-based models: schedules, reproducible draws, event logs, spaces and choices (decision 0004);
- surveys and social networks (planned).

It builds on four sibling libraries:

- **tyche**, for probability, estimation, regression and time series;
- **calx**, for numerical methods: linear algebra, optimisation, random streams;
- **morphe**, for exact and combinatorial structures: graphs, orders;
- **physis**, through tyche.

Out of scope:

| Topic | Where it belongs |
|---|---|
| General probability, estimation, testing, regression (linear and generalised), time series, causal graphs | tyche |
| Floating-point numerical methods and devices | calx |
| Exact algebraic, combinatorial and topological structures | morphe |
| Particular economies, polities, societies and their data (countries, currencies, elections, censuses) | the domains of The Science Infrastructures |

The boundary with tyche is [decision 0002](development/decisions/0002-scope-boundaries-and-dependencies.md): a general statistical model or procedure, used across the sciences, is tyche's; a method whose formulation is that of agents, choices, institutions, welfare and distribution, or whose established use is in the social sciences, is nomos's, built on tyche's.

## 2. Invariants

1. **Dependencies.** nomos depends only on tyche, calx, morphe and physis among project libraries, and never on the infrastructure's domains or applications.
2. **No particular objects.** Nothing in nomos names an object in nature or a particular society. Units, groups, items, players, applicants and institutions are numbered; their meaning is the caller's.
3. **Reuse.** What a sibling library holds is used through its public headers, never written again here. A general method nomos needs that names no social object belongs in the sibling where it fits.
4. **Cited sources.** Every implemented index, rule, estimator or algorithm cites its source; equation, theorem or section numbers are given only where they are certain.
5. **Conventions stated.** Every method states its conventions where they vary in the literature: population or sample divisors, ties, truncated ballots, strictness, degrees of freedom.

## 3. Design Philosophy

nomos shares its governing principle with its sibling repositories: **modular, extensible, revisable**.

- **Status vocabulary:**
  - `planned`: named, no content yet;
  - `seed`: an idea, possibly wrong;
  - `draft`: being worked out;
  - `stable`: tested and relied upon;
  - `deprecated`: superseded, with a pointer to its replacement.
- **Axioms as checks.** A social-scientific method is usually characterised by the properties it satisfies (the principle of transfers, time reversal, stability, Condorcet consistency). Each is checked numerically where it holds, and shown to fail where it is known to fail.
- **Exact where the structure is discrete.** Matchings, profiles and pure equilibria are computed exactly; floating point enters only through payoffs, amounts and estimates.

## 4. Repository Layout

```
nomos/
├── CMakeLists.txt   the project: one target nomos::<module> per module; the siblings added once
├── measurement/ · indices/ · games/ · matching/ · choice/ · econometrics/    the modules written
└── development/     decisions · ideas · release notes
```

Each module directory holds `CMakeLists.txt`, `README.md` (concepts first, then the content, then the checks), `include/nomos/<module>/` and `tests/` (check programs on calx's `verification`, registered with CTest by `calx_add_check`). Module directories are created with their first content; until then a module is a `planned` row in the README. The modules and their order are in [decision 0003](development/decisions/0003-the-first-modules.md).

## 5. Working Procedures

### 5.1 Module Contract

Each module has:

- a `README.md` giving its concepts, content, checks, dependencies, references and open questions;
- public headers in `include/nomos/<module>/`, header-only C++20 for now;
- checks in `tests/`;
- optionally, Python bindings in `bindings/python/`.

### 5.2 Decisions and Ideas

- Decisions are numbered records in `development/decisions/`. Records are never rewritten; a later record supersedes an earlier one.
- Ideas are files named `development/ideas/YYYY-MM-DD-slug.md`, with status `seed`, `promoted` or `abandoned`.
- **Release notes.** At each commit, add to `development/release-notes/YYYY-MM-DD.md` (one file per day, one short section per commit) a concise summary in plain prose of what is new, what one can now do and what works better, not a list of changes (the git log keeps those); not between commits. This convention is shared with the sibling repositories.

### 5.3 Change Discipline

- Before restructuring or changing an invariant (§2), propose the change and wait for the user's agreement.
- Keep each change to one concept.

### 5.4 Shared Conventions

calx, morphe, physis, tyche, nomos and the infrastructure keep these conventions identical:

- writing style;
- status vocabulary;
- decision, idea and release-note conventions.

Change them in every repository of [The Science Infrastructures](https://github.com/huangwanhong-maker/The-Science-Infrastructures) together: the foundation libraries and the infrastructure's repository.

### 5.5 Working Inside the Infrastructure

When you edit nomos inside the infrastructure's `foundations/nomos`:

1. commit in nomos first;
2. then commit the updated pointer in the infrastructure.

When a change needs a sibling to change as well, commit them in dependency order: calx and morphe first, then physis, then tyche, then nomos.

## 6. Conventions

- **Weights.** A weight is a frequency or sampling weight; an index of a weighted distribution is that of the population the weights describe, so integer weights equal replication.
- **Divisors.** Population divisors (W) for indices of a distribution; the degrees of freedom of each estimator are stated where it is defined.
- **Orientation.** Lists of preferences and ballots run best first. A matrix of a two-player game gives the row player's payoffs in `a` and the column player's in `b`.
- **Writing style.**
  - Headings are declarative, numbered and academic, never questions.
  - Concepts come before structure.
  - Label asides with *Rationale.*, *Remark.* or *Caveat.*
  - British spelling; English throughout.

## 7. Technology

- **Language: C++.** This follows ADR 0003 of The Science Infrastructures.
- **Provisional:** C++20, header-only targets `nomos::<module>`; CMake 3.20 or later; checks are programs on calx's `verification`, registered with CTest. On the development machine the compiler is GCC in WSL.
- **Dependency resolution.** Inside the infrastructure, nomos uses the targets of `foundations/calx`, `foundations/morphe`, `foundations/physis` and `foundations/tyche`, adding each that no larger project has added. Never nest submodules.
- **Licences.** nomos is MIT-licensed. Check every dependency before adopting it.

## 8. Open Decisions

1. **The boundary of causal inference.** Closed (decision 0007): causal graphs, instruments, differences in differences, regression discontinuity and synthetic control are peira's (`quasi`, `graphs`), panels tyche's (`panel`).
2. **Exact games.** Equilibria are computed in floating point; exact rational support enumeration and Lemke–Howson on morphe's rationals are possible when exact payoffs matter.
3. **Large markets and profiles.** The algorithms are for the sizes of research examples; large instances (national matching markets, Kemeny rankings of many candidates) would need better data structures or heuristics.

## 9. Current State

As of 2026-10-25, the first six modules are written ([decision 0003](development/decisions/0003-the-first-modules.md)): `measurement`, `indices`, `games`, `matching`, `choice` and `econometrics`, each with its checks. `agents` followed on 2026-11-04 ([decision 0004](development/decisions/0004-agent-based-models.md)). `psychometrics`, `surveys`, `networks` and `mechanisms` are planned.
