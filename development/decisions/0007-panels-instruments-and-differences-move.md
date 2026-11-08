# 0007. Panels, Instruments and Differences in Differences Move Out

- **Status:** Accepted, provisionally (2026-11-08). Amends decision 0002 (the econometrics module's scope).
- **Decided by:** the infrastructure's ADR 0082 (stage 2), on the user's go-ahead.

## 1. Decision

From `econometrics`:
- **To tyche:** panels and the absorption of fixed effects (`tyche::panel`; tyche's decision 0014), and the robust covariances' pieces (`tyche::regression`'s `sandwich.hpp`).
- **To peira:** instrumental variables and differences in differences (`peira::quasi`; peira's decision 0004).

Each went with its checks.

nomos keeps the econometrics of choices and counts that the social sciences use: average marginal effects of logit and probit, the multinomial logit, and Poisson pseudo-maximum likelihood with two sets of absorbed effects (the gravity equation's). `common.hpp` names what those read from tyche.

## 2. Rationale

Panels are inference within a model, which is tyche's. Instruments and differences in differences identify effects by assumptions on the outcome, which is peira's: they serve clinical, ecological and urban studies as much as economic ones.

## 3. Consequences

- nomos's checks fall from 277 to 242, the 35 moved checks passing where they now stand.
- Callers moved with them: the API's toolbox and the engines that compile it.
