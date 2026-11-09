# 0009. Choice Among Alternatives

- **Status:** Accepted, provisionally (2026-11-09).
- **Decided by:** the infrastructure's ADR 0088 (the urban priorities' stage 4, demand), at the user's word ("continue the next feature").

## 1. Decision

`econometrics` gains the models of choice among alternatives described by their attributes, on one long form of choice situations (frequency weights, decision makers):
- `conditional_logit.hpp`:
  - McFadden's conditional logit by Newton's method, with model, robust or clustered covariance;
  - willingness to pay as −β_k/β_p, with its delta-method standard error;
- `nested_logit.hpp`: the nested logit with a dissimilarity for each nest (McFadden 1978);
- `mixed_logit.hpp`: the mixed logit with normal coefficients, simulated with Halton draws over each decision maker's panel, its gradient analytic (Revelt and Train 1998; Train 2009).

## 2. Rationale

Preferences revealed by choices (a route, a destination, a mode, a home) and stated in surveys are measured by these models across the social sciences. The multinomial logit of characteristics of the chooser was already here; these are its counterpart for the alternatives' attributes.

## 3. Verification

`check_nomos_choice_models` (11):
- **Against statsmodels and numpy:**
  - the conditional logit against statsmodels' ConditionalLogit, β̂ to 10⁻¹⁴ and standard errors to 10⁻⁶;
  - willingness to pay and its standard error against numpy.
- **Exact properties:**
  - weights as replication;
  - clustered covariance equal to the robust one times G/(G − 1);
  - the nested logit at λ = 1 equal to the conditional logit;
  - the mixed logit without random columns equal to it.
- **Recovery** of λ = 0.5 and of a normal price coefficient (−1, 0.5²) within four standard errors.

nomos's checks rise from 250 to 261.
