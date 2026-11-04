# measurement: Inequality, Poverty, Concentration and Specialisation

*Status: draft (2026-10-25), [decision 0003](../development/decisions/0003-the-first-modules.md).* Scales, reliability and item response theory are planned in `psychometrics`.

## 1. Concepts

**Distributions.** A distribution is a list of amounts xᵢ held by units with weights wᵢ ≥ 0, frequency or sampling weights; none given means one each. Every index is a functional of the weighted empirical distribution, so it is unchanged when the population is replicated or a unit of weight two is split in two (Dalton's principle of population), and integer weights are replication. Population divisors are used throughout.

**Inequality** is measured relative to the mean, so every index is unchanged when all amounts are scaled, and falls under a transfer from a richer to a poorer unit (the principle of transfers; Pigou 1912, Dalton 1920).

- The *Gini coefficient* is half the relative mean difference, Σᵢⱼ wᵢwⱼ|xᵢ − xⱼ| / (2W²μ), computed in O(n log n) as 2 cov(x, F̂)/μ with the midpoint cumulative weights F̂ (Lerman and Yitzhaki 1984, 1989); the two agree exactly for any weights and ties. One holder of everything among n gives (n − 1)/n.
- The *Lorenz curve* L(p) is the share held by the poorest fraction p (Gastwirth 1971); for a discrete distribution it is a polygon, so interpolation is exact, and G = 1 − 2∫L. The *generalised Lorenz curve* is μL(p) (Shorrocks 1983). *Quantile shares* and the *Palma ratio*, the share of the richest tenth over that of the poorest four tenths (Palma 2011; Cobham and Sumner 2013), are read from it.
- The *generalised entropy family* GE(α) = (E[r^α] − 1)/(α(α − 1)), r = x/μ, with Theil's L = GE(0) (the mean logarithmic deviation) and Theil's T = GE(1) as its limits (Theil 1967; Shorrocks 1980); GE(2) is half the squared coefficient of variation. It is evaluated with expm1 so that it is continuous and accurate at α = 0 and 1.
- *Atkinson's index* A(ε) = 1 − (E[r^{1−ε}])^{1/(1−ε)}, one minus the equally distributed equivalent over the mean, ε the aversion to inequality (Atkinson 1970); an ordinal transform of GE(1 − ε).

**Poverty** below a line z: the Foster–Greer–Thorbecke family FGT(α) = (1/W) Σ_{x<z} w ((z − x)/z)^α, the headcount ratio, the poverty gap and the severity at α = 0, 1, 2 (Foster, Greer and Thorbecke 1984), and the Watts index (Watts 1968).

**Concentration** of a total among holders, by their shares: the Herfindahl–Hirschman index Σs² (Hirschman 1945; Herfindahl 1950) and its normalised form, the concentration ratio CR_k, the entropy index −Σ s ln s (Theil 1967), and Hannah and Kay's numbers equivalents (Σ s^α)^{1/(1−α)}, which unite the others (Hannah and Kay 1977).

**Specialisation.** In a table of amounts by unit and category, the *location quotient* (x_rc/x_r·)/(x_·c/x_··), which is *Balassa's index* of revealed comparative advantage (Isard 1960; Balassa 1965), and its symmetric form (Laursen 2015); the *index of dissimilarity* ½Σ|a/Σa − b/Σb| (Duncan and Duncan 1955); for flows in two directions by category, the *Grubel–Lloyd index* 1 − |a − b|/(a + b), its aggregate, and the aggregate adjusted for the overall imbalance (Grubel and Lloyd 1975).

**Decompositions.** Every generalised entropy index splits over a partition into groups as GE = Σ_g v_g^{1−α} s_g^α GE_g + GE_between, v_g and s_g the population and amount shares (Bourguignon 1979; Shorrocks 1980, which shows these are the only additively decomposable relative indices); FGT(α) = Σ_g v_g FGT_g (Foster, Greer and Thorbecke 1984). The Gini coefficient does not decompose so.

## 2. Content

| Header | Content |
|---|---|
| `inequality.hpp` | `mean(x, w)`, `gini(x, w)`, `lorenz(x, w)` → `LorenzCurve{p, L, mean; at(q), generalised_at(q)}`, `gini(LorenzCurve)`, `quantile_share(c, p0, p1)`, `top_share(c, q)`, `palma_ratio(c)` and `palma_ratio(x, w)`, `generalised_entropy(x, α, w)`, `theil_t`, `theil_l`, `mean_log_deviation`, `atkinson(x, ε, w)`, `coefficient_of_variation(x, w)` |
| `poverty.hpp` | `fgt(x, z, α, w)`, `headcount_ratio`, `poverty_gap_index`, `poverty_severity`, `watts_index` |
| `concentration.hpp` | `shares(x)`, `herfindahl_hirschman(x)`, `herfindahl_hirschman_normalised(x)`, `concentration_ratio(x, k)`, `entropy_index(x)`, `hannah_kay(x, α)` |
| `specialisation.hpp` | `location_quotients(Rectangular)`, `balassa_index(Rectangular)`, `symmetric_location_quotients(Rectangular)`, `dissimilarity(a, b)`, `grubel_lloyd(a, b)` for one category and for spans, `grubel_lloyd_adjusted(a, b)` |
| `decomposition.hpp` | `decompose_generalised_entropy(x, group, α, w)` → `EntropyDecomposition{total, within, between, population_share, amount_share, mean, index, weight}`; `decompose_fgt(x, group, z, α, w)` → `PovertyDecomposition{total, population_share, index, contribution}` |

Weights are `std::span<const double>` (empty for ones); tables are `calx::linalg::Rectangular`. Invalid input (negative weights, zero mean, zeros where a logarithm needs positive amounts) throws `std::invalid_argument`.

## 3. Checks

- `tests/check_inequality.cpp` (35 checks): equal amounts give zero everywhere; one holder of everything among ten gives a Gini coefficient of 0.9, Theil's T of ln 10 and GE(2) of 4.5; the amounts 1, …, 5 have G = 4/15; 1, …, 10 have a Palma ratio of 1 and a middle half holding 27.5/55; the Gini coefficient by covariance against the double sum and against 1 − 2∫L with random weights and ties; integer weights as replication, three copies of the population and scaling by 7.3 change no index; a unit of zero weight counts for nothing; Atkinson against GE(1 − ε) and against the power mean, A(1) against the geometric mean, CV² = 2GE(2), GE continuous at α = 0 and 1 from both sides; a transfer from the richest to the poorest lowers every index; the exponential distribution at 10⁶ midpoint quantiles against its closed forms (G = ½, T = 1 − γ, L = γ, GE(2) = ½, A(½) = 1 − π/4, A(1) = 1 − e^{−γ}, its Lorenz curve) to 10⁻⁶.
- `tests/check_poverty_concentration.cpp` (35 checks): FGT(0, 1, 2, 3) and Watts on 1, …, 10 by hand, the line itself not poor, weights as replication, monotonicity, a transfer among the poor leaving the gap unchanged and lowering severity and Watts; H = 0.365 for shares 50, 30, 15, 5 per cent, CR₂ = 0.8, equal holders (H = 1/n, E = ln n, Hannah–Kay n), a single holder, Hannah–Kay decreasing in α, equal to 1/H at 2 and continuous at 1; location quotients by hand, their two weighted means of one on a random table, NaN for an empty column; dissimilarity at its ends; Grubel–Lloyd by category, aggregate and adjusted, by hand.
- `tests/check_decomposition.cpp` (10 checks): within + between = total for six values of α with an empty group; the weights of Theil's L and T; no between part when groups share a distribution, no within part when groups are equal inside; FGT contributions adding up.

## 4. Dependencies

calx `linalg` (`Rectangular`). In the checks, calx `random` and `verification`.

## 5. References

- Atkinson, A. B. (1970), On the measurement of inequality, *J. Econ. Theory* 2, 244–263.
- Balassa, B. (1965), Trade liberalisation and "revealed" comparative advantage, *Manchester School* 33, 99–123.
- Bourguignon, F. (1979), Decomposable income inequality measures, *Econometrica* 47, 901–920.
- Cobham, A. and Sumner, A. (2013), Is it all about the tails? The Palma measure of income inequality, Center for Global Development, Working Paper 343.
- Cowell, F. A. (2011), *Measuring Inequality*, 3rd ed., Oxford University Press.
- Dalton, H. (1920), The measurement of the inequality of incomes, *Econ. J.* 30, 348–361.
- Deaton, A. (1997), *The Analysis of Household Surveys*, Johns Hopkins University Press.
- Duncan, O. D. and Duncan, B. (1955), A methodological analysis of segregation indexes, *Am. Sociol. Rev.* 20, 210–217.
- Foster, J., Greer, J. and Thorbecke, E. (1984), A class of decomposable poverty measures, *Econometrica* 52, 761–766.
- Gastwirth, J. L. (1971), A general definition of the Lorenz curve, *Econometrica* 39, 1037–1039.
- Gini, C. (1912), *Variabilità e mutabilità*, Bologna.
- Grubel, H. G. and Lloyd, P. J. (1975), *Intra-Industry Trade*, Macmillan.
- Hannah, L. and Kay, J. A. (1977), *Concentration in Modern Industry*, Macmillan.
- Haughton, J. and Khandker, S. R. (2009), *Handbook on Poverty and Inequality*, World Bank.
- Herfindahl, O. C. (1950), *Concentration in the Steel Industry*, PhD thesis, Columbia University.
- Hirschman, A. O. (1945), *National Power and the Structure of Foreign Trade*, University of California Press; (1964), The paternity of an index, *Am. Econ. Rev.* 54, 761–762.
- Isard, W. (1960), *Methods of Regional Analysis*, MIT Press.
- Laursen, K. (2015), Revealed comparative advantage and the alternatives as measures of international specialization, *Eurasian Bus. Rev.* 5, 99–115.
- Lerman, R. I. and Yitzhaki, S. (1984), A note on the calculation and interpretation of the Gini index, *Econ. Lett.* 15, 363–368; (1989), Improving the accuracy of estimates of Gini coefficients, *J. Econometrics* 42, 43–47.
- Palma, J. G. (2011), Homogeneous middles vs. heterogeneous tails, *Dev. Change* 42, 87–153.
- Pigou, A. C. (1912), *Wealth and Welfare*, Macmillan.
- Shorrocks, A. F. (1980), The class of additively decomposable inequality measures, *Econometrica* 48, 613–625; (1983), Ranking income distributions, *Economica* 50, 3–17.
- Theil, H. (1967), *Economics and Information Theory*, North-Holland.
- Watts, H. W. (1968), An economic definition of poverty, in D. P. Moynihan (ed.), *On Understanding Poverty*, Basic Books.

## 6. Open Questions

1. **Sampling inference** for the indices (standard errors of the Gini coefficient by linearisation or the bootstrap) belongs with `surveys`.
2. **Poverty indices of the Sen family** and the decomposition of the Gini coefficient by sources (Lerman and Yitzhaki 1985), when a domain needs them.
