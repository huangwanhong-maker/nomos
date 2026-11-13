"""Reference values for nomos's agreement coefficients and the AHP, from independent implementations.

statsmodels 0.14 (statsmodels.stats.inter_rater: cohens_kappa, fleiss_kappa), krippendorff 0.9 (Castro 2017, the package
`krippendorff`, after Krippendorff 2011), and numpy's eigensolver for the AHP's Perron vector. The values printed are copied, with
this script named, into measurement/tests/check_agreement.cpp and appraisal/tests/check_appraisal.cpp. Generated 2026-11-13.

Run:  PYTHONPATH=<a folder holding statsmodels, patsy and krippendorff> .venv/Scripts/python foundations/nomos/appraisal/tests/fixtures/make_appraisal_references.py
"""
import krippendorff
import numpy as np
from statsmodels.stats.inter_rater import cohens_kappa, fleiss_kappa

# Two coders, twenty units, three ordered categories.
a = [0, 1, 2, 1, 0, 2, 2, 1, 0, 1, 2, 0, 1, 1, 2, 0, 0, 2, 1, 1]
b = [0, 1, 2, 2, 0, 1, 2, 1, 1, 1, 2, 0, 0, 1, 2, 0, 1, 2, 1, 2]
table = np.zeros((3, 3))
for x, y in zip(a, b):
    table[x, y] += 1
print("cohen nominal", repr(cohens_kappa(table).kappa))
print("cohen linear", repr(cohens_kappa(table, wt="linear").kappa))
print("cohen quadratic", repr(cohens_kappa(table, wt="quadratic").kappa))

# Fleiss: ten units, five coders, four categories.
counts = np.array([[5, 0, 0, 0], [2, 2, 1, 0], [0, 4, 1, 0], [1, 1, 1, 2], [0, 0, 5, 0], [3, 2, 0, 0], [0, 0, 0, 5], [1, 3, 1, 0], [0, 1, 4, 0], [2, 0, 0, 3]])
print("fleiss", repr(fleiss_kappa(counts)))

# Krippendorff's example: four coders, twelve units, values 1–5, missing values (nan).
n = np.nan
data = np.array([[1, 2, 3, 3, 2, 1, 4, 1, 2, n, n, n],
                 [1, 2, 3, 3, 2, 2, 4, 1, 2, 5, n, 3],
                 [n, 3, 3, 3, 2, 3, 4, 2, 2, 5, 1, n],
                 [1, 2, 3, 3, 2, 4, 4, 1, 2, 5, 1, n]])
for level in ("nominal", "ordinal", "interval", "ratio"):
    print("krippendorff", level, repr(krippendorff.alpha(reliability_data=data, level_of_measurement=level)))

# The AHP: Saaty's 3 × 3 and an inconsistent 4 × 4.
for P in (np.array([[1, 3, 5], [1 / 3, 1, 3], [1 / 5, 1 / 3, 1]]), np.array([[1, 2, 4, 3], [1 / 2, 1, 3, 1 / 2], [1 / 4, 1 / 3, 1, 1 / 5], [1 / 3, 2, 5, 1]])):
    vals, vecs = np.linalg.eig(P)
    k = int(np.argmax(vals.real))
    w = np.abs(vecs[:, k].real)
    w /= w.sum()
    print("ahp", len(P), repr(vals[k].real), [repr(x) for x in w])
