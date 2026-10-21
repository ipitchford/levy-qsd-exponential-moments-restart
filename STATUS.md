# Status

Version: 1.0.0-candidate (21 October 2026). DOI: [10.5281/zenodo.23545672](https://doi.org/10.5281/zenodo.23545672). Scholarly creator: Anonymous.

**Unrefereed research candidate.** The manuscript has not been peer reviewed by a journal or by named human specialists. It was checked by internal AI-assisted review in several roles, an AI-assisted adversarial review by a different model (OpenAI gpt-5.6-sol) and a verification re-review (see [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md)). That is not peer review, not independent reproduction and not formal verification.

Object: real Lévy processes with a Gaussian component and a Cramér root, killed on leaving (0, L); the conjecture of Garcia, Bermolen, Jonckheere and Shneer, *Efficiency of parallel and restart exploration strategies in model-free stochastic simulations*, Stochastic Systems 16(3), 264–288 (2026), doi:10.1287/stsy.2025.0108 (arXiv:2503.03565v3, paragraph after Corollary 3.15), and their restart theorem (arXiv Theorem 3.14, journal Theorem 4).

* **Proved in writing:** Theorem 2.2 (exact limits of all exponential moments of the interval QSD); Theorem 2.3 (endpoint masses of the centred QSD); Theorem 2.6 (restart rate, uniform relative-error saturation, what is not uniform); Corollary 2.7 and Remark 2.8 (corrected budget range of the source's restart theorem); Proposition 2.9 (ratio-one reading false; sign of C* − 1 for one-sided jumps); Theorem 8.1 and Corollary 8.2 (spectrally negative processes without a Gaussian part, under scale-function positivity (Q), automatic for absolutely continuous laws); Theorems 9.1 and 9.2 (effective length and 1/L corrections in the source model, taking the standard QSD-density regularity (F1) as given).
* **Certified:** the enclosure C* ∈ [0.85014864899051539225, 0.85014864899051539226] and R* for the source model (outward-rounded interval arithmetic, `e4_ratio_one.py`). Exact rational lattice conjugacy and exact Brownian identities (`verify.py`).
* **Numerical:** 70–85-digit constants by two formulas; finite-L tables from the exact four-mode solution; the lattice diagnostic; scale-function computations in Section 8; second-order finite-L constants; condition (Q) for compound-Poisson-plus-drift models.
* **Exploratory:** Conjecture 9.9 (first-order finite-L formulas for the whole class).
* **Open:** Fleming–Viot or other approximate QSD sampling, restart overhead and discrete monitoring; random walks; two-sided pure-jump processes; Conjecture 9.9.

No proof assistant was used. The novelty statement is bounded by the searches summarised in [RESEARCH_GATES.json](RESEARCH_GATES.json); it does not establish priority.
