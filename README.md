# Exact exponential moments of interval quasi-stationary laws for Lévy processes, and the saturation of idealised QSD restart

Author: Anonymous. Manuscript: `levy_qsd_restart.pdf` (source `levy_qsd_restart.tex`).

## In plain English

Suppose a simulated random process drifts away from a high target level `L`, and we want to see it reach `L` within a fixed budget of simulated time. A simple trick is *restart*: whenever a run falls below a floor, throw it away and start again from a fresh point drawn from the process's *quasi-stationary distribution* (QSD), the law of the process conditioned to have stayed between the floor and the target for a long time. Garcia, Bermolen, Jonckheere and Shneer (Stochastic Systems, 2026) showed that the time to success is then exactly exponential, and that its rate is governed by an exponential moment of the QSD. They conjectured how fast that moment grows.

This paper computes that moment exactly, for every exponent, under the source's assumptions (a Lévy process with a Gaussian component and a Cramér root), and draws three practical consequences:

* the success rate of idealised QSD restart is `r_L ~ R_* exp(-eta L)`, with an explicit constant `R_*`, so the budget needed for success probability `p` is about `-log(1-p) exp(eta L) / R_*`;
* the curve `1 - exp(-R_* B exp(-eta L))` is accurate at every budget `B` in relative error;
* the linear budget regime claimed in the source's restart theorem (Theorem 3.14 of arXiv v3, Theorem 4 of the journal version) is `B = O(exp(eta L))`, not `B = o(exp(lambda_* L))`; an exact Brownian example shows the original range is too wide.

The conjecture is true on the logarithmic scale and up to constants, but its "ratio one, exactly as for Brownian motion" reading is false: the constant is `C_* = 0.85014864899...` in the source's own example, below one for every spectrally negative process with jumps and above one for every spectrally positive one.

## What is proved, certified, numerical or exploratory

| Statement | Where | Status |
|---|---|---|
| Exact endpoint masses of the centred QSD (core result) | Thm 2.3, §4 | proof |
| All exponential moments of the interval QSD, constants `C(theta)`, `C_*` | Thm 2.2, §5 | proof |
| Restart constant `R_*`; uniform relative-error saturation | Thm 2.6, §6 | proof |
| Corrected budget range of the source's Theorem 3.14 | Cor 2.7, Remark 2.8 | proof |
| Ratio-one reading false; sign of `C_* - 1` for one-sided jumps | Prop 2.9, §7 | proof |
| `C_*` of the source model in a rigorous enclosure, `C_* < 1` | Prop 2.9(iv), `e4_ratio_one.py` | certified (interval arithmetic) |
| Spectrally negative processes without Gaussian part (`C_* = eta^2 psi''(eta)/(2 kappa_inf) <= 1`) | Thm 8.1, Cor 8.2 | proof, under a scale-function positivity condition (Q) that holds when the laws of `X_t` are absolutely continuous; for compound-Poisson-plus-drift models (Q) is only checked numerically |
| Effective length and first-order `1/L` correction in the source model: `L^2 alpha_L = (pi^2 v/2)(L/(L+delta))^2 + O(L^-2)`, `C_L/C_*, R_L/R_* = 1 + g/L + O(L^-2)`, `g = (pi^2/2)(delta_- - delta_+) = 0.4974` | Thms 9.1, 9.2 | proof conditional on an unproved regularity property (F1) of the QSD density; second-order constants numerical |
| Same formulas for the whole class (supported by spectrally negative scale-function models outside the proof class) | Conj. 9.9 | exploratory |
| Finite-`L` tables, lattice diagnostic, figures | §10, App. A | numerical |

The proofs do not depend on any computation. The computations check transcription, constants and finite-`L` behaviour. The subcritical moment limit, the sine law and the `L^-2` killing-rate expansion are recovered or standard and are labelled as such in the paper.

## Replay

Requirements: Python 3.13 with NumPy, SciPy, SymPy, mpmath and Matplotlib ; a TeX distribution with `pdflatex`.

```bash
python3 verify_all.py            # reruns every check into a temporary folder and compares data outputs byte for byte
python3 verify_all.py --figures  # also regenerates and compares the figures
pdflatex levy_qsd_restart.tex && pdflatex levy_qsd_restart.tex
```

Individual scripts accept `--out DIR` (default `checks/`):

| Script | Output | What it checks |
|---|---|---|
| `verify.py` | `verification.json`, two CSVs, Figs 2 and 3 | Brownian identities and counterexample, exact rational lattice conjugacy, lattice diagnostic |
| `double_exponential_prefactor.py` | `prefactor_decimal.txt` | 70-digit factorisation; two formulas for `C_*` |
| `manuscript_constant_check.py` | `manuscript_constant_check.json` | 85-digit constants; exact finite-interval solution of the source model |
| `e4_ratio_one.py` | `e4_ratio_one.json` | certified enclosure of `C_*`, `R_*` |
| `e1_effective_length.py` | `e1_effective_length.json` | high-precision finite-`L` solution and `1/L` coefficients |
| `e2_spectrally_negative.py` | `e2_spectrally_negative.json` | scale-function computations for spectrally negative models |
| `fig_boundary_layers.py` | `boundary_layers.pdf` | Figure 1 |

All outputs are deterministic: no random numbers, no network, no time stamps or version strings, fixed rounding, and round-off-level residuals recorded as thresholds. `SHA256SUMS` lists the hashes of every file in the package; `MANIFEST.md` describes each file.

## Scope

The restart results assume exact, free, independent draws from the QSD and continuous monitoring of exits. Fleming–Viot implementation, discrete monitoring and the random-walk case of the source's theorem are discussed heuristically or left open. Theorem and lemma numbers of the source refer to arXiv:2503.03565v3; the journal version (DOI 10.1287/stsy.2025.0108) numbers them consecutively.

## Use of AI tools

Large language model-based AI tools were used extensively to develop the proofs, computations and text. Data and code DOI: [10.5281/zenodo.23545672](https://doi.org/10.5281/zenodo.23545672) (version 1.0.0-candidate, Evidence Press release page: https://evidencepress.org/releases/levy-qsd-exponential-moments-restart/). Status: unrefereed research candidate; see STATUS.md, ASSURANCE.md and AI_INDEX.md.
