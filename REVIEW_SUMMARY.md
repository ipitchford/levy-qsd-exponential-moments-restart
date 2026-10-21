# Review summary (public, rating-free)

This release was examined only by AI-assisted review. That is not peer review, independent reproduction or formal verification. The original reports are kept private. This summary records the substantive objections and how the revised manuscript deals with each one. Section, theorem and table numbers refer to the released manuscript.

## What was done

1. **First round, on the 20 October 2026 draft (15 pages).**
   - An internal AI-assisted review with five roles: an editor, methodology, domain, perspective and devil's advocate. Each role worked independently of the others.
   - Each role also ran its own checks. The constants were recomputed to 40–60 digits, an independent finite-difference QSD solver was written, and the Brownian counterexample was checked against the source's text.
   - Separately, an AI-assisted adversarial review by a different model (OpenAI gpt-5.6-sol) attempted to refute the core argument.
2. **First revision.** The manuscript was rewritten and four extensions were added. These are the material now in §§8–9, Theorem 2.6(iii)–(iv) and Proposition 2.9.
3. **Verification round, on the revised text.**
   - An internal AI-assisted verification re-review: a traceability check of every item plus a targeted audit of the new proofs, with independent scale-function and determinant computations.
   - A second adversarial audit of the new or rewritten material only, again by OpenAI gpt-5.6-sol.
   - A second, smaller revision followed.
4. **Replay.** The check scripts were replayed in a separate copy on the same machine, and every data output was byte-identical. Only the interval-certificate script changed after that replay; it was rerun, and its output is byte-identical across two runs.

No reviewer found a mathematical error in the core chain. That chain runs from the Esscher-centred QSD to the diffusive sine limit, then to the endpoint Laplace masses, the three moment regimes and the restart rate. The objections below therefore concern faithfulness to the source, scope, labelling, positioning and readability, together with the conditions attached to the new extensions.

## First-round objections and how the revision deals with them

| Objection | Disposition in the released manuscript |
|---|---|
| The source's conjecture and restart theorem were paraphrased, not quoted. The theorem numbers followed the journal version, which cannot be matched against the open preprint | §1.2 quotes verbatim, from arXiv:2503.03565v3, the conjecture, Theorem 3.14 with its display (17), Corollary 3.15 and Lemmas 5.1, 6.1 and 6.2. Table 1 maps the notation and gives both numberings. The journal number of Lemma 6.1 is marked as unverified |
| "Grow as … exactly as for Brownian motion" can be read in several ways, and the reading decides whether the paper proves or partly refutes the conjecture | §1.3, Table 2: true on the logarithmic scale; true up to constants; a ratio limit exists with an explicit `C*`; the ratio-one reading is false. Proposition 2.9 shows `C* = 1` for Brownian motion, `C* < 1` for every spectrally negative process with jumps, and `C* > 1` for every spectrally positive one. Part (iv) encloses `C* = 0.8501…` in the source's own model |
| The Brownian counterexample had to be checked against the source's exact hypotheses, the corrected range was not stated, and the effect on the source's Brownian corollary was not drawn | Corollary 2.7: the source's display (17) holds if and only if `B = O(e^{ηL})`, and the ratio tends to 0 once `B e^{-ηL} → ∞`. Remark 2.8 quotes the faulty step of the source's Proof 6.4 and notes that the example also refutes Corollary 3.15. The source's exact restart law and lower budget range are unaffected. Random walks are not addressed |
| "Uniformly accurate" could be read as absolute or as relative error | Theorem 2.6(iii) proves uniform relative error, at most `\|a_L − 1\|` at every positive budget. Theorem 2.6(iv) shows that the failure probability is not uniformly approximated. Absolute statements are now labelled "absolute" |
| Much of the machinery is classical: the subcritical regime recovers Kyprianou–Palmowski, and it was unclear what the exact constant adds to a cheaper two-sided bound | §1.5 explains what the ratio limit adds: `R*`, the budget rule, the exact corrected range, and the sign of `C* − 1`. §1.6 and Table 3 separate the core results from consequences and from recovered results. The sine bulk, the `L^-2` eigenvalue correction and the renewal shape of the layers are stated not to be new |
| Literature was missing: two-boundary problems (the Lotov school, Kadankova–Veraverbeke), classical fluctuation theory, QSD background, restart and resetting, splitting and Fleming–Viot | §1.6 adds these, with every reference checked against a primary record. The expansion of `κ_L` and the exit split `Q_L → 1/2` are described as closely related to, and plausibly implied by, the Lotov-school expansions, and neither is claimed. Several of those works could be read only through reviews (see Remaining limits) |
| The Gaussian component does most of the technical work, and the source also covers pure-jump processes and random walks | Theorem 8.1 removes the Gaussian component for spectrally negative processes under a scale-function positivity condition (Q). §12 lists the four places where `σ > 0` is used for two-sided jumps |
| Hard to read for non-specialists: no heuristic, symbol clashes, an orientation easily confused | §1.1 states the problem in plain terms. §1.4 gives the heuristic and Figure 1. Table 4 assigns each symbol one meaning, and §1.7 gives a roadmap. Remark 2.4 explains which ladder sits at which end; swapping them would give 1.2926 instead of 0.8501 |
| Some proof steps were compressed: moving starting points, full support, class-D stopping, interval exit versus first passage | Lemma 3.1 now has four explicit steps. Remark 2.1 proves full support. Lemma 4.4 uses bounded optional stopping, and finiteness follows rather than being assumed. Lemma 4.3 states the comparison explicitly |
| The finite-`L` error of about `0.5/L` was visible in the tables but not reported | Table 6 reports it, with Richardson extrapolation. §9 identifies the first-order terms in the source model (Theorems 9.1–9.2), with conditions given below |
| No guidance for practitioners, and the cost of sampling the QSD was not sized | §11 and Table 7 compare plain simulation, QSD restart and splitting-type methods, and footnote (b) notes that the last are not directly comparable. A budget rule is given. Fleming–Viot burn-in and discretisation are discussed and labelled heuristic |
| Workflow statements appeared in the abstract, and the numerical claims could not be checked from the text alone | These are removed. An AI-use statement and a data and code statement have been added. Every number is regenerated by deterministic scripts (Appendix A, `verify_all.py`) |

## Verification-round objections

In the verification round, no defect was found in Theorem 2.3, Theorem 2.6(iii)–(iv), the one-sided sign results of Proposition 2.9, Corollary 2.7 or the source-proof diagnosis in Remark 2.8. The independent computations agreed with the paper: `δ`, `g`, the values of `L²α_L` in Table 6, the spectrally negative constant (including a new compound-Poisson model) and a spectrally positive model with `C* > 1`.

| Objection | Disposition |
|---|---|
| Hypothesis (Q) is not proved for the pure-jump example `ψ(s) = s − 2s/(1+s)`. Grid checks cannot establish it for every large `L`, yet `C* = √2/2` was stated as exact | Resolved by relabelling, not by proof. §2.4, §8 and Table 3 now say that (Q) is verified only numerically for this model and that `C* = 1/√2` holds *if (Q) holds*. The abstract states the extension "under a scale-function positivity condition" |
| The effective-length theorems depend on the regularity fact (F1), but the theorem statements did not say so, and "standard regularity" was neither proved nor cited | Resolved by labelling. Theorems 9.1 and 9.2 now begin "Assume (F1)" and carry "conditional on (F1)" in their titles; the abstract and Remark 9.7 match. (F1) is still neither proved nor tied to a specific cited theorem |
| The 20-digit "certified" enclosure of `C*` could not be checked from the text | Appendix B now gives the bisection polynomial, the rational isolating interval and outward-rounded enclosures of every intermediate quantity. The script prints directed-rounding endpoints. The four-point factorisation check is described as a consistency test, not part of the certificate |
| In Step 4 of the proof of Theorem 8.1, the Brownian case (no jumps) needed separate treatment | Resolved: Step 4 now treats `Π = 0` separately |
| In Step 3 of the proof of Theorem 8.1, the identity `∫₀^∞ e^{sw} Π̄(w) dw = s^{-1}∫(e^{s\|x\|} − 1) Π(dx) < ∞` can fail near zero for unbounded-variation jumps. Only the tail beyond `R` is needed | **Not changed** in the released text. The reviewer judged the argument repairable by keeping the truncation, and the conclusion is unaffected |
| Positioning: Kadankova–Veraverbeke's model is not the source's model, and the scope of the Asmussen–Avram–Pistorius identities was misdescribed | Both sentences are corrected in §1.6. The Asmussen–Avram–Pistorius description follows that paper's title and subject; the paper was not reread in full |
| A sentence in Remark 9.7 referred to an earlier version of a table | Deleted |

## Remaining limits

- **Restart model.** The restart results assume exact, free, independent QSD draws and continuous monitoring. Fleming–Viot or other approximate sampling, restart overhead and discrete monitoring are treated only heuristically.
- **Scope.** There are no results for random walks, including the random-walk case of the source's restart theorem, or for two-sided pure-jump processes.
- **Conditional results.** The pure-jump spectrally negative example depends on (Q), which is checked only numerically. The finite-`L` theorems (§9) depend on (F1), which is not proved. The second-order finite-`L` constants are numerical only. Conjecture 9.9, on the general effective length, is exploratory.
- **The interval certificate** relies on the closed form for `a±`. That form is checked numerically to about `1e-58`, but has not been verified symbolically.
- **Sources read only in part.** Lotov (1979), Lotov (2016) and Bertoin (1997) could not be opened; Lotov–Khodzhibaev (1998/99) was read only through reviews. The journal numbering of the source's Lemma 6.1 was not re-checked, because the publisher page was behind a bot check. The original wording of Kolb–Savov was not re-read in this cycle.
- **Prior-art search.** The search was bounded, so no historical-priority claim is made. The source's authors have not been contacted about the budget-range correction.
