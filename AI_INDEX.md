# Agent navigation

## Identity

- Release slug: `levy-qsd-exponential-moments-restart`; version `1.0.0-candidate`; DOI 10.5281/zenodo.23545672 (concept DOI 10.5281/zenodo.23545671)
- Scholarly creator: Anonymous. Status: unrefereed candidate ([STATUS.md](STATUS.md))
- Manuscript: [PDF](levy_qsd_restart.pdf), [LaTeX source](levy_qsd_restart.tex); human summary and claim-status table: [README.md](README.md); file descriptions: [MANIFEST.md](MANIFEST.md)
- Machine-readable claims: [CLAIMS.json](CLAIMS.json); frozen target with hypotheses, notation and fingerprints: [target.yaml](target.yaml)
- Assurance limits: [ASSURANCE.md](ASSURANCE.md); internal review summary: [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md)
- Provenance: [PROVENANCE.md](PROVENANCE.md); prior art and limitation triage: [RESEARCH_GATES.json](RESEARCH_GATES.json)
- Citation: [CITATION.cff](CITATION.cff), [Zenodo metadata](.zenodo.json); licences: [LICENSES.md](LICENSES.md)

## Scientific map

Object: a real Lévy process X with Gaussian coefficient σ > 0, cumulant ψ finite on [0, s_max), Cramér root λ* (ψ(λ*) = 0) and η = argmin ψ on (0, λ*); ν_L is its quasi-stationary distribution on (0, L) and M_L(θ) = ∫ e^{θy} ν_L(dy). Source: Garcia, Bermolen, Jonckheere and Shneer (Stochastic Systems 2026; arXiv:2503.03565v3), who conjectured M_L(λ*) grows like e^{(λ*−η)L}, "exactly the case for the linear Brownian motion".

- **Answer** (Theorem 2.2, proof): M_L(θ) → K_ℓ(η−θ)/K_ℓ(η) for θ < η; M_L(η)/L² → 1/K_ℓ(η); e^{−(θ−η)L} M_L(θ) → C(θ) = K_r(θ−η)/K_ℓ(η) for θ > η. So M_L(λ*) ∼ C* e^{(λ*−η)L}, C* = G₊(λ*−η)/G₋(η).
- **Engine** (Theorem 2.3, proof): endpoint masses of the Esscher-centred QSD, L²A_L(−b) → K_ℓ(b), L²e^{−bL}A_L(b) → K_r(b); left layer = ascending ladder renewal function, right layer = descending.
- **Ratio one fails** (Proposition 2.9, proof; certified value): C* = 1 Brownian; C* < 1 spectrally negative with jumps; C* > 1 spectrally positive with jumps; source model C* = 0.85014864899… ([checks/e4_ratio_one.json](checks/e4_ratio_one.json)).
- **Restart** (Theorem 2.6, Corollary 2.7, proof): r_L ∼ R* e^{−ηL}; sup_B |p_L(B)/(1 − exp(−R* B e^{−ηL})) − 1| ≤ |a_L − 1|; the source restart theorem's linear regime is B = O(e^{ηL}), not B = o(e^{λ*L}).
- **Extensions**: spectrally negative without Gaussian part under (Q) (Theorem 8.1); effective length and 1/L corrections in the source model (Theorems 9.1–9.2, given (F1)); Conjecture 9.9 exploratory.

Do not read: "certified" as covering the general theorems; numerical constants as enclosures; the restart results as statements about Fleming–Viot or discretely monitored implementations; any claim for random walks or two-sided pure-jump processes.

## Replay

Python 3.13 with numpy, scipy, sympy, mpmath, matplotlib ([requirements.txt](requirements.txt); tested versions in [ENVIRONMENT.txt](ENVIRONMENT.txt)):

```sh
python3 -B e4_ratio_one.py --out /tmp/c        # certified enclosure of C* and R*
python3 -B -O negative_controls.py             # swapped orientation, ratio one and tampering are rejected
python3 -B verify_all.py                       # rerun every check; byte-for-byte comparison with checks/
shasum -a 256 -c SHA256SUMS
```

[verify_all.py](verify_all.py) reruns every script in a temporary directory. Scripts: [verify.py](verify.py), [double_exponential_prefactor.py](double_exponential_prefactor.py), [manuscript_constant_check.py](manuscript_constant_check.py), [e4_ratio_one.py](e4_ratio_one.py), [e1_effective_length.py](e1_effective_length.py), [e2_spectrally_negative.py](e2_spectrally_negative.py), [fig_boundary_layers.py](fig_boundary_layers.py); outputs and figures: [checks/](checks/).

## Useful next research question

Quantify implementable QSD restart (Fleming–Viot approximation error, restart overhead, discrete monitoring), or extend the endpoint-mass theorem to random walks and two-sided pure-jump processes.
