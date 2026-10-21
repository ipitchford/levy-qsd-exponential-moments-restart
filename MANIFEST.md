# Manifest

| File | Description |
|---|---|
| `levy_qsd_restart.pdf` | Manuscript (35 pages), built with pdflatex (two passes) |
| `levy_qsd_restart.tex` | LaTeX source; figures are read from `checks/` |
| `README.md` | Plain-English summary, status of each claim, replay instructions |
| `target.yaml` | Machine-readable statement of hypotheses, claims and labels |
| `verify_all.py` | Reruns every script into a temporary folder; compares outputs byte for byte with `checks/` |
| `verify.py` | Brownian identities and counterexample; rational lattice conjugacy; lattice diagnostic; Figures 2-3 |
| `double_exponential_prefactor.py` | 70-digit decimal factorisation and two formulas for `C_*` (standard library only) |
| `manuscript_constant_check.py` | 85-digit constants; exact four-mode finite-interval solution of the source model (Table 6) |
| `e1_effective_length.py` | 70-digit finite-interval solver, effective length and 1/L coefficients (Section 9, Table 5) |
| `e2_spectrally_negative.py` | Scale-function computations for spectrally negative models (Section 8) |
| `e4_ratio_one.py` | Interval-arithmetic enclosure of `C_*`, `R_*` in the source model (Prop. 2.9(iv)) |
| `fig_boundary_layers.py` | Figure 1 |
| `checks/*.json, *.csv, *.txt` | Deterministic outputs of the scripts above |
| `checks/*.pdf, *.png` | Figures |
| `negative_controls.py` | Semantic negative controls (standard library): swapped ladder orientation, ratio one, tampering, linear budget regime |
| `STATUS.md`, `ASSURANCE.md`, `PROVENANCE.md`, `AI_INDEX.md`, `CLAIMS.json`, `RESEARCH_GATES.json`, `REVIEW_SUMMARY.md` | Release records: status, assurance vector, provenance, agent index, machine-readable claims, dated prior-art and limitation record, rating-free review summary |
| `LICENSE`, `LICENSE-CODE`, `LICENSES.md`, `PUBLIC_DOMAIN.md`, `CITATION.cff`, `.zenodo.json` | Licences (CC0-1.0 prose and data, MIT code) and citation metadata |
| `ENVIRONMENT.txt`, `requirements.txt`, `.github/workflows/replay.yml` | Tested environment, pinned packages and the quick CI replay |
| `SHA256SUMS` | SHA-256 of every other file in the package |
