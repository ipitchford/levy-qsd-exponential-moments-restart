# Assurance

| Dimension | State | Note |
|---|---|---|
| Written proofs | present | Complete proofs of every theorem in the manuscript; not machine-checked. The proofs do not depend on any computation |
| Certified numerics | passed | `e4_ratio_one.py`: outward-rounded interval enclosure of C* and R* in the source model (Proposition 2.9(iv)) |
| Exact checks | passed | `verify.py`: exact Brownian identities and the Brownian counterexample to the source's budget range; exact rational lattice conjugacy |
| High-precision numerics | passed | `double_exponential_prefactor.py` (standard library, 70 digits, two formulas for C*), `manuscript_constant_check.py` (85 digits; exact finite-interval solution), `e1_effective_length.py`, `e2_spectrally_negative.py` |
| Producer replay | passed | `verify_all.py` reruns every script into a temporary folder and compares each data output byte for byte. A separate clean replay in an independent copy reproduced all eight data outputs; the released scripts and outputs were compared with that copy by SHA-256 and only changed checks were rerun after DOI insertion |
| Negative controls | passed | `negative_controls.py`: the swapped ladder orientation (C* ≈ 1.2926) and the ratio-one value C* = 1 lie outside the certified enclosure; a one-digit change in a saved output is detected; the linear budget regime is rejected beyond e^{ηL} |
| Independent rerun | not assessed | No unaffiliated rerun |
| Independent reimplementation | partial | Two independent formulas for every constant and an exact finite-interval solution inside the producer workflow; implementation diversity, not external reproduction |
| Formal verification | not assessed | No proof assistant; interval arithmetic is rigorous numerics, not formal verification |
| Specialist review | not assessed | AI-assisted internal and adversarial review is not external specialist review |
| Editorial peer review | not assessed | None |
| Prior art | bounded | Dated bounded search summarised in [RESEARCH_GATES.json](RESEARCH_GATES.json); no priority claim |

Trust boundary: the general statements rest on the written proofs. The certified enclosure covers the stated constants of one model; the Python interpreter, mpmath interval arithmetic and the written reduction from the paper's definitions to the computed quantities are in its trust boundary. Values labelled numerical are deterministic high-precision computations, not enclosures; Conjecture 9.9 is exploratory. The restart results assume exact, free QSD draws and continuous monitoring.
