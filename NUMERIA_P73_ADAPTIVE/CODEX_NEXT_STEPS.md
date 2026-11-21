# Codex continuation · NUMERIA P73

Branch: `numeria-p73-random-adaptive`, based on P72 branch `numeria-p72-zero-side-codex`.

1. **Reproduce** `certificate_20_t46.6.json` exactly with `certify_adaptive.py` and the committed `arb_zeros_20.json`; match `cell_coverage_sha256` and rational margin.
2. **Reauthenticate** the zeros: run P72 `export_zeros_flint.py` with python-flint 0.9.0 / 192 bits, compare bytes and JSON SHA-256 to the committed P72 zero file. Also independently audit FLINT documentation and zero indexing.
3. **Review mathematics separately**: Suzuki (1.3), Platt–Trudgian RH finite height, Bellotti–Wong bound and Stieltjes integration, curvature interpolation, upper bound for complex zero terms. A passing workflow alone does not discharge this literature review.
4. Compare uniform/random/adaptive sampling using multiple seeds and budgets, not one seed as universal evidence. Keep all stochastic data outside the proof trace.
5. Add 50, 100, 500, then 3,000 Arb zeros if needed; attempt `t=47.48` with **strictly positive** cell-by-cell lower bound, document negative lower-bound cells as inconclusive.
6. Keep PRs #71 and #72 untouched; preserve checksums and do not merge before code and analytic review.

A passing positive rational certificate only proves a conditional mathematical implication until external sources and zero provenance are independently reviewed. Never describe a new finite interval as a solution of RH.
