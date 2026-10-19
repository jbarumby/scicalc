# NUMERIA–RIEMANN — audit notes for adendum v2 (2026-10-18)

## Scope and provenance

These notes summarize the supplied `NUMERIA_RIEMANN_AUDIT_PACKAGE_ADENDO_v2_2026-10-18.zip`. The package report labels itself a partial reproducible audit. The package is not vendored here; these notes preserve its principal findings and limitations without importing its separate scripts or dependencies into the isolated prototype.

The figures below are **reported by the supplied package**. They have not all been independently rerun in this repository.

## Findings worth retaining

| Area | Reported result | Correct interpretation |
|---|---|---|
| File integrity | 48/48 SHA-256 checks for the v1 files; 14/14 upstream Lean repository checks | Integrity/provenance checks only; not compilation or theorem verification |
| Goel–Murty Theorem 4.1 | 1,720/1,720 finite cases matched | Finite consistency checks, not a proof for all parameters or of Hardy–Littlewood |
| Mertens function | Dirichlet recurrence reports M(10^6)=212, M(10^7)=1037, M(10^8)=1928 | Useful independent algorithmic cross-check; values do not imply an asymptotic theorem |
| Multiple zeta values | 1,521 (r,k) pairs; no grid/conjecture mismatches reported; 11,700 sign checks agreed at 100 vs 160 dps | Numerical evidence for the tested grid. Sign-change counts can miss even-multiplicity roots or pairs within one cell; this does not prove the conjecture or imply RH |
| Weil matrix sweep | 20 configurations; 11 negative computed minimum eigenvalues; 0 negative cases resolved beyond estimated quadrature spread | The computed signs are numerically unresolved; no positivity or non-positivity conclusion is justified |
| Spectral benchmark | 120 ordinates split 80/40; Riemann–von Mangoldt inverse RMSE on test set about 0.959, lower than the tested oscillator models | Descriptive comparison only. The oscillator fit is not a derived spectral correspondence; only 39 test gaps makes spacing statistics low-powered |
| SOS / operator bound | Exact symbolic SOS identities and a positive robust coefficient are reported; min–max lower bound follows for the stated operator/domain assumptions | The bound concerns the modeled Schrödinger operator. The package does not establish a rigorous map from its spectrum to zeta zeros |
| Lean formalization | Fresh clone and repository manifest checks reported as 14 PASS | `lake build` and `#print axioms` were not run. File checks do not certify the Lean declarations |

## Priority improvements

1. **Separate numerical residuals from certification.** Record precision, algorithm, truncation, and error estimates separately. A small residual or a negative eigenvalue smaller than quadrature spread must be marked unresolved.
2. **Make benchmark runs automatic and retain their artifacts.** The GitHub Actions workflow now has a benchmark job that publishes CSVs, environment metadata, and a SHA-256 manifest. Artifacts are diagnostics, not proof certificates.
3. **Test mathematical conventions explicitly.** The explicit-formula implementation now distinguishes positive ordinates (which represent a conjugate pair) from a list of individual complex zeros. Supplying both conjugates must not double-count the pair.
4. **Avoid overinterpreting fitted spectra.** Keep train/test splits and fixed parameters visible; do not describe affine calibration or finite-difference eigenvalues as a theorem connecting an operator to zeta zeros.
5. **Formal verification remains a separate gate.** Before describing Lean results as compiled verification, run the pinned toolchain's `lake build` and inspect `#print axioms` outputs.

## Status

This is a research audit summary, not a proof of the Riemann Hypothesis, global Weil positivity, Hardy–Littlewood, or a Hilbert–Pólya construction. The broader package's experiments are kept separate from this repository's smaller explicit-formula prototype until their inputs, environments, and scripts can be independently reproduced.
