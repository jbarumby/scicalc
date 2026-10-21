# NUMERIA–RIEMANN: finite Weil interval audit, v2 enclosures

**Date:** 2026-10-21. **Status:** certificates for the *finite constrained matrix only*, under the stated parameters. Not global Weil positivity, not RH.

Parameters: L = 0.45, R = 160, 192-bit Arb, active prime powers {2}, `--archimedean adaptive`.

## Result

| N | Projected dim. | Status | LDL | Rump-type margin (δ) | Min. eigenvalue of midpoint (non-rigorous) |
|---|---:|---|---|---|---|
| 6  | 4 | certified | pass | pass (δ ≈ 0.0127) | ≈ 0.133 |
| 8  | 6 | certified | pass | pass (δ ≈ 0.0279) | ≈ 0.105 |
| 10 | 8 | certified | pass | pass (δ ≈ 0.0518) | ≈ 0.103 |

N = 10 interval-LDL pivot lower bounds: 0.2298; 0.8884; 0.5022; 1.1256; 1.1088; 1.1704; 1.6935; 1.4274. Raw reports and per-run diagnostics are in `results/`.

## Why earlier runs were inconclusive

The post-correction runs (2026-10-20) were inconclusive because the enclosures were too wide, not because a negative direction was found:

1. **Range sums** (`width(I) * f(I)`) enclose each factor separately (dependency problem). Entry radii were 0.1–1; they are now ~1e-12 with `acb.integral` (mean radius 0.077 → 1e-13 on the N = 6 archimedean block).
2. Digamma series truncation contributed a radius of order `r²/(M³)` at large r.
3. Structurally zero entries (i − j odd) received tail radius needlessly; they are now exactly 0.
4. Interval LDL accumulates wrapping error pivot by pivot. The Rump-type test (`M − ‖rad‖∞ I` positive definite ⇒ every symmetric matrix inside the enclosure is positive definite) avoids that.

Old and new enclosures are consistent: every legacy range-sum entry overlaps the adaptive one, and all legacy pivot lower bounds (e.g. 0.2015; 0.8598; 0.2929; 1.0805 for N = 10, step 0.005) lie below the new ones.

## Rigorous export and exact re-verification

- **Export.** Decimal bounds are rounded *outward* with exact rational arithmetic (lower: floor, upper: ceil). The previous export used `str(40)`, which rounds to nearest and, in a concrete example, put a "lower bound" *above* the true bound. Balls also carry exact `mid_exact`/`rad_exact` as `[mantissa, exponent]`. A randomized test (3000 balls, tiny and huge radii, negative midpoints) found 0 violations after fixing a context-rounding bug in my first implementation. Keys listed in `display_only_keys` (`interval_ldl_pivots`) are display strings and are not rigorous.
- **Exact verifier.** `python -m numeria_riemann.verify_certificate REPORT.json` reads only the exported decimal endpoints as exact rationals (no Arb, no floats) and redoes the Rump-type test with an exact rational LDL of `M − δI`. It validates the final step (exported enclosure ⇒ positivity), not the Arb construction of the enclosures.
- **Reproduction at another precision.** Runs at 256 bits with integral tolerance 2⁻⁷⁰ certify again, pass the exact verifier, and every pair of enclosures intersects the 192-bit run (`--compare`).

| N | 192 bits: δ / min shifted pivot | 256 bits: δ / min shifted pivot | Reports consistent |
|---|---|---|---|
| 6  | 0.012676 / 0.21892 | 0.012676 / 0.21892 | yes |
| 8  | 0.027929 / 0.20366 | 0.027929 / 0.20366 | yes |
| 10 | 0.051784 / 0.17981 | 0.051784 / 0.17981 | yes |

## Independent validation performed

- Three archimedean entries agree with `mpmath` quadrature to ~1e-21 (inside the Arb radii ~1e-13).
- `_overlap` equals a direct `mpmath` quadrature of `(1/L)∫ sin(w_i(x+a)) sin(w_j x) dx` (≤ 1e-14), and the zero-shift overlap contains the identity Gram matrix (regression for the earlier factor-2 error).
- The pole-constraint basis annihilates `∫ g_k(x) e^{±(x−L)/2} dx` numerically (≤ 1e-14).
- **End to end** (`NUMERIA_SLOW=1`): for random vectors in the constrained space, `cᵀAc` plus the archimedean tail beyond R equals `Σ_zeros |F_c(γ)|²` over the first 150–300 real zeta zeros plus a density tail, to ~1e-8. This is numerical evidence that normalization, signs, the prime term and the constraints are mutually consistent. It is not rigorous, and it assumes the zeros are the real ones.

## Open caveats

- The derivation of the quadratic form from the Guinand–Weil formula, with every sign and normalization factor mapped to code and tests, is in [WEIL_FORM_DERIVATION_PT.md](WEIL_FORM_DERIVATION_PT.md) (Portuguese). It relies on a cited theorem whose exact hypotheses still need a human cross-read.
- The analytic tail constants are now proved in [TAIL_PROOF.md](TAIL_PROOF.md) (Lemmas 1–4, Proposition 5); the proof is new and should be reviewed independently. The explicit-formula normalization is validated numerically (zeta-zero check), not proved.
- The tail radius still dominates the remaining enclosure widths (∝ R⁻³); raise `--R` for more margin at larger N.
- Parity-zero entries rely on `Re(F_i · conj F_j) = 0` for odd i − j, which follows from the Fourier-amplitude phase and is inherited from the earlier code.
- A certificate here says nothing about other L, other bases, or the global criterion.
