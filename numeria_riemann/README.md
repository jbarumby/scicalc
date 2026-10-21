# NUMERIA-RIEMANN: integrated numerical research module

This package is isolated from the SciCalc application. It adds reproducible
numerical experiments and optional cross-checks without changing calculator
files or claiming to prove the Riemann Hypothesis.

## Integrated components

- `explicit_formula.py`: von Mangoldt function, finite exact inclusive
  Chebyshev ψ baseline, and a truncated explicit-formula diagnostic.
- `backend_compare.py`: zero-residual reports with mpmath and optional
  cross-backend ζ evaluation using python-flint.
- `number_theory_backends.py`: optional SymPy checks for von Mangoldt values,
  prime lists, and Möbius values; optional primesieve prime-list comparison.
- `benchmark.py`: CSV outputs, environment/version metadata, and SHA-256
  manifest for repeatable runs.
- `test_*.py`: unit tests run by GitHub Actions.
- `REFERENCES_AND_PROVENANCE.md`: upstream references and licensing boundaries.

## Install and test

Python 3.10+:

```bash
python -m pip install -r numeria_riemann/requirements.txt sympy
python -m unittest numeria_riemann.test_explicit_formula numeria_riemann.test_backend_compare numeria_riemann.test_number_theory_backends numeria_riemann.test_flint_integration -v
```

The optional independent backend can be installed with:

```bash
python -m pip install python-flint sympy
python -m unittest numeria_riemann.test_explicit_formula numeria_riemann.test_backend_compare numeria_riemann.test_number_theory_backends numeria_riemann.test_flint_integration -v
```

## Run a reproducible benchmark

```bash
python -m numeria_riemann.benchmark --out numeria_riemann/results/run1 --zeros 20 --dps 40 --x 10.5
```

The output folder contains `zero_residuals.csv`, `explicit_formula_convergence.csv`,
`metadata.json`, and `SHA256SUMS.json`. Additional optional checks are available via
`compare_prime_lists(limit)` and `compare_mobius_values(limit)` in
`number_theory_backends.py`. Do not commit generated results as
validated evidence without recording the exact runtime and reviewing outputs.

## Scientific status and limitations

- mpmath zeros are numerical approximations; small `|zeta(rho)|` is a
  residual check, not a proof that a zero is certified or isolated.
- python-flint/Arb can provide ball enclosures. Midpoint agreement between
  libraries alone does not constitute a rigorous certificate.
- The truncated explicit formula can oscillate and is not a certified
  evaluation. At prime powers the usual formula uses a midpoint convention,
  while `psi_exact` reports the inclusive staircase value.
- No code from the Riemann-Zeta-Zero-Finding-Suite is copied: its licensing
  notice reserves software artifacts absent separate written permission.
- The package does not prove RH, global Weil positivity, or a Hilbert–Pólya
  operator construction.

See [REFERENCES_AND_PROVENANCE.md](REFERENCES_AND_PROVENANCE.md) for sources,
license checks, and what remains to be independently validated.

## Latest audit adendum

The key findings from the supplied 2026-10-18 NUMERIA–RIEMANN audit package are summarized in [AUDIT_ADENDUM_2026-10-18.md](AUDIT_ADENDUM_2026-10-18.md). Its reported finite checks are not a proof of RH or related conjectures. In particular, the Weil-matrix sweep did not resolve any negative minimum eigenvalue beyond estimated quadrature spread, and the Lean source-integrity checks did not include `lake build` or `#print axioms`.

GitHub Actions also runs an exploratory benchmark and publishes its CSVs, metadata, and SHA-256 manifest as a downloadable workflow artifact. These artifacts record the computation; they are not mathematical certificates.

## Weil matrix: interval audit (experimental)

`weil_interval_audit.py` is a separate Arb-based path for the finite sine-basis matrix. It encloses the finite archimedean integral with Arb certified adaptive integration (`--archimedean adaptive`, default; the older interval range sums remain as `--archimedean rangesum` for cross-checks), bounds the remaining integral tail analytically, computes prime-overlap terms and the two pole constraints with Arb intervals, then tries two positivity tests on a parameterization of the exact constrained subspace: interval LDL and a Rump-type margin test (`M - ||rad||∞·I` positive definite).

Run the reference case:

```bash
python -m pip install python-flint==0.9.0
python -m numeria_riemann.weil_interval_audit --L 0.45 --N 6 --R 160 --precision 192 --out weil-interval-audit.json
```

The report status is either `certified_positive_definite_finite_matrix` (interval LDL or the Rump-type margin test passed) or `inconclusive_interval_ldl`. The latter is not evidence of a negative eigenvalue. Even a successful result certifies only the selected finite-dimensional constrained matrix, not positivity of the full Weil criterion and not RH. The workflow publishes the report as an artifact. Review the assumptions and interval bounds before interpreting any certificate. The reference case's initial successful interval certificate and its pivot lower bounds are recorded in [WEIL_INTERVAL_AUDIT_RESULTS_2026-10-19.md](WEIL_INTERVAL_AUDIT_RESULTS_2026-10-19.md). Results of the v2 enclosures, their independent validation, and open caveats are in [WEIL_INTERVAL_AUDIT_V2_2026-10-21.md](WEIL_INTERVAL_AUDIT_V2_2026-10-21.md). The analytic tail bound is proved in [TAIL_PROOF.md](TAIL_PROOF.md) ([PT](TAIL_PROOF_PT.md)), and the derivation of the quadratic form in [WEIL_FORM_DERIVATION_PT.md](WEIL_FORM_DERIVATION_PT.md). Reports export outward-rounded bounds; `python -m numeria_riemann.verify_certificate REPORT.json` re-checks the final step in exact rational arithmetic. The 2026-10-19 certificates remain superseded.
