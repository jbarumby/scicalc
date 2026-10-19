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
