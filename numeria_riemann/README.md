# NUMERIA-RIEMANN: integrated numerical research module

This package is isolated from the SciCalc application. It adds reproducible
numerical experiments and optional cross-checks without changing calculator
files or claiming to prove the Riemann Hypothesis.

## Integrated components

- `explicit_formula.py`: von Mangoldt function, finite exact inclusive
  Chebyshev ψ baseline, and a truncated explicit-formula diagnostic.
- `backend_compare.py`: zero-residual reports with mpmath and optional
  cross-backend ζ evaluation using python-flint.
- `benchmark.py`: CSV outputs, environment/version metadata, and SHA-256
  manifest for repeatable runs.
- `test_*.py`: unit tests run by GitHub Actions.
- `REFERENCES_AND_PROVENANCE.md`: upstream references and licensing boundaries.

## Install and test

Python 3.10+:

```bash
python -m pip install -r numeria_riemann/requirements.txt
python -m unittest numeria_riemann.test_explicit_formula numeria_riemann.test_backend_compare -v
```

The optional independent backend can be installed with:

```bash
python -m pip install python-flint
python -m unittest numeria_riemann.test_explicit_formula numeria_riemann.test_backend_compare -v
```

## Run a reproducible benchmark

```bash
python -m numeria_riemann.benchmark --out numeria_riemann/results/run1 --zeros 20 --dps 40 --x 10.5
```

The output folder contains `zero_residuals.csv`, `explicit_formula_convergence.csv`,
`metadata.json`, and `SHA256SUMS.json`. Do not commit generated results as
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
