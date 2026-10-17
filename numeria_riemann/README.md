# NUMERIA-RIEMANN: explicit-formula prototype

This is an isolated research prototype inside the SciCalc repository. It does
not modify existing calculator files or claim to prove the Riemann Hypothesis.

## Purpose

- Compute the von Mangoldt function and the exact finite baseline
  psi(x) = sum_{n <= x} Lambda(n).
- Evaluate a finite truncation of the classical explicit formula using
  ordinates of non-trivial zeta zeros.
- Report truncation error against the inclusive exact baseline.

## Run

Requires Python 3.10+ and mpmath:

    python -m pip install mpmath
    python -m unittest numeria_riemann.test_explicit_formula

Example:

    from numeria_riemann.explicit_formula import convergence_report
    for row in convergence_report(10.5, (1, 5, 10, 20), dps=30):
        print(row)

## Scientific limitations

The finite zero sum is a diagnostic, not a certified evaluation. Its behavior
depends on truncation and can be oscillatory. The standard explicit formula
uses a midpoint convention at prime powers; psi_exact instead returns the
inclusive staircase value. This prototype does not yet implement smoothing,
rigorous interval arithmetic, or an independently certified zero list.

## Next validation step

Run the unit tests, compare multiple truncation levels at points away from
prime powers, and record dependency versions, inputs, output CSV, precision,
and a hash manifest before treating results as benchmark evidence.
