# NUMERIA–RIEMANN: interval audit of the finite Weil matrix

**Reference configuration:** L = 0.45, N = 6, R = 160, requested step = 0.02, 8,000 integration cells, 128 terms in the digamma series, Arb precision = 192 bits.

## Result

The GitHub Actions run for commit `97e832b61b3479b1c261e54d6dc81f2484ea6cbc` reported:

- Status: `certified_positive_definite_finite_matrix`
- Active prime powers: `n = 2`
- Dimension after the two exact pole constraints: 4
- Interval LDL pivot lower bounds:
  1. 0.2149252225123455007633539770168471790632
  2. 0.6591348761638142069767114402974921544929
  3. 0.3643216324439690834660359698579271340626
  4. 0.9835229029954672328707609368191653524151

All four pivot intervals have strictly positive lower bounds. For the interval-enclosed symmetric matrix and the exact constrained subspace represented by the interval basis, this certifies positive definiteness of this **finite-dimensional matrix**, rather than merely reporting a floating-point eigenvalue.

## Error-control path

1. The finite archimedean integral on [0,R] is enclosed by Arb interval range sums. For each cell I, the integral is contained in width(I) times the interval enclosure of the integrand over I.
2. The real digamma factor is evaluated using reflection/duplication identities and a positive series. Its omitted series tail is enclosed by the integral test plus the first omitted term.
3. The remaining archimedean integral over [R, infinity) is enclosed by an analytic bound based on Binet's formula and the decay bound for the Fourier transform of the sine basis.
4. Prime-power overlap integrals and the two pole-constraint vectors are evaluated with Arb interval arithmetic. The constrained subspace is parameterized by solving the first two constraint columns; the resulting basis coefficients remain interval enclosures.
5. Interval LDL^T is attempted on the projected matrix. A certificate is emitted only if every pivot's lower bound is strictly positive. An inconclusive pivot must not be interpreted as a negative eigenvalue.

## Scope and caveats

This result applies only to the specified finite sine-basis matrix, support, prime cutoff, and the stated form of the Weil quadratic form. It does **not** prove positivity for every admissible test function, does not establish the global Weil criterion, and does not prove or refute the Riemann Hypothesis.

The implementation and its analytic inequalities should receive independent mathematical review before this result is treated as a publication-grade certificate. The reproducible workflow reruns the reference case and fails if the positive interval LDL condition is no longer met.

## Reproduction

```bash
python -m pip install python-flint==0.9.0
python -m numeria_riemann.weil_interval_audit --L 0.45 --N 6 --R 160 --step 0.02 --series-terms 128 --precision 192 --out weil-interval-audit.json --require-positive
```

Workflow artifact: [download the interval audit report](https://github.com/jbarumby/scicalc/actions/runs/39180103550/artifacts/12155479271).
