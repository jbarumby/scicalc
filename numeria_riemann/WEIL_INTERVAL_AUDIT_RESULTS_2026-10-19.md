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

## Controlled dimension increase: N = 8

The workflow run [#39337029515](https://github.com/jbarumby/scicalc/actions/runs/39337029515) tested the same parameters with N = 8 (projected dimension 6) at two resolutions:

| Configuration | Status | Interval LDL pivot lower bounds |
|---|---|---|
| step = 0.02, 128 series terms, 192 bits, 8,000 cells | inconclusive | 0.2149252225; 0.6591348762; 0.3643216324; 0.983522902995; -0.2238731886 |
| step = 0.01, 256 series terms, 256 bits, 16,000 cells | certified positive definite finite matrix | 0.2263857073; 0.7849070623; 0.4694446513; 1.0920803137; 0.6082769443; 1.2479011409 |

The baseline run became inconclusive at the fifth pivot. This is **not** a certified negative direction: interval LDL returns “inconclusive” whenever a pivot interval is not provably positive. Refining the mesh and increasing precision/series terms produced six strictly positive pivot lower bounds. This comparison is evidence of sensitivity of the *certificate procedure* to enclosure width, not evidence for or against global Weil positivity.

Artifact: [N=8 and subsequent interval reports](https://github.com/jbarumby/scicalc/actions/runs/39337029515/artifacts/12226434916).

A further controlled step to N = 10 has been added to the workflow, with both the baseline and refined configurations. Its result must be recorded only after the corresponding CI run finishes.

## Controlled dimension increase: N = 10

The follow-up workflow [#39337323734](https://github.com/jbarumby/scicalc/actions/runs/39337323734) completed successfully. For N = 10 (projected dimension 8):

| Configuration | Status | Interval LDL pivot lower bounds |
|---|---|---|
| step = 0.02, 128 series terms, 192 bits, 8,000 cells | inconclusive | first five computed lower bounds: 0.2149252225; 0.6591348762; 0.3643216324; 0.9835229030; -0.2238731886 |
| step = 0.01, 256 series terms, 256 bits, 16,000 cells | certified positive definite finite matrix | 0.2263857073; 0.7849070623; 0.4694446513; 1.0920803137; 0.6082769443; 1.2479011409; 0.5969515469; 1.4503530679 |

The refined N = 10 run therefore certifies positive definiteness of the finite projected matrix of dimension 8. The baseline's inconclusive result is not a negative-eigenvalue certificate. The first five baseline pivot lower bounds repeat those seen at N = 8; the extra coordinates only appear after the earlier pivots.

Artifact (N=8 and N=10, baseline and refined): [download from GitHub Actions](https://github.com/jbarumby/scicalc/actions/runs/39337323734/artifacts/12227530705).

This is a second finite-dimensional certificate under the same support/cutoff settings, not evidence that the full Weil quadratic form is positive on all admissible test functions.

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
