# NUMERIA–RIEMANN P77: rigorous K1 and B1 reconstruction

**Scientific state:** implementation prepared, rigorous evaluation in FLINT/Arb **pending CI replay**. The independent nonrigorous diagnostic with `mpmath` (45 dps) gave positive margins for all 16 rows, minimum about `8.571548196e-5` for the first row. Those numbers are **not certificates** until the Arb job passes.

**Continuity:** based on PR #76, preserving PR #75, PR #74 and parallel `audit/p75-journal-only-baseline`. Nothing is merged automatically.

Primary source: [Chourasiya–Simonič v2](https://arxiv.org/html/2507.15184v2), Theorem 1 and eqs. (10), (35)–(40), Corollary 5, Lemmas 5 and 6. The numerical quadrature is *not* a proof; only `acb.integral` enclosures are accepted.

## Rigorous quadrature and explicit tail

Equation (37) defines an integral

\[
I=\int_0^\infty e^{-u}u^{1/6}\left(1+\frac{\log(1+cu)}{2L}\right)^{4/3}(C_1\sqrt{cu}+C_2)^{2/3}\,du,
\quad c=3/(4\kappa),\quad L=\log T_0.
\]

The substitution `u=x^6` gives, on `[0,64]`, the integral over **`x ∈ [0,2]`** of

`6*x^6*exp(-x^6)*(1+log(1+c*x^6)/(2*L))^(4/3)*(C1*sqrt(c)*x^3+C2)^(2/3)`.

This removes the fractional power at the endpoint x=0 and is analytic near the real interval (with strictly positive power bases on the real segment). The integrator uses `acb.integral` with branch-aware `analytic=True` in the callback for both log and fractional powers. If a nonfinite or nonreal enclosure occurs, the gate fails closed.

For `u≥64≥1`, `log(1+cu)≤cu`, `1+au≤(1+a)u`, and since the bracket is `≥1` its exponent `4/3` can be majorized by `2`. Also `(C1√(cu)+C2)^(2/3) ≤ (1+C1√c+C2) u^(1/3)` by `y^(2/3)≤1+y`. Thus the original integrand is bounded above by

`D * exp(-u) * u^(5/2) ≤ D * exp(-u) * u^3`, where `D=(1+c/(2L))^2*(1+C1*sqrt(c)+C2)`.

Hence **without numerical truncation assumptions**:

`0 ≤ ∫_64^∞ integrand du ≤ D*exp(-64)*(64^3+3*64^2+6*64+6)`.

All these quantities are represented as rigorous Arb intervals with outward rounding.

## Exact choice of A0

The original table prints `log10(A0)` to 2 decimals, which does not specify the optimizer used by the authors. We instead **define our own admissible parameter** for every row, `A0=10^(the exact printed decimal)`; the theorem permits any `A0≥10^10`. It is rigorously checked that `10^10 < A0 < T0=3·10^12`. This avoids assuming anything about the missing digits of an optimizer.

For every row, the code rebuilds the Corollary 5 fourth-moment majorants `a1(A0),a2(A0)`, Lemma 5 `C1,C2`, Lemma 6 `C3,C4`, then `K1,K2,K3,L0,L1,L2`, the Theorem 1 `A1,A2`, and Corollary 1 `B1`. It compares exact published rational decimals to the **upper Arb endpoint** `B1_upper`, failing if `inf(B1_published-B1_upper) <= 0`.

The audit also needs a separate review of the hypotheses and external theorems; even a passing B1 check does **not** prove RH globally or automatically validate all downstream P75 conclusions.

## Reproduction

```sh
python -m pip install python-flint==0.9.0 mpmath==1.3.0 scipy
python -m unittest NUMERIA_P77_B1.test_math -v
python -m unittest NUMERIA_P77_B1.test_flint -v
python -O -m unittest NUMERIA_P77_B1.test_flint -v
python -m NUMERIA_P77_B1.audit_b1 --precision 192 --out /tmp/P77_B1_192.json
python -m NUMERIA_P77_B1.audit_b1 --precision 256 --out /tmp/P77_B1_256.json
```

The JSON certificates contain exact dyadic upper bounds and strict lower margins for each printed B1. The CI job must fail closed if any margin is nonpositive. `mpmath` is for diagnostics only; publication of any positive theorem is contingent upon the Arb execution **and** correct use of source hypotheses.