# NUMERIA–RIEMANN — corrected finite Weil interval audit

**Audit date:** 2026-10-20  
**Status:** exploratory certificate attempt; no positive-definiteness certificate after independent normalization review.

## 1. Formula review: prime-overlap normalization

The Fourier amplitudes and the pole-constraint vectors are consistent with the unit-normalized sine basis
\[
g_k(x)=\frac{\sin(w_kx)}{\sqrt L},\qquad 0\le x\le 2L,\quad w_k=\frac{k\pi}{2L}.
\]
Indeed, \(L^{-1}\int_0^{2L}\sin^2(w_kx)\,dx=1\). For shift \(a=\log n\), the normalized overlap is
\[
C_{ij}(a)=\frac1L\int_0^{2L-a}\sin(w_i(x+a))\sin(w_jx)\,dx.
\]
Product-to-sum contributes a factor \(1/2\), so if `first` and `second` are the two cosine integrals, the overlap is \((\text{first}-\text{second})/(2L)\), not \((\text{first}-\text{second})/(4L)\). The old code therefore halved the prime contribution while the explicit-formula coefficient remained \(-2\log p/\sqrt n\).

The corrected implementation uses `(first - second) / (2 * L)`. A regression test checks that the zero-shift overlap matrix contains the identity Gram matrix. The older positive certificates are superseded.

## 2. Corrected CI sweep

Workflow run [#39345184865](https://github.com/jbarumby/scicalc/actions/runs/39345184865) completed with all four jobs successful, including the interval tests and all six matrix runs. The downloadable artifact is [here](https://github.com/jbarumby/scicalc/actions/runs/39345184865/artifacts/12229824716).

Parameters common to the sweep: \(L=0.45\), \(R=160\), active prime-power set \(\{2\}\). Baseline: step 0.02, 128 digamma-series terms, 192-bit Arb, 8,000 cells. Refined: step 0.01, 256 terms, 256-bit Arb, 16,000 cells.

| N | Projected dimension | Resolution | Status | LDL pivot lower bounds returned |
|---|---:|---|---|---|
| 6 | 4 | baseline | inconclusive | 0.1371285697; 0.7033747368; -0.8535492277 |
| 6 | 4 | refined | inconclusive | 0.1485890546; 0.8291504881; -0.4839289451 |
| 8 | 6 | baseline | inconclusive | 0.1371285697; 0.7033747368; -0.8535492277 |
| 8 | 6 | refined | inconclusive | 0.1485890546; 0.8291504881; -0.4839289451 |
| 10 | 8 | baseline | inconclusive | 0.1371285697; 0.7033747368; -0.8535492277 |
| 10 | 8 | refined | inconclusive | 0.1485890546; 0.8291504881; -0.4839289451 |

The LDL routine stops at the third pivot because its lower bound is not strictly positive. A negative lower bound alone does **not** certify a negative pivot: the upper bound is needed to establish that. The first three pivot bounds are identical across N in this parameterization because the leading projected coordinates use the same first five sine modes and the same two constraint columns. Refinement improves the third lower bound but does not certify positive definiteness.

**Correct conclusion:** no positive-definiteness certificate was obtained after correcting the normalization. These computations neither prove nor disprove positivity of the full Weil form or the Riemann Hypothesis.

## 3. Error and tail enclosure

The finite archimedean integral on \([0,R]\) is enclosed by Arb interval range sums. The digamma series tail is bounded by the integral test plus the first omitted term. The infinite archimedean tail uses
\[
|F_k(r)|\le\frac{4w_k}{\sqrt L\,r^2},\qquad
|\Omega(r)|\le \log r+\log\pi+\frac{4}{3r},
\]
and therefore
\[
|G_{ij}^{\rm tail}|\le
\frac{16w_iw_j}{\pi L}
\left[
\frac{\log R+\frac13+\log\pi}{3R^3}
+\frac{1}{3R^4}
\right].
\]
For the absolute \(\Omega\) bound, the corrected implementation now requires \(R\ge32\). This range can be justified from
\[
\Omega(r)=-\gamma-3\log2-\log\pi+
\sum_{m\ge0}\frac{r^2}{(m+\frac12)((m+\frac12)^2+r^2)}
-\frac{\pi}{2\cosh(\pi r)}
\]
and the integral-test bounds
\[
\frac12\log(1+4r^2)\le S(r)\le
\frac12\log(1+4r^2)+\frac{2r^2}{r^2+\frac14}.
\]
For \(r\ge32\), these bounds imply \(\Omega(r)>0\) and the stated, deliberately loose upper bound. All reported runs use \(R=160\), so they satisfy this cutoff requirement.

## 4. Reproduction

The CI pins `python-flint==0.9.0` and Python 3.11. It runs:

```bash
python -m unittest numeria_riemann.test_weil_interval_audit -v
python -m numeria_riemann.weil_interval_audit --L 0.45 --N 6 --R 160 --step 0.02 --series-terms 128 --precision 192 --out N6.json
python -m numeria_riemann.weil_interval_audit --L 0.45 --N 6 --R 160 --step 0.01 --series-terms 256 --precision 256 --out N6-fine.json
python -m numeria_riemann.weil_interval_audit --L 0.45 --N 8 --R 160 --step 0.02 --series-terms 128 --precision 192 --out N8.json
python -m numeria_riemann.weil_interval_audit --L 0.45 --N 8 --R 160 --step 0.01 --series-terms 256 --precision 256 --out N8-fine.json
python -m numeria_riemann.weil_interval_audit --L 0.45 --N 10 --R 160 --step 0.02 --series-terms 128 --precision 192 --out N10.json
python -m numeria_riemann.weil_interval_audit --L 0.45 --N 10 --R 160 --step 0.01 --series-terms 256 --precision 256 --out N10-fine.json
```

The correction is isolated on branch `numeria-riemann/overlap-normalization-fix-ci4`; no branch was merged into `main`. The original PR #1 remains unmerged. Before publication, the exact full pivot intervals and a second independent implementation of the prime-overlap formula should be retained alongside the JSON artifact.


## Final guarded rerun

The implementation now explicitly rejects cutoffs below R=32 for the stated absolute Omega tail bound, and a regression test checks this guard. The complete workflow was rerun at commit `7629b8c87cbcafefc05929ce5d9119f95fa27378`: [run #39345501884](https://github.com/jbarumby/scicalc/actions/runs/39345501884). All four jobs passed, including the new cutoff-guard test and all six N=6/N=8/N=10 runs. The pivot lower bounds and inconclusive status are unchanged from the table above. Final artifact: [N=6 through N=10 JSON reports](https://github.com/jbarumby/scicalc/actions/runs/39345501884/artifacts/12231150827).
