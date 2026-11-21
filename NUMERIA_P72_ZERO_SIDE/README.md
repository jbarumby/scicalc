# NUMERIA–RIEMANN · P72 — Zero-side extension, no extra prime powers

**Date:** 2026-11-21. **Baseline:** [PR #71](https://github.com/jbarumby/scicalc/pull/71), whose certified arithmetic coverage ends at `t=log(10^8)`; and the local 3-zero rational audit through `t=34.5`.

## Status: conditional calculation, not yet a new rigorously certified 20-zero theorem

The **exact-rational verifier** and its continuity/tail arguments have run. The **20 input zero enclosures have NOT been validated locally by Arb**: `preview_zeros.json` contains high-precision mpmath approximations surrounded by proposed intervals, not certified intervals. Therefore every positive result below is **conditional on those intervals containing the corresponding zeta zero ordinates**. The program explicitly reports `certified_positive: false` for these preview calculations.

The separate `export_zeros_flint.py` script calls the proof-producing `acb.zeta_zeros()` and exports exact *dyadic* interval endpoints. Running it with `python-flint==0.9.0` and rerunning the rational verifier on its output is the next mandatory gate. This environment lacks python-flint and cannot download it. Until that gate passes, do not claim that the 20-zero result is certified.

## Executed tests

| Zeros retained | Interval for t | Grid points | Lower bound for min S_N on grid | Curvature loss | Tail upper bound | **Conditional** uniform lower |
|---:|---:|---:|---:|---:|---:|---:|
| 3 | 18–34.5 | 3,301 | 0.000516561422957 | 0.00001875 | 0.0000516926271243 | **+0.000446118795833** |
| 10 | 18–40 | 4,401 | 0.00531855386899 | 0.0000625 | 0.000808608664018 | **+0.00444744520497** |
| 20 | 18–40 | 4,401 | 0.00741908812479 | 0.000125 | 0.000808608664018 | **+0.00648547946077** |
| 20 | 18–44.35 | 5,271 | 0.00741908812479 | 0.000125 | 0.00711752311308 | **+0.000176565011706** |
| 20 | 18–44.40 | 5,281 | 0.00741908812479 | 0.000125 | 0.00729770406840 | **−0.00000361594361** |
| 20 | 18–47.48 | 5,897 | 0.00741908812479 | 0.000125 | 0.0340407993797 | **−0.0267467112549** |

The final two rows are **inconclusive**, not counterexamples to Suzuki positivity. The positive `t=44.35` row proves an exact implication from input zero enclosures; it will become a certified statement only after the Arb zero generator and full replay pass. No prime-power enumeration is used here.

All grid lines are evaluated with fixed-point integer interval arithmetic at scale 10^75, a rational Machin interval for pi, range reduction to [-4,4], Taylor cosine with verified remainder <10^-50, and outward rounding at every operation. Each grid has a SHA-256 digest in its JSON file.

## Tail proof (the key improvement)

Let `H=3*10^12` and `N(u)` be the Riemann–von Mangoldt zero-counting function, counting zeros of positive imaginary part with multiplicity. Bellotti–Wong (2025) give, for `u>=e`, 

`|N(u) - u/(2*pi)*log(u/(2*pi*e))| <= 0.10076*log(u)+0.24460*log(log(u))+8.08344.`

Define `M(u)=u/(2*pi)*log(u/(2*pi*e))`, `E(u)=(0.10076+0.24460)*log(u)+8.08344`; we used `log(log(u)) <= log(u)`.

**Retaining the boundary term** in the Stieltjes integration gives, exactly,

`sum_{Im(rho)>H} 1/Im(rho)^2 = 2*integral_H^infty N(u)/u^3 du - N(H)/H^2`.

Using `N(u)<=M(u)+E(u)` and `N(H)>=M(H)-E(H)`,

`tail <= (log(H/(2*pi))+1)/(2*pi*H) + [0.34536*(2*log(H)+1/2)+2*8.08344]/H^2`.

Since `log H<30`, `pi>3`, `log(2*pi)>1`, and the bracketed constant is below `38`, we obtain

`sum_{Im(rho)>H} 1/Im(rho)^2 < 5/H + 38/H^2`.

Suzuki (2023), formula (1.3), is an **unconditional** expansion over the possibly complex zero ordinates. Platt–Trudgian rigorously verified RH below `H`, so every omitted zero below `H` contributes **nonnegatively**. Above `H`, `|Im(gamma)|<1/2`, `|gamma|>=Im(rho)`, and

`sum_{|Re(gamma)|>H} |(1-cos(gamma*t))/gamma^2| <= 2*(1+cosh(t/2))*(5/H + 38/H^2)`.

For `0<=t<=T`, the right side is strictly below `(3+exp(T/2))*(5/H+38/H^2)`; we enclose `exp(T/2)` from above by a **rational** Taylor-series/geometric-tail bound and round it upwards to the next 10^-9.

## Curvature interpolation

For the `N` retained positive zero ordinates `gamma_j`,

`S_N(t)=2*sum_{j=1}^N (1-cos(gamma_j*t))/gamma_j^2`, so `S_N''(t)=2*sum_j cos(gamma_j*t)<=2N` for real ordinates.

For any cell `[a,b]` of length `h`, the comparison of `S_N` with its secant line implies

`S_N(t) >= min(S_N(a),S_N(b)) - N*h*h/4`.

This is the certified **entire interval** bound used in the program; grid sampling alone would not suffice.

## Reproduction

With `python-flint==0.9.0` available:

```bash
python export_zeros_flint.py --n 20 --precision 192 --out arb_zeros_20.json
python verify_spectral_rational.py --zeros arb_zeros_20.json --n 10 --to 40 --out arb_10_t40.json
python verify_spectral_rational.py --zeros arb_zeros_20.json --n 20 --to 44.35 --out arb_20_t44.35.json
python verify_spectral_rational.py --zeros arb_zeros_20.json --n 20 --to 47.48 --out arb_20_t47.48.json
python -m unittest test_p72 -v
python -O -m unittest test_p72 -v
```

The first two positive runs with Arb-generated intervals should give `certified_positive: true` **if** their computed uniform bounds remain positive. This result is not asserted before running those commands. The `t=47.48` run is currently inconclusive because the bound for zeros above H is larger than the retained 20-zero contribution.

For a GitHub Actions replay, the package also contains `github_workflow_numeria_p72.yml` (copy it to `.github/workflows/numeria-p72-zero-side.yml`). No file in the actual GitHub repository has been modified by this local preparation.

## Sources and remaining gates

- Suzuki 2023, eq. (1.3): https://doi.org/10.1112/jlms.12785
- Platt–Trudgian, rigorous RH up to height 3e12: https://doi.org/10.1112/blms.12460
- Bellotti–Wong 2025, explicit N(T): https://doi.org/10.1090/mcom/4133
- FLINT verified zero isolation + Turing counting: https://flintlib.org/doc/acb_dirichlet/zeta-zeros.html
- python-flint `acb.zeta_zeros`: https://python-flint.readthedocs.io/en/latest/acb.html
- Independent published zero-side benchmark (3000 zeros, `t<=47.48`): https://zenodo.org/records/23481975

**Gate A:** independently generate and validate FLINT/Arb zeros 1–20, recording exact dyadic intervals, version and hash.  
**Gate B:** rerun rational verifier in two environments; compare SHA-256 and margins.  
**Gate C:** to approach `t=47.48`, retain substantially more positive ordinates (e.g. thousands) and optimize the high-zero tail; 20 zero terms plus this error bound cannot close it.  
**Gate D:** formalize the analytic zero-counting and boundary-term lemmas separately from the numerical replay.

A finite positive window is never a proof of global RH, nor is this certificate for the full Weil positivity form.
