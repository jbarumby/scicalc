# NUMERIA–RIEMANN · P73 — randomized discovery and deterministic zero-side certification

**Working base:** [P72](https://github.com/jbarumby/scicalc/pull/72), based on [P71](https://github.com/jbarumby/scicalc/pull/71).  
**Date:** 2026-11-21. **Status:** computable *finite-window* certificate **conditional on verified external lemmas and provenance of 20 FLINT/Arb zero enclosures**; *not* proof of global RH.

## Core finding

The P72 verifier uses a **single worst-case tail correction at final T**, subtracting it from *all* values of the first 20 zero terms. P73 instead certifies each finite interval with a **local, monotone upper tail bound** and adaptive rational subdivision. Holding the same zeros and Platt–Trudgian height `H=3×10¹²`, its exact-rational certificate reaches:

`18 <= t <= 46.6`, with **565 covering intervals**, **566** rational endpoint evaluations (see JSON), strict lower bound **+0.0000419846671888** (rounded display).

The rigorous implication is this: *assuming the imported Suzuki spectral identity, Platt–Trudgian finite RH, Bellotti–Wong counting bound, and that the 20 supplied dyadic intervals enclose the first 20 positive zero ordinates*, the rational cell calculations imply `Psi(t) > 0` on the stated interval. Link this to earlier finite certificates to obtain coverage starting at t=0 only after verifying endpoints and overlap. The current code does not reprove these external theorems or authenticate arbitrary input zero JSON files by itself.

### Mathematics of the local-cell proof

For certified first `N` positive zero ordinates (`γ_j` real), define

`S_N(t) = 2 Σ_{j=1}^N (1 - cos(γ_j t))/γ_j²`.

The unconditional spectral expression of Suzuki is a sum over all zeros, which may not lie on the critical line. Platt–Trudgian show that all zeros with `|Im(ρ)| <= H` are on the critical line; therefore omitted terms from that finite height range are nonnegative, but zeros above `H` require absolute bounds. P72 combines Suzuki's identity and Bellotti–Wong to give

`Psi(t) >= S_N(t) - U(t)` for `t>=0`,  
`U(t) = (3 + exp(t/2)) (5/H + 38/H²)`.

For any cell `[a,b]` of length `h=b-a`, positivity follows from the *rational* lower bound

`L[a,b] = min(S_N^lo(a), S_N^lo(b)) - N*h²/4 - U^up(b)`.

Here `U^up(b)` is computed using an exact-rational majorant for `exp` at `ceil(b/bucket)*bucket` (bucket `1/20` by default). Since `|S_N''|<=2N`, the secant-interpolation error is at most `N*h²/4`. A cell is accepted iff `L[a,b]>0`; otherwise it is **bisected deterministically**, until a limit or a failure. Accepted cells are sorted, their union is checked to cover the entire interval without gaps, and a canonical exact-rational SHA-256 transcript is emitted. No randomized quantity enters that test.

The previous *global-tail* margin at `t=44.35` was around `0.000176565`; the P73 margin at `46.6` is positive. The uniform minimum in P73 is a *certificate-method margin*, not the actual minimum of `Psi`.

### Distinct non-proof discovery experiment

`experiment.py` samples `t` values in `44.35 <= t <= 47.48`, using the 20 zero midpoints in IEEE-754 double precision and a monotone tail proxy. It compares equal budgets of uniform sampling, seeded random sampling (seed 73), and a width-weighted adaptive search favoring small proxy values.

At **512 samples each** (see `discovery_512.json`):

| Method | Score evaluations | Negative bound-proxy samples | First negative sample (approx.) |
|---|---:|---:|---:|
| Uniform grid | 512 | 56 | 46.6347 |
| Seeded random | 512 | 48 | 46.6434 |
| Adaptive heuristic | 512 | 80 | 46.6364 |

`score(t) < 0` only says the **20-zero lower-bound method is inconclusive at t**. It is *not* an instance of `Psi(t)<0` and is *not* a counterexample to RH. The experiment is limited to one seed and one budget; it does **not** establish that the adaptive search is universally faster. The randomized strategy is a discovery method, not a proof or a stochastic generation of actual primes. Statistical models of random primes or random matrices are not substituted for zeta zeros in the certificate.

## Reproduction

Execute from repository root. The P72 JSON is pinned by SHA-256 and should first be **regenerated with FLINT/Arb** and matched byte-for-byte in CI. Merely copying a file and writing `zero_certification: arb_enclosures` is not proof of provenance.

```bash
python -m pip install python-flint==0.9.0 mpmath==1.3.0
python NUMERIA_P72_ZERO_SIDE/export_zeros_flint.py --n 20 --precision 192 --out /tmp/p73_zeros_replay.json
cmp /tmp/p73_zeros_replay.json NUMERIA_P72_ZERO_SIDE/arb_zeros_20.json

python NUMERIA_P73_ADAPTIVE/certify_adaptive.py \
  --zeros NUMERIA_P72_ZERO_SIDE/arb_zeros_20.json \
  --expected-zeros-sha256 c5993b75a6db7803251655955859cd935f4381b2c0bc20fd77b9f0258bdeddab \
  --to 46.6 --out /tmp/p73_replayed.json
cmp /tmp/p73_replayed.json NUMERIA_P73_ADAPTIVE/certificate_20_t46.6.json

python NUMERIA_P73_ADAPTIVE/experiment.py \
  --zeros NUMERIA_P72_ZERO_SIDE/arb_zeros_20.json \
  --budget 512 --seed 73 --out /tmp/p73_discovery.json
(cd NUMERIA_P73_ADAPTIVE && python -m unittest test_p73 -v)
(cd NUMERIA_P73_ADAPTIVE && python -O -m unittest test_p73 -v)
```

The certificate can be regenerated with `python -O` and compared byte-for-byte. The workflow automates these gates, including checking the 20 zero enclosures using the proof-producing FLINT routine before accepting the certificate.

## Next scientific gates

1. Review **Suzuki eq. (1.3)** and verify normalization, zero-pairing/symmetries, and convergence; don't let `certified_positive` be inferred from a JSON metadata string.
2. Independently check the **Bellotti–Wong** constants and the retained `-N(H)/H²` boundary term; cross-check **Platt–Trudgian** `H=3×10¹²` theorem assumptions.
3. Replay P73 in a second environment at 192/256 bits. 256-bit intervals will not produce byte-identical files; independently compare containment and the sign of the certificate.
4. If extending beyond approximately `46.63` stalls, regenerate additional *real* zeros by Arb (50–3,000); use local-cell certificates and optimize the tail. Do not interpret a negative lower bound as a negative Suzuki function.
5. Consider stronger zero-tail estimates and adaptive precision. Randomized prime simulations (Cramér/Granville) or CUE random-matrix models (Keating–Snaith) may help **heuristic prioritization**, but cannot replace certified zeta zero ordinates.

## Literature

- Masatoshi Suzuki, *Aspects of the screw function corresponding to the Riemann zeta-function*, J. London Math. Soc. (2023), [doi:10.1112/jlms.12785](https://doi.org/10.1112/jlms.12785).
- D. Platt and T. Trudgian, *The Riemann hypothesis is true up to 3·10¹²*, Bull. London Math. Soc. (2021), [doi:10.1112/blms.12460](https://doi.org/10.1112/blms.12460).
- C. Bellotti and P.-J. Wong, *Improved estimates for the argument and zero-counting function of the Riemann zeta-function*, Math. Comp. (2025), [doi:10.1090/mcom/4133](https://doi.org/10.1090/mcom/4133). Use the **published** `8.08344`, not the changed constant in early arXiv abstracts.
- J. P. Keating and N. C. Snaith, *Random Matrix Theory and ζ(1/2+it)*, Comm. Math. Phys. 214 (2000), [doi:10.1007/s002200000261](https://doi.org/10.1007/s002200000261), **heuristic background only**.
- A. Granville, *Harald Cramér and the Distribution of Prime Numbers*, Scand. Actuarial J. (1995), [article PDF](https://dms.umontreal.ca/~andrew/PDF/cramer.pdf), **random prime models; not a proof input**.

**Critical limitation:** Finite verification, even very large, cannot prove global RH. A finite positive `Psi(t)` interval is much weaker than positivity of Suzuki's full Hermitian kernel for all parameters.
