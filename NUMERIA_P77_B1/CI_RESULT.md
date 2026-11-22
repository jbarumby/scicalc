# P77 — Arb B1 certificate results (GitHub Actions)

**Run:** https://github.com/jbarumby/scicalc/actions/runs/45021396110  
**Artifact:** `NUMERIA-P77-B1-certificates`, GitHub artifact ID `14754026929` (`b1-192.json`, `b1-256.json`).  
**Run state:** completed, success. Executed 2026-11-22.

## Rigorous theorem-to-table gate

At **192 and 256 bits** all 16 rows passed the strict Arb inequality

`inf(B1_published - B1_computed_Arb_upper) > 0`.

The smallest certified positive lower margin occurs in row 1 (`σ∈[1/2,17/32]`):

- 192-bit exact dyadic lower margin: `16813900017557291983815347421280044986879157385642227 × 2^-187`;
- 256-bit exact dyadic lower margin: `1240646842019279607076415312874891156210052344595007339131031016833389825 × 2^-253`;
- decimal orientation only: approximately **8.571548196 × 10^-5**.

The CI independently checked the SHA-256 source manifest, ran four algebraic tests, performed the 16-row rigorous certificate at both precisions and ran three additional Arb regressions in normal Python and with `-O`. The certificate JSONs were preserved as a downloadable workflow artifact.

### Reproduction method

The original improper integral in Eq. (37) was split at `u=64`. The substitution `u=x^6` regularized the finite part `x∈[0,2]`, which was evaluated using `acb.integral` with analytic branch checking. The infinite tail is bounded by the explicitly proved elementary majorant in `README.md`. Values of `A0` were specified as exact `10^(printed log10 A0)`, not taken from unknown rounding of an optimizer.

### Boundaries of what this certifies

This closes the **numerical coefficient reconstruction of B1** for the 16 rows from Chourasiya–Simonič v2, complementing the B2/B3 gates in P76. It does not independently prove every theorem quoted by that preprint or all analytic dependencies used by P75, and does not itself rerun the P75 t≤79 certificate with an alternate tail source. It is **not** proof of the global Riemann hypothesis. Keep PR #77 draft pending mathematical peer review.
