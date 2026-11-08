# NUMERIA ES8-P100 — determinant decision with tighter Arb integration

Independent direct P92-normalized 2x2 matrix via ES8-P99's
`certifica_matriz_2x2_P98.py` and `W12_arb_direto_P98.py`.

- eps = 1/2048 (best previous choice)
- bits = 128 and 256 in two separate jobs
- `--tol 1e-8` (P99 1e-7)
- `--eval-limit 60000` (P99 15000)
- `--depth-limit 40` (P99 24)
- `python-flint==0.9.0`, 35-minute job upper bound

The workflow preserves complete JSON outputs even for correctly
reported `INCONCLUSIVE_DETERMINANT_NOT_POSITIVE` results. The
`auditar_uma_execucao.py` recomputes the Sylvester criterion with exact
rational arithmetic. It rejects missing or nonfinite enclosures.

**Pass/fail policy:** A green CI status means the computation and audit
completed, **not** that the determinant was certified. A validated
`INTERVAL_POSITIVE_DEFINITE_CONDITIONAL` state requires W11_lower>0,
W22_lower>0, and a positive exact lower bound on the determinant.

Results require review of formula identification and nested `acb.integral`
assumptions. This is about one finite 2x2 Weil subspace and **not** RH.

The 50-mode experiment is preserved and is not changed in this branch.
