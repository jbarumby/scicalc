# Codex handoff — NUMERIA P77

1. Start from PR #77, branch `numeria-p77-cs-v2-b1-quadrature`, based on PR #76; **do not overwrite** the P75 or parallel `audit/p75-journal-only-baseline`.
2. Inspect the original Chourasiya–Simonič v2 eqs. (10), (35)–(40); independently verify every A1/A2 factor, every C1…C4 term and the exact parameter choices `A0=10^printed_log10A0`.
3. Run FLINT 0.9.0 192/256 bits and verify all 16 `B1_published > B1_Arb_upper` inequalities. Inspect quadrature convergence and the explicit `u>=64` majorant. If CI fails, publish the failure and fix only the proved bug; do not force success.
4. Replay all of P76 B2/B3, P75 500-zero t<=79 and `audit/p75-journal-only-baseline`. Be explicit whether the table's 16 coefficients have mathematically audited provenance and whether all external lemmas hold.
5. Only then consider extending to t>79, in an independent PR with SHA256, JSON intervals and reproducible CI. Never assert a global RH proof.