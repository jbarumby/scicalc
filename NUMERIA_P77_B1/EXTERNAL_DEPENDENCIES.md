# P77 — External analytic dependencies and P78 priorities

## Independent coefficient audit versus upstream theorem

The GitHub Actions execution https://github.com/jbarumby/scicalc/actions/runs/45021396110 verified the 16 printed B1 coefficient comparisons under FLINT/Arb at 192 and 256 bits. See CI_RESULT.md and the archived interval JSONs.

This verifies the reconstruction of the numerical coefficient formulas from Chourasiya–Simonič v2, **not** an independent proof of every analytic inequality used to establish the underlying zero-density theorem.

## Specific unresolved source dependency: Lemma 6

Source: https://arxiv.org/html/2507.15184v2 (Lemma 6 and its footnote, around the definitions of C3 and C4). The paper uses the parameter b4 = 0.54 in a Dirichlet-coefficient bound and credits a crucial inequality to a private communication from Olivier Ramaré; it mentions anticipated joint work with Sebastián Zúñiga Alterman.

The P77 arithmetic accepts the displayed Lemma 6 estimate and faithfully reconstructs C3, C4, K3, L1, L2, A1, A2 and B1. It does not independently prove the private-communication inequality, nor has a published source with precisely matching hypotheses yet been authenticated. This is an **open dependency gate**, not evidence that the stated inequality is false.

## Priority P78

1. Locate and inspect a public proof that supplies precisely the Lemma 6 bound behind b4 = 0.54, including the same definitions, domain and constants. Do not substitute a superficially similar inequality.
2. If no proof can be authenticated, derive an independently documented conservative alternative, with explicit mathematical hypotheses and rational/Arb error controls; reconstruct all 16 B1 upper bounds and re-execute the P75 window from t=18 to t=79.
3. Independently audit the chain Suzuki representation / Platt–Trudgian finite-height RH / Bellotti–Wong total count / Chourasiya–Simonič v2 density assumptions, including endpoint multiplicities, zero-counting conventions and symmetry factors.
4. Keep the parallel journal-only baseline in audit/p75-journal-only-baseline separate; preserve PRs #71–#77 and original certificates and do not merge them automatically.
5. Do not attempt to declare any extension beyond t=79 until the source-level gate and finite-window replay have been reproduced.

**Scientific scope:** the finite-window positivity checks and the B1 coefficient inequalities do not prove the global Riemann hypothesis.
