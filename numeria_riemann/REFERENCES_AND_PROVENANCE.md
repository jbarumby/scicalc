# External references and integration boundaries

This module is an original, isolated integration layer. It does not copy
source code from the projects below. Recheck upstream license files and
versions before adapting any implementation.

## Numerical backends

- **mpmath** — https://github.com/mpmath/mpmath
  Used by the existing prototype for arbitrary-precision evaluation and
  exploratory computation of zeta zeros. Upstream documentation describes
  `zetazero`, `zeta`, `siegelz`, and `nzeros`.
- **python-flint** — https://github.com/flintlib/python-flint
  Optional second numerical implementation. Its Arb/Acb ball arithmetic can
  provide enclosures; a midpoint comparison alone is not a certificate.
  Upstream repository reports the MIT license; preserve notices if code is
  ever adapted or copied.

## Independent zero-finding implementation reviewed

- **Riemann-Zeta-Zero-Finding-Suite** —
  https://github.com/error-wtf/Riemann-Zeta-Zero-Finding-Suite
  Reviewed for architecture and validation ideas (separating approximate
  Riemann–Siegel main-sum evaluation from complete evaluations and interval
  diagnostics). Its repository licensing notice says software artifacts are
  ALL RIGHTS RESERVED absent separate written permission. Therefore no code
  from it is copied into NUMERIA. Any future code reuse requires permission.
- **zeta-lab** — https://github.com/teal-sea/zeta-lab
  Candidate for future independent audit of explicit-formula workflows.
  Verify current source licenses and validate claims independently before
  reusing code or treating any stated mathematical result as established.

## Interpretation rules

1. Numerical residuals, cross-backend agreement, and finite truncations are
   evidence about computations, not proof that a computed value is an exact
   zero.
2. A certified zero claim needs a mathematically justified enclosure and a
   validated root-count/isolation argument.
3. Finite explicit-formula sums can oscillate and need careful handling at
   prime powers. The current baseline uses the inclusive staircase convention.
4. This repository does not establish the Riemann Hypothesis, global Weil
   positivity, or a Hilbert–Pólya operator construction.
5. Record software versions, precision, inputs, outputs, and hashes for
   benchmark runs. Never describe a test as executed unless a runtime result
   exists.
