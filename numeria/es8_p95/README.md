# NUMERIA ES8-P95 — real Arb pilot for W12

Research-only branch; no change to main. Inherits P94 outer GL12,
Cauchy/Taylor bounds and original interval Legendre nodes.

The **push workflow** runs Python 3.11 and python-flint 0.9.0 in GitHub
Actions. It first checks that the 12-point rule encloses polynomial moments
through order 23; then runs exactly **one** Fourier segment (segment 0,
r in [0,2.5]), for the center pair (1,2).

**A successful pilot is NOT a certificate for W12**, which needs all 240
segments, together with the P89 prime/pole terms and infinite Fourier tail.
The workflow never labels a partial run as a full W12 proof.

The original source is kept in p94_arb_fourier_gl.py and
originais/gauss_legendre_intervalar_flint.py. See workflow logs and its
ES8P95_W12_PILOT artifact (result and log on success).

Next gate after a green pilot: eight shards covering all 240 segments,
hash-checked aggregation and cross comparison with the P89 interval.
