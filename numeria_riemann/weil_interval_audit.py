"""Finite-dimensional NUMERIA–RIEMANN Weil interval audit.

Uses python-flint/Arb for enclosures. A positive interval LDL result concerns
only the stated finite constrained matrix, not global Weil positivity or RH.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from flint import arb, ctx

def _a(value):
    """Parse a user decimal as an Arb ball, avoiding binary-float constants."""
    return arb(str(value))


def _sinc(x):
    """Enclose sin(x)/x, including intervals containing zero."""
    if not x.contains(0):
        return x.sin() / x
    if not x.abs_upper() < arb("0.1"):
        raise ValueError("interval crossing sinc's removable singularity is too wide")
    total = arb(1)
    term = arb(1)
    for n in range(1, 9):
        term = term * x * x / ((2 * n) * (2 * n + 1))
        total += term if n % 2 == 0 else -term
    # For |x| <= 0.1, the alternating-series remainder after x^16/17!
    # is bounded by 0.1^18/19! < 1e-34.
    return total + arb("0 +/- 1e-34")


def _omega(r, terms):
    """Enclose Re digamma(1/4+i*r/2)-log(pi) for r>=0.

    Reflection and duplication give
      Re psi(1/4+i*r/2) = Re psi(1/2+i*r) - log(2) - pi/(2 cosh(pi*r)).
    Also Re psi(1/2+i*r) = -gamma - 2 log(2) + sum_{n>=0}
      r^2 / ((n+1/2)*((n+1/2)^2+r^2)).
    The summand decreases in n. The tail sum is between its integral
    0.5*log(1+r^2/(M+1/2)^2) and that integral plus the first omitted term.
    """
    pi = arb.pi()
    r2 = r * r
    total = -arb.const_euler() - 3 * arb(2).log() - pi.log()
    for n in range(terms):
        y = arb(n) + arb("0.5")
        total += r2 / (y * (y * y + r2))
    y = arb(terms) + arb("0.5")
    integral_tail = arb("0.5") * (1 + r2 / (y * y)).log()
    first_omitted = r2 / (y * (y * y + r2))
    # Enclose [integral_tail, integral_tail + first_omitted].
    sum_tail = integral_tail + first_omitted / 2 + arb(0, (first_omitted / 2).upper())
    q = (-pi * r).exp()
    sech_term = pi * q / (1 + q * q)
    return total + sum_tail - sech_term


def _fourier_amplitude(r, L, k):
    """Real amplitude in F_k(r)=i*exp(-i*k*pi/2)*amplitude."""
    pi = arb.pi()
    w = arb(k) * pi / (2 * L)
    return 2 * L.sqrt() * w / (w + r) * _sinc(L * (r - w))


def _archimedean_matrix(L, N, R, step, terms):
    """Rigorous entrywise enclosure of (1/pi) integral_0^R Re(F_i conj(F_j))*Omega dr."""
    count = int(math.ceil(float(R.upper()) / float(step.lower())))
    width = R / count
    pi = arb.pi()
    G = [[arb(0) for _ in range(N)] for _ in range(N)]
    for cell in range(count):
        center = (arb(2 * cell + 1) * R) / (2 * count)
        r = center + arb(0, width / 2)
        omega = _omega(r, terms)
        amplitudes = [_fourier_amplitude(r, L, k) for k in range(1, N + 1)]
        for i in range(N):
            for j in range(i, N):
                difference = (i + 1) - (j + 1)
                phase = 0 if difference % 2 else (-1) ** (difference // 2)
                if phase == 0:
                    contribution = arb(0)
                else:
                    integrand = phase * amplitudes[i] * amplitudes[j] * omega / pi
                    contribution = width * integrand
                G[i][j] += contribution
                if i != j:
                    G[j][i] += contribution
    # Conservative analytic tail bound, r >= max(1, 2*w_max):
    # |F_k(r)| <= 4*w_k/(sqrt(L)*r^2),
    # |Omega(r)| <= 8 + log(3*r) (see TAIL_PROOF.md).
    # Integrating (8+log(3*r))/r^4 from R to infinity gives
    # (8+log(3*R)+1/3)/(3*R^3).
    wmax = arb(N) * pi / (2 * L)
