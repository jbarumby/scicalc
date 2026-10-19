"""Exploratory explicit-formula diagnostics for NUMERIA-RIEMANN.

This isolated module provides an exact finite prime-power baseline and a
truncated diagnostic version of Chebyshev's explicit formula. A finite zero
sum is not an exact identity: results depend on truncation and smoothing.

Dependency: mpmath.
"""
from __future__ import annotations

from math import floor, isfinite
from typing import Iterable, Sequence

import mpmath as mp


def von_mangoldt(n: int) -> float:
    """Return Lambda(n): log(p) for a prime power, otherwise zero."""
    if not isinstance(n, int) or isinstance(n, bool) or n < 1:
        raise ValueError("n must be a positive integer")
    if n == 1:
        return 0.0
    remaining = n
    prime_factor = None
    divisor = 2
    while divisor * divisor <= remaining:
        if remaining % divisor == 0:
            if prime_factor is not None:
                return 0.0
            prime_factor = divisor
            while remaining % divisor == 0:
                remaining //= divisor
        divisor = 3 if divisor == 2 else divisor + 2
    if remaining > 1:
        if prime_factor is not None:
            return 0.0
        prime_factor = remaining
    if prime_factor is not None:
        return float(mp.log(prime_factor))
    return 0.0


def psi_exact(x: float) -> float:
    """Compute psi(x)=sum_{n<=x} Lambda(n) directly for finite x >= 1."""
    if not isfinite(float(x)) or x < 1:
        raise ValueError("x must be finite and >= 1")
    total = mp.mpf("0")
    for n in range(2, floor(x) + 1):
        total += von_mangoldt(n)
    return float(total)


def first_zero_ordinates(count: int, dps: int = 30) -> list[mp.mpf]:
    """Return the first positive ordinates of non-trivial zeta zeros."""
    if not isinstance(count, int) or isinstance(count, bool) or count < 1:
        raise ValueError("count must be a positive integer")
    if not isinstance(dps, int) or dps < 15:
        raise ValueError("dps must be an integer >= 15")
    with mp.workdps(dps):
        return [mp.im(mp.zetazero(k)) for k in range(1, count + 1)]


def psi_from_zeros_truncated(
    x: float,
    zeros: Sequence[complex | float],
    *,
    dps: int = 30,
    assume_critical_line: bool = True,
) -> float:
    """Evaluate a finite explicit-formula diagnostic at x > 1.

    If assume_critical_line is true, each input is a positive ordinate gamma
    and rho=1/2+i*gamma is assumed; each represents a conjugate pair.
    Otherwise each input is an individual complex zero rho. Supply both
    members of each conjugate pair; their real contributions are summed
    without an extra factor of two.

    This finite truncation is not certified and may converge slowly. At prime
    powers the standard psi_0 formula uses the midpoint of the jump, while
    psi_exact returns the inclusive staircase value.
    """
    if not isfinite(float(x)) or x <= 1:
        raise ValueError("x must be finite and > 1")
    if not isinstance(dps, int) or dps < 15:
        raise ValueError("dps must be an integer >= 15")
    if len(zeros) == 0:
        raise ValueError("provide at least one zero or ordinate")

    with mp.workdps(dps):
        xx = mp.mpf(x)
        oscillation = mp.mpf("0")
        for item in zeros:
            if assume_critical_line:
                rho = mp.mpf("0.5") + 1j * mp.mpf(item)
                if abs(rho) == 0:
                    raise ValueError("zero location cannot be 0")
                # One positive ordinate represents the conjugate pair.
                oscillation += 2 * mp.re(mp.power(xx, rho) / rho)
            else:
                rho = mp.mpc(item)
                if abs(rho) == 0:
                    raise ValueError("zero location cannot be 0")
                # Here callers supply individual complex zeros. To reproduce
                # the real sum, include both members of each conjugate pair.
                # Do not double each term: that would count each pair twice.
                oscillation += mp.re(mp.power(xx, rho) / rho)
        correction = mp.log(2 * mp.pi) + mp.mpf("0.5") * mp.log1p(-1 / (xx * xx))
        return float(xx - oscillation - correction)


def convergence_report(
    x: float,
    zero_counts: Iterable[int] = (1, 5, 10, 20, 50),
    *,
    dps: int = 30,
) -> list[dict[str, float | int | str]]:
    """Compare truncations with exact psi(x), exposing the error explicitly."""
    counts = list(zero_counts)
    if not counts or any(not isinstance(n, int) or n < 1 for n in counts):
        raise ValueError("zero_counts must contain positive integers")
    ordinates = first_zero_ordinates(max(counts), dps=dps)
    exact = psi_exact(x)
    rows = []
    for count in counts:
        estimate = psi_from_zeros_truncated(x, ordinates[:count], dps=dps)
        rows.append({
            "x": float(x),
            "zero_count": count,
            "psi_exact_inclusive": exact,
            "psi_truncated_diagnostic": estimate,
            "absolute_error_vs_inclusive": abs(estimate - exact),
            "status": "diagnostic_only_not_certified",
        })
    return rows
