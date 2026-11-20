"""Independent Arb replay of the finite ES8-P249 prime operator and its symbol."""

from __future__ import annotations

import json
from fractions import Fraction
from operator import index

import numpy as np
import flint
from flint import arb, ctx

from .prime_operator import P236_HALF_WIDTH, P236_PRIME_POWERS, PrimeShiftOperator


def _dyadic(point):
    """Serialize an exact finite Arb endpoint without a decimal conversion."""
    if not point.is_finite() or not point.is_exact():
        raise ArithmeticError("expected a finite, exact Arb endpoint")
    mantissa, exponent = point.man_exp()
    return [str(mantissa), str(exponent)]


def _arb_prime_terms(powers):
    terms = []
    for power, prime in powers:
        power, prime = index(power), index(prime)
        if power < 2 or prime < 2:
            raise ValueError("prime and prime power must be at least two")
        terms.append((2 * arb(prime).log() / arb(power).sqrt(), arb(power).log()))
    return terms


def _real_binary64_vector(vector):
    vector = np.asarray(vector)
    if vector.ndim != 1 or vector.size == 0 or np.iscomplexobj(vector):
        raise ValueError("expected a nonempty, one-dimensional real vector")
    vector = vector.astype(np.float64)
    if not np.isfinite(vector).all():
        raise ValueError("expected finite binary64 coordinates")
    return vector


def prime_shift_lags_arb(cells, half_width=P236_HALF_WIDTH, powers=P236_PRIME_POWERS):
    """Enclose each Toeplitz lag by the geometric lengths of shifted overlaps.

    The action uses the unsigned sum of the two symmetric shifts for every
    supplied (prime power, prime) pair, on normalized cell indicators.
    """
    cells, half_width = index(cells), Fraction(half_width)
    if cells < 1 or half_width <= 0:
        raise ValueError("cell count and half width must be positive")
    width = arb(2 * half_width.numerator) / (half_width.denominator * cells)
    terms = _arb_prime_terms(powers)
    zero = arb(0)
    lags = []
    for lag in range(cells):
        separation = lag * width
        entry = arb(0)
        for amplitude, shift in terms:
            overlap_right = (width - abs(separation - shift)).max(zero)
            overlap_left = (width - abs(separation + shift)).max(zero)
            entry += amplitude * (overlap_right + overlap_left) / (2 * width)
        lags.append(entry)
    return lags


def prime_shift_matvec_arb(vector, half_width=P236_HALF_WIDTH, powers=P236_PRIME_POWERS):
    """Enclose the exact finite action on the binary64 input coordinates."""
    vector = _real_binary64_vector(vector)
    # arb(float) preserves the original dyadic value at precision >= 53.
    with ctx.workprec(max(ctx.prec, 53)):
        lags = prime_shift_lags_arb(vector.size, half_width, powers)
        coefficients = [arb(float(coefficient)) for coefficient in vector]
        return [
            sum((lags[abs(row - column)] * value for column, value in enumerate(coefficients)), arb(0))
            for row in range(vector.size)
        ]


def audit_fft_vector(vector, *, half_width=P236_HALF_WIDTH, powers=P236_PRIME_POWERS, precision_bits=128):
    """Enclose the error of one actual binary64 FFT action against the Arb replay."""
    precision_bits = index(precision_bits)
    if precision_bits < 53:
        raise ValueError("at least 53 Arb bits are needed for exact binary64 inputs")
    vector = _real_binary64_vector(vector)
    powers = tuple((index(power), index(prime)) for power, prime in powers)
    approximate = PrimeShiftOperator(vector.size, half_width=half_width, powers=powers).matvec(vector)
    if not np.isfinite(approximate).all():
        raise ArithmeticError("nonfinite FFT output")
    with ctx.workprec(precision_bits):
        exact = prime_shift_matvec_arb(vector, half_width, powers)
        errors = [(arb(float(value)) - enclosure).abs_upper() for value, enclosure in zip(approximate, exact, strict=True)]
        return {
            "cells": vector.size,
            "half_width": str(Fraction(half_width)),
            "prime_powers": powers,
            "precision_bits": precision_bits,
            "input_vector_dyadic": [_dyadic(arb(float(value))) for value in vector],
            "max_abs_error_upper_dyadic": _dyadic(max(errors)),
            "coordinates": [
                {
                    "fft_output_dyadic": _dyadic(arb(float(approximation))),
                    "exact_lower_dyadic": _dyadic(enclosure.lower()),
                    "exact_upper_dyadic": _dyadic(enclosure.upper()),
                    "abs_error_upper_dyadic": _dyadic(error),
                }
                for approximation, enclosure, error in zip(approximate, exact, errors, strict=True)
            ],
        }


def _cosine_taylor(cell_index, subdivisions, terms, radius, remainder):
    center = arb(2 * cell_index + 1) / (2 * subdivisions)
    symbol, derivative = arb(0), arb(0)
    for amplitude, frequency in terms:
        sine, cosine = (center * frequency).sin_cos()
        symbol += amplitude * cosine
        derivative -= amplitude * frequency * sine
    return symbol, radius * abs(derivative) + remainder


def certify_cancellation_shells(*, powers=P236_PRIME_POWERS, endpoint=504, subdivisions_per_unit=32, precision_bits=128):
    """Bound the cosine symbol on dyadic shells using signed Taylor coefficients.

    At a subcell center c and radius r, the bound is
    |S(c)| + r |S'(c)| + r²/2 sum_q a_q log(q)². The first two terms are
    accumulated with their signs before taking magnitudes. The final remainder
    follows from the global bound on the second derivative of the finite sum.
    """
    endpoint, subdivisions_per_unit = index(endpoint), index(subdivisions_per_unit)
    precision_bits = index(precision_bits)
    if endpoint < 1:
        raise ValueError("frequency endpoint must be a positive integer")
    if subdivisions_per_unit < 1:
        raise ValueError("subdivisions per unit must be a positive integer")
    if precision_bits < 53:
        raise ValueError("use at least 53 Arb bits")
    powers = tuple((index(power), index(prime)) for power, prime in powers)
    with ctx.workprec(precision_bits):
        terms = _arb_prime_terms(powers)
        triangle = sum((amplitude for amplitude, _ in terms), arb(0))
        if not triangle > 0:
            raise ValueError("at least one positive prime amplitude is required")
        curvature = sum((amplitude * frequency ** 2 for amplitude, frequency in terms), arb(0))
        radius = arb(1) / (2 * subdivisions_per_unit)
        remainder = radius ** 2 * curvature / 2
        shells = []
        lower, upper = 0, 1
        while lower < endpoint:
            first = lower * subdivisions_per_unit
            stop = upper * subdivisions_per_unit
            taylor_upper = arb(0)
            for cell_index in range(first, stop):
                symbol, variation = _cosine_taylor(cell_index, subdivisions_per_unit, terms, radius, remainder)
                bound = (abs(symbol) + variation).upper()
                taylor_upper = max(taylor_upper, bound)
            envelope_upper = min(taylor_upper, triangle.upper())
            ratio_upper = min((envelope_upper / triangle).upper(), arb(1))
            shells.append({
                "range": [lower, upper],
                "subcells": stop - first,
                "envelope_upper_dyadic": _dyadic(envelope_upper),
                "ratio_upper_dyadic": _dyadic(ratio_upper),
                "strict_reduction": bool(ratio_upper < 1),
            })
            lower, upper = upper, min(2 * upper, endpoint)
        return {
            "precision_bits": precision_bits,
            "prime_powers": powers,
            "subdivisions_per_unit": subdivisions_per_unit,
            "triangle_lower_dyadic": _dyadic(triangle.lower()),
            "triangle_upper_dyadic": _dyadic(triangle.upper()),
            "shells": shells,
        }


def certify_p236_tail_start(*, subdivisions_per_unit=128, precision_bits=128):
    """Certify a finite signed prime estimate for the unchanged P236 beta.

    B(t) = log(t/(2*pi)) - 1/t increases for t > 0. On each subcell,
    B(left) minus beta minus a signed Taylor upper bound on S proves
    a lower bound for B(t) - beta - S(t) across that subcell.
    """
    subdivisions_per_unit = index(subdivisions_per_unit)
    precision_bits = index(precision_bits)
    if subdivisions_per_unit < 1 or precision_bits < 53:
        raise ValueError("use positive subdivisions and at least 53 Arb bits")
    endpoint = 504
    with ctx.workprec(precision_bits):
        terms = _arb_prime_terms(P236_PRIME_POWERS)
        triangle = sum((amplitude for amplitude, _ in terms), arb(0))
        curvature = sum((amplitude * frequency ** 2 for amplitude, frequency in terms), arb(0))
        radius = arb(1) / (2 * subdivisions_per_unit)
        remainder = radius ** 2 * curvature / 2
        two_pi = 2 * arb.pi()
        beta = (arb(endpoint) / two_pi).log() - arb(1) / endpoint - triangle
        start = endpoint
        certified_lower = None
        first_inconclusive_cell = None
        for left_integer in range(endpoint - 1, 0, -1):
            unit_bounds = []
            first = left_integer * subdivisions_per_unit
            for cell_index in range(first + subdivisions_per_unit - 1, first - 1, -1):
                left = arb(cell_index) / subdivisions_per_unit
                symbol, variation = _cosine_taylor(cell_index, subdivisions_per_unit, terms, radius, remainder)
                bound = ((left / two_pi).log() - 1 / left - beta - symbol - variation).lower()
                if bound <= 0:
                    first_inconclusive_cell = {
                        "interval_numerators": [cell_index, cell_index + 1],
                        "lower_bound_dyadic": _dyadic(bound),
                    }
                    break
                unit_bounds.append(bound)
            if first_inconclusive_cell is not None:
                break
            start = left_integer
            unit_lower = min(unit_bounds)
            certified_lower = unit_lower if certified_lower is None else min(certified_lower, unit_lower)
        return {
            "precision_bits": precision_bits,
            "subdivisions_per_unit": subdivisions_per_unit,
            "subcell_width": str(Fraction(1, subdivisions_per_unit)),
            "prime_powers": P236_PRIME_POWERS,
            "original_cutoff": endpoint,
            "certified_range": [start, endpoint],
            "certified_subcells": (endpoint - start) * subdivisions_per_unit,
            "beta_lower_dyadic": _dyadic(beta.lower()),
            "beta_upper_dyadic": _dyadic(beta.upper()),
            "strict_gap_lower_dyadic": _dyadic(certified_lower) if certified_lower is not None else None,
            "first_inconclusive_cell": first_inconclusive_cell,
        }


def run_audit():
    """Generate the reproducible default ES8-P249 audit with no wall-clock fields."""
    vector = [((3 * index * index + index) % 23 - 11) / 8 for index in range(65)]
    return {
        "python_flint_version": flint.__version__,
        "numpy_version": np.__version__,
        "fft_vector": audit_fft_vector(vector),
        "cancellation": certify_cancellation_shells(),
        "p236_tail_start": certify_p236_tail_start(),
        "p236_tail_coarser": [certify_p236_tail_start(subdivisions_per_unit=n) for n in (32, 64)],
    }


if __name__ == "__main__":
    print(json.dumps(run_audit(), indent=2))
