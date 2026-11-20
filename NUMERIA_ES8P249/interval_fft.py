"""Outward-rounded circulant convolution using FLINT complex balls."""

from __future__ import annotations

import operator
from collections.abc import Iterable

from flint import acb, arb, ctx


def _fft_in_place(values: list[acb], roots: tuple[acb, ...]) -> None:
    """Apply a radix-2 DFT with the supplied powers of its root of unity."""
    size = len(values)
    reversed_index = 0
    for index in range(1, size):
        bit = size // 2
        while reversed_index & bit:
            reversed_index ^= bit
            bit //= 2
        reversed_index ^= bit
        if index < reversed_index:
            values[index], values[reversed_index] = values[reversed_index], values[index]

    span = 2
    while span <= size:
        half_span = span // 2
        root_step = size // span
        for start in range(0, size, span):
            for offset in range(half_span):
                left = start + offset
                right = left + half_span
                first = values[left]
                second = roots[offset * root_step] * values[right]
                values[left] = first + second
                values[right] = first - second
        span *= 2


def toeplitz_matvec_arb(
    lags: Iterable, vector: Iterable, *, precision_bits: int = 128
) -> tuple[acb, ...]:
    """Enclose the action of the symmetric matrix A[i,j] = lags[abs(i-j)].

    Inputs can be real or complex values accepted by ``acb``, including balls.
    Negative lags use the same coefficients as positive lags. The two forward
    transforms and one inverse transform use an embedding of length at least
    ``2 * len(lags) - 1`` and require no matrix entries beyond the lags.
    """
    precision_bits = operator.index(precision_bits)
    if precision_bits < 2:
        raise ValueError("precision_bits must be at least 2")

    with ctx.workprec(precision_bits):
        coefficients = tuple(acb(value) for value in lags)
        inputs = tuple(acb(value) for value in vector)
        dimension = len(coefficients)
        if dimension == 0 or len(inputs) != dimension:
            raise ValueError("lags and vector must be nonempty and have equal lengths")
        if any(not value.is_finite() for value in (*coefficients, *inputs)):
            raise ValueError("lags and vector must be finite")

        size = 1 << (2 * dimension - 2).bit_length()
        zero = acb(0)
        embedding = [*coefficients, *([zero] * (size - 2 * dimension + 1)), *reversed(coefficients[1:])]
        padded_vector = [*inputs, *([zero] * (size - dimension))]

        angle_step = -2 * arb.pi() / size
        forward_roots = tuple(acb(0, k * angle_step).exp() for k in range(size // 2))
        inverse_roots = tuple(root.conjugate() for root in forward_roots)
        _fft_in_place(embedding, forward_roots)
        _fft_in_place(padded_vector, forward_roots)

        product = [a * b for a, b in zip(embedding, padded_vector, strict=True)]
        _fft_in_place(product, inverse_roots)
        return tuple(value / size for value in product[:dimension])
