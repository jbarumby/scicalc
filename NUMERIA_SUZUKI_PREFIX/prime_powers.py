"""Enumerate prime powers in increasing order with a segmented sieve."""

from heapq import merge
from math import isqrt
from operator import index


def _clear_multiples(flags, low, high, prime):
    first = max(prime * prime, ((low + prime - 1) // prime) * prime)
    if first <= high:
        flags[first - low :: prime] = bytes((high - first) // prime + 1)


def _base_primes(limit):
    flags = bytearray(b"\1") * (limit - 1)
    for candidate in range(2, isqrt(limit) + 1):
        if flags[candidate - 2]:
            _clear_multiples(flags, 2, limit, candidate)
    return [value for value, flag in enumerate(flags, 2) if flag]


def _prime_segments(limit, block_size, base_primes):
    for low in range(2, limit + 1, block_size):
        high = min(low + block_size - 1, limit)
        flags = bytearray(b"\1") * (high - low + 1)
        for prime in base_primes:
            if prime * prime > high:
                break
            _clear_multiples(flags, low, high, prime)
        offset = flags.find(b"\1")
        while offset != -1:
            prime = low + offset
            yield prime, prime
            offset = flags.find(b"\1", offset + 1)


def iter_prime_powers(limit, block_size=2**20):
    """Return a stream of exact pairs (q, p) for q = p**k <= limit, k >= 1.

    Composite integers are marked using primes at most isqrt(limit). Only
    one segment, the base primes, and the powers with k >= 2 are retained.
    """
    limit = index(limit)
    block_size = index(block_size)
    if limit < 2:
        raise ValueError("the sieve limit must be at least 2")
    if block_size < 1:
        raise ValueError("the sieve block must contain at least one integer")

    base_primes = _base_primes(isqrt(limit))
    higher_powers = []
    for prime in base_primes:
        power = prime * prime
        while power <= limit:
            higher_powers.append((power, prime))
            power *= prime
    higher_powers.sort()
    return merge(_prime_segments(limit, block_size, base_primes), higher_powers)
