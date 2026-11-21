"""Compare every prime-power pair with an independent, unsegmented sieve."""

import argparse
from heapq import merge
from itertools import chain, compress
from math import isqrt

from .prime_powers import iter_prime_powers


DEFAULT_LIMIT = 10**8


def _monolithic_prime_powers(limit):
    """Use an odd-only sieve for primes and trial division for higher powers."""
    flags = bytearray(b"\1") * ((limit + 1) // 2)
    flags[0] = 0
    root = isqrt(limit)
    for candidate in range(3, root + 1, 2):
        if flags[candidate // 2]:
            first = candidate**2 // 2
            flags[first::candidate] = bytes((len(flags) - 1 - first) // candidate + 1)

    higher = []
    for prime in range(2, root + 1):
        if any(prime % divisor == 0 for divisor in range(2, isqrt(prime) + 1)):
            continue
        exponent = 2
        while (power := prime**exponent) <= limit:
            higher.append((power, prime))
            exponent += 1
    primes = chain((2,), compress(range(1, limit + 1, 2), flags))
    return merge(((prime, prime) for prime in primes), sorted(higher))


def audit(limit):
    segmented = iter_prime_powers(limit)
    monolithic = _monolithic_prime_powers(limit)
    count = prime_count = 0
    for count, (actual, expected) in enumerate(zip(segmented, monolithic, strict=True), 1):
        if actual != expected:
            raise ArithmeticError(f"pair {count} differs ({actual} versus {expected})")
        q, prime = actual
        prime_count += q == prime
    if limit == DEFAULT_LIMIT and (count, prime_count) != (5_762_859, 5_761_455):
        raise ArithmeticError(f"unexpected prime-power and prime counts ({count}, {prime_count})")
    return count, prime_count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=DEFAULT_LIMIT)
    args = parser.parse_args()
    count, prime_count = audit(args.limit)
    print(f"{count} matching pairs through {args.limit}, including {prime_count} primes")


if __name__ == "__main__":
    main()
