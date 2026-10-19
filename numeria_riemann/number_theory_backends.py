"""Optional independent checks for prime and Möbius computations.

These helpers compare outputs from separate libraries; they do not establish
any conjecture. SymPy is optional, and primesieve is optional.
"""
from __future__ import annotations

import hashlib
import json


def sympy_von_mangoldt(n: int) -> float:
    """Independent von Mangoldt evaluation using SymPy factorization."""
    if not isinstance(n, int) or isinstance(n, bool) or n < 1:
        raise ValueError("n must be a positive integer")
    try:
        from sympy import factorint
    except ImportError as exc:
        raise ImportError("install optional dependency: pip install sympy") from exc
    factors = factorint(n)
    if n == 1 or len(factors) != 1:
        return 0.0
    prime, exponent = next(iter(factors.items()))
    if prime ** exponent != n:
        return 0.0
    from math import log
    return log(int(prime))


def sympy_primes(limit: int) -> list[int]:
    """Return primes <= limit using SymPy's prime-range implementation."""
    if not isinstance(limit, int) or isinstance(limit, bool) or limit < 0:
        raise ValueError("limit must be a non-negative integer")
    try:
        from sympy import primerange
    except ImportError as exc:
        raise ImportError("install optional dependency: pip install sympy") from exc
    return list(primerange(2, limit + 1))


def primesieve_primes(limit: int) -> list[int]:
    """Return primes <= limit using the optional primesieve Python binding."""
    if not isinstance(limit, int) or isinstance(limit, bool) or limit < 0:
        raise ValueError("limit must be a non-negative integer")
    try:
        import primesieve
    except ImportError as exc:
        raise ImportError(
            "primesieve is optional; install its Python binding separately"
        ) from exc
    return list(primesieve.primes(limit))


def compare_prime_lists(limit: int, *, include_primesieve: bool = True) -> dict:
    """Compare prime lists and report counts, equality, and SHA-256 hashes."""
    reference = sympy_primes(limit)
    digest = lambda values: hashlib.sha256(
        ",".join(map(str, values)).encode("ascii")
    ).hexdigest()
    result = {
        "limit_inclusive": limit,
        "sympy_count": len(reference),
        "sympy_sha256": digest(reference),
        "status": "crosscheck_only_not_proof",
        "comparisons": {},
    }
    if include_primesieve:
        try:
            candidate = primesieve_primes(limit)
            result["comparisons"]["primesieve"] = {
                "count": len(candidate),
                "sha256": digest(candidate),
                "matches_sympy": candidate == reference,
            }
        except ImportError:
            result["comparisons"]["primesieve"] = {
                "available": False,
                "status": "optional_dependency_not_installed",
            }
    return result


def compare_mobius_values(limit: int) -> dict:
    """Cross-check Möbius values 1..limit using SymPy and return a hash."""
    if not isinstance(limit, int) or isinstance(limit, bool) or limit < 1:
        raise ValueError("limit must be a positive integer")
    try:
        from sympy import mobius
    except ImportError as exc:
        raise ImportError("install optional dependency: pip install sympy") from exc
    values = [int(mobius(n)) for n in range(1, limit + 1)]
    encoded = ",".join(map(str, values)).encode("ascii")
    return {
        "limit_inclusive": limit,
        "count": len(values),
        "sha256": hashlib.sha256(encoded).hexdigest(),
        "sum": sum(values),
        "status": "sympy_reference_values_not_conjecture_evidence",
    }
