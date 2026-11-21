"""Isolate zeta zeros in Arb and audit the ordinates by Hardy signs and a count."""

from hashlib import sha256
import json

from flint import acb, arb, ctx

from NUMERIA_ES8P259.suzuki import _dyadic


def hardy_z(t):
    """Enclose real Hardy Z(t) using the gamma factor and the zeta function."""
    half_t = t / 2
    theta = acb(arb(1) / 4, half_t).lgamma().imag - half_t * arb.pi().log()
    zeta = acb(arb(1) / 2, t).zeta()
    return (acb(0, theta).exp() * zeta).real


def _count_height(last, count):
    height = arb(last.upper().floor() + 1)
    if height.is_exact() and last < height and height.zeta_nzeros() == count:
        return height

    next_lower = acb.zeta_zero(count + 1).imag.lower()
    height = ((last.upper() + next_lower) / 2).mid()
    if not (last < height < next_lower) or not height.zeta_nzeros() == count:
        raise ArithmeticError(f"cannot certify exactly {count} zeros below the cutoff")
    return height


def _audit_signs(ordinates):
    minimum = None
    with ctx.workprec(2 * ctx.prec):
        for index, interval in enumerate(ordinates, 1):
            lower, upper = hardy_z(interval.lower()), hardy_z(interval.upper())
            if not (lower > 0 and upper < 0 or lower < 0 and upper > 0):
                raise ArithmeticError(f"Hardy Z endpoint signs are undecided at zero {index}")
            margin = min(abs(lower).lower(), abs(upper).lower())
            if minimum is None or margin < minimum:
                minimum = margin
    return minimum


def _bounds(interval):
    return [_dyadic(interval.lower()), _dyadic(interval.upper())]


def verified_ordinates(count):
    """Return disjoint ordinate balls and an independently checked zero count.

    Opposite Hardy signs prove existence in every padded ball. The total zero
    count through a separating height then certifies simplicity and that no
    preceding zero, including an off-line zero, is missing from the list.
    """
    if type(count) is not int or count < 1:
        raise ValueError("the number of zeros must be a positive integer")
    half = arb(1) / 2
    padding_exponent = 16 - ctx.prec
    padding = arb(0, arb(2)**padding_exponent)
    ordinates = []
    for index, zero in enumerate(acb.zeta_zeros(1, count), 1):
        ordinate = zero.imag + padding
        if not zero.real == half or not ordinate.is_finite() or not ordinate > 0:
            raise ArithmeticError(f"invalid critical-line root enclosure at zero {index}")
        if ordinates and not ordinates[-1] < ordinate:
            raise ArithmeticError(f"overlapping or unordered enclosures at zero {index}")
        ordinates.append(ordinate)
    if len(ordinates) != count:
        raise ArithmeticError(f"expected {count} enclosures, received {len(ordinates)}")

    height = _count_height(ordinates[-1], count)
    sign_margin = _audit_signs(ordinates)
    digest = sha256()
    for index, ordinate in enumerate(ordinates, 1):
        payload = json.dumps([index, _bounds(ordinate)], separators=(",", ":"))
        digest.update((payload + "\n").encode("ascii"))
    integer_height = height.unique_fmpz()
    return ordinates, {
        "zero_count": count,
        "count_at_height": count,
        "count_height": int(integer_height) if integer_height is not None else None,
        "count_height_dyadic": _dyadic(height),
        "first_ordinate_bounds_dyadic": _bounds(ordinates[0]),
        "last_ordinate_bounds_dyadic": _bounds(ordinates[-1]),
        "ordered_ordinate_bounds_sha256": digest.hexdigest(),
        "hardy_endpoint_sign_verified_count": count,
        "hardy_endpoint_sign_minimum_dyadic": _dyadic(sign_margin),
        "outward_padding_dyadic": ["1", padding_exponent],
    }
