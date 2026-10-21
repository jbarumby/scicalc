
"""Exact decimal validation and quadrature mesh counting."""

from decimal import Decimal, InvalidOperation, ROUND_CEILING
import math


def validate_decimal(value):
    """Validate a finite, positive decimal input."""
    if isinstance(value, bool):
        raise ValueError("Boolean values are not valid")

    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError("Value must be finite")

    try:
        number = Decimal(str(value))
    except (InvalidOperation, ValueError, TypeError) as exc:
        raise ValueError("Invalid decimal value") from exc

    if not number.is_finite() or number <= 0:
        raise ValueError("Value must be finite and positive")

    return str(number)


def exact_cell_count(R_value, step_value):
    """Compute ceil(R / step) using exact decimal arithmetic."""
    R = Decimal(validate_decimal(R_value))
    step = Decimal(validate_decimal(step_value))

    # Decimal integer arithmetic avoids rounded division.
    numerator, denominator = R.as_integer_ratio()
    step_num, step_den = step.as_integer_ratio()

    a = numerator * step_den
    b = denominator * step_num

    return (a + b - 1) // b
