"""Exact constants for the separated-strip bound with the raised floor."""

from fractions import Fraction

from NUMERIA_ES8P257.p257_strip_taper import (
    HIGH_REMAINDER, KAPPA_UP, PI_LO_TARGET, arctanh_lower, prime_7_constants,
)


OMEGA = 2520
LOSS = 1416


def certify_strip_constants():
    """Verify the scalar Schur bound for the two-strip estimate."""
    prime_7_constants()
    high_lower = arctanh_lower(Fraction(2, 3), 8) + Fraction(1, 630) - Fraction(1, 4)
    used_high = Fraction(51, 40)
    if not high_lower > used_high + HIGH_REMAINDER:
        raise ArithmeticError("the raised-floor model lacks the required high reserve")
    loss = (OMEGA + 1) * KAPPA_UP * used_high / (PI_LO_TARGET * (used_high - KAPPA_UP))
    if not loss < LOSS:
        raise ArithmeticError("the separated-strip estimate exceeds the claimed loss")
    return {
        "scope": "conditional on an auxiliary coercivity margin for the raised floor",
        "four_prime_beta_exact": "1/4",
        "absorption_frequency": OMEGA,
        "high_reserve_strict_lower_exact": str(high_lower),
        "retained_high_reserve_exact": str(HIGH_REMAINDER),
        "schur_high_consumption_exact": str(used_high),
        "strip_loss_factor_exact_upper": str(loss),
        "strip_loss_factor_integer_upper": LOSS,
    }
