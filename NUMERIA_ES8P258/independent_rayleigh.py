"""Exact Rayleigh replay for dyadic NUMERIA matrix enclosures."""

from fractions import Fraction
from operator import index

from numeria.p237.p237_core import as_endpoint_interval


def _decode(interval):
    lower, upper, exponent = as_endpoint_interval(interval)
    scale = Fraction(2) ** exponent
    return lower * scale, upper * scale


def quadratic_interval(upper_triangle, witness):
    """Enclose the form on a nonzero sparse exact witness."""
    witness = tuple(Fraction(value) for value in witness)
    nonzero = {i: value for i, value in enumerate(witness) if value}
    if not nonzero:
        raise ValueError("a Rayleigh witness must be nonzero")
    lower = upper = Fraction(0)
    previous = (-1, -1)
    handled = 0
    for i, j, endpoints in upper_triangle:
        i, j = index(i), index(j)
        if not 0 <= i <= j < len(witness) or (i, j) <= previous:
            raise ValueError("invalid or unordered upper triangle")
        previous = i, j
        if i not in nonzero or j not in nonzero:
            continue
        lo, hi = _decode(endpoints)
        coefficient = nonzero[i] * nonzero[j]
        if i != j:
            coefficient *= 2
        lower += coefficient * (lo if coefficient > 0 else hi)
        upper += coefficient * (hi if coefficient > 0 else lo)
        handled += 1
    if handled != len(nonzero) * (len(nonzero) + 1) // 2:
        raise ValueError("matrix enclosure does not cover the witness support")
    return lower, upper


def replay_rayleigh(witness, components, *, beta=None, operator_errors=None):
    """Replay the Rayleigh enclosure using exact fractions and error budgets."""
    witness = tuple(Fraction(value) for value in witness)
    norm_squared = sum((value * value for value in witness), Fraction(0))
    if not norm_squared or not components:
        raise ValueError("expected a nonzero witness and at least one matrix box")
    lower = upper = Fraction(0)
    component_results = {}
    for name, (upper_triangle, multiplier) in components.items():
        multiplier = Fraction(multiplier)
        lo, hi = quadratic_interval(upper_triangle, witness)
        scaled_lo, scaled_hi = sorted((multiplier * lo, multiplier * hi))
        lower += scaled_lo
        upper += scaled_hi
        component_results[name] = {
            "scalar_multiplier_exact": str(multiplier),
            "quadratic_interval_exact": [str(scaled_lo), str(scaled_hi)],
        }
    beta_results = None
    if beta is not None:
        beta_lo, beta_hi = _decode(beta)
        lower += beta_lo * norm_squared
        upper += beta_hi * norm_squared
        beta_results = [str(beta_lo), str(beta_hi)]
    errors = {name: Fraction(value) for name, value in (operator_errors or {}).items()}
    if any(value < 0 for value in errors.values()):
        raise ValueError("operator norm bounds must be nonnegative")
    total_error = sum(errors.values(), Fraction(0))
    lower -= total_error * norm_squared
    upper += total_error * norm_squared
    sign = "negative" if upper < 0 else "positive" if lower > 0 else "inconclusive"
    return {
        "schema": "NUMERIA.ES8P258.EXACT_RATIONAL_RAYLEIGH.1",
        "dimension": len(witness),
        "nonzero_witness_exact": [[i, str(value)] for i, value in enumerate(witness) if value],
        "norm_squared_exact": str(norm_squared),
        "matrix_components": component_results,
        "beta_interval_exact": beta_results,
        "operator_norm_errors_exact": {name: str(value) for name, value in errors.items()},
        "total_operator_norm_error_exact": str(total_error),
        "unnormalized_interval_exact": [str(lower), str(upper)],
        "rayleigh_interval_exact": [str(lower / norm_squared), str(upper / norm_squared)],
        "witness_strict_sign": sign,
    }
