"""Directed Binet tail bounds and exact Schur reserve for the V.39 symbol."""

import argparse
from collections import Counter
import csv
from fractions import Fraction
from hashlib import sha256
import json
from math import ceil, lcm
from pathlib import Path

from flint import acb, arb, ctx


BINET_SOURCE = "https://dlmf.nist.gov/5.9.E15"
DIGAMMA_SERIES_SOURCE = "https://dlmf.nist.gov/5.7.E6"
PRIME_POWERS = ((2, 2), (3, 3), (4, 2), (5, 5), (7, 7), (8, 2), (9, 3), (11, 11), (13, 13))
DEFAULT_ARCHIVE = Path(__file__).resolve().parents[1] / "NUMERIA_P93V39"


def _fraction(value):
    if type(value) not in (int, str, Fraction):
        raise ValueError("expected an exact rational, not a floating-point value")
    try:
        return Fraction(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError("invalid exact rational") from exc


def _ball(value):
    value = Fraction(value)
    return arb(value.numerator) / value.denominator


def _dyadic(point):
    mantissa, exponent = point.man_exp()
    return [str(mantissa), int(exponent)]


def _dyadic_fraction(point):
    mantissa, exponent = point.man_exp()
    return int(mantissa) * Fraction(2)**int(exponent)


def prime_amplitude():
    """V.39's nine terms, with their historical factor of two."""
    return sum((2 * arb(p).log() / arb(q).sqrt() for q, p in PRIME_POWERS), arb(0))


def binet_majorants(cutoff):
    """Return archimedean anchor, signed exponential loss, and compact upper error.

    Let a=1/4, y=|t|/2, z=a+iy, and v=Re(u^2+z^2). The real
    reciprocal is v/(v^2+(2ay)^2). Its absolute value is at most
    1/(4ay). For 0<=u<=y/2 it is nonpositive, and its negative
    is at most 1/(3*y^2/4-a^2). Use the total thermal integral
    1/24 on this compact interval. Above y/2 integrate the
    exponential majorant of u/(exp(2*pi*u)-1).
    """
    cutoff = _fraction(cutoff)
    if cutoff <= 0 or 3 * cutoff**2 <= 1:
        raise ValueError("the Binet sign split requires a positive cutoff with 3*cutoff^2 > 1")

    a = arb(1) / 4
    y = _ball(cutoff / 2)
    pi = arb.pi()
    radius_squared = a**2 + y**2
    elementary = radius_squared.log() / 2 - a / (2 * radius_squared) - pi.log()
    exponential = (-pi * y).exp()
    loss = exponential / (1 - exponential) * (1 + 1 / (pi * y)) / (8 * pi * a)
    compact_upper = 1 / (12 * (3 * y**2 / 4 - a**2))
    return elementary, loss, compact_upper


def audit(cutoff=4500, beta=2, precision=192, no_deficit_cutoff=10**6):
    """Certify the high-frequency upper loss and the eventual zero deficit."""
    if type(precision) is not int or precision < 96:
        raise ValueError("precision must be an integer of at least 96 bits")
    cutoff, beta, no_deficit_cutoff = map(_fraction, (cutoff, beta, no_deficit_cutoff))
    if no_deficit_cutoff < cutoff:
        raise ValueError("the eventual-zero threshold must be above the loss cutoff")
    with ctx.workprec(precision):
        elementary, loss, compact_upper = binet_majorants(cutoff)
        floor = (elementary - loss.upper()).lower()
        pi = arb.pi()
        mass = prime_amplitude()
        if not mass < _ball(Fraction(9944, 1000)):
            raise ArithmeticError("the nine-term historical mass bound did not close")
        # log|1/4+it/2| <= log(|t|/2) + 1/(8*t^2), and the pole term
        # -Re(1/(2z)) is nonpositive. All error terms decrease with |t|.
        unit_error = 1 / (8 * _ball(cutoff)**2)
        reference = mass - _ball(beta) - (2 * pi).log() + unit_error
        original_remainder_upper = 1 / (3 * _ball(cutoff))
        original_constant = (reference + original_remainder_upper).upper()
        refined_constant = (reference + compact_upper.upper() + loss.upper()).upper()
        scale = 10**4
        refined_grid = Fraction(ceil(_dyadic_fraction(refined_constant) * scale), scale)

        # Independently evaluate the digamma at the cutoff. Its real part
        # is increasing in |t| by the summands of the Euler series.
        probe = acb(arb(1) / 4, _ball(cutoff / 2)).digamma().real - pi.log()
        if not probe.is_finite() or not floor.is_finite() or not probe > floor:
            raise ArithmeticError("the independent Arb anchor did not exceed the Binet floor")

        # The elementary anchor increases. The exponential loss decreases
        # because each factor q/(1-q), 1+1/(pi*y) decreases with y.
        high_anchor, high_error, _ = binet_majorants(no_deficit_cutoff)
        high_gap = (high_anchor - high_error.upper() - mass - _ball(beta)).lower()
        if not high_gap > 0:
            raise ArithmeticError("the requested eventual-zero deficit threshold did not close")
        if not original_constant < 8 or not _ball(cutoff).log() + 8 > 0:
            raise ArithmeticError("the historical upper-loss function did not close")
        if not _ball(cutoff).log() + _ball(refined_grid) > 0:
            raise ArithmeticError("the refined upper-loss function must be nonnegative")

        return {
            "schema": "NUMERIA.P93V40.BINET.1",
            "archimedean_symbol": "Re psi(1/4+i*t/2)-log(pi)",
            "historical_symbol": "archimedean_symbol - 2*sum_{q} Lambda(q)*cos(t*log(q))/sqrt(q)",
            "prime_powers": [list(pair) for pair in PRIME_POWERS],
            "cutoff_exact": str(cutoff),
            "beta_exact": str(beta),
            "precision_bits": precision,
            "prime_amplitude_upper_dyadic": _dyadic(mass.upper()),
            "uniform_archimedean_binet_floor_lower_dyadic": _dyadic(floor),
            "uniform_archimedean_digamma_floor_lower_dyadic": _dyadic(probe.lower()),
            "binet_exponential_loss_upper_dyadic": _dyadic(loss.upper()),
            "historical_upper_loss_constant": "8",
            "historical_budget_derived_constant_upper_dyadic": _dyadic(original_constant),
            "refined_upper_loss_constant_exact": str(refined_grid),
            "refined_budget_derived_constant_upper_dyadic": _dyadic(refined_constant),
            "no_deficit_from_exact": str(no_deficit_cutoff),
            "no_deficit_above_beta_lower_dyadic": _dyadic(high_gap),
            "references": [BINET_SOURCE, DIGAMMA_SERIES_SOURCE],
        }


def _components(csv_path, cutoff):
    components = []
    previous = -1
    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        if next(reader, None) != ["left16", "right16", "deficit_micro"]:
            raise ValueError("unexpected header in deficit table")
        for start, end, deficit in (tuple(map(int, row)) for row in reader):
            if not (start >= previous and end == start + 1 and deficit > 0):
                raise ValueError("unsorted or invalid deficit interval")
            previous = end
            if start < 16 * cutoff:
                continue
            if components and components[-1][1] == start:
                l, _, mass, height = components[-1]
                components[-1] = (l, end, mass + deficit, max(height, deficit))
            else:
                components.append((start, end, deficit, deficit))
    return components


def replay_reserve(archive=DEFAULT_ARCHIVE):
    """Rebuild the two adaptive Schur bounds from inherited CSV enclosures.

    The rounded sinc bound uses, on each pair of bands, a<=1387/1000
    and pi>333/106. The latter, the kernel grid, and the archived
    multipliers and weights are kept as stated in the V.39 proof.
    """
    archive = Path(archive)
    settings_bytes = (archive / "hybrid39.json").read_bytes()
    settings = json.loads(settings_bytes)
    csv_path = archive / "deficit_2000000.csv"
    components = _components(csv_path, 4500)
    weights, factors = settings["weights"], settings["ks"]
    expected = [2 if mass < 16 * 10**6 else 14 for _, _, mass, _ in components]
    if factors != expected or any(type(k) is not int for k in factors) or len(weights) != len(components):
        raise ValueError("component mass or multiplier mismatch")
    if any(type(w) is not int or w <= 0 for w in weights):
        raise ValueError("Schur weights must be positive integers")
    with ctx.workprec(96):
        if not arb.pi() > _ball(Fraction(333, 106)):
            raise ArithmeticError("the Schur lower bound on pi did not close")

    grid = 2**20
    den = lcm(*(k - 1 for k in factors))
    cap = (1387 * grid + 999) // 1000
    weighted_components = [
        (l, r, mass * w * k * (den // (k - 1)))
        for (l, r, mass, _), w, k in zip(components, weights, factors, strict=True)
    ]
    greatest_num, greatest_den = 0, 1
    for (l, r, _, _), w in zip(components, weights, strict=True):
        row_sum = 0
        for start, end, weighted_mass in weighted_components:
            gap = max(0, start - r, l - end)
            difference_kernel = min(cap, (16 * grid + gap - 1) // gap) if gap else cap
            sum_kernel = (16 * grid + l + start - 1) // (l + start)
            row_sum += (difference_kernel + sum_kernel) * weighted_mass
        if row_sum * greatest_den > greatest_num * w:
            greatest_num, greatest_den = row_sum, w
    c = Fraction(106 * greatest_num, 333 * den * grid * 16 * 10**6 * greatest_den)
    d = max(Fraction(height * k, 10**6 * (k - 1))
            for (_, _, _, height), k in zip(components, factors, strict=True))
    lam = Fraction(settings["lambda"])
    if not 0 < lam < 1:
        raise ValueError("the harmonic parameter must lie in (0,1)")
    s = lam * d / ((1 - lam) * c)
    f_1 = c * d / ((1 - lam) * c + lam * d)
    if any(calculated != Fraction(settings[key]) for key, calculated in
           {"c": c, "d": d, "s": s, "f1": f_1}.items()):
        raise ArithmeticError("the reconstructed Schur bounds differ from V.39")
    if not (0 < c < 2 and 0 < f_1 < 2):
        raise ArithmeticError("the archived penalties have exhausted the reference")
    return {
        "basis": "V.39 deficit CSV and exact weighted row sums",
        "deficit_csv_sha256": sha256(csv_path.read_bytes()).hexdigest(),
        "configuration_sha256": sha256(settings_bytes).hexdigest(),
        "components": len(components),
        "multipliers": {str(k): n for k, n in sorted(Counter(factors).items())},
        "C_exact": str(c),
        "D_exact": str(d),
        "lambda_exact": str(lam),
        "s_exact": str(s),
        "F_1_exact": str(f_1),
        "adaptive_spectral_reserve_exact": str(2 - c),
        "harmonic_spectral_reserve_exact": str(2 - f_1),
        "maximum_halfwidth_exact": "1387/1000",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--precision", type=int, default=192)
    parser.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    payload = json.dumps({
        "analytic": audit(precision=args.precision),
        "schur": replay_reserve(args.archive),
    }, indent=2) + "\n"
    if args.out:
        args.out.write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
