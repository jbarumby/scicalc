"""Spectral diagnostics for exact dyadic V.40 witnesses; no full zero sum."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

from flint import acb, arb, arb_mat, ctx

from .assembly import T, WIDTH, _bessels, band_gram
from .interval_ldl import dyadic


def bounds(ball):
    return [dyadic(ball.lower()), dyadic(ball.upper())]


def load_witness(record):
    detail = record["preconditioned_congruence"]
    support = detail["unproved_point_pivot"]
    raw = detail["rayleigh_check"]["witness_dyadic"]
    truncated = [arb(int(m)) * arb(2) ** e for m, e in raw]
    if not all(x.is_exact() for x in truncated) or any(not x.is_zero() for x in truncated[support:]):
        raise ValueError("witness is not exact or exceeds the stated support")
    values = truncated[:support]
    norm = sum((v * v for v in values), arb(0)).sqrt()
    if not norm > 0:
        raise ValueError("witness has zero norm")
    return [v / norm for v in values]


def diagnostic(path, bits, zero_count):
    raw = path.read_bytes()
    data = json.loads(raw)
    if data["schema"] != "NUMERIA.P93V40.PREFIXES.1":
        raise ValueError("expected a certified V.40 prefix containing exact witnesses")
    if bits < 768 or zero_count < 1:
        raise ValueError("require 768 or more working bits and a positive zero count")
    results = []
    with ctx.workprec(bits):
        a = arb(WIDTH.numerator) / WIDTH.denominator
        delta = 2 * a - arb(16).log()
        if not delta > 0:
            raise ArithmeticError("the support must cross the prime-16 threshold")
        kappa = arb(2).log() / 4
        zeros = acb.zeta_zeros(1, zero_count)
        if any(z.real != arb(1) / 2 or not z.imag > 0 for z in zeros):
            raise ArithmeticError("a reported zero could not be identified on the positive critical line")
        for row in data["sectors"]:
            parity = row["parity"]
            detail = row["preconditioned_congruence"]
            if not detail or "witness_dyadic" not in detail.get("rayleigh_check", {}):
                continue
            values = load_witness(row)
            support = len(values)
            n = range(parity, 2 * support, 2)
            fweights = [(-1) ** k * c * (2 * a * (2 * m + 1)).sqrt() for k, (m, c) in enumerate(zip(n, values))]

            def ft(t):
                # Actual transform has the common unit phase i**parity.
                if t.is_zero():
                    return fweights[0] if parity == 0 else arb(0)
                js = _bessels(a * t, n[-1])
                return sum((c * js[m] for c, m in zip(fweights, n)), arb(0))

            samples = [abs(ft(z.imag)) for z in zeros]
            concentration = band_gram([(T, 1)], support)[parity]
            column = arb_mat(support, 1, values)
            tail_mass = 1 - (column.transpose() * (concentration * column))[0, 0]
            physical_samples = {str(t): bounds(abs(ft(arb(t)))) for t in (0, 4, 10, 30, 70, 100, 500)}
            boundary_values = [c * ((2 * m + 1) / (2 * a)).sqrt() for m, c in zip(n, values)]
            endpoint = abs(sum(boundary_values, arb(0)))
            derivative_cap = sum((abs(c) * m**2 for m, c in zip(n, boundary_values)), arb(0)) / a
            # Markov bounds each Legendre derivative. By parity, both
            # strips have the same supremum and each has length delta.
            edge = endpoint + delta * derivative_cap
            strip_cost = (2 * kappa * delta * edge**2).upper()
            results.append({
                "parity": parity,
                "support_modes": support,
                "fourier_absolute_values_at_selected_ordinary_frequencies": physical_samples,
                "absolute_value_at_first_zero": bounds(samples[0]),
                "maximum_sample_upper_dyadic": dyadic(max(b.upper() for b in samples)),
                "two_sided_partial_zero_sum_dyadic": bounds(2 * sum((b**2 for b in samples), arb(0))),
                "spectral_l2_mass_above_cutoff_dyadic": bounds(tail_mass),
                "absolute_value_at_support_endpoint": bounds(endpoint),
                "prime16_absolute_cost_upper_dyadic": dyadic(strip_cost),
            })
    if not results:
        raise ValueError("the prefix replay contains no exact witness")
    return {
        "schema": "NUMERIA.P93V40.WITNESS_SPECTRAL_DIAGNOSTIC.1",
        "input_sha256": sha256(raw).hexdigest(),
        "precision_bits": bits,
        "normalization": "L2-normalized witness in orthonormal Legendre basis, Fourier exp(itx) without prefactor",
        "zero_samples_per_parity": zero_count,
        "frequency_cutoff": T,
        "scope": "finite diagnostics for the witnesses, not a certificate for the full sum over zeros",
        "sectors": results,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("proof", type=Path)
    parser.add_argument("--bits", type=int, default=768)
    parser.add_argument("--zeros", type=int, default=30)
    args = parser.parse_args()
    print(json.dumps(diagnostic(args.proof, args.bits, args.zeros), indent=2))


if __name__ == "__main__":
    main()
