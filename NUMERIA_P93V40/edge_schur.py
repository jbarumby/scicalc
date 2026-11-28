"""Replay the structured prime-16 Schur reduction in directed Arb."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from flint import arb, arb_mat, ctx

from .assembly import ARCHIVE, LOCALIZATION, T, WIDTH, assemble
from .edge_direction import directional_cross_majorant
from .edge_schur_algebra import finite_rank_schur
from .edge_sinc import endpoint_coefficients, endpoint_norm_squares
from .prefix import _head_digest
from .precondition import certify_by_congruence
from .transfer import RHO, WIDTH_CAP, ball, count_mass, dyadic, low_symbol_cap, pole_audit


OMEGA = 1_600_000
ETA = Fraction(8, 15)
YOUNG = Fraction(1, 16)
COMPLEMENT_SIZES = (200, 400, 600, 1000, 2000, 3000, 3309, 3400, 3500)
SCALAR_GATE_PROBES = (Fraction("4e-72"), Fraction("2e-72"))
FINITE_REPLAY = Path(__file__).with_name("edge_schur_finite_200.json")


def edge_constants():
    a, eta = ball(WIDTH), ball(ETA)
    delta = 2 * a - arb(16).log()
    if not (delta > 0 and 2 * a < arb(17).log()):
        raise ArithmeticError("support width does not cross only the prime power 16")
    cap = ball(WIDTH_CAP)
    if not (a < cap and 2 * cap * cap.exp() < 15 and 2 * cap * cap.sinh() < 9):
        raise ArithmeticError("the support lies outside the pole and deficit bound hypotheses")
    kappa = arb(2).log() / 4
    omega = arb(OMEGA)
    far_lower = (omega / (2 * arb.pi())).log() - 2 / omega**2 - 1 / (3 * omega)
    far_lower -= count_mass() + 2
    far_slack = (far_lower - kappa * (1 + 1 / eta)).lower()
    if not far_slack > 0:
        raise ArithmeticError("the far spectrum does not absorb the strip high-frequency cost")
    # Within each parity the two squared endpoint evaluations coincide.
    coefficient = (4 * kappa * (1 + eta) * delta).upper()
    scalar = (4 * kappa * (1 + eta) * delta**3 * omega**3 / (9 * arb.pi())).upper()
    return coefficient, scalar, far_slack


def _block(matrix, rows, columns):
    return arb_mat([[matrix[i, j] for j in columns] for i in rows])


def _proof(matrix):
    report = certify_by_congruence(matrix, [matrix.nrows()], 0, absolute_bound=True)
    if report["status"] != "certified":
        return {"status": report["status"], "diagnostic": report}
    return {"status": report["status"], **report["prefixes"][-1]}


def _upper(value):
    return dyadic(value.upper())


def _lower(value):
    return dyadic(value.lower())


def _interval(value):
    return [_lower(value), _upper(value)]


def _constants_record(coefficient, scalar, far_slack):
    return {
        "omega": OMEGA,
        "strip_young_parameter": str(ETA),
        "width_exact": str(WIDTH),
        "edge_coefficient_upper_dyadic": dyadic(coefficient),
        "edge_scalar_remainder_upper_dyadic": dyadic(scalar),
        "far_reserve_after_absorption_lower_dyadic": dyadic(far_slack),
    }


def finite_certificate(modes, split, bits):
    if not 0 < split < modes <= LOCALIZATION:
        raise ValueError("require a nonempty split inside the localized polynomial span")
    with ctx.workprec(bits):
        sectors = assemble(modes, bits)
        return finite_from_assembly(sectors, modes, split, bits)


def finite_from_assembly(sectors, modes, split, bits):
    with ctx.workprec(bits):
        gamma, epsilon, far_slack = edge_constants()
        endpoints = endpoint_coefficients(modes, width=WIDTH, omega=OMEGA)
        infinite_split = 3309
        budgets = {row["parity"]: row for row in complement_budget()["results"]
                   if row["modes"] == infinite_split}
        records = []
        first, second = range(split), range(split, modes)
        for parity, (source, endpoint) in enumerate(zip(sectors, endpoints, strict=True)):
            base = arb_mat(source["head"])
            for i in range(modes):
                base[i, i] -= epsilon
            a, c = _block(base, first, first), _block(base, first, second)
            d = _block(base, second, second)
            u = arb_mat(split, 1, endpoint[:split])
            v = arb_mat(modes-split, 1, endpoint[split:])
            tail_proof = _proof(d)
            row = {"parity": parity, "base_matrix_sha256": _head_digest(base),
                   "finite_complement": tail_proof}
            budget = budgets[parity]
            cross_weight = _decode(budget["cross_majorant_multiplier_upper_dyadic"])
            edge_weight = _decode(budget["required_head_endpoint_weight_upper_dyadic"])
            endpoint_vector = arb_mat(modes, 1, endpoint)
            directed_head = (base - directional_cross_majorant(modes, infinite_split, parity)*cross_weight
                             - (endpoint_vector*endpoint_vector.transpose())*edge_weight)
            row["infinite_directional_head_prefix"] = {
                "infinite_split_per_parity": infinite_split,
                "corrected_matrix_sha256": _head_digest(directed_head),
                **_proof(directed_head),
            }
            if modes >= 125:
                trial = SCALAR_GATE_PROBES[parity]
                probe = certify_by_congruence(base, [modes], trial)
                row["scalar_gate_probe"] = {
                    "head_trial_exact": str(trial),
                    "status": probe["status"],
                    **({"support_modes": probe["unproved_point_pivot"],
                        "rayleigh_check": probe["rayleigh_check"]}
                       if "rayleigh_check" in probe else {}),
                }
            if tail_proof["status"] == "certified":
                result = finite_rank_schur(a, c, d, u, v, gamma)
                effective = result["effective_endpoint"]
                schur_proof = _proof(result["schur"])
                row.update({
                    "tail_loading_interval_dyadic": _interval(result["tail_loading"]),
                    "rank_denominator_lower_dyadic": _lower(result["denominator"]),
                    "schur_rank_multiplier_upper_dyadic": _upper(result["multiplier"]),
                    "effective_endpoint_norm_squared_upper_dyadic":
                        _upper((effective.transpose() * effective)[0, 0]),
                    "corrected_schur_sha256": _head_digest(result["schur"]),
                    "corrected_schur": schur_proof,
                })
                if schur_proof["status"] == "certified":
                    elimination_norm = sum((x**2 for x in result["corrected_elimination"].entries()),
                                           arb(0)).sqrt().upper()
                    tail_floor = (_decode(tail_proof["absolute_margin_lower_dyadic"])
                                  * result["denominator"]).lower()
                    schur_floor = _decode(schur_proof["absolute_margin_lower_dyadic"])
                    form_floor = (min(tail_floor, schur_floor)/(1+elimination_norm)**2).lower()
                    if not form_floor > 0:
                        raise ArithmeticError("the directed full finite gap is inconclusive")
                    row["full_finite_form_gap_lower_dyadic"] = dyadic(form_floor)
            records.append(row)
    return {
        "schema": "NUMERIA.P93V40.EDGE_SCHUR_FINITE.1",
        "scope": "full Weil form restricted to the stated finite polynomial span",
        "modes_per_parity": modes, "split_per_parity": split, "precision_bits": bits,
        **_constants_record(gamma, epsilon, far_slack),
        "sectors": records,
        "certified": all(row.get("corrected_schur", {}).get("status") == "certified"
                         for row in records),
    }


def complement_budget(bits=512):
    """Bound the infinite tail and its cross block, conditional on a corrected head."""
    with ctx.workprec(bits):
        gamma, epsilon, far_slack = edge_constants()
        c = ball(json.loads((ARCHIVE / "hybrid39.json").read_text())["c"])
        cap = arb(18)
        low_symbol_cap(count_mass())
        rho = ball(RHO)
        d_rho = rho - 1 / rho
        z = ball(WIDTH) * T * d_rho
        half = (z / 2).exp()
        factor = 2 * (half - 1 + arb.pi() * (z.exp()-half) / (2*z))
        factor /= d_rho * (1-rho**-4)
        endpoints = endpoint_coefficients(max(COMPLEMENT_SIZES), width=WIDTH, omega=OMEGA)
        norms = endpoint_norm_squares(width=WIDTH, omega=OMEGA)
        young = ball(YOUNG)
        rows = []
        for n in COMPLEMENT_SIZES:
            for parity, (vector, norm) in enumerate(zip(endpoints, norms, strict=True)):
                b = (norm - sum((x**2 for x in vector[:n]), arb(0))).upper()
                if not b >= 0:
                    raise ArithmeticError("negative upper enclosure for the endpoint tail norm")
                g = min(arb(1), (factor * rho**(-2*(2*n+parity))).upper())
                pole_cross, pole_tail = map(ball, pole_audit(n))
                d = (2 - c - cap*g - pole_tail - epsilon).lower()
                tail_after_edge = (d - gamma*b).lower()
                row = {
                    "modes": n, "parity": parity, "localization": min(n, LOCALIZATION),
                    "low_band_tail_concentration_upper_dyadic": dyadic(g),
                    "endpoint_tail_norm_squared_upper_dyadic": dyadic(b),
                    "tail_reserve_lower_dyadic": dyadic(d),
                    "corrected_tail_reserve_lower_dyadic": dyadic(tail_after_edge),
                    "tail_certified_by_bound": bool(tail_after_edge > 0),
                }
                if tail_after_edge > 0:
                    base_cross_bound = (cap*g.sqrt() + pole_cross).upper()
                    endpoint_weight = (gamma*d*(1+young)/tail_after_edge).upper()
                    cross_weight = 1/d + gamma*(1+1/young)*b/(d*tail_after_edge)
                    row.update({
                        "base_cross_norm_upper_dyadic": dyadic(base_cross_bound),
                        "cross_majorant_multiplier_upper_dyadic": _upper(cross_weight),
                        "required_corrected_head_gap_upper_dyadic":
                            _upper(cross_weight*base_cross_bound**2),
                        "required_head_endpoint_weight_upper_dyadic": dyadic(endpoint_weight),
                    })
                rows.append(row)
    return {
        "schema": "NUMERIA.P93V40.EDGE_SCHUR_COMPLEMENT.1", "precision_bits": bits,
        **_constants_record(gamma, epsilon, far_slack),
        "scope": "conditional infinite Schur test; no positivity assertion for unassembled heads",
        "hybrid_theta_exact": "0", "rho_exact": str(RHO),
        "schur_young_parameter": str(YOUNG), "results": rows,
    }


def _decode(pair):
    mantissa, exponent = pair
    return arb(int(mantissa)) * arb(2)**exponent


def directional_witness_report(path=FINITE_REPLAY, bits=768):
    """Test the infinite directional Schur bound on reproducible finite witnesses."""
    raw = path.read_bytes()
    source = json.loads(raw)
    if source["schema"] != "NUMERIA.P93V40.EDGE_SCHUR_FINITE.1":
        raise ValueError("expected a physical finite Schur replay")
    split = 3309
    budgets = {row["parity"]: row for row in complement_budget()["results"]
               if row["modes"] == split}
    records = []
    with ctx.workprec(bits):
        for sector in source["sectors"]:
            parity = sector["parity"]
            probe = sector["scalar_gate_probe"]
            result = probe.get("rayleigh_check", {})
            if not result.get("proven_below_requested_margin"):
                continue
            support = probe["support_modes"]
            coordinates = [_decode(pair) for pair in result["witness_dyadic"]]
            if not (all(x.is_exact() for x in coordinates)
                    and all(x.is_zero() for x in coordinates[support:])):
                raise ArithmeticError("the scalar-gate witness has the wrong support")
            vector = arb_mat(support, 1, coordinates[:support])
            norm = (vector.transpose()*vector)[0, 0]
            if not norm > 0:
                raise ArithmeticError("a zero witness is invalid")
            lower, upper = map(_decode, result["rayleigh_numerator_dyadic"])
            q = lower.union(upper) / norm
            bound = directional_cross_majorant(support, split, parity)
            cross = (vector.transpose()*bound*vector)[0, 0] / norm
            u = arb_mat(support, 1, endpoint_coefficients(support, width=WIDTH, omega=OMEGA)[parity])
            endpoint = (u.transpose()*vector)[0, 0]**2 / norm
            budget = budgets[parity]
            cross_cost = _decode(budget["cross_majorant_multiplier_upper_dyadic"]) * cross
            edge_cost = _decode(budget["required_head_endpoint_weight_upper_dyadic"]) * endpoint
            scalar_gate = _decode(budget["required_corrected_head_gap_upper_dyadic"])
            slack = (q-cross_cost-edge_cost).lower()
            records.append({
                "parity": parity, "witness_support": support,
                "base_rayleigh_interval_dyadic": _interval(q),
                "directional_cross_cost_upper_dyadic": _upper(cross_cost),
                "effective_border_cost_upper_dyadic": _upper(edge_cost),
                "directional_schur_slack_lower_dyadic": dyadic(slack),
                "fails_scalar_sufficient_threshold": bool(q < scalar_gate),
                "passes_directional_sufficient_test_on_this_vector": bool(slack > 0),
            })
    if not records:
        raise ValueError("expected at least one dyadic scalar-gate witness")
    return {
        "schema": "NUMERIA.P93V40.EDGE_SCHUR_DIRECTIONS.1",
        "source_sha256": sha256(raw).hexdigest(),
        "precision_bits": bits, "infinite_split_per_parity": split,
        "scope": "validated directions only, not a full-head or operator certificate",
        "sectors": records,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--modes", type=int, default=200)
    parser.add_argument("--split", type=int, default=100)
    parser.add_argument("--bits", type=int, default=768)
    output_kind = parser.add_mutually_exclusive_group()
    output_kind.add_argument("--complement", action="store_true")
    output_kind.add_argument("--directional", action="store_true")
    parser.add_argument("--head-proof", type=Path, default=FINITE_REPLAY)
    io = parser.add_mutually_exclusive_group(required=True)
    io.add_argument("--out", type=Path)
    io.add_argument("--verify", type=Path)
    args = parser.parse_args()
    expected = json.loads(args.verify.read_text()) if args.verify else None
    if expected:
        if expected["schema"] == "NUMERIA.P93V40.EDGE_SCHUR_COMPLEMENT.1":
            report = complement_budget(expected["precision_bits"])
        elif expected["schema"] == "NUMERIA.P93V40.EDGE_SCHUR_DIRECTIONS.1":
            report = directional_witness_report(args.head_proof, expected["precision_bits"])
        elif expected["schema"] == "NUMERIA.P93V40.EDGE_SCHUR_FINITE.1":
            report = finite_certificate(expected["modes_per_parity"],
                                        expected["split_per_parity"], expected["precision_bits"])
        else:
            raise ValueError("unknown Schur certificate schema")
        if report != expected:
            raise ValueError("the directed Schur replay differs from the certificate")
        print(f"Verified {args.verify}")
    else:
        if args.complement:
            report = complement_budget(args.bits)
        elif args.directional:
            report = directional_witness_report(args.head_proof, args.bits)
        else:
            report = finite_certificate(args.modes, args.split, args.bits)
        args.out.write_text(json.dumps(report, indent=2) + "\n")
        print(report.get("certified", report["schema"]))


if __name__ == "__main__":
    main()
