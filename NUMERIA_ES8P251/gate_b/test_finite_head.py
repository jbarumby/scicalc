"""Independent and adversarial checks for the ES8-P251 Gate B replay."""

from copy import deepcopy
from fractions import Fraction
import gzip
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from flint import acb, arb, arb_mat, ctx

from NUMERIA_ES8P251.gate_b.finite_head import aggregate_chunks
from NUMERIA_ES8P251.gate_b.finite_head_certify import (
    assemble_finite_head,
    certify_center_weyl,
    gl60_error_budget,
)
from NUMERIA_ES8P251.gate_b.finite_head_quadrature import chunk_report, quadrature_chunk, spherical_bessel_j_sequence
from NUMERIA_ES8P251.gate_b.historical_audit import P245_ROOTS, P245_ROOTS_SHA, authenticated_bytes
from NUMERIA_ES8P251.gate_b.mpfr_audit import _read_trace


def decode_dyadic(pair):
    mantissa, exponent = (int(value) for value in pair)
    return mantissa * Fraction(2) ** exponent


def _write_chunk(folder, name, report):
    path = Path(folder) / name
    path.write_bytes(gzip.compress(json.dumps(report).encode("utf-8"), mtime=0))
    return path


def _amplitudes(z, length):
    """Real parity amplitudes of the first four Fourier–Legendre modes."""
    sine, cosine = z.sin_cos()
    return (
        length.sqrt() * sine / z,
        (3 * length).sqrt() * (sine / z**2 - cosine / z),
        -(5 * length).sqrt() * ((3 / z**3 - 1 / z) * sine - 3 * cosine / z**2),
        -(7 * length).sqrt() * ((15 / z**4 - 6 / z**2) * sine + (1 / z - 15 / z**3) * cosine),
    )


def _scalar_panel(panel):
    """Direct GL60 sum from the rational P245B root brackets and weights."""
    length = arb(161) / 100
    even, odd = arb_mat(2, 2), arb_mat(2, 2)
    brackets = json.loads(authenticated_bytes(P245_ROOTS, P245_ROOTS_SHA))
    denominator = arb(2) ** brackets["root_denominator_power_of_two"]
    degree = brackets["root_degree"]
    for root in brackets["roots"]:
        left, right = [arb(int(root[f"{key}_numerator"])) / denominator for key in ("lower", "upper")]
        relative_node = (left + right) / 2 + arb(0, (right - left) / 2)
        # At a root of P_n, (1-x²) P_n'(x) = n P_(n-1)(x).
        weight = 2 * (1 - relative_node**2) / (degree * relative_node.legendre_p(degree - 1)) ** 2
        t = (2 * panel + 1 + relative_node) / 4
        amplitude = _amplitudes(length * t / 2, length)
        # Same beta(504) from P237, rewritten to remove the prime mass and pi.
        multiplier = acb(arb(1) / 4, t / 2).digamma().real - arb(252).log() + arb(1) / 504
        for power, prime in ((2, 2), (3, 3), (4, 2), (5, 5)):
            phase = (t * arb(power).log()) / 2
            multiplier += (4 * arb(prime).log() / arb(power).sqrt()) * phase.sin() ** 2
        scaled = multiplier * weight / (4 * arb.pi())
        for parity, matrix in enumerate((even, odd)):
            for i in range(2):
                for j in range(2):
                    matrix[i, j] += scaled * amplitude[2 * i + parity] * amplitude[2 * j + parity]
    return even, odd


def _pole_integrals():
    """Compute pole coefficients by polynomial moments of exp(x/2)."""
    length = arb(161) / 100
    z = length / 4
    e_plus, e_minus = z.exp(), (-z).exp()
    powers = [(e_plus - e_minus) / z]
    for degree in (1, 2, 3):
        powers.append((e_plus - (-1) ** degree * e_minus - degree * powers[-1]) / z)
    legendre = (powers[0], powers[1], (3 * powers[2] - powers[0]) / 2, (5 * powers[3] - 3 * powers[1]) / 2)
    return [legendre[n] * ((2 * n + 1) * length).sqrt() / 2 for n in range(4)]


class FiniteQuadratureTests(unittest.TestCase):
    def test_first_and_last_panels_match_independent_scalar_reference(self):
        for panel in (0, 617):
            with self.subTest(panel=panel):
                obtained = quadrature_chunk(panel, panel + 1, modes=2)
                with ctx.workprec(384):
                    expected = _scalar_panel(panel)
                    for parity, (actual, reference) in enumerate(zip(obtained, expected, strict=True)):
                        for i in range(2):
                            for j in range(2):
                                self.assertTrue(actual[i, j].overlaps(reference[i, j]), (parity, i, j))

    def test_adjacent_chunks_cover_seam_and_grid_validation_rejects_outside(self):
        first = quadrature_chunk(616, 617, modes=1)
        second = quadrature_chunk(617, 618, modes=1)
        joined = quadrature_chunk(616, 618, modes=1)
        with ctx.workprec(256):
            for left, right, complete in zip(first, second, joined, strict=True):
                self.assertTrue((left[0, 0] + right[0, 0]).overlaps(complete[0, 0]))
        for start, stop in ((-1, 0), (617, 617), (618, 619)):
            with self.subTest(start=start, stop=stop):
                with self.assertRaises(ValueError):
                    quadrature_chunk(start, stop, modes=1)

    def test_serialized_chunk_retains_exact_coverage_and_declares_unapplied_terms(self):
        report = chunk_report(616, 618, modes=2)
        self.assertEqual(report["panel_range"], [616, 618])
        self.assertEqual(report["cutoff"], 309)
        self.assertEqual(report["beta_cutoff"], 504)
        self.assertEqual(report["GL60_nodes"], 120)
        self.assertFalse(report["beta_identity_poles_and_analytic_remainder_added"])
        self.assertEqual(set(report["sectors"]), {"even", "odd"})
        for entries in report["sectors"].values():
            self.assertEqual([(i, j) for i, j, _ in entries], [(0, 0), (0, 1), (1, 1)])
            for _, _, endpoints in entries:
                lo, hi = map(decode_dyadic, endpoints)
                self.assertLessEqual(lo, hi)

    def test_downward_recurrence_encloses_independent_native_bessel_values(self):
        with ctx.workprec(700):
            for argument in (arb(1) / 100, arb(161 * 309) / 200):
                sequence = spherical_bessel_j_sequence(argument, 517)
                self.assertEqual(len(sequence), 518)
                for degree in (0, 1, 52, 258, 516, 517):
                    expected = (arb.pi() / (2 * argument)).sqrt() * argument.bessel_j(arb(2 * degree + 1) / 2)
                    self.assertTrue(sequence[degree].overlaps(expected), degree)

    def test_matrix_free_actions_overlap_the_dense_product(self):
        vectors = ([1, -2, 3], [-1, 1, 0])
        dense = quadrature_chunk(616, 618, modes=3, batch_panels=1)
        actions = quadrature_chunk(616, 618, modes=3, vectors=vectors)
        for matrix, vector, action in zip(dense, vectors, actions, strict=True):
            reference = matrix * arb_mat(3, 1, vector)
            self.assertEqual(action.ncols(), 1)
            for direct, expected in zip(action.entries(), reference.entries(), strict=True):
                self.assertTrue(direct.overlaps(expected))
        with self.assertRaises(ValueError):
            quadrature_chunk(616, 617, modes=3, vectors=([1], [1]))


class FiniteCertificateTests(unittest.TestCase):
    def test_remainder_uses_panel_count_and_archived_bernstein_gl60_constant(self):
        for sector, norms in (("even", (4, 6, 8)), ("odd", (6, 8, 10))):
            with self.subTest(sector=sector):
                bound = gl60_error_budget(sector, modes=3)
                factor = Fraction(618 * 50, 3**120)
                self.assertEqual(bound["panels"], 618)
                self.assertEqual(bound["panel_factor"], Fraction(50, 3**120))
                self.assertEqual(bound["all_panels_factor"], factor)
                self.assertEqual(bound["max_entry_bound"], factor * max(norms) ** 2)
                self.assertEqual(bound["row_sum_bound"], factor * max(norms) * sum(norms))

    def test_assembly_uses_fixed_beta_and_opposite_pole_signs_once(self):
        with ctx.workprec(320):
            quadrature = arb_mat([[arb(1) / 8, arb(1) / 32], [arb(1) / 32, -arb(1) / 64]])
            mass = sum(2 * arb(p).log() / arb(q).sqrt() for q, p in ((2, 2), (3, 3), (4, 2), (5, 5)))
            beta = (arb(504) / (2 * arb.pi())).log() - arb(1) / 504 - mass
            poles = _pole_integrals()
            for sector, parity, sign in (("even", 0, 1), ("odd", 1, -1)):
                with self.subTest(sector=sector):
                    result, _ = assemble_finite_head(quadrature, sector)
                    for i in range(2):
                        for j in range(2):
                            expected = quadrature[i, j] + (beta if i == j else arb(0))
                            expected += sign * 2 * poles[2 * i + parity] * poles[2 * j + parity]
                            self.assertTrue(result[i, j].overlaps(expected), (i, j))

    def test_center_weyl_certifies_only_when_both_independent_gates_pass(self):
        symmetric = arb_mat([[2, 1], [1, 2]])
        success = certify_center_weyl(symmetric, analytic_row_error=Fraction(1, 4))
        self.assertTrue(success["certified_positive"])
        self.assertEqual(success["status"], "CERTIFIED_POSITIVE")
        self.assertEqual(success["positive_pivot_count"], 2)
        self.assertEqual(decode_dyadic(success["perturbation_upper_dyadic"]), Fraction(1, 4))
        self.assertTrue(0 < decode_dyadic(success["margin_lower_dyadic"]) <= Fraction(3, 4))

        bad_center = certify_center_weyl(arb_mat([[1, 2], [2, 1]]))
        self.assertFalse(bad_center["certified_positive"])
        self.assertEqual(bad_center["status"], "MIDPOINT_LDL_INCONCLUSIVE")
        self.assertEqual(bad_center["positive_pivot_count"], 1)

        large_error = certify_center_weyl(symmetric, analytic_row_error=Fraction(1))
        self.assertFalse(large_error["certified_positive"])
        self.assertEqual(large_error["status"], "WEYL_INCONCLUSIVE")
        self.assertLessEqual(decode_dyadic(large_error["margin_lower_dyadic"]), 0)

        wide_box = certify_center_weyl(arb_mat([[arb(2, 1), 1], [1, 2]]))
        self.assertFalse(wide_box["certified_positive"])
        self.assertEqual(wide_box["status"], "WEYL_INCONCLUSIVE")
        self.assertGreaterEqual(decode_dyadic(wide_box["perturbation_upper_dyadic"]), 1)

    def test_center_weyl_rejects_nonsymmetric_boxes(self):
        with self.assertRaises(ValueError):
            certify_center_weyl(arb_mat([[2, 0], [1, 2]]))

    def test_mpfr_trace_rejects_invalid_pivots_and_inverse(self):
        header = "# index low_mantissa low_exp high_mantissa high_exp\n"
        correct = "0 1 0 2 0\n# inverse_frobenius_squared_upper 3 0\n"
        self.assertEqual(_read_trace((header + correct).encode(), 1), (["1", "0"], ["3", "0"]))
        for bad in (correct.replace("0 1", "1 1", 1), correct.replace("0 1", "0 0", 1),
                    correct.replace("0 1 0 2", "0 3 0 2", 1), correct.replace("upper 3", "upper 0")):
            with self.subTest(trace=bad), self.assertRaises((ArithmeticError, ValueError)):
                _read_trace((header + bad).encode(), 1)


class FiniteCoverageTests(unittest.TestCase):
    def test_historical_root_digest_must_match(self):
        with self.assertRaises(ArithmeticError):
            authenticated_bytes(P245_ROOTS, "0" * 64)

    def test_chunk_aggregation_encloses_exact_dyadic_sums(self):
        reports = [chunk_report(start, start + 1, cutoff=1, modes=2) for start in (0, 1)]
        with TemporaryDirectory() as folder, ctx.workprec(320):
            paths = [_write_chunk(folder, f"{index}.json.gz", report) for index, report in enumerate(reports)]
            obtained, records = aggregate_chunks(paths, cutoff=1, modes=2)
            self.assertEqual(len(records), 2)
            for sector, matrix in zip(("even", "odd"), obtained, strict=True):
                for left, right in zip(*(report["sectors"][sector] for report in reports), strict=True):
                    i, j = left[:2]
                    expected_lo, expected_hi = (sum(map(decode_dyadic, endpoints)) for endpoints in zip(left[2], right[2]))
                    self.assertLessEqual(decode_dyadic(matrix[i, j].lower().man_exp()), expected_lo)
                    self.assertGreaterEqual(decode_dyadic(matrix[i, j].upper().man_exp()), expected_hi)
                    self.assertEqual(matrix[i, j].mid(), matrix[j, i].mid())

    def test_chunk_aggregation_rejects_incomplete_or_mislabelled_inputs(self):
        left, right = [chunk_report(start, start + 1, cutoff=1, modes=2) for start in (0, 1)]
        wrong_beta = deepcopy(left)
        wrong_beta["beta_cutoff"] = 309
        wrong_order = deepcopy(left)
        wrong_order["sectors"]["even"].reverse()
        with TemporaryDirectory() as folder:
            first, second, bad_beta, bad_order = [
                _write_chunk(folder, f"{index}.json.gz", report)
                for index, report in enumerate((left, right, wrong_beta, wrong_order))
            ]
            for inputs in ([first], [first, first, second], [bad_beta, second], [bad_order, second]):
                with self.subTest(inputs=[path.name for path in inputs]):
                    with self.assertRaises(ValueError):
                        aggregate_chunks(inputs, cutoff=1, modes=2)


if __name__ == "__main__":
    unittest.main()
