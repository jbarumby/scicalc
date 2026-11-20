"""Independent dense and interval checks for the finite P249 operators."""

import cmath
from fractions import Fraction
import math
import unittest

from flint import acb, arb, ctx, fmpq
import numpy as np

from NUMERIA_ES8P249.interval_audit import (
    audit_constant_rayleigh,
    audit_fft_vector,
    certify_cancellation_shells,
    certify_p236_tail_start,
    prime_shift_lags_arb,
    prime_shift_matvec_arb,
)
from NUMERIA_ES8P249.interval_fft import toeplitz_matvec_arb
from NUMERIA_ES8P249.legendre_tail import modal_tail_bounds
from NUMERIA_ES8P249.p236_replay import panel_cell_lags, replay_panel
from NUMERIA_ES8P249.prime_operator import (
    P236_HALF_WIDTH,
    P236_PRIME_POWERS,
    PrimeShiftOperator,
    ToeplitzFFT,
    dyadic_prime_components,
    fourier_cell_lags,
    prime_symbol,
)


def dense_shift_reference(cells, half_width, powers):
    """Integrate shifted indicator products by intersecting endpoints."""
    edges = np.linspace(-float(half_width), float(half_width), cells + 1)
    width = 2 * float(half_width) / cells
    result = np.zeros((cells, cells))
    for i in range(cells):
        for j in range(cells):
            for power, prime in powers:
                frequency = math.log(power)
                coefficient = math.log(prime) / math.sqrt(power)
                for shift in (frequency, -frequency):
                    left = max(edges[i], edges[j] + shift)
                    right = min(edges[i + 1], edges[j + 1] + shift)
                    result[i, j] += coefficient * max(0.0, right - left) / width
    return result


def arb_from_fraction(value):
    return arb(fmpq(value.numerator, value.denominator))


def decode_dyadic(pair):
    mantissa, exponent = map(int, pair)
    return Fraction(mantissa) * Fraction(2) ** exponent


def symbol_reference(frequency):
    return math.fsum(
        2 * math.log(prime) / math.sqrt(power) * math.cos(frequency * math.log(power))
        for power, prime in ((2, 2), (3, 3), (4, 2), (5, 5))
    )


def p236_panel_nodes(panel):
    for index in range(60):
        root, weight = arb.legendre_p_root(60, index, weight=True)
        yield (arb(2 * panel + 1) + root) / 4, weight / (4 * arb.pi())


def p236_arb_terms(frequency):
    symbol = arb(0)
    correction = arb(0)
    for power, prime in ((2, 2), (3, 3), (4, 2), (5, 5)):
        coefficient = 2 * arb(prime).log() / arb(power).sqrt()
        phase = frequency * arb(power).log()
        symbol += coefficient * phase.cos()
        correction += 2 * coefficient * (phase / 2).sin() ** 2
    return symbol, correction


class MatrixFreeTests(unittest.TestCase):
    def test_fft_matches_dense_reference_for_real_and_complex_vectors(self):
        for size in (1, 2, 5, 13, 32):
            with self.subTest(size=size):
                lags = np.array([(-1) ** k / (k + 1) for k in range(size)])
                operator = ToeplitzFFT(lags)
                dense = np.array([[lags[abs(i - j)] for j in range(size)] for i in range(size)])
                real = np.linspace(-0.7, 1.3, size)
                complex_vector = real + 1j * real[::-1]
                for vector in (real, complex_vector):
                    np.testing.assert_allclose(operator.matvec(vector), dense @ vector, atol=4e-15)

    def test_end_impulse_keeps_linear_convolution(self):
        lags = np.array([2.0, -1.0, 0.0, 0.0, 0.25])
        impulse = np.zeros(len(lags))
        impulse[-1] = 1
        np.testing.assert_allclose(ToeplitzFFT(lags).matvec(impulse), lags[::-1], atol=1e-15)

    def test_prime_shift_matches_spatial_intersections(self):
        cases = (
            (1, P236_HALF_WIDTH, P236_PRIME_POWERS),
            (7, Fraction(7, 20), ((2, 2), (3, 3))),
            (12, P236_HALF_WIDTH, P236_PRIME_POWERS),
            (19, Fraction(1, 4), P236_PRIME_POWERS),
            (23, Fraction(1), ((4, 2),)),
        )
        for cells, half_width, powers in cases:
            with self.subTest(cells=cells):
                dense = dense_shift_reference(cells, half_width, powers)
                vector = np.array([(-1) ** k * (k + 2) / (k + 3) for k in range(cells)])
                operator = PrimeShiftOperator(cells, half_width=half_width, powers=powers)
                np.testing.assert_allclose(operator.lags, dense[:, 0], rtol=0, atol=1e-14)
                np.testing.assert_allclose(operator.matvec(vector), dense @ vector, rtol=0, atol=2e-14)

    def test_fourier_cell_lags_match_products_of_analytic_integrals(self):
        cells = 9
        half_width = Fraction(13, 17)
        frequencies = [0.0, 0.75, 4.0, 13.5]
        weights = [1.5, -0.4, -1.1, 0.7]
        edges = np.linspace(-float(half_width), float(half_width), cells + 1)
        width = edges[1] - edges[0]
        gram = np.zeros((cells, cells))
        for frequency, weight in zip(frequencies, weights):
            if frequency == 0:
                values = np.full(cells, math.sqrt(width))
            else:
                primitive = [cmath.exp(1j * frequency * x) / (1j * frequency) for x in edges]
                values = np.diff(primitive) / math.sqrt(width)
            gram += weight * np.outer(values, np.conjugate(values)).real
        lags = fourier_cell_lags(cells, half_width, frequencies, weights)
        vector = np.linspace(-1, 2, cells)
        np.testing.assert_allclose(lags, gram[:, 0], atol=3e-15)
        np.testing.assert_allclose(ToeplitzFFT(lags).matvec(vector), gram @ vector, atol=1e-14)


class OscillationTests(unittest.TestCase):
    def test_prime_symbol_preserves_sign_and_exact_triangle_at_zero(self):
        points = [0.0, 0.125, 3.25, 79.785, 504.0]
        expected = [symbol_reference(t) for t in points]
        self.assertLess(min(expected), 0)
        np.testing.assert_allclose(prime_symbol(points), expected, atol=1.5e-13)

    def test_dyadic_components_telescope_and_respect_taylor_remainder(self):
        points = np.array([0, 0.5, 0.999, 1, 4.625, 19.571, 251.317, 504])
        truth = np.array([symbol_reference(t) for t in points])
        second_derivative_bound = math.fsum(
            2 * math.log(prime) / math.sqrt(power) * math.log(power) ** 2
            for power, prime in P236_PRIME_POWERS
        )
        for max_level in (0, 1, 4, 8):
            with self.subTest(max_level=max_level):
                components = dyadic_prime_components(points, max_level)
                self.assertEqual(len(components), max_level + 2)
                np.testing.assert_allclose(np.sum(components, axis=0), truth, atol=1.5e-13)
                for level in range(max_level + 1):
                    approximation = np.sum(components[:level + 1], axis=0)
                    bound = second_derivative_bound * 2.0 ** (-2 * level) / 8
                    self.assertLessEqual(float(np.max(np.abs(approximation - truth))), bound + 1.5e-13)
                self.assertGreater(float(np.max(np.abs(components[-1]))), 0)

    def test_shell_envelopes_cover_signed_values_and_register_cancellation(self):
        report = certify_cancellation_shells(endpoint=8, subdivisions_per_unit=16)
        shells = report["shells"]
        self.assertEqual([shell["range"] for shell in shells], [[0, 1], [1, 2], [2, 4], [4, 8]])
        self.assertFalse(shells[0]["strict_reduction"])
        self.assertTrue(any(shell["strict_reduction"] for shell in shells[1:]))
        self.assertGreaterEqual(float(decode_dyadic(shells[0]["envelope_upper_dyadic"])), symbol_reference(0))
        for shell in shells:
            start, stop = shell["range"]
            envelope = float(decode_dyadic(shell["envelope_upper_dyadic"]))
            sampled = max(abs(symbol_reference(t)) for t in np.linspace(start, stop, 129))
            self.assertLessEqual(sampled, envelope + 3e-15)
            if shell["strict_reduction"]:
                self.assertLess(decode_dyadic(shell["ratio_upper_dyadic"]), 1)

    def test_signed_tail_certificate_keeps_beta_and_covers_the_whole_range(self):
        certificate = certify_p236_tail_start()
        start, stop = map(Fraction, certificate["certified_range"])
        subdivisions = certificate["subdivisions_per_unit"]
        self.assertEqual(subdivisions, 128)
        self.assertEqual(start, 309)
        self.assertEqual(stop, certificate["original_cutoff"])
        self.assertEqual(stop, 504)
        self.assertEqual(certificate["certified_subcells"], (stop - start) * subdivisions)
        self.assertEqual(certificate["certified_subcells"], 24960)
        strict_gap = decode_dyadic(certificate["strict_gap_lower_dyadic"])
        self.assertGreater(strict_gap, 0)
        inconclusive = certificate["first_inconclusive_cell"]
        left, right = (Fraction(x, subdivisions) for x in inconclusive["interval_numerators"])
        self.assertEqual(right - left, Fraction(1, subdivisions))
        self.assertLessEqual(right, start)
        self.assertLessEqual(decode_dyadic(inconclusive["lower_bound_dyadic"]), 0)
        with ctx.workprec(256):
            beta = (arb(504) / (2 * arb.pi())).log() - arb(1) / 504
            beta -= sum(2 * arb(prime).log() / arb(power).sqrt() for power, prime in P236_PRIME_POWERS)
            self.assertLessEqual(arb_from_fraction(decode_dyadic(certificate["beta_lower_dyadic"])), beta.lower())
            self.assertLessEqual(beta.upper(), arb_from_fraction(decode_dyadic(certificate["beta_upper_dyadic"])))
            # Check independently at both edges and an interior point using
            # A-S(t)=sum 2*a*sin²(t*log(q)/2), with no subtraction of S(t).
            for point in (start, (start + stop) / 2, stop):
                t = arb_from_fraction(point)
                _, correction = p236_arb_terms(t)
                gap = (t / 504).log() + arb(1) / 504 - 1 / t + correction
                self.assertLessEqual(arb_from_fraction(strict_gap), gap.lower())


class ArbAuditTests(unittest.TestCase):
    def test_p248_constant_rayleigh_is_independent_of_grid_size(self):
        one_cell = audit_constant_rayleigh(1)
        lower_reference = Fraction('1.05798304576502875357447158043965')
        upper_reference = Fraction('1.05798304576502875357447158043966')
        for cells in (7, 128):
            with self.subTest(cells=cells):
                report = audit_constant_rayleigh(cells)
                self.assertEqual(report['direct_overlap_dyadic'], one_cell['direct_overlap_dyadic'])
                lower, upper = map(decode_dyadic, report['cell_fft_dyadic'])
                self.assertLess(lower_reference, lower)
                self.assertLess(upper, upper_reference)
                self.assertLess(decode_dyadic(report['abs_difference_upper_dyadic']), Fraction(2) ** -140)

    def test_outward_fft_contains_exact_dense_complex_actions(self):
        dimension = 7
        radius = arb(1) / 64
        lag_centers = [arb((-1) ** k * (k + 1)) / 16 for k in range(dimension)]
        vector_centers = [acb(arb(k - 2) / 8, arb(5 - k) / 16) for k in range(dimension)]
        lags = [arb(center, radius) for center in lag_centers]
        vector = [acb(arb(center.real, radius), arb(center.imag, radius)) for center in vector_centers]
        enclosures = toeplitz_matvec_arb(lags, vector)
        for direction in (-1, 0, 1):
            exact_lags = [center + direction * radius for center in lag_centers]
            exact_vector = [center + acb(direction * radius, -direction * radius) for center in vector_centers]
            for row, enclosure in enumerate(enclosures):
                direct = sum(
                    (exact_lags[abs(row - col)] * value for col, value in enumerate(exact_vector)), acb(0)
                )
                self.assertTrue(enclosure.contains(direct))

    def test_exact_one_cell_ball_and_actual_binary_fft_error_are_enclosed(self):
        vector = [0.3]
        report = audit_fft_vector(vector, precision_bits=128)
        coordinate = report["coordinates"][0]
        with ctx.workprec(256):
            width = 2 * arb_from_fraction(P236_HALF_WIDTH)
            exact_lag = sum(
                2 * arb(prime).log() / arb(power).sqrt() * (1 - arb(power).log() / width)
                for power, prime in ((2, 2), (3, 3), (4, 2), (5, 5))
            )
            exact = exact_lag * arb_from_fraction(Fraction(vector[0]))
            lower = arb_from_fraction(decode_dyadic(coordinate["exact_lower_dyadic"]))
            upper = arb_from_fraction(decode_dyadic(coordinate["exact_upper_dyadic"]))
            self.assertLessEqual(lower, exact.lower())
            self.assertLessEqual(exact.upper(), upper)
            actual = float(PrimeShiftOperator(1).matvec(vector)[0])
            error = abs(exact - arb_from_fraction(Fraction(actual)))
            bound = arb_from_fraction(decode_dyadic(coordinate["abs_error_upper_dyadic"]))
            self.assertLessEqual(error.upper(), bound)
            self.assertEqual(coordinate["abs_error_upper_dyadic"], report["max_abs_error_upper_dyadic"])

    def test_support_cutoff_is_exact_in_arb(self):
        half_width = Fraction(1, 4)
        lags = prime_shift_lags_arb(7, half_width=half_width)
        applied = prime_shift_matvec_arb([0.5, 0.25, 1, 2, 3, 4, 5], half_width=half_width)
        self.assertTrue(all(lag.is_zero() for lag in lags))
        self.assertTrue(all(value.is_zero() for value in applied))

    def test_p236_signed_cell_lag_matches_independent_complex_integrals(self):
        cells = 11
        panel = 503
        lags = panel_cell_lags(panel, cells, bits=160)
        with ctx.workprec(256):
            half_width = arb(161) / 200
            width = 2 * half_width / cells
            left = -half_width
            reference = [arb(0), arb(0)]
            for t, weight in p236_panel_nodes(panel):
                symbol, _ = p236_arb_terms(t)
                values = []
                for cell in (0, 1):
                    a = left + cell * width
                    b = a + width
                    values.append(
                        (acb(0, t * b).exp() - acb(0, t * a).exp()) / acb(0, t * width.sqrt())
                    )
                for lag in (0, 1):
                    reference[lag] += weight * symbol * (values[lag] * values[0].conjugate()).real
            for lag in (0, 1):
                self.assertTrue(lags[lag].overlaps(reference[lag]))

    def test_legendre_panel_replay_matches_closed_forms_in_both_parities(self):
        for panel, degree, companion in ((0, 0, 2), (1007, 1, 3)):
            with self.subTest(panel=panel, degree=degree):
                replay = replay_panel(panel, 17, (degree, companion), bits=160)
                entries = {tuple(entry["degrees"]): entry for entry in replay["entries"]}
                entry = entries[(degree, degree)]
                with ctx.workprec(256):
                    half_width = arb(161) / 200
                    prime, full = arb(0), arb(0)
                    for t, weight in p236_panel_nodes(panel):
                        z = half_width * t
                        if degree == 0:
                            amplitude = (2 * half_width).sqrt() * z.sin() / z
                        else:
                            amplitude = (6 * half_width).sqrt() * (z.sin() / z ** 2 - z.cos() / z)
                        symbol, correction = p236_arb_terms(t)
                        baseline = acb(arb(1) / 4, t / 2).digamma().real - arb(252).log() + arb(1) / 504
                        prime += weight * amplitude ** 2 * symbol
                        full += weight * amplitude ** 2 * (baseline + correction)
                    self.assertTrue(entry["legendre_prime"].overlaps(prime))
                    self.assertTrue(entry["p236_full_quadrature"].overlaps(full))
                    if degree == 0:
                        self.assertTrue(entry["grid_bound"].is_zero())
                        self.assertTrue(entry["grid_residual"].contains(0))
                    else:
                        self.assertLess(abs(entry["grid_residual"]).upper(), entry["grid_bound"].lower())


class ModalTailTests(unittest.TestCase):
    def test_reduced_cutoff_and_archived_cutoff_meet_rational_budgets(self):
        coupling_budget = Fraction(1, 10**180)
        tail_budget = Fraction(1, 10**190)
        for cutoff, modes in ((309, 259), (504, 377)):
            for parity in (0, 1):
                with self.subTest(cutoff=cutoff, modes=modes, parity=parity):
                    coupling_sq, deviation = modal_tail_bounds(cutoff, modes, parity)
                    self.assertIsInstance(coupling_sq, Fraction)
                    self.assertIsInstance(deviation, Fraction)
                    self.assertTrue(0 < coupling_sq < coupling_budget)
                    self.assertTrue(0 < deviation < tail_budget)

    def test_deviation_bound_with_258_modes_does_not_meet_budget_for_both_parities(self):
        tail_budget = Fraction(1, 10**190)
        deviations = [modal_tail_bounds(309, 258, parity)[1] for parity in (0, 1)]
        self.assertGreaterEqual(max(deviations), tail_budget)


if __name__ == "__main__":
    unittest.main()
