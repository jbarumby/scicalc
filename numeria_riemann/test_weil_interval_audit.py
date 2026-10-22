"""Testes da auditoria intervalar finita de Weil.

Rápidos por padrão (~30 s). A checagem ponta a ponta contra zeros reais da zeta
é lenta e só roda com NUMERIA_SLOW=1 (precisa de mpmath).
"""
import os
import random
import unittest

from flint import acb, arb, ctx

from numeria_riemann import audit_diagnostics as diag
from numeria_riemann import verify_certificate as ver
from numeria_riemann import weil_interval_audit as v2
from numeria_riemann.weil_interval_audit import _omega, _sinc, run

try:
    import mpmath as mp
except ImportError:  # pragma: no cover
    mp = None

CERT = "certified_positive_definite_finite_matrix"


class _Precision(unittest.TestCase):
    def setUp(self):
        self._old = ctx.prec
        ctx.prec = 192

    def tearDown(self):
        ctx.prec = self._old


class WeilIntervalAuditTests(unittest.TestCase):
    def test_sinc_interval_contains_zero_singularity(self):
        value = _sinc(arb(0, "0.01"))
        self.assertTrue(value.is_finite())
        self.assertTrue(value.lower() < arb("1.001"))
        self.assertTrue(value.upper() > arb("0.999"))

    def test_omega_interval_is_finite(self):
        value = _omega(arb("1 +/- 0.001"), terms=32)
        self.assertTrue(value.is_finite())

    def test_legacy_rangesum_returns_explicit_status(self):
        report = run(L_value=0.45, N=4, R_value=40, step_value=0.05, terms=32,
                     precision=128, archimedean="rangesum")
        self.assertIn(report["status"], {CERT, "inconclusive_interval_ldl"})
        self.assertEqual(report["projected_dimension"], 2)
        self.assertEqual(report["parameters"]["quadrature_cells"], 800)
        pivots = report["interval_ldl_pivots"]
        self.assertGreaterEqual(len(pivots), 1)
        self.assertLessEqual(len(pivots), 2)
        if report["status"] == CERT:
            self.assertEqual(len(pivots), 2)
        self.assertIn(2, report["active_prime_powers"])
        self.assertIn("not global Weil positivity", report["claim_scope"])


class WeilV2Tests(_Precision):
    def test_phase_is_integer_and_zero_for_odd_difference(self):
        self.assertEqual((v2._phase(1, 2), v2._phase(4, 1)), (0, 0))
        self.assertEqual((v2._phase(3, 3), v2._phase(1, 3), v2._phase(5, 1)), (1, -1, 1))
        self.assertTrue(all(isinstance(v2._phase(i, j), int) for i in range(1, 7) for j in range(1, 7)))

    def test_zero_shift_overlap_is_identity_gram_matrix(self):
        # regressão do fator 2 (normalização do overlap), ver CORRECTED_RESULTS_2026-10-20
        L = arb("0.45")
        for i in range(1, 5):
            for j in range(1, 5):
                value = v2._overlap(arb(0), L, i, j)
                self.assertTrue(value.contains(1 if i == j else 0), (i, j, value))

    def test_rump_rejects_when_margin_too_small(self):
        A = [[arb(1) + arb(0, "0.6") if i == j else arb(0, "0.6") for j in range(3)] for i in range(3)]
        self.assertFalse(v2._rump_certificate(A)[0])

    def test_adaptive_integral_consistent_with_legacy_rangesum(self):
        N, L, R = 4, arb("0.45"), arb("40")
        new = v2._archimedean_matrix_adaptive(L, N, R, 40)
        old, _ = v2._archimedean_matrix_rangesum(L, N, R, "0.02", "40", 128)
        for i in range(N):
            for j in range(N):
                if v2._phase(i + 1, j + 1) == 0:
                    self.assertEqual((float(new[i][j].lower()), float(new[i][j].upper())), (0.0, 0.0))
                    continue
                a, b = old[i][j], new[i][j]
                self.assertLessEqual(abs(float(a.mid()) - float(b.mid())), float(a.rad()) + float(b.rad()))
                self.assertLess(float(b.rad()), 1e-8)

    def test_full_run_small_case_reports_both_certificates(self):
        rep = run(N=6, precision=192)
        self.assertEqual(rep["status"], CERT)
        self.assertEqual(len(rep["interval_ldl_pivot_upper_bounds"]), rep["projected_dimension"])
        self.assertEqual(diag.check_consistency(rep), [])

    def test_diagnostics_handles_nan_and_flags_inconsistency(self):
        bad = {"status": diag.CERT, "projected_dimension": 3,
               "interval_ldl_pivot_lower_bounds": ["0.5", "NaN", "-1"]}
        self.assertTrue(diag.check_consistency(bad))
        self.assertIn("inválido", diag.summarize(bad))


class WeilExportAndVerifierTests(_Precision):
    def test_export_is_outward_rounded_on_random_balls(self):
        from decimal import Decimal
        from fractions import Fraction
        random.seed(5)
        for _ in range(400):
            mid = arb(random.uniform(-1, 1)) * arb(10) ** random.randint(-80, 40)
            rad = arb(10) ** random.randint(-120, 10) * random.random()
            x = arb(mid.mid(), rad.upper()) if random.random() < 0.9 else mid.mid()
            d = v2._ball_json(x, 40)
            self.assertLessEqual(Fraction(Decimal(d["lower"])), v2._exact_fraction(x.lower()))
            self.assertGreaterEqual(Fraction(Decimal(d["upper"])), v2._exact_fraction(x.upper()))
            m = int(d["mid_exact"][0]) * arb(2) ** int(d["mid_exact"][1])
            r = int(d["rad_exact"][0]) * arb(2) ** int(d["rad_exact"][1])
            self.assertTrue(arb(m.mid(), r.upper()).contains(x.mid()))

    def test_directed_decimal_never_crosses_the_bound(self):
        # regressão: str(40) arredondava o limite inferior PARA CIMA
        lo = arb(2).sqrt().lower()
        self.assertTrue(arb(v2._lower_decimal(lo)) <= lo or not (arb(v2._lower_decimal(lo)) > lo))
        from fractions import Fraction
        from decimal import Decimal
        self.assertLessEqual(Fraction(Decimal(v2._lower_decimal(lo))), v2._exact_fraction(lo))
        self.assertGreaterEqual(Fraction(Decimal(v2._upper_decimal(lo))), v2._exact_fraction(lo))

    def test_exact_verifier_accepts_run_and_rejects_perturbations(self):
        report = run(N=6, precision=192)
        result = ver.verify(report)
        self.assertTrue(result["verified"])
        self.assertTrue(report["status"] == CERT)
        # alargar todos os enclosures até a margem sumir: o verificador tem que rejeitar
        import copy
        from decimal import Decimal
        wide = copy.deepcopy(report)
        for row in wide["projected_matrix_enclosures"]:
            for e in row:
                e["lower"] = str(Decimal(e["lower"]) - Decimal("0.5"))
                e["upper"] = str(Decimal(e["upper"]) + Decimal("0.5"))
        self.assertFalse(ver.verify(wide)["verified"])

    def test_exact_ldl_detects_indefinite_matrix(self):
        from fractions import Fraction as F
        _, fail = ver.exact_ldl_pivots([[F(1), F(2)], [F(2), F(1)]])
        self.assertEqual(fail, 1)
        pivots, fail = ver.exact_ldl_pivots([[F(2), F(1)], [F(1), F(2)]])
        self.assertIsNone(fail)
        self.assertEqual(pivots, [F(2), F(3, 2)])

    def test_compare_detects_disjoint_reports(self):
        a = {"projected_matrix_enclosures": [[{"lower": "0", "upper": "1"}]]}
        b = {"projected_matrix_enclosures": [[{"lower": "2", "upper": "3"}]]}
        self.assertTrue(ver.compare(a, a))
        self.assertFalse(ver.compare(a, b))
        details = ver.compare_details(a, b)
        self.assertFalse(details["enclosures_intersect"])
        self.assertEqual(details["max_midpoint_difference"], 2.0)
        self.assertEqual(details["max_entry_width_first"], 1.0)


@unittest.skipIf(mp is None, "mpmath não instalado")
class WeilMpmathCrossChecks(_Precision):
    def test_omega_analytic_matches_mpmath(self):
        mp.mp.dps = 25
        for rv in ("0.01", "0.5", "3", "40", "160"):
            ball = v2._omega_analytic(acb(rv)).real
            ref = mp.re(mp.digamma(mp.mpf("0.25") + 0.5j * mp.mpf(rv))) - mp.log(mp.pi)
            self.assertLess(abs(mp.mpf(ball.mid().str(30, radius=False)) - ref), mp.mpf("1e-20"), rv)
            self.assertLess(float(ball.rad()), 1e-30)

    def test_prime_overlap_matches_direct_quadrature(self):
        mp.mp.dps = 20
        L, a = arb("0.45"), arb(2).log()
        Lm, am = mp.mpf("0.45"), mp.log(2)
        w = lambda k: k * mp.pi / (2 * Lm)

        def directed(k, m):
            return mp.quad(lambda x: mp.sin(w(k) * (x + am)) * mp.sin(w(m) * x), [0, 2 * Lm - am]) / Lm

        for (i, j) in [(1, 1), (1, 2), (2, 5), (4, 3)]:
            ref = (directed(i, j) + directed(j, i)) / 2
            got = v2._overlap(a, L, i, j)
            self.assertLess(abs(mp.mpf(got.mid().str(20, radius=False)) - ref), mp.mpf("1e-14"), (i, j))

    def test_pole_constraint_vectors_match_direct_integrals(self):
        # h(±i/2) = ∫_0^{2L} g_k(x) e^{±(x-L)/2} dx, com g_k = sin(w_k x)/sqrt(L)
        mp.mp.dps = 20
        Lm = mp.mpf("0.45")
        L = arb("0.45")
        B = v2._constraint_nullspace(L, 6)
        pi = arb.pi()
        for sign in (+1, -1):
            for col in range(4):
                total = mp.mpf(0)
                for k in range(1, 7):
                    coeff = mp.mpf(B[k - 1][col].mid().str(20, radius=False))
                    w = k * mp.pi / (2 * Lm)
                    total += coeff * mp.quad(lambda x: mp.sin(w * x) / mp.sqrt(Lm) * mp.exp(sign * (x - Lm) / 2), [0, 2 * Lm])
                self.assertLess(abs(total), mp.mpf("1e-14"), (sign, col))

    def test_tail_proof_lemma1_fourier_transform_matches_closed_form(self):
        # TAIL_PROOF.md, Lema 1: F_k(r) = i e^{-ik pi/2} A_k(r)
        mp.mp.dps = 20
        Lm = mp.mpf("0.45")
        for k, r in [(1, "0.7"), (3, "5.3"), (4, "50")]:
            r = mp.mpf(r)
            w = k * mp.pi / (2 * Lm)
            direct = mp.quad(lambda x: mp.sin(w * x) / mp.sqrt(Lm) * mp.exp(1j * r * (x - Lm)), [0, 2 * Lm])
            amp = 2 * mp.sqrt(Lm) * w / (w + r) * mp.sin(Lm * (r - w)) / (Lm * (r - w))
            self.assertLess(abs(direct - 1j * mp.exp(-1j * k * mp.pi / 2) * amp), mp.mpf("1e-15"))

    def test_tail_proof_lemma4_omega_bounds(self):
        # TAIL_PROOF.md, Lema 4: -3.14 <= Omega(r) <= log r - 0.99 para r >= 1
        mp.mp.dps = 20
        om = lambda r: mp.re(mp.digamma(mp.mpf("0.25") + 0.5j * r)) - mp.log(mp.pi)
        for rv in ["1", "1.01", "1.5", "2", "3", "5", "10", "30", "100", "1e3", "1e5", "1e8"]:
            r = mp.mpf(rv)
            self.assertGreaterEqual(om(r), mp.mpf("-3.14"), rv)
            self.assertLessEqual(om(r), mp.log(r) - mp.mpf("0.99"), rv)

    def test_tail_proof_proposition5_closed_form_integral_and_code_radius(self):
        mp.mp.dps = 20
        R = mp.mpf(160)
        numeric = mp.quad(lambda r: (8 + mp.log(3 * r)) / r ** 4, [R, 1000, mp.inf])
        closed = (8 + mp.log(3 * R) + mp.mpf(1) / 3) / (3 * R ** 3)
        self.assertLess(abs(numeric - closed), mp.mpf("1e-18"))
        tails = v2._tail_bound(arb("0.45"), 6, arb("160"))
        wi = 3 * mp.pi / (2 * mp.mpf("0.45"))
        expected = 16 * wi * wi / (mp.pi * mp.mpf("0.45")) * closed
        self.assertLess(abs(mp.mpf(tails[2][2].mid().str(20, radius=False)) - expected), mp.mpf("1e-15"))
        self.assertEqual(float(tails[0][1].upper()), 0.0)  # i-j ímpar: cauda exatamente 0

    def test_normalization_parseval_links_fk_to_g0_log_pi_term(self):
        # (1/2π)∫_R |F_k|² dr = ∫ h_k² dt = 1  (Gram = I); com h par: (1/π)∫_0^∞ A_k² = 1.
        # É o que faz o termo -g(0) log π da fórmula explícita virar o "-log π" em Ω.
        mp.mp.dps = 20
        Lm = mp.mpf("0.45")

        def amp(r, k):
            w = k * mp.pi / (2 * Lm)
            x = Lm * (r - w)
            return 2 * mp.sqrt(Lm) * w / (w + r) * (mp.sin(x) / x if x != 0 else 1)

        for k in (1, 3, 6):
            pts = [mp.mpf(t) * 160 / 40 for t in range(41)]
            head = mp.quad(lambda r: amp(r, k) ** 2, pts) / mp.pi
            tail = mp.quad(lambda r: amp(r, k) ** 2, [160, 1000, 10000, mp.inf]) / mp.pi
            self.assertLess(abs(head + tail - 1), mp.mpf("1e-7"), k)

    def test_rump_certificate_is_sound_on_random_perturbations(self):
        random.seed(1)
        mp.mp.dps = 30
        n = 4
        base = [[arb("1.0") if i == j else arb("0.1") for j in range(n)] for i in range(n)]
        A = [[base[i][j] + arb(0, "0.05") for j in range(n)] for i in range(n)]
        ok, delta, _, _ = v2._rump_certificate(A)
        self.assertTrue(ok)
        self.assertLessEqual(float(delta), 0.2 + 1e-6)
        for _ in range(200):
            E = mp.matrix(n, n)
            for i in range(n):
                for j in range(i, n):
                    E[i, j] = E[j, i] = mp.mpf(random.uniform(-0.05, 0.05))
            M = mp.matrix([[float(base[i][j].mid()) for j in range(n)] for i in range(n)]) + E
            self.assertGreater(min(mp.eigsy(M, eigvals_only=True)), 0)


@unittest.skipUnless(mp is not None and os.environ.get("NUMERIA_SLOW") == "1",
                     "lento: defina NUMERIA_SLOW=1 (precisa de mpmath)")
class WeilEndToEndZeros(_Precision):
    def test_quadratic_form_equals_sum_over_real_zeta_zeros(self):
        """Cross-check numérico condicionado a zeros conhecidos na linha crítica.

        NÃO verifica zeros fora da linha crítica nem prova a fórmula para zeros
        complexos. A identidade |F_c(gamma)|^2 usa gamma real. A cauda da soma
        de zeros é estimada, não certificada. Este teste é evidência numérica
        de consistência da implementação para os zeros consultados, não uma
        validação independente do critério global de Weil nem da RH.
        """
        mp.mp.dps = 20
        N = 6
        L, R = arb("0.45"), arb("160")
        G = v2._archimedean_matrix_adaptive(L, N, R, 40)
        P, _ = v2._prime_matrix(L, N)
        W = [[G[i][j] + P[i][j] for j in range(N)] for i in range(N)]
        B = v2._constraint_nullspace(L, N)
        A = v2._project(B, W)
        Bm = [[mp.mpf(B[k][i].mid().str(25, radius=False)) for i in range(N - 2)] for k in range(N)]
        Am = mp.matrix([[mp.mpf(A[i][j].mid().str(25, radius=False)) for j in range(N - 2)] for i in range(N - 2)])
        Lm = mp.mpf("0.45")
        w = lambda k: k * mp.pi / (2 * Lm)

        def amp(r, k):
            x = Lm * (r - w(k))
            return 2 * mp.sqrt(Lm) * w(k) / (w(k) + r) * (mp.sin(x) / x if x != 0 else 1)

        def fc2(r, c):
            return sum(c[i - 1] * c[j - 1] * v2._phase(i, j) * amp(r, i) * amp(r, j)
                       for i in range(1, N + 1) for j in range(1, N + 1) if v2._phase(i, j))

        om = lambda r: mp.re(mp.digamma(mp.mpf("0.25") + 0.5j * r)) - mp.log(mp.pi)
        gam = [mp.im(mp.zetazero(n)) for n in range(1, 151)]
        random.seed(3)
        cr = [mp.mpf(random.uniform(-1, 1)) for _ in range(N - 2)]
        c = [sum(Bm[k][i] * cr[i] for i in range(N - 2)) for k in range(N)]
        q = (mp.matrix(cr).T * Am * mp.matrix(cr))[0]
        q += mp.quad(lambda r: fc2(r, c) * om(r) / mp.pi, [160, 400, 1000, 4000, 20000])
        s = 2 * sum(fc2(g, c) for g in gam)
        s += 2 * mp.quad(lambda r: mp.log(r / (2 * mp.pi)) / (2 * mp.pi) * fc2(r, c), [gam[-1], 1000, 4000, 20000])
        self.assertLess(abs(q - s), mp.mpf("1e-6"))


if __name__ == "__main__":
    unittest.main()
