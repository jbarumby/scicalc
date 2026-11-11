"""Testes de regressão para o laboratório de cálculo avançado do SciCalc."""
import importlib.util
import json
import pathlib
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("scicalc_pro_core", ROOT / "scicalc-pro-core.py")
core = importlib.util.module_from_spec(spec)
spec.loader.exec_module(core)

def solve(**request):
    return json.loads(core.scicalc_pro_solve(json.dumps(request)))

class AdvancedCalculusTests(unittest.TestCase):
    def test_non_elementary_antiderivative(self):
        result = solve(action="integral", expr="exp(-x^2)")
        self.assertTrue(result["ok"])
        self.assertIn("erf", result["result_text"])
        self.assertTrue(result["constant"])

    def test_definite_fresnel_integral(self):
        result = solve(action="integral", expr="sin(x^2)", lower="0", upper="1")
        self.assertTrue(result["ok"])
        self.assertNotEqual(result["status"], "Integral formal")

    def test_infinite_bounds(self):
        result = solve(action="integral", expr="exp(-x)", lower="0", upper="inf")
        self.assertEqual(result["result_text"], "1")

    def test_removable_zero_zero_limit(self):
        result = solve(action="limit", expr="(exp(x)-1-x)/x^2", target="0", direction="both")
        self.assertEqual(result["result_text"], "1/2")

    def test_one_sided_disagreement(self):
        result = solve(action="limit", expr="1/x", target="0", direction="both")
        self.assertEqual(result["status"], "Não existe (bilateral)")
        self.assertIn("left_tex", result)
        self.assertIn("right_tex", result)

    def test_limit_at_infinity(self):
        result = solve(action="limit", expr="x^(1/x)", target="inf", direction="both")
        self.assertEqual(result["result_text"], "1")

    def test_higher_derivative(self):
        result = solve(action="derivative", expr="x^2*exp(x)", order=3, point="0")
        self.assertEqual(result["point_result_tex"], "6")

    def test_trigonometric_limit_has_actual_algebra(self):
        r=solve(action="limit",expr="(1-cos(x))/x^2",target="0",direction="both")
        self.assertEqual(r["result_text"], "1/2")
        steps=r["steps"]
        self.assertGreaterEqual(len(steps),5)
        full_tex=" ".join(x.get("tex","") for x in steps)
        self.assertIn(r"\sin^2",full_tex)
        self.assertIn(r"\frac12",full_tex)
        self.assertIn(r"\boxed{\frac12}",full_tex)
        self.assertNotIn("não foi possível gerar",full_tex.lower())

    def test_limit_hopital_two_steps_verified(self):
        r=solve(action="limit",expr="(exp(x)-1-x)/x^2",target="0",direction="both")
        self.assertEqual(r["result_text"],"1/2")
        titles=[x["title"] for x in r["steps"]]
        self.assertTrue(any("L’Hôpital" in t for t in titles))
        self.assertGreaterEqual(len(titles),4)
        self.assertTrue(any("quociente" in t.lower() for t in titles))

    def test_limit_polynomial_cancellation_steps(self):
        r=solve(action="limit",expr="(x^2-9)/(x-3)",target="3",direction="both")
        self.assertEqual(r["result_text"],"6")
        titles=[x["title"] for x in r["steps"]]
        self.assertIn("Cancelar fatores comuns",titles)

    def test_limit_at_infinity_has_power_division(self):
        r=solve(action="limit",expr="(2*x+1)/x",target="inf",direction="both")
        self.assertEqual(r["result_text"],"2")
        self.assertTrue(any("maior potência" in x["title"] for x in r["steps"]))

    def test_derivative_product_has_successive_derivatives(self):
        r=solve(action="derivative",expr="x^2*exp(x)",order=3)
        self.assertTrue(any("regra do produto" in x["title"].lower() for x in r["steps"]))
        self.assertTrue(any("3ª vez" in x["title"] for x in r["steps"]))
        self.assertTrue(any("Confirmar" in x["title"] for x in r["steps"]))

    def test_integral_by_parts_has_intermediate_identity(self):
        r=solve(action="integral",expr="x*exp(x)")
        self.assertEqual(r["verified"],True)
        self.assertTrue(any("por partes" in x["title"].lower() for x in r["steps"]))
        joined=" ".join(x.get("tex","") for x in r["steps"])
        self.assertIn(r"\int xe^x",joined)

    def test_integral_first_step_latex_command_is_separated(self):
        """Regression: '\\int' + 'e^{...}' must not form invalid '\\inte'."""
        for expr in ("exp(-x^2)", "exp(x)", "x^2", "sin(x)", "cos(x)"):
            with self.subTest(expr=expr):
                response=solve(action="integral",expr=expr)
                first=response["steps"][0]["tex"]
                self.assertTrue(first.startswith(r"\int\,"),first)
                self.assertNotIn(r"\inte",first)
                self.assertIn(r"\,dx",first)

    def test_symbolic_primitives_do_not_claim_numeric_approximation(self):
        r=solve(action="integral",expr="exp(-x^2)")
        self.assertIsNone(r.get("approx"))
        self.assertEqual(r.get("engine_revision"),"20261111c")
        self.assertIn("erf",r["result_text"])

    def test_numeric_constant_has_actual_approximation(self):
        r=solve(action="integral",expr="exp(-x)",lower="0",upper="1")
        self.assertIsNotNone(r.get("approx"))
        self.assertFalse(any(c.isalpha() for c in r["approx"] if c not in 'eE'))

    def test_bounded_integral_first_step_latex(self):
        response=solve(action="integral",expr="exp(-x)",lower="0",upper="1")
        first=response["steps"][0]["tex"]
        self.assertTrue(first.startswith(r"\int_{0}^{1}\,"),first)
        self.assertIn(r"e^{- x}",first)
        self.assertNotIn(r"\inte",first)

    def test_gaussian_integral_uses_function_erf(self):
        r=solve(action="integral",expr="exp(-x^2)")
        self.assertTrue(any("erf" in x.get("tex","") for x in r["steps"]))

    def test_unrecognized_method_does_not_pretend_to_prove(self):
        r=solve(action="limit",expr="x^(1/x)",target="inf",direction="both")
        self.assertEqual(r["result_text"],"1")
        self.assertTrue(any("não foi possível gerar" in x.get("note","").lower()
            for x in r["steps"]))

    def test_invalid_input_rejected(self):
        with self.assertRaises(ValueError):
            solve(action="derivative", expr="__import__('os').system('x')", order=1)

    def test_no_false_numeric_integral_over_pole(self):
        result = solve(action="integral", expr="1/x", lower="-1", upper="1")
        self.assertNotEqual(result.get("approx"), "0")

if __name__ == "__main__":
    unittest.main()
