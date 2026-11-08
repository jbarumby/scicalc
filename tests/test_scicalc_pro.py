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

    def test_invalid_input_rejected(self):
        with self.assertRaises(ValueError):
            solve(action="derivative", expr="__import__('os').system('x')", order=1)

    def test_no_false_numeric_integral_over_pole(self):
        result = solve(action="integral", expr="1/x", lower="-1", upper="1")
        self.assertNotEqual(result.get("approx"), "0")

if __name__ == "__main__":
    unittest.main()
