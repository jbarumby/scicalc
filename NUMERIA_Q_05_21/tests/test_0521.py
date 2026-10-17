import unittest
from classical_reference import modular_inverse_euclid
from generators import generate_minterm, generate_anf, is_prime
from evaluator import verify_exact
from numq_model import Gate, validate_gate

class TestNumeriaQ0521(unittest.TestCase):
    def test_primes(self):
        self.assertTrue(is_prime(17)); self.assertFalse(is_prime(15))
    def test_classical_inverse(self):
        for x in range(1,7): self.assertEqual((x*modular_inverse_euclid(x,7))%7,1)
    def test_noninvertible_raises(self):
        with self.assertRaises(ValueError): modular_inverse_euclid(0,7)
    def test_minterm_exact(self):
        self.assertTrue(verify_exact(generate_minterm(7),7)["passed"])
    def test_anf_exact(self):
        self.assertTrue(verify_exact(generate_anf(7),7)["passed"])
    def test_ancillas_clean(self):
        self.assertTrue(verify_exact(generate_anf(11),11)["ancillas_clean"])
    def test_invalid_gate(self):
        with self.assertRaises(ValueError): validate_gate(Gate("BAD",(0,)),2)
    def test_reject_nonprime(self):
        with self.assertRaises(ValueError): generate_anf(9)
if __name__=="__main__": unittest.main()
