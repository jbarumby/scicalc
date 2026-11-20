"""Replay the arithmetic inputs to the local Weil positivity argument."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from numeria.p244.package.p244_analytic_gates import exp_lower

from NUMERIA_ES8P250.adaptive_signed_tail import run as replay_signed_tail
from NUMERIA_ES8P251.gate_b.finite_head import certify_aggregate
from NUMERIA_ES8P251.p251_transfer_gate import PINNED_P250_SHA256


ROOT = Path(__file__).resolve().parents[2]
GATE_B_DATA = ROOT / "NUMERIA_ES8P251/gate_b/data"
GATE_B_CERT_SHA = "86d9c56d8bcd7308120163a3bcd450558eca02b30d7348d050fee2847526f3b2"
GATE_B_MATRIX_SHA = "8854f0ff5c7c70aca84b3941ddda63e7f5c3cf4046515dfaa5edc6e0c983c4ae"
HALF_WIDTH = Fraction(161, 200)
CUTOFF = 309
BASE_MARGIN = Fraction(16186, 10**22)
INCREMENT = Fraction(5, 10**22)


def authenticated_input(path, expected_digest):
    raw = path.read_bytes()
    if sha256(raw).hexdigest() != expected_digest:
        raise ArithmeticError(f"SHA-256 mismatch in {path.name}")
    return raw


def check_model_floor(report):
    """Verify that a Gate B replay certifies the stated fixed-beta model."""
    parameters = {
        "schema": "NUMERIA.ES8P251.GATE_B.FINITE_HEAD_CERTIFICATE.1",
        "input_GL60_aggregate_sha256": GATE_B_MATRIX_SHA,
        "L": str(HALF_WIDTH),
        "cutoff": CUTOFF,
        "beta_reference_cutoff": 504,
        "modes_per_sector": 259,
        "panels_per_sector": 618,
        "both_finite_heads_positive": True,
    }
    for name, expected in parameters.items():
        value = report.get(name)
        if type(value) is not type(expected) or value != expected:
            raise ArithmeticError(f"incompatible Gate B {name}")
    sectors = report["sectors"]
    if (set(sectors) != {"even", "odd"}
            or any(sector["certified_positive"] is not True for sector in sectors.values())):
        raise ArithmeticError("both Gate B sectors must be positive")
    transfer = report["schur_transfer"]
    if transfer["certified_for_reduced_model"] is not True:
        raise ArithmeticError("Gate B Schur transfer is inconclusive")
    margin = Fraction(transfer["model_floor_lower_exact"])
    if not margin > BASE_MARGIN:
        raise ArithmeticError("Gate B margin is below the bound used in the bridge")
    return margin


def local_extension(symbol_lower, *, half_width_increment=INCREMENT):
    """Check the rational inequalities in the one-sided dilation proof."""
    increment = Fraction(half_width_increment)
    symbol_lower = Fraction(symbol_lower)
    if increment < 0:
        raise ValueError("the local extension is one-sided")

    half_width_right = HALF_WIDTH + increment
    # exp(2 L0)>5 and exp(2 L)<=exp(5/3)<7 since e<3 and 3^5<7^3.
    if not (exp_lower(2 * HALF_WIDTH, 8) > 5 and half_width_right <= Fraction(5, 6)
            and 3**5 < 7**3):
        raise ArithmeticError("the proof's prime-power window is not established")

    log2, log3, log5 = Fraction(7, 10), Fraction(11, 10), Fraction(81, 50)
    sqrt2, sqrt3, sqrt5 = Fraction(7, 5), Fraction(17, 10), Fraction(11, 5)
    primes = ((2, log2, sqrt2), (3, log3, sqrt3), (5, log5, sqrt5))
    if not all(exp_lower(log_bound, 8) > prime and sqrt_bound**2 < prime
               for prime, log_bound, sqrt_bound in primes):
        raise ArithmeticError("prime slope input inequalities do not hold")
    prime_slope = 2 * (log2**2 / sqrt2 + log3**2 / sqrt3 + log2**2 + log5**2 / sqrt5)
    slope_cap = Fraction(11, 2)
    if not prime_slope < slope_cap:
        raise ArithmeticError("prime slope budget exceeded")

    global_slope_cap = 9 + slope_cap
    common_band_slope_cap = max(2 * global_slope_cap, 3 + CUTOFF * slope_cap)
    # The pole derivative has norm strictly below 24 for L<=1.
    total_lipschitz_cap = common_band_slope_cap / HALF_WIDTH + 24
    perturbation_factor = 2140
    if not total_lipschitz_cap < perturbation_factor:
        raise ArithmeticError("common-band and pole perturbation budget exceeded")
    sliver_lower = symbol_lower - global_slope_cap * CUTOFF * increment / HALF_WIDTH
    if not sliver_lower > 0:
        raise ArithmeticError("positivity on the newly added frequency band is inconclusive")
    loss = perturbation_factor * increment
    floor = BASE_MARGIN - loss
    if not floor > 0:
        raise ArithmeticError("model coercivity after dilation is inconclusive")

    return {
        "half_width_left_exact": str(HALF_WIDTH),
        "half_width_right_exact": str(half_width_right),
        "half_width_increment_exact": str(increment),
        "prime_powers": [2, 3, 4, 5],
        "prime_slope_upper_exact": str(prime_slope),
        "symbol_derivative_upper_exact": str(global_slope_cap),
        "common_band_and_pole_derivative_upper_exact": str(total_lipschitz_cap),
        "common_band_and_pole_loss_upper_exact": str(loss),
        "new_band_symbol_lower_exact": str(sliver_lower),
        "model_floor_lower_exact": str(floor),
    }


def run_audit():
    matrix_path = GATE_B_DATA / "P251_GATE_B_GL60.json.gz"
    stored = json.loads(authenticated_input(GATE_B_DATA / "P251_GATE_B_CERT.json", GATE_B_CERT_SHA))
    authenticated_input(matrix_path, GATE_B_MATRIX_SHA)
    check_model_floor(stored)
    check_model_floor(certify_aggregate(matrix_path))

    tail = replay_signed_tail()
    raw = (json.dumps(tail, indent=2, ensure_ascii=False) + "\n").encode("utf-8")
    if sha256(raw).hexdigest() != PINNED_P250_SHA256:
        raise ArithmeticError("replayed P250 certificate did not match the pinned digest")
    if (tail["L"] != str(HALF_WIDTH) or tail["beta_frozen_at"] != 504
            or tail["certified_positive_frequency_range"] != [CUTOFF, "infinity"]):
        raise ArithmeticError("the P250 tail does not match Gate B")
    tail_lower = Fraction(int(tail["min_lower_dyadic_numerator"]), 2**tail["interval_scale_power_of_two"])
    if not tail_lower > Fraction(1794, 10**7):
        raise ArithmeticError("the P250 lower bound at the seam is insufficient")

    return {
        "schema": "NUMERIA.ES8P251.GATE_C.WEIL_BRIDGE.1",
        "gate_b_GL60_sha256": GATE_B_MATRIX_SHA,
        "gate_b_certificate_sha256": GATE_B_CERT_SHA,
        "p250_recomputed_sha256": PINNED_P250_SHA256,
        "cutoff": CUTOFF,
        "beta_frozen_at": 504,
        "fixed_window_model_floor_lower_exact": str(BASE_MARGIN),
        "p250_tail_symbol_lower_exact": str(tail_lower),
        "local_extension": local_extension(Fraction(1794, 10**7)),
        "test_function_class": "complex C_c^infty(R) with support in [-L,L]",
        "analytic_inputs": [
            "Guinand–Weil explicit formula, Fourier inversion and Plancherel",
            "Binet digamma lower bound and trigamma series",
            "Legendre–Bessel identities and Gauss–Legendre Bernstein error bound",
        ],
        "cofinal_extension_established": False,
        "Riemann_hypothesis_proved": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True, help="write the exact arithmetic certificate")
    args = parser.parse_args()
    output = json.dumps(run_audit(), indent=2, ensure_ascii=False) + "\n"
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(output, encoding="utf-8")


if __name__ == "__main__":
    main()
