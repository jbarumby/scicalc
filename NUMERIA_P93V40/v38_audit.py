"""Replay the attached V.38 A1/A4 Schur tables with integers and Arb."""

import argparse
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
import json
from math import factorial
from pathlib import Path

from flint import acb, arb, ctx
import numpy as np

from .spectral_replay import CSV_PATH, CSV_SHA256, SCALE, TUPLE, NineTermSymbol, dyadic, read_rows, source_cells


INPUTS = Path(__file__).with_name("v38_inputs")
ATTACHMENT = INPUTS / "P93V38_A1_A4_resultados.json"
ATTACHMENT_SHA256 = "7a2282afd44635232cff2249a22d0f5f3becaf70e203e95961d12cf0ef83d1d6"
OPT_SHA256 = "a96514ec7e7b06eac39ca6c2d055b3979ad8416209d84ad994695e425d686d2b"
WIDTH_CAP = Fraction(1387, 1000)
GRID = 2**20
BETA_GRID = (Fraction(3, 2), Fraction(13, 8), Fraction(7, 4), Fraction(15, 8), Fraction(2))
CUTOFFS = (1600, 3100, 4500, 6400, 10000)
POLE_EPS_DENOMINATOR = 2**2000
HEAD_MARGIN = Fraction(1, 10**60)
TARGET = Fraction(9, 10) * HEAD_MARGIN


def _require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def _read_json(path, expected_sha=None):
    payload = path.read_bytes()
    digest = sha256(payload).hexdigest()
    if expected_sha is not None and digest != expected_sha:
        raise ValueError(f"SHA-256 mismatch for {path.name}")
    return json.loads(payload), digest


def _ball(value):
    exact = Fraction(value)
    return arb(exact.numerator) / exact.denominator


def components(rows, beta, cutoff):
    """Lower the beta=2 cap in exact millionths, then group adjacent cells."""
    reduction = (2 - beta) * 10**6
    if reduction < 0 or reduction.denominator != 1:
        raise ValueError("beta must be on the grid of millionths and at most two")
    result = []
    for left, right, deficit in rows:
        deficit -= reduction.numerator
        if left < 16 * cutoff or deficit <= 0:
            continue
        if result and result[-1][1] == left:
            first, _, mass, height = result[-1]
            result[-1] = first, right, mass + deficit, max(height, deficit)
        else:
            result.append((left, right, deficit, deficit))
    return result


def _kernel(left, right):
    distance = np.maximum(1, np.maximum(left[:, None] - right[None, :], left[None, :] - right[:, None]))
    total = left[:, None] + left[None, :]
    cap = (1387 * GRID + 999) // 1000
    grid16 = 16 * GRID
    return np.minimum(cap, (grid16 + distance - 1) // distance) + (grid16 + total - 1) // total


def schur_bound(kernel, masses, weights):
    """Compute max row ratio exactly, switching products to Python integers if needed."""
    if len(weights) != len(masses) or any(type(w) is not int or w <= 0 for w in weights):
        raise ValueError("invalid positive integral Schur weights")
    weighted = [int(m) * w for m, w in zip(masses, weights, strict=True)]
    if sum(weighted) * int(kernel.max()) < 2**63:
        sums = kernel @ np.array(weighted, dtype=np.int64)
    else:
        sums = kernel.astype(object) @ np.array(weighted, dtype=object)
    denominator = 333 * GRID * 16_000_000
    return max(Fraction(106 * int(s), denominator * w) for s, w in zip(sums, weights, strict=True))


class _DirectDigammaSymbol(NineTermSymbol):
    def arch_lower(self, left):
        # Euler's digamma series gives monotonicity of Re psi at fixed positive
        # real part and increasing nonnegative imaginary part, term by term.
        return (acb(arb(1) / 4, left / 2).digamma().real - self.log_pi).lower()


def early_spectrum():
    """Verify the range preceding V.40's archived spectral replay at 4500."""
    digest = sha256()
    accepted = [0, 0]
    minima = [None, None]
    with ctx.workprec(128):
        symbol = _DirectDigammaSymbol()

        def verify(left, right, cap, value):
            margin = symbol.cell_margin(arb(left) / SCALE, arb(right) / SCALE, cap, False)
            if margin > 0:
                group = int(value != 0)
                accepted[group] += 1
                if minima[group] is None or margin < minima[group]:
                    minima[group] = margin
                digest.update(TUPLE.pack(left, right, value))
                return
            split = (left + right) // 2
            _require(split != left, f"early spectral enclosure did not close at {left}/{SCALE}")
            verify(left, split, cap, value)
            verify(split, right, cap, value)

        for left, right, value in source_cells(CSV_PATH, 1600, 4500):
            verify(left, right, arb(value) / 10**6, value)
        return {
            "start": 1600, "stop": 4500, "working_precision_bits": 128,
            "omitted_accepted_leaves": accepted[0], "listed_accepted_leaves": accepted[1],
            "omitted_margin_lower_dyadic": dyadic(minima[0]),
            "listed_margin_lower_dyadic": dyadic(minima[1]),
            "ordered_leaf_partition_sha256": digest.hexdigest(),
        }


@lru_cache(maxsize=None)
def _pole_ok(modes):
    value = factorial(2 * modes)
    return value > 48 * POLE_EPS_DENOMINATOR and value**2 > 72 * POLE_EPS_DENOMINATOR


class _GateReplay:
    """Apply the historical Schur formula with exact constants and Arb exp."""

    def __init__(self):
        self._concentrations = {}
        self.sector_tests = 0
        self.display_underflow = {"g": 0, "b2": 0}

    def _concentration(self, cutoff, rho, modes, parity):
        key = cutoff, rho, modes, parity
        if key not in self._concentrations:
            prefactor = (8 * WIDTH_CAP * cutoff / 3) / ((1 - 1 / rho)**2 * (1 - rho**-4))
            exponent = WIDTH_CAP * cutoff * (rho - 1 / rho)
            degree = 2 * (2 * modes + parity)
            # Evaluate the product in logarithmic form to avoid underflow and
            # avoid large numerator/denominator powers in the rational replay.
            self._concentrations[key] = _ball(prefactor) * (_ball(exponent) - degree * _ball(rho).log()).exp()
        return self._concentrations[key]

    def check(self, cutoff, rho, modes, parity, penalty, mixing):
        g = self._concentration(cutoff, rho, modes, parity)
        pole = arb(1) / POLE_EPS_DENOMINATOR
        b2 = 2 * _ball(mixing**2) * g + 2 * pole**2
        reserve = _ball(2 - penalty - TARGET) - 50 * g - pole
        slack = _ball(HEAD_MARGIN - TARGET) * reserve - b2
        if reserve > 0 and slack > 0:
            return True, g, b2
        _require(reserve < 0 or slack < 0, f"indeterminate gate at T={cutoff}, N={modes}, parity={parity}")
        return False, g, b2

    def first(self, cutoff, rho, penalty, mixing):
        edge = int(WIDTH_CAP * cutoff / 2)
        start, stop = edge // 25, (edge + 1000) // 25
        # Both margins increase with N and the even sector has the larger g.
        if not _pole_ok(25 * stop) or not self.check(cutoff, rho, 25 * stop, 0, penalty, mixing)[0]:
            return None
        while start < stop:
            middle = (start + stop) // 2
            if _pole_ok(25 * middle) and self.check(cutoff, rho, 25 * middle, 0, penalty, mixing)[0]:
                stop = middle
            else:
                start = middle + 1
        return 25 * start

    def mode(self, cutoff, mode, info, penalty, mixing):
        actual = self.first(cutoff, Fraction(5, 4), penalty, mixing)
        _require(actual == info["N"], f"different first grid dimension in {mode}")
        for row in info["scan"]:
            modes = row["N"]
            _require(row["pole_ok"] is _pole_ok(modes), f"pole disagreement in {mode} N={modes}")
            _require([sector["parity"] for sector in row["sectors"]] == [0, 1], "invalid gate parities")
            for sector in row["sectors"]:
                passes, g, b2 = self.check(cutoff, Fraction(5, 4), modes, sector["parity"], penalty, mixing)
                _require(sector["positive"] is passes, f"gate sign disagreement in {mode} N={modes}")
                self.sector_tests += 1
                for label, exact in (("g", g), ("b2", b2)):
                    if sector[label] == 0 and exact > 0:
                        self.display_underflow[label] += 1
        return actual


def _penalties(c, d, lam, k):
    s = lam * d / ((1 - lam) * c)
    f1 = c * s / (lam * (1 + s))
    eta = k * c / (k - 1)
    return s, f1, eta


def run(attachment=ATTACHMENT):
    data, source_sha = _read_json(attachment, ATTACHMENT_SHA256)
    archived, opt_sha = _read_json(CSV_PATH.with_name("opt38.json"), OPT_SHA256)
    perron_weights, weights_sha = _read_json(INPUTS / "perron_weights.json")
    by_profile = {(item["beta"], item["T"]): item for item in archived}
    rows = tuple(read_rows(CSV_PATH))
    a1 = data["A1"]
    expected_profiles = [(str(beta), cutoff) for beta in BETA_GRID for cutoff in CUTOFFS]
    _require([(p["beta"], p["T"]) for p in a1["uniform_25"]] == expected_profiles, "wrong profile grid")
    checked_perron = {}
    derived = {}
    for profile in a1["uniform_25"]:
        beta, cutoff = profile["beta"], profile["T"]
        stored = by_profile[beta, cutoff]
        bands = components(rows, Fraction(beta), cutoff)
        _require(len(bands) == profile["ncomp"] == stored["components"], f"wrong component count at {(beta, cutoff)}")
        left, right, mass, height = np.array(bands, dtype=np.int64).T
        kernel = _kernel(left, right)
        c = schur_bound(kernel, mass, stored["weights"])
        uniform = schur_bound(kernel, mass, [1] * len(mass))
        d = Fraction(int(max(height)), 10**6)
        _require(c == Fraction(profile["c_archived"]) == Fraction(stored["c"]), "archived Schur bound differs")
        _require(uniform == Fraction(profile["c_uniform"]) and d == Fraction(stored["d"]), "uniform bound or d differs")
        _require(float(uniform / c) == profile["ratio"], "archived displayed ratio differs")
        if beta == "2" and cutoff in (4500, 6400):
            old = a1[f"perron_T{cutoff}"]
            _require(old["ncomp"] == len(bands) and old["T"] == cutoff, "wrong Perron profile")
            _require(Fraction(old["c_archived"]) == c and Fraction(old["c_uniform"]) == uniform, "Perron profile header differs")
            results = {}
            for scale, info in old["perron"].items():
                vector = perron_weights[str(cutoff)][scale]
                bound = schur_bound(kernel, mass, vector)
                _require(bound == Fraction(info["c"]) and min(vector) == int(scale) and max(vector) == info["wmax"], "Perron replay differs")
                results[scale] = {"c_exact": str(bound), "maximum_weight": max(vector)}
            checked_perron[str(cutoff)] = results
            derived[cutoff] = c, d, min(Fraction(record["c_exact"]) for record in results.values())

    directed_gates = []
    with ctx.workprec(192):
        symbol = _DirectDigammaSymbol()
        _require(arb.pi() > _ball("333/106"), "the rational pi bound failed")
        low_symbol = {}
        for cutoff in derived:
            min_value = (arb(1) / 4).digamma() - symbol.log_pi - symbol.mass - 2
            max_value = acb(arb(1) / 4, arb(cutoff) / 2).digamma().real - symbol.log_pi + symbol.mass - 2
            absolute = max((-min_value).upper(), max_value.upper())
            _require(absolute < 50, "the absolute low symbol bound failed")
            low_symbol[str(cutoff)] = dyadic(absolute)
        replay = _GateReplay()
        gates = [(data["A4"], derived[4500][0])]
        for cutoff, (c, d, best) in sorted(derived.items()):
            packet = a1[f"gate_T{cutoff}"]
            _require(Fraction(packet["c_used"]) == best and Fraction(packet["d"]) == d, "selected A1 bounds differ")
            gates.extend((packet[label]["gate"], bound) for label, bound in (("archived", c), ("best", best)))
        for gate, c in gates:
            cutoff, k, lam = gate["T"], gate["k"], Fraction(gate["lam"])
            d = derived[cutoff][1]
            _require(Fraction(gate["c"]) == c and Fraction(gate["d"]) == d, "Schur parameters differ")
            s, f1, eta = _penalties(c, d, lam, k)
            _require(tuple(Fraction(gate[field]) for field in ("s", "f1", "eta")) == (s, f1, eta), "harmonic or localized penalty differs")
            localized = f"localized_M600_k{k}"
            modes = gate["modes"]
            directed_gates.append({
                "tag": gate["tag"], "cutoff": cutoff,
                "harmonic_first_dimension": replay.mode(cutoff, "harmonic", modes["harmonic"], f1, 50 + d / (1 - lam)),
                "localized_first_dimension": replay.mode(cutoff, localized, modes[localized], eta, 50),
            })
        grid_results = 0
        for cutoff, (c, d, best) in derived.items():
            packet = a1["rho_lambda_scan"][f"T{cutoff}"]
            _require(tuple(Fraction(packet[field]) for field in ("c_arch", "d", "c_best")) == (c, d, best), "grid parameters differ")
            for row in packet["grid"]:
                rho, lam = map(Fraction, (row["rho"], row["lam"]))
                for field, bound in (("N_archived", c), ("N_best", best)):
                    _, f1, eta = _penalties(bound, d, lam, packet["k"])
                    expected = row[field]
                    _require(replay.first(cutoff, rho, f1, 50 + d / (1 - lam)) == expected["harmonic"], "harmonic grid differs")
                    _require(replay.first(cutoff, rho, eta, 50) == expected["localized"], "localized grid differs")
                    grid_results += 2
    return {
        "schema": "NUMERIA.P93V40.V38_A1_A4_REPLAY.1",
        "attachment_sha256": source_sha, "deficit_csv_sha256": CSV_SHA256,
        "archived_opt38_sha256": opt_sha, "perron_weights_sha256": weights_sha,
        "profiles_reproduced": len(expected_profiles), "perron_reproduced": checked_perron,
        "earlier_spectral_interval": early_spectrum(),
        "low_absolute_symbol_upper_dyadic": low_symbol,
        "directed_schur_checks": directed_gates,
        "individual_sector_signs_reproduced": replay.sector_tests,
        "rho_lambda_grid_dimensions_reproduced": grid_results,
        "zero_float_displays_with_strictly_positive_mathematical_value": replay.display_underflow,
        "scope": "Schur transfer conditions only. The complete matrix heads are not certified by these tables.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    io = parser.add_mutually_exclusive_group(required=True)
    io.add_argument("--out", type=Path)
    io.add_argument("--verify", type=Path)
    args = parser.parse_args()
    result = run()
    if args.verify is not None:
        if result != json.loads(args.verify.read_text(encoding="utf-8")):
            raise ValueError("the V.38 A1/A4 replay differs")
        print(f"Verified {args.verify}")
    else:
        args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
