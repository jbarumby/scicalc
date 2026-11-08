#!/usr/bin/env python3
"""ES8-P101: calculate ONE direct Arb Weil matrix entry, fail closed.

Based on ES8-P99's `compute_entry`, with no P89/P95 matrix data as input.
W22 equals W11 EXACTLY by translation invariance for the chosen two
identical bumps, so only W11 and W12 are required.
"""
import argparse
import json
import sys
import time
from pathlib import Path

from flint import ctx

ROOT = Path(__file__).resolve().parent.parent / "es8_p99"
sys.path.insert(0, str(ROOT))
from certifica_matriz_2x2_P98 import compute_entry, make_integrator


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--entry", choices=["W11", "W12"], required=True)
    p.add_argument("--prec", type=int, default=128)
    p.add_argument("--tol", default="1e-8")
    p.add_argument("--eps-den", type=int, default=2048)
    p.add_argument("--eval-limit", type=int, default=60000)
    p.add_argument("--depth-limit", type=int, default=40)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()

    if a.prec < 128 or a.eps_den != 2048 or a.tol != "1e-8":
        p.error("ES8-P101 fixes bits>=128, epsilon=1/2048 and tol=1e-8")
    started = time.monotonic()
    ctx.prec = a.prec
    n, m = (1, 1) if a.entry == "W11" else (1, 2)
    print("START", a.entry, "PREC", a.prec, "TOL", a.tol, flush=True)
    try:
        calc = make_integrator(a.tol, a.eval_limit, a.depth_limit, a.eps_den)
        entry = compute_entry(n, m, calc)
        result = {
            "status": "FINITE_INTERVAL_ENTRY_CONDITIONAL_ON_FORMULA",
            "entry_name": a.entry,
            "normalization": "P92",
            "precision_bits": a.prec,
            "epsilon_denominator": a.eps_den,
            "tolerance": a.tol,
            "eval_limit": a.eval_limit,
            "depth_limit": a.depth_limit,
            "entry": entry,
            "elapsed_seconds": time.monotonic() - started,
            "review_warning": "Nested Arb integral and explicit formula need independent proof; no RH claim."
        }
    except Exception as ex:
        result = {
            "status": "INCONCLUSIVE",
            "entry_name": a.entry,
            "reason": f"{type(ex).__name__}: {ex}",
            "elapsed_seconds": time.monotonic() - started
        }
    a.output.write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print("FINISHED", a.entry, result["status"],
          "ELAPSED", result["elapsed_seconds"], flush=True)
    if result["status"] == "INCONCLUSIVE":
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
