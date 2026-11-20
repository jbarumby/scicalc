#!/usr/bin/env python3
"""Informative Monte Carlo estimate for hddivsteps walk length; not a proof."""
import collections
import json
import random
import sys
import time

P = 2**256 - 2**32 - 977


def steps(x: int) -> int:
    delta, f, g, n = 1, P, x, 0
    while g:
        if g & 1 and delta > 0:
            delta, f, g = 2 - delta, g, (g - f) >> 1
        elif g & 1:
            delta, g = 2 + delta, (g + f) >> 1
        else:
            delta, g = 2 + delta, g >> 1
        n += 1
    return n


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 100_000
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 20261018
    rng = random.Random(seed)
    started = time.monotonic()
    counts = collections.Counter(steps(rng.randrange(1, P)) for _ in range(n))
    thresholds = (500, 510, 520, 525, 530, 533, 540, 550, 590)
    print(json.dumps({
        "n": n,
        "seed": seed,
        "min_steps": min(counts),
        "max_steps": max(counts),
        "mean_steps": sum(k * v for k, v in counts.items()) / n,
        "count_steps_greater_than": {str(t): sum(v for k, v in counts.items() if k > t) for t in thresholds},
        "elapsed_seconds": round(time.monotonic() - started, 3),
        "interpretation": "sampling only; not a universal bound or circuit failure diagnosis"
    }, indent=2))
