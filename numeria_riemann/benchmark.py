"""Reproducible benchmark runner for the isolated NUMERIA-RIEMANN package.

Usage:
  python -m numeria_riemann.benchmark --out numeria_riemann/results/run1 --zeros 20

The output includes CSV tables, environment metadata, and SHA-256 hashes.
It does not claim certified zeros or prove the Riemann Hypothesis.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import platform
from pathlib import Path
from datetime import datetime, timezone

import mpmath as mp

from numeria_riemann.backend_compare import zero_residual_report
from numeria_riemann.explicit_formula import convergence_report


def _write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run(out: Path, zero_count: int, dps: int, x: float) -> dict:
    if zero_count < 1:
        raise ValueError("zero_count must be >= 1")
    if dps < 15:
        raise ValueError("dps must be >= 15")
    out.mkdir(parents=True, exist_ok=True)

    residuals = zero_residual_report(zero_count, dps=dps)
    counts = sorted(set(n for n in (1, 5, 10, 20, 50, zero_count) if n <= zero_count))
    convergence = convergence_report(x, counts, dps=dps)
    residual_path = out / "zero_residuals.csv"
    convergence_path = out / "explicit_formula_convergence.csv"
    _write_csv(residual_path, residuals)
    _write_csv(convergence_path, convergence)

    try:
        flint_version = importlib.metadata.version("python-flint")
    except importlib.metadata.PackageNotFoundError:
        flint_version = None
    metadata = {
        "project": "NUMERIA-RIEMANN",
        "status": "exploratory_numerical_benchmark_not_certified",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "mpmath": importlib.metadata.version("mpmath"),
        "python_flint": flint_version,
        "zero_count": zero_count,
        "decimal_precision": dps,
        "explicit_formula_x": x,
        "outputs": {},
        "limitations": [
            "mpmath zero residuals are numerical diagnostics, not zero certificates",
            "finite explicit-formula truncations may oscillate",
            "this benchmark does not prove the Riemann Hypothesis",
        ],
    }
    for path in (residual_path, convergence_path):
        metadata["outputs"][path.name] = {"sha256": _sha256(path), "bytes": path.stat().st_size}
    metadata_path = out / "metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    manifest = {
        path.name: _sha256(path)
        for path in sorted(out.iterdir())
        if path.is_file() and path.name != "SHA256SUMS.json"
    }
    (out / "SHA256SUMS.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("numeria_riemann/results/latest"))
    parser.add_argument("--zeros", type=int, default=20)
    parser.add_argument("--dps", type=int, default=40)
    parser.add_argument("--x", type=float, default=10.5)
    args = parser.parse_args()
    print(json.dumps(run(args.out, args.zeros, args.dps, args.x), indent=2))


if __name__ == "__main__":
    main()
