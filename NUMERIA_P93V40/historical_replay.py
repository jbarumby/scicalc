"""Reproduce the V39 archive's historical 12-mode certificates in a scratch directory."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
from shutil import copytree
import subprocess
import sys
from tempfile import TemporaryDirectory


SOURCE = Path(__file__).resolve().parents[1] / "NUMERIA_P93V39"
MANIFEST = "SHA256SUMS.txt"


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def audit(source=SOURCE):
    manifest = source / MANIFEST
    expected = {}
    for line in manifest.read_text(encoding="ascii").splitlines():
        key, name = line.split(maxsplit=1)
        if Path(name).name != name or name in expected:
            raise ValueError("invalid or repeated manifest entry")
        if digest(source / name) != key:
            raise ValueError(f"V39 SHA-256 mismatch: {name}")
        expected[name] = key
    if len(expected) != 35:
        raise ValueError("the V39 manifest must describe all 35 historical files")

    def run(args, cwd):
        return subprocess.run([sys.executable, *args], cwd=cwd,
                              stdout=subprocess.PIPE, check=True).stdout

    with TemporaryDirectory(prefix="v39-certificate-replay-") as temp:
        scratch = Path(temp) / "archive"
        copytree(source, scratch)
        reference = expected["replay39.json"]
        for options in ([], ["-O"]):
            if sha256(run([*options, "replay39.py"], scratch)).hexdigest() != reference:
                raise ArithmeticError("the historical adaptive Schur ledger differed")
        for name, output in (("congruence39.py", "congruence39.json"),
                             ("schur_adaptive39.py", "schur_adaptive39.json")):
            run([name], scratch)
            if digest(scratch / output) != expected[output]:
                raise ArithmeticError(f"historical replay differs: {output}")

        outcomes = {}
        width = json.loads((source / "moments_adaptive_binet_12.json").read_bytes())["width"]
        for script, tag in (("moments39.py", "adaptive"),
                            ("moments39_binet.py", "adaptive_binet")):
            run([script, "--N", "12", "--adaptive", "--width", width], scratch)
            names = [f"matrices_{tag}_12_p{p}.json" for p in (0, 1)]
            proof = f"moments_{tag}_12.json"
            names.append(proof)
            hashes = {name: digest(scratch / name) for name in names}
            if any(key != expected[name] for name, key in hashes.items()):
                raise ArithmeticError(f"historical {tag} matrix replay changed")
            records = json.loads((scratch / proof).read_bytes())["results"]
            gates = [r for r in records if r["N"] == 12 and r["h"] == "1/" + str(10**60)]
            if len(gates) != 2:
                raise ArithmeticError("expected 12 modes in both parities")
            outcomes[tag] = {
                "bytes_match": True,
                "hashes": hashes,
                "parities": [{key: row[key] for key in ("p", "positive", "pivots")} for row in gates],
            }

    return {
        "schema": "NUMERIA.P93V40.HISTORICAL_V39_REPLAY.1",
        "manifest_sha256": digest(manifest),
        "manifest_entries": len(expected),
        "width_exact": width,
        "schur_and_congruence_replay_match": True,
        "schur_exact_replay_sha256": reference,
        "adaptive": outcomes,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    io = parser.add_mutually_exclusive_group(required=True)
    io.add_argument("--out", type=Path)
    io.add_argument("--verify", type=Path)
    args = parser.parse_args()
    result = audit()
    if args.verify is not None:
        if json.loads(args.verify.read_bytes()) != result:
            raise ValueError("the historical V39 audit differs")
        print(f"Verified {args.verify}")
    else:
        args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
