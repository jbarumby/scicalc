"""Authenticate Alpöge's original source and inspect the MIT challenge export."""

import argparse
from collections import Counter
from hashlib import file_digest, sha256
import json
from pathlib import Path, PurePosixPath
import re
import tarfile


ARCHIVE_SHA256 = "9631b4f020199e580ba86178cb56c6e9b09f690c2369a80dc4334aa93f49c7fc"
CHALLENGE_SHA256 = "bd484ac9886a96f26b1d9a7a59ebd45b7a70d102c015d4f31f9302bbdd6ab2b9"
KEY_FILES = {
    "lean-toolchain",
    "lakefile.toml",
    "lake-manifest.json",
    "src/Main.lean",
    "src/Statements.lean",
    "src/Supplements/Kernels/MainTheorems.lean",
    "src/Supplements/Kernels/MainTheoremsAll.lean",
    "src/Formalism/GammaFactorNonvanishing/RhOfRight.lean",
    "src/Formalism/GammaFactorNonvanishing/CritLine.lean",
    "src/Development/AutomorphicArithmetic/CharacterBoxes/EveryCharacter.lean",
}
SOURCE_PATTERNS = {
    "placeholder_terms": re.compile(rb"\b(?:sorry|admit|sorryAx)\b"),
    "axiom_commands": re.compile(rb"^\s*(?:(?:private|protected|public)\s+)?axiom\s+", re.M),
    "native_decide_terms": re.compile(rb"\b(?:native_decide|ofReduceBool)\b"),
    "addDecl_terms": re.compile(rb"\baddDecl\b"),
}
STATEMENT_TYPES = {
    "riemannHypothesis": "RiemannHypothesis",
    "riemannZeta_eq_zero_re": "RiemannZetaZerosOnLine",
    "LFunction_ne_zero_of_half_lt_re": "DirichletLNonvanishingRightOfHalf",
    "LFunction_eq_zero_re": "DirichletLZerosOnLine",
    "LFunction_eq_zero_re_of_isPrimitive": "PrimitiveDirichletLZerosOnLine",
}


def authenticate(path, expected):
    with path.open("rb") as stream:
        digest = file_digest(stream, "sha256").hexdigest()
    if digest != expected:
        raise ValueError(f"SHA-256 mismatch for {path}: {digest}")
    return {"sha256": digest, "bytes": path.stat().st_size}


def audit_source(path):
    archive = authenticate(path, ARCHIVE_SHA256)
    sources = Counter()
    matched = Counter()
    examples = {}
    selected = {}
    line_total = 0
    byte_total = 0
    with tarfile.open(path, mode="r|gz") as stream:
        for member in stream:
            name = PurePosixPath(member.name)
            if not name.parts or name.parts[0] != "lean" or ".." in name.parts:
                raise ValueError(f"unexpected archive path: {member.name}")
            if member.isdir():
                continue
            if not member.isfile():
                raise ValueError(f"unexpected archive member type: {member.name}")
            relative = name.relative_to("lean")
            is_lean = len(relative.parts) > 1 and relative.parts[0] == "src" and relative.suffix == ".lean"
            if not is_lean and str(relative) not in KEY_FILES:
                continue
            with stream.extractfile(member) as source:
                raw = source.read()
            if str(relative) in KEY_FILES:
                selected[str(relative)] = sha256(raw).hexdigest()
            if not is_lean:
                continue
            sources[relative.parts[1]] += 1
            line_total += raw.count(b"\n") + int(bool(raw) and not raw.endswith(b"\n"))
            byte_total += len(raw)
            for label, pattern in SOURCE_PATTERNS.items():
                if pattern.search(raw):
                    matched[label] += 1
                    examples.setdefault(label, str(relative))
    if selected.keys() != KEY_FILES:
        raise ValueError("missing key files")
    return {
        **archive,
        "lean_source_files": sum(sources.values()),
        "lean_source_lines": line_total,
        "lean_source_bytes": byte_total,
        "files_by_top_level": dict(sorted(sources.items())),
        "text_matching_file_counts": {label: matched[label] for label in SOURCE_PATTERNS},
        "first_text_matches": examples,
        "key_file_sha256": dict(sorted(selected.items())),
    }


def audit_challenge(path):
    """Confirm the exact constant types in the authenticated challenge data."""
    audit = authenticate(path, CHALLENGE_SHA256)
    targets = set(STATEMENT_TYPES) | set(STATEMENT_TYPES.values()) | {"sorryAx"}
    names = {}
    decls = {}
    name_pattern = re.compile(
        rb'"str":"(' + b"|".join(name.encode("ascii") for name in sorted(targets)) + rb')"'
    )
    with path.open("rb") as stream:
        for line in stream:
            if name_pattern.search(line):
                record = json.loads(line)
                entry = record.get("str", {})
                if entry.get("pre") == 0:
                    names[entry["str"]] = record["in"]
            elif line.startswith((b'{"thm":', b'{"def":', b'{"axiom":')):
                record = json.loads(line)
                kind, decl = next(iter(record.items()))
                if decl["name"] in names.values():
                    decls[decl["name"]] = (kind, decl["type"])

    type_symbols = {}
    type_ids = {decls[names[theorem]][1] for theorem in STATEMENT_TYPES}
    with path.open("rb") as stream:
        for line in stream:
            if not line.startswith(b'{"const":'):
                continue
            expr = json.loads(line)
            if expr["ie"] in type_ids:
                type_symbols[expr["ie"]] = expr["const"]

    for theorem, proposition in STATEMENT_TYPES.items():
        kind, type_id = decls[names[theorem]]
        if kind != "thm" or type_symbols.get(type_id) != {"name": names[proposition], "us": []}:
            raise ValueError(f"unexpected type for {theorem}")
        if decls[names[proposition]][0] != "def":
            raise ValueError(f"expected a definition for {proposition}")

    if decls[names["sorryAx"]][0] != "axiom":
        raise ValueError("expected the challenge placeholder axiom")
    return {
        **audit,
        "statement_types": STATEMENT_TYPES,
        "challenge_declares_sorryAx": True,
        "scope": "constant/type inspection only; this is the challenge, not the solution",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--challenge", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = {"source": audit_source(args.archive)}
    if args.challenge is not None:
        report["challenge"] = audit_challenge(args.challenge)
    args.out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
