#!/usr/bin/env python3
"""P256: source-pinned map of the RH -> EveryCharacter proof route.

This checks exact source bytes and selected declarations/calls only. It does
NOT invoke Lean, replay the proof, validate imported theorem bodies, or prove RH.
"""
from __future__ import annotations
import argparse
from hashlib import sha256
import json
from pathlib import Path

SOURCE_COMMIT = "e50d09817ad81469f737e8f637561e6ccf86c670"
SOURCE_REPO = "https://github.com/Deicyde/riemann-hypothesis"
SOURCE_HASHES = {
    "EveryCharacter.lean": "31f6e49936226cfc3d71b15685347bd2874db0aef384e5cf4fea7d09ba68dc71",
    "MainConclusion.lean": "a5f71e107087e37b2c598709795b0f56c0a6741cf1e8ab99ee1f97e2fce0fe34",
    "VwinOfRow.lean": "d2df45a5cf337190ae6f885e8e43056ae4e804ae8da250a2ab0abd64b1ca7080",
    "EveryCharacterOfLiteratureOfVWin.lean": "246b74eeb03dae250c45b46c155b3590cb53c2bac86f5e35f4b74338db487068",
    "CritLine.lean": "5b1552ac13442150f6900f37b9b44eff1491fad1bebe2ec07a959eecb9bcfe02",
    "RhOfRight.lean": "7e2c0b519ee420192cad8fb281c2d2bfebae2ba553342fd579c81c2aa0bd3dc5",
    "MainTheoremsAll.lean": "a209a6a2f6ea1cc4abe865996b909bfeef391adcb89d3a948858874cb5956fec",
}
REQUIRED = {
    "EveryCharacter.lean": (
        "def UpperTheory.EveryCharacter : Prop :=",
        "DirichletCharacter.LFunction χ s ≠ 0",
    ),
    "MainConclusion.lean": (
        "theorem UpperTheory.FinalTheorems.everyCharacter_of_row",
        "UpperTheory.Gal.Int.galoisLocalGlobalStrong_recParameter",
        "theorem UpperTheory.FinalTheorems.everyCharacter : UpperTheory.EveryCharacter",
    ),
    "VwinOfRow.lean": (
        "theorem UpperTheory.One.vwin_of_row",
        "UpperTheory.GalStar.Row0b",
        "UpperTheory.Top.VWin",
        "automorphic_window_of_inputs",
    ),
    "EveryCharacterOfLiteratureOfVWin.lean": (
        "every_character_of_literature_of_vWin",
        "every_character_of_theta_integrals_local_matching",
        "UpperTheory.Top.VWin",
    ),
    "CritLine.lean": (
        "theorem Formalism.VW.rh_of_every",
        "Formalism.VW.rh_of_right",
        "DirichletCharacter ℂ 1",
    ),
    "RhOfRight.lean": (
        "theorem Formalism.VW.rh_of_right",
        "riemannZeta_one_sub",
        "RiemannHypothesis",
    ),
    "MainTheoremsAll.lean": (
        "import Supplements.Kernels.MainTheoremsCritical",
        "import Development.AutomorphicArithmetic.AutomorphicVarieties.SmallClassNumbers.CorGKOf",
    ),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def audit(source_dir: Path, p255_dir: Path) -> dict:
    from NUMERIA_ES8P255.audit_statements import check_source

    source_text: dict[str, str] = {}
    verified: dict[str, str] = {}
    for name, expected in SOURCE_HASHES.items():
        raw = (source_dir / name).read_bytes()
        actual = sha256(raw).hexdigest()
        require(actual == expected, f"{name}: source SHA-256 mismatch")
        source_text[name] = raw.decode("utf-8")
        verified[name] = actual

    # Reuse P255's exact source pins for the published theorem names and entry point.
    p255_text: dict[str, str] = {}
    for name in ("Main.lean", "MainTheorems.lean"):
        p255_text[name] = check_source(name, (p255_dir / name).read_bytes())

    checks = [
        ("top_level_RH", p255_text["Main.lean"], "theorem riemannHypothesis : RiemannHypothesis := rhOfAll"),
        ("allChi_entry", p255_text["MainTheorems.lean"], "theorem allChi : UpperTheory.EveryCharacter := UpperTheory.FinalTheorems.everyCharacter"),
        ("RH_reduction", p255_text["MainTheorems.lean"], "theorem rhOfAll : RiemannHypothesis := Formalism.VW.rh_of_every allChi"),
    ]
    for name, text in source_text.items():
        checks.extend((f"{name}:{i}", text, snippet) for i, snippet in enumerate(REQUIRED[name]))
    for label, text, snippet in checks:
        require(snippet in text, f"expected source link missing: {label}")

    return {
        "schema": "NUMERIA.ES8P256.RH_DEPENDENCY_MAP.1",
        "source_repository": SOURCE_REPO,
        "source_commit": SOURCE_COMMIT,
        "source_sha256_verified": verified,
        "p255_headline_sources_verified_by_pinned_checker": ["Main.lean", "MainTheorems.lean"],
        "checked_links": [label for label, _, _ in checks],
        "formal_input": "EveryCharacter: nonvanishing of every Dirichlet L-function in Re(s)>1/2, away from the trivial pole",
        "formal_reduction": "allChi -> rhOfAll -> rh_of_every -> rh_of_right; the trivial character modulo 1 gives zeta nonvanishing to the right, then the functional equation reflects zeros",
        "concentrated_upstream_input": "galoisLocalGlobalStrong_recParameter supplies Row0b; vwin_of_row yields VWin; the theta-integral/local-matching theorem yields EveryCharacter",
        "lean_kernel_invoked": False,
        "full_solution_replayed": False,
        "imported_theorem_bodies_independently_verified": False,
        "interpretation": "source-pinned dependency map only; not an independent proof or kernel replay",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=Path(__file__).with_name("lean_sources"))
    parser.add_argument("--p255-dir", type=Path, default=Path(__file__).parents[1] / "NUMERIA_ES8P255" / "fixtures")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.source_dir, args.p255_dir)
    args.out.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"P256 verified {len(result['checked_links'])} pinned source links; no Lean kernel invoked")


if __name__ == "__main__":
    main()
