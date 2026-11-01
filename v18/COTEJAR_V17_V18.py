#!/usr/bin/env python3
"""Cotejo fail-closed V17/V18: hashes, evidencias, replay e adulteracoes reais.
Nao certifica as hipoteses analiticas nem a Hipotese de Riemann.
"""
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SHA = {
    "v17/RESULTADO_E8K_V17.json": "3e57f1c79a0f07850db70e0513438b920ff5204781cc7d81ba16048eaa2023e3",
    "v17/PARECER_TECNICO_V17.md": "f1517c4db25fb0d9d918e852be02f0c037842d95e665a525793ecaa8732fc67e",
    "v17/AUDITAR_E8K_V17_NORMALIZACAO.py": "8d2fa7cdfc2281da9bde56dce665ce115124e6359165996d05fe91e50bcb13e4",
    "v18/RESULTADO_E8K_V18.json": "93efb06d5a2909b48c901b3e29943b8981c9958322cab2958e52bc8da6500050",
    "v18/AUDITAR_E8K_V18_TESTEMUNHAS.py": "0c78d6211ce145d99a808e58532b73c2c3efcdac41b38b3a6dbf3833be3ca",
    "NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip": "4acf99a2625a262b670ee73ed104b970da0951ffd98250cdd254f387110b02a7",
}
MAP = {
    "v17_report_sha256": "v17/RESULTADO_E8K_V17.json",
    "v17_parecer_sha256": "v17/PARECER_TECNICO_V17.md",
    "v17_source_sha256": "v17/AUDITAR_E8K_V17_NORMALIZACAO.py",
}

def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()

def validate():
    for rel, expected in SHA.items():
        p = ROOT / rel
        assert p.is_file(), "Ausente: " + rel
        assert digest(p) == expected, "Hash alterado: " + rel
    manifest = ROOT / "v18/SHA256SUMS_V18.txt"
    lines = manifest.read_text(encoding="utf-8").splitlines()
    expected_rel = set(SHA) - {"NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip"}
    seen = set()
    for line in lines:
        h, relative = line.split("  ", 1)
        p = (manifest.parent / relative).resolve()
        rel = p.relative_to(ROOT.resolve()).as_posix()
        assert rel in expected_rel and rel not in seen, "Linha/arquivo inesperado no manifesto: " + rel
        assert h == digest(p), "SHA256SUMS divergiu: " + rel
        seen.add(rel)
    assert seen == expected_rel, "Manifesto incompleto"
    v17 = json.loads((ROOT / "v17/RESULTADO_E8K_V17.json").read_text(encoding="utf-8"))
    v18 = json.loads((ROOT / "v18/RESULTADO_E8K_V18.json").read_text(encoding="utf-8"))
    assert v17["schema"] == "NUMERIA.E8K.V17.EXPLICIT_FORMULA_SPECIALIZATION.1"
    assert v18["schema"] == "NUMERIA.E8K.V18.WITNESS_INDEPENDENT_REPLAY.1"
    assert v18["status"] == "V17_WITNESSES_REPLAYED_MUTATIONS_DETECTED"
    assert len(v18["gates"]) == 14 and all(v["status"] == "PASS" for v in v18["gates"].values())
    assert len(v18["mutations"]) == 13 and all(m["detectada"] is True for m in v18["mutations"])
    assert all(v18["v17_artifact_sha256"][key] == SHA[rel] for key, rel in MAP.items())
    assert len(v18["decorative_cases_v17_historical_not_evidence"]) == 1
    assert v18["decorative_cases_v17_historical_not_evidence"][0]["k"] == 5
    rs = v18["resonance_independent_v18"]
    assert len(rs["amostras"]) == 10 and all(x["concorda_float"] is True for x in rs["amostras"])
    assert all(x["modulo_quadrado_analitico_exato"] == "9/20" for x in rs["amostras"])
    assert v18["not_certified"], "Certificacoes pendentes omitidas"

def run(out, v17_report=None, e7=True):
    return subprocess.run([
        sys.executable, str(ROOT / "v18/AUDITAR_E8K_V18_TESTEMUNHAS.py"),
        "--v17-report", str(v17_report or ROOT / "v17/RESULTADO_E8K_V17.json"),
        "--v17-parecer", str(ROOT / "v17/PARECER_TECNICO_V17.md"),
        "--v17-source", str(ROOT / "v17/AUDITAR_E8K_V17_NORMALIZACAO.py"),
        *(["--e7", str(ROOT / "NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip")] if e7 else []),
        "--out", str(out), "--strict"
    ], capture_output=True, text=True, timeout=90)

def main():
    try:
        validate()
        with tempfile.TemporaryDirectory() as temp:
            d = Path(temp)
            baseline = run(d / "base.json")
            assert baseline.returncode == 0, baseline.stderr + baseline.stdout
            assert (d / "base.json").read_bytes() == (ROOT / "v18/RESULTADO_E8K_V18.json").read_bytes(), "Cotejo byte a byte divergente"
            second = run(d / "second.json")
            assert second.returncode == 0
            assert (d / "base.json").read_bytes() == (d / "second.json").read_bytes(), "Nao deterministico"
            without_e7 = run(d / "missing-e7.json", e7=False)
            assert without_e7.returncode == 2, "Sem E7 deveria reprovar"
            assert json.loads((d / "missing-e7.json").read_text())["status"] == "AUDIT_FAILED"
            v17 = json.loads((ROOT / "v17/RESULTADO_E8K_V17.json").read_text())
            for label, change in (
                ("prime_factor", lambda x: x["fourier_convention"].__setitem__("prime_term", "-1*(log 2)/sqrt 2*C_ij(log 2)")),
                ("support_witness", lambda x: x["support_witness"].__setitem__("strict_inequality", False)),
                ("source_hash", lambda x: x["source_hashes"].__setitem__("e7_zip_sha256", "0"*64)),
            ):
                m = json.loads(json.dumps(v17))
                change(m)
                p = d / (label + ".json")
                p.write_text(json.dumps(m), encoding="utf-8")
                result = run(d / (label + "-out.json"), v17_report=p)
                assert result.returncode == 2, "Adulteracao nao detectada: " + label
                assert json.loads((d / (label + "-out.json")).read_text())["status"] == "AUDIT_FAILED"
        print("PASS: 6 hashes + manifesto, 14 gates, 13 mutacoes, 10 ressonancias, replay deterministico, 4 cenarios negativos")
        print("AVISO: testes de software; cotas de quadratura e cauda continuam condicionais.")
        return 0
    except Exception as exc:
        print("FAIL:", exc, file=sys.stderr)
        return 2

if __name__ == "__main__":
    sys.exit(main())
