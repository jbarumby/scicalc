#!/usr/bin/env python3
"""E8-K V18: cotejo independente, fail-closed, executado sobre arquivos do GitHub."""
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SHA = {
    "v17/RESULTADO_E8K_V17.json": "671daf0fb8b7304cd0722205db1b9b4dab3479f9cbb7019c4b2a66b777230d89",
    "v17/PARECER_TECNICO_V17.md": "58a45c21c29f02cf98b5a40c20686c87ae58c889354fc16de58e704a4f599276",
    "v17/AUDITAR_E8K_V17_NORMALIZACAO.py": "09053afc435e409f8254330896759b913254b38cd9c5a5c39088bb6b1880898a",
    "v18/AUDITAR_E8K_V18_TESTEMUNHAS.py": "0c78d6211ce145d99a808e58532b73c2c3efcdac41b38b3a6dbf3836833be3ca",
    "NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip": "4acf99a2625a262b670ee73ed104b970da0951ffd98250cdd254f387110b02a7",
}
GENERATED_SHA = "ee1ad12b3b04cefc995f627fd4cee9c73abb25a11180f5e8c13405ac2c402484"
ARTIFACT_MAP = {
    "v17_report_sha256": "v17/RESULTADO_E8K_V17.json",
    "v17_parecer_sha256": "v17/PARECER_TECNICO_V17.md",
    "v17_source_sha256": "v17/AUDITAR_E8K_V17_NORMALIZACAO.py",
}

def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()

def run(out, report=None, e7=True):
    cmd = [
        sys.executable, str(ROOT / "v18/AUDITAR_E8K_V18_TESTEMUNHAS.py"),
        "--v17-report", str(report or ROOT / "v17/RESULTADO_E8K_V17.json"),
        "--v17-parecer", str(ROOT / "v17/PARECER_TECNICO_V17.md"),
        "--v17-source", str(ROOT / "v17/AUDITAR_E8K_V17_NORMALIZACAO.py"),
        *(["--e7", str(ROOT / "NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip")] if e7 else []),
        "--out", str(out), "--strict",
    ]
    return subprocess.run(cmd, capture_output=True, text=True, timeout=90)

def assert_evidence(report):
    assert report["schema"] == "NUMERIA.E8K.V18.WITNESS_INDEPENDENT_REPLAY.1"
    assert report["status"] == "V17_WITNESSES_REPLAYED_MUTATIONS_DETECTED"
    gates = report["gates"]
    assert len(gates) == 14 and all(v["status"] == "PASS" for v in gates.values())
    assert len(report["mutations"]) == 13 and all(m["detectada"] is True for m in report["mutations"])
    assert len(report["decorative_cases_v17_historical_not_evidence"]) == 1
    assert report["decorative_cases_v17_historical_not_evidence"][0]["k"] == 5
    samples = report["resonance_independent_v18"]["amostras"]
    assert len(samples) == 10 and all(x["concorda_float"] is True for x in samples)
    assert all(x["modulo_quadrado_analitico_exato"] == "9/20" for x in samples)
    assert all(report["v17_artifact_sha256"][key] == SHA[value] for key, value in ARTIFACT_MAP.items())
    assert report["not_certified"], "Limites cientificos omitidos"

def main():
    try:
        for path, expected in SHA.items():
            full = ROOT / path
            assert full.is_file() and digest(full) == expected, "Falha SHA-256: " + path
        original = json.loads((ROOT / "v17/RESULTADO_E8K_V17.json").read_text())
        assert original["schema"] == "NUMERIA.E8K.V17.EXPLICIT_FORMULA_SPECIALIZATION.1"
        snapshot = json.loads((ROOT / "v18/RESULTADO_E8K_V18.json").read_text())
        assert_evidence(snapshot)
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            first = directory / "v18-first.json"
            second = directory / "v18-second.json"
            one = run(first)
            assert one.returncode == 0, one.stderr + one.stdout
            report = json.loads(first.read_text())
            assert_evidence(report)
            assert digest(first) == GENERATED_SHA, "Digest do replay nao reproduzido"
            two = run(second)
            assert two.returncode == 0 and first.read_bytes() == second.read_bytes(), "JSON nao deterministico"
            missing = directory / "sem-e7.json"
            assert run(missing, e7=False).returncode == 2, "E7 ausente nao reprovou"
            assert json.loads(missing.read_text())["status"] == "AUDIT_FAILED"
            for tag, mutation in (
                ("fator_primo", lambda j: j["fourier_convention"].__setitem__("prime_term", "-1*(log 2)/sqrt 2*C_ij(log 2)")),
                ("testemunha_suporte", lambda j: j["support_witness"].__setitem__("strict_inequality", False)),
                ("hash_e7", lambda j: j["source_hashes"].__setitem__("e7_zip_sha256", "0"*64)),
            ):
                altered = json.loads(json.dumps(original))
                mutation(altered)
                input_path = directory / (tag + "-input.json")
                output_path = directory / (tag + "-result.json")
                input_path.write_text(json.dumps(altered))
                assert run(output_path, report=input_path).returncode == 2, tag + " aceito incorretamente"
                assert json.loads(output_path.read_text())["status"] == "AUDIT_FAILED"
        print("PASS: 5 hashes, 14/14 gates, 13/13 mutacoes, 10 ressonancias, 2 replays identicos, 4 testes negativos")
        print("NAO CERTIFICADO: cotas de resto, cauda, validade global Weil ou RH")
        return 0
    except Exception as exc:
        print("FAIL:", exc, file=sys.stderr)
        return 2

if __name__ == "__main__":
    sys.exit(main())
