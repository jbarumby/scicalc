#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""NUMERIA P78-S2 — Varredura heuristica de consumidores de B1 (downstream).

Percorre o repo (*.py e *.json) e lista, por arquivo, tudo que pode precisar
de re-certificacao apos a re-impressao S2 dos b1 das linhas 10-16:

  - referencias a b1 (row.b1, "b1" como chave, identificador B1);
  - ocorrencias dos valores ANTIGOS impressos do P77
    (17.21, 19.31, 21.77, 24.83, 29.01, 35.35, 46.06) — tabelas,
    certificados e caches que ainda carregam os B1 antigos;
  - ocorrencias dos valores NOVOS S2 (17.49, 19.66/19.67, 22.24/22.25, ...);
  - imports/referencias a modulos NUMERIA_P77* e NUMERIA_P78*.

HEURISTICO (nivel de texto): serve para planejar a re-certificacao
downstream e nao decide nada sozinho; a classificacao sugere o papel de
cada arquivo (TABLE = fonte da verdade; certificados = regenerar;
prototipos S1/S2 = hits esperados). A decisao final e do mantenedor.

Uso:
  python3 p78_s2_downstream_scan.py                  (repo auto-detectado)
  python3 p78_s2_downstream_scan.py --root CAMINHO   (raiz arbitraria)
  python3 p78_s2_downstream_scan.py --saida out.json
"""
import argparse
import json
import re
from pathlib import Path

SCHEMA = "NUMERIA.P78.S2.DOWNSTREAM.SCAN"
OLD_B1 = {"17.21": 10, "19.31": 11, "21.77": 12, "24.83": 13,
          "29.01": 14, "35.35": 15, "46.06": 16}
NEW_B1 = ["17.49", "19.66", "19.67", "22.24", "22.25",
          "25.48", "29.89", "36.59", "51.65"]
PROTOTYPE_HINTS = ("sensitivity", "resolution", "recert",
                   "downstream_scan", "table_patch", "make_s2_patch")


def find_repo():
    here = Path(__file__).resolve().parent
    for cand in [here] + list(here.parents):
        if (cand / "NUMERIA_P76_COEFFICIENTS").is_dir():
            return cand
    for cand in (Path("/home/z/my-project/scicalc-p78"), Path.cwd()):
        if cand is not None and (cand / "NUMERIA_P76_COEFFICIENTS").is_dir():
            return cand
    return None


def classify(rel):
    low = rel.lower()
    if any(h in low for h in PROTOTYPE_HINTS):
        return "prototipo_S1_S2 (hits esperados)"
    if "p76" in low and "table" in low:
        return "TABLE_fonte_da_verdade (atualizar b1 linhas 10-16)"
    if "certificate" in low or "cert" in low:
        return "certificado (regenerar apos re-certificacao)"
    if "p77" in low or "p78" in low:
        return "etapa_P77_P78 (revisar)"
    return "revisar"


def scan(root):
    pats = {
        "b1_attr": re.compile(r"\.b1\b"),
        "b1_key": re.compile(r"[\"']b1[\"']"),
        "B1_ident": re.compile(r"(?<![A-Za-z0-9_])B1(?![A-Za-z0-9_])"),
        "ref_P77": re.compile(r"NUMERIA_P77"),
        "ref_P78": re.compile(r"NUMERIA_P78"),
    }
    old_pats = {v: re.compile(r"(?<![\d.])" + re.escape(v) + r"(?![\d])")
                for v in OLD_B1}
    new_pats = {v: re.compile(r"(?<![\d.])" + re.escape(v) + r"(?![\d])")
                for v in NEW_B1}

    files_out = []
    paths = sorted(root.rglob("*.py")) + sorted(root.rglob("*.json"))
    for f in paths:
        try:
            text = f.read_text(errors="replace")
        except OSError:
            continue
        rel = str(f.relative_to(root))
        refs = {name: [] for name in pats}
        old_hits, new_hits = {}, {}
        for i, line in enumerate(text.splitlines(), 1):
            for name, pat in pats.items():
                if pat.search(line):
                    refs[name].append(i)
            for v, pat in old_pats.items():
                if pat.search(line):
                    old_hits.setdefault(str(OLD_B1[v]), []).append(i)
            for v, pat in new_pats.items():
                if pat.search(line):
                    new_hits.setdefault(v, []).append(i)
        if not (any(refs.values()) or old_hits or new_hits):
            continue
        rec = {
            "file": rel,
            "classification": classify(rel),
            "b1_references": {k: v for k, v in refs.items() if v},
        }
        if old_hits:
            rec["old_B1_values_still_present"] = old_hits
        if new_hits:
            rec["new_B1_values_present"] = new_hits
        files_out.append(rec)

    review = [r["file"] for r in files_out
              if not r["classification"].startswith("prototipo")]
    summary = {
        "root": str(root),
        "files_scanned": len(paths),
        "files_with_hits": len(files_out),
        "files_to_review": review,
        "table_updated_hint": (
            "a TABLE (fonte da verdade) deve receber os novos b1 das linhas "
            "10-16; certificados e caches com os valores antigos devem ser "
            "regenerados; os demais arquivos listados entram na re-certificacao "
            "downstream, no padrao do pipeline NUMERIA."),
        "verdict": ("HITs encontrados — revisar arquivos em files_to_review"
                    if files_out else "nenhum hit — nada a revisar"),
    }
    return {"schema": SCHEMA, "summary": summary, "files": files_out}


def main():
    ap = argparse.ArgumentParser(
        description="Varredura heuristica de consumidores de B1 (downstream)")
    ap.add_argument("--root", default=None, help="raiz do repo (default: auto-detectar)")
    ap.add_argument("--saida", default=None, help="arquivo de saida JSON")
    args = ap.parse_args()
    root = Path(args.root) if args.root else find_repo()
    if root is None:
        raise SystemExit("repo nao encontrado; use --root CAMINHO")
    doc = scan(root)
    print(json.dumps(doc["summary"], indent=1))
    for r in doc["files"]:
        extra = ""
        if r.get("old_B1_values_still_present"):
            extra += " ANTIGOS:" + ",".join(sorted(r["old_B1_values_still_present"]))
        if r.get("new_B1_values_present"):
            extra += " NOVOS:" + ",".join(sorted(r["new_B1_values_present"]))
        print(f"  {r['file']:<60} [{r['classification']}]{extra}")
    dest = (Path(args.saida) if args.saida
            else root.parent / "p78_s2_downstream_scan.json")
    dest.write_text(json.dumps(doc, indent=1) + "\n")
    print(f"[ok] varredura gravada em {dest}")


if __name__ == "__main__":
    main()