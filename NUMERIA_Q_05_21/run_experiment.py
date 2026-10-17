#!/usr/bin/env python3
"""Executa benchmark de inversão modular e grava CSV/JSON."""
from __future__ import annotations
import argparse, csv, json
from pathlib import Path
from generators import is_prime, generate_minterm, generate_anf
from evaluator import verify_exact, resource_metrics, semantic_hash

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--primes", nargs="+", type=int, default=[7,11,13,17])
    ap.add_argument("--out", default="resultados_05_21")
    args = ap.parse_args(); out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    rows=[]; summary=[]
    for p in args.primes:
        if not is_prime(p): raise SystemExit(f"{p} não é primo")
        for gen in (generate_minterm, generate_anf):
            circuit=gen(p); verification=verify_exact(circuit,p)
            metrics=resource_metrics(circuit)
            row={"p":p,"architecture":circuit.name,**metrics,**verification,"semantic_hash":semantic_hash(p,circuit.n_input)}
            rows.append(row); summary.append({"p":p,"architecture":circuit.name,"passed":verification["passed"],"checks":verification.get("checks",0),**metrics})
            if not verification["passed"]: raise SystemExit(f"Falha de equivalência: p={p}, arquitetura={circuit.name}: {verification}")
    fields=sorted({k for row in rows for k in row if not isinstance(row.get(k), (dict,list))})
    with (out/"comparacao_arquiteturas.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
        for row in rows: w.writerow({k:row.get(k,"") for k in fields})
    (out/"resumo_05_21.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding="utf-8")
    print(f"OK: {len(rows)} circuitos avaliados; relatório em {out}")
if __name__=="__main__": main()
