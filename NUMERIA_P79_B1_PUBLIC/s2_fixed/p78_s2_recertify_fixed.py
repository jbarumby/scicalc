#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""NUMERIA P78-S2 (prototipo) — Re-certificacao dos B1 re-impressos com b4 publico.

Etapa de adocao do kit S2. Valida os novos b1 das linhas do certificado P78
(10-16) contra B1_upper recalculado em aritmetica de intervalos (K1/K2
rigorosos + algebra b4-dependente), com gate certificado do P77:

  (printed - B1_upper(b4 publico)).lower() > 0

MODOS
  PRE-ADOCAO  (--proposed PATCH.json): certifica os valores propostos do patch
              SEM editar a TABLE (variante minimal ou robust). Deve passar
              ANTES de qualquer edicao da TABLE.
  POS-ADOCAO  (sem --proposed): certifica o que esta na TABLE. Se a TABLE
              ainda tiver os valores antigos do P77, o gate DEVE falhar —
              isso confirma que a re-impressao e obrigatoria.

Varredura dyadic: --dyadic 192,256,384 repete todo o calculo nas precisoes
indicadas (K1/K2 recompute por precisao) e o gate tem de passar em todas.
Cada linha ainda passa pela checagem do modelo exato (e verificado contra
e = (2*s2-1)/(2-s2), s2 = (sigma_start+1)/32) e pela comparacao com a
expectativa float do patch (aviso se divergir > 1e-6 relativo).

Saida: p78_s2_recert_certificate.json em repo.parent, com margens por
precisao, SHA256 das entradas (table.py, audit_b23.py, certificate.json,
patch) e veredito. Prototipo single-run: o certificado OFICIAL continua
sendo o do pipeline NUMERIA (dyadic completo + CI + SHA256); este script
produz numeros/gates no mesmo padrao para integrar a esse pipeline.
Se validate_table() tiver checagens internas sobre b1, atualize-as junto
com os valores.

Uso:
  python3 p78_s2_recertify.py --proposed p78_s2_table_patch.json
  python3 p78_s2_recertify.py --proposed p78_s2_table_patch.json --variant minimal
  python3 p78_s2_recertify.py --proposed p78_s2_table_patch.json --dyadic 192,256,384
  python3 p78_s2_recertify.py                      (pos-adocao, TABLE)
"""
import argparse
import hashlib
import json
import math
import platform
import sys
from fractions import Fraction as Fr
from pathlib import Path

SCHEMA = "NUMERIA.P78.S2.RECERT.PROTOTYPE"
B54 = Fr(27, 50)          # b4 interno de referencia do P77 (0.54 exato)


def find_repo():
    here = Path(__file__).resolve().parent
    for cand in [here] + list(here.parents):
        if (cand / "NUMERIA_P76_COEFFICIENTS").is_dir():
            return cand
    for cand in (Path("/home/z/my-project/scicalc-p78"), Path.cwd()):
        if cand is not None and (cand / "NUMERIA_P76_COEFFICIENTS").is_dir():
            return cand
    return None


def e_closed(sigma_start):
    s2 = (sigma_start + 1) / 32.0
    return (2.0 * s2 - 1.0) / (2.0 - s2)


def sha256_file(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(
        description="Re-certificacao S2 dos B1 re-impressos (prototipo)")
    ap.add_argument("--repo", default=None, help="raiz do repositório NUMERIA")
    ap.add_argument("--out", default=None, help="arquivo JSON de saída")
    ap.add_argument("--proposed", metavar="PATCH_JSON", default=None,
                    help="patch JSON: certifica os valores propostos SEM editar a TABLE")
    ap.add_argument("--variant", choices=("minimal", "robust"), default=None,
                    help="variante do patch (default: default_variant do patch)")
    ap.add_argument("--dyadic", metavar="PRECS", default="256",
                    help="precisoes da varredura dyadic, ex.: 192,256,384 (default: 256)")
    args = ap.parse_args()

    repo = Path(args.repo).resolve() if args.repo else find_repo()
    if repo is None or not (repo / "NUMERIA_P76_COEFFICIENTS").is_dir():
        raise SystemExit("repo NUMERIA nao encontrado (NUMERIA_P76_COEFFICIENTS "
                         "ausente); rode dentro do repo")
    sys.path.insert(0, str(repo))
    try:
        import flint as _flint
        from flint import acb, arb, ctx
    except ImportError as exc:
        raise SystemExit(f"python-flint ausente: {exc}")
    from NUMERIA_P76_COEFFICIENTS.audit_b23 import constants, q, dyadic
    from NUMERIA_P76_COEFFICIENTS.table import TABLE, validate_table

    T0 = 3 * 10 ** 12

    # ------------------------------------------------------------------
    # Maquinaria de intervalos (mesmas definicoes do S1/p78_b1_sensitivity)
    # ------------------------------------------------------------------
    def rat(n, d):
        return arb(n) / arb(d)

    def pw(x, e):
        if not x > 0:
            raise ArithmeticError("non-positive real power base")
        return (x.log() * e).exp()

    def row_parts(row, cs):
        d, kap, delta = q(row.d), q(row.kappa), q(row.delta)
        L, pi, m0 = cs["logh"], cs["pi"], cs["m0"]
        A0 = (arb(10).log() * q(row.log10_a0)).exp()
        La = A0.log()
        half = (A0 / 2).log()
        a1 = A0 / (2 * pi * pi) * half ** 4 + q("34.842") * A0 * half ** 3 + arb("1.3897e16")
        r = 1 - arb.const_log2() / La
        a2 = r ** 4 / (2 * pi * pi) + q("34.842") / La * r ** 3 + arb("9.5079e15") / (A0 * La ** 4)
        c = 3 / (4 * kap)
        return d, kap, delta, L, pi, m0, A0, a1, a2, c

    def k1_k2(row, cs, parts):
        d, kap, delta, L, pi, m0, A0, a1, a2, c = parts

        def f(x, analytic):
            v = x ** 6
            A = acb(1) + (acb(1) + acb(c) * v).log(analytic=analytic) / (2 * acb(L))
            B = acb(cs["c1"] * c.sqrt()) * x ** 3 + acb(cs["c2"])
            return 6 * v * (-v).exp() * A.pow(acb(4) / 3, analytic=analytic) \
                * B.pow(acb(2) / 3, analytic=analytic)

        I = acb.integral(f, 0, 2, rel_tol=arb('1e-25'), abs_tol=arb('1e-25'),
                         eval_limit=30000, depth_limit=30)
        if not I.is_finite() or not I.imag.contains(0) or not I.real.lower() > 0:
            raise ArithmeticError("quadrature failed or inconclusive")
        U = arb(64)
        D = (1 + c / (2 * L)) ** 2 * (1 + cs["c1"] * c.sqrt() + cs["c2"])
        tail = (D * (-U).exp() * (U ** 3 + 3 * U ** 2 + 6 * U + 6)).upper()
        lower, upper = I.real.lower(), (I.real.upper() + tail).upper()
        F = pw(2 * a2 * c.sqrt(), rat(1, 3))
        correction = (pw(2 * a1, rat(1, 3)) * pw(cs["c2"] + A0 * cs["c1"] / T0, rat(2, 3))
                      / (pw(arb(T0), rat(1, 3)) * pw(L, rat(4, 3))))
        k1lo = (F * lower + correction).lower()
        k1 = (F * upper + correction).upper()
        if not k1lo > 0:
            raise ArithmeticError("K1 lower bound inconclusive")
        two13 = pw(arb(2), rat(1, 3))
        k2 = two13 * c.sqrt() * pi.sqrt() / 2 + two13 * (4 * kap / 3) * pw(A0 / arb(T0), arb(3))
        return k1, k2

    def b1_upper(row, cs, parts, k1, k2, b4):
        d, kap, delta, L, pi, m0, A0, a1, a2, c = parts
        s1, s2 = rat(row.sigma_start, 32), rat(row.sigma_start + 1, 32)
        c3 = pi * m0 * b4 / (2 * delta) * (2 * delta * arb.const_euler() / L).exp()
        c4 = b4 / (5 * delta * delta.exp()) + q(".605") * (-2 * delta).exp() / L
        k3 = c3 + c4 * (pi.sqrt() / (2 * (2 * kap).sqrt()) + pi * m0 / T0)
        l0 = (kap * (1 + (1 / (4 * arb(T0)) * (1 + 2 * delta / L)) ** 2)).exp() / (1 - arb(1) / T0)
        l1 = max(arb(1), pw(k1, 3 * delta / L)) * max(arb(1), pw(k3, -3 * delta / L)) \
            * (3 * delta + 4 * delta * L.log() / L).exp()
        l2 = max(arb(1), pw(k2, 3 * delta / L)) * max(arb(1), pw(k3, -3 * delta / L)) \
            * (3 * delta).exp()
        common = max(arb(1), pw(l0, -2 * d / L)) * max(arb(1), pw(k3, -3 * d / L))
        common *= max(pw(l0, 2 / (2 - s1)), pw(l0, 2 / (2 - s2)))
        common *= max(pw(k3, (2 * s1 - 1) / (2 - s1)), pw(k3, (2 * s2 - 1) / (2 - s2)))
        A1 = l1 * (d * (3 + 4 * L.log() / L)).exp() * common
        A1 *= max(arb(1), pw(k1, 3 * d / L))
        A1 *= max(pw(k1, 3 * (1 - s1) / (2 - s1)), pw(k1, 3 * (1 - s2) / (2 - s2)))
        A2 = l2 * (3 * d).exp() * common
        A2 *= max(arb(1), pw(k2, 3 * d / L))
        A2 *= max(pw(k2, 3 * (1 - s1) / (2 - s1)), pw(k2, 3 * (1 - s2) / (2 - s2)))
        B1 = (A1 + A2 * pw(L, -6 * (1 - s2) / (2 - s2))) / (2 * pi * d)
        return B1.upper(), k3

    def make_constants(prec):
        with ctx.workprec(prec):
            cs = constants()
            cs["c1"] = 6 / cs["pi"] ** 2 + q("1.048") / cs["logh"]
            cs["c2"] = (cs["pi"] * cs["m0"] * q(".62") / cs["logh"]
                        + 6 * cs["m0"] / (cs["pi"] * T0)
                        + cs["pi"] * cs["m0"] * q("1.048") / (arb(T0) * cs["logh"]))
            return cs

    def certify_full(row, cs, prec, bpub, printed_str, old_str=None):
        with ctx.workprec(prec):
            # Both non-dyadic input rationals must be constructed at target precision.
            b54a = arb(B54.numerator) / arb(B54.denominator)
            bpuba = arb(bpub.numerator) / arb(bpub.denominator)
            parts = row_parts(row, cs)
            k1, k2 = k1_k2(row, cs, parts)
            up54, k3_54 = b1_upper(row, cs, parts, k1, k2, b54a)
            up78, k3_78 = b1_upper(row, cs, parts, k1, k2, bpuba)
            printed = q(printed_str)
            m = (printed - up78).lower()
            # When b4 is unchanged, each logarithm in the ratio is log(1)=0:
            # the original code raised ZeroDivisionError (0/0) in lines 1..9.
            if bpub == B54:
                e_checked = None
                e_check_mode = "unchanged_b4_no_ratio"
            else:
                k_ratio = float(k3_78)/float(k3_54)
                b_ratio = float(up78)/float(up54)
                if not (k_ratio > 1.0 and b_ratio > 1.0):
                    raise ArithmeticError("unexpected nonpositive B1/K3 sensitivity")
                e_checked = math.log(b_ratio)/math.log(k_ratio)
                e_check_mode = "float_ratio_diagnostic_only"
            rec = {
                "workprec": prec,
                "k1": float(k1), "k2": float(k2),
                "k3": float(k3_78), "k3_at_0_54": float(k3_54),
                "B1_upper_0_54_upper": float(up54),
                "up_lo": float(up78.lower()), "up_hi": float(up78.upper()),
                "margin_lower": float(m),
                "gate_pass": bool(m > 0),
                "B1_upper_dyadic": dyadic(up78.upper()),
                "margin_lower_dyadic": dyadic(m),
                "e_exponent_checked": e_checked,
                "e_exponent_mode": e_check_mode,
            }
            if old_str is not None:
                rec["margin_table_old_lower"] = float((q(old_str) - up78).lower())
            return rec

    # ------------------------------------------------------------------
    # Plano de certificacao
    # ------------------------------------------------------------------
    cert_path = repo / "NUMERIA_P78_RAMARE_PUBLIC" / "certificate.json"
    cert = json.loads(cert_path.read_text())
    p78_b4 = {r["row"]: Fr(str(r["b4"])) for r in cert["rows"]}
    validate_table()

    patch = None
    if args.proposed:
        patch = json.loads(Path(args.proposed).read_text())
        if patch.get("schema") != "NUMERIA.P78.S2.TABLE.PATCH":
            raise SystemExit(f"patch com schema inesperado: {patch.get('schema')}")
    variant = args.variant or (patch or {}).get("default_variant", "robust")
    key = "b1_new_" + variant
    precs = sorted({int(p) for p in args.dyadic.split(",") if p.strip()})
    if not precs or min(precs) < 128:
        raise SystemExit("informe precisões >= 128 bits, por exemplo 192,256,384")

    plan = []
    if patch:
        for pr in patch["rows"]:
            idx = int(pr["row"])
            if idx not in p78_b4:
                raise SystemExit(f"linha {idx} do patch ausente no certificado P78")
            if key not in pr:
                raise SystemExit(f"campo {key} ausente na linha {idx} do patch")
            if Fr(pr["b4_public"]) != p78_b4[idx] or Fr(pr["delta"]) != Fr(TABLE[idx-1].delta):
                raise SystemExit(f"patch com b4 ou delta inconsistente na linha {idx}")
            plan.append((idx, str(pr[key]), "PROPOSED"))
    else:
        for idx in sorted(p78_b4):
            plan.append((idx, TABLE[idx - 1].b1, "TABLE"))

    rows_out = []
    for idx, printed_str, origem in plan:
        row = TABLE[idx - 1]
        bpub = p78_b4[idx]
        old_str = TABLE[idx - 1].b1 if origem == "PROPOSED" else None
        per_prec = {}
        for prec in precs:
            cs = make_constants(prec)
            per_prec[str(prec)] = certify_full(row, cs, prec, bpub, printed_str, old_str)
        rec0 = per_prec[str(precs[0])]

        exp_up = exp_ok = None
        if patch:
            pr = next(p for p in patch["rows"] if int(p["row"]) == idx)
            exp_up = pr.get("expected_upper_public_float")
            if exp_up is not None:
                exp_ok = abs(rec0["up_hi"] - float(exp_up)) / float(exp_up) < 1e-6
                if not exp_ok:
                    print(f"[aviso] linha {idx}: B1_upper (intervalo) {rec0['up_hi']:.12f} "
                          f"diverge da expectativa float do patch {float(exp_up):.12f}")

        e_c = e_closed(row.sigma_start)
        model_chk = rec0["e_exponent_checked"]
        model_ok = ((bpub == B54 and model_chk is None) or
                    (bpub != B54 and model_chk is not None and abs(model_chk-e_c) < 1e-9))
        # Compare exact dyadic Arb endpoints with 0.001 as a rational, never floats.
        def exact_dyad(pair):
            return Fr(int(pair[0])) * (Fr(2) ** int(pair[1]))
        robust_ok = (origem != "PROPOSED" or variant != "robust" or
                     all(exact_dyad(v["margin_lower_dyadic"]) >= Fr("0.001")
                         for v in per_prec.values()))
        rows_out.append({
            "row": idx, "delta": str(row.delta),
            "sigma_start": row.sigma_start,
            "b4_public": str(bpub),
            "printed": printed_str,
            "printed_variant": variant if origem == "PROPOSED" else "TABLE",
            "printed_source": origem,
            "table_b1_current": str(TABLE[idx - 1].b1),
            "table_already_matches": bool(str(TABLE[idx - 1].b1) == printed_str),
            "dyadic": {k: {"workprec": v["workprec"],
                           "k1": v["k1"], "k2": v["k2"],
                           "k3": v["k3"], "k3_at_0_54": v["k3_at_0_54"],
                           "B1_upper_lower": v["up_lo"], "B1_upper_upper": v["up_hi"],
                           "B1_upper_0_54_upper": v["B1_upper_0_54_upper"],
                           "margin_lower": v["margin_lower"],  # float diagnostic only
                           "B1_upper_endpoint_dyadic": v["B1_upper_dyadic"],
                           "margin_lower_dyadic": v["margin_lower_dyadic"],
                           "gate_pass": v["gate_pass"],
                           **({"margin_table_old_lower": v["margin_table_old_lower"]}
                              if "margin_table_old_lower" in v else {})}
                       for k, v in per_prec.items()},
            "gate_pass_all_precisions": all(v["gate_pass"] for v in per_prec.values()) and robust_ok,
            "robust_gate": robust_ok,
            "e_exponent_closed": e_c, "e_exponent_checked": model_chk,
            "model_consistent": bool(model_ok),
            "expected_upper_public_float": exp_up,
            "expected_upper_consistent": exp_ok,
        })

    fails = [r["row"] for r in rows_out if not r["gate_pass_all_precisions"]]
    summary = {
        "mode": "PROPOSED_PRE_ADOPTION" if patch else "TABLE_POST_ADOPTION",
        "variant": variant,
        "dyadic_precisions": precs,
        "rows_checked": [r["row"] for r in rows_out],
        "gate_fail_rows": fails,
        "all_gates_pass": not fails,
        "model_consistent_all": all(r["model_consistent"] for r in rows_out),
        "min_margin_lower": min(v for r in rows_out
                                for v in r["dyadic"].values()
                                if isinstance(v, dict) and "margin_lower" in v),
        "verdict": ("todos os gates passam em todas as precisoes; integrar ao "
                    "certificado oficial dyadic + CI + SHA256 e recertificar "
                    "downstream (p78_s2_downstream_scan.py)."
                    if not fails else
                    "GATE FALHOU nas linhas " + str(fails) +
                    (" — se a TABLE ainda tem os valores antigos do P77, use "
                     "--proposed (pre-adocao)." if not patch else "")),
    }

    sha = {}
    for label, p in [("table_py", repo / "NUMERIA_P76_COEFFICIENTS" / "table.py"),
                     ("audit_b23_py", repo / "NUMERIA_P76_COEFFICIENTS" / "audit_b23.py"),
                     ("certificate_json", cert_path),
                     ("patch_json", Path(args.proposed) if args.proposed else None)]:
        if p is not None and p.exists():
            sha[label] = {"path": str(p), "sha256": sha256_file(p)}

    doc = {
        "schema": SCHEMA,
        "mode": summary["mode"],
        "variant": variant,
        "dyadic_precisions": precs,
        "python": platform.python_version(),
        "flint_version": getattr(_flint, "__version__", "unknown"),
        "sha256_inputs": sha,
        "summary": summary,
        "rows": rows_out,
        "notes": [
            "Gate do P77: (printed - B1_upper(b4 publico)).lower() > 0, por precisao.",
            "Variante robust exige margem absoluta >= 1e-3 por construção do patch; "
            "a checagem final da margem robusta e feita no certificado oficial.",
            "Prototipo single-run: numeros no padrao do P77; certificado oficial e o "
            "do pipeline NUMERIA (dyadic completo + CI + SHA256).",
            "Apos adocao: recertificar downstream com p78_s2_downstream_scan.py.",
        ],
    }
    print(json.dumps(doc, indent=1))
    print()
    print(f"{'row':>3} {'printed':>8} {'prec':>5} {'B1_up_lower':>14} "
          f"{'margem_lower':>14} {'gate':>6}")
    for r in rows_out:
        for k, v in sorted(r["dyadic"].items()):
            print(f"{r['row']:>3} {r['printed']:>8} {k:>5} "
                  f"{v['B1_upper_lower']:>14.9f} {v['margin_lower']:>14.9f} "
                  f"{'PASS' if v['gate_pass'] else 'FAIL':>6}")
    print("VEREDITO:", summary["verdict"])

    dest = None
    candidates = ([Path(args.out)] if args.out else
                  [Path(__file__).resolve().parent / "p78_s2_recert_certificate.json",
                   repo.parent / "p78_s2_recert_certificate.json"])
    for cand in candidates:
        try:
            cand.parent.mkdir(parents=True, exist_ok=True)
            cand.write_text(json.dumps(doc, indent=1) + "\n")
            dest = cand
            break
        except OSError:
            continue
    if dest is None:
        raise SystemExit("não foi possível gravar relatório de diagnóstico")
    print(f"[ok] certificado-prototipo gravado em {dest}")
    if fails:
        raise SystemExit("S2 GATE FAIL: linhas " + str(fails))


if __name__ == "__main__":
    main()